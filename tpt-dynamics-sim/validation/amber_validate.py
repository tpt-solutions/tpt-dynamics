"""Check the converted Amber parameters against OpenMM (ADR-010).

For a residue sequence, the Rust example `ffenergy` builds the chain from the
converted ffparams file, and OpenMM builds the same chain from the user's XML
with the same positions. Both potentials and forces must agree.

Two modes:

- default (bonded and Lennard-Jones): electrostatics are off on both sides. The
  charges are zeroed in OpenMM, and the engine's `None` setting adds no Coulomb.
  OpenMM's no-cutoff mode is used, so no reaction-field terms appear.
- `--pme TOL ORDER`: electrostatics are on, with PME, in a periodic box. The engine
  reports the Ewald splitting alpha and FFT grid it chose, and OpenMM is given the
  same values. The B-spline order must be 5, which is OpenMM's fixed order.

Impropers are included in both modes; the converter writes them as periodic improper types.
The user's XML and the converted file are never committed. Run from the repo root:

    python tpt-dynamics-sim/validation/amber_validate.py \
        --ffxml <path>/protein.ff14SB.xml --residues ACE,ALA,NME [--pme 1e-5 5]
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import openmm as mm
from openmm import unit

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import amber_convert as ac  # noqa: E402

REPO = HERE.parent.parent
BOND_LENGTH = 0.15
BOX_NO_PME = 40.0
BOX_PME = 9.0  # must exceed 2 * (cutoff + skin) = 8.4 nm
CUTOFF = 4.0
GENERATED = HERE / "generated"


def run_ffenergy(params_path, residues, box, pme=None):
    """Runs the Rust example. Returns the keyed header values, positions and forces."""
    cmd = [
        "cargo", "run", "--release", "-q", "-p", "tpt-dynamics-sim", "--example", "ffenergy", "--",
        str(params_path), ",".join(residues), str(BOND_LENGTH), str(box),
    ]
    if pme is not None:
        cmd += [str(pme[0]), str(pme[1])]
    lines = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, check=True).stdout.splitlines()

    header = {}
    i = 0
    while lines[i] not in ("forces",) and lines[i].split()[0] in ("potential", "atoms", "pme_alpha", "pme_grid"):
        key, *rest = lines[i].split()
        if key == "pme_grid":
            header[key] = tuple(int(x) for x in rest)
        elif key == "atoms":
            header[key] = int(rest[0])
        else:
            header[key] = float(rest[0])
        i += 1
    n = header["atoms"]
    positions = [tuple(float(x) for x in line.split()) for line in lines[i : i + n]]
    i += n
    assert lines[i] == "forces", "unexpected ffenergy output"
    forces = [tuple(float(x) for x in line.split()) for line in lines[i + 1 : i + 1 + n]]
    return header, positions, forces


def load_templates(ffxml_paths):
    types, res_defs = {}, {}
    for p in ffxml_paths:
        t, r, _ = ac.parse_ffxml(p)
        types.update(t)
        res_defs.update(r)
    return types, res_defs


def openmm_reference(ffxml_paths, residues, positions, box=None, pme=None, lj=True):
    types, res_defs = load_templates(ffxml_paths)
    atoms, bonds, _ = ac.chain_topology(res_defs, residues)
    system = ac.build_openmm_system(ffxml_paths, res_defs, residues, atoms, bonds, types)

    # All torsions are kept, impropers included: the converter writes both kinds.

    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    if pme is None:
        # Electrostatics off: zero the charges, no-cutoff mode.
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        for k in range(nb.getNumParticles()):
            _, sigma, eps = nb.getParticleParameters(k)
            nb.setParticleParameters(k, 0.0, sigma, eps if lj else 0.0)
        for k in range(nb.getNumExceptions()):
            i, j, _, sigma, eps = nb.getExceptionParameters(k)
            nb.setExceptionParameters(k, i, j, 0.0, sigma, eps if lj else 0.0)
    else:
        alpha, (nx, ny, nz) = pme
        nb.setNonbondedMethod(mm.NonbondedForce.PME)
        nb.setCutoffDistance(CUTOFF * unit.nanometer)
        nb.setPMEParameters(alpha / unit.nanometer, nx, ny, nz)
        nb.setUseDispersionCorrection(False)
        if not lj:
            for k in range(nb.getNumParticles()):
                q, sigma, _ = nb.getParticleParameters(k)
                nb.setParticleParameters(k, q, sigma, 0.0)
            for k in range(nb.getNumExceptions()):
                i, j, qq, sigma, _ = nb.getExceptionParameters(k)
                nb.setExceptionParameters(k, i, j, qq, sigma, 0.0)
        system.setDefaultPeriodicBoxVectors(mm.Vec3(box, 0, 0), mm.Vec3(0, box, 0), mm.Vec3(0, 0, box))

    context = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    context.setPositions([mm.Vec3(*p) for p in positions])
    state = context.getState(getEnergy=True, getForces=True)
    potential = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = [(f.x, f.y, f.z) for f in state.getForces().value_in_unit(unit.kilojoule_per_mole / unit.nanometer)]
    return potential, forces, len(atoms)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ffxml", nargs="+", required=True, type=pathlib.Path)
    ap.add_argument("--residues", required=True)
    ap.add_argument("--name", default="amber-validation")
    ap.add_argument("--pme", nargs=2, metavar=("TOL", "ORDER"), help="check electrostatics with PME")
    ap.add_argument("--tolerance", type=float, default=1e-9)
    ap.add_argument("--no-lj", action="store_true",
                    help="switch Lennard-Jones off on both sides, to check electrostatics alone")
    args = ap.parse_args(argv)
    residues = [r.strip() for r in args.residues.split(",") if r.strip()]

    GENERATED.mkdir(parents=True, exist_ok=True)
    params = GENERATED / f"{args.name}.toml"
    info = ac.convert(args.ffxml, residues, params, args.name)
    print("converted:", info)
    lj = not args.no_lj
    molecule = "-".join(residues)  # the converter writes the whole sequence as one molecule
    if not lj:
        import re
        text = params.read_text(encoding="utf-8")
        params = GENERATED / f"{args.name}-nolj.toml"
        params.write_text(re.sub(r"epsilon = [^\n]+", "epsilon = 0.0", text), encoding="utf-8")

    if args.pme is None:
        box = BOX_NO_PME
        header, positions, rust_f = run_ffenergy(params, [molecule], box)
        om_e, om_f, n = openmm_reference(args.ffxml, residues, positions, lj=lj)
        mode = "bonded + Lennard-Jones (electrostatics off)"
    else:
        tol, order = float(args.pme[0]), int(args.pme[1])
        if order != 5:
            raise SystemExit("OpenMM's B-spline order is 5; use --pme TOL 5")
        box = BOX_PME
        header, positions, rust_f = run_ffenergy(params, [molecule], box, pme=(tol, order))
        alpha, grid = header["pme_alpha"], header["pme_grid"]
        print(f"pme: alpha {alpha:.6f} 1/nm, grid {grid}")
        om_e, om_f, n = openmm_reference(args.ffxml, residues, positions, box=box, pme=(alpha, grid), lj=lj)
        mode = f"PME (tolerance {tol:g}, order {order})" + ("" if lj else ", Lennard-Jones off")
    rust_e = header["potential"]
    assert n == len(positions) == len(rust_f) == len(om_f)

    rel_e = abs(rust_e - om_e) / max(abs(om_e), 1.0)
    scale = max(max(abs(c) for f in om_f for c in f), 1.0)
    worst_f = max(max(abs(a - b) for a, b in zip(rf, of)) for rf, of in zip(rust_f, om_f))
    rel_f = worst_f / scale
    print(f"mode: {mode}")
    print(f"atoms: {n}")
    print(f"potential: rust {rust_e:.9f}  openmm {om_e:.9f}  relative {rel_e:.3e}")
    print(f"forces: largest {scale:.3f}  worst component difference {worst_f:.3e}  relative {rel_f:.3e}")
    ok = rel_e < args.tolerance and rel_f < args.tolerance
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
