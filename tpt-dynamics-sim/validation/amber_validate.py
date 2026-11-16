"""Check the converted Amber parameters against OpenMM (ADR-010, version 1).

For a residue sequence:

1. The Rust example `ffenergy` builds the chain from the converted ffparams file
   and prints its potential, positions and forces (Lennard-Jones and bonded terms).
2. OpenMM builds the same chain from the user's XML, with the same positions.
   Improper torsions are removed and charges are set to zero, so only the terms
   the converter covers are compared.
3. Potential and forces must agree to a tight relative tolerance.

The user's XML and the converted file are never committed. Run from the repository root:

    python tpt-dynamics-sim/validation/amber_validate.py \
        --ffxml <path>/protein.ff14SB.xml --residues ACE,ALA,NME
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
BOX = 40.0
GENERATED = HERE / "generated"


def run_ffenergy(params_path, residues):
    cmd = [
        "cargo", "run", "--release", "-q", "-p", "tpt-dynamics-sim", "--example", "ffenergy", "--",
        str(params_path), ",".join(residues), str(BOND_LENGTH), str(BOX),
    ]
    out = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, check=True).stdout.splitlines()
    potential = float(out[0].split()[1])
    n = int(out[1].split()[1])
    positions = [tuple(float(x) for x in line.split()) for line in out[2 : 2 + n]]
    assert out[2 + n] == "forces", "unexpected ffenergy output"
    forces = [tuple(float(x) for x in line.split()) for line in out[3 + n : 3 + 2 * n]]
    return potential, positions, forces


def openmm_reference(ffxml_paths, residues, positions, lj=True):
    types, res_defs, _ = {}, {}, None
    for p in ffxml_paths:
        t, r, _ = ac.parse_ffxml(p)
        types.update(t)
        res_defs.update(r)
    atoms, bonds, _ = ac.chain_topology(res_defs, residues)
    system = ac.build_openmm_system(ffxml_paths, res_defs, residues, atoms, bonds, types)

    graph = {(min(i, j), max(i, j)) for i, j in bonds}
    # Keep only proper torsions, matching what the converter writes.
    old_index = next(i for i, f in enumerate(system.getForces()) if isinstance(f, mm.PeriodicTorsionForce))
    old = system.getForce(old_index)
    new = mm.PeriodicTorsionForce()
    for k in range(old.getNumTorsions()):
        i, j, kk, l_, n, phase, kf = old.getTorsionParameters(k)
        if ac.is_proper((i, j, kk, l_), graph):
            new.addTorsion(i, j, kk, l_, n, phase, kf)
    system.removeForce(old_index)
    system.addForce(new)

    # Electrostatics are off on both sides for this check: the engine's `None` setting
    # has no unscreened Coulomb mode to match OpenMM's no-cutoff sum, and PME would need
    # matched Ewald parameters. So the charges are zeroed here, and the comparison covers
    # the bonded terms and Lennard-Jones. OpenMM's no-cutoff mode is used so that no
    # reaction-field terms appear.
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    for k in range(nb.getNumParticles()):
        _, sigma, eps = nb.getParticleParameters(k)
        nb.setParticleParameters(k, 0.0, sigma, eps if lj else 0.0)
    for k in range(nb.getNumExceptions()):
        i, j, _, sigma, eps = nb.getExceptionParameters(k)
        nb.setExceptionParameters(k, i, j, 0.0, sigma, eps if lj else 0.0)

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
    ap.add_argument("--tolerance", type=float, default=1e-9)
    args = ap.parse_args(argv)
    residues = [r.strip() for r in args.residues.split(",") if r.strip()]

    GENERATED.mkdir(parents=True, exist_ok=True)
    params = GENERATED / f"{args.name}.toml"
    info = ac.convert(args.ffxml, residues, params, args.name)
    print("converted:", info)

    rust_e, positions, rust_f = run_ffenergy(params, residues)
    om_e, om_f, n = openmm_reference(args.ffxml, residues, positions)
    assert n == len(positions) == len(rust_f) == len(om_f)

    rel_e = abs(rust_e - om_e) / max(abs(om_e), 1.0)
    scale = max(max(abs(c) for f in om_f for c in f), 1.0)
    worst_f = max(max(abs(a - b) for a, b in zip(rf, of)) for rf, of in zip(rust_f, om_f))
    rel_f = worst_f / scale
    print(f"atoms: {n}")
    print(f"potential: rust {rust_e:.9f}  openmm {om_e:.9f}  relative {rel_e:.3e}")
    print(f"forces: largest {scale:.3f}  worst component difference {worst_f:.3e}  relative {rel_f:.3e}")
    ok = rel_e < args.tolerance and rel_f < args.tolerance
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
