"""Independent OpenMM reference for residue-chain energies.

The chain geometry, the bonded topology (bonds, angles, proper dihedrals, 1-4
pairs, exclusions) and the parameter lookup are implemented here, separately
from the Rust code, and the energy is evaluated by OpenMM. The Rust test
tests/openmm_reference.rs builds the same chains with ParamSet::assemble_chains
and checks that positions and total potential agree.

Inputs are the synthetic test parameters (data/synthetic-test.toml, ours) and
the residue templates (validation/chain-residues.toml, ours). No third-party
parameter data is read or written.

Output: validation/fixtures/<name>.txt, a plain-text file with the energy and the
positions, so the Rust test needs no extra dependency:

    potential <kJ/mol>
    components <bond> <angle> <torsion> <nonbonded>
    atoms <n>
    x y z           (one line per atom, nm)

Run from the repository root:

    python tpt-dynamics-sim/validation/openmm_chain_reference.py
"""

from __future__ import annotations

import math
import pathlib
import sys
import tomllib

import openmm as mm

HERE = pathlib.Path(__file__).resolve().parent
SYNTHETIC = HERE.parent / "data" / "synthetic-test.toml"
TEMPLATES = HERE / "chain-residues.toml"
FIXTURES = HERE / "fixtures"

BOND_LENGTH = 0.15  # nm, matches the Rust test
ORIGIN = (0.5, 0.5, 0.5)
BOX = 12.0  # nm; the cutoff is well below half of it, so no periodic images matter
CUTOFF = 4.0  # nm, no shift and no switch, to match Modifier::None

# name -> residue list for each fixture
CASES = {
    "rsa5": ["RSA"] * 5,
    "lip3": ["LIP"] * 3,
    "mixed": ["RSA", "RSA", "LIP"],
}


def load_inputs():
    synth = tomllib.loads(SYNTHETIC.read_text(encoding="utf-8"))
    templates = tomllib.loads(TEMPLATES.read_text(encoding="utf-8"))
    atomtypes = {a["name"]: a for a in synth["atomtype"]}
    molecules = {m["name"]: m for m in templates["molecule"]}
    return synth, atomtypes, molecules


def build_chain(molecules, residues):
    """Places residues head to tail, using the same rule as the Rust builder:
    each residue after the first is translated so that its head lies
    BOND_LENGTH along +x from the previous tail."""
    positions, types, bonds = [], [], []
    link = None  # (global tail index, tail world position)
    for k, name in enumerate(residues):
        m = molecules[name]
        base = len(positions)
        if link is None:
            shift = ORIGIN
        else:
            tail_index, tail_world = link
            head_ref = m["atoms"][m["head"]]["position"]
            shift = tuple(
                tail_world[i] + (BOND_LENGTH if i == 0 else 0.0) - head_ref[i] for i in range(3)
            )
            bonds.append((tail_index, base + m["head"]))
        for atom in m["atoms"]:
            positions.append(tuple(atom["position"][i] + shift[i] for i in range(3)))
            types.append(atom["type"])
        bonds.extend((base + i, base + j) for i, j in m.get("bonds", []))
        if k + 1 < len(residues):
            t = m["tail"]
            link = (base + t, positions[base + t])
    return positions, types, bonds


def topology(n, bonds):
    """Angles, proper dihedrals, 1-4 pairs and exclusions from the bond graph."""
    adj = [set() for _ in range(n)]
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)

    # Shortest path lengths up to 3 bonds.
    def distances(src):
        dist = {src: 0}
        frontier = [src]
        for d in (1, 2, 3):
            nxt = []
            for a in frontier:
                for b in adj[a]:
                    if b not in dist:
                        dist[b] = d
                        nxt.append(b)
            frontier = nxt
        return dist

    exclusions, pairs14 = [], []
    for i in range(n):
        for j, d in distances(i).items():
            if j > i:
                exclusions.append((i, j))
                if d == 3:
                    pairs14.append((i, j))

    angles = []
    for c in range(n):
        nb = sorted(adj[c])
        for x in range(len(nb)):
            for y in range(x + 1, len(nb)):
                angles.append((nb[x], c, nb[y]))

    dihedrals, seen = [], set()
    for j, k in bonds:
        for a, b in ((j, k), (k, j)):
            for i in adj[a]:
                if i == b:
                    continue
                for l_ in adj[b]:
                    if l_ in (a, i):
                        continue
                    path = (i, a, b, l_)
                    key = min(path, path[::-1])
                    if key not in seen:
                        seen.add(key)
                        dihedrals.append(path)
    return exclusions, pairs14, angles, dihedrals


def match(pattern_types, types):
    """Forward or reverse match with '*' as a wildcard."""
    def ok(p, t):
        return all(a == "*" or a == b for a, b in zip(p, t))
    return ok(pattern_types, types) or ok(pattern_types, types[::-1])


def build_system(synth, atomtypes, types, bonds, exclusions, pairs14, angles, dihedrals):
    fudge_lj = float(synth["fudge_lj"])
    system = mm.System()
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
    nb.setCutoffDistance(CUTOFF)
    for t in types:
        at = atomtypes[t]
        system.addParticle(at["mass"])
        nb.addParticle(0.0, at["sigma"], at["epsilon"])

    pair14 = set(pairs14)
    for i, j in exclusions:
        ti, tj = atomtypes[types[i]], atomtypes[types[j]]
        sigma = 0.5 * (ti["sigma"] + tj["sigma"])
        eps = fudge_lj * math.sqrt(ti["epsilon"] * tj["epsilon"]) if (i, j) in pair14 else 0.0
        nb.addException(i, j, 0.0, sigma, eps)
    nb.setForceGroup(3)
    system.addForce(nb)

    bond_types = synth["bondtype"]
    bf = mm.HarmonicBondForce()
    for i, j in bonds:
        entry = next(b for b in bond_types if sorted(b["types"]) == sorted([types[i], types[j]]))
        bf.addBond(i, j, entry["r0"], entry["k"])
    bf.setForceGroup(0)
    system.addForce(bf)

    af = mm.HarmonicAngleForce()
    for i, j, k in angles:
        entry = next(a for a in synth["angletype"] if match(a["types"], [types[i], types[j], types[k]]))
        af.addAngle(i, j, k, math.radians(entry["theta0"]), entry["k"])
    af.setForceGroup(1)
    system.addForce(af)

    tf = mm.PeriodicTorsionForce()
    for i, j, k, l_ in dihedrals:
        entry = next(
            d for d in synth["dihedraltype"] if match(d["types"], [types[i], types[j], types[k], types[l_]])
        )
        for term in entry["terms"]:
            tf.addTorsion(i, j, k, l_, term["n"], math.radians(term["phase"]), term["k"])
    tf.setForceGroup(2)
    system.addForce(tf)
    return system


def energies(system, positions):
    platform = mm.Platform.getPlatformByName("Reference")
    context = mm.Context(system, mm.VerletIntegrator(0.001), platform)
    context.setPositions([mm.Vec3(*p) for p in positions])
    out = {}
    total = context.getState(getEnergy=True).getPotentialEnergy()
    out["potential"] = total.value_in_unit(mm.unit.kilojoule_per_mole)
    for name, group in (("bond", 0), ("angle", 1), ("torsion", 2), ("nonbonded", 3)):
        e = context.getState(getEnergy=True, groups={group}).getPotentialEnergy()
        out[name] = e.value_in_unit(mm.unit.kilojoule_per_mole)
    return out


def write_fixture(name, positions, e):
    FIXTURES.mkdir(parents=True, exist_ok=True)
    lines = [
        f"potential {e['potential']!r}",
        f"components {e['bond']!r} {e['angle']!r} {e['torsion']!r} {e['nonbonded']!r}",
        f"atoms {len(positions)}",
    ]
    lines += [f"{x!r} {y!r} {z!r}" for x, y, z in positions]
    path = FIXTURES / f"{name}.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    synth, atomtypes, molecules = load_inputs()
    for name, residues in CASES.items():
        positions, types, bonds = build_chain(molecules, residues)
        exclusions, pairs14, angles, dihedrals = topology(len(positions), bonds)
        system = build_system(synth, atomtypes, types, bonds, exclusions, pairs14, angles, dihedrals)
        e = energies(system, positions)
        path = write_fixture(name, positions, e)
        print(
            f"{name}: atoms={len(positions)} bonds={len(bonds)} angles={len(angles)} "
            f"dihedrals={len(dihedrals)} pairs14={len(pairs14)} "
            f"potential={e['potential']:.9f} -> {path.relative_to(HERE.parent.parent)}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
