"""Convert user-supplied Amber force-field XML into tpt-dynamics ffparams TOML.

This is the first version of the converter described in ADR-010. It covers
bonds, angles, proper torsions and Lennard-Jones. Impropers and Urey-Bradley
terms are not converted yet, and the converter says so in its output.

Division of work:

- The residue templates (atoms, charges, bonds, external bonds) and the atom
  types are read from the XML the user passes in.
- OpenMM builds the System for the requested residue sequence from that XML.
  OpenMM resolves Amber's wildcard and ordering rules, so this tool does not
  reimplement them.
- The bond, angle, torsion and Lennard-Jones terms are read back from the System
  and written as explicit entries, one per type combination.

Nothing from the input files is bundled with the repository. The output belongs
to the user and goes to a path they choose (validation/generated/ is git-ignored).

Usage (from the repository root):

    python tpt-dynamics-sim/validation/amber_convert.py \
        --ffxml <path>/protein.ff14SB.xml --residues ACE,ALA,NME \
        --out tpt-dynamics-sim/validation/generated/dipeptide.toml

The geometry written into the templates is a deterministic tree embedding, not a
real structure. It is only used to give the energy evaluation positions, and
the validation step does not rely on it being realistic.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import sys
import xml.etree.ElementTree as ET

import openmm as mm
from openmm import app, unit

BOND_LENGTH = 0.15  # nm; used for the tree embedding and for head-to-tail links
TETRAHEDRAL = math.radians(109.5)


# ---------------------------------------------------------------- XML parsing

def parse_ffxml(path: pathlib.Path):
    root = ET.parse(path).getroot()
    types = {}
    for t in root.iter("Type"):
        types[t.attrib["name"]] = {
            "mass": float(t.attrib["mass"]),
            "element": t.attrib["element"],
        }
    residues = {}
    for res in root.iter("Residue"):
        name = res.attrib["name"]
        atoms = [
            {"name": a.attrib["name"], "type": a.attrib["type"], "charge": float(a.attrib["charge"])}
            for a in res.findall("Atom")
        ]
        bonds = []
        for b in res.findall("Bond"):
            bonds.append((b.attrib["atomName1"], b.attrib["atomName2"]))
        external = [e.attrib["atomName"] for e in res.findall("ExternalBond")]
        residues[name] = {"atoms": atoms, "bonds": bonds, "external": external}
    nb = root.find("NonbondedForce")
    if nb is None:
        raise SystemExit(f"{path}: no NonbondedForce section")
    scales = {
        "fudge_lj": float(nb.attrib.get("lj14scale", "0.5")),
        "fudge_qq": float(nb.attrib.get("coulomb14scale", "0.8333333333")),
    }
    return types, residues, scales


# ------------------------------------------------------------ residue chains

def residue_atoms(residues, name):
    res = residues[name]
    names = [a["name"] for a in res["atoms"]]
    index = {n: i for i, n in enumerate(names)}
    bonds = [(index[a], index[b]) for a, b in res["bonds"]]
    return res["atoms"], bonds, index


def head_and_tail(residues, name):
    """Head is the external atom named N, tail the one named C (Amber convention)."""
    ext = residues[name]["external"]
    head = "N" if "N" in ext else None
    tail = "C" if "C" in ext else None
    return head, tail


def chain_topology(residues, sequence):
    """Global atom list, bonds (including head-to-tail links) and residue offsets."""
    atoms, bonds, offsets = [], [], []
    link = None
    for name in sequence:
        res_atoms, res_bonds, index = residue_atoms(residues, name)
        base = len(atoms)
        offsets.append(base)
        atoms.extend((name, a) for a in res_atoms)
        bonds.extend((base + i, base + j) for i, j in res_bonds)
        head, tail = head_and_tail(residues, name)
        if link is not None:
            if head is None:
                raise SystemExit(f"residue {name} follows another residue but has no N external bond")
            bonds.append((link, base + index[head]))
        link = None
        if tail is not None:
            link = base + index[tail]
    return atoms, bonds, offsets


def embed(n, bonds):
    """Deterministic tree embedding with bond length BOND_LENGTH and tetrahedral-ish
    angles. Used only to give the energy evaluation positions."""
    adj = [[] for _ in range(n)]
    for i, j in bonds:
        adj[i].append(j)
        adj[j].append(i)
    pos = [None] * n
    parent = [None] * n
    pos[0] = (0.0, 0.0, 0.0)
    order = [0]
    seen = {0}
    k = 0
    while k < len(order):
        p = order[k]
        k += 1
        child_no = 0
        for c in adj[p]:
            if c in seen:
                continue
            seen.add(c)
            parent[c] = p
            pos[c] = place(pos, parent, p, child_no)
            child_no += 1
            order.append(c)
    if any(q is None for q in pos):
        raise SystemExit("residue graph is not connected; cannot embed")
    return pos


def _norm(v):
    m = math.sqrt(sum(x * x for x in v))
    return tuple(x / m for x in v)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def place(pos, parent, p, child_no):
    g = parent[p]
    if g is None:
        dirs = [(1.0, 0.0, 0.0), (-0.5, math.sqrt(3) / 2, 0.0), (-0.5, -math.sqrt(3) / 2, 0.0), (0.0, 0.0, 1.0)]
        d = dirs[child_no % len(dirs)]
        return tuple(pos[p][i] + BOND_LENGTH * d[i] for i in range(3))
    e = _norm(tuple(g_i - p_i for g_i, p_i in zip(pos[g], pos[p])))  # p -> grandparent
    axis = _norm(_cross(e, (0.0, 0.0, 1.0))) if abs(e[2]) < 0.9 else _norm(_cross(e, (1.0, 0.0, 0.0)))
    perp = _cross(e, axis)
    phi = child_no * 2.0944  # spread siblings 120 degrees apart around the bond
    w = tuple(math.cos(phi) * axis[i] + math.sin(phi) * perp[i] for i in range(3))
    dvec = tuple(math.cos(TETRAHEDRAL) * e[i] + math.sin(TETRAHEDRAL) * w[i] for i in range(3))
    return tuple(pos[p][i] + BOND_LENGTH * dvec[i] for i in range(3))


# ------------------------------------------------------------- OpenMM System

def build_openmm_system(ffxml_paths, residues, sequence, atoms, bonds, types):
    ff = app.ForceField(*[str(p) for p in ffxml_paths])
    top = app.Topology()
    chain = top.addChain()
    res_objs = []
    for name in sequence:
        res_objs.append(top.addResidue(name, chain))
    atom_objs = []
    # Atoms are added residue by residue, in the same order as `atoms`.
    for r_i, name in enumerate(sequence):
        for a in residues[name]["atoms"]:
            elem = app.element.get_by_symbol(types[a["type"]]["element"])
            atom_objs.append(top.addAtom(a["name"], elem, res_objs[r_i]))
    for i, j in bonds:
        top.addBond(atom_objs[i], atom_objs[j])
    system = ff.createSystem(
        top,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=4.0 * unit.nanometer,
        constraints=None,
        removeCMMotion=False,
    )
    return system


def openmm_terms(system):
    """Read bonds, angles, torsions and per-particle Lennard-Jones from the System."""
    bonds, angles, torsions, nb_particles = [], [], [], {}
    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            for k in range(force.getNumBonds()):
                i, j, r0, kk = force.getBondParameters(k)
                bonds.append((i, j, r0.value_in_unit(unit.nanometer), kk.value_in_unit(unit.kilojoule_per_mole / unit.nanometer**2)))
        elif isinstance(force, mm.HarmonicAngleForce):
            for k in range(force.getNumAngles()):
                i, j, kk_, theta0, kf = force.getAngleParameters(k)
                angles.append((i, j, kk_, theta0.value_in_unit(unit.radian), kf.value_in_unit(unit.kilojoule_per_mole / unit.radian**2)))
        elif isinstance(force, mm.PeriodicTorsionForce):
            for k in range(force.getNumTorsions()):
                i, j, kk_, l_, n, phase, kf = force.getTorsionParameters(k)
                torsions.append((i, j, kk_, l_, n, phase.value_in_unit(unit.radian), kf.value_in_unit(unit.kilojoule_per_mole)))
        elif isinstance(force, mm.NonbondedForce):
            for k in range(force.getNumParticles()):
                charge, sigma, eps = force.getParticleParameters(k)
                nb_particles[k] = (sigma.value_in_unit(unit.nanometer), eps.value_in_unit(unit.kilojoule_per_mole))
    return bonds, angles, torsions, nb_particles


def is_proper(quad, graph_bonds):
    a, b, c, d = quad
    return (min(a, b), max(a, b)) in graph_bonds and (min(b, c), max(b, c)) in graph_bonds and (min(c, d), max(c, d)) in graph_bonds


# ----------------------------------------------------------------- TOML output

def fmt(x):
    return repr(float(x))


def write_toml(path, name, types, used_types, sigma_eps, bond_types, angle_types, dihedral_types, molecules, scales):
    lines = [
        "# Generated by tpt-dynamics-sim/validation/amber_convert.py (ADR-010, version 1).",
        "# Derived from the user's own Amber force-field XML. Do not commit this file.",
        "# Not included in v1: impropers, Urey-Bradley terms, electrostatics are passed through as charges.",
        "",
        f'name = "{name}"',
        "format_version = 1",
        f"fudge_lj = {fmt(scales['fudge_lj'])}",
        f"fudge_qq = {fmt(scales['fudge_qq'])}",
        "",
    ]
    for t in sorted(used_types):
        sigma, eps = sigma_eps[t]
        lines += [
            "[[atomtype]]",
            f'name = "{t}"',
            f"mass = {fmt(types[t]['mass'])}",
            f"sigma = {fmt(sigma)}",
            f"epsilon = {fmt(eps)}",
            "",
        ]
    for key, (k, r0) in sorted(bond_types.items()):
        lines += ["[[bondtype]]", f"types = [{', '.join(chr(34) + t + chr(34) for t in key)}]", f"k = {fmt(k)}", f"r0 = {fmt(r0)}", ""]
    for key, (k, theta0) in sorted(angle_types.items()):
        lines += ["[[angletype]]", f"types = [{', '.join(chr(34) + t + chr(34) for t in key)}]", f"k = {fmt(k)}", f"theta0 = {fmt(math.degrees(theta0))}", ""]
    for key, terms in sorted(dihedral_types.items()):
        term_text = ", ".join(
            f"{{ k = {fmt(kf)}, n = {n}, phase = {fmt(math.degrees(phase))} }}" for kf, n, phase in terms
        )
        lines += ["[[dihedraltype]]", f"types = [{', '.join(chr(34) + t + chr(34) for t in key)}]", f"terms = [{term_text}]", ""]
    for mol in molecules:
        lines += ["[[molecule]]", f'name = "{mol["name"]}"', "atoms = ["]
        for a in mol["atoms"]:
            x, y, z = a["position"]
            lines.append(
                f'  {{ name = "{a["name"]}", type = "{a["type"]}", charge = {fmt(a["charge"])}, '
                f"position = [{fmt(x)}, {fmt(y)}, {fmt(z)}] }},"
            )
        lines.append("]")
        bond_text = ", ".join(f"[{i}, {j}]" for i, j in mol["bonds"])
        lines.append(f"bonds = [{bond_text}]")
        if mol.get("head") is not None:
            lines.append(f"head = {mol['head']}")
        if mol.get("tail") is not None:
            lines.append(f"tail = {mol['tail']}")
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ------------------------------------------------------------------- driver

def convert(ffxml_paths, sequence, out_path, name):
    types, residues, scales = None, None, None
    types, residues, scales = {}, {}, {}
    for p in ffxml_paths:
        t, r, s = parse_ffxml(p)
        types.update(t)
        residues.update(r)
        scales = s
    for res in sequence:
        if res not in residues:
            raise SystemExit(f"residue {res} is not in the given XML files")

    atoms, bonds, offsets = chain_topology(residues, sequence)
    positions = embed(len(atoms), bonds)
    system = build_openmm_system(ffxml_paths, residues, sequence, atoms, bonds, types)
    om_bonds, om_angles, om_torsions, nb = openmm_terms(system)

    graph = {(min(i, j), max(i, j)) for i, j in bonds}
    # Atom types per global index, from the residue templates.
    atom_type = [a["type"] for _, a in atoms]

    # Bonds: canonical by type pair.
    bond_types = {}
    for i, j, r0, k in om_bonds:
        key = tuple(sorted((atom_type[i], atom_type[j])))
        val = (k, r0)
        if key in bond_types and not all(abs(x - y) < 1e-9 for x, y in zip(bond_types[key], val)):
            raise SystemExit(f"bond type {key} has conflicting parameters")
        bond_types[key] = val

    # Angles: canonical by type triple.
    angle_types = {}
    for i, j, k_, theta0, k in om_angles:
        fwd = (atom_type[i], atom_type[j], atom_type[k_])
        key = min(fwd, fwd[::-1])
        val = (k, theta0)
        if key in angle_types and not all(abs(x - y) < 1e-9 for x, y in zip(angle_types[key], val)):
            raise SystemExit(f"angle type {key} has conflicting parameters")
        angle_types[key] = val

    # Proper torsions: group the OpenMM terms by atom quadruple, skip impropers.
    by_quad = {}
    skipped = 0
    for i, j, k_, l_, n, phase, kf in om_torsions:
        quad = (i, j, k_, l_)
        if not is_proper(quad, graph):
            skipped += 1
            continue
        key_atoms = min(quad, quad[::-1])
        by_quad.setdefault(key_atoms, []).append((kf, n, phase))
    dihedral_types = {}
    for quad, terms in by_quad.items():
        tkey = tuple(atom_type[a] for a in quad)
        tkey = min(tkey, tkey[::-1])
        terms = sorted(terms)
        if tkey in dihedral_types and dihedral_types[tkey] != terms:
            raise SystemExit(f"dihedral type {tkey} has conflicting terms")
        dihedral_types[tkey] = terms
    # Our topology enumerates every 3-bond path. A path with no OpenMM torsion
    # gets a zero-force entry, so the engine can still find it.
    for i, j, k_, l_ in enumerate_paths(graph, len(atoms)):
        tkey = tuple(atom_type[a] for a in (i, j, k_, l_))
        tkey = min(tkey, tkey[::-1])
        dihedral_types.setdefault(tkey, [(0.0, 1, 0.0)])

    sigma_eps = {}
    for g, (sigma, eps) in nb.items():
        t = atom_type[g]
        if t in sigma_eps and not all(abs(x - y) < 1e-9 for x, y in zip(sigma_eps[t], (sigma, eps))):
            raise SystemExit(f"atom type {t} has conflicting Lennard-Jones parameters")
        sigma_eps[t] = (sigma, eps)

    used_types = set(atom_type)
    molecules = []
    for r_i, res in enumerate(sequence):
        base = offsets[r_i]
        res_atoms, res_bonds, index = residue_atoms(residues, res)
        head, tail = head_and_tail(residues, res)
        mol_atoms = []
        for local, a in enumerate(res_atoms):
            g = base + local
            mol_atoms.append(
                {"name": a["name"], "type": a["type"], "charge": a["charge"], "position": positions[g]}
            )
        molecules.append(
            {
                "name": res,
                "atoms": mol_atoms,
                "bonds": res_bonds,
                "head": index[head] if head else None,
                "tail": index[tail] if tail else None,
            }
        )
    write_toml(out_path, name, types, used_types, sigma_eps, bond_types, angle_types, dihedral_types, molecules, scales)
    return {
        "atoms": len(atoms),
        "bonds": len(bonds),
        "bond_types": len(bond_types),
        "angle_types": len(angle_types),
        "dihedral_types": len(dihedral_types),
        "skipped_impropers": skipped,
        "path": str(out_path),
    }


def enumerate_paths(graph, n):
    adj = [set() for _ in range(n)]
    for i, j in graph:
        adj[i].add(j)
        adj[j].add(i)
    seen = set()
    for j, k in graph:
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
                        yield path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ffxml", nargs="+", required=True, type=pathlib.Path, help="Amber force-field XML file(s)")
    ap.add_argument("--residues", required=True, help="comma-separated residue sequence, e.g. ACE,ALA,NME")
    ap.add_argument("--out", required=True, type=pathlib.Path, help="output ffparams TOML path")
    ap.add_argument("--name", default="amber-converted", help="name stored in the output")
    args = ap.parse_args(argv)
    sequence = [r.strip() for r in args.residues.split(",") if r.strip()]
    info = convert(args.ffxml, sequence, args.out, args.name)
    for k, v in info.items():
        print(f"{k}: {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
