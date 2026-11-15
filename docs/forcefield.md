# Force-field parameter sets

## Policy

- tpt-dynamics **does not bundle** any third-party force-field parameter set.
- Users supply their own parameter files in the project's TOML format (see the
  module docs in `tpt-dynamics-sim/src/ffparams.rs`). The only file shipped is
  `tpt-dynamics-sim/data/synthetic-test.toml`, written for the tests.
- Loading a set is the user's act. **Whoever loads a parameter set must follow that
  set's licence and citation terms.** The project does not grant rights to any set.

## Licence status of common sets (checked 2026-11-16, re-check before use)

| Set | Status | Bundled? |
|-----|--------|----------|
| OpenFF Sage / Parsley (small molecules) | [openff-forcefields](https://github.com/openforcefield/openff-forcefields) states CC-BY-4.0. Format is SMIRNOFF `.offxml`, not this project's format. | No |
| AMBER | [AmberTools](https://www.ambermd.org/GetAmber.php) is free; components are mostly GPL. Terms for the parameter files were not verified against official docs. | No |
| CHARMM36 | No explicit licence found. The lab asks for citation. | No |
| GROMOS 54A7 | Licence not found. Check `forcefield.itp` in a local GROMACS install. | No |

Revisit this table before any change to the policy above.

## Why

Parameter sets are data, and each has its own terms. Bundling one would make the
project's distributions (crates, wheels) carry those terms. Keeping sets out keeps
the project's MIT OR Apache-2.0 licence unaffected.

## Residue chains (proteins, nucleic acids, lipid tails)

`ParamSet::assemble_chains` builds linear chains of residue templates. A
template may name a `head` atom (bonds to the previous residue) and a `tail`
atom (bonds to the next one). Each link is a bond between the previous tail and
the next head. Angles, dihedrals and 1-4 pairs come from the whole bond graph, so
they span residue boundaries. Several chains can share one system, for example
a membrane of lipid chains.

Geometry rules:

- Each residue after the first is translated so that its head lies
  `bond_length` nm along +x from the previous tail. Templates keep their
  reference geometry otherwise.
- A link whose neighbouring angles are nearly straight (`|cos θ| > 0.996`) is
  refused. The angle force is singular there, so the energy is not smooth.
  Templates need a bent geometry at their head and tail.

What is and is not established:

- The tests in `tpt-dynamics-sim/tests/residue_chain.rs` check topology counts
  for a linear chain, link lengths, residue numbering, separate chains, the
  error paths, and energy conservation under NVE (relative drift
  −6.9e-4 over 4000 steps at dt = 0.0005, reduced units). Finite-difference
  checks agree with the analytic forces to about 1e-6 for one, two and five
  residues.
- The residue templates in those tests use the synthetic test parameters. They
  are not a protein, lipid or nucleic-acid force field, and nothing has been
  validated against reference energies or forces.
- Supporting real biomolecules needs a parameter set chosen and supplied by the
  user (none is bundled, per the policy above) and a reference to validate it
  against.
