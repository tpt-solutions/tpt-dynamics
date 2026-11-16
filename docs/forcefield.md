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

## Parameter-set recommendation (2026-11-16)

Recommended target: the AMBER family, **ff14SB** for proteins, **OL15** for DNA,
**OL3** for RNA and **lipid17** for lipids. It covers all three kinds of molecule
with one family, and OpenMM ships the same parameters as XML, so an independent
reference is available locally.

Licence status, as far as it could be checked from the sources found:

- A 2023 J. Chem. Inf. Model. paper on AmberTools states that the Amber force
  fields are in the public domain and are distributed with AmberTools. This is
  a secondary source about the family as a class, not a per-file licence. The
  search results reported it from AmberTools' own licence pages on ambermd.org.
  The header of each file must be checked before any file is bundled.
- No licence was verified for CHARMM36 (the forum thread and the
  Mackerell-lab notes describe informal permission, not a licence), GROMOS 54A7
  or OPLS-AA. These are not recommended as the first target.

Policy stays as before until the file-level check is done: no parameter data is
bundled. Users obtain the files from the source and point the converter at them.

## OpenMM cross-check

`validation/openmm_chain_reference.py` builds residue chains with its own
geometry and topology code and evaluates them with OpenMM 8.6 (Reference
platform). `tests/openmm_reference.rs` rebuilds the same chains with
`assemble_chains` and checks:

- positions agree to 1e-12 nm;
- the total potential agrees to 1e-9 relative, for a five-residue chain of four
  carbon residues, a three-residue chain of six-carbon units, and a mixed chain.

Scope: this uses the synthetic test parameters, so it checks the implementation
(bonds, angles, proper dihedrals, exclusions, 1-4 scaling, Lennard-Jones mixing,
cutoff without shift). It does not check any parameter set. The test was checked
for sensitivity: a 1 kJ/mol change in the reference makes it fail.

Forces are also compared (1e-9 relative to the largest force), so the gradients
agree as well as the energies.

Not yet covered: electrostatics, impropers, Urey-Bradley terms, and periodic images.

### Licence check, second pass (2026-11-16)

- The AMBER parameter files are distributed with AmberTools. The AmberTools
  page says its components are mostly GPL, with some public-domain and other
  open-source parts. It does not name a licence for the parameter files.
- The AmberModels page gives no licence or redistribution terms for them.
- The `openmmforcefields` repository is MIT licensed. That licence covers its
  code, not the AMBER parameter data it converts.

So the public-domain statement for the AMBER files is still secondary. It must be
confirmed from the header or licence notice of the files the user obtains from
AmberTools before anything is bundled.

## Amber conversion, version 1 (2026-11-16)

`validation/amber_convert.py` converts the user's own Amber XML, and
`validation/amber_validate.py` checks the result against OpenMM. For ACE-ALA-NME
from ff14SB, bonded and Lennard-Jones terms agree with OpenMM exactly. Impropers
and electrostatics are not yet covered (see ADR-010).

The check also found that `Electrostatics::None` had applied Coulomb to 1-4 pairs
only. That is fixed. The fix is in the engine, and a regression test covers it.

Electrostatics are now checked too. With PME on both sides, using the engine's
Ewald parameters, the potential and forces agree with OpenMM to about 1e-10
(electrostatics alone) and 1e-13 (all terms). Impropers are the remaining gap.

Impropers are checked too (periodic form, as in Amber). For ACE-ALA-NME,
ACE-ALA-ALA-NME and ACE-GLY-SER-NME from ff14SB, all terms agree with OpenMM to
about 1e-14 in potential and forces. Urey-Bradley terms are not yet converted.
