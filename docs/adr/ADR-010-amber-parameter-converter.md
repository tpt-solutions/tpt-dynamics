# ADR-010: Converting user-supplied Amber parameters (proposed)

- **Status:** Accepted. Bonds, angles, proper and improper torsions, Lennard-Jones and electrostatics (PME) are implemented and checked against OpenMM for ACE-ALA-NME, ACE-ALA-ALA-NME and ACE-GLY-SER-NME (see "Results" below). Urey-Bradley terms are not converted.
- **Date:** 2026-11-16
- **Related:** ADR-005 (licence policy), ADR-009, `docs/forcefield.md`

## Context

The project reads its own parameter format (`ffparams` TOML). No published
parameter set is bundled. The recommended target for proteins, nucleic acids and
lipids is the AMBER family (ff14SB, OL15, OL3, lipid17). Its licence terms are not
yet confirmed from the primary files (see `docs/forcefield.md`).

The OpenMM cross-check already validates the engine's implementation with
synthetic parameters, down to forces (1e-9 relative). What is missing is a way to
run real AMBER parameters through the engine without shipping them.

## Decision (recommended, not yet approved)

1. **No parameter data in the repository.** The tool converts files the user
   already has, from their own AmberTools or OpenMM install, into `ffparams` TOML
   on their machine. The output stays with the user.
2. **The converter is ours**, Python with OpenMM as an optional dependency. It is
   MIT OR Apache-2.0 like the rest of the code.
3. **OpenMM resolves the parameters at conversion time.** Bond, angle, torsion and
   nonbonded terms are read from the OpenMM `System` built from the user's XML.
   The converter does not reimplement Amber's wildcard and ordering rules. Each
   torsion term that appears in the system is written as an explicit entry, so
   the engine's matching has nothing to resolve. The matching problem is therefore
   removed from the converter, and an unexpected mismatch shows up in validation.
4. **Validation**: a capped dipeptide (ACE-ALA-NME) and a short peptide, built
   with OpenMM's own topology code. Compare positions and energy components
   against OpenMM with the converted parameters.
5. **Electrostatics** are validated separately, after the converter works. The
   engine's PME must be matched to OpenMM's grid and tolerance settings first.

## Consequences

- Repository stays free of third-party parameter data, whatever the licence turns
  out to be.
- Users need OpenMM installed to convert. This is acceptable, and it is a stated
  requirement.
- Explicit per-torsion entries make the parameter files larger. That is acceptable
  for a validation and research workflow.
- Impropers, Urey-Bradley and PME are not covered by the first version.

## Open questions for the owner

- Confirm the licence check. Until the AMBER file headers are read, the converter
  is written but never used on real files in the repo's own tests.
- Confirm OpenMM as the conversion-time dependency (it is installed locally for
  development).
- Confirm that impropers and Urey-Bradley terms come in the second version.

## Version 1 result (2026-11-16)

- Tool: `tpt-dynamics-sim/validation/amber_convert.py`. Check: `amber_validate.py`.
- Test system: ACE-ALA-NME (22 atoms) from the user's `protein.ff14SB.xml`.
  Converted: 11 bond types, 21 angle types, 26 dihedral types. Improper torsions skipped: 4.
- Electrostatics off on both sides. Potential identical. Forces agree to 1e-15 relative.
- The geometry is an embedding, not a real structure, so the energies are large. The
  check tests the implementation, not the chemistry.

Pending:

1. Impropers. Amber's improper rules are wildcard-based and need OpenMM's
   resolution, which the converter does not reimplement yet.
2. Electrostatics. The engine has no unscreened Coulomb mode. The check needs
   either a PME comparison with matched Ewald parameters, or a plain-Coulomb
   mode in the engine.
3. Licence. The Amber files' terms are still unconfirmed from their own headers.

Finding during the check: `Electrostatics::None` applied Coulomb to 1-4 pairs but
not to other pairs. This is fixed, and the regression test is in place. Earlier
simulations with `None` and non-zero charges were affected.

## Electrostatics check (2026-11-16)

`amber_validate.py --pme 1e-5 5` runs the ACE-ALA-NME dipeptide with PME. The
engine reports its Ewald alpha (0.780853 1/nm) and grid (128^3), and OpenMM is
given the same values in a 9 nm periodic box. The B-spline order is 5 on both sides.

- Potential, all terms: relative difference 9e-14. Forces: 1.4e-14.
- Electrostatics alone (`--no-lj`): potential relative difference 9e-11, forces 3e-10.

This covers the 1-4 Coulomb scaling, the exclusion corrections and the PME
reciprocal sum. Impropers are still the only term not checked.

## Results (all terms)

Impropers are now converted. The engine accepts `[[periodicimpropertype]]`, with
the same `k(1 + cos(n phi - phase))` form Amber uses. An improper matched by a
periodic type uses it. Otherwise it falls back to the harmonic `impropertype`.
Each residue sequence is written as one molecule, so impropers can span residues.

Checked against OpenMM with the same positions, all terms on, PME with the engine's
own Ewald parameters:

| Sequence | Atoms | Potential relative | Forces relative |
|---|---|---|---|
| ACE-ALA-NME | 22 | 8.6e-14 | 1.4e-15 |
| ACE-ALA-ALA-NME | 32 | 6.9e-15 | 1.3e-15 |
| ACE-GLY-SER-NME | 30 | 1.2e-14 | 1.1e-15 |

Electrostatics alone (Lennard-Jones off) for ACE-ALA-NME: potential 9e-11, forces 5e-11.

Not yet covered: Urey-Bradley terms (Amber uses them for a few angle types), and
residue-template reuse for long chains. The converter currently writes one molecule
per sequence, which is fine for validation and for small systems.
