# ADR-010: Converting user-supplied Amber parameters (proposed)

- **Status:** Proposed. No code yet.
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
