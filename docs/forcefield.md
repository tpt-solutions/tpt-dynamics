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
