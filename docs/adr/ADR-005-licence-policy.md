# ADR-005: Licence policy (no Apache-2.0-only dependencies)

- **Status:** Accepted
- **Date:** 2026-11-16

## Context

The spec's `deny.toml` allows plain `Apache-2.0`. That would admit Apache-2.0-only
crates such as `arrow` and `parquet`, which the project's MIT-first policy bans.
The spec also uses the `unlicensed` and `copyleft` keys, which cargo-deny has
removed.

## Decision

- Allow `MIT`, `BSD-2-Clause`, `BSD-3-Clause`, `Zlib`, `ISC` and `Unicode-3.0`.
- Do **not** allow bare `Apache-2.0`. A dual `MIT OR Apache-2.0` crate passes
  through its MIT branch. An Apache-2.0-only crate is rejected.
- Reject the copyleft family by not allowing it (no explicit list is required).
- Ban by name the crates listed in `deny.toml`, including `arrow` and `parquet`.

This is the deviation from the spec. It was tested: a scratch copy that depends on
`arrow` fails `cargo deny check` with `arrow` reported as banned and the licence
rejected.

## Consequences

- `tpt-strata-parquet` cannot be used (see ADR-004).
- The `*-sys` ban cannot be expressed as a glob in cargo-deny. A CI job
  (`zero-ffi`) fails if any `*-sys` crate is in the normal dependency tree.

## Amendment (2026-11-16): `Apache-2.0 WITH LLVM-exception`

`target-lexicon` (a build dependency of PyO3's build configuration) is licensed
`Apache-2.0 WITH LLVM-exception`. The LLVM exception permits linking and
redistribution without the Apache-2.0 restrictions that matter here, so this exact
SPDX expression is added to the allow list by name. Bare `Apache-2.0` remains not
allowed, and any other `WITH` exception needs a new amendment.

## Amendment (2026-11-16, later): `target-lexicon` replaced; exception removed

The 2026-11-16 exception above is withdrawn. `target-lexicon` is now patched to
`target-lexicon-mit`, an MIT-licensed, clean-room crate with the same name and API
subset that PyO3's build configuration uses. Cargo's `[patch.crates-io]` in the
workspace `Cargo.toml` selects it. `Apache-2.0 WITH LLVM-exception` has been removed
from the `deny.toml` allow list, so the allow list is MIT/BSD/Zlib/ISC/Unicode-3.0
only again. The patch is a local path until the crate is pushed; it then moves to a
git `rev`.
