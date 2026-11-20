# Full-tree licence review

- **Date:** 2026-11-20
- **Commit:** `88c407a` plus the uncommitted working tree (the tree was dirty when reviewed)
- **Tools:** cargo 1.100.0, cargo-deny 0.20.2
- **Policy:** [ADR-005](../adr/ADR-005-licence-policy.md), `deny.toml` (allow list: MIT, BSD-2-Clause, BSD-3-Clause, Zlib, ISC, Unicode-3.0, and `Apache-2.0 WITH LLVM-exception` by name)
- **Scope:** every third-party package in `cargo metadata --all-features`: 116 packages. Each is classed as `normal` (reached from a workspace member through normal dependencies), `build` (reached only through build dependencies), or `dev` (dev-dependency only, not shipped). Workspace members are excluded.

Counts by scope: 113 normal, 3 build, 0 dev. The workspace members' own licences are `MIT OR Apache-2.0` and are not reviewed here.

## Result of `cargo deny check`

```text
advisories ok, bans ok, licenses ok, sources ok
```

Warnings (no failure): `unmatched-source` for `allow-git` (no git dependencies); `license-not-encountered` for `ISC` (allowed, not used); `duplicate` for `syn` 2.0.119 and 3.0.10 (bans `multiple-versions = "warn"`; both are `MIT OR Apache-2.0`).

## Changes since the 2026-11-16 review

- The Apache-2.0-only dev crate `normalize-line-endings` is no longer in `Cargo.lock`. The CLI test harness (`tpt-dynamics-cli/tests/cli.rs`) is written locally and no longer depends on `assert_cmd` or `predicates`. The earlier finding 1 is resolved.
- There are no dev-only packages in the current graph, so the earlier gap in `cargo deny` coverage no longer applies.
- `target-lexicon` is now the MIT-licensed build dependency (`0.13.6`). The earlier `Apache-2.0 WITH LLVM-exception` entry is gone.
- `snap` (BSD-3-Clause, SNAPPY) and `flate2` with its backends `miniz_oxide`, `adler2`, `crc32fast` and `simd-adler32` (GZIP) are pulled in by the `parquet` feature of `tpt-dynamics-data` (`parquet = ["dep:snap", "dep:flate2"]`). They are new since the last review. All pass the allow list.

## Headline findings

1. **`windows-sys` 0.61.2 is a `*-sys` crate in the normal tree on Windows targets (decision needed).** It is a raw Windows API binding (`MIT OR Apache-2.0`). It reaches `tpt-dynamics-cli` through `clap` -> `anstream` -> `anstyle-query` / `anstyle-wincon` (terminal colour output), on Windows targets only. The `zero-ffi` CI job runs on Linux, where `cargo tree` does not resolve Windows-only dependencies, so it does not see this crate. Confirmed by `cargo tree -i windows-sys --target x86_64-pc-windows-msvc -e normal`. Owner decision: record a named exception for Windows system bindings in ADR-005 or a new ADR, or remove the path.
2. **Non-SPDX expression `MIT/Apache-2.0`** on `matrixmultiply` 0.3.11 and `rawpointer` 0.2.1. cargo-deny accepts it, but the expression is not valid SPDX. Upstream fix needed; no local action.
3. **Reliance on the MIT branch is universal.** Every `X OR MIT` or `MIT OR Apache-2.0` expression passes only because its MIT branch is used, since bare `Apache-2.0` is not allowed. This is the policy as written in ADR-005, not a defect.
4. **`0BSD OR MIT OR Apache-2.0`** (`adler2`) and **`MIT OR Zlib OR Apache-2.0`** (`miniz_oxide`) pass through their MIT branch. `0BSD` is not on the allow list, so it is not relied on.
5. **`syn` duplicate** (2.0.119 and 3.0.10). Both are `MIT OR Apache-2.0`. Warning only.

## Table: every distinct licence expression

| Licence expression | Count | Allow-list status | Crates (scope when not normal) |
|---|---|---|---|
| `MIT OR Apache-2.0` | 94 | Allowed via the MIT branch | anstream, anstyle, anstyle-parse, anstyle-query, anstyle-wincon, cfg-if, clap, clap_builder, clap_complete, clap_derive, clap_lex, clap_mangen, colorchoice, crc32fast, crossbeam-deque, crossbeam-epoch, crossbeam-utils, either, flate2, half, hashbrown, heck, is_terminal_polyfill, libc, log, ndarray, num-complex, num-integer, num-traits, once_cell, pyo3 and pyo3 crates, quote, proc-macro2, rayon, rayon-core, roff, serde family, syn 2.0.119 and 3.0.10, toml family, windows-link, windows-sys (Windows only, finding 1), and others (full list: `cargo metadata --all-features`) |
| `Apache-2.0 OR MIT` | 7 | Allowed via the MIT branch | autocfg (build), equivalent, indexmap, portable-atomic, portable-atomic-util, rustc-hash, utf8parse |
| `MIT` | 6 | Allowed (plain) | crunchy, libm, simd-adler32, strsim, target-lexicon (build), winnow |
| `BSD-2-Clause OR Apache-2.0 OR MIT` | 2 | Allowed via the BSD-2-Clause or MIT branch | zerocopy, zerocopy-derive |
| `MIT/Apache-2.0` | 2 | Allowed via the MIT branch (non-SPDX `/` separator, finding 2) | matrixmultiply, rawpointer |
| `(MIT OR Apache-2.0) AND Unicode-3.0` | 1 | Allowed (both terms are allowed) | unicode-ident |
| `0BSD OR MIT OR Apache-2.0` | 1 | Allowed via the MIT branch | adler2 |
| `MIT OR Zlib OR Apache-2.0` | 1 | Allowed via the MIT or Zlib branch | miniz_oxide |
| `BSD-2-Clause` | 1 | Allowed (plain) | numpy (`tpt-dynamics-py`, `publish = false`) |
| `BSD-3-Clause` | 1 | Allowed (plain) | snap (parquet feature) |

## Added 2026-11-20: memory-mapped reader (ADR-011)

- `memmap2` 0.9.12, `MIT OR Apache-2.0`, optional feature `mmap` of `tpt-dynamics-data` only. Normal dependencies: none on Windows; `libc` on Unix (already in the tree). Passes the allow list via the MIT branch.

## Dependency scope summary

- `normal`: shipped in the library or CLI, or used by the Python extension module.
- `build`: build-dependency only (`target-lexicon`, `autocfg`, and PyO3 build helpers).
- `dev`: none in the current graph.

## Method

1. `cargo metadata --format-version 1 --all-features` (all features, all targets, host-independent).
2. Licence strings are the `license` field of each manifest. Each is matched against ADR-005.
3. Scope: reachability from the workspace members. `normal` follows normal edges only. `build` is reached through normal and build edges, and is not reached by normal edges alone. `dev` is whatever remains reachable through the resolve graph.
4. `cargo deny check` and `cargo tree -i windows-sys --target x86_64-pc-windows-msvc -e normal` were run for the zero-FFI check.
