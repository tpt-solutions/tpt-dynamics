# Full-tree licence review

- **Date:** 2026-11-16
- **Commit:** `0bca58a` plus the uncommitted working tree (the tree was dirty when reviewed)
- **Tools:** cargo 1.100.0, cargo-deny 0.20.2
- **Policy:** [ADR-005](../adr/ADR-005-licence-policy.md), `deny.toml` (allow list: MIT, BSD-2-Clause, BSD-3-Clause, Zlib, ISC, Unicode-3.0, and `Apache-2.0 WITH LLVM-exception` by name)
- **Scope:** every third-party package in `cargo metadata --all-features`: 125 packages. Each is classed as `normal` (reached from a workspace member through normal dependencies), `build` (build-dependency only), or `dev` (dev-dependency only, not shipped). Workspace members are excluded.

Counts by scope: 107 normal, 3 build, 15 dev. The workspace members' own
licences are `MIT OR Apache-2.0` and are not reviewed here.

## Result of `cargo deny check`

```text
advisories ok, bans ok, licenses ok, sources ok
```

Warnings (no failure): `unmatched-source` for `allow-git` (no git dependencies);
`license-not-encountered` for `BSD-3-Clause`, `ISC` and `Zlib` (allowed, not used);
`duplicate` for `syn` 2.0.119 and 3.0.10 (bans `multiple-versions = "warn"`; both are `MIT OR Apache-2.0`).

## Headline findings

1. **Apache-2.0-only crate in the tree (dev-only): `normalize-line-endings` 0.3.0.** It is pulled by `predicates` 3.1.4 -> `assert_cmd` 2.2.2 -> `tpt-dynamics-cli` dev-dependencies. ADR-005 rejects bare `Apache-2.0`. `cargo deny check` passes because cargo-deny 0.20.2 with the current config does not list dev-dependencies at all (`cargo deny list` shows 113 packages; `cargo metadata` shows 125). Checked with a scratch config that sets `[graph] exclude-dev = false`: the result did not change. The cause is not confirmed. **Decision needed:** replace or drop the dev-dependency (for example, hand-written assertions), or record a dev-only exception in ADR-005 and make the gate cover dev-dependencies.
2. **Dev-only crates are invisible to the gate.** All 15 `dev` packages are unchecked by `cargo deny`: `assert_cmd`, `predicates`, `predicates-core`, `predicates-tree`, `bstr`, `difflib`, `float-cmp`, `regex`, `regex-automata`, `regex-syntax`, `aho-corasick`, `memchr`, `termtree`, `wait-timeout`, `normalize-line-endings`. Their licences were checked by hand in this review; the full list is in the table below.
3. **Non-SPDX expression `MIT/Apache-2.0`** on `matrixmultiply` 0.3.11, `rawpointer` 0.2.1 and (dev) `wait-timeout` 0.2.1. cargo-deny accepted it, but the expression is not valid SPDX. Upstream fix needed; no local action.
4. **`Unlicense OR MIT`** on `aho-corasick` and `memchr` (both dev-only). `Unlicense` is not on the allow list. It is accepted only through the MIT branch. Acceptable under ADR-005 wording, but it is a dual-licence reliance on the MIT branch and should be recorded.
5. **Reliance on the MIT branch is universal.** Every `X OR MIT` or `MIT OR Apache-2.0` expression passes only because its MIT branch is used, since bare `Apache-2.0` is not allowed. This is the policy as written in ADR-005, not a defect. The non-trivial cases are listed in the table: `Apache-2.0 OR MIT` (7), `MIT/Apache-2.0` (3), `Unlicense OR MIT` (2), and `BSD-2-Clause OR Apache-2.0 OR MIT` (zerocopy, which uses the BSD-2-Clause or MIT branch).
6. **`windows-sys` 0.61.2 is a `*-sys` crate in the normal tree on Windows targets.** It is a raw Windows API binding (`MIT OR Apache-2.0`). It reaches `tpt-dynamics-cli` through `clap` -> `anstream` -> `anstyle-query` / `anstyle-wincon` (terminal colour output), on Windows targets only. The `zero-ffi` CI job runs on Linux, where `cargo tree` does not resolve Windows-only dependencies, so it does not see this crate. Decision needed: record a named exception for Windows system bindings, or remove the path.
7. **`target-lexicon`** (upstream `Apache-2.0 WITH LLVM-exception`, build scope via `pyo3-build-config`) is replaced by the MIT-licensed `target-lexicon-mit` through `[patch.crates-io]`. The exception has been removed from `deny.toml` (ADR-005 amendment, later 2026-11-16). Resolved.
8. **`numpy` 0.30.0** is `BSD-2-Clause` (allowed, normal scope of `tpt-dynamics-py`, which is `publish = false`).
9. **Cargo-deny allow entries not used:** `BSD-3-Clause`, `ISC`, `Zlib`. Harmless. Keep for planned crates or remove when next tidying `deny.toml`.

## Table: every distinct licence expression

| Licence expression | Count | Allow-list status | Crates (scope when not normal) |
|---|---|---|---|
| `MIT OR Apache-2.0` | 100 | Allowed via the MIT branch | anstream@1.0.0, anstyle@1.0.14, anstyle-parse@1.0.0, anstyle-query@1.1.5, anstyle-wincon@3.0.11, assert_cmd@2.2.2 (dev), bstr@1.13.1 (dev), cfg-if@1.0.5, clap@4.6.7, clap_builder@4.6.7, clap_complete@4.6.11, clap_derive@4.6.7, clap_lex@1.1.1, clap_mangen@0.3.3, colorchoice@1.0.5, crossbeam-deque@0.8.8, crossbeam-epoch@0.9.21, crossbeam-utils@0.8.23, either@1.19.0, half@2.7.1, hashbrown@0.17.1, heck@0.5.0, is_terminal_polyfill@1.70.2, libc@0.2.191, log@0.4.35, ndarray@0.17.2, num-complex@0.4.6, num-integer@0.1.47, num-traits@0.2.19, once_cell@1.21.4, once_cell_polyfill@1.70.2, predicates@3.1.4 (dev), predicates-core@1.0.10 (dev), predicates-tree@1.0.13 (dev), proc-macro2@1.0.107, pyo3@0.30.1, pyo3-build-config@0.30.1 (build), pyo3-ffi@0.30.1, pyo3-macros@0.30.1, pyo3-macros-backend@0.30.1, quote@1.0.47, rayon@1.12.0, rayon-core@1.13.0, regex@1.13.1 (dev), regex-automata@0.4.18 (dev), regex-syntax@0.8.11 (dev), roff@1.1.1, serde@1.0.229, serde_core@1.0.229, serde_derive@1.0.229, serde_spanned@1.1.3, syn@2.0.119, syn@3.0.10, toml@1.2.1+spec-1.1.0, toml_datetime@1.1.3+spec-1.1.0, toml_parser@1.1.5+spec-1.1.0, toml_writer@1.2.0+spec-1.1.0, tpt-simd@0.1.1, tpt-simd-aligned@0.1.1, tpt-simd-audio@0.1.1, tpt-simd-bits@0.1.1, tpt-simd-blas@0.1.1, tpt-simd-blend@0.1.1, tpt-simd-butterfly@0.1.1, tpt-simd-compare@0.1.1, tpt-simd-complex@0.1.1, tpt-simd-convolve@0.1.1, tpt-simd-core@0.1.1, tpt-simd-dct@0.1.1, tpt-simd-distance@0.1.1, tpt-simd-dot@0.1.1, tpt-simd-fixed@0.1.1, tpt-simd-gather@0.1.1, tpt-simd-geom@0.1.1, tpt-simd-half@0.1.1, tpt-simd-hash@0.1.1, tpt-simd-horizontal@0.1.1, tpt-simd-interpolate@0.1.1, tpt-simd-math@0.1.1, tpt-simd-matrix@0.1.1, tpt-simd-mul@0.1.1, tpt-simd-ntt@0.1.1, tpt-simd-permute@0.1.1, tpt-simd-pixel@0.1.1, tpt-simd-reduce@0.1.1, tpt-simd-rng@0.1.1, tpt-simd-rounding@0.1.1, tpt-simd-saturate@0.1.1, tpt-simd-scatter@0.1.1, tpt-simd-search@0.1.1, tpt-simd-select@0.1.1, tpt-simd-shift@0.1.1, tpt-simd-sort@0.1.1, tpt-simd-sparse@0.1.1, tpt-simd-text@0.1.1, tpt-simd-transpose@0.1.1, tpt-simd-vector@0.1.1, tpt-simd-window@0.1.1, windows-link@0.2.1, windows-sys@0.61.2 |
| `Apache-2.0 OR MIT` | 7 | Allowed via the MIT branch | autocfg@1.5.1 (build), equivalent@1.0.3, indexmap@2.14.2, portable-atomic@1.15.0, portable-atomic-util@0.2.8, rustc-hash@2.1.3, utf8parse@0.2.2 |
| `MIT` | 7 | Allowed (plain) | crunchy@0.2.4, difflib@0.4.0 (dev), float-cmp@0.10.0 (dev), libm@0.2.16, strsim@0.11.1, termtree@0.5.1 (dev), winnow@1.0.5 |
| `MIT/Apache-2.0` | 3 | Allowed via the MIT branch (non-SPDX `/` separator, see finding 3) | matrixmultiply@0.3.11, rawpointer@0.2.1, wait-timeout@0.2.1 (dev) |
| `BSD-2-Clause OR Apache-2.0 OR MIT` | 2 | Allowed via the BSD-2-Clause or MIT branch | zerocopy@0.8.63, zerocopy-derive@0.8.63 |
| `Unlicense OR MIT` | 2 | Allowed via the MIT branch (`Unlicense` is not on the allow list, see finding 4) | aho-corasick@1.1.5 (dev), memchr@2.8.3 (dev) |
| `(MIT OR Apache-2.0) AND Unicode-3.0` | 1 | Allowed (both terms are allowed) | unicode-ident@1.0.27 |
| `Apache-2.0` | 1 | NOT ALLOWED under ADR-005 (bare Apache-2.0). Dev-only, so cargo deny does not report it (finding 1) | normalize-line-endings@0.3.0 (dev) |
| `BSD-2-Clause` | 1 | Allowed (plain) | numpy@0.30.0 |

## Dependency scope summary

- `normal`: shipped in the library or CLI, or used by the Python extension module.
- `build`: build-dependency only (`target-lexicon`, `pyo3-build-config`, related PyO3 build helpers).
- `dev`: `dev-dependencies` only. Not shipped, but not checked by `cargo deny` (finding 2).

## Method

1. `cargo metadata --format-version 1 --all-features` (all features, all targets, host-independent).
2. Licence strings are the `license` field of each manifest. Each is matched against ADR-005.
3. Scope: breadth-first reachability from the workspace members over `normal`, `build`, and `dev` edges in the resolve graph.
4. `cargo deny list -l crate` was compared to the metadata list to find the dev-dependency gap.

Regenerate by repeating the method above; the policy text is in ADR-005.
