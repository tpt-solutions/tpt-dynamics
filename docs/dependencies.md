# Dependency audit

Every dependency is recorded here with its licence and the decision taken. The
rules are in [ADR-005](adr/ADR-005-licence-policy.md). `cargo deny check` enforces
them in CI.

Checked on 2026-11-16 from the crate tarballs on crates.io.

## Accepted

| Crate | Version checked | Licence | Notes |
|-------|-----------------|---------|-------|
| `libm` | 0.2 | MIT OR Apache-2.0 | Used by `tpt-dynamics-core` for `sqrt` in `no_std` builds. |

## Candidates (Phase 0 audit)

| Crate | Version checked | Licence | Decision |
|-------|-----------------|---------|----------|
| `tpt-math` | 0.1.0 | MIT OR Apache-2.0 | Candidate for FFT. Optional deps are many; enable only the feature needed. MSRV 1.95. Not yet added. |
| `tpt-simd` | 0.1.1 | MIT OR Apache-2.0 | Candidate for SIMD kernels (Phase 3). Pulls 39 `tpt-simd-*` sub-crates, all MIT OR Apache-2.0 (checked). **MSRV 1.95**. Adopted in `tpt-dynamics-compute`; the workspace MSRV is now 1.95. |
| `tpt-strata` | 0.1.0 | MIT OR Apache-2.0 | Candidate for the data engine. **No dependencies.** Adopt (ADR-004). Not yet added. |
| `tpt-strata-parquet` | 0.1.0 | MIT OR Apache-2.0 | **Rejected.** Depends on `arrow` 60 (Apache-2.0 only) and `parquet` 60. |
| `snap` | 1.x | BSD-3-Clause | Adopted for the `parquet` feature (SNAPPY pages). Pure Rust. |
| `flate2` | 1.x | MIT OR Apache-2.0 | Adopted for the `parquet` feature (GZIP pages), `rust_backend` only. |
| `glam` | 0.34.3 | MIT OR Apache-2.0 | Deferred. The local `Vec3` is used for now. Revisit if SIMD layout needs it. |
| `clap` | 4.6.7 | MIT OR Apache-2.0 | Planned for the CLI (Phase 4). |
| `pyo3` | 0.30.1 | MIT OR Apache-2.0 | Planned for Python bindings (Phase 4). |
| `rayon` | 1.12.0 | MIT OR Apache-2.0 | Planned for shared-memory parallelism. |
| `maturin` | 1.16.0 | MIT OR Apache-2.0 | Build tool for Python wheels (not a Rust dependency). |
| `arrow` | 60.1.0 | Apache-2.0 | **Rejected** (Apache-2.0 only; banned in `deny.toml`). |

## Notes

- The dual-licence check is on the crate's `license` field from its manifest.
  Each dependency's own transitive tree must also be checked by `cargo deny`.
- `tpt-simd` sets the workspace MSRV (1.95). `tpt-math` (1.85) is not yet added;
  it needs no reconciliation unless adopted.
- **FFT decision (2026-11-16):** the 3D FFT for PME is the in-house radix-2 implementation in `tpt-dynamics-compute/src/fft.rs`. `tpt-math` is not adopted for FFT: the in-house code has no dependencies, is covered by tests, and avoids a new transitive tree that would need its own licence review. Revisit only if non-power-of-two lengths are needed.

## Adopted during Phases 1–8 (checked 2026-11-16)

| Crate | Licence | Used by | Notes |
|-------|---------|---------|-------|
| `tpt-simd` 0.1.1 (+ 39 sub-crates) | MIT OR Apache-2.0 | compute | Sets workspace MSRV 1.95. |
| `rayon` 1.x | MIT OR Apache-2.0 | data, sim, compute | Parallel chunks and neighbour lists. |
| `libm` 0.2 | MIT OR Apache-2.0 | core | `no_std` maths. |
| `toml` 1.2.1, `serde` 1.x | MIT OR Apache-2.0 | sim, cli | Config files. |
| `clap`, `clap_complete`, `clap_mangen` | MIT OR Apache-2.0 | cli | Subcommands, completions, man pages. |
| `log` | MIT OR Apache-2.0 | cli | Logger is a small in-crate stderr implementation. |
| `assert_cmd`, `predicates` | MIT OR Apache-2.0 | cli (dev) | Integration tests. |
| `pyo3` 0.30.1, `pyo3-ffi`, `pyo3-macros`, `pyo3-build-config` | MIT OR Apache-2.0 | py | Requires ADR-007 (unsafe exception, generated code only). |
| `numpy` 0.30.0 | BSD-2-Clause | py | Allowed by ADR-005. |
| `ndarray` 0.17.2 and its tree (`matrixmultiply`, `num-*`, `rawpointer`, `autocfg`, …) | MIT OR Apache-2.0 | py (via numpy) | |
| `target-lexicon` 0.13.6 | MIT (patched to `target-lexicon-mit`) | py (build dep of pyo3-build-config) | Resolved: upstream Apache-2.0 WITH LLVM-exception crate replaced via `[patch.crates-io]` (ADR-005 amendment, later 2026-11-16). |

Python side (not a Cargo dependency): maturin, pytest, numpy (development only).
