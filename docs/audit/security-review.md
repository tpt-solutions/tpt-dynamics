# Security review: parser attack surface and resource limits

- **Date:** 2026-11-16
- **Scope:** untrusted-input readers in `tpt-dynamics-data` and
  `tpt-dynamics-sim`, plus the trajectory reader if present. Fuzz coverage is
  in `tpt-dynamics-data/tests/fuzz_parsers.rs`; see `docs/fuzzing.md`.
- **Method:** code reading of each reader, a grep for allocations driven by
  file-declared lengths (`Vec::with_capacity`, `vec![0; n]`, `read_exact`,
  `take`), and the randomised tests above.
- **Scope note:** findings in files outside `tpt-dynamics-data/src/{arrow_ipc,jagged}.rs`
  are reported only. They were not modified.

## Attack surface

| Reader | Crate / file | Input | Entry point | Status |
|---|---|---|---|---|
| CSV | data / `csv.rs` | any byte stream | `CsvReader::new` | Line length capped by `max_line_bytes` (default 1 MiB). Fuzzed. |
| Columnar binary | data / `columnar.rs` | byte stream | `read_table` | Length fields checked, reads grow with real data. Fuzzed. |
| Histogram binary | data / `serialise.rs` | `&[u8]` | `from_bytes` | Counts checked against remaining input before allocation. CRC-32 verified first. Fuzzed. |
| Histogram JSON | data / `serialise.rs`, `json.rs` | `&str` | `from_json` | Nesting depth capped (`MAX_DEPTH`). Input size not capped (caller's responsibility). Fuzzed. |
| Arrow-IPC-style | data / `arrow_ipc.rs` (new) | `&[u8]` | `read_all` | Name, column, and row limits. Every length checked against remaining bytes. Fuzzed. |
| PDB | sim / `pdb.rs` | `BufRead` | `read_pdb` | Uses `BufRead::lines()`, so one line is unbounded. **Finding S-1.** |
| GRO | sim / `gro.rs` | `BufRead` | `read_gro` | Same as PDB (**S-1**). Atom count capped at pre-allocation (`min(1 << 20)`). |
| TOML config | sim / `config.rs` | `&str` | `from_toml_str` | Parsed by the `toml` crate. No size cap. **Finding S-3.** |
| Trajectory | data or sim (in progress) | file | not yet on `master` | Not reviewed. Must follow the limits in this document. |
| Checkpoint | none | n/a | Phase 5 | Not implemented. |

## Resource limits: allocations driven by file-declared lengths

| Site | Driven by | Bound | Verdict |
|---|---|---|---|
| `columnar.rs` `vec![0u8; len]` (column name) | u16 length field | 65,535 bytes | Bounded. |
| `columnar.rs` `Vec::with_capacity(ncols.min(1024))` | u32 column count | 1024 | Bounded. The loop still runs `ncols` times, but each iteration consumes input, so it stops at EOF. |
| `columnar.rs` `take(nbytes).read_to_end` | `nrows * width` | Grows with bytes actually read; `checked_mul` used | Bounded by input size. |
| `serialise.rs` `Dec::len(elem_bytes)` | u64 count | Checked `n <= remaining / elem_bytes` before any allocation | Bounded by input size. |
| `arrow_ipc.rs` | u16 ncols, u16 name length, u32 nrows | 1024 columns, 256 bytes, 2^24 rows; `take` checks remaining bytes before every allocation | Bounded. Tested with `u32::MAX` rows. |
| `csv.rs` line buffer | physical line | `max_line_bytes` | Bounded. |
| `sim/gro.rs` `with_capacity(count.min(1 << 20))` | atom count in header | 2^20 | Bounded. |
| `json.rs` recursion | nesting depth | `MAX_DEPTH` | Bounded. |

No unbounded allocation driven by a file length was found in the data crate.

## Findings

| ID | Severity | Location | Description | Recommendation |
|---|---|---|---|---|
| S-1 | Medium (memory DoS) | `tpt-dynamics-sim/src/pdb.rs:49`, `gro.rs:45` (`BufRead::lines`) | A single very long line (no newline) is read fully into memory before parsing. The CSV reader already caps this. | Replace `lines()` with a bounded line reader (for example, the approach in `csv.rs::read_line` with a `max_line_bytes`). Owner: sim crate. |
| S-2 | Low | `tpt-dynamics-sim/src/gro.rs` | The atom count is parsed from text. The pre-allocation is capped, but the loop trusts the count until EOF. | Acceptable. Optionally reject counts above a configured maximum. |
| S-3 | Low | `tpt-dynamics-sim/src/config.rs::from_toml_str` | No input size cap. The `toml` crate's depth handling was not verified in this review. | Cap input length (for example, 1 MiB) before parsing. |
| S-4 | Info | `bindings/tpt-dynamics-py` | `unsafe` is allowed for generated glue, and ADR-007 is referenced but does not exist. | Write ADR-007 before release (see `docs/audit/unsafe-audit.md`). |
| S-5 | Info | `tpt-dynamics-data/src/columnar.rs` | Corrupt files can yield `Truncated` or `ChecksumMismatch` only after reading the column. Memory is bounded by the file, so this is not a vulnerability. | None. |
| S-6 | Info | Trajectory reader | Not on `master` at review time. | Apply the same rules: length checks against remaining bytes, caps on counts and line lengths, and fuzz tests before merge. |

## Fixes made

None in existing files. The new files `arrow_ipc.rs` and `jagged.rs` were written
with the limits above. Findings S-1 to S-3 belong to the sim crate and are
reported only, as the scope required.

## Verification

- `cargo test -p tpt-dynamics-data --test fuzz_parsers` passes with the
  default iteration count, and with `FUZZ_ITERS` set higher (see
  `docs/fuzzing.md`).
- The truncation tests check every prefix of a valid input for the columnar and
  Arrow-IPC readers. The serialise tests do the same for each histogram type.
