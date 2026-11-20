# ADR-011: Zero-copy memory-mapped reader (unsafe exception)

- **Status:** Accepted
- **Date:** 2026-11-20

## Context

`todo.md` (Phase 1, I/O) has an open item: "Zero-copy / memory-mapped readers where
possible". The columnar format (version 1) is read by copying every column into a
`Vec`. For large detector files this doubles memory use and forces a full read before
any value is seen.

Memory mapping is the standard way to avoid that copy. In Rust, the mapping call
itself is `unsafe`, because the mapping is only sound while the file is not modified by
another party. ADR-006 forbids unsafe code by default and requires an ADR for any
exception.

## Decision

1. **Crate and block.** The exception is in `tpt-dynamics-data`, module
   `src/mapped.rs`, in one function: `map_file`, which calls `memmap2::Mmap::map`. It
   carries `#[allow(unsafe_code)]` and a `// SAFETY:` comment.
2. **Lint change.** The workspace lint `unsafe_code` changes from `forbid` to `deny`,
   because `forbid` cannot be overridden locally. The crates that previously relied on
   the workspace `forbid` all have `#![forbid(unsafe_code)]` in their crate root
   (core, compute, sim, cli, embedded-demo), so their policy is unchanged. The data
   crate root has `#![deny(unsafe_code)]`, and only `mapped.rs` may override it.
3. **Invariant (the safety condition).** The mapping is read-only and is only read
   through shared slices. Soundness requires that no process truncates or rewrites the
   file while it is mapped. This is a caller obligation. The reader verifies every
   header and column CRC when it opens the file, but a CRC check cannot make a file that
   is changing under the mapping safe.
4. **Caller policy.** `MappedTable::open` may only be pointed at complete files that are
   not being written. Files written by `write_table` should be closed before they are
   read through this module. The reader does not check for concurrent writers.
5. **Dependency.** `memmap2` 0.9.12 (MIT OR Apache-2.0), behind the optional feature
   `mmap`. It has no new normal dependencies on Windows. On Unix it depends on `libc`,
   which is already in the tree. The crate itself contains the platform unsafe code.
   Its licence and FFI status are recorded in `docs/audit/licence-review.md`.
6. **No format change.** Values are decoded from the mapped bytes with
   `from_le_bytes` as they are read. Typed slices are not reinterpreted in place, so
   the format does not need alignment padding and version 1 is unchanged.

## Why safe code is insufficient

Without `unsafe`, a file can only be read into owned memory. A mapping gives a view of
the file's pages with no copy, and the process reads them lazily. Safe Rust has no
mapping primitive in the standard library, so the mapping call is the minimum unsafe
surface.

## Consequences

- The data crate is no longer `forbid(unsafe_code)`. It is `deny`, with one excepted
  function. Phase 8's unsafe-code audit must review this block (`docs/audit/unsafe-audit.md`).
- The feature `mmap` is off by default. CI does not yet build or test it. Before it
  is relied on in a release, add `--features mmap` to the clippy and test jobs.
- Mapped data cannot outlive the `MappedTable` that owns the mapping. The API returns
  iterators that borrow from it, so the borrow checker enforces the lifetime.
