<!--
Pull requests are not accepted at this time (see CONTRIBUTING.md). Please open an
issue instead. This template applies once external pull requests are opened.
-->

## Summary

<!-- What changes, and why. -->

## Checklist

- [ ] `cargo fmt --all -- --check` passes
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` passes
- [ ] `cargo test --workspace` passes
- [ ] `cargo deny check` passes (no new licence, ban or source issues)
- [ ] No `*-sys`, `bindgen`, `cc` or `cmake` dependency added (zero-FFI policy)
- [ ] New or changed dependencies are MIT, BSD-2/3, Zlib, ISC, Unicode or `MIT OR Apache-2.0`
- [ ] `CHANGELOG.md` updated under `[Unreleased]`
- [ ] `todo.md` updated if a checklist item changed state
