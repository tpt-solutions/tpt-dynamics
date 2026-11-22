# Release checklist

This checklist is the gate sequence for a tpt-dynamics release (v0.1.0 and later).
Every gate must pass on a clean checkout of the commit being released. Run the
commands from the repository root.

> **Publishing is not approved.** The publish steps in section 5 are marked
> **DO NOT RUN until approved**. Only the dry runs in section 4 may be run before
> approval. A dry run never uploads anything, but it still reads the crates.io index.

## 1. Quality gates (CI and local)

| # | Gate | Command | CI job |
|---|------|---------|--------|
| 1 | Formatting | `cargo fmt --all -- --check` | `fmt` |
| 2 | Clippy, warnings denied | `cargo clippy --workspace --all-targets -- -D warnings` | `clippy` |
| 3 | Tests, all OS and arch entries | `cargo test --workspace` | `test` |
| 4 | Rustdoc, warnings denied | `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps` | `doc` |
| 5 | Dependency gate | `cargo deny check` (see section 2) | `deny` |
| 6 | Zero-FFI: no `*-sys` crate in the normal tree | see section 3 | `zero-ffi` |
| 7 | `no_std` build of the core | `cargo build -p tpt-dynamics-core --no-default-features --target thumbv7em-none-eabihf` | `no-std` |
| 8 | MSRV (1.95) | `cargo +1.95 check --workspace --all-targets` | `msrv` |
| 9 | Embedded demo | `cargo build --manifest-path embedded-demo/Cargo.toml --target thumbv7em-none-eabihf`, then `cargo test --manifest-path embedded-demo/Cargo.toml` | `embedded-demo` (build only) |
| 10 | Benchmark regression (25% limit) | `cargo run --release -p tpt-dynamics-data --example bench_json -- --check docs/benchmarks/baseline.json` | `bench-regression` (PRs) |
| 11 | Python wheels: abi3 build for each target | `maturin build --release --out dist` in `bindings/tpt-dynamics-py` | `wheels` (`build`) |
| 12 | Python tests against the built wheel | `python -m pytest bindings/tpt-dynamics-py/tests -q` | `wheels` (`test`) |
| 13 | Source distribution | `maturin sdist --out dist` in `bindings/tpt-dynamics-py` | `wheels` (`sdist`) |
| 14 | Documentation index builds and matches the committed `docs/site` | `python docs/site/build.py`, then `git diff --exit-code docs/site` | `docs-site` |

Gate 3 runs on these runners: `ubuntu-latest` (x86_64), `ubuntu-24.04-arm` (aarch64),
`windows-latest` (x86_64), `windows-11-arm` (aarch64, experimental: failures do not
block), and `macos-latest` (aarch64). The job name is `test (<os>, <arch>)`.
Branch protection must use the names as they now appear.

Before tagging, confirm that the release commit is green for every CI job, including
the experimental `windows-11-arm` entry (record any failure in the release notes).

## 2. cargo-deny pre-release gate

1. `cargo deny check` must report `advisories ok, bans ok, licenses ok, sources ok`.
   Warnings are reviewed, not ignored. At the time of writing the expected warnings
   are: an unused source allowance (`allow-git`), unused licence allowances (`BSD-3-Clause`,
   `ISC`, `Zlib`), and the duplicate `syn` 2.x/3.x entries (bans `multiple-versions = "warn"`).
2. Known limitation: cargo-deny 0.20.2 with the current `deny.toml` does not list
   dev-dependencies of workspace members (for example `assert_cmd`, `predicates`,
   `normalize-line-endings`). Those crates therefore pass the gate unchecked. The
   full-tree review in section 3 covers them manually until this is resolved.
3. Run `cargo deny list -l crate` and compare the output with
   `docs/audit/licence-review.md`. Any new crate, licence or version change needs a
   licence review entry in that file before the release.

## 3. Full-tree licence review procedure

The review covers every package in `cargo metadata`, including dev-dependencies and
build-dependencies, because `cargo deny` does not currently report dev-only packages.

1. Regenerate the package list and licence expressions:
   `cargo metadata --format-version 1 --all-features > <scratch>/meta.json`.
   Group the third-party packages by their `license` field. Record each package's
   dependency kinds (normal, build, dev) from the resolve graph.
2. For each distinct expression, check it against ADR-005:
   * A plain expression from the allow list (`MIT`, `BSD-2-Clause`, `BSD-3-Clause`,
     `Zlib`, `ISC`, `Unicode-3.0`) passes.
   * An `OR` expression passes if at least one branch is on the allow list. Record which
     branch is relied on. Note any expression that relies on the MIT branch alone.
   * `Apache-2.0` alone is rejected (ADR-005). `Apache-2.0 WITH LLVM-exception` is
     allowed by name only (ADR-005 amendment).
   * Any GPL, LGPL, AGPL or MPL term is rejected.
   * A non-SPDX form (for example `MIT/Apache-2.0`) is flagged for a manifest fix upstream.
3. Check the zero-FFI rules: no crate named `*-sys` in the **normal** tree, and no
   `cc`, `cmake` or `bindgen` in any normal or build tree. Note that `windows-sys` is a
   raw Windows API binding and is present through the Windows target only, so the
   Linux-hosted `zero-ffi` job does not see it. Decide and record that exception in
   ADR-005 or a new ADR before release.
4. Write the results to `docs/audit/licence-review.md` with the date, the commit,
   and the tool versions. Every flagged item needs an owner decision recorded there.

## 4. cargo publish --dry-run (dry run only)

A dry run packages the crate, verifies that the package builds from the packaged
files, and stops before the upload. It never publishes. Run it for each crate in the
publish order in section 5.

1. Commit the release state first. A dry run refuses a dirty tree. Use
   `--allow-dirty` only for a local preview on uncommitted work, and never for the
   release record.
2. Per crate: `cargo publish --dry-run -p <crate>`. Record the final line of output.
3. Expected behaviour: `tpt-dynamics-core` has no internal dependencies, so its dry run
   can pass before anything is published. The other crates depend on internal crates
   that are not yet on crates.io. Their dry runs resolve those dependencies against the
   registry and will fail until the upstream crates are live. Do not work around that
   by editing manifests. Record the failure and re-run after each upstream publish, or
   rely on the `cargo package --list` check (below) for contents only.
4. Content check for every crate: `cargo package --list -p <crate>`. Confirm that
   `README.md`, `LICENSE-MIT`, `LICENSE-APACHE` and `NOTICE` (where the crate ships
   them) are present, and that no scratch files, `target/` or fixtures are included.
5. `tpt-dynamics-py` has `publish = false`. Do not run a dry run for it on crates.io.
   Its distribution is through PyPI wheels (section 5, step 6).

## 5. Publish order

**DO NOT RUN until approved.** Every step below uploads to a public registry and
cannot be undone (a published version can be yanked but never replaced).

1. `tpt-dynamics-core` (no internal dependencies). **DO NOT RUN until approved.**
2. `tpt-dynamics-data` and `tpt-dynamics-compute`. Both depend only on external
   crates, so they can go in either order. **DO NOT RUN until approved.**
3. `tpt-dynamics-sim` (depends on core and compute). **DO NOT RUN until approved.**
4. `tpt-dynamics-cli` (depends on core, data and sim). **DO NOT RUN until approved.**
5. Wait for each crate to appear on the index before publishing the next.
6. PyPI: the `tpt-dynamics` wheels and sdist from the wheels workflow. Upload is a
   separate approval. **DO NOT RUN until approved.**
7. Tag the release (`v0.1.0`), publish the release notes, and update `CHANGELOG.md`.
   **DO NOT RUN until approved.**

## 6. Sign-off record

| Item | Result | Date | By |
|------|--------|------|----|
| Gates 1 to 14 green on release commit | | | |
| `cargo deny check` summary recorded | | | |
| Licence review updated (`docs/audit/licence-review.md`) | | | |
| Dry runs recorded (section 4) | | | |
| Publish approval for steps in section 5 | | | |
