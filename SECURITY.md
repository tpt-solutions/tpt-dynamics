# Security policy

## Supported versions

tpt-dynamics is pre-release. Only the latest commit on `master` receives security
fixes until v0.1.0 is tagged. After v0.1.0, the latest minor release receives fixes.

## Reporting a vulnerability

Please **do not** open a public issue for a security problem.

Report it privately through GitHub's private vulnerability reporting for this
repository (the **Security** tab, then **Report a vulnerability**). Include:

* the crate, version or commit, and operating system;
* a minimal reproduction, or a description of the input that triggers the problem
  (for example a malformed PDB, GRO, CSV, TOML or binary file);
* the impact you expect (crash, memory unsafety, incorrect results, data leak).

We aim to acknowledge a report within 7 days and to agree a disclosure timeline
with the reporter. Fixes are released with a changelog entry that credits the
reporter, unless they ask not to be named.

## Scope

In scope: parsers and readers (PDB, GRO, TOML, CSV, JSON, binary columnar and
trajectory formats, Arrow-IPC-style reader), the `tptdyn` command line, the Python
bindings, and any `unsafe` code (see `docs/audit/unsafe-audit.md` and ADR-006/007).

Out of scope: denial of service from deliberately enormous inputs that the
documented memory bounds do not promise to handle, and issues in third-party
crates that have already been reported upstream (report those upstream and tell us
so we can update the pin).

## Dependency advisories

`cargo deny check advisories` runs in CI. Advisories are reviewed, not silently
ignored; any ignore entry must name the advisory and the reason.
