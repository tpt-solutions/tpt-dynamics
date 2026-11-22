# Security Policy

## Supported versions

tpt-dynamics is pre-release (`0.1.0-dev`). Security fixes are applied to the
`master` branch. Once v0.1.0 is published, the latest 0.1.x release will be
supported.

## Reporting a vulnerability

Please do **not** open a public issue for a security problem.

Report it privately through GitHub's private vulnerability reporting:
<https://github.com/tpt-solutions/tpt-dynamics/security/advisories/new>

Include:

- the affected crate(s) and version or commit,
- a description of the problem and its impact,
- steps or a minimal input that reproduce it (for parsers, a malformed file
  or byte sequence is usually enough),
- whether you would like to be credited in the advisory.

We aim to acknowledge a report within 5 working days and to agree a fix and
disclosure timeline with the reporter. Reports that involve untrusted input
to the parsers (PDB, GRO, CSV, Parquet, checkpoint, trajectory and network
data sources) are treated as high priority.

## Scope notes

- The crates are written in safe Rust. The one reviewed `unsafe` exception is
  the memory-mapped reader behind the `mmap` feature of `tpt-dynamics-data`
  (see ADR-011 in `docs/adr/`).
- The network data source reads plain HTTP only; TLS is not implemented, so
  do not send sensitive data over it.
