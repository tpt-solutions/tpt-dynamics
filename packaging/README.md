# Packaging templates

These files are **drafts**. They are not submitted to Homebrew, winget or any other
index. Every `PLACEHOLDER` value must be replaced with the real release version and
SHA-256 digests (taken from the `.sha256` files attached to the GitHub Release) before
a manifest is used. Do not submit them before the first release exists.

- `homebrew/tptdyn.rb`: Homebrew formula template for the macOS (Apple silicon) and
  Linux archives.
- `winget/`: winget manifest template (schema 1.6.0) for the Windows x64 zip.

The install scripts at the repository root (`install.sh`, `install.ps1`) are the
supported install paths until these are published.
