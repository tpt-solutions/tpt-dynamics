#!/usr/bin/env sh
# Install the prebuilt tptdyn command-line tool from a GitHub Release.
#
#   curl -fsSL https://raw.githubusercontent.com/tpt-solutions/tpt-dynamics/master/install.sh | sh
#
# Environment variables:
#   TPTDYN_VERSION  release version without the leading v (default: the version below)
#   PREFIX          install prefix; the binary goes to $PREFIX/bin (default: $HOME/.local)
#
# The script downloads the archive for this OS and CPU, checks its SHA-256 against the
# .sha256 file published with the release, and installs the binary. It works only
# after the first GitHub Release exists. Until then it fails with a 404.
#
# Supported: Linux (x86_64, aarch64) and macOS (Apple silicon). Windows: use install.ps1.

set -eu

REPO="tpt-solutions/tpt-dynamics"
# Keep this equal to workspace.package.version in Cargo.toml. Update it with the version.
VERSION="${TPTDYN_VERSION:-0.1.0-dev}"
PREFIX="${PREFIX:-$HOME/.local}"

say() { printf '%s\n' "$*"; }
die() { printf 'install.sh: error: %s\n' "$*" >&2; exit 1; }

need() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

detect_target() {
    os=$(uname -s)
    arch=$(uname -m)
    case "$os" in
        Linux)
            case "$arch" in
                x86_64|amd64) echo "x86_64-unknown-linux-gnu" ;;
                aarch64|arm64) echo "aarch64-unknown-linux-gnu" ;;
                *) die "no prebuilt binary for Linux/$arch" ;;
            esac
            ;;
        Darwin)
            case "$arch" in
                arm64|aarch64) echo "aarch64-apple-darwin" ;;
                *) die "no prebuilt binary for macOS/$arch (Intel Macs are not built; build from source)" ;;
            esac
            ;;
        *)
            die "unsupported operating system: $os (on Windows use install.ps1)"
            ;;
    esac
}

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        die "need sha256sum or shasum to verify the download"
    fi
}

need curl
need tar
need uname
need awk

TARGET=$(detect_target)
NAME="tptdyn-v${VERSION}-${TARGET}"
BASE="https://github.com/${REPO}/releases/download/v${VERSION}"
ARCHIVE="${NAME}.tar.gz"

TMP=$(mktemp -d "${TMPDIR:-/tmp}/tptdyn-install.XXXXXX")
trap 'rm -rf "$TMP"' EXIT INT TERM

say "Downloading ${ARCHIVE} (v${VERSION}, ${TARGET})"
curl -fsSL -o "$TMP/$ARCHIVE" "$BASE/$ARCHIVE" \
    || die "download failed: $BASE/$ARCHIVE (is the release published?)"
curl -fsSL -o "$TMP/$ARCHIVE.sha256" "$BASE/$ARCHIVE.sha256" \
    || die "download failed: $BASE/$ARCHIVE.sha256"

expected=$(awk '{print $1}' "$TMP/$ARCHIVE.sha256")
actual=$(sha256_of "$TMP/$ARCHIVE")
if [ "$expected" != "$actual" ]; then
    die "checksum mismatch for $ARCHIVE (expected $expected, got $actual); nothing was installed"
fi
say "SHA-256 verified"

tar -C "$TMP" -xzf "$TMP/$ARCHIVE"
BIN="$TMP/$NAME/tptdyn"
[ -f "$BIN" ] || die "archive does not contain $NAME/tptdyn"

DEST="$PREFIX/bin"
mkdir -p "$DEST"
install -m 0755 "$BIN" "$DEST/tptdyn" 2>/dev/null || {
    cp "$BIN" "$DEST/tptdyn"
    chmod 0755 "$DEST/tptdyn"
}

say "Installed $DEST/tptdyn"
case ":$PATH:" in
    *":$DEST:"*) ;;
    *) say "Note: $DEST is not on PATH. Add it to your shell profile." ;;
esac
"$DEST/tptdyn" --version || true
