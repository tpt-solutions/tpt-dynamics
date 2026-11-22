# TEMPLATE. Not for submission. Replace every PLACEHOLDER with the release values.
# Version and SHA-256 digests come from the tptdyn-v<version>-<target>.tar.gz.sha256
# files attached to the GitHub Release.

class Tptdyn < Formula
  desc "Command-line tool for tpt-dynamics: molecular dynamics and HEP-style data analysis"
  homepage "https://github.com/tpt-solutions/tpt-dynamics"
  version "PLACEHOLDER_VERSION"
  license any_of: ["MIT", "Apache-2.0"]

  on_macos do
    on_arm do
      url "https://github.com/tpt-solutions/tpt-dynamics/releases/download/v#{version}/tptdyn-v#{version}-aarch64-apple-darwin.tar.gz"
      sha256 "PLACEHOLDER_SHA256_AARCH64_APPLE_DARWIN"
    end
  end

  on_linux do
    on_intel do
      url "https://github.com/tpt-solutions/tpt-dynamics/releases/download/v#{version}/tptdyn-v#{version}-x86_64-unknown-linux-gnu.tar.gz"
      sha256 "PLACEHOLDER_SHA256_X86_64_UNKNOWN_LINUX_GNU"
    end
    on_arm do
      url "https://github.com/tpt-solutions/tpt-dynamics/releases/download/v#{version}/tptdyn-v#{version}-aarch64-unknown-linux-gnu.tar.gz"
      sha256 "PLACEHOLDER_SHA256_AARCH64_UNKNOWN_LINUX_GNU"
    end
  end

  def install
    bin.install "tptdyn"
    doc.install "README.md", "NOTICE", "LICENSE-MIT", "LICENSE-APACHE"
  end

  test do
    assert_match "tptdyn", shell_output("#{bin}/tptdyn --version")
  end
end
