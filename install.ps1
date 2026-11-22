# Install the prebuilt tptdyn command-line tool from a GitHub Release (Windows x64).
#
#   irm https://raw.githubusercontent.com/tpt-solutions/tpt-dynamics/master/install.ps1 | iex
#
# Environment variables:
#   TPTDYN_VERSION  release version without the leading v (default: the version below)
#   PREFIX          install prefix; the binary goes to $PREFIX\bin (default: $HOME\.local)
#
# The script downloads the zip archive, checks its SHA-256 against the .sha256 file
# published with the release, and installs tptdyn.exe. It works only after the first
# GitHub Release exists. Until then it fails with a 404.

$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'

$Repo = 'tpt-solutions/tpt-dynamics'
# Keep this equal to workspace.package.version in Cargo.toml. Update it with the version.
$Version = if ($env:TPTDYN_VERSION) { $env:TPTDYN_VERSION } else { '0.1.0-dev' }
$Prefix = if ($env:PREFIX) { $env:PREFIX } else { Join-Path $HOME '.local' }

if (-not [Environment]::Is64BitOperatingSystem) {
    throw 'install.ps1: only 64-bit Windows is supported.'
}
if ($env:PROCESSOR_ARCHITECTURE -eq 'ARM64') {
    throw 'install.ps1: no prebuilt Windows ARM64 binary is published. Build from source.'
}

$Target = 'x86_64-pc-windows-msvc'
$Name = "tptdyn-v$Version-$Target"
$Base = "https://github.com/$Repo/releases/download/v$Version"
$Archive = "$Name.zip"

$Tmp = Join-Path ([IO.Path]::GetTempPath()) ("tptdyn-install-" + [Guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $Tmp | Out-Null
try {
    Write-Host "Downloading $Archive (v$Version, $Target)"
    $ArchivePath = Join-Path $Tmp $Archive
    $ShaPath = "$ArchivePath.sha256"
    Invoke-WebRequest -Uri "$Base/$Archive" -OutFile $ArchivePath -UseBasicParsing
    Invoke-WebRequest -Uri "$Base/$Archive.sha256" -OutFile $ShaPath -UseBasicParsing

    $Expected = ((Get-Content -Path $ShaPath -TotalCount 1) -split '\s+')[0].ToLower()
    $Actual = (Get-FileHash -Algorithm SHA256 -Path $ArchivePath).Hash.ToLower()
    if ($Expected -ne $Actual) {
        throw "checksum mismatch for $Archive (expected $Expected, got $Actual); nothing was installed"
    }
    Write-Host 'SHA-256 verified'

    $Extract = Join-Path $Tmp 'extract'
    Expand-Archive -Path $ArchivePath -DestinationPath $Extract
    $Exe = Get-ChildItem -Path $Extract -Recurse -Filter 'tptdyn.exe' | Select-Object -First 1
    if (-not $Exe) { throw "archive does not contain tptdyn.exe" }

    $Dest = Join-Path $Prefix 'bin'
    New-Item -ItemType Directory -Force -Path $Dest | Out-Null
    Copy-Item -Path $Exe.FullName -Destination (Join-Path $Dest 'tptdyn.exe') -Force

    Write-Host "Installed $(Join-Path $Dest 'tptdyn.exe')"
    if (-not (($env:Path -split ';') -contains $Dest)) {
        Write-Host "Note: $Dest is not on PATH. Add it with: [Environment]::SetEnvironmentVariable('Path', `$env:Path + ';$Dest', 'User')"
    }
    & (Join-Path $Dest 'tptdyn.exe') --version
}
finally {
    Remove-Item -Recurse -Force -Path $Tmp -ErrorAction SilentlyContinue
}
