$ErrorActionPreference = 'Stop'

$BuildDir = if ($env:BUILD_DIR) { $env:BUILD_DIR } else { Join-Path $env:RUNNER_TEMP 'openbabel-build' }
$DepsDir = if ($env:DEPS_DIR) { $env:DEPS_DIR } else { Join-Path $env:RUNNER_TEMP 'msvc-dependencies' }
$VCRedist = if ($env:VCREDIST) { $env:VCREDIST } else { throw 'VCREDIST is not set.' }
$SourceDir = if ($env:GITHUB_WORKSPACE) { $env:GITHUB_WORKSPACE } else { (Get-Location).Path }
$InstallerScript = Join-Path $SourceDir 'packaging\windows\installer.nsi'

if (-not (Test-Path $InstallerScript)) {
    throw "Installer script was not found: $InstallerScript"
}

$Matches = Select-String -Path (Join-Path $SourceDir 'CMakeLists.txt') -Pattern 'set\(BABEL_(MAJ|MIN|PATCH)_VER\s+([0-9]+)\)' -AllMatches
$Version = ($Matches.Matches | ForEach-Object { $_.Groups[2].Value }) -join '.'
if ([string]::IsNullOrWhiteSpace($Version)) {
    throw 'Could not determine Open Babel version.'
}

$DepsDllDir = Join-Path $DepsDir 'libs-common\x64'
if (-not (Test-Path $DepsDllDir)) {
    throw "MSVC dependency directory was not found: $DepsDllDir"
}

$ReleaseDir = Join-Path $BuildDir 'bin\Release'
foreach ($Required in @('obgui.exe', 'openbabel-3.dll')) {
    if (-not (Test-Path (Join-Path $ReleaseDir $Required))) {
        throw "Required installer payload is missing: $Required"
    }
}

$Output = Join-Path $SourceDir "OpenBabel-$Version-x64.exe"

makensis.exe `
    "/DSourceDir=$SourceDir" `
    "/DBuildDir=$BuildDir" `
    "/DDepsDir=$DepsDir" `
    "/DVCRedist=$VCRedist" `
    '/DArch=x64' `
    "/DOBVersion=$Version" `
    "/DmyOutFile=$Output" `
    $InstallerScript

if ($LASTEXITCODE -ne 0) {
    throw "makensis.exe failed with exit code $LASTEXITCODE"
}

if (-not (Test-Path $Output)) {
    throw "Installer was not created: $Output"
}

"INSTALLER=$Output" >> $env:GITHUB_ENV
"VERSION=$Version" >> $env:GITHUB_ENV
Write-Host "Installer: $Output"
