$ErrorActionPreference = 'Stop'

$BuildDir = if ($env:BUILD_DIR) { $env:BUILD_DIR } else { Join-Path $env:RUNNER_TEMP 'openbabel-build' }
$DepsDir = if ($env:DEPS_DIR) { $env:DEPS_DIR } else { Join-Path $env:RUNNER_TEMP 'msvc-dependencies' }
$VCRedist = if ($env:VCREDIST) { $env:VCREDIST } else { throw 'VCREDIST is not set.' }
$WxLibDir = if ($env:WX_LIB_DIR) { $env:WX_LIB_DIR } else { throw 'WX_LIB_DIR is not set.' }
$SourceDir = if ($env:GITHUB_WORKSPACE) { $env:GITHUB_WORKSPACE } else { (Get-Location).Path }
$InstallerScript = Join-Path $SourceDir 'packaging\windows\installer.nsi'

if (-not (Test-Path $InstallerScript)) { throw "Installer script was not found: $InstallerScript" }

$Matches = Select-String -Path (Join-Path $SourceDir 'CMakeLists.txt') -Pattern 'set\(BABEL_(MAJ|MIN|PATCH)_VER\s+([0-9]+)\)' -AllMatches
$Version = ($Matches.Matches | ForEach-Object { $_.Groups[2].Value }) -join '.'
if ([string]::IsNullOrWhiteSpace($Version)) { throw 'Could not determine Open Babel version.' }

$DepsDllDir = Join-Path $DepsDir 'libs-common\x64'
if (-not (Test-Path $DepsDllDir)) { throw "MSVC dependency directory was not found: $DepsDllDir" }

# wxWidgets is built as a shared library. Stage its complete x64 runtime set
# next to the other installer DLL dependencies so installer.nsi can package it.
$WxDlls = Get-ChildItem $WxLibDir -Filter 'wx*.dll' -File
if ($WxDlls.Count -eq 0) { throw "No wxWidgets runtime DLLs found in $WxLibDir" }
foreach ($Required in @('wxbase30u_vc_custom.dll', 'wxmsw30u_core_vc_custom.dll', 'wxmsw30u_adv_vc_custom.dll')) {
    if (-not ($WxDlls.Name -contains $Required)) { throw "Required wxWidgets runtime DLL is missing: $Required" }
}
foreach ($Dll in $WxDlls) {
    Copy-Item $Dll.FullName $DepsDllDir -Force
    Write-Host "Staging wxWidgets runtime: $($Dll.Name)"
}

$ReleaseDir = Join-Path $BuildDir 'bin\Release'
foreach ($Required in @('obgui.exe', 'openbabel-3.dll')) {
    if (-not (Test-Path (Join-Path $ReleaseDir $Required))) { throw "Required installer payload is missing: $Required" }
}

$Output = Join-Path $SourceDir "OpenBabel-$Version-x64.exe"
makensis.exe "/DSourceDir=$SourceDir" "/DBuildDir=$BuildDir" "/DDepsDir=$DepsDir" "/DVCRedist=$VCRedist" '/DArch=x64' "/DOBVersion=$Version" "/DmyOutFile=$Output" $InstallerScript
if ($LASTEXITCODE -ne 0) { throw "makensis.exe failed with exit code $LASTEXITCODE" }
if (-not (Test-Path $Output)) { throw "Installer was not created: $Output" }

"INSTALLER=$Output" >> $env:GITHUB_ENV
"VERSION=$Version" >> $env:GITHUB_ENV
Write-Host "Installer: $Output"
