$ErrorActionPreference = 'Stop'

$BuildDir = if ($env:BUILD_DIR) { $env:BUILD_DIR } else { Join-Path $env:RUNNER_TEMP 'openbabel-build' }
$DepsDir = if ($env:DEPS_DIR) { $env:DEPS_DIR } else { Join-Path $env:RUNNER_TEMP 'msvc-dependencies' }
$MaintenanceDir = if ($env:MAINTENANCE_DIR) { $env:MAINTENANCE_DIR } else { Join-Path $env:RUNNER_TEMP 'openbabel-maintenance' }
$WxLibDir = if ($env:WX_LIB_DIR) { $env:WX_LIB_DIR } else { throw 'WX_LIB_DIR is not set.' }
$SourceDir = if ($env:GITHUB_WORKSPACE) { $env:GITHUB_WORKSPACE } else { (Get-Location).Path }

$DistDir = Join-Path $MaintenanceDir 'windows\for_dist'
$Nsi = Join-Path $DistDir 'NSISScriptToCreateInstaller.nsi'

$Matches = Select-String -Path (Join-Path $SourceDir 'CMakeLists.txt') -Pattern 'set\(BABEL_(MAJ|MIN|PATCH)_VER\s+([0-9]+)\)' -AllMatches
$Version = ($Matches.Matches | ForEach-Object { $_.Groups[2].Value }) -join '.'
if ([string]::IsNullOrWhiteSpace($Version)) {
    throw 'Could not determine Open Babel version.'
}

# Keep the maintenance repository's NSIS script as the packaging contract for
# now. Patch only the options that differ from this x64 native-only build.
$Text = Get-Content $Nsi -Raw
$Text = $Text -replace '!define OBVersion\s+[^\r\n]+', "!define OBVersion $Version"

foreach ($Pattern in @(
    '(?m)^\s*File \$\{SourceDir\}\\scripts\\java\\openbabel\.jar\s*\r?\n',
    '(?m)^\s*File \$\{BuildDir\}\\bin\\Release\\openbabel_java\.dll\s*\r?\n',
    '(?m)^\s*File \$\{BuildDir\}\\bin\\Release\\openbabel_csharp\.dll\s*\r?\n',
    '(?m)^\s*File \$\{BuildDir\}\\bin\\Release\\OBDotNet\.dll\s*\r?\n'
)) {
    $Text = $Text -replace $Pattern, ''
}

$Text = $Text -replace '(?s)\s*StrCmp \$\{Arch\} "i386" 0 archIs64\s*File vc_redist\.x86\.exe\s*Goto done\s*archIs64:\s*File vc_redist\.x64\.exe\s*done:', ([Environment]::NewLine + '  File vc_redist.x64.exe')
$Text = $Text -replace '(?s)\s*StrCmp \$\{Arch\} "i386" 0 vcredist_for_x64\s*ExecWait ''"\$INSTDIR/vc_redist\.x86\.exe" /quiet''\s*Goto vcredist_done\s*vcredist_for_x64:\s*ExecWait ''"\$INSTDIR/vc_redist\.x64\.exe" /quiet''\s*vcredist_done:', ([Environment]::NewLine + '  ExecWait ''"$INSTDIR/vc_redist.x64.exe" /quiet''')

Set-Content -Path $Nsi -Value $Text -NoNewline

# wxWidgets is dynamically linked. Reuse the legacy NSIS script's existing
# dependency collection by staging its complete runtime set in DepsDir.
$DepsDllDir = Join-Path $DepsDir 'libs-common\x64'
$WxDlls = Get-ChildItem $WxLibDir -Filter 'wx*.dll' -File
if ($WxDlls.Count -eq 0) {
    throw "Could not locate wxWidgets runtime DLLs in $WxLibDir"
}
foreach ($Required in @('wxbase30u_vc_custom.dll', 'wxmsw30u_core_vc_custom.dll')) {
    if (-not ($WxDlls.Name -contains $Required)) {
        throw "Required wxWidgets runtime DLL is missing: $Required"
    }
}
foreach ($Dll in $WxDlls) {
    Copy-Item $Dll.FullName $DepsDllDir -Force
}

$ReleaseDir = Join-Path $BuildDir 'bin\Release'
foreach ($Required in @('obgui.exe', 'openbabel-3.dll')) {
    if (-not (Test-Path (Join-Path $ReleaseDir $Required))) {
        throw "Required installer payload is missing: $Required"
    }
}

$Output = Join-Path $SourceDir "OpenBabel-$Version-x64.exe"
Push-Location $DistDir
try {
    makensis.exe `
        "/DSourceDir=$SourceDir" `
        "/DBuildDir=$BuildDir" `
        '/DArch=x64' `
        "/DDepsDir=$DepsDir" `
        "/DmyOutFile=$Output" `
        'NSISScriptToCreateInstaller.nsi'
    if ($LASTEXITCODE -ne 0) {
        throw "makensis.exe failed with exit code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}

if (-not (Test-Path $Output)) {
    throw "Installer was not created: $Output"
}

"INSTALLER=$Output" >> $env:GITHUB_ENV
"VERSION=$Version" >> $env:GITHUB_ENV
Write-Host "Installer: $Output"
