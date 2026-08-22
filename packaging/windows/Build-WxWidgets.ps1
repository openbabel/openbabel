$ErrorActionPreference = 'Stop'

$RunnerTemp = if ($env:RUNNER_TEMP) { $env:RUNNER_TEMP } else { Join-Path $env:TEMP 'openbabel-windows' }
$WxRoot = Join-Path $RunnerTemp 'wxWidgets'
$WxVersion = 'v3.0.5'

if (Test-Path $WxRoot) {
    Remove-Item $WxRoot -Recurse -Force
}

git clone --depth 1 --branch $WxVersion https://github.com/wxWidgets/wxWidgets.git $WxRoot
if (-not (Test-Path "$WxRoot\include\wx\wx.h")) {
    throw 'wxWidgets source checkout is incomplete.'
}

Copy-Item "$WxRoot\include\wx\msw\setup0.h" "$WxRoot\include\wx\msw\setup.h" -Force

# Keep the proven CI build path: initialize the MSVC environment in cmd.exe and
# invoke nmake there. Calling vcvarsall from PowerShell can produce misleading
# "network path was not found" errors because of command-line quoting/context.
$VcVars = Join-Path $env:VCINSTALLDIR 'Auxiliary\Build\vcvarsall.bat'
if (-not (Test-Path $VcVars)) {
    throw "vcvarsall.bat was not found: $VcVars"
}

$WxBuildDir = Join-Path $WxRoot 'build\msw'
cmd.exe /d /s /c "call \"$VcVars\" x64 && cd /d \"$WxBuildDir\" && nmake /f makefile.vc BUILD=release SHARED=1 TARGET_CPU=X64"
if ($LASTEXITCODE -ne 0) {
    throw "wxWidgets build failed with exit code $LASTEXITCODE"
}

$Lib = Get-ChildItem "$WxRoot\lib" -Recurse -Filter 'wxmsw*.lib' -File |
    Where-Object { $_.FullName -match 'vc.*x64.*dll|vc.*dll.*x64|x64.*dll' } |
    Select-Object -First 1
if ($null -eq $Lib) {
    Get-ChildItem "$WxRoot\lib" -Recurse -File | Select-Object FullName, Length
    throw 'Could not locate the wxWidgets x64 release library.'
}

$WxLibDir = $Lib.Directory.FullName
$WxDlls = Get-ChildItem $WxLibDir -Filter 'wx*.dll' -File
if ($WxDlls.Count -eq 0) {
    throw "Could not locate wxWidgets runtime DLLs in $WxLibDir"
}

foreach ($Required in @('wxbase30u_vc_custom.dll', 'wxmsw30u_core_vc_custom.dll')) {
    if (-not ($WxDlls.Name -contains $Required)) {
        throw "Required wxWidgets runtime DLL is missing: $Required"
    }
}

"WXWIN=$WxRoot" >> $env:GITHUB_ENV
"WX_LIB_DIR=$WxLibDir" >> $env:GITHUB_ENV

Write-Host "wxWidgets root: $WxRoot"
Write-Host "wxWidgets library: $WxLibDir"
$WxDlls | Select-Object Name, Length
