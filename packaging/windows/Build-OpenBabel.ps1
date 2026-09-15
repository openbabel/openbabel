$ErrorActionPreference = 'Stop'

$BuildDir = if ($env:BUILD_DIR) { $env:BUILD_DIR } else { Join-Path $env:RUNNER_TEMP 'openbabel-build' }

cmake --build $BuildDir --config Release --parallel
if ($LASTEXITCODE -ne 0) {
    throw "Open Babel build failed with exit code $LASTEXITCODE"
}

$Gui = Get-ChildItem $BuildDir -Filter 'obgui.exe' -Recurse | Select-Object -First 1
if ($null -eq $Gui) {
    throw 'obgui.exe was not built.'
}

$Native = Join-Path $BuildDir 'bin\Release\openbabel-3.dll'
if (-not (Test-Path $Native)) {
    throw "Open Babel runtime was not built: $Native"
}

Write-Host "Built GUI: $($Gui.FullName)"
Write-Host "Built runtime: $Native"
