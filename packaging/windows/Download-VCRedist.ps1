$ErrorActionPreference = 'Stop'

$OutputDir = if ($env:RUNNER_TEMP) {
    Join-Path $env:RUNNER_TEMP 'openbabel-windows'
} elseif ($env:BUILD_DIR) {
    Join-Path $env:BUILD_DIR 'packaging'
} else {
    throw 'Neither RUNNER_TEMP nor BUILD_DIR is available.'
}

New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null
$Output = Join-Path $OutputDir 'vc_redist.x64.exe'
$Uri = 'https://aka.ms/vc14/vc_redist.x64.exe'

Invoke-WebRequest -Uri $Uri -OutFile $Output

if (-not (Test-Path $Output)) {
    throw "VC++ redistributable was not downloaded: $Output"
}

"VCREDIST=$Output" >> $env:GITHUB_ENV
Write-Host "VC++ redistributable: $Output"
