$ErrorActionPreference = 'Stop'

$MaintenanceDir = if ($env:MAINTENANCE_DIR) {
    $env:MAINTENANCE_DIR
} else {
    throw 'MAINTENANCE_DIR is not set.'
}

$DistDir = Join-Path $MaintenanceDir 'windows\for_dist'
if (-not (Test-Path $DistDir)) {
    throw "Windows distribution directory was not found: $DistDir"
}

$Output = Join-Path $DistDir 'vc_redist.x64.exe'
$Uri = 'https://aka.ms/vc14/vc_redist.x64.exe'

Invoke-WebRequest -Uri $Uri -OutFile $Output

if (-not (Test-Path $Output)) {
    throw "VC++ redistributable was not downloaded: $Output"
}

Write-Host "VC++ redistributable: $Output"
