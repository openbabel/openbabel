$ErrorActionPreference = 'Stop'

$RunnerTemp = if ($env:RUNNER_TEMP) { $env:RUNNER_TEMP } else { Join-Path $env:TEMP 'openbabel-windows' }
$MaintenanceDir = Join-Path $RunnerTemp 'openbabel-maintenance'
$Commit = 'aa42db57cf9bae8dacf44b99e958d9f2a6a578e3'

if (Test-Path $MaintenanceDir) {
    Remove-Item $MaintenanceDir -Recurse -Force
}

git clone --depth 1 https://github.com/openbabel/maintenance.git $MaintenanceDir
git -C $MaintenanceDir fetch --depth 1 origin $Commit
git -C $MaintenanceDir checkout $Commit

$DistDir = Join-Path $MaintenanceDir 'windows\for_dist'
if (-not (Test-Path (Join-Path $DistDir 'NSISScriptToCreateInstaller.nsi'))) {
    throw 'Legacy NSIS installer script was not found.'
}

"MAINTENANCE_DIR=$MaintenanceDir" >> $env:GITHUB_ENV
Write-Host "Maintenance assets: $MaintenanceDir"
