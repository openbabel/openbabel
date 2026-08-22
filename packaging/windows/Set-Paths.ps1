$ErrorActionPreference = 'Stop'

if (-not $env:RUNNER_TEMP) {
    throw 'RUNNER_TEMP is not set. Set it when running this script outside GitHub Actions.'
}

"BUILD_DIR=$env:RUNNER_TEMP\openbabel-build" >> $env:GITHUB_ENV
"DEPS_DIR=$env:RUNNER_TEMP\msvc-dependencies" >> $env:GITHUB_ENV
"MAINTENANCE_DIR=$env:RUNNER_TEMP\openbabel-maintenance" >> $env:GITHUB_ENV
"WXWIN=$env:RUNNER_TEMP\wxWidgets" >> $env:GITHUB_ENV

Write-Host "BUILD_DIR=$env:RUNNER_TEMP\openbabel-build"
Write-Host "DEPS_DIR=$env:RUNNER_TEMP\msvc-dependencies"
Write-Host "MAINTENANCE_DIR=$env:RUNNER_TEMP\openbabel-maintenance"
Write-Host "WXWIN=$env:RUNNER_TEMP\wxWidgets"
