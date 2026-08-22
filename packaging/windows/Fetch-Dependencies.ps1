$ErrorActionPreference = 'Stop'

$Root = if ($PSScriptRoot) { Split-Path -Parent (Split-Path -Parent $PSScriptRoot) } else { $env:GITHUB_WORKSPACE }
$RunnerTemp = if ($env:RUNNER_TEMP) { $env:RUNNER_TEMP } else { Join-Path $Root 'build-windows' }

$DepsDir = Join-Path $RunnerTemp 'msvc-dependencies'
$ZlibRepo = Join-Path $RunnerTemp 'msvc-dependencies-zlib'

$DepsCommit = '15f6b268d2a0c0e4605f15848e4050a1127fc0a2'
$ZlibCommit = 'aeba12e73e82acf1ae94fbe275c35846e7ec4dfb'

function Checkout-Repository([string]$Url, [string]$Directory, [string]$Commit) {
    if (Test-Path $Directory) {
        Remove-Item $Directory -Recurse -Force
    }
    git clone --depth 1 $Url $Directory
    git -C $Directory fetch --depth 1 origin $Commit
    git -C $Directory checkout $Commit
}

Checkout-Repository 'https://github.com/openbabel/msvc-dependencies.git' $DepsDir $DepsCommit

# These headers shadow the MSVC headers with incompatible legacy versions.
Remove-Item (Join-Path $DepsDir 'include\stdint.h') -Force -ErrorAction SilentlyContinue
Remove-Item (Join-Path $DepsDir 'include\inttypes.h') -Force -ErrorAction SilentlyContinue

# The current dependency snapshot does not contain the zlib import library used
# by this Open Babel Windows build. Restore it from the known-good legacy commit.
Checkout-Repository 'https://github.com/openbabel/msvc-dependencies.git' $ZlibRepo $ZlibCommit

$TargetDir = Join-Path $DepsDir 'libs-common\x64'
$ZlibLib = Join-Path $ZlibRepo 'libs-common\x64\zlib1.lib'
$ZlibDll = Join-Path $ZlibRepo 'libs-common\x64\zlib1.dll'

foreach ($File in @($ZlibLib, $ZlibDll)) {
    if (-not (Test-Path $File)) {
        throw "Required zlib file is missing: $File"
    }
    Copy-Item $File $TargetDir -Force
}

Remove-Item $ZlibRepo -Recurse -Force

"DEPS_DIR=$DepsDir" >> $env:GITHUB_ENV
Write-Host "MSVC dependencies: $DepsDir"
Write-Host "Using zlib from commit $ZlibCommit"
Get-ChildItem $TargetDir -Filter 'zlib1.*' | Select-Object Name, Length
