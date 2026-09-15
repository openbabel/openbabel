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
    if ($LASTEXITCODE -ne 0) {
        throw "git clone failed for $Url"
    }

    git -C $Directory fetch --depth 1 origin $Commit
    if ($LASTEXITCODE -ne 0) {
        throw "git fetch failed for commit $Commit"
    }

    git -C $Directory checkout $Commit
    if ($LASTEXITCODE -ne 0) {
        throw "git checkout failed for commit $Commit"
    }

    $Head = git -C $Directory rev-parse HEAD
    if ($LASTEXITCODE -ne 0 -or $Head -ne $Commit) {
        throw "Checked out commit $Head does not match requested commit $Commit"
    }
}

Checkout-Repository 'https://github.com/openbabel/msvc-dependencies.git' $DepsDir $DepsCommit

# These headers shadow the MSVC headers with incompatible legacy versions.
Remove-Item (Join-Path $DepsDir 'include\stdint.h') -Force -ErrorAction SilentlyContinue
Remove-Item (Join-Path $DepsDir 'include\inttypes.h') -Force -ErrorAction SilentlyContinue

# Restore only the import library needed by the Open Babel link step. Do not
# restore the legacy DLL: libpng16.dll must use the matching current zlib DLL.
Checkout-Repository 'https://github.com/openbabel/msvc-dependencies.git' $ZlibRepo $ZlibCommit

$TargetDir = Join-Path $DepsDir 'libs-common\x64'
$ZlibLib = Join-Path $ZlibRepo 'libs-common\x64\zlib1.lib'
if (-not (Test-Path $ZlibLib)) {
    throw "Required zlib import library is missing: $ZlibLib"
}
Copy-Item $ZlibLib $TargetDir -Force
Remove-Item $ZlibRepo -Recurse -Force

$ZlibRuntime = Join-Path $TargetDir 'zlib1.dll'
if (-not (Test-Path $ZlibRuntime)) {
    throw "Current dependency snapshot lacks zlib runtime: $ZlibRuntime"
}

"DEPS_DIR=$DepsDir" >> $env:GITHUB_ENV
Write-Host "MSVC dependencies: $DepsDir"
Write-Host "Using zlib import library from commit $ZlibCommit"
Get-ChildItem $TargetDir -Filter 'zlib1.*' | Select-Object Name, Length
