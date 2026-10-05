$ErrorActionPreference = 'Stop'

$BuildDir = if ($env:BUILD_DIR) { $env:BUILD_DIR } else { Join-Path $env:RUNNER_TEMP 'openbabel-build' }
$DepsDir = if ($env:DEPS_DIR) { $env:DEPS_DIR } else { Join-Path $env:RUNNER_TEMP 'msvc-dependencies' }
$Workspace = if ($env:GITHUB_WORKSPACE) { $env:GITHUB_WORKSPACE } else { (Get-Location).Path }

$Swig = (Get-Command swig.exe -ErrorAction Stop).Source

# Open Babel's MSVC prebuilt-binary mode expects the system InChI header and
# library. The dependency snapshot already provides libinchi.lib/libinchi.dll.
$InchiIncludeDir = Join-Path $Workspace 'include\inchi'
$InchiHeader = Join-Path $InchiIncludeDir 'inchi_api.h'
$InchiLibrary = Join-Path $DepsDir 'libs-common\x64\libinchi.lib'
$InchiRuntime = Join-Path $DepsDir 'libs-common\x64\libinchi.dll'
if (-not (Test-Path $InchiHeader)) { throw "InChI API header was not found: $InchiHeader" }
if (-not (Test-Path $InchiLibrary)) { throw "InChI import library was not found: $InchiLibrary" }
if (-not (Test-Path $InchiRuntime)) { throw "InChI runtime DLL was not found: $InchiRuntime" }

cmake -S $Workspace -B $BuildDir `
  -G 'Visual Studio 17 2022' -A x64 `
  -DCMAKE_INSTALL_PREFIX="$BuildDir\install" `
  -DOB_USE_PREBUILT_BINARIES=ON `
  -DMINIMAL_BUILD=OFF `
  -DBUILD_GUI=ON `
  -DwxWidgets_ROOT_DIR="$env:WXWIN" `
  -DwxWidgets_LIB_DIR="$env:WX_LIB_DIR" `
  -DwxWidgets_CONFIGURATION=mswu `
  -DwxWidgets_USE_STATIC=OFF `
  -DLIBXML2_LIBRARIES="$DepsDir\libs-common\x64\libxml2.lib" `
  -DLIBXML2_INCLUDE_DIR="$DepsDir\include" `
  -DCAIRO_INCLUDE_DIRS="$DepsDir\include\cairo" `
  -DCAIRO_LIBRARIES="$DepsDir\libs-common\x64\cairo.lib" `
  -DZLIB_LIBRARY="$DepsDir\libs-common\x64\zlib1.lib" `
  -DZLIB_INCLUDE_DIR="$DepsDir\include" `
  -DXDR_LIBRARY="$DepsDir\libs-common\x64\xdr.lib" `
  -DXDR_INCLUDE_DIR="$DepsDir\include" `
  -DEIGEN3_INCLUDE_DIR="$DepsDir\include" `
  -DJSON_LIBRARY="$DepsDir\libs-vs12\x64\jsoncpp.lib" `
  -DINCHI_INCLUDE_DIR="$InchiIncludeDir" `
  -DINCHI_LIBRARY="$InchiLibrary" `
  -DRUN_SWIG=ON `
  -DSWIG_EXECUTABLE="$Swig" `
  -DJAVA_BINDINGS=OFF `
  -DCSHARP_BINDINGS=OFF `
  -DWITH_INCHI=ON `
  -DOPENBABEL_USE_SYSTEM_INCHI=ON `
  -DENABLE_TESTS=OFF

if ($LASTEXITCODE -ne 0) {
    throw "CMake configure failed with exit code $LASTEXITCODE"
}

Get-Content "$BuildDir\CMakeCache.txt" |
    Select-String 'wxWidgets|WX|ZLIB|INCHI|JAVA|CSHARP|BUILD_GUI|MINIMAL_BUILD'
