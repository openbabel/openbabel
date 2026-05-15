#!/usr/bin/env bash
# Windows bootstrap for cibuildwheel runs. The windows-latest runner ships
# with chocolatey, CMake, and Visual Studio build tools; we add swig and
# vcpkg-built dependencies (eigen, libxml2, zlib) so the CMake configure
# finds everything it needs.
set -ev

choco install -y --no-progress swig

# vcpkg ships pre-installed at C:\vcpkg on GitHub-hosted Windows runners.
if [ -d "/c/vcpkg" ]; then
  (cd /c/vcpkg && git pull origin master)
  /c/vcpkg/vcpkg.exe install eigen3:x64-windows libxml2:x64-windows zlib:x64-windows
fi

swig -version
cmake --version
