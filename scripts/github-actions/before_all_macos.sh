#!/usr/bin/env bash
# macOS bootstrap for cibuildwheel runs. The macos-15-intel and macos-latest
# (arm64) runners already have CMake; we just need swig and eigen on top.
set -ev

brew install swig eigen
# libxml2 and zlib ship with macOS but pkg-config may need them surfaced.
brew install libxml2 || true

swig -version
cmake --version
