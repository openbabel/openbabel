#!/usr/bin/env bash
# Runs inside the manylinux container before any wheel is built.
# manylinux_2_28 is based on AlmaLinux 8, so we get dnf/yum.
set -ev

if command -v dnf >/dev/null 2>&1; then
  dnf -y install swig perl eigen3-devel libxml2-devel zlib-devel
elif command -v yum >/dev/null 2>&1; then
  yum -y install swig perl eigen3-devel libxml2-devel zlib-devel
else
  echo "No supported package manager (dnf/yum) found in manylinux container." >&2
  exit 1
fi

swig -version
cmake --version
