#!/usr/bin/env bash
# cibuildwheel's repair step for Linux. We unpack the freshly-built wheel
# so auditwheel can see libopenbabel and the plugin .so files (which live
# inside the wheel at openbabel/ and openbabel/plugins/) when resolving
# dependent ELF symbols.
set -ev

WHEEL="$1"
DEST_DIR="$2"

WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

# Extract so auditwheel can find sibling shared objects on LD_LIBRARY_PATH.
unzip -q "$WHEEL" -d "$WORK_DIR"

export LD_LIBRARY_PATH="$WORK_DIR/openbabel:$WORK_DIR/openbabel/plugins${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

auditwheel repair --plat "$AUDITWHEEL_PLAT" -w "$DEST_DIR" "$WHEEL"
