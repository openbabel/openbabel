#!/usr/bin/env bash

set -x

mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../
make
# only install if tests passed
make test && sudo make install
