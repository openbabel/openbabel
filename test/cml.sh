#!/bin/sh
# test Open Babel with CML, if the CML test suite exists

if [ -f cmltest/test.sh ]; then
    cd cmltest
    source test.sh
    cd ..
fi
