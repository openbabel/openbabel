#!/bin/sh
# test Open Babel with CML, if the CML test suite exists

echo
echo "Testing CML support..."
if [ -f cmltest/test.sh ]; then
    cd cmltest
    source test.sh
    cd ..
else
    echo "CML tests not found. Skipping test."
    exit 77
fi
