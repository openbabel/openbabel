#!/bin/sh
# test Open Babel with CML, if the CML test suite exists

builddir=`pwd`
export builddir
if `env | grep ^srcdir > /dev/null 2>&1`; then
   cmltestdir=$srcdir/cmltest
else
   srcdir=`pwd`
   export srcdir
   cmltestdir=./cmltest
   export cmltestdir
fi

echo
echo "# Testing CML support..."
if [ -f $cmltestdir/test.sh ]; then
    echo "# This test will take quite some time!"
    echo "# and it currently always passes, so check the results file"
    (cd $cmltestdir; source test.sh 2>/dev/null)
else
    echo "1..0 # skipping - CML test set not found"
fi
