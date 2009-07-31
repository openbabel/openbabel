#!/bin/sh

# Generate test results

TESTS="ffghemical ffuff ffmmff94 formula formalcharge rings smarts strip"

unset BABEL_LIBDIR
unset BABEL_DATADIR
if [ -d ../src/formats/.libs ]; then
    if [ "x${BABEL_LIBDIR}" = "x" ]; then
	BABEL_LIBDIR="`pwd`/../src/formats/.libs:`pwd`/../src/formats/xml/.libs"
	export BABEL_LIBDIR
    fi
    if [ "x${BABEL_DATADIR}" = "x" ]; then
	if [ "x${top_srcdir}" != "x" ]; then
	    BABEL_DATADIR="${top_srcdir}/data"
	else
	    BABEL_DATADIR="`pwd`/../data"
	fi
	export BABEL_DATADIR
    fi
fi

for test in ${TESTS}; do
 echo "Generating results for test $test"
 ./${test} -g
done
