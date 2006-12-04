#!/bin/sh

# Run "prove" on all Perl programs

TESTS="aromatic.pl atom bond conversion data format"
TESTS="${TESTS} formula formalcharge internalcoord"
TESTS="${TESTS} iterators"
TESTS="${TESTS} invalidsmarts invalidsmiles math"
TESTS="${TESTS} mol residue rings"
TESTS="${TESTS} smarts smilesmatch unitcell"
TESTS="${TESTS} cml.sh test-set.sh"
PROVE=prove

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

${PROVE} "$@" ${TESTS}
