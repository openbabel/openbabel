#!/bin/sh

# Run "prove" on all Perl programs

TESTS="aromatest.pl atom bond conversion data format"
TESTS="${TESTS} internalcoord matrix mol residue rings"
TESTS="${TESTS} smarts unitcell"
TESTS="${TESTS} cml.sh test-set.sh"
PROVE=prove

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

prove ${TESTS}
