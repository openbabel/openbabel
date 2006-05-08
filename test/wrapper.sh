#!/bin/sh

# Run "prove" on all Perl programs

TESTS="matrix.pl unitcell.pl rings.pl smarts.pl aromatest.pl cml.pl test-set.pl"
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
