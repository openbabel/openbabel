#!/bin/sh

# Run "prove" on all Perl programs
# Passes along any arguments to the prove tool
# so it can also be used for debugging or running in random order, etc.

# wrapper.sh --shuffle --debug
# wrapper.sh --verbose # display full output of tests while running

TESTS="aromatic.pl atom bond cansmi cmlreadfile conversion data"
TESTS="${TESTS} ffghemical ffmmff94 ffuff"
TESTS="${TESTS} format formula formalcharge inchi.pl inchi2.pl"
TESTS="${TESTS} internalcoord iterators" 
TESTS="${TESTS} invalidsmarts invalidsmiles"
TESTS="${TESTS} logp_psa math"
TESTS="${TESTS} mol residue rings"
TESTS="${TESTS} smarts smilesmatch unitcell"
TESTS="${TESTS} cml.sh test-set.sh"
PROVE=prove

echo "top srcdir: .${topsrcdir}."

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
