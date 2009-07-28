#!/bin/sh

# Run "prove" on all Perl programs
# Passes along any arguments to the prove tool
# so it can also be used for debugging or running in random order, etc.

# wrapper.sh --shuffle --debug
# wrapper.sh --verbose # display full output of tests while running

TESTS="aromatic.pl atom bond cansmi cmlreadfile conversion data"
#TESTS="${TESTS} ffghemical ffmmff94 ffuff"
TESTS="${TESTS} format formula formalcharge"
TESTS="${TESTS} internalcoord iterators" 
TESTS="${TESTS} invalidsmarts invalidsmiles"
TESTS="${TESTS} logp_psa math"
TESTS="${TESTS} mol phmodel residue rings"
TESTS="${TESTS} smarts smilesmatch strip unitcell"
TESTS="${TESTS} zipstream cml.sh test-set.sh"
if [ "x${srcdir}" != "x" ]; then
  TESTS="${TESTS} ${srcdir}/inchi.pl ${srcdir}/inchi2.pl"
else
  TESTS="${TESTS} inchi.pl inchi2.pl"
fi
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
	      if [ "x${srcdir}" != "x" ]; then
	          BABEL_DATADIR="${srcdir}/../data"
		  TESTDATADIR="${srcdir}/files"
	      else
	          BABEL_DATADIR="`pwd`/../data"
	          BABEL_DATADIR="`pwd`/files"
	      fi
	      export BABEL_DATADIR
	      export TESTDATADIR
    fi
fi

${PROVE} "$@" ${TESTS}
