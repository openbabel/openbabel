#!/bin/sh
# test OB with CML
# $1 is inputfile root
# $2 is inputfile suffix
# $3 is inputfile type
# e.g. cycleformats.sh foo mol mdl

# Make sure we know the absolute path of the programs we're trying to 
# call, in case the test is called from somewhere eles
if `env | grep ^builddir > /dev/null 2>&1`; then
   BABEL=$builddir/../tools/babel
else
   builddir=..
   BABEL=../../tools/babel
fi

# CML1 output
$BABEL -i$3 $1.$2 -ocml $1.1.cml -x1v || exit 1;
$builddir/roundtrip $1.$2 $1.1.cml || exit 1;
# CML2 output
$BABEL -i$3 $1.$2 -ocml $1.2.cml -x2v || exit 1;
$builddir/roundtrip $1.$2 $1.2.cml || exit 1;
# CML1+array output
$BABEL -i$3 $1.$2 -ocml $1.a1.cml -xa1v || exit 1;
$builddir/roundtrip $1.$2 $1.a1.cml || exit 1;
# CML2+array output
$BABEL -i$3 $1.$2 -ocml $1.a2.cml -xa2v || exit 1;
$builddir/roundtrip $1.$2 $1.a2.cml || exit 1;

# roundtrip to MOL; should be identical
$BABEL -icml $1.1.cml  -o$3 $1.1.$2 -x2v || exit 1;
$builddir/roundtrip $1.1.cml $1.1.$2 || exit 1;
$BABEL -icml $1.2.cml  -o$3 $1.2.$2 -x2v || exit 1;
$builddir/roundtrip $1.2.cml $1.2.$2 || exit 1;
$BABEL -icml $1.a1.cml -o$3 $1.a1.$2 -x2v || exit 1;
$builddir/roundtrip $1.a1.cml $1.a1.$2 || exit 1;
$BABEL -icml $1.a2.cml -o$3 $1.a2.$2 -x2v || exit 1;
$builddir/roundtrip $1.a2.cml $1.a2.$2 || exit 1;

# And check to make sure the four $2 files are the same !
$builddir/roundtrip $1.1.$2 $1.2.$2 || exit 1;
$builddir/roundtrip $1.1.$2 $1.a1.$2 || exit 1;
$builddir/roundtrip $1.1.$2 $1.a2.$2 || exit 1;
