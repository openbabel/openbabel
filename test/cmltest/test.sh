#!/bin/sh
# test Open Babel with CML

# Make sure we know the absolute path of the programs we're trying to 
# call, in case the test is called from somewhere else
unset BABEL_LIBDIR
if `env | grep ^builddir > /dev/null 2>&1`; then
   BABEL=$builddir/../tools/babel
   BABEL_LIBDIR="$builddir/../src/formats/.libs:$builddir/../src/formats/xml/.libs"
   export BABEL_LIBDIR
else
   builddir=..
   BABEL=../../tools/babel
   BABEL_LIBDIR="../../src/formats/.libs:../../src/formats/xml/.libs"
   export BABEL_LIBDIR
fi

# check to see if we have CML support or bail!
cmltest=`${BABEL} -Hcml | grep "not recognized"`
if [ "x${cmltest}" != "x" ]; then
    echo "1..0 # skipping - CML format not loaded"
    exit
fi

ROUNDTRIP=${builddir}/roundtrip

# test input
# CML2 with array
echo "1..6"
echo "# CML2 with array"
if (${BABEL} -icml cs2a.cml -omdl cs2a.mol && ${ROUNDTRIP} cs2a.cml cs2a.mol) then
    echo "ok 1"
else
    echo "not ok 1"
fi

#  3D molecules in SDF 
#  CML2 with XML version
echo "# CML2 with XML version"
if (${BABEL} -isdf 3d.head.sdf -ocml 3d.head.2.cml -x2v && ${ROUNDTRIP} 3d.head.sdf 3d.head.2.cml) then
    echo "ok 2"
else
    echo "not ok 2"
fi

#  CML1 with DOCTYPE
echo "# CML1 with DOCTYPE"
if (${BABEL} -isdf cs2a.mol -ocml cs2a.mol.cml -x1d && ${ROUNDTRIP} cs2a.mol cs2a.mol.cml); then
    echo "ok 3"
else
    echo "not ok 3"
fi

#  CML2 arrays with namespaces (large)
echo "# CML2 arrays with namespaces"
# PR#1486678
if (${BABEL} -isdf 3d.head.sdf -ocml 3d.head.2an.cml -x2an && ${ROUNDTRIP} 3d.head.sdf 3d.head.2an.cml); then
    echo "ok 4"
else
    echo "not ok 4"
fi

#  roundtripping; arguments are fileroot; input format; input suffix
#  2d MDL to CML and back again through all main variants
echo "# Roundtripping from 2D MDL Molfile to CML and back"
if ./cycleformats.sh nsc2dmol mol mdl; then
    echo "ok 5"
else
    echo "not ok 5"
fi

#  3d MDL to CML and back again through all main variants
echo "# Roundtripping from 3D MDL Molfile to CML and back"
if ./cycleformats.sh nsc3dmol mol mdl; then
    echo "ok 6"
else
    echo "not ok 6"
fi
