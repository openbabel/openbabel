#!/bin/sh
# test Open Babel with CML

# Make sure we know the absolute path of the programs we're trying to 
# call, in case the test is called from somewhere eles
unset BABEL_LIBDIR
if `env | grep ^builddir > /dev/null 2>&1`; then
   BABEL=$builddir/../src/babel
   BABEL_LIBDIR="$builddir/../src/formats/.libs:$builddir/../src/formats/xml/.libs"
   export BABEL_LIBDIR
else
   builddir=..
   BABEL=../../src/babel
   BABEL_LIBDIR="../../src/formats/.libs:../../src/formats/xml/.libs"
   export BABEL_LIBDIR
fi

# check to see if we have CML support or bail!
cmltest=`${BABEL} -Hcml | grep "not recognized"`
if [ "x${cmltest}" != "x" ]; then
    echo "CML format not loaded. Skipping tests."
    exit 77
fi

ROUNDTRIP=${builddir}/roundtrip

# test input
# CML2 with array
echo "CML2 with array"
${BABEL} -icml cs2a.cml -omdl cs2a.mol 
${ROUNDTRIP} cs2a.cml cs2a.mol
#  3D molecules in SDF 
#  CML2 with XML version
echo "CML2 with XML version"
${BABEL} -isdf 3d.head.sdf -ocml 3d.head.2.cml -x2v 
${ROUNDTRIP} 3d.head.sdf 3d.head.2.cml
#  CML1 with DOCTYPE
echo "CML1 with DOCTYPE"
${BABEL} -isdf cs2a.mol -ocml cs2a.mol.cml -x1d
${ROUNDTRIP} cs2a.mol cs2a.mol.cml
#  CML2 arrays with namespaces (large)
echo "CML2 arrays with namespaces"
${BABEL} -isdf 3d.head.sdf -ocml 3d.head.2an.cml -x2an 
${ROUNDTRIP} 3d.head.sdf 3d.head.2an.cml

#  roundtripping; arguments are fileroot; input format; input suffix
#  2d MDL to CML and back again through all main variants
./roundtrip.sh nsc2dmol mol mdl
#  3d MDL to CML and back again through all main variants
./roundtrip.sh nsc3dmol mol mdl

