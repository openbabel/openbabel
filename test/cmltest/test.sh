#!/bin/sh
# test Open Babel with CML

BABEL=../../src/babel

# test input
# CML2 with array
echo "CML2 with array"
$BABEL -icml cs2a.cml -omdl cs2a.mol 
../roundtrip cs2a.cml cs2a.mol
#  3D molecules in SDF 
#  CML2 with XML version
echo "CML2 with XML version"
$BABEL -isdf 3d.head.sdf -ocml 3d.head.2.cml -x2v 
../roundtrip 3d.head.sdf 3d.head.2.cml
#  CML1 with DOCTYPE
echo "CML1 with DOCTYPE"
$BABEL -isdf cs2a.mol -ocml cs2a.mol.cml -x1d
../roundtrip cs2a.mol cs2a.mol.cml
#  CML2 arrays with namespaces (large)
echo "CML2 arrays with namespaces"
$BABEL -isdf 3d.head.sdf -ocml 3d.head.2an.cml -x2an 
../roundtrip 3d.head.sdf 3d.head.2an.cml

#  roundtripping; arguments are fileroot; input format; input suffix
#  2d MDL to CML and back again through all main variants
./roundtrip.sh nsc2dmol mol mdl
#  3d MDL to CML and back again through all main variants
./roundtrip.sh nsc3dmol mol mdl

