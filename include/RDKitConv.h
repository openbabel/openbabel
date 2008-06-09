#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <graphmol/RWMol.h>
#include <graphmol/Atom.h>

///Convert OpenBabel OBMol to and from RGKit molecules
RDKit::RWMol OBMolToRWMol(OpenBabel::OBMol* pOBMol);

//! \file RDKitConv.h
//! \brief Allow conversion from OBMol to RDKit RWMol