/**********************************************************************
op2.cpp - Plugin to adds 2D coordinates using RDKit routines

Copyright (C) 2007 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
**********************************************************************

This code calls C++ routines in RDKit which are
  Copyright (C) 2003-2006 Rational Discovery LLC
    BSD license

***********************************************************************/

#include <openbabel/babelconfig.h>
#include <iostream>
#include<openbabel/op.h>
#include<openbabel/mol.h>
#include <RDKitConv.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <Geometry/point.h>
#include <GraphMol/conformer.h>
#include <GraphMol/molops.h>

#ifndef OBERROR
 #define OBERROR
#endif

namespace OpenBabel
{
  class Op2D : public OBOp //was OBERROR when with OpenBabelDLL
{
public:
  Op2D(const char* ID) : OBOp(ID, false){};
  const char* Description()
  {
    return "Generate 2D coordinates\n"
      "Uses RDKit http://www.rdkit.org";
  }
  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }

  virtual bool Do(OBBase* pOb, OpMap*, const char* OptionText);
};

Op2D theOp2D("2D"); //Global instance

/////////////////////////////////////////////////////////////////
bool Op2D::Do(OBBase* pOb, OpMap*, const char* OptionText)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  //Do the conversion
  try
  {
    RDKit::RWMol RDMol = OBMolToRWMol(pmol);
    RDKit::MolOps::sanitizeMol(RDMol); //initializes various internl parameters

    unsigned int ConformerID = RDDepict::compute2DCoords(RDMol);
    RDKit::Conformer confmer = RDMol.getConformer(ConformerID);
    for(int i=0; i<confmer.getNumAtoms(); ++i)
    {
      //transfer coordinates from the RDKit conformer to equivalent atoms in the OBMol
      RDGeom::Point3D atompos = confmer.getAtomPos(i);
      OBAtom* obat = pmol->GetAtom(i+1);
      obat->SetVector(atompos.x, atompos.y, 0);
    }
  }
  catch(...)
  {
    obErrorLog.ThrowError(__FUNCTION__, "Op2D failed with an exception, probably in RDKit Code" , obError);
    return false;
  }
  pmol->SetDimension(2); //No longer without coordinates!
  return true;
}

}//namespace
