/**********************************************************************
gen3d.cpp - A OBOp for generation of 3D coordinates (wrapper for OBBuilder)

Copyright (C) 2006-2007 by Tim Vandermeersch
          (C) 2007 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/builder.h>
#include <openbabel/forcefield.h>

namespace OpenBabel
{

class OpGen3D : public OBOp
{
public:
  OpGen3D(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return "Generate 3D coordinates"; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
};

/////////////////////////////////////////////////////////////////
OpGen3D theOpGen3D("gen3D"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpGen3D::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  OBBuilder builder;
  builder.Build(*pmol);
  pmol->SetDimension(3);

  OBForceField* pFF = OBForceField::FindForceField("MMFF94");
  if (!pFF)
    return true;

  pmol->AddHydrogens(false, true); // Add some hydrogens before running MMFF
  if (!pFF->Setup(*pmol)) {
    pFF = OBForceField::FindForceField("UFF");
    if (!pFF || !pFF->Setup(*pmol)) return true; // can't use either MMFF94 or UFF
  }

  // Since we only want a rough geometry, use distance cutoffs for VDW, Electrostatics
  pFF->EnableCutOff(true);
  pFF->SetVDWCutOff(10.0);
  pFF->SetElectrostaticCutOff(20.0);
  pFF->SetUpdateFrequency(10); // update non-bonded distances

  pFF->SteepestDescent(250, 1.0e-4);
  pFF->WeightedRotorSearch(200, 25);
  pFF->ConjugateGradients(250, 1.0e-6);
  pFF->UpdateCoordinates(*pmol);

  return true;
}
}//namespace
