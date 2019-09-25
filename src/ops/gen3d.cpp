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
#include <openbabel/distgeom.h>
#include <openbabel/forcefield.h>

#include <cstdlib> // needed for strtol and gcc 4.8

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

  // As with gen2D, we need to perceive the stereo if coming from 0D.
  // Otherwise, unspecified cis/trans stereobonds become specified.
  if (pmol->GetDimension() == 0) {
    pmol->UnsetFlag(OB_CHIRALITY_MOL);
    StereoFrom0D(pmol);
  }

  // 1 is best quality, slowest
  // 2 is good quality, slow
  // 3 is balance   (FF cleanup + FastRotorSearch)
  // 4 is fast      (OBBuilder + FF cleanup)
  // 5 is fastest   (only OBBuilder)
  // 6 is DistanceGeometry
  int speed;
  bool useDistGeom = false;

  // first try converting OptionText to an integer
  char *endptr;
  speed = strtol(OptionText, &endptr, 10);
  if (endptr == OptionText) { // not a number
    speed = 3; // we'll default to balanced
    // but let's also check if it's words like "fast" or "best"
    if (strncasecmp(OptionText, "fastest", 7) == 0)
      speed = 5;
    else if (strncasecmp(OptionText, "fast", 4) == 0) // already matched fastest
      speed = 4;
    else if (strncasecmp(OptionText, "med", 3) == 0) // or medium
      speed = 3;
    else if ( (strncasecmp(OptionText, "slowest", 7) == 0)
             || (strncasecmp(OptionText, "best", 4) == 0) )
      speed = 1;
    else if ( (strncasecmp(OptionText, "slow", 4) == 0)
              || (strncasecmp(OptionText, "better", 6) == 0) )
      speed = 2;
    else if ( (strncasecmp(OptionText, "dist", 4) == 0)
               || (strncasecmp(OptionText, "dg", 2) == 0) ) {
      useDistGeom = true;
    }
  }

  // Give some limits so we can use switch statements
  if (speed < 1)
    speed = 1;
  else if (speed > 5)
    speed = 5;

  // This is done for all speed levels (i.e., create the structure)
  OBBuilder builder;
  bool attemptBuild = !useDistGeom;
  if (attemptBuild && !builder.Build(*pmol) ) {
    std::cerr << "Warning: Stereochemistry is wrong, using the distance geometry method instead" << std::endl;
    useDistGeom = true;
  }

#ifdef HAVE_EIGEN
  OBDistanceGeometry dg;
  if (useDistGeom) {
    dg.Setup(*pmol, attemptBuild); // use the bond lengths and angles if we ran the builder
    dg.GetGeometry(*pmol); // ensured to have correct stereo
  }
#endif

  // rule-based builder worked
  pmol->SetDimension(3);
  pmol->AddHydrogens(false, false); // Add some hydrogens before running MMFF

  if (speed == 5)
    return true; // done

  // All other speed levels do some FF cleanup
  // Try MMFF94 first and UFF if that doesn't work
  OBForceField* pFF = OBForceField::FindForceField("MMFF94");
  if (!pFF)
    return true;
  if (!pFF->Setup(*pmol)) {
    pFF = OBForceField::FindForceField("UFF");
    if (!pFF || !pFF->Setup(*pmol)) return true; // can't use either MMFF94 or UFF
  }

  // Since we only want a rough geometry, use distance cutoffs for VDW, Electrostatics
  pFF->EnableCutOff(true);
  pFF->SetVDWCutOff(10.0);
  pFF->SetElectrostaticCutOff(20.0);
  pFF->SetUpdateFrequency(10); // update non-bonded distances infrequently

  // How many cleanup cycles?
  int iterations = 250;
  switch (speed) {
  case 1:
    iterations = 500;
    break;
  case 2:
    iterations = 250;
    break;
  case 3:
  case 4:
  default:
    iterations = 100;
  }

  // Initial cleanup for every level
  pFF->ConjugateGradients(iterations, 1.0e-4);

  if (speed == 4) {
    pFF->UpdateCoordinates(*pmol);
    return true; // no conformer searching
  }

  switch(speed) {
  case 1:
    pFF->WeightedRotorSearch(250, 10); // maybe based on # of rotatable bonds?
    break;
  case 2:
    pFF->FastRotorSearch(true); // permute central rotors
    break;
  case 3:
  default:
    pFF->FastRotorSearch(false); // only one permutation
  }

  // Final cleanup and copy the new coordinates back
  pFF->ConjugateGradients(iterations, 1.0e-6);
  pFF->UpdateCoordinates(*pmol);

  return true;
}
}//namespace
