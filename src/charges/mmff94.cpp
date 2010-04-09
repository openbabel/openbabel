/**********************************************************************
mmff94.cpp - A OBPartialCharge model to assign MMFF94 partial charges

Copyright (C) 2010 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/partialcharge.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

namespace OpenBabel
{

class MMFF94Charges : public OBPartialCharge
{
public:
  MMFF94Charges(const char* ID) : OBPartialCharge(ID, false){};
  const char* Description(){ return "   Assign MMFF94 partial charges"; }

  /// \return whether partial charges were successfully assigned to this molecule
  bool AssignPartialCharges(OBMol &mol);
};

/////////////////////////////////////////////////////////////////
MMFF94Charges theMMFF94Charges("mmff94"); //Global instance

/////////////////////////////////////////////////////////////////

  bool MMFF94Charges::AssignPartialCharges(OBMol &mol)
  {
    mol.SetPartialChargesPerceived();

    // Annotate that partial charges come from MMFF94
    OBPairData *dp = new OBPairData;
    dp->SetAttribute("PartialCharges");
    dp->SetValue("MMFF94");
    dp->SetOrigin(perceived);
    mol.SetData(dp);

    OBForceField* pFF = OBForceField::FindForceField("MMFF94");
    if (!pFF || !pFF->Setup(mol))
      return false;

    pFF->GetPartialCharges(mol);
    FOR_ATOMS_OF_MOL(atom, mol) {
      OBPairData *chg = (OpenBabel::OBPairData*) atom->GetData("FFPartialCharge");
      if (chg)
        atom->SetPartialCharge(atof(chg->GetValue().c_str()));
    }
  }

}//namespace
