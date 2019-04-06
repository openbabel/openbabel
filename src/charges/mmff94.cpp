/**********************************************************************
mmff94.cpp - A OBChargeModel to assign MMFF94 partial charges

Copyright (C) 2010 by Geoffrey R. Hutchison

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
#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/obiter.h>
#include <openbabel/generic.h>

#include <cstdlib>

namespace OpenBabel
{

class MMFF94Charges : public OBChargeModel
{
public:
  MMFF94Charges(const char* ID) : OBChargeModel(ID, false){};
  const char* Description(){ return "   Assign MMFF94 partial charges"; }

  /// \return whether partial charges were successfully assigned to this molecule
  bool ComputeCharges(OBMol &mol);

  double DipoleScalingFactor() { return 3.8558; } // fit from regression across MMFF94 validation set
};

/////////////////////////////////////////////////////////////////
MMFF94Charges theMMFF94Charges("mmff94"); //Global instance

/////////////////////////////////////////////////////////////////

  bool MMFF94Charges::ComputeCharges(OBMol &mol)
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
    m_partialCharges.clear();
    m_partialCharges.reserve(mol.NumAtoms());
    m_formalCharges.clear();
    m_formalCharges.reserve(mol.NumAtoms());
    FOR_ATOMS_OF_MOL(atom, mol) {
      OBPairData *chg = (OpenBabel::OBPairData*) atom->GetData("FFPartialCharge");
      if (chg)
        atom->SetPartialCharge(atof(chg->GetValue().c_str()));
      m_partialCharges.push_back(atom->GetPartialCharge());
      m_formalCharges.push_back(atom->GetFormalCharge());
    }
    return true;
  }

}//namespace
