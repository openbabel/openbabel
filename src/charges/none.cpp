/**********************************************************************
none.cpp - A OBChargeModel to clear all partial charges

Copyright (C) 2011 by Geoffrey R. Hutchison
(Based on a suggestion from Scott Brozell)

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
#include <openbabel/obiter.h>
#include <openbabel/atom.h>

namespace OpenBabel
{

class NoCharges : public OBChargeModel
{
public:
  NoCharges(const char* ID) : OBChargeModel(ID, false){};
  const char* Description(){ return "Clear all partial charges"; }

  /// \return whether partial charges were successfully assigned to this molecule
  bool ComputeCharges(OBMol &mol);
};

/////////////////////////////////////////////////////////////////
NoCharges theNoCharges("none"); //Global instance

/////////////////////////////////////////////////////////////////

  bool NoCharges::ComputeCharges(OBMol &mol)
  {
    mol.SetPartialChargesPerceived();

    FOR_ATOMS_OF_MOL(atom, mol) {
      atom->SetPartialCharge(0.0);
    }

    OBChargeModel::FillChargeVectors(mol);

    return true;
  }

}//namespace
