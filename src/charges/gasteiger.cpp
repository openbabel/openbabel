/**********************************************************************
gasteiger.cpp - A OBPartialCharge model to handle Gasteiger sigma charges

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
#include <openbabel/molchrg.h>

namespace OpenBabel
{

class GasteigerCharges : public OBPartialCharge
{
public:
  GasteigerCharges(const char* ID) : OBPartialCharge(ID, false){};
  const char* Description(){ return "Assign Gasteiger-Marsili sigma partial charges"; }

  /// \return whether partial charges were successfully assigned to this molecule
  bool AssignPartialCharges(OBMol &mol);
};

/////////////////////////////////////////////////////////////////
GasteigerCharges theGasteigerCharges("gasteiger"); //Global instance

/////////////////////////////////////////////////////////////////

  bool GasteigerCharges::AssignPartialCharges(OBMol &mol)
  {
    mol.SetPartialChargesPerceived();

    OBGastChrg gc;
    return gc.AssignPartialCharges(mol);
  }

}//namespace
