/**********************************************************************
Copyright (C) 2019 by NextMove Software

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

#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <openbabel/kekulize.h>


OpenBabel::OBAtom *NMOBMolNewAtom(OpenBabel::OBMol *mol, unsigned int elem)
{
  OpenBabel::OBAtom *result = mol->NewAtom();
  result->SetAtomicNum(elem);
  return result;
}


OpenBabel::OBBond *NMOBMolNewBond(OpenBabel::OBMol *mol,
                                  OpenBabel::OBAtom *beg,
                                  OpenBabel::OBAtom *end,
                                  unsigned int order, bool arom)
{
  if (!mol->AddBond(beg->GetIdx(),end->GetIdx(),order))
    return (OpenBabel::OBBond*)0;
  OpenBabel::OBBond *bptr = mol->GetBond(mol->NumBonds()-1);
  if (arom)
    bptr->SetAromatic();
  return bptr;
}


void NMOBAtomSetAromatic(OpenBabel::OBAtom *atm, bool arom)
{
  OpenBabel::OBMol *mol = (OpenBabel::OBMol*)atm->GetParent();
  if (mol && !mol->HasAromaticPerceived())
    mol->SetAromaticPerceived();

  if (arom)
    atm->SetAromatic();
  else
    atm->UnsetAromatic();
}


bool NMOBSanitizeMol(OpenBabel::OBMol *mol)
{
  if (!OBKekulize(mol))
    return false;
  mol->SetAromaticPerceived(false);
  return true;
}

