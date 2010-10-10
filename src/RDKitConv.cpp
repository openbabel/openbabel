/**********************************************************************
RDKitConv.cpp - Convert OpenBabel OBMol to RGKit molecules

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
#include <RDKitConv.h>

using OpenBabel::OBMolAtomIter;
using OpenBabel::OBMolBondIter;

RDKit::RWMol OBMolToRWMol(OpenBabel::OBMol* pOBMol);

RDKit::RWMol OBMolToRWMol(OpenBabel::OBMol* pOBMol)
{
  RDKit::RWMol RDMol;

  FOR_ATOMS_OF_MOL(a,pOBMol)
  {
    RDMol.addAtom();
    RDKit::Atom* pRDAtom = RDMol.getActiveAtom();
    pRDAtom->setAtomicNum(a->GetAtomicNum());
    pRDAtom->setFormalCharge(a->GetFormalCharge());
  }
  FOR_BONDS_OF_MOL(b,pOBMol)
  {
    //bond order >3 needs doing properly
    //assume RDKit atom indices start at 0
    RDMol.addBond(b->GetBeginAtomIdx()-1, b->GetEndAtomIdx()-1, (RDKit::Bond::BondType)b->GetBO());
  }
  std::string msg("RWMol made from ");
  if(pOBMol->GetTitle())
    msg += pOBMol->GetTitle();
  else
    msg += "OBMol";
  OpenBabel::obErrorLog.ThrowError(__FUNCTION__, msg, OpenBabel::obInfo);
  return RDMol;
}

