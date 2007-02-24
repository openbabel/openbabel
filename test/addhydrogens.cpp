/**********************************************************************
addhydrogen.cpp - Unit tests for OBMol::AddHydrogens

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <stdio.h>
#include <iostream>

using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      cout << "Usage: addhydrogens" << endl;
      cout << " Unit tests for OBMol::AddHydrogens " << endl;
      return(-1);
    }

  cout << "# Unit tests for OBMol::AddHydrogens \n";

  // the number of tests for "prove"
  cout << "1..9\n";

  cout << "ok 1\n"; // for loading tests

  OBMol obMol;
  int nbonds;
  OBConversion obConversion;
  obConversion.SetInFormat("smi");
  cout << "ok 2\n";

  obConversion.ReadString(&obMol, "CC");
  cout << "ok 3\n";

  if (obMol.NumAtoms() == 2) {
    cout << "ok 4\n";
  } else {
    cout << "not ok 4\n";
  }

  obMol.AddHydrogens();
  
  if (obMol.NumAtoms() == 8) {
    cout << "ok 5\n";
  } else {
    cout << "not ok 5\n";
  }

  if (obMol.NumBonds() == 7) {
    cout << "ok 6\n";
  } else {
    cout << "not ok 6\n";
  }

  nbonds = obMol.NumBonds();
  FOR_BONDS_OF_MOL (b, obMol)
    nbonds--;
    
  if (!nbonds) {
    cout << "ok 7\n";
  } else {
    cout << "not ok 7\n";
  }
 
  obMol.DeleteHydrogens();
  
  nbonds = obMol.NumBonds();
  FOR_BONDS_OF_MOL (b, obMol)
    nbonds--;
  
  if (!nbonds) {
    cout << "ok 8\n";
  } else {
    cout << "not ok 8\n";
  }

  OBAtom *extraAtom = obMol.NewAtom();
  extraAtom->SetAtomicNum(6);
  OBBond *extraBond = obMol.NewBond();
  extraBond->SetBegin(obMol.GetAtom(1));
  extraBond->SetEnd(extraAtom);
  
  nbonds = obMol.NumBonds();
  FOR_BONDS_OF_MOL (b, obMol)
    nbonds--;
  
  if (!nbonds) {
    cout << "ok 9\n";
  } else {
    cout << "not ok 9\n";
  }

  return(0);
}


