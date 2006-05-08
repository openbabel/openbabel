/**********************************************************************
mol.cpp - Unit tests for Open Babel OBMol class

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

#include "babelconfig.h"
#include "mol.h"
#include "obconversion.h"

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
      cout << "Usage: mol" << endl;
      cout << " Unit tests for OBMol " << endl;
      return(-1);
    }

  cout << "# Unit tests for OBMol \n";

  // the number of tests for "prove"
  cout << "1..7\n";

  cout << "ok 1\n"; // for loading tests

  OBMol emptyMol, testMol1;
  cout << "ok 2\n"; // ctor works

  testMol1.ReserveAtoms(-1);
  testMol1.ReserveAtoms(0);
  testMol1.ReserveAtoms(2);
  cout << "ok 3\n";

  // atom component tests

  if (testMol1.NumAtoms() == 0) {
    cout << "ok 4\n";
  } else {
    cout << "not ok 4\n";
  }

  testMol1.NewAtom();
  if (testMol1.NumAtoms() == 1) {
    cout << "ok 5\n";
  } else {
    cout << "not ok 5\n";
  }

  testMol1.NewAtom();
  testMol1.AddBond(1, 2, 1);
  if (testMol1.NumBonds() == 1) {
    cout << "ok 6\n";
  } else {
    cout << "not ok 6\n";
  }

  testMol1.Clear();
  if (testMol1.NumAtoms() == 0) {
    cout << "ok 7\n";
  } else {
    cout << "not ok 7\n";
  }

  return(0);
}
