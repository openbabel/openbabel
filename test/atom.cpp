/**********************************************************************
atom.cpp - Unit tests for Open Babel OBAtom class

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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
      cout << "Usage: atom" << endl;
      cout << " Unit tests for OBAtom " << endl;
      return(-1);
    }

  cout << "# Unit tests for OBAtom \n";

  // the number of tests for "prove"
  cout << "1..6\n";

  cout << "ok 1\n"; // for loading tests

  OBAtom emptyAtom, testAtom1, testAtom2;
  cout << "ok 2\n"; // constructor works OK

  // OBAtom isolation tests (no connection to residue, bond, molecule...)

  // Beat up on SetIdx
  testAtom1.SetIdx(0);
  cout <<  testAtom1.GetIdx() << "\n";
  testAtom1.SetIdx(-1);
  cout <<  testAtom1.GetIdx() << "\n";
  testAtom1.SetIdx(1);
  cout <<  "ok 3\n";

  // Beat up on atomic number
  testAtom1.SetAtomicNum(0);
  cout <<  testAtom1.GetAtomicNum() << "\n";
  testAtom1.SetAtomicNum(-1);
  cout <<  testAtom1.GetAtomicNum() << "\n";
  testAtom1.SetAtomicNum(200);
  cout <<  testAtom1.GetAtomicNum() << "\n";
  testAtom1.SetAtomicNum(300);
  cout <<  testAtom1.GetAtomicNum() << "\n";
  testAtom1.SetAtomicNum(1);
  cout <<  "ok 4\n";

  // PR#1373650
  double *coordPtr = testAtom1.GetCoordinate();
  cout <<  "ok 5\n";

  testAtom1.SetCoordPtr(&coordPtr);
  cout <<  "ok 6\n";

  return(0);
}

