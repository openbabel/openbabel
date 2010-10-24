/**********************************************************************
bond.cpp - Unit tests for Open Babel OBBond class

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
      cout << "Usage: bond" << endl;
      cout << " Unit tests for OBBond " << endl;
      return(-1);
    }

  cout << "# Unit tests for OBBond \n";

  // the number of tests for "prove"
  cout << "1..3\n";

  cout << "ok 1\n"; // for loading tests

  // OBBond isolation tests (no connection to residue, molecule...)

  OBAtom emptyAtom, begin1, end1;
  OBBond emptyBond, bond1;
  cout << "ok 2\n"; // constructors work OK

  bond1.SetBegin(&begin1);
  bond1.SetEnd(&end1);
  cout << "ok 3\n";

  return(0);
}

