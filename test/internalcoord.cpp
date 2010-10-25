/**********************************************************************
internalcoord.cpp - Unit tests for Open Babel OBInternalCoord class

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
      cout << "Usage: internalcoord " << endl;
      cout << " Unit tests for internal coordinate conversion " << endl;
      return(-1);
    }

  cout << "# Unit tests for z-matrix conversion \n";

  // the number of tests for "prove"
  cout << "1..2\n";

  cout << "ok 1\n"; // for loading tests

  // OBInternalCoord isolation tests

  OBInternalCoord emptyIC, testIC1, testIC2, testIC3, testIC4;
  cout << "ok 2\n"; // ctor works OK

  return(0);
}

