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
#include <openbabel/internalcoord.h>

#include <stdio.h>
#include <iostream>

using namespace std;
using namespace OpenBabel;

int internalcoord(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
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

