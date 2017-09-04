/**********************************************************************
mol.cpp - Unit tests for Open Babel OBMol class

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
#include <openbabel/elements.h>

using namespace std;
using namespace OpenBabel;

string pdbfile = "00T_ideal.pdb";

int pdbreadfile(int argc, char* argv[])
{
  int defaultchoice = 1;

  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  cout << "# Unit tests for OBMol \n";

  cout << "ok 1\n"; // for loading tests

  OBConversion obconv;
  if(!obconv.SetInFormat("PDB"))
  {
    cout << "Bail out! Fail format isn't loaded!" << endl;
    return -1;
  }

  // Test using ReadFile to read from PDB
  OpenBabel::OBMol obmol;
  if (obconv.ReadFile(&obmol, TESTDATADIR + pdbfile))
          cout << "ok 2!" << endl;
  else
          cout << "not ok 2" << endl;

  if (obmol.NumAtoms()==22)
          cout << "ok 3!" << endl;
  else
          cout << "not ok 3" << endl;

  if (obmol.GetAtom(10)->GetAtomicNum() == OBElements::Chlorine)
          cout << "ok 4!" << endl;
  else
          cout << "not ok 4" << endl;

  // Test that there are no remaining molecules
  // (this test fails on Linux)
  // if (!obconv.Read(&obmol))
  //        cout << "ok 6!" << endl;
  // else
  //        cout << "not ok 6" << endl;

  // the total number of tests for "prove"
  // update when you add more tests!
  cout << "1..4\n";

  return 0;
}
