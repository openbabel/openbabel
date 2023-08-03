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
#include <openbabel/atom.h>
#include <cstdlib>

#include <stdio.h>
#include <iostream>
#include <openbabel/elements.h>

using namespace std;
using namespace OpenBabel;


int pdbreadfile(int argc, char* argv[])
{
  int defaultchoice = 1;

  int choice = defaultchoice;

  string pdbfile;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  switch(choice) {
  case 1:
    pdbfile = "00T_ideal.pdb";
    break;
  case 2:
    pdbfile = "00T_nonstandard.pdb";
    break;
  case 3:
    pdbfile = "00T_ideal_het.pdb";
    break;
  case 4:
    pdbfile = "00T_nonstandard_het.pdb";
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
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

  if (obmol.GetAtom(6)->GetAtomicNum() == OBElements::Nitrogen)
          cout << "ok 5!" << endl;
  else
          cout << "not ok 5" << endl;

  if (obmol.GetAtom(12)->GetAtomicNum() == OBElements::Hydrogen)
          cout << "ok 6!" << endl;
  else
          cout << "not ok 6" << endl;

  if (obmol.GetAtom(13)->GetAtomicNum() == OBElements::Hydrogen)
          cout << "ok 7!" << endl;
  else
          cout << "not ok 7" << endl;

  if (obmol.GetAtom(14)->GetAtomicNum() == OBElements::Hydrogen)
          cout << "ok 8!" << endl;
  else
          cout << "not ok 8" << endl;

  // Test that there are no remaining molecules
  // (this test fails on Linux)
  // if (!obconv.Read(&obmol))
  //        cout << "ok 6!" << endl;
  // else
  //        cout << "not ok 6" << endl;

  // the total number of tests for "prove"
  // update when you add more tests!
  cout << "1..6\n";

  return 0;
}
