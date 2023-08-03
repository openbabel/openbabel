 /**********************************************************************
 cansmi.cpp - Test Canonical SMILES generation -- write and then read

 This file is part of the Open Babel project.
 For more information, see <http://openbabel.org/>

 Copyright (C) 2001-2007 Geoffrey R. Hutchison

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

#include <fstream>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <cstdlib>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string btestdatadir = TESTDATADIR;
  string bsmilestypes_file = btestdatadir + "nci.smi";
#else
   string bsmilestypes_file = "nci.smi";
#endif

int cansmi(int argc, char* argv[])
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

  unsigned int currentMol = 0;
  OBMol mol, mol2;
  string buffer;

  cout << endl << "# Testing Canonical SMILES Generation ...  \n";
  
  std::ifstream mifs(bsmilestypes_file.c_str());
  if (!mifs)
    {
      cout << "Bail out! Cannot read test data " << bsmilestypes_file << endl;
      return -1; // test failed
    }

  OBConversion conv(&mifs, &cout);
  if (! conv.SetInAndOutFormats("SMI","CAN"))
    {
      cout << "Bail out! SMILES format is not loaded" << endl;
      return -1;
    }

  switch(choice) {
  case 1:
    //read in molecules (as SMI), write as CANSMI, read in again
    while (getline(mifs, buffer))
      {
        mol.Clear();
        if (!conv.ReadString(&mol, buffer)) {
          cout << "not ok " << ++currentMol << " # SMILES read failed"
               << " buffer was " << buffer << "\n";
          continue;
        }
        if (mol.Empty())
          continue;

        buffer = conv.WriteString(&mol);

        mol2.Clear();
        if (!conv.ReadString(&mol2, buffer)) {
          cout << "not ok " << ++currentMol << " # SMARTS did not match"
               << " for molecule " << buffer << "\n";
          continue;
        }

        // Now make sure the molecules are roughly equivalent
        if (mol.NumHvyAtoms() == mol2.NumHvyAtoms())
          cout << "ok " << ++currentMol << " # number of heavy atoms match\n";
        else
          cout << "not ok " << ++currentMol << " # number of heavy atoms wrong"
               << " for molecule " << buffer << "\n";
      }
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  // Passed Test
  return 0;
}
