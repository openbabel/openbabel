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

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string smilestypes_file = testdatadir + "nci.smi";
#else
   string smilestypes_file = "nci.smi";
#endif

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  if (argc != 1)
    {
      cout << "Usage: cansmih\n";
      cout << "   Tests Open Babel canonical SMILES generation." << endl;
      return 0;
    }
  
  cout << endl << "# Testing Canonical SMILES Generation ...  \n";
  
  std::ifstream mifs(smilestypes_file.c_str());
  if (!mifs)
    {
      cout << "Bail out! Cannot read test data " << smilestypes_file << endl;
      return -1; // test failed
    }

  OBConversion conv(&mifs, &cout);
  if (! conv.SetInAndOutFormats("SMI","CAN"))
    {
      cout << "Bail out! SMILES format is not loaded" << endl;
      return -1;
    }

  unsigned int currentMol = 0;
  OBMol mol, mol2;
  string buffer;

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

  // output the number of tests run
  cout << "1.." << currentMol << endl;

  // Passed Test
  return 0;
}
