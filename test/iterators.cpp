 /**********************************************************************
 iterators.cpp - tests for iterators

 This file is part of the Open Babel project.
 For more information, see <http://openbabel.org/>

 Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
 Some portions Copyright (C) 2001-2005 Geoffrey R. Hutchison

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
#include <cstdlib>

#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>
#include <openbabel/obutil.h>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string jtestdatadir = TESTDATADIR;
  string jsmilestypes_file = jtestdatadir + "attype.00.smi";
#else
   string jsmilestypes_file = "attype.00.smi";
#endif

int iterators(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }


  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  cout << endl << "# Testing iterators...  \n";
  
  std::ifstream mifs;
  if (!SafeOpen(mifs, jsmilestypes_file.c_str()))
    {
      cout << "Bail out! Cannot read test data " << jsmilestypes_file << endl;
      return -1; // test failed
    }

  OBConversion conv(&mifs, &cout);
  if (! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cout << "Bail out! SMILES format is not loaded" << endl;
      return -1;
    }

  OBMol mol;

  unsigned int currentTest = 0;
  unsigned int counter = 0;

  // run through atom and bond iterators
  while(mifs.peek() != EOF && mifs.good())
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;

      counter = 0;
      FOR_ATOMS_OF_MOL(atom, mol)
        ++counter;
      if (counter != mol.NumAtoms())
        cout << "not ok " << ++currentTest 
             << " # atom iterator failed: expected " << mol.NumAtoms()
             << " but found " << counter << '\n';
      else
        cout << "ok " << ++currentTest << '\n';

      counter = 0;
      FOR_DFS_OF_MOL(atom, mol)
        ++counter;
      if (counter != mol.NumAtoms())
        cout << "not ok " << ++currentTest 
             << " # depth-first atom iterator failed: expected " 
             << mol.NumAtoms() << " but found " << counter << '\n';
      else
        cout << "ok " << ++currentTest << '\n';

      counter = 0;
      FOR_BFS_OF_MOL(atom, mol)
        ++counter;
      if (counter != mol.NumAtoms())
        cout << "not ok " << ++currentTest 
             << " # breadth-first atom iterator failed: expected " 
             << mol.NumAtoms() << " but found " << counter << '\n';
      else
        cout << "ok " << ++currentTest << '\n';

      counter = 0;
      FOR_BONDS_OF_MOL(bond, mol)
        ++counter;
      if (counter != mol.NumBonds())
        cout << "not ok " << ++currentTest 
             << " # bond iterator failed: expected " << mol.NumBonds()
             << " but found " << counter << '\n';
      else
        cout << "ok " << ++currentTest << '\n';

      // test ring iterators: PR#2815025
      counter = 0;
      OBRing *prevRing = 0;
      if (mol.GetSSSR().size() != 0) {
        FOR_RINGS_OF_MOL(ring, mol) {
          if (prevRing != 0 && &*ring == prevRing) {
            break; // OOPS, same ring in a row?!
          }
          prevRing = &*ring;
          ++counter;
        }
        if (counter != mol.GetSSSR().size()) // number of rings
          cout << "not ok " << ++currentTest
               << " # ring iterator failed: expected " << mol.GetSSSR().size()
               << " but found " << counter << '\n';
        else
          cout << "ok " << ++currentTest << '\n';
      } // if (have rings)
    }

  // output the number of tests run
  cout << "1.." << currentTest << endl;

  // Passed Test
  return 0;
}
