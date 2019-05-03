/**********************************************************************
invalid-smiles.cpp - Test SMILES pattern parsing (rejecting invalid patterns)

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

Some portions Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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

#include <cstdlib>
#include <fstream>
#include <openbabel/obutil.h>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string itestdatadir = TESTDATADIR;
  string iinvalid_file = itestdatadir + "invalid-smiles.txt";
  string irandom1_file = itestdatadir + "random";
  string irandom2_file = itestdatadir + "random2";
  string irandom3_file = itestdatadir + "random3";
#else
  string iinvalid_file = "files/invalid-smiles.txt";
  string irandom1_file = "files/random";
  string irandom2_file = "files/random2";
  string irandom3_file = "files/random3";
#endif

int invalidsmiles(int argc, char* argv[])
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

  if (argc != 1)
    {
      cout << "Usage: invalid-smiles" << endl;
      cout << "   Tests Open Babel SMILES parsing - rejecting invalid patterns."
           << endl;
      return 0;
    }

  cout << "# Testing invalid SMILES parsing..." << endl;

  // make sure to kill off all error reporting
  obErrorLog.StopLogging();

  std::ifstream mifs;
  if (!SafeOpen(mifs, iinvalid_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << iinvalid_file << endl;
      return -1; // test failed
    }

  unsigned int currentTest = 0;
  OBMol mol;
  OBConversion conv(&mifs, &cout);
  if (! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cout << "Bail out! SMILES format is not loaded" << endl;
      return -1;
    }

  while(mifs.good())
    {
      mol.Clear();
      if (conv.Read(&mol))
        cout << "not ok " << ++currentTest << " # invalid SMILES was parsed "
             << " but should have been rejected" << endl;
      else
        cout << "ok " << ++currentTest << " # invalid SMILES rejected" << endl;
    }

  mifs.close();
  mifs.clear();

  // A known bug in 2.1.1 -- we hang on random SMILES garbage. This will be
  // fixed in the next release. PR#1730132
//   // random file#1
//   if (!SafeOpen(mifs, random1_file.c_str()))
//     {
//       cout << "Bail out! Cannot read file " << random1_file << endl;
//       return -1; // test failed
//     }

//   mol.Clear();
//   if (conv.Read(&mol, &mifs))
//     cout << "not ok " << ++currentTest << " # random data was parsed "
//          << " but should have been rejected\n";
//   else
//     cout << "ok " << ++currentTest << " # random data 1\n";
  
//   mifs.close();
//   mifs.clear();

//   // random2
//   if (!SafeOpen(mifs, random2_file.c_str()))
//     {
//       cout << "Bail out! Cannot read file " << random2_file << endl;
//       return -1; // test failed
//     }

//   mol.Clear();
//   if (conv.Read(&mol, &mifs))
//     cout << "not ok " << ++currentTest << " # random data #2 was parsed "
//          << " but should have been rejected\n";
//   else
//     cout << "ok " << ++currentTest << " # random data 2\n";

//   mifs.close();
//   mifs.clear();

//   // random3
//   if (!SafeOpen(mifs, random3_file.c_str()))
//     {
//       cout << "Bail out! Cannot read file " << random3_file << endl;
//       return -1; // test failed
//     }

//   mol.Clear();
//   if (conv.Read(&mol, &mifs))
//     cout << "not ok " << ++currentTest << " # random data #3 was parsed "
//          << " but should have been rejected\n";
//   else
//     cout << "ok " << ++currentTest << " # random data 3\n";

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}
