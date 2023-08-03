 /**********************************************************************
 smilesmatch.cpp - Test SMARTS matching (i.e., SMILES matching themselves)

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
#include <openbabel/obconversion.h>
#include <openbabel/obutil.h>
#include <openbabel/parsmart.h>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string otestdatadir = TESTDATADIR;
  string osmilestypes_file = otestdatadir + "nci.smi";
#else
   string osmilestypes_file = "nci.smi";
#endif

int smilesmatch(int argc, char* argv[])
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

  
  cout << endl << "# Testing SMILES self-matching using SMARTS...  \n";
  
  std::ifstream mifs;
  if (!SafeOpen(mifs, osmilestypes_file.c_str()))
    {
      cout << "Bail out! Cannot read test data " << osmilestypes_file << endl;
      return -1; // test failed
    }

  OBConversion conv(&mifs, &cout);
  if (! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cout << "Bail out! SMILES format is not loaded" << endl;
      return -1;
    }

  unsigned int currentMol = 0;
  OBSmartsPattern smarts;
  OBMol mol;
  string buffer;

  //read in molecules and see if their SMARTS matches themselves
  while (getline(mifs, buffer))
    {
      mol.Clear();
      conv.ReadString(&mol, buffer);
      if (mol.Empty())
        continue;

      // trim off any title, etc.
      string::size_type pos = buffer.find_first_of(" \t\n\r");
      if (pos != string::npos)
        buffer.erase(pos);

      pos = buffer.find_first_of('.');
      if (pos != string::npos)
        continue;

      smarts.Init(buffer);
      if (smarts.Match(mol))
        cout << "ok " << ++currentMol << " # SMARTS matched the"
             << " SMILES molecule\n";
      else
        cout << "not ok " << ++currentMol << " # SMARTS did not match"
             << " for molecule " << buffer << "\n";
    }

  // output the number of tests run
  cout << "1.." << currentMol << endl;

  // Passed Test
  return 0;
}
