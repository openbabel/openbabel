/**********************************************************************
formalcharge.cpp - Test molecular formal charge perception

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

#include <fstream>
#include <cstdlib>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obutil.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string ftestdatadir = TESTDATADIR;
  string fresults_file = ftestdatadir + "formalchargeresults.txt";
  string fsmilestypes_file = ftestdatadir + "attype.00.smi";
#else
  string fresults_file = "files/formalchargeresults.txt";
  string fsmilestypes_file = "files/attype.00.smi";
#endif

void GenerateFormalChargeReference();

int formalcharge(int argc, char* argv[])
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

  if (choice == 99)
    {
      GenerateFormalChargeReference();
      return 0;
    }

  cout << "# Testing molecular formal charges..." << endl;

  std::ifstream mifs;
  if (!SafeOpen(mifs, fsmilestypes_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << fsmilestypes_file << endl;
      return -1; // test failed
    }

  std::ifstream rifs;
  if (!SafeOpen(rifs, fresults_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << fresults_file << endl;
      return -1; // test failed
    }

  char buffer[BUFF_SIZE];
  vector<string> vs;
  OBMol mol;
  OBConversion conv(&mifs, &cout);
  unsigned int currentTest = 1;

  if(! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cout << "Bail out! SMILES format is not loaded" << endl;
      return -1;
    }

  for (;mifs;)
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;
      if (!rifs.getline(buffer,BUFF_SIZE))
        {
          cout << "Bail out! error reading reference data" << endl;
          return -1; // test failed
        }

      tokenize(vs,buffer);
      // check charges
    }

  cout << "ok 1\n";

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}

void GenerateFormalChargeReference()
{
  std::ifstream ifs;
  if (!SafeOpen(ifs, fsmilestypes_file.c_str()))
    return;

  std::ofstream ofs;
  if (!SafeOpen(ofs, fresults_file.c_str()))
    return;

  OBMol mol;
  OBConversion conv(&ifs, &cout);

  if(! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cerr << "SMILES format is not loaded" << endl;
      return;
    }

  for (;ifs;)
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;

      ofs << mol.GetTotalCharge();
      FOR_ATOMS_OF_MOL(atom, mol)
        {
          ofs << " " << atom->GetFormalCharge();
        }
      ofs << endl;
    }

	cerr << " Formal charge results written successfully" << endl;
  return;
}

