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

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obutil.h>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string results_file = testdatadir + "formalchargeresults.txt";
  string smilestypes_file = testdatadir + "attype.00.smi";
#else
  string results_file = "files/formalchargeresults.txt";
  string smilestypes_file = "files/attype.00.smi";
#endif

void GenerateFormalChargeReference();

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  if (argc != 1)
    {
      if (strncmp(argv[1], "-g", 2))
        {
          cout << "Usage: formalcharge" << endl;
          cout << "   Tests Open Babel molecular formal charge perception." 
               << endl;
          return 0;
        }
      else
        {
          GenerateFormalChargeReference();
          return 0;
        }
    }

  cout << "# Testing molecular formal charges..." << endl;

  std::ifstream mifs;
  if (!SafeOpen(mifs, smilestypes_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << smilestypes_file << endl;
      return -1; // test failed
    }

  std::ifstream rifs;
  if (!SafeOpen(rifs, results_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << results_file << endl;
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
  if (!SafeOpen(ifs, smilestypes_file.c_str()))
    return;

  std::ofstream ofs;
  if (!SafeOpen(ofs, results_file.c_str()))
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

