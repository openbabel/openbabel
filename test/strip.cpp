/**********************************************************************
strip.cpp - Test OBMol::StripSalts and friends

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

Copyright (C) 2009 Geoffrey R. Hutchison

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

namespace OpenBabel
{
  bool SafeOpen(std::ifstream &fs, const char *filename);
  bool SafeOpen(std::ofstream &fs, const char *filename);
}

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string results_file = testdatadir + "stripresults.txt";
  string smilestypes_file = testdatadir + "attype.00.smi";
#else
  string results_file = "files/stripresults.txt";
  string smilestypes_file = "files/attype.00.smi";
#endif

void GenerateReference();

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      if (strncmp(argv[1], "-g", 2))
        {
          cout << "Usage: strip" << endl;
          cout << "   Tests Open Babel stripping and splitting functions." << endl;
          return 0;
        }
      else
        {
          GenerateReference();
          return 0;
        }
    }

  cout << "# Testing molecular stripping/splitting..." << endl;

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
  unsigned int currentTest = 0;
  // double mass;

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
      if (vs.size() != 2)
        {
          cout << "Bail out! Reference data has incorrect format" << endl;
          return -1; // test failed
        }

      if (mol.StripSalts(2)) {
        if (atoi(vs[0].c_str()) != mol.NumAtoms() || atoi(vs[1].c_str()) != mol.NumBonds())
        {
          cout << "not ok " << ++currentTest << " # molecular strip incorrect "
               << " for molecule " << mol.GetTitle() << "\n";
        }
      }
      else
        cout << "ok " << ++currentTest << " # one component molecule\n";

    }

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}

void GenerateReference()
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

      // if we've stripped components, output the new number of atoms
      if (mol.StripSalts(2))
        ofs << mol.NumAtoms() << " " << mol.NumBonds() << endl;
      else
        ofs << "0 0" << endl; // doesn't matter, as long as it's 2 components
    }

	cerr << " Strip results written successfully" << endl;
  return;
}
