/**********************************************************************
ffmmff94.cpp - Test energy and gradients for MMFF94 force field

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

Some portions Copyright (C) 2008 Geoffrey R. Hutchison
 
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

#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/obutil.h>

using namespace std;
using namespace OpenBabel;

#ifndef TESTDATADIR
#define TESTDATADIR="files/";
#endif

  unsigned int currentTest = 0;

void GenerateEnergies(string molecules_file, string results_file)
{
  std::ifstream ifs;
  if (!SafeOpen(ifs, molecules_file.c_str()))
    return;

  std::ofstream ofs;
  if (!SafeOpen(ofs, results_file.c_str()))
    return;

  OBMol mol;
  OBConversion conv(&ifs, &cout);
  char buffer[BUFF_SIZE];
  
  if(! conv.SetInAndOutFormats("SDF","SDF"))
    {
      cerr << "SDF format is not loaded" << endl;
      return;
    }

  OBForceField* pFF = OBForceField::FindForceField("MMFF94");
  OB_REQUIRE(pFF != NULL);

  pFF->SetLogFile(&cout);
  pFF->SetLogLevel(OBFF_LOGLVL_NONE);

  for (;ifs;)
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;

      if (!pFF->Setup(mol)) {
        cerr << "Could not setup force field on molecule: " << mol.GetTitle() << endl;
        return;
      }
      
      // Don't compute gradients
      sprintf(buffer, "%15.5f\n", pFF->Energy(false));
      ofs << buffer;
    }

	cerr << " MMFF94 force field energies written successfully" << endl;
  return;
}

void TestFile(string filename, string results_file)
{
  std::ifstream mifs;
  if (!SafeOpen(mifs, filename.c_str()))
    {
      cout << "Bail out! Cannot read file " << filename << endl;
      return;
    }

  std::ifstream rifs;
  if (!SafeOpen(rifs, results_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << results_file << endl;
      return;
    }

  char buffer[BUFF_SIZE];
  vector<string> vs;
  OBMol mol;
  OBConversion conv(&mifs, &cout);

  if(! conv.SetInAndOutFormats("SDF","SDF"))
    {
      cout << "Bail out! SDF format is not loaded" << endl;
      return;
    }
    
  OBForceField* pFF = OBForceField::FindForceField("MMFF94");
  OB_REQUIRE(pFF != NULL);

  pFF->SetLogFile(&cout);
  pFF->SetLogLevel(OBFF_LOGLVL_NONE);

  double energy;
  while(mifs)
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;
      if (!rifs.getline(buffer,BUFF_SIZE))
        {
          cout << "Bail out! error reading reference data" << endl;
          return;
        }
        
      if (!pFF->Setup(mol)) {
        cout << "Bail out! could not setup force field on " << mol.GetTitle() << endl;
        return;
      }

      // compare the calculated energy to our reference data
      energy = pFF->Energy(false);
      if ( fabs(atof(buffer) - energy ) > 1.0e-3)
        {
          cout << "not ok " << ++currentTest << " # calculated energy incorrect "
               << " for molecule " << mol.GetTitle() << "\n";
          cout << "# Expected " << buffer << " found " <<
            energy << "\n";
        }
      else
        cout << "ok " << ++currentTest << " # energy \n";

      // check that gradients validate too
      if (!pFF->ValidateGradients())
        {
          cout << "not ok " << ++currentTest << " # gradients do not validate "
               << " for molecule " << mol.GetTitle() << "\n";
        }
      else
        cout << "ok " << ++currentTest << " # gradients \n";
    }
} // end TestFile

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

  string testdatadir = TESTDATADIR;
  string molecules_file = testdatadir + "forcefield.sdf";
  string results_file = testdatadir + "mmff94results.txt";

  if (argc > 1)
    {
      if (strncmp(argv[1], "-g", 2))
        { // Get the filenames from the command-line
          molecules_file = argv[1];
          results_file = argv[2];
          TestFile( string(argv[1]), string(argv[2]));
          cout << "1.." << currentTest << endl;
          return 0;
        }
      else
        {
          if (argc > 3) {
            molecules_file = argv[2];
            results_file = argv[3];
          }
          GenerateEnergies(molecules_file, results_file);
          return 0;
        }
    }

  cout << "# Testing MMFF94 Force Field..." << endl;
  TestFile(testdatadir + "forcefield.sdf", testdatadir + "mmff94results.txt");
  TestFile(testdatadir + "more-mmff94.sdf", testdatadir + "more-mmff94results.txt"); // provided by Paolo Tosco

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}

