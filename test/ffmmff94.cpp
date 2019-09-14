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

  int currentTest = 0;

void GenerateEnergies(string molecules_file, string results_file, string method, double epsilon = 1.0)
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

  OBForceField* pFF = OBForceField::FindForceField(method);
  OB_REQUIRE(pFF != NULL);

  pFF->SetLogFile(&cout);
  pFF->SetLogLevel(OBFF_LOGLVL_NONE);
  pFF->SetDielectricConstant(epsilon);

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

void TestFile(string filename, string results_file, string method, double epsilon = 1.0)
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

  OBForceField* pFF = OBForceField::FindForceField(method);
  OB_REQUIRE(pFF != NULL);

  pFF->SetLogFile(&cout);
  pFF->SetLogLevel(OBFF_LOGLVL_NONE);
  pFF->SetDielectricConstant(epsilon);

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

int ffmmff94(int argc, char* argv[])
{
  int defaultchoice = 1;

  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  string testdatadir = TESTDATADIR;

  if (choice == 99)
    {
      GenerateEnergies(testdatadir + "forcefield.sdf", testdatadir + "mmff94results.txt", "MMFF94");
      GenerateEnergies(testdatadir + "more-mmff94.sdf", testdatadir + "more-mmff94results.txt", "MMFF94"); // provided by Paolo Tosco
      GenerateEnergies(testdatadir + "forcefield.sdf", testdatadir + "mmff94sresults.txt", "MMFF94s");
      GenerateEnergies(testdatadir + "more-mmff94.sdf", testdatadir + "more-mmff94sresults.txt", "MMFF94s"); // ditto
      GenerateEnergies(testdatadir + "forcefield.sdf", testdatadir + "mmff94e4results.txt", "MMFF94", 4.0);
      GenerateEnergies(testdatadir + "more-mmff94.sdf", testdatadir + "more-mmff94e4results.txt", "MMFF94", 4.0); // provided by Paolo Tosco

      return 0;
    }

  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  cout << "# Testing MMFF94 Force Field..." << endl;
  switch(choice) {
  case 1:
    TestFile(testdatadir + "forcefield.sdf", testdatadir + "mmff94results.txt", "MMFF94");
    break;
  case 2:
    TestFile(testdatadir + "more-mmff94.sdf", testdatadir + "more-mmff94results.txt", "MMFF94"); // provided by Paolo Tosco
    break;
  case 3:
    TestFile(testdatadir + "forcefield.sdf", testdatadir + "mmff94sresults.txt", "MMFF94s");
    break;
  case 4:
    TestFile(testdatadir + "more-mmff94.sdf", testdatadir + "more-mmff94sresults.txt", "MMFF94s"); // ditto
    break;
  case 5:
    TestFile(testdatadir + "forcefield.sdf", testdatadir + "mmff94e4results.txt", "MMFF94", 4.0);
    break;
  case 6:
    TestFile(testdatadir + "more-mmff94.sdf", testdatadir + "more-mmff94e4sresults.txt", "MMFF94", 4.0);
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  // Passed tests
  return 0;
}
