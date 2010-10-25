/**********************************************************************
ffuff.cpp - Test energy and gradients for UFF force field

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

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/obutil.h>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string results_file = testdatadir + "uffresults.txt";
  string molecules_file = testdatadir + "forcefield.sdf";
#else
  string results_file = "files/uffresults.txt";
  string molecules_file = "files/forcefield.sdf";
#endif

void GenerateEnergies();

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
          cout << "Usage: ffuff" << endl;
          cout << "   Tests Open Babel UFF force field implementation." << endl;
          return 0;
        }
      else
        {
          GenerateEnergies();
          return 0;
        }
    }

  cout << "# Testing UFF Force Field..." << endl;

  std::ifstream mifs;
  if (!SafeOpen(mifs, molecules_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << molecules_file << endl;
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

  if(! conv.SetInAndOutFormats("SDF","SDF"))
    {
      cout << "Bail out! SDF format is not loaded" << endl;
      return -1; // test failed
    }
    
  OBForceField* pFF = OBForceField::FindForceField("UFF");

  if (pFF == NULL) {
    cerr << "Bail out! Cannot load force field!" << endl;
    return -1; // test failed
  }

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
          return -1; // test failed
        }
        
      if (!pFF->Setup(mol)) {
        cout << "Bail out! could not setup force field on " << mol.GetTitle() << endl;
        return -1; // test failed
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

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}

void GenerateEnergies()
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

  OBForceField* pFF = OBForceField::FindForceField("UFF");

  if (pFF == NULL) {
    cerr << "Cannot load force field!" << endl;
    return;
  }

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

	cerr << " UFF force field energies written successfully" << endl;
  return;
}
