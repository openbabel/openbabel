/**********************************************************************
formula.cpp - Test molecular formula, weight & exact mass

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

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string gtestdatadir = TESTDATADIR;
  string gresults_file = gtestdatadir + "formularesults.txt";
  string gsmilestypes_file = gtestdatadir + "attype.00.smi";
#else
  string gresults_file = "files/formularesults.txt";
  string gsmilestypes_file = "files/attype.00.smi";
#endif

void GenerateFormulaReference();

int formula(int argc, char* argv[])
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

  bool check = true;
  string filename;
  if (choice == 99)
    {
      GenerateFormulaReference();
      return 0;
    }

  cout << "# Testing molecular formulas..." << endl;

  std::ifstream mifs;
  char buffer[BUFF_SIZE];
  vector<string> vs;
  OBMol mol;
  OBConversion conv(&mifs, &cout);
  unsigned int currentTest = 0;

  if (check) {
    filename = gsmilestypes_file;
  }

  if (!SafeOpen(mifs, filename.c_str()))
    {
      cout << "Bail out! Cannot read file " << filename << endl;
      return -1; // test failed
    }

  OBFormat *format = conv.FormatFromExt(filename.c_str());

  std::ifstream rifs;
  if (check) {
    if (!SafeOpen(rifs, gresults_file.c_str()))
      {
        cout << "Bail out! Cannot read file " << gresults_file << endl;
        return -1; // test failed
      }
  }

  if(! conv.SetInFormat(format))
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

      if (!check) {
        // just give the molecular formula
        cout << "Formula: " << mol.GetFormula() << '\n';
        continue;
      }

      if (!rifs.getline(buffer,BUFF_SIZE))
        {
          cout << "Bail out! error reading reference data" << endl;
          return -1; // test failed
        }

      tokenize(vs,buffer);
      if (vs.size() != 3)
        {
          cout << "Bail out! Reference data has incorrect format" << endl;
          return -1; // test failed
        }

      if (vs[0] != mol.GetFormula())
        {
          cout << "not ok " << ++currentTest << " # molecular formula incorrect"
               << " for molecule " << mol.GetTitle() << "\n";
        }
      else
        cout << "ok " << ++currentTest << " # molecular formula\n";

      if ( fabs(atof(vs[1].c_str()) - mol.GetMolWt() ) > 1.0e-3)
        {
          cout << "not ok " << ++currentTest << " # molecular weight incorrect"
               << " for molecule " << mol.GetTitle() << "\n";
          cout << "# Expected " << atof(vs[1].c_str()) << " found " <<
            mol.GetMolWt() << "\n";
        }
      else
        cout << "ok " << ++currentTest << " # molecular weight\n";

      if ( fabs(atof(vs[2].c_str()) - mol.GetExactMass() ) > 1.0e-3)
        {
          cout << "not ok " << ++currentTest << " # exact mass incorrect"
               << " for molecule " << mol.GetTitle() << "\n";
          cout << "# Expected " << atof(vs[2].c_str()) << " found " <<
            mol.GetExactMass() << "\n";
        }
      else
        cout << "ok " << ++currentTest << " # molecular exact mass\n";


      // now after adding explicit hydrogens -- should be identical
      //  since we'll add hydrogens that were implicit before

      // PR#1485580
      mol.AddHydrogens();

      if (vs[0] != mol.GetFormula())
        {
          cout << "not ok " << ++currentTest << " # molecular formula incorrect"
               << " for hydrogen-added molecule " << mol.GetTitle() << "\n";
        }
      else
        cout << "ok " << ++currentTest << " # molecular hydrogen-added formula\n";

      if ( fabs(atof(vs[1].c_str()) - mol.GetMolWt() ) > 1.0e-3)
        {
          cout << "not ok " << ++currentTest << " # molecular weight incorrect"
               << " for hydrogen-added molecule " << mol.GetTitle() << "\n";
          cout << "# Expected " << atof(vs[1].c_str()) << " found " <<
            mol.GetMolWt() << "\n";
          cout << "# Difference " << fabs(atof(vs[1].c_str()) - mol.GetMolWt())
               << "\n";
        }
      else
        cout << "ok " << ++currentTest << " # molecule + hydrogens weight\n";

      if ( fabs(atof(vs[2].c_str()) - mol.GetExactMass() ) > 1.0e-3)
        {
          cout << "not ok " << ++currentTest << " # exact mass incorrect"
               << " for hydrogen-added molecule " << mol.GetTitle() << "\n";
          cout << "# Expected " << atof(vs[2].c_str()) << " found " <<
            mol.GetExactMass() << "\n";
          cout << "# Difference " << fabs(atof(vs[2].c_str()) - mol.GetExactMass())
               << "\n";
        }
      else
        cout << "ok " << ++currentTest << " # molecular exact mass"
             << " after hydrogen addition\n";

    }

  // ---------------------------------------------------------------------------
  // Regression test for issue #3008: molreport / GetFormula dropped an explicit
  // hydrogen that is bonded to no heavy atom (isolated, or bonded only to other
  // hydrogens). Such a hydrogen was counted by nobody, so it vanished from the
  // formula even though it remained in the atom list and the molecular weight.
  // ---------------------------------------------------------------------------
  {
    struct { const char *smiles; const char *formula; } issue3008[] = {
      { "[N][H][H]",            "H2N"     }, // stray H must not be dropped
      { "[H].[O][H]",           "H2O"     }, // isolated H plus an O-H fragment
      { "[He].[H].[H].[Co][C]", "CH2CoHe" }, // isolated H with a bonded fragment
      // Controls that must stay correct after the fix:
      { "[H][O][H]",            "H2O"     }, // both H bonded to a heavy atom
      { "[H].[O].[H]",          "H2O"     }, // no bonds: implicit-H path unused
      { "[H][H]",               "H2"      }  // no heavy atoms present
    };
    OBConversion smiconv;
    if (!smiconv.SetInFormat("smi"))
      cout << "not ok " << ++currentTest << " # SMILES format not loaded\n";
    else
      for (size_t i = 0; i < sizeof(issue3008) / sizeof(issue3008[0]); ++i)
        {
          OBMol cmol;
          smiconv.ReadString(&cmol, issue3008[i].smiles);
          string got = cmol.GetFormula();
          if (got != issue3008[i].formula)
            cout << "not ok " << ++currentTest << " # formula \"" << got
                 << "\" != \"" << issue3008[i].formula << "\" for "
                 << issue3008[i].smiles << " (issue #3008)\n";
          else
            cout << "ok " << ++currentTest << " # formula \"" << got
                 << "\" for " << issue3008[i].smiles << " (issue #3008)\n";
        }
  }

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}

void GenerateFormulaReference()
{
  std::ifstream ifs;
  if (!SafeOpen(ifs, gsmilestypes_file.c_str()))
    return;

  std::ofstream ofs;
  if (!SafeOpen(ofs, gresults_file.c_str()))
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

      //write out formula, molecular weight and exact mass
      ofs << mol.GetFormula() << " " << mol.GetMolWt() << " "
          << mol.GetExactMass() << endl;
    }

	cerr << " Molecular formula results written successfully" << endl;
  return;
}
