 /**********************************************************************
 smartstest.cpp - Test SMARTS algorithms and atom typing.

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

void GenerateSmartsReference();

#ifdef TESTDATADIR
  string mtestdatadir = TESTDATADIR;
  string msmarts_file = mtestdatadir + "smartstest.txt";
  string mresults_file = mtestdatadir + "smartsresults.txt";
  string msmilestypes_file = mtestdatadir + "attype.00.smi";
#else
   string msmarts_file = "files/smartstest.txt";
   string mresults_file = "files/smartsresults.txt";
   string msmilestypes_file = "files/attype.00.smi";
#endif

int smartstest(int argc, char* argv[])
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

    if (choice==99)
    {
          GenerateSmartsReference();
          return 0;
    }
  
  cout << endl << "# Testing SMARTS...  \n";

  std::ifstream ifs;
  if (!SafeOpen(ifs, msmarts_file.c_str()))
    {
      cout << "Bail out! Cannot read " << msmarts_file << endl;
      return -1; // test failed
    }

  //read in the SMARTS test patterns
  char buffer[BUFF_SIZE];
  vector<OBSmartsPattern*> vsp;
  for (;ifs.getline(buffer,BUFF_SIZE);)
    {
      if (buffer[0] == '#') // skip comment line
        continue;

      OBSmartsPattern *sp = new OBSmartsPattern;

      if (sp->Init(buffer))
        vsp.push_back(sp);
      else
        delete sp;
    }
  ifs.close();

  std::ifstream rifs;
  if (!SafeOpen(rifs, mresults_file.c_str()))
    {
      cout << "Bail out! Cannot read in results file " << mresults_file << endl;
      return -1; // test failed
    }
  unsigned int npats;
  rifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%d %*s",&npats);

  //make sure the number of SMARTS patterns is the same as in the
  //reference data
  if (npats != vsp.size())
    {
      cout << "Bail out! Correct number of patterns not read in" <<
        "Read in " << vsp.size() << " expected " << npats << endl;
      return -1; // test failed
    }

  std::ifstream mifs;
  if (!SafeOpen(mifs, msmilestypes_file.c_str()))
    {
      cout << "Bail out! Cannot read atom types " << msmilestypes_file << endl;
      return -1; // test failed
    }

  unsigned int k;
  unsigned int res_line = 0;
  OBMol mol;
  vector<string> vs;
  vector<OBSmartsPattern*>::iterator i;
  vector<vector<int> > mlist;
  unsigned int currentMol = 0; // each molecule is a separate test
  bool molPassed = true;

  OBConversion conv(&mifs, &cout);
  if (! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cout << "Bail out! SMILES format is not loaded" << endl;
      return -1;
    }

  //read in molecules, match SMARTS, and compare results to reference data
  for (;mifs;)
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;

      currentMol++;
      molPassed = true;

      for (i = vsp.begin();i != vsp.end();++i)
        {
          if (!rifs.getline(buffer,BUFF_SIZE))
            {
              cout << "Bail out! Error reading reference data" << endl;
              return -1; // test failed
            }
          res_line++;

          tokenize(vs,buffer);
          (*i)->Match(mol);
          mlist = (*i)->GetMapList();
          if (mlist.size() != vs.size())
            {
              cout << "not ok " << currentMol
                   << " # number of matches different than reference\n";
              cout << "# Expected " << vs.size() << " matches, found "
                   << mlist.size() << "\n";
              cout << "# Error with molecule " << mol.GetTitle();
              cout << "#  on pattern " << (*i)->GetSMARTS() << "\n";
              if (mlist.size())
                cout << "# First match: atom #" << mlist[0][0] << "\n";
              molPassed = false;
              continue;
            }

          if (mlist.size())
            {
              for (k = 0;k < vs.size();k++)
                {
                  if (atoi(vs[k].c_str()) != mlist[k][0])
                    {
                      cout << "not ok " << currentMol
                           << "# matching atom numbers different than reference\n";
                      cout << "# Expected " << vs[k] << " but found "
                           << mlist[k][0] << "\n";
                      cout << "# Molecule: " << mol.GetTitle() << "\n";
                      cout << "# Pattern: " << (*i)->GetSMARTS() << "\n";
                      molPassed = false;
                      break;
                    }
                }
              if (k != vs.size())
                {
                  molPassed = false;
                  continue;
                }
            }
        }
      if (molPassed)
        cout << "ok " << currentMol << " # molecule passed tests\n";
    }

  // output the number of tests run
  cout << "1.." << currentMol << endl;

  // Passed Test
  return 0;
}

void GenerateSmartsReference()
{
  std::ifstream ifs;
  if (!SafeOpen(ifs,msmarts_file.c_str()))
    return;

  char buffer[BUFF_SIZE];
  vector<OBSmartsPattern*> vsp;
  for (;ifs.getline(buffer,BUFF_SIZE);)
    {
      if (buffer[0] == '#') // skip comment line
        continue;

      OBSmartsPattern *sp = new OBSmartsPattern;

      if (sp->Init(buffer))
        vsp.push_back(sp);
      else
        delete sp;
    }

  std::ofstream ofs;
  if (!SafeOpen(ofs, mresults_file.c_str()))
    return;

  ofs << vsp.size() << " patterns" << endl;
  std::ifstream mifs;
  if (!SafeOpen(mifs, msmilestypes_file.c_str()))
    return;

  vector<int> vm;
  vector<vector<int> > mlist;
  vector<vector<int> >::iterator j;
  vector<OBSmartsPattern*>::iterator i;
  OBMol mol;
  OBConversion conv(&mifs, &cout);

  if(! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cerr << "SMILES format is not loaded" << endl;
      return;
    }

  for (;mifs;)
    {
      mol.Clear();
      conv.Read(&mol);

      if (mol.Empty())
        continue;
      for (i = vsp.begin();i != vsp.end();++i)
        {
          (*i)->Match(mol);
          mlist = (*i)->GetMapList();
          for (j = mlist.begin();j != mlist.end();++j)
            {
              sprintf(buffer,"%3d",*(j->begin()));
              ofs << buffer;
            }
          ofs << endl;
        }
    }


  cerr << " SMARTS test results written successfully" << endl;
  return;
}

