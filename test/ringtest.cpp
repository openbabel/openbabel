/**********************************************************************
ringtest.cpp - Test ring perception algorithms.

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

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obutil.h>

using namespace std;
using namespace OpenBabel;

void GenerateRingReference();

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string results_file = testdatadir + "ringresults.txt";
  string smilestypes_file = testdatadir + "attype.00.smi";
#else
  string results_file = "files/ringresults.txt";
  string smilestypes_file = "files/attype.00.smi";
#endif

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
          cout << "Usage: ringtest" << endl;
          cout << "   Tests Open Babel ring perception testing." << endl;
          return 0;
        }
      else
        {
          GenerateRingReference();
          return 0;
        }
    }

  cout << "# Testing ring perception..." << endl;

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

  unsigned int size;
  OBBond *bond;
  OBAtom *atom;
  int count;
  char buffer[BUFF_SIZE];
  vector<string> vs;
  vector<OBRing*> vr;
  vector<bool> vb;
  vector<int> vi;
  OBMol mol;
  vector<string>::iterator i;
  vector<OBEdgeBase*>::iterator j;
  vector<OBNodeBase*>::iterator k;
  vector<OBRing*>::iterator m;
  OBConversion conv(&mifs, &cout);
  unsigned int currentTest = 0;

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

      vb.clear();
      vb.resize(mol.NumBonds(),false);
      //check ring bonds
      tokenize(vs,buffer);
      for (i = vs.begin();i != vs.end();i++)
        vb[atoi(i->c_str())] = true;

      for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
        {
          if (vb[bond->GetIdx()] != bond->IsInRing())
            {
              cout << "not ok " << ++currentTest
                   << " # ring bond data different than reference\n";
              cout << "# Molecule: " << mol.GetTitle() << "\n";
            }
          else
            cout << "ok " << ++currentTest
                 << " # correct ring bond data\n";
        }

      vr = mol.GetSSSR();
      if (!rifs.getline(buffer,BUFF_SIZE))
        {
          cout << "Bail out! error reading reference data\n";
          return -1; // test failed
        }
      sscanf(buffer,"%d",&size);
      if (vr.size() != size) //check SSSR size
        {
          cout << "not ok " << ++currentTest 
               << " # SSSR size different than reference\n";
          cout << "# Molecule: " << mol.GetTitle() << "\n";
        }
      else
        cout << "ok " << ++currentTest
             << " # SSSR size matches reference\n";

      if (!rifs.getline(buffer,BUFF_SIZE))
        {
          cout << "Bail out! error reading reference data" << endl;
          return -1; // test failed
        }

      tokenize(vs,buffer);
      i = vs.begin();
      for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k))
        {
          if (i == vs.end())
            {
              cout << "not ok " << ++currentTest << " # error in SSSR count\n";
              cout << "# Molecule: " << mol.GetTitle() << "\n";
            }
          else
            cout << "ok " << ++currentTest << " # correct SSSR count\n";

          count = 0;
          for (m = vr.begin();m != vr.end();m++)
            if ((*m)->_pathset[atom->GetIdx()])
              count++;

          if (atoi(i->c_str()) != count)
            {
              cout << "not ok " << ++currentTest << "# ring membership test failed\n";
              cout << "# Molecule: " << mol.GetTitle() << "\n";
            }
          else
            cout << "ok " << ++currentTest << " # ring membership passed\n";

          i++;
        }


    }

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}

void GenerateRingReference()
{
  std::ifstream ifs;

  if (!SafeOpen(ifs,smilestypes_file.c_str()))
    return;

  std::ofstream ofs;
  if (!SafeOpen(ofs,results_file.c_str()))
    return;

  int count;
  OBAtom *atom;
  OBBond *bond;
  char buffer[BUFF_SIZE];
  vector<OBRing*> vr;
  vector<OBEdgeBase*>::iterator i;
  vector<OBNodeBase*>::iterator j;
  vector<OBRing*>::iterator k;
  OBMol mol;
  OBConversion conv(&ifs, &cout);

  if(! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cerr << " SMILES format is not loaded" << endl;
      return;
    }

  for (;ifs;)
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;

      //write out ring bonds
      for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
        if (bond->IsInRing())
          {
            sprintf(buffer,"%3d",bond->GetIdx());
            ofs << buffer;
          }
      ofs << endl;

      vr = mol.GetSSSR();
      //write the total number of rings
      ofs << vr.size() << endl;

      //write the number of rings that each atom is a member of
      for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
        {
          count = 0;
          for (k = vr.begin();k != vr.end();k++)
            if ((*k)->_pathset[atom->GetIdx()])
              count++;

          sprintf(buffer,"%3d",count);
          ofs << buffer;
        }
      ofs << endl;

    }

  cerr << " Ring perception test results written successfully" << endl;
  return;
}

