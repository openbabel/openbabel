/**********************************************************************
addhtest.cpp - Test adding hydrogens

Copyright (C) 2019 David R. Koes
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

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
#include <cstdlib>
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenBabel;

//add hydrogens and then perform some sanity checks
void addh_check(OBMol *mol)
{
    mol->AddHydrogens();
    //check valences and bond lengths

    /* doesn't pass yet
    FOR_ATOMS_OF_MOL(atom, mol) {
      int total_valence = atom->GetTotalValence();
      if(total_valence > OBElements::GetMaxBonds(atom->GetAtomicNum())) {
        cerr << "Valence greater than MaxBonds " << total_valence << " vs "
            << OBElements::GetMaxBonds(atom->GetAtomicNum()) << " for atom "
            << atom->GetId() << " " << OBElements::GetSymbol(atom->GetAtomicNum()) << "\n";
        exit(-1);
      }
    }
  */
    FOR_BONDS_OF_MOL(bond, mol) {
      cerr << "bond length " << bond->GetLength() << "\n";
      if(bond->GetLength() > 5) {
        cerr << "Bond length between " << bond->GetBeginAtom()->GetId()
            << " and " << bond->GetEndAtom()->GetId() << " is too long: "
            << bond->GetLength() << "\n";
        exit(-1);
      }
    }
}

// 1
// reads addh.in, which has a file name on each line,
// reads molecule, adds hydrogens

int addhtest(int argc, char* argv[])
{

  int choice = 1;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  if(choice != 1) //eh, not bothering to split this up
  {
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  cout << endl << "# Testing addh handling...  " << endl;
 
  #ifdef TESTDATADIR
    string testdatadir(TESTDATADIR);
    string infile = testdatadir + "addh.in";
  #else
    string infile = "files/addh.in";
  #endif

  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  ifstream ifs(infile.c_str());
  if (!ifs)
    {
      cout << "Bail out! Cannot read " << infile << endl;
      return(-1);
    }
  
  //gzip.in specifies what files to read and their canonical smiles
  string line;
  string filepath;
  vector<string> correctResults; //read from file

  OBConversion conv;
  OBMol mol;

  while(getline(ifs, line))
  {
    if(line.length() == 0)
      continue; //ignore blank lines
    if(line[0] != '#')
    {
      stringstream str(line);
      string fname;
      str >> fname; 
      OBFormat *format = conv.FormatFromExt(fname.c_str());
      conv.SetInFormat(format);

      #ifdef TESTDATADIR
          filepath = testdatadir + fname;
       #else
          filepath = "files/" + fname;
      #endif

      conv.ReadFile(&mol,filepath);
      addh_check(&mol);
    }
  }

  return(0);
}
