/**********************************************************************
aromatest.cpp - Test Open Babel aromaticity perception

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>
#include <openbabel/elements.h>

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenBabel;

// Molecules aromatized in error (in the course of porting from SMARTS patterns)
void NegativeTestCases(int &molCount, unsigned int &testCount)
{
  // Check that not every atom is aromatic (i.e. negative test cases)
  const char* smiles[] = { "c1ccc2[N+]=c3ccccc3=Nc2c1", // N radical found in eMolecules
                           "N1S[SH+]C=C1",
                           "S1C=[NH+]=[NH+]=C1",
                           "C1(N23)=CC=CC2=CC=CC3=CC=C1", // pyrido[2,1,6-de]quinolizine - no atom is Daylight aromatic
                           0 };
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");

  for (int i = 0; smiles[i]; ++i) {
    mol.Clear();
    conv.ReadString(&mol, smiles[i]);
    molCount++;
    bool found_non_aromatic = false;
    FOR_ATOMS_OF_MOL(atom, mol) {
      if (atom->GetAtomicNum() == OBElements::Hydrogen)
        continue;
      if (!atom->IsAromatic()) {
        found_non_aromatic = true;
        break;
      }
    }
    if (found_non_aromatic)
      cout << "ok " << ++testCount << "\n";
    else
      cout << "not ok " << ++testCount << " # all atoms are aromatic in molecule " << molCount << " "
      << mol.GetTitle() << "\n";
  }
}

int aromatest(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  cout << endl << "# Testing aromaticity perception...  " << endl;
 
  #ifdef TESTDATADIR
    string testdatadir(TESTDATADIR);
    string filename = testdatadir + "aromatics.smi";
  #else
    string filename = "files/aromatics.smi";
  #endif

  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  ifstream ifs(filename.c_str());
  if (!ifs)
    {
      cout << "Bail out! Cannot read input file!" << endl;
      return(-1);
    }
  
  OBConversion conv(&ifs, &cout);
  OBFormat* pFormat;
  
  pFormat = conv.FormatFromExt("aromatics.smi");
  if ( pFormat == NULL )
    {
      cout << "Bail out! Cannot read file format!" << endl;
      return(-1);
    }
  
  // Finally, we can do some work!
  OBMol mol;
  
  unsigned int testCount = 0;

  if (! conv.SetInAndOutFormats(pFormat, pFormat))
    {
      cout << "Bail out! File format isn't loaded" << endl;
      return (-1);
    }
  
  int molCount;

  switch(choice) {
  case 1:
    molCount = 0;
    while(ifs.peek() != EOF && ifs.good())
      {
        mol.Clear();
        conv.Read(&mol);
        molCount++;
        for (int N = 0; N < 2; ++N)
        {
          if (N == 0)
            mol.AddHydrogens();
          else
            mol.DeleteHydrogens();
          FOR_ATOMS_OF_MOL(atom, mol)
          {
            if (atom->GetAtomicNum() == OBElements::Hydrogen)
              continue;

            if (atom->IsAromatic())
              cout << "ok " << ++testCount << "\n";
            else
            {
              cout << "not ok " << ++testCount << " # atom isn't aromatic!\n";
              cout << "# atom idx " << atom->GetIdx()
                << " in molecule " << molCount << " "
                << mol.GetTitle() << "\n";
            }
          }
        }
      } // while reading molecules

    NegativeTestCases(molCount, testCount);

    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return(0);
}
