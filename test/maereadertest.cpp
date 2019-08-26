/**********************************************************************
maereadertest.cpp - Unit tests for the Maestro format reader

Copyright (C) 2019 by Schrodinger Inc.

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

#include "obtest.h"
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/generic.h>

#include <string>
#include <algorithm>

using namespace std;
using namespace OpenBabel;

void testMaeReader()
{
  OBConversion conv;
  OBMol mol;
  conv.SetOutFormat("smi");

  conv.SetInFormat("mae");
  conv.ReadFile(&mol, OBTestUtil::GetFilename("maereader.mae"));
  string mae_smi = conv.WriteString(&mol);

  conv.SetInFormat("maegz");
  conv.ReadFile(&mol, OBTestUtil::GetFilename("maereader.maegz"));
  string maegz_smi = conv.WriteString(&mol);

  OB_ASSERT(mae_smi == maegz_smi);

  // Erase any potential newlines
  mae_smi.erase(remove(mae_smi.begin(), mae_smi.end(), '\n'), mae_smi.end());
  mae_smi.erase(remove(mae_smi.begin(), mae_smi.end(), '\r'), mae_smi.end());

  const string known_smi = "C([N+](=O)[O-])[N+](=O)[O-]\t2:Acids";
  OB_COMPARE(mae_smi, known_smi);
}


void testMaeWriter()
{
  OBConversion conv;

  OBMol mol;
  conv.SetInFormat("mae");
  conv.SetOutFormat("mae");
  conv.ReadFile(&mol, OBTestUtil::GetFilename("maereader.mae"));
  string mae_file_txt = conv.WriteString(&mol);

  conv.SetInFormat("sdf");
  conv.ReadFile(&mol, OBTestUtil::GetFilename("gaff.sdf"));
  string mae_txt = conv.WriteString(&mol);

  OBMol mae_mol;
  conv.SetInFormat("mae");
  conv.ReadString(&mae_mol, mae_txt);

  OB_COMPARE(mae_mol.NumAtoms(), 39);

  // Verify that reading from sequential strings works
  conv.ReadString(&mae_mol, mae_file_txt);
  OB_COMPARE(mae_mol.NumAtoms(), 9);


}

int maereadertest(int argc, char* argv[])
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

  switch(choice) {
  case 1:
    testMaeReader();
    break;
  case 2:
    testMaeWriter();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return(0);
}
