/**********************************************************************
carspacegrouptest.cpp - Unit tests for to check if space group is being handled
properly in .car format.

Copyright (C) 2013 by Schrodinger Inc.
 
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
#include <openbabel/generic.h>
#include <openbabel/obconversion.h>

#include <string>

using namespace std;
using namespace OpenBabel;

std::string static GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;
  return path;
}

void testSpaceGroupWithSpace()
{
  OBConversion conv;
  OBMol mol;
  conv.SetInFormat("car");
  conv.ReadFile(&mol, GetFilename("test3.car"));
  OBUnitCell* pUC = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
  OB_ASSERT( pUC->GetSpaceGroupName() == "P 4" );
}


void testSpaceGroupWithoutParentheses()
{
  OBConversion conv;
  OBMol mol;
  conv.SetInFormat("car");
  conv.ReadFile(&mol, GetFilename("test2.car"));
  OBUnitCell* pUC = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
  OB_ASSERT( pUC->GetSpaceGroupName() == "P4" );
}

void testSpaceGroupWithParentheses()
{
  OBConversion conv;
  OBMol mol;
  conv.SetInFormat("car");
  conv.ReadFile(&mol, GetFilename("test1.car"));
  OBUnitCell* pUC = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
  OB_ASSERT( pUC->GetSpaceGroupName() == "P4" );
  
}

void testDefaultSpaceGroup()
{
  OBConversion conv;
  OBMol mol;
  conv.SetInFormat("car");
  conv.ReadFile(&mol, GetFilename("monoclinic.car"));
  OBUnitCell* pUC = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
  OB_ASSERT( pUC->GetSpaceGroupName() == "" );
}

int carspacegrouptest(int argc, char* argv[])
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
    testDefaultSpaceGroup();
    break;
  case 2:
    testSpaceGroupWithParentheses();
    break;
  case 3:
    testSpaceGroupWithoutParentheses();
    break;
  case 4:
    testSpaceGroupWithSpace();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return(0);
}
