/**********************************************************************
cifspacegrouptest.cpp - Unit tests for to check if space group is being handled
properly in .cif format.

Copyright (C) 2016 by Schrodinger Inc.
 
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

#include <string>
#include <algorithm>

using namespace std;
using namespace OpenBabel;

std::string static GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;
  return path;
}

void testSpaceGroupUniqueTransformations()
{
  // See  https://github.com/openbabel/openbabel/pull/260
  // also https://github.com/openbabel/openbabel/pull/255
  OBConversion conv;
  OBMol mol;
  conv.SetInFormat("cif");
  conv.ReadFile(&mol, GetFilename("test01.cif"));
  OBUnitCell* pUC = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
  const SpaceGroup* pSG = pUC->GetSpaceGroup();
  SpaceGroup* sg = new SpaceGroup(*pSG);
  pSG = SpaceGroup::Find(sg);
  OB_ASSERT( pSG != NULL );

  // Check also for errors and warnings
  string summary = obErrorLog.GetMessageSummary();
  OB_ASSERT( summary.find("error") == string::npos);
  OB_ASSERT( summary.find("warning") == string::npos);

  OB_ASSERT( pSG->GetId() == 64 );
}

void testSpaceGroupClean()
{
  // See https://github.com/openbabel/openbabel/pull/254
  OBConversion conv;
  OBMol mol;
  conv.SetInFormat("cif");
  conv.ReadFile(&mol, GetFilename("test02.cif"));
  OBUnitCell* pUC = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
  const SpaceGroup* pSG = pUC->GetSpaceGroup();
  SpaceGroup* sg = new SpaceGroup(*pSG);
  pSG = SpaceGroup::Find(sg);
  OB_ASSERT( pSG != NULL );

  // Check also for errors and warnings
  string summary = obErrorLog.GetMessageSummary();
  OB_ASSERT( summary.find("error") == string::npos);
  OB_ASSERT( summary.find("warning") == string::npos);

  OB_ASSERT( pSG->GetId() == 166 );
}

void testSpaceGroupTransformations()
{
  // See https://github.com/openbabel/openbabel/pull/254
  SpaceGroup group;
  vector<string> trans;
  vector<string> trans_exp;
  vector<string> trans_got;

  // Same transformation
  trans_exp.push_back("x,y,z");
  trans.push_back("x,y,z");
  trans.push_back(" x , y , z ");
  trans.push_back("x,y+1,z");
  trans.push_back("x,y-1,z-1");
  trans.push_back("x,y-1,z");
  // Same transformation
  trans_exp.push_back("x,-y,z");
  trans.push_back("x,-y,z");
  trans.push_back("x,1-y,z");
  trans.push_back("x,-y+1,z");
  // Same transformation
  trans_exp.push_back("x,-y,-z");
  trans.push_back("x,-y,-z");
  trans.push_back("x,1-y,1-z");
  trans.push_back("x,-y+1,-z+1");
  // Same transformation
  trans_exp.push_back("x,-y,-y-z");
  trans.push_back("x,-y,-z-y");
  trans.push_back("x,1-y,-y+1-z");
  trans.push_back("x,-y+1,-y-z+1");
  // Same transformation
  trans_exp.push_back("x,-y,1/6-y-z");
  trans.push_back("x,-y,-z-y+1/6");
  trans.push_back("x,-y,-z-y+7/6");
  trans.push_back("x,1-y,7/6-y+1-z");
  trans.push_back("x,-y+1,-y-z+1+1/6");
  trans.push_back("x,-y+1,-y+1/6-z+1");
  // Same transformation
  trans_exp.push_back("x,3/4-y+z,5/6-y-z");
  trans.push_back("x,3/4-y+z,-1/6-z-y");
  trans.push_back("x,-y+3/4+z,-z-y-1/6");
  trans.push_back("x,z-y+3/4,-z-y-7/6");
  trans.push_back("x,1+z+3/4-y,-7/6-y-z");
  trans.push_back("X , 3 / 4 - Y + 1 + z , - Y - Z - 1 / 6 ");
  trans.push_back("x,z-y+1+3/4,-y-1/6-z+1");

  vector<string>::const_iterator i, iend;
  iend = trans.end();
  for (i = trans.begin(); i != iend; ++i)
    group.AddTransform(i->c_str());

  // Loop over symmetry operators
  transform3dIterator ti;
  const transform3d *t = group.BeginTransform(ti);
  while(t){
    trans_got.push_back(t->DescribeAsString());
    //cout << t->DescribeAsString() << "\n";
    t = group.NextTransform(ti);
  }

  OB_ASSERT( trans_exp.size() == trans_got.size() );
  OB_ASSERT( equal(trans_exp.begin(), trans_exp.end(), trans_got.begin()) );
}

int cifspacegrouptest(int argc, char* argv[])
{
  int defaultchoice = 3;

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
    testSpaceGroupUniqueTransformations();
    break;
  case 2:
    testSpaceGroupClean();
    break;
  case 3:
    testSpaceGroupTransformations();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return(0);
}
