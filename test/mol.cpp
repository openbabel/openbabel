/**********************************************************************
mol.cpp - Unit tests for Open Babel OBMol class

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

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <cstdlib>

#include <stdio.h>
#include <iostream>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string ktestdatadir = TESTDATADIR;
  string kd2file = ktestdatadir + "test2d.xyz";
  string kd3file = ktestdatadir + "test3d.xyz";
#else
  string kd2file = "files/test2d.xyz";
  string kd3file = "files/test3d.xyz";
#endif

int mol(int argc, char* argv[])
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

  cout << "# Unit tests for OBMol \n";

  cout << "ok 1\n"; // for loading tests

  OBMol emptyMol, testMol1;
  cout << "ok 2\n"; // ctor works

  testMol1.ReserveAtoms(-1);
  testMol1.ReserveAtoms(0);
  testMol1.ReserveAtoms(2);
  cout << "ok 3\n";

  // atom component tests

  if (testMol1.NumAtoms() == 0) {
    cout << "ok 4\n";
  } else {
    cout << "not ok 4\n";
  }

  testMol1.NewAtom();
  if (testMol1.NumAtoms() == 1) {
    cout << "ok 5\n";
  } else {
    cout << "not ok 5\n";
  }

  testMol1.NewAtom();
  testMol1.AddBond(1, 2, 1);
  if (testMol1.NumBonds() == 1) {
    cout << "ok 6\n";
  } else {
    cout << "not ok 6\n";
  }

  testMol1.Clear();
  if (testMol1.NumAtoms() == 0) {
    cout << "ok 7\n";
  } else {
    cout << "not ok 7\n";
  }

  ifstream ifs1(kd3file.c_str());
  if (!ifs1)
    {
      cout << "Bail out! Cannot read input file!" << endl;
      return(-1);
    }
  OBConversion conv(&ifs1, &cout);
  OBFormat* pFormat;
  
  pFormat = conv.FindFormat("XYZ");
  if ( pFormat == NULL )
    {
      cout << "Bail out! Cannot read file format!" << endl;
      return(-1);
    }
  if (! conv.SetInAndOutFormats(pFormat, pFormat))
    {
      cout << "Bail out! File format isn't loaded" << endl;
      return (-1);
    }

  OBMol testMol2D, testMol3D;
  if (conv.Read(&testMol3D))
    cout << "ok 8\n";
  else
    cout << "not ok 8\n";
  testMol3D.Center();
  
  // test bond insertion (PR#1665649)
  OBMol doubleBondMol;
  OBAtom *a1, *a2;
  OBBond *b;
  doubleBondMol.BeginModify();
  a1 = doubleBondMol.NewAtom();
  a1->SetVector(0.0, 0.0, 0.0);
  a1->SetAtomicNum(6);
  a2 = doubleBondMol.NewAtom();
  a2->SetVector(1.6, 0.0, 0.0);
  a2->SetAtomicNum(6);
  b = doubleBondMol.NewBond();
  b->SetBegin(a1);
  b->SetEnd(a2);
  a1->AddBond(b);
  a2->AddBond(b);
  doubleBondMol.EndModify();
  cout << "ok 9" << endl;

  // test AddHydrogens
  OBMol testMolH;
  testMolH.BeginModify();
  OBAtom *testAtom = testMolH.NewAtom();
  testAtom->SetVector(0.5f, 0.5f, 0.5f);
  testAtom->SetAtomicNum(6);
  testAtom->SetImplicitHCount(4);
  testMolH.EndModify();
  testMolH.AddHydrogens();
  if (testMolH.NumAtoms() == 5) {
    cout << "ok 10" << endl;
  } else {
    cout << "not ok 10" << endl;
  }

  // test AddHydrogens (pr #1665519)
  OBMol testMolH2;
  OBAtom *testAtom2 = testMolH2.NewAtom();
  testAtom2->SetVector(0.5f, 0.5f, 0.5f);
  testAtom2->SetAtomicNum(6);
  testAtom2->SetImplicitHCount(4);
  testMolH2.AddHydrogens();
  if (testMolH2.NumAtoms() == 5) {
    cout << "ok 11" << endl;
  } else {
    cout << "not ok 11 # hydrogen additions" << endl;
  }
  
  // Attempt to write an empty InChI (PR#2864334)
  pFormat = conv.FindFormat("InChI");
  if ( pFormat != NULL && conv.SetOutFormat(pFormat))
    {
      if (conv.Write(&emptyMol))
        cout << "ok 12" << endl;
      else
        cout << "not ok 12 # failed empty InChI" << endl;
    }

  OBMol testMolFormula;
  string formula("C6");
  testMolFormula.SetFormula(formula);
  if ( testMolFormula.GetFormula() == formula ) {
     cout << "ok 13" << endl;
  } else {
    cout << "not ok 13 # SetFormula "<< endl;
  }
  // Reset the formula to test for a double delete error
  testMolFormula.SetFormula(formula);
  
  // Test molecular formulas with large atomic numbers
  OBMol testLgAtNo;
  testLgAtNo.BeginModify();
  OBAtom *lgAtom = testLgAtNo.NewAtom();
  lgAtom->SetAtomicNum(118);
  // Undefined atomic numbers should be ignored with an obWarning instead of segfault
  lgAtom = testLgAtNo.NewAtom();
  lgAtom->SetAtomicNum(200);
  lgAtom = testLgAtNo.NewAtom();
  lgAtom->SetAtomicNum(1);
  lgAtom->SetIsotope(2);
  testLgAtNo.EndModify();
  if ( testLgAtNo.GetFormula() == "DOg" ) {
    cout << "ok 14" << endl;
  } else {
    cout << "not ok 14" << endl;
  }
  

  double dihedral = CalcTorsionAngle(vector3(-1., -1.,  0.),
                                     vector3(-1.,  0.,  0.),
                                     vector3( 1.,  0.,  0.),
                                     vector3( 1.,  1.,  0.));

  double dihedral_error = fabs(dihedral) - 180.0;

  if (fabs(dihedral_error) < 0.001) {
      std::cout << "ok 15 " << dihedral_error << std::endl;
  } else {

      std::cout << "not ok 15 # CalcTorsionAngle " << dihedral << "!= 180.0" << std::endl;
  }

  cout << "1..15\n"; // total number of tests for Perl's "prove" tool
  return(0);
}
