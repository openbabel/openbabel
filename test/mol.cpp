/**********************************************************************
mol.cpp - Unit tests for Open Babel OBMol class

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

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

#include <stdio.h>
#include <iostream>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string d2file = testdatadir + "test2d.xyz";
  string d3file = testdatadir + "test3d.xyz";
#else
  string d2file = "files/test2d.xyz";
  string d3file = "files/test3d.xyz";
#endif

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      cout << "Usage: mol" << endl;
      cout << " Unit tests for OBMol " << endl;
      return(-1);
    }

  cout << "# Unit tests for OBMol \n";

  // the number of tests for "prove"
  cout << "1..9\n";

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

  ifstream ifs1(d3file.c_str());
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
  //  testMol3D.Center();
  
  // test AddHydrogens
  OBMol testMolH;
  testMolH.BeginModify();
  OBAtom *testAtom = testMolH.NewAtom();
  testAtom->SetVector(0.5f, 0.5f, 0.5f);
  testAtom->SetAtomicNum(6);
  testMolH.EndModify();
  testMolH.AddHydrogens();
  if (testMolH.NumAtoms() == 5) {
    cout << "ok 9" << endl;;
  } else {
    cout << "not ok 9" << endl;;
  }
  
   // test AddHydrogens
  OBMol testMolH2;
  OBAtom *testAtom2 = testMolH2.NewAtom();
  testAtom2->SetVector(0.5f, 0.5f, 0.5f);
  testAtom2->SetAtomicNum(6);
  testMolH2.AddHydrogens();
  if (testMolH2.NumAtoms() == 5) {
    cout << "ok 10" << endl;;
  } else {
    cout << "not ok 10" << endl;;
  }
  
  return(0);
}
