/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

void GenerateSmartsReference();
bool TestSmarts(void);
void GenerateRingReference();
bool TestRings(void);
bool TestMatrixAlgebra(void);
bool TestUnitCell(void);

using namespace std;
using namespace OpenBabel;



int main(int argc,char *argv[])
{
  if (argc != 1) {
    cout << "Usage: obtest" << endl;
    return 0;
  }

  //GenerateSmartsReference();
  //GenerateRingReference();

  bool alltests = true;

  cout << endl << "Testing Matrix Algebra..." << endl;
  if (!TestMatrixAlgebra()) {
    cout << "ERROR: ***MATRIX test failed***" << endl;
    alltests = false;
  }
  else
    cout << "MATRIX test passed" << endl;

  cout << endl << "Testing SMARTS..." << endl;
  if (!TestSmarts()) {
    cout << "ERROR: ***SMARTS test failed***" << endl;
    alltests = false;
  }
  else
    cout << "SMARTS test passed" << endl;
    
  cout << endl << "Testing RINGS..." << endl;
  if (!TestRings()) {
    cout << "ERROR: ***RING test failed***" << endl;
    alltests = false;
  }
  else
    cout << "RING test passed" << endl;

  cout << endl << "Testing datatype conversions ..." << endl;
  cout << " DATATYPE: unit cell vectors";
  if (!TestUnitCell()) {
    cout << endl << "ERROR: *** Unit Cell test failed***" << endl;
  }
  else
    cout << " test passed " << endl;

  if (alltests) {
    cout << endl << "All tests passed" << endl;
    return 0;
  }
  else {
    cout << endl << "Some test(s) failed" << endl;
    return -1;
  }
  return 0;
}
