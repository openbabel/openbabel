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

bool TestUnitCell(void);

using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
  if (argc != 1) {
    cout << "Usage: obtest" << endl;
    cout << "   Tests Open Babel generic datatype conversions." << endl;
    return 0;
  }

  // Right now, only the unit cell conversion tests, so we fail if that fails
  // Otherwise, we should decide if there should be unique tests for unit-tests
  //  or keep one tests for data type conversions and fail if any break

  cout << endl << "Testing datatype conversions ..." << endl;
  cout << " DATATYPE: unit cell vectors";
  if (!TestUnitCell()) {
    cout << endl << "ERROR: *** Unit Cell test failed***" << endl;
    return -1;
  }
  else
    cout << " test passed " << endl;

  return 0;
}
