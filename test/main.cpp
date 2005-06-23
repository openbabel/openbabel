/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 Geoffrey R. Hutchison
 
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

#include "mol.h"

bool TestUnitCell(void);

using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
    if (argc != 1)
    {
        cout << "Usage: obtest" << endl;
        cout << "   Tests Open Babel generic datatype conversions." << endl;
        return 0;
    }

    // Right now, only the unit cell conversion tests, so we fail if that fails
    // Otherwise, we should decide if there should be unique tests for unit-tests
    //  or keep one tests for data type conversions and fail if any break

    cout << endl << "Testing datatype conversions ..." << endl;
    cout << " DATATYPE: unit cell vectors";
    if (!TestUnitCell())
    {
        cout << endl << "ERROR: *** Unit Cell test failed***" << endl;
        return -1;
    }
    else
        cout << " test passed " << endl;

    // Need to test:
    // atomic spin (i.e. radicals)
    // stereochemistry
    // bond types
    // atom types
    // fractional coordinate conversion
    // internal -> cartesian conversion
    // cartesian -> internal conversion
    // residue information
    // isotopes
    // chain information
    // chain perception
    // bond type perception
    // aromatic perception
    // ring perception
    // chirality perception
    // compress / uncompress
    // NumRotors

    // Builder utils
    // H->Methyl
    // Add hydrogens
    // SetHybAndGeom
    // SetLength
    // SetTorsion (obrotate)
    // Center
    // ToIntertialFrame
    // Rotate
    // Kekulize
    // Delete hydrogens
    // Add polar hydrogens
    // Add pH hydrogens
    // strip salts
    // obfit

    return 0;
}
