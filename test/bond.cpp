/**********************************************************************
bond.cpp - Unit tests for Open Babel OBBond class

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

#include <stdio.h>
#include <iostream>

#include <iostream>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
OBMol * mol = new OBMol();
OBConversion obconversion;
obconversion.SetInFormat("smiles");
//obconversion.ReadString(mol, "C1(N=C(NC(C)C)N=C(N=1)OC)NC(C)C");
obconversion.ReadString(mol,
"N12[C@@H]([C@@H](NC([C@@H](c3ccsc3)C(=O)O)=O)C2=O)SC(C)(C)[C@@-]1C(=O)O");
cout << "formula "<< mol->GetFormula() << "\n";
cout << "Has3D: " << mol->Has3D() << "\n";
OBBuilder builder;
builder.Build(*mol);
cout << "Has3D: " << mol->Has3D() << "\n";
delete(mol);
return 0;
}