/**********************************************************************
Copyright (C) 2003 Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "version.h"
#include "data.h"

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
  io_type inFile1 = UNDEFINED, inFile2 = UNDEFINED;
  OBFileFormat fileFormat;

  if (argc != 3) {
    cout << "Usage: roundtest <file1> <file2>" << endl;
    return(-1);
  }

  if (extab.CanReadExtension(argv[1]))
    inFile1 = extab.FilenameToType(argv[1]);
  else
    {
      cerr << argv[0] << ": Cannot read file-1 format!" << endl;
      return(-1);
    }
  if (extab.CanReadExtension(argv[2]))
    inFile2 = extab.FilenameToType(argv[2]);
  else
    {
      cerr << argv[0] << ": Cannot read file-2 format!" << endl;
      return(-1);
    }

  // Finally, we can do some work!
  OBMol mol(inFile1, UNDEFINED);
  OBMol mol2(inFile2, UNDEFINED);
  ifstream inFileStream1(argv[1]);
  ifstream inFileStream2(argv[2]);
  OBAtom *atom1, *atom2;

  if (!inFileStream1)
    {
      cerr << argv[0] << ": Cannot read input file-1!" << endl;
      return(-1);
    }
  else if (!inFileStream2)
    {
      cerr << argv[0] << ": Cannot read input file-2!" << endl;
      return(-1);
    }

  fileFormat.ReadMolecule(inFileStream1, mol, argv[1]);
  fileFormat.ReadMolecule(inFileStream2, mol2, argv[2]);

  if (mol.NumAtoms() != mol2.NumAtoms())
    {
      cout << " ** ERROR ** Number of atoms differ: " << mol.NumAtoms()
	   << " and " << mol2.NumAtoms() << endl;
      return(-1);
    }

  for(int i = 1;i <= mol.NumAtoms(); i++)
  {
    atom1 = mol.GetAtom(i);
    atom2 = mol2.GetAtom(i);

    if (atom1->GetAtomicNum() != atom2->GetAtomicNum())
      {
	cout << " ** ERROR ** Elements for atom " << i << " differ: " <<
	     atom1->GetAtomicNum() << " and " << atom2->GetAtomicNum() << endl;
	return(-1);
      }
    if ((atom1->GetX()-atom2->GetX()>1e-4) || 
	(atom1->GetY()-atom2->GetY()>1e-4) ||
	(atom1->GetZ()-atom2->GetZ()>1e-4))
      {
      	cout << " ** ERROR ** Coordinates for atom " << i << " differ." << endl;
	return(-1);
      }
  }
  
  return(0);
}
