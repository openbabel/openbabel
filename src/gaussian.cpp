/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2004 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

using namespace std;

namespace OpenBabel {

// Gaussian output is only tested for Gaussian98. If
// you encounter problems or success with earlier version
// please report it. 2002-07-23, Michael Banck<mbanck@gmx.net>
bool WriteGaussianCart(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  int charge = 0;
  unsigned int multiplicity = 0;
  char buffer[BUFF_SIZE];
  
  ofs << "%cmem=20000000" << endl << '\045';
  ofs << "#Put Keywords Here, check Charge and Multiplicity" << endl << endl;
  ofs << "XX " << mol.GetTitle() << endl << endl;

  OBAtom *atom;
  string str,str1;
  sprintf(buffer,"  %d  %d", mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity());
  ofs << buffer << endl;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%-3s      %10.5f      %10.5f      %10.5f ",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX(), atom->GetY(), atom->GetZ());
    ofs << buffer << endl;
  }
  // file should end with a blank line
  ofs << endl;
  return(true);
}

}
