/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

namespace OpenBabel {

bool WriteGaussianCart(ostream &ofs,OEMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  
  ofs << "%cmem=20000000" << endl << '\045';
  ofs << "#Put Keywords Here" << endl << endl;
  ttab.SetFromType("INT"); ttab.SetToType("XYZ");

  OEAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%-3s0      x%-5d     y%-5d     z%-5d ",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    i, i, i);
    ofs << buffer << endl;
  }

  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"x%-4d %10.5f",
	    i,
	    atom->GetX());
    ofs << buffer << endl;
  }
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"y%-4d %10.5f",
	    i,
	    atom->GetY());
    ofs << buffer << endl;
  }
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"z%-4d %10.5f",
	    i,
	    atom->GetZ());
    ofs << buffer << endl;
  }

  return(true);
}

}
