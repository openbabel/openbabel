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

using namespace std;

namespace OpenBabel {
	
bool WriteTinker(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;
  
  sprintf(buffer,"%6d %-20s   MM2 parameters",mol.NumAtoms(),mol.GetTitle());
  ofs << buffer << endl;

  ttab.SetFromType("INT");

  OBAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    str = atom->GetType();
    ttab.SetToType("MM2"); ttab.Translate(str1,str);
    sprintf(buffer,"%6d %2s  %12.6f%12.6f%12.6f %5d",
	    i,
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ(),
	    atoi((char*)str1.c_str()));
    ofs << buffer;
    
    for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j))
      {
	sprintf(buffer,"%6d", (bond->GetNbrAtom(atom))->GetIdx());
	ofs << buffer;
      }

    ofs << endl;
  }

  return(true);
}

}
