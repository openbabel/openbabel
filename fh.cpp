/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

namespace OpenBabel
{

bool WriteFenskeZmat(ostream &ofs,OEMol &mol)
{
  OEAtom *atom,*a,*b,*c;
  char type[10],buffer[BUFF_SIZE];
  vector<OENodeBase*>::iterator i;
  
  vector<OEInternalCoord*> vic;
  vic.push_back((OEInternalCoord*)NULL);
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    vic.push_back(new OEInternalCoord);

  CartesianToInternal(vic,mol);

  ofs << endl << mol.NumAtoms() << endl;
    
  float r,w,t;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  { 
    a = vic[atom->GetIdx()]->_a;
    b = vic[atom->GetIdx()]->_b;
    c = vic[atom->GetIdx()]->_c;
    r = vic[atom->GetIdx()]->_dst;
    w = vic[atom->GetIdx()]->_ang;
    t = vic[atom->GetIdx()]->_tor;
    strcpy(type,etab.GetSymbol(atom->GetAtomicNum()));

    if (atom->GetIdx() == 1)
    {
      sprintf(buffer,"%-2s  1",type); 
      ofs << buffer << endl;
      continue;
    }

    if (atom->GetIdx() == 2)
    {
      sprintf(buffer,"%-2s%3d%6.3f",
	      type,
	      a->GetIdx(),r);
      ofs << buffer << endl;
      continue;
    }

    if (atom->GetIdx() == 3)
    {
      sprintf(buffer,"%-2s%3d%6.3f%3d%8.3f",
	      type,
	      a->GetIdx(),r, b->GetIdx(),w);
      ofs << buffer << endl;
      continue;
    }

    if (t < 0) t += 360;
    sprintf(buffer,"%-2s%3d%6.3f%3d%8.3f%3d%6.1f",
	    type,
	    a->GetIdx(),r,b->GetIdx(),w,c->GetIdx(),t);
    ofs << buffer << endl;
  }

  return(true);
}

}
