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

namespace OpenEye
{

bool WriteCSSR(ostream &ofs,OEMol &mol)
{ 
  char buffer[BUFF_SIZE];

  sprintf(buffer,
	  " REFERENCE STRUCTURE = 00000   A,B,C =  %6.3f  %6.3f  %6.3f",
	  1.0,1.0,1.0); 
  ofs << buffer << endl;
  sprintf(buffer,
	  "   ALPHA,BETA,GAMMA =  90.000  90.000  90.000    SPGR =    P1");
  ofs << buffer << endl;
  sprintf(buffer,"%4d\n",mol.NumAtoms());
  ofs << buffer << endl;

  OEAtom *atom,*nbr;
  vector<OEAtom*>::iterator i;
  vector<OEBond*>::iterator j;
  vector<int> vtmp(106,0);

  for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    //assign_pdb_number(pdb_types,atom->GetIdx());
    vtmp[atom->GetAtomicNum()]++;
    sprintf(buffer," %3d%2s%-3d  %8.4f  %8.4f  %8.4f ",
	    atom->GetIdx(),
	    etab.GetSymbol(atom->GetAtomicNum()),
	    vtmp[atom->GetAtomicNum()],
	    atom->x(),
	    atom->y(),
	    atom->z());
    ofs << buffer;
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      {
	sprintf(buffer,"%4d",nbr->GetIdx());
	ofs << buffer;
      }
    ofs << endl;
  }

  return(true);
}

}
