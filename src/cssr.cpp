/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

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

namespace OpenBabel
{

bool WriteCSSR(ostream &ofs,OBMol &mol)
{ 
  char buffer[BUFF_SIZE];

  if (!mol.HasData(obUnitCell))
    {
      sprintf(buffer,
	      " REFERENCE STRUCTURE = 00000   A,B,C =%8.3f%8.3f%8.3f",
	      1.0,1.0,1.0); 
      ofs << buffer << endl;
      sprintf(buffer,
	      "   ALPHA,BETA,GAMMA =%8.3f%8.3f%8.3f    SPGR =    P1"
	      , 90.0f, 90.0f, 90.0f);
      ofs << buffer << endl;
    }
  else
    {
      OBUnitCell *uc = (OBUnitCell*)mol.GetData(obUnitCell);
      sprintf(buffer,
	      " REFERENCE STRUCTURE = 00000   A,B,C =%8.3f%8.3f%8.3f",
	      uc->GetA(), uc->GetB(), uc->GetC()); 
      ofs << buffer << endl;
      sprintf(buffer,
	      "   ALPHA,BETA,GAMMA =%8.3f%8.3f%8.3f    SPGR =    P1",
	      uc->GetAlpha() , uc->GetBeta(), uc->GetGamma());
      ofs << buffer << endl;
    }

  sprintf(buffer,"%4d   1 %s",mol.NumAtoms(), mol.GetTitle());
  ofs << buffer << endl << endl;

  OBAtom *atom,*nbr;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j;
  vector<int> vtmp(106,0);
  int bonds;

  for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    //assign_pdb_number(pdb_types,atom->GetIdx());
    vtmp[atom->GetAtomicNum()]++;
    sprintf(buffer,"%4d%2s%-3d  %9.5f %9.5f %9.5f ",
	    atom->GetIdx(),
	    etab.GetSymbol(atom->GetAtomicNum()),
	    vtmp[atom->GetAtomicNum()],
	    atom->x(),
	    atom->y(),
	    atom->z());
    ofs << buffer;
    bonds = 0;
    for (nbr = atom->BeginNbrAtom(j); nbr; nbr = atom->NextNbrAtom(j))
      {
	if (bonds > 8) break;
	sprintf(buffer,"%4d",nbr->GetIdx());
	ofs << buffer;
	bonds++;
      }
    for (; bonds < 8; bonds ++)
      {
    	sprintf(buffer,"%4d",0);
    	ofs << buffer;
      }
    sprintf(buffer," %7.3f%4d", atom->GetPartialCharge(), 1);
    ofs << buffer << endl;
  }

  return(true);
}

}
