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

using namespace std;

namespace OpenBabel
{

bool ReadFeat(istream &ifs,OBMol &mol, const char *title)
{
  char buffer[BUFF_SIZE];
  int i,natoms;

  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%d",&natoms);

  mol.ReserveAtoms(natoms);
  mol.BeginModify();

  if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
  mol.SetTitle(buffer);

  float x,y,z;
  char type[20];
  OBAtom *atom;
  for (i = 0; i < natoms;i++)
  {
    if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
    sscanf(buffer,"%s %f %f %f",
	   type,
	   &x,
	   &y,
	   &z);
    CleanAtomType(type);
    atom = mol.NewAtom();
    atom->SetVector(x,y,z);
    atom->SetAtomicNum(etab.GetAtomicNum(type));
  }

  mol.EndModify();
  return(true);
}

bool WriteFeat(ostream &ofs,OBMol &mol)
{ 
  char buffer[BUFF_SIZE];
  
  ofs << mol.NumAtoms() << endl;
  ofs << mol.GetTitle() << endl;

  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    sprintf(buffer,"%-3s %8.5f  %8.5f  %8.5f ",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->x(),
	    atom->y(),
	    atom->z());
    ofs << buffer << endl;
  }

  return(true);
}

}
