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

bool ReadUnichem(istream &ifs,OBMol &mol,char *title)
{
  int i;
  int natoms;
  char buffer[BUFF_SIZE];

  ifs.getline(buffer,BUFF_SIZE);
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%d", &natoms);
  if (!natoms) return(false);

  mol.ReserveAtoms(natoms);

  string str;
  float x,y,z;
  OBAtom *atom;
  vector<string> vs;

  for (i = 1; i <= natoms; i ++)
  {
    if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
    tokenize(vs,buffer);
    if (vs.size() != 4) return(false);
    atom = mol.NewAtom();
    x = atof((char*)vs[1].c_str());
    y = atof((char*)vs[2].c_str());
    z = atof((char*)vs[3].c_str());
    atom->SetVector(x,y,z); //set coordinates

    //set atomic number
    atom->SetAtomicNum(atoi((char*)vs[0].c_str()));
  }
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  mol.SetTitle(title);
  return(true);
}

bool WriteUnichem(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];

  ofs << mol.GetTitle() << endl;
  ofs << mol.NumAtoms() << endl;

  OBAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%3d%15.5f%15.5f%15.5f",
	    atom->GetAtomicNum(),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer << endl;
  }

  return(true);
}

}
