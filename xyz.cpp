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

namespace OpenEye {

bool ReadXYZ(istream &ifs,OEMol &mol,char *title)
{
  int i;
  int natoms;
  char buffer[BUFF_SIZE];

  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%d", &natoms);
  if (!natoms) return(false);

  mol.ReserveAtoms(natoms);
  ttab.SetFromType("XYZ");

  string str;
  float x,y,z;
  OEAtom *atom;
  vector<string> vs;

  ifs.getline(buffer,BUFF_SIZE);
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
    atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
    //set type
    ttab.SetToType("INT"); ttab.Translate(str,vs[0]); 
    atom->SetType(str);
  }
  mol.ConnectTheDots();
  mol.SetTitle(title);
  return(true);
}

bool WriteXYZ(ostream &ofs,OEMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  
  sprintf(buffer,"%d", mol.NumAtoms());
  ofs << buffer << endl;
  sprintf(buffer,"%s\t%15.7f", mol.GetTitle(), mol.GetEnergy());
  ofs << buffer << endl;
  ttab.SetFromType("INT"); ttab.SetToType("XYZ");

  OEAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%3s%15.5f%15.5f%15.5f",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer << endl;
  }

  return(true);
}

}
