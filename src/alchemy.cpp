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

namespace OpenBabel {

bool ReadAlchemy(istream &ifs,OBMol &mol,char *title)
{
  int i;
  int natoms,nbonds;
  char buffer[BUFF_SIZE];

  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%d %*s %d", &natoms, &nbonds);
  if (!natoms) return(false);

  mol.ReserveAtoms(natoms);
  ttab.SetFromType("ALC");

  string str;
  float x,y,z;
  OBAtom *atom;
  vector<string> vs;

  for (i = 1; i <= natoms; i ++)
  {
    if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
    tokenize(vs,buffer);
    if (vs.size() != 6) return(false);
    atom = mol.NewAtom();
    x = atof((char*)vs[2].c_str());
    y = atof((char*)vs[3].c_str());
    z = atof((char*)vs[4].c_str());
    atom->SetVector(x,y,z); //set coordinates

    //set atomic number
    ttab.SetToType("ATN"); ttab.Translate(str,vs[1]);
    atom->SetAtomicNum(atoi(str.c_str()));
    //set type
    ttab.SetToType("INT"); ttab.Translate(str,vs[1]); 
    atom->SetType(str);
  }

  char bobuf[100];
  string bostr;
  int bgn,end,order;

  for (i = 0; i < nbonds; i++)
  {
    if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
    sscanf(buffer,"%*d%d%d%s",&bgn,&end,bobuf);
    bostr = bobuf; order = 1;
    if      (bostr == "DOUBLE")   order = 2;
    else if (bostr == "TRIPLE")   order = 3;
    else if (bostr == "AROMATIC") order = 5;
    mol.AddBond(bgn,end,order);
  }

  return(true);
}

bool WriteAlchemy(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  char bond_string[10];
  
  sprintf(buffer,"%5d ATOMS, %5d BONDS,     0 CHARGES",
	  mol.NumAtoms(),
	  mol.NumBonds());
  ofs << buffer << endl;
  ttab.SetFromType("INT"); ttab.SetToType("ALC");

  OBAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    str = atom->GetType();
    ttab.Translate(str1,str);
    sprintf(buffer,"%5d %-6s%8.4f %8.4f %8.4f     0.0000",
	    i,
	    (char*)str1.c_str(),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer << endl;
  }

  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;

  for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
  {
    switch(bond->GetBO())
    {
    case 1 :  strcpy(bond_string,"SINGLE");  break;
    case 2 :  strcpy(bond_string,"DOUBLE");  break;
    case 3 :  strcpy(bond_string,"TRIPLE");  break;
    case 5 :  strcpy(bond_string,"AROMATIC"); break;
    default : strcpy(bond_string,"SINGLE");
    }
    sprintf(buffer,"%5d  %4d  %4d  %s",
	    bond->GetIdx()+1,
	    bond->GetBeginAtomIdx(),
	    bond->GetEndAtomIdx(),
	    bond_string);
    ofs << buffer << endl;
  }
  return(true);
}

}

