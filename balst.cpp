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

bool ReadBallAndStick(istream &ifs,OEMol &mol,char *title)
{
  int i,natoms;
  char buffer[BUFF_SIZE];  

  if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
  if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
  sscanf(buffer,"%d",&natoms);
  mol.ReserveAtoms(natoms);

  float x,y,z;
  OEAtom *atom;
  vector<string> vs;
  vector<string>::iterator j;

  for (i = 1; i <= natoms;i ++)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
      tokenize(vs,buffer);
      if (vs.size() < 4) return(false);
      if (vs[0].size() > 1) vs[0][1] = tolower(vs[0][1]);
      atom = mol.NewAtom();
      x = atof((char*)vs[1].c_str());
      y = atof((char*)vs[2].c_str());
      z = atof((char*)vs[3].c_str());
      atom->SetVector(x,y,z); //set coordinates
      atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
      atom->SetType(vs[0]);

      for (j = vs.begin()+4;j != vs.end();j++)
	mol.AddBond(atom->GetIdx(),atoi((char*)j->c_str()),1);
    }

  //result = assign_radii(mol);
  //result = assign_types(mol);
  //result = build_connection_table(mol);
  //assign_bond_order(mol);
  return(true);
}

bool WriteBallAndStick(ostream &ofs,OEMol &mol)
{ 
  char tmptype[10];
  char buffer[BUFF_SIZE];
  
  if (strlen(mol.GetTitle()) > 0)  ofs << mol.GetTitle() << endl; 
  else                             ofs << "Untitled" << endl;
  
  sprintf(buffer,"%d",mol.NumAtoms()); ofs << buffer << endl;
  
  OEAtom *atom,*nbr;
  vector<OENodeBase*>::iterator i;
  vector<OEEdgeBase*>::iterator j;

  for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      strcpy(tmptype,etab.GetSymbol(atom->GetAtomicNum()));
      if (strlen(tmptype) > 1) tmptype[1] = toupper(tmptype[1]);
      sprintf(buffer,"%-3s %8.4f  %8.4f  %8.4f",
	      tmptype,
	      atom->GetX(),
	      atom->GetY(),
	      atom->GetZ());
      ofs << buffer; 
      for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
	{
	  sprintf(buffer,"%6d",nbr->GetIdx());
	  ofs << buffer;
	}
      ofs << endl;
    }

  return(true);
}

}
