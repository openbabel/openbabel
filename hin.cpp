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

bool ReadHIN(istream &ifs,OEMol &mol,char *title)
{ 
  // Right now only read in the first molecule
  int i;
  int max, bo;
  char buffer[BUFF_SIZE];
  string str,str1;
  float x,y,z;
  OEAtom *atom;
  vector<string> vs;

  ttab.SetFromType("XYZ");
  while (strstr(buffer,"mol") == NULL)
    ifs.getline(buffer, BUFF_SIZE);
  ifs.getline(buffer, BUFF_SIZE);
  
  mol.BeginModify();
  while (strstr(buffer,"enmol") == NULL)
    {
      tokenize(vs,buffer); // Don't really know how long it'll be
      if (vs.size() <= 11) break;
      atom = mol.NewAtom();
      atom->SetAtomicNum(etab.GetAtomicNum(vs[3].c_str()));
      x = atof((char*)vs[6].c_str());
      y = atof((char*)vs[7].c_str());
      z = atof((char*)vs[8].c_str());
      atom->SetVector(x,y,z);
      ttab.SetToType("INT"); ttab.Translate(str,vs[3]); 
      atom->SetType(str);
      
      max = 11 + 2 * atoi((char *)vs[10].c_str());
      for (i = 11; i < max; i+=2)
	{
	  switch(((char*)vs[i+1].c_str())[0]) // First char in next token
	    {
	    case 's': bo = 1; break;
	    case 'd': bo = 2; break;
	    case 't': bo = 3; break;
	    case 'a': bo = 5; break;
	    default : bo = 1; break;
	    }
	  mol.AddBond(mol.NumAtoms(), atoi((char *)vs[i].c_str()), bo);
	}
      ifs.getline(buffer, BUFF_SIZE);
    }
  mol.EndModify();

  mol.SetTitle(title);
  return(true);
}

bool WriteHIN(ostream &ofs,OEMol &mol)
{
  unsigned int i, file_num = 1;
  string str,str1;
  char buffer[BUFF_SIZE];
  OEAtom *atom;
  OEBond *bond;
  vector<OEBond*>::iterator j;
  char bond_char;

  ttab.SetFromType("INT"); ttab.SetToType("XYZ");

  ofs << "mol " << file_num << " " << mol.GetTitle() << endl;;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"atom %d - %-3s **  - %8.5f %8.5f  %8.5f  %8.5f %d ",
	    i,
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetPartialCharge(),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ(),
	    atom->GetValence());
    ofs << buffer;
    for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j))
    {
      switch(bond->GetBO())
      {
      case 1 : bond_char = 's'; break;
      case 2 : bond_char = 'd'; break;
      case 3 : bond_char = 't'; break;
      case 5 : bond_char = 'a'; break;
      default: bond_char = 's'; break;
      }
      sprintf(buffer,"%d %c ", (bond->GetNbrAtom(atom))->GetIdx(), bond_char);
      ofs << buffer;
    }
    ofs << endl;
  }
  ofs << "endmol " << file_num << endl;
  return(true);
}

}
