/**********************************************************************
Copyright (C) 2000 by Geoffrey Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989

namespace OpenBabel {

bool ReadDMol(istream &ifs,OBMol &mol,char *title)
{
  char buffer[BUFF_SIZE];
  string str;
  float x,y,z;
  OBAtom *atom;
  vector<string> vs;

  ttab.SetFromType("XYZ");

  ifs.getline(buffer,BUFF_SIZE);
  while (strstr(buffer,"$coordinates") == NULL)
    ifs.getline(buffer,BUFF_SIZE);

  //  ifs.getline(buffer,BUFF_SIZE);
  while (strstr(buffer,"$end") == NULL)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) break;
      tokenize(vs,buffer);
      if (vs.size() != 4) break;
      atom = mol.NewAtom();
      //set atomic number
      atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
      //set type
      ttab.SetToType("INT"); ttab.Translate(str,vs[0]); 
      atom->SetType(str);
      x = atof((char*)vs[1].c_str()) * BOHR_TO_ANGSTROM;
      y = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
      z = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
      atom->SetVector(x,y,z); //set coordinates
    }
  mol.ConnectTheDots();
  mol.SetTitle(title);
  return(true);
}

bool WriteDMol(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  
  ofs << "$coordinates" << endl;

  OBAtom *atom;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%-3s% 27.14f% 20.14f% 20.14f",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX() * ANGSTROM_TO_BOHR,
	    atom->GetY() * ANGSTROM_TO_BOHR,
	    atom->GetZ() * ANGSTROM_TO_BOHR);
    ofs << buffer << endl;
  }

  ofs << "$end" << endl;

  return(true);
}

}
