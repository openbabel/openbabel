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

using namespace std;

namespace OpenBabel {

#define BOHR_TO_ANGSTROM 0.529177

bool ReadMPQC(istream &ifs,OBMol &mol,char *title)
{
  char buffer[BUFF_SIZE];
  string str,str1;
  float x,y,z;
  OBAtom *atom;
  vector<string> vs;
  bool bohr = true;

  mol.BeginModify();
  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"<Molecule>:") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  mol.BeginModify();
	  while	(strstr(buffer,"atoms geometry") == NULL)
	    {
	      if (strstr(buffer,"angstrom") != NULL)
		bohr = false;
	      ifs.getline(buffer,BUFF_SIZE);
	    }
	  ifs.getline(buffer,BUFF_SIZE); // Now we're on the atoms
	  tokenize(vs,buffer);
	  while (vs.size() == 6)
	    {
	      if (bohr)
		{
		  x = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
		  y = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
		  z = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
		}
	      else {
		  x = atof((char*)vs[2].c_str());
		  y = atof((char*)vs[3].c_str());
		  z = atof((char*)vs[4].c_str());
		}
	      atom = mol.NewAtom();
	      atom->SetVector(x,y,z);
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
    }
  mol.EndModify();

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  mol.SetTitle(title);
  
  return(true);
}

}
