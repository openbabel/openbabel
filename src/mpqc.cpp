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

#define BOHR_TO_ANGSTROM 0.529177249

bool ReadMPQC(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  string str,str1;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;
  bool bohr = true;
  int i = 0;

  mol.BeginModify();
  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"<Molecule>:") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  while	(strstr(buffer,"geometry") == NULL)
	    {
	      if (strstr(buffer,"angstrom") != NULL)
		bohr = false;
	      if (!ifs.getline(buffer,BUFF_SIZE))
		return(false);
	    }
	  ifs.getline(buffer,BUFF_SIZE); // Now we're on the atoms
	  tokenize(vs,buffer);
	  while (vs.size() == 6)
	    {
	      if (bohr)
		{
		  x = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
		  y = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
		  z = atof((char*)vs[5].c_str()) * BOHR_TO_ANGSTROM;
		}
	      else {
		  x = atof((char*)vs[3].c_str());
		  y = atof((char*)vs[4].c_str());
		  z = atof((char*)vs[5].c_str());
		}
	      atom = mol.NewAtom();
	      atom->SetVector(x,y,z);
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
    }

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  mol.SetTitle(title);
  
  return(true);
}

}
