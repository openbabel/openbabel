/**********************************************************************
Copyright (C) 2001-2003 by Geoffrey R. Hutchison

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

bool ReadNWChem(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  string str;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;

  mol.BeginModify();
  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"Output coordinates") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  mol.BeginModify();
	  ifs.getline(buffer,BUFF_SIZE);	// blank
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);	// ---- ----- ----
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 6)
	    {
	      atom = mol.NewAtom();
	      x = atof((char*)vs[3].c_str());
	      y = atof((char*)vs[4].c_str());
	      z = atof((char*)vs[5].c_str());
	      atom->SetVector(x,y,z); //set coordinates

	      //set atomic number
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	} // if "output coordinates"
    } // while 
  mol.EndModify();
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  mol.SetTitle(title);
  return(true);
}

bool WriteNWChem(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  
  ofs << "start molecule" << endl << endl;
  ofs << "title " << endl << " " << mol.GetTitle() << endl << endl;

  ofs << "geometry units angstroms print xyz autosym" << endl;

  OBAtom *atom;
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

  ofs << "end" << endl;

  return(true);
}

}
