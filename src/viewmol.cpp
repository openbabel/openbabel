/**********************************************************************
Copyright (C) 2002 by Geoffrey Hutchison
Based on code Copyright (C) 1999 Jörg-Rüdiger Hill

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
#define ANGSTROM_TO_BOHR 1.889725989

bool ReadViewMol(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  OBAtom *atom;
  double x,y,z, border;
  double factor = 1.0;
  int bgn, end, order;
  vector<string> vs;
  bool foundTitle = false;
  bool foundBonds = false;

  mol.BeginModify();

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if (strstr(buffer,"$title") != NULL)
	{
	  if (!ifs.getline(buffer,BUFF_SIZE)) return (false);
	  mol.SetTitle(buffer);
	  foundTitle = true;
	}
      else if (strstr(buffer,"$coord") != NULL)
	{
	  tokenize(vs,buffer);
	  if (vs.size() == 2)
	    factor = atof((char*)vs[1].c_str()); // conversion to angstrom
	  while (ifs.getline(buffer,BUFF_SIZE))
	    {
	      if (buffer[0] == '$') break;
	      tokenize(vs,buffer);
	      if (vs.size() != 4) return(false);
	      atom = mol.NewAtom();
	      x = atof((char*)vs[0].c_str()) * factor;
	      y = atof((char*)vs[1].c_str()) * factor;
	      z = atof((char*)vs[2].c_str()) * factor;
	      atom->SetVector(x,y,z); //set coordinates
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[3].c_str()));
	    }
	}
      else if (strstr(buffer,"$bonds") != NULL) 
	{
	  foundBonds = true;
	  while (ifs.getline(buffer,BUFF_SIZE))
	    {
	      if (buffer[0] == '$') break;
	      sscanf(buffer,"%d %d %lf",&bgn,&end, &border);
	      if (border > 1.0)
		order = int(border);
	      else
		order = 1;
	      mol.AddBond(bgn+1,end+1,order);
	    }
	}
      else if (strstr(buffer,"$end") != NULL)
	break;
    } // while
  
  mol.EndModify();

  if (!foundTitle)
    mol.SetTitle(title);
  if (!foundBonds)
    {
      mol.ConnectTheDots();
      mol.PerceiveBondOrders();
    }
  return(true);
}

bool WriteViewMol(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  
  if (strlen(mol.GetTitle()) > 0)
    ofs << "$title" << endl << mol.GetTitle() << endl;

  ofs << "$coord 1.0" << endl;

  OBAtom *atom;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%22.14f%22.14f%22.14f %s",
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ(),
	    etab.GetSymbol(atom->GetAtomicNum()));
    ofs << buffer << endl;
  }

  ofs << "$end" << endl;

  return(true);

}

}
