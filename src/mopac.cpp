/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 Geoffrey R. Hutchison

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

bool ReadMOPAC(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  string str,str1;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;
  vector<double> charges;
  bool hasPartialCharges = false;
  double energy;

  mol.BeginModify();
  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"CARTESIAN COORDINATES") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  mol.BeginModify();
	  ifs.getline(buffer,BUFF_SIZE);	// blank
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);	// blank
	  ifs.getline(buffer,BUFF_SIZE);	  
	  tokenize(vs,buffer);
	  while (vs.size() == 5)
	    {
	      atom = mol.NewAtom();
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));
	      x = atof((char*)vs[2].c_str());
	      y = atof((char*)vs[3].c_str());
	      z = atof((char*)vs[4].c_str());
	      atom->SetVector(x,y,z);

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      else if(strstr(buffer,"FINAL HEAT") != NULL)
	{
	  sscanf(buffer,"%*s%*s%*s%*s%*s%lf",&energy);
	  mol.SetEnergy(energy);
	}
      else if(strstr(buffer,"NET ATOMIC CHARGES") != NULL)
	{
	  hasPartialCharges = true;
	  charges.clear();
	  ifs.getline(buffer,BUFF_SIZE);	// blank
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 4)
	    {
	      atom = mol.GetAtom(atoi(vs[0].c_str()));
	      atom->SetPartialCharge(atof(vs[2].c_str()));
	      charges.push_back(atof(vs[2].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
    }
  mol.EndModify();
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  if (hasPartialCharges)
    {
      mol.SetPartialChargesPerceived();
      for (unsigned int i = 1; i <= mol.NumAtoms(); i++)
	{
	  atom = mol.GetAtom(i);
	  atom->SetPartialCharge(charges[i-1]);
	}
    }
  mol.SetTitle(title);
  return(true);
}

bool ReadMOPACCartesian(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  string str;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;

  ifs.getline(buffer,BUFF_SIZE); // keywords
  ifs.getline(buffer,BUFF_SIZE); // filename
  ifs.getline(buffer,BUFF_SIZE); // title (currently ignored)

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      tokenize(vs,buffer);
      if (vs.size() == 0) break;
      else if (vs.size() < 7) return false;
      atom = mol.NewAtom();
      x = atof((char*)vs[1].c_str());
      y = atof((char*)vs[3].c_str());
      z = atof((char*)vs[5].c_str());
      atom->SetVector(x,y,z); //set coordinates

      //set atomic number
      atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
    }

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  mol.SetTitle(title);

  return(true);
}

bool WriteMOPACCartesian(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  
  ofs << "PUT KEYWORDS HERE" << endl;
  ofs << endl;
  ofs << mol.GetTitle() << endl;

  OBAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%-3s%8.5f 1 %8.5f 1 %8.5f 1",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer << endl;
  }

  return(true);
}

}
