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

namespace OpenBabel {

#define BOHR_TO_ANGSTROM 0.529177

bool ReadGAMESS(istream &ifs,OEMol &mol,char *title)
{
  char buffer[BUFF_SIZE];
  string str,str1;
  float x,y,z;
  OEAtom *atom;
  vector<string> vs;
  bool hasPartialCharges = false;

  ttab.SetFromType("ATN");
  ttab.SetToType("INT");
  mol.BeginModify();
  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"ATOMIC                      COORDINATES (BOHR)") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  mol.BeginModify();
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 5)
	    {
	      atom = mol.NewAtom();
	      atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
	      x = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
	      y = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
	      z = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
	      atom->SetVector(x,y,z);
	      vs[1].erase(vs[1].size() - 2, 2);
	      ttab.Translate(str,vs[1]); 
	      atom->SetType(str);

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      else if(strstr(buffer,"COORDINATES OF ALL ATOMS ARE (ANGS)") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  mol.BeginModify();
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);	// ---------------
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 5)
	    {
	      atom = mol.NewAtom();
	      atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
	      x = atof((char*)vs[2].c_str());
	      y = atof((char*)vs[3].c_str());
	      z = atof((char*)vs[4].c_str());
	      atom->SetVector(x,y,z);
	      vs[1].erase(vs[1].size() - 2, 2);
	      ttab.Translate(str,vs[1]); 
	      atom->SetType(str);

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      else if(strstr(buffer,"MOPAC CHARGES") != NULL)
	{
	  hasPartialCharges = true;
	  ifs.getline(buffer,BUFF_SIZE);	// ---------------
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 4)
	    {
	      atom = mol.GetAtom(atoi(vs[0].c_str()));
	      atom->SetPartialCharge(atof(vs[2].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
    }
  mol.EndModify();
  mol.ConnectTheDots();

  if (hasPartialCharges)
    mol.SetPartialChargesPerceived();
  mol.SetTitle(title);
  return(true);
}

bool WriteGAMESS(ostream &ofs,OEMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];

  ofs << " $CONTRL COORD=CART UNITS=ANGS $END" << endl;
  ofs << " $DATA" << endl;
  ofs << mol.GetTitle() << endl;
  ofs << "Put symmetry info here" << endl << endl;
  
  OEAtom *atom;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%-3s %4d.0    %8.5f  %8.5f  %8.5f ",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetAtomicNum(),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer << endl;
  }

  ofs << " $END" << endl << endl << endl;
  return(true);
}

}
