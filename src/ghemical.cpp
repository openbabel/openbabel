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

bool ReadGhemical(istream &ifs,OBMol &mol, const char *title)
{
  int i;
  int natoms, nbonds;
  char buffer[BUFF_SIZE];
  string str,str1;
  float x,y,z;
  OBAtom *atom;
  vector<string> vs;
  char bobuf[100];
  string bostr;
  int bgn,end,order;
  bool hasPartialCharges = false;

  mol.BeginModify();

  // Get !Header line with version number
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%*s %*s %d", &i);
  if (!i || i > 100) return false;

  // Get !Info line with number of coord sets
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%*s %d", &i);
  if (!i || i != 1) return false;

  // Get !Atoms line with number
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%*s %d", &natoms);
  if (!natoms) return(false);

  for (i = 1; i <= natoms; i ++)
  {
    if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
    tokenize(vs,buffer);
    if (vs.size() != 2) return(false);
    atom = mol.NewAtom();
    atom->SetAtomicNum(atoi(vs[1].c_str()));
  }

  // Get !Bonds line with number
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%*s %d", &nbonds);
  if (nbonds != 0)
    for (i = 0; i < nbonds; i++)
      {
	if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
	if (!sscanf(buffer,"%d%d%s",&bgn,&end,bobuf)) return (false);
	bostr = bobuf; order = 1;
	if      (bostr == "D")   order = 2;
	else if (bostr == "T")   order = 3;
	else if (bostr == "C")   order = 5; // Conjugated ~= Aromatic
	mol.AddBond(bgn+1,end+1,order);
      }

  // Get !Coord line
  ifs.getline(buffer,BUFF_SIZE);
  for (i = 1; i <= natoms; i ++)
  {
    if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
    tokenize(vs,buffer);
    if (vs.size() != 4) return(false);
    atom = mol.GetAtom(i);
    x = 10.0*atof((char*)vs[1].c_str());
    y = 10.0*atof((char*)vs[2].c_str());
    z = 10.0*atof((char*)vs[3].c_str());
    atom->SetVector(x,y,z); //set coordinates
  }

  if (ifs.getline(buffer,BUFF_SIZE) && strstr(buffer, "!Charges") != NULL)
    {
      hasPartialCharges = true;
      for (i = 1; i <= natoms; i ++)
	{
	  if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
	  tokenize(vs,buffer);
	  if (vs.size() != 2) return(false);
	  atom = mol.GetAtom(i);
	  atom->SetPartialCharge(atof((char*)vs[1].c_str()));
	}
    }

  mol.EndModify();
  if (hasPartialCharges)
    mol.SetPartialChargesPerceived();
  mol.SetTitle(title);
  return(true);
}

bool WriteGhemical(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  char bond_char;
  
  // Ghemical header -- here "version 1.0" format
  ofs << "!Header mm1gp 100" << endl;

  // Number of coordinate sets
  ofs << "!Info 1" << endl;

  // Atom definitions
  sprintf(buffer,"!Atoms %d", mol.NumAtoms());
  ofs << buffer << endl;

  OBAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    ofs << (i - 1) << " " << atom->GetAtomicNum() << endl;
  }

  // Bond definitions
  sprintf(buffer, "!Bonds %d", mol.NumBonds());
  ofs << buffer << endl;

  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;

  for (bond = mol.BeginBond(j); bond; bond = mol.NextBond(j))
    {
      switch(bond->GetBO())
	{
	case 1 :  bond_char = 'S';  break;
	case 2 :  bond_char = 'D';  break;
	case 3 :  bond_char = 'T';  break;
	case 4 :  bond_char = 'C';  break;
	case 5 :  bond_char = 'C';  break;
	default : bond_char = 'S';
	}
      sprintf(buffer,"%d %d %c",
	      bond->GetBeginAtomIdx()-1,
	      bond->GetEndAtomIdx()-1,
	      bond_char);
      ofs << buffer << endl;
    }

  // Coordinate sets (here only 1)
  ofs << "!Coord" << endl;

  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%d %f %f %f",
	    i-1,
	    atom->GetX()/10.0,
	    atom->GetY()/10.0,
	    atom->GetZ()/10.0);

    ofs << buffer << endl;
  }

  // Calculated charges
  ofs << "!Charges" << endl;

  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%d %f",
	    i-1,
	    atom->GetPartialCharge());

    ofs << buffer << endl;
  }

  ofs << "!End" << endl;

  return(true);

}

}
