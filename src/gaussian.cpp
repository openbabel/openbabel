/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2004 by Geoffrey R. Hutchison

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

// Gaussian input is only tested for Gaussian98. If
// you encounter problems or success with earlier version
// please report it. 2002-07-23, Michael Banck<mbanck@gmx.net>
bool WriteGaussianCart(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  int charge = 0;
  unsigned int multiplicity = 0;
  char buffer[BUFF_SIZE];
  
  ofs << "%" << endl << '\045';
  ofs << "#" << endl << endl;
  ofs << " " << mol.GetTitle() << endl << endl;

  OBAtom *atom;
  string str,str1;
  sprintf(buffer,"%d  %d", mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity());
  ofs << buffer << endl;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    if (atom->GetIsotope() == 0)
      sprintf(buffer,"%-3s      %10.5f      %10.5f      %10.5f ",
	      etab.GetSymbol(atom->GetAtomicNum()),
	      atom->GetX(), atom->GetY(), atom->GetZ());
    else
      sprintf(buffer,"%-3s(Iso=%d) %10.5f      %10.5f      %10.5f ",
	      etab.GetSymbol(atom->GetAtomicNum()),
	      atom->GetX(), atom->GetY(), atom->GetZ());
    ofs << buffer << endl;
  }
  // file should end with a blank line
  ofs << endl;
  return(true);
}

// Reading Gaussian output has been tested for G98 and G03 to some degree
// If you have problems (or examples of older output), please contact
// the openbabel-discuss@lists.sourceforge.net mailing list and/or post a bug
bool ReadGaussianOut(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  string str,str1;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;
  int index, charge = 0;
  unsigned int spin = 1;
  bool hasPartialCharges = false;

  mol.BeginModify();

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if (strstr(buffer,"Symbolic Z-matrix:") != NULL)
	{
 	  ifs.getline(buffer,BUFF_SIZE); // Charge = # Multiplicty = #
	  tokenize(vs, buffer, " \t\n");
	  if (vs.size() == 6)
	    {
	      charge = atoi(vs[2].c_str());
	      spin = atoi(vs[5].c_str());
	    }

 	  ifs.getline(buffer,BUFF_SIZE);
	}
      else if(strstr(buffer,"Coordinates (Angstroms)") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  mol.BeginModify();
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);	// ---------------
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 6)
	    {
	      atom = mol.NewAtom();
	      atom->SetAtomicNum(atoi((char*)vs[1].c_str()));
	      x = atof((char*)vs[3].c_str());
	      y = atof((char*)vs[4].c_str());
	      z = atof((char*)vs[5].c_str());
	      atom->SetVector(x,y,z);

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      else if(strstr(buffer,"Total atomic charges") != NULL ||
	      strstr(buffer,"Mulliken atomic charges") != NULL)
	{
	  hasPartialCharges = true;
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() >= 3 && 
		 strstr(buffer,"Sum of ") == NULL)
	    {
	      atom = mol.GetAtom(atoi(vs[0].c_str()));
	      atom->SetPartialCharge(atof(vs[2].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
    } // end while
  mol.EndModify();
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  if (hasPartialCharges)
    mol.SetPartialChargesPerceived();
  mol.SetTotalCharge(charge);
  mol.SetTotalSpinMultiplicity(spin);

  mol.SetTitle(title);
  return(true);
}

}
