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

bool WriteQChem(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  OBAtom *atom;
  
  ofs << "$comment" << endl;
  ofs << mol.GetTitle() << endl;
  ofs << "$end" << endl;
  ofs << endl << "$molecule" << endl << "0 1" << endl;

  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    ofs << atom->GetAtomicNum() << " "
	<< atom->GetX() << " " << atom->GetY() << " " << atom->GetZ() << endl;
  }
  ofs << "$end" << endl;
  ofs << endl << "$rem" << endl << "$end" << endl;

  return(true);
}

bool ReadQChem(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  string str,str1;
  double x,y,z;
  OBAtom *atom;
  OBInternalCoord *coord;
  vector<string> vs;
  vector<OBInternalCoord *> internals; // If we get a z-matrix
  int index;
  bool hasPartialCharges = false;

  mol.BeginModify();

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"Standard Nuclear Orientation") != NULL)
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
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));
	      x = atof((char*)vs[2].c_str());
	      y = atof((char*)vs[3].c_str());
	      z = atof((char*)vs[4].c_str());
	      atom->SetVector(x,y,z);

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      else if(strstr(buffer,"Mulliken Net Atomic Charges") != NULL)
	{
	  hasPartialCharges = true;
	  ifs.getline(buffer,BUFF_SIZE);	// (blank)
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);	// -----------------
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() >= 3)
	    {
	      atom = mol.GetAtom(atoi(vs[0].c_str()));
	      atom->SetPartialCharge(atof(vs[2].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      // This is still not working completely, so it is ignored for now
      // In principle, the final geometry in cartesians is exactly the same
      // as this geometry, so we shouldn't lose any information
//      else if(strstr(buffer,"OPTIMIZATION CONVERGED") != NULL)
//	{
// 	  // mol.EndModify();
// 	  mol.Clear();
// 	  mol.BeginModify();
// 	  // Unfortunately this will come in a Z-matrix
// 	  ifs.getline(buffer,BUFF_SIZE); // *****
// 	  ifs.getline(buffer,BUFF_SIZE); // blank
// 	  ifs.getline(buffer,BUFF_SIZE); // Z-matrix Print:
// 	  ifs.getline(buffer,BUFF_SIZE); // $molecule
// 	  ifs.getline(buffer,BUFF_SIZE); // Charge,Spin
	  
// 	  ifs.getline(buffer,BUFF_SIZE);
// 	  index = 1;
// 	  while (strstr(buffer, "$end") == NULL)
// 	    {
// 	      tokenize(vs,buffer);
// 	      atom = mol.NewAtom();
// 	      atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));

// 	      tokenize(vs,buffer);
// 	      coord = new OBInternalCoord();
// 	      if (index > 1)
// 		{
// 		  coord->_a = mol.GetAtom(atoi(vs[2].c_str()));
// 		  coord->_dst = atof(vs[3].c_str());
// 		}
// 	      if (index > 2)
// 		{
// 		  coord->_b = mol.GetAtom(atoi(vs[4].c_str()));
// 		  coord->_ang = atof(vs[5].c_str()) * DEG_TO_RAD;
// 		}
// 	      if (index > 3)
// 		{
// 		  coord->_c = mol.GetAtom(atoi(vs[6].c_str()));
// 		  coord->_tor = atof(vs[7].c_str()) * DEG_TO_RAD;
// 		}
// 	      index++;
// 	      ifs.getline(buffer,BUFF_SIZE);
// 	      internals.push_back(coord);
// 	    }
// 	  InternalToCartesian(internals, mol);
//	}
    } // while
  mol.EndModify();
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  if (hasPartialCharges)
    mol.SetPartialChargesPerceived();

  mol.SetTitle(title);
  return(true);
}

} // namespace
