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


bool ReadAmberPrep(istream &ifs,OBMol &mol,char *title)
{
  char buffer[BUFF_SIZE];
  string str,str1;
  OBAtom *atom;
  OBInternalCoord *coord;
  vector<string> vs;
  vector<OBInternalCoord*> internals;

  ttab.SetFromType("XYZ");
  mol.BeginModify();

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      tokenize(vs,buffer);
      if (vs.size() > 8)
	{
	  atom = mol.NewAtom();
	  coord = new OBInternalCoord();
	  if (mol.NumAtoms() > 1)
	    coord->_a = mol.GetAtom(atoi(vs[4].c_str()));
	  if (mol.NumAtoms() > 2)
	    coord->_b = mol.GetAtom(atoi(vs[5].c_str()));
	  if (mol.NumAtoms() > 3)
	    coord->_c = mol.GetAtom(atoi(vs[6].c_str()));
	  coord->_dst = atof(vs[7].c_str());
	  coord->_ang = atof(vs[8].c_str());
	  coord->_tor = atof(vs[9].c_str());
	  internals.push_back(coord);

	  atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));
	  ttab.SetToType("INT"); ttab.Translate(str,vs[1]); 
	  atom->SetType(str);

	  if (!ifs.getline(buffer,BUFF_SIZE)) break;
	  tokenize(vs,buffer);
	}
    }
  InternalToCartesian(internals,mol);
  mol.EndModify();

  mol.ConnectTheDots();
  mol.SetTitle(title);
  return(true);
}

} // namespace
