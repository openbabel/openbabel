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

bool ReadBiosymCAR(istream &ifs,OEMol &mol,char *title)
{
  char buffer[BUFF_SIZE];
  string str;
  float x,y,z;
  OEAtom *atom;
  vector<string> vs;

  ttab.SetFromType("XYZ");
  mol.BeginModify();

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"PBC") != NULL)
	{
	  if(strstr(buffer,"ON") != NULL)
	    {
	      ifs.getline(buffer,BUFF_SIZE);
	      ifs.getline(buffer,BUFF_SIZE);
	      ifs.getline(buffer,BUFF_SIZE);
	    }
	  else
	    {
	      ifs.getline(buffer,BUFF_SIZE);
	      ifs.getline(buffer,BUFF_SIZE);
	    }
	  break;
	}
    }

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"end") != NULL)
	break;

      atom = mol.NewAtom();

      tokenize(vs,buffer);
      ttab.SetToType("ATN"); ttab.Translate(str,vs[7]);
      atom->SetAtomicNum(atoi(str.c_str()));
      x = atof((char*)vs[1].c_str());
      y = atof((char*)vs[2].c_str());
      z = atof((char*)vs[3].c_str());
      atom->SetVector(x,y,z);
      ttab.SetToType("INT"); ttab.Translate(str,vs[7]); 
      atom->SetType(str);
    }
  mol.EndModify();

  mol.ConnectTheDots();
  mol.SetTitle(title);
  return(true);
}

} // namespace
