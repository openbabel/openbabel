/**********************************************************************
Copyright (C) 1998-2003 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "math/matrix3x3.h"

using namespace std;

namespace OpenBabel
{

  // ShelX homepage http://shelx.uni-ac.gwdg.de/SHELX/
  //   and 	    http://www.msg.ucsf.edu/local/programs/shelxl/SHELX_97.html

bool ReadShelX(istream &ifs,OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  int natoms;
  double A,B,C,Alpha,Beta,Gamma;
  matrix3x3 m;
  
  ifs.getline(buffer,BUFF_SIZE); mol.SetTitle(buffer);
  ifs.getline(buffer,BUFF_SIZE); sscanf(buffer,"%d",&natoms);
  
  while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"CELL",4));

  if (!EQn(buffer,"CELL",4)) return(false);
  vector<string> vs;
  tokenize(vs,buffer," \n\t,");
  if (vs.size() != 7) return(false);

  //parse cell values
  A = atof((char*)vs[1].c_str());
  B = atof((char*)vs[2].c_str());
  C = atof((char*)vs[3].c_str());
  Alpha = atof((char*)vs[4].c_str());
  Beta  = atof((char*)vs[5].c_str());
  Gamma = atof((char*)vs[6].c_str());

  OBUnitCell *uc = new OBUnitCell;
  uc->SetData(A, B, C, Alpha, Beta, Gamma);
  mol.SetData(uc);
  m = uc->GetOrthoMatrix();

  int i;
  double x,y,z;
  char type[10];
  OBAtom *atom;
  vector3 v;

  for (i = 1; i <= natoms;i++)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
      tokenize(vs,buffer," \n\t,");
      if (vs.size() < 4) return(false);
      atom = mol.NewAtom();

      x = atof((char*)vs[1].c_str());
      y = atof((char*)vs[2].c_str());
      z = atof((char*)vs[3].c_str());
      v.Set(x,y,z); v *= m;

      strcpy(type,vs[0].c_str());
      atom->SetAtomicNum(etab.GetAtomicNum(type));
      atom->SetVector(v);
    }

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  return(true);
}

} // end namespace
