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
#include <ctype.h>

using namespace std;

namespace OpenBabel {

bool ReadJaguar(istream &ifs,OBMol &mol, const char *title)
{
  char buffer[BUFF_SIZE];
  string str,str1;
  double x,y,z;
  unsigned int i;
  OBAtom *atom;
  vector<string> vs;

  mol.BeginModify();
  while (ifs.getline(buffer,BUFF_SIZE))
  {
    if (strstr(buffer,"Input geometry:") != NULL
        || strstr(buffer,"symmetrized geometry:") != NULL
        || strstr(buffer,"new geometry:") != NULL
        || strstr(buffer,"final geometry:") != NULL)
    {
      // mol.EndModify();
      mol.Clear();
      mol.BeginModify();
      ifs.getline(buffer,BUFF_SIZE);  //     angstroms
      ifs.getline(buffer,BUFF_SIZE);  // column headings
      ifs.getline(buffer,BUFF_SIZE);
      tokenize(vs,buffer);
      while (vs.size() == 4)
      {
        atom = mol.NewAtom();
        str = vs[0]; // Separate out the Symbol# into just Symbol ...
        for (i = 0;i < str.size();i++)
          if (isdigit(str[i]))
            str[i] = '\0';

        atom->SetAtomicNum(etab.GetAtomicNum(str.c_str()));
        x = atof((char*)vs[1].c_str());
        y = atof((char*)vs[2].c_str());
        z = atof((char*)vs[3].c_str());
        atom->SetVector(x,y,z);

        if (!ifs.getline(buffer,BUFF_SIZE)) break;
        tokenize(vs,buffer);
      }
    }
    if (strstr(buffer, "Atomic charges from electrostatic potential") != NULL)
    { 
      mol.SetAutomaticPartialCharge(false);
      unsigned int chgcount=0;
      while (chgcount<mol.NumAtoms())
      {
        ifs.getline(buffer,BUFF_SIZE);  // blank line
        ifs.getline(buffer,BUFF_SIZE);  // header line
        ifs.getline(buffer,BUFF_SIZE);  // data line
        tokenize(vs,buffer);
        for (vector<string>::size_type icount=1;icount<vs.size();++icount)
        {
          chgcount=chgcount+1;
          mol.GetAtom(chgcount)->SetPartialCharge(atof((char*)vs[icount].c_str()));
        }
      }
    }
  }
  mol.EndModify();

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  mol.SetTitle(title);
  return(true);
}

bool WriteJaguar(ostream &ofs,OBMol &mol)
{
  unsigned int i;
  char buffer[BUFF_SIZE];
  OBAtom *atom;

  ofs << mol.GetTitle() << endl << endl;
  ofs << "&gen" << endl;
  ofs << "&" << endl;
  ofs << "&zmat" << endl;

  for (i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"  %s%d   %12.7f  %12.7f  %12.7f",
            etab.GetSymbol(atom->GetAtomicNum()), i,
            atom->GetX(),
            atom->GetY(),
            atom->GetZ());
    ofs << buffer << endl;
  }

  ofs << "&" << endl;
  return(true);
}

} // Namespace
