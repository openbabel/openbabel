/**********************************************************************
Copyright (C) 2003 by Geoffrey R. Hutchison

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
  bool SafeOpen(std::ifstream &fs, char *filename);
  bool SafeOpen(std::ofstream &fs, char *filename);
}

using namespace std;
using namespace OpenBabel;

bool TestUnitCell()
{
  double a, b, c, alpha, beta, gamma;
  vector3 v1, v2, v3;
  double x = 0.0, y = 0.0, z = 0.0;
  char buffer[BUFF_SIZE];
  std::ifstream ifs;
  OBUnitCell cell, cell2;
  vector<vector3> v3Return;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string unitcell_file = testdatadir + "unitcell.txt";
#else
  string unitcell_file = "unitcell.txt";
#endif

  if (!SafeOpen(ifs, (char*)unitcell_file.c_str())) return(false);
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%lf %lf %lf",&x, &y, &z);
  v1.Set(x, y, z);
  
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%lf %lf %lf",&x, &y, &z);
  v2.Set(x, y, z);

  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%lf %lf %lf",&x, &y, &z);
  v3.Set(x, y, z);

  cell.SetData(v1, v2, v3);
  a = cell.GetA();
  b = cell.GetB();
  c = cell.GetC();
  alpha = cell.GetAlpha();
  beta = cell.GetBeta();
  gamma = cell.GetGamma();

  v3Return = cell.GetCellVectors();

  cell2.SetData(v1, v2, v3);
  a = cell2.GetA();
  b = cell2.GetB();
  c = cell2.GetC();
  alpha = cell2.GetAlpha();
  beta = cell2.GetBeta();
  gamma = cell2.GetGamma();

  if (cell.GetA() != cell2.GetA() || cell.GetB() != cell2.GetB() ||
      cell.GetC() != cell2.GetC() || cell.GetAlpha() != cell2.GetAlpha() ||
      cell.GetBeta() != cell2.GetBeta() || cell.GetGamma() != cell2.GetGamma())
    return(false);

  return(true);
}
