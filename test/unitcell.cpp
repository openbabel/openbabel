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

int main(int argc,char *argv[])
{
  float a, b, c, alpha, beta, gamma;
  vector3 v1, v2, v3;
  float x, y, z;
  char buffer[BUFF_SIZE];
  std::ifstream ifs;

  if (!SafeOpen(ifs,"unitcell.txt")) return(false);
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%f %f %f",&x, &y, &z);
  v1.Set(x, y, z);
  
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%f %f %f",&x, &y, &z);
  v2.Set(x, y, z);

  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%f %f %f",&x, &y, &z);
  v3.Set(x, y, z);

  a = v1.length();
  b = v2.length();
  c = v3.length();
  cout << " a: " << a << " b: " << b << " c: " << c << endl;
  alpha = vectorAngle(v2, v3);
  beta = vectorAngle(v1, v3);
  gamma = vectorAngle(v1, v2);
  cout << " alpha : " << alpha << " beta: " << beta << " gamma: " << gamma << endl;;

  v1.Set(a, 0.0f, 0.0f);
  v2.Set(b*cos(DEG_TO_RAD*gamma), b*sin(DEG_TO_RAD*gamma), 0.0f);
  v3.Set(c*cos(DEG_TO_RAD*beta)*sin(DEG_TO_RAD*alpha),
	 c*sin(DEG_TO_RAD*beta)*cos(DEG_TO_RAD*alpha),
	 c*sin(DEG_TO_RAD*beta)*sin(DEG_TO_RAD*alpha));
  cout << " v1: " << v1 << endl;
  cout << " v2: " << v2 << endl;
  cout << " v3: " << v3 << endl;

  a = v1.length();
  b = v2.length();
  c = v3.length();
  cout << " a: " << a << " b: " << b << " c: " << c << endl;
  alpha = vectorAngle(v2, v3);
  beta = vectorAngle(v1, v3);
  gamma = vectorAngle(v1, v2);
  cout << " alpha : " << alpha << " beta: " << beta << " gamma: " << gamma << endl;

  return(0);
}
