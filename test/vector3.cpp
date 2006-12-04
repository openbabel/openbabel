/**********************************************************************
vector3.cpp - Unit tests for vector3 manipulations.

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include "babelconfig.h"
#include "math/vector3.h"
#include "obutil.h"

#include <iostream>

using namespace std;
using namespace OpenBabel;

OBRandom randomizer;

vector3 randomVector(void)
{
  vector3 A(randomizer.NextFloat() * 2.0 - 1.0,
            randomizer.NextFloat() * 2.0 - 1.0,
            randomizer.NextFloat() * 2.0 - 1.0);

  return A;
}

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);
  
  if (argc != 1)
    {
      cout << "Usage: vector3" << endl;
      cout << "   Tests Open Babel 3D vector manipulations." << endl;
      return 0;
    }
  
  cout << "# Testing vector3 .." << endl;
  
  randomizer.TimeSeed();
  
  bool passedAll = true;
  unsigned int testCount = 0;
  vector3 test1, test2, test3;
  for (unsigned int i = 0; i < 500; ++i)
    {
      test1 = randomVector();
      test2 = test1.normalize();
      if (IsApprox(test2.length(), 1.0, 1.0e-6))
        cout << "ok " << ++testCount << "\n";
      else
        cout << "not ok " << ++testCount << "\n";

      test1 = randomVector();
      test1.createOrthoVector(test3);
      if (IsNegligible(dot(test1, test3), 1.0, 1.0e-6))
        cout << "ok " << ++testCount << "\n";
      else
        cout << "not ok " << ++testCount << "# " << dot(test1, test3) 
             << "\n";
    }

  return 0;
}
