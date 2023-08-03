/**********************************************************************
multicml.cpp - Test handling of cml files with multiple molecules.

Copyright (C) 2015 David R. Koes
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

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
#include <cstdlib>
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenBabel;

//The XML reader was reading up to </cml>.
//When Read was called on this, it would not return false.

// 1
// reads in molecules from a file
int multicmltest(int argc, char* argv[])
{

  int choice = 1;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  if(choice != 1) //eh, not bothering to split this up
  {
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  cout << endl << "# Testing handling of multi-molecule cml files...  " << endl;
 
  #ifdef TESTDATADIR
    string testdatadir(TESTDATADIR);
    string infile = testdatadir + "c3.cml";
  #else
    string infile = "files/c3.cml";
  #endif

  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  ifstream ifs(infile.c_str());
  if (!ifs)
    {
      cout << "Bail out! Cannot read c3.cml!" << endl;
      return(-1);
    }
  
  OBMol mol;
  OBConversion conv(&ifs);
  conv.SetInFormat("CML");
  int cnt = 0;
  while(conv.Read(&mol))
  {
    cnt++;
  }

  if(cnt != 3)
  {
    cout << "Wrong number of molecules found in file " << cnt << "\n";
    return -1;
  }
  return(0);
}
