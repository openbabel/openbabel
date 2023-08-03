 /**********************************************************************
 smartsparse.cpp - Test SMARTS parsing for valid patterns

 This file is part of the Open Babel project.
 For more information, see <http://openbabel.org/>

 Copyright (C) 2010 by Geoffrey R. Hutchison

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

#include <openbabel/babelconfig.h>

#include <fstream>
#include <iostream>

#include <openbabel/parsmart.h>
#include <openbabel/obutil.h>

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif


using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string ntestdatadir = TESTDATADIR;
  string nsmarts_file = ntestdatadir + "validsmarts.txt";
#else
   string nsmarts_file = "files/validsmarts.txt";
#endif

int smartsparse(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }
  
  cout << "# Testing SMARTS Parsing...  \n";

  std::ifstream ifs;
  if (!SafeOpen(ifs, nsmarts_file.c_str()))
    {
      cout << "Bail out! Cannot read " << nsmarts_file << endl;
      return -1; // test failed
    }

  //read in the SMARTS test patterns
  char buffer[BUFF_SIZE];
  OBSmartsPattern sp;
  unsigned int patterns = 0;
  for (;ifs.getline(buffer,BUFF_SIZE);)
    {
      if (buffer[0] == '#') // skip comment line
        continue;
      
      if (sp.Init(buffer))
        cout << "ok " << ++patterns << endl;
      else
        cout << "not ok " << ++patterns << " failed on " << buffer << endl;
    }
  ifs.close();
  // output the number of tests run
  cout << "1.." << patterns << endl;

  // Passed Test
  return 0;
}
