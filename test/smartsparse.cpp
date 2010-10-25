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

#include <openbabel/parsmart.h>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string smarts_file = testdatadir + "validsmarts.txt";
#else
   string smarts_file = "files/validsmarts.txt";
#endif

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);
  
  cout << "# Testing SMARTS Parsing...  \n";

  std::ifstream ifs;
  if (!SafeOpen(ifs, smarts_file.c_str()))
    {
      cout << "Bail out! Cannot read " << smarts_file << endl;
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
