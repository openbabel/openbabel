/**********************************************************************
invalid-smarts.cpp - Test SMARTS pattern parsing (rejecting invalid patterns)

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

Some portions Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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

#include <openbabel/mol.h>
#include <openbabel/obutil.h>
#include <openbabel/oberror.h>
#include <openbabel/parsmart.h>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string htestdatadir = TESTDATADIR;
  string hinvalid_file = htestdatadir + "invalid-smarts.txt";
  string hrandom1_file = htestdatadir + "random";
  string hrandom2_file = htestdatadir + "random2";
  string hrandom3_file = htestdatadir + "random3";
#else
  string hinvalid_file = "files/invalid-smarts.txt";
  string hrandom1_file = "files/random";
  string hrandom2_file = "files/random2";
  string hrandom3_file = "files/random3";
#endif

int invalidsmarts(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  cout << "# Testing invalid SMARTS parsing..." << endl;

  // make sure to kill off all error reporting
  obErrorLog.StopLogging();

  unsigned int currentTest = 0;
  OBSmartsPattern smarts;
  string pattern, buffer;

  std::ifstream mifs;
  if (!SafeOpen(mifs, hinvalid_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << hinvalid_file << endl;
      return -1; // test failed
    }

  while (getline(mifs, pattern))
    {
      // each line is a potential test pattern

      if (smarts.Init(pattern))
        cout << "not ok " << ++currentTest << " # pattern was parsed "
             << " but should have been rejected\n";
      else
        cout << "ok " << ++currentTest << "\n";
    }

  mifs.close();
  mifs.clear();

  // random file#1
  if (!SafeOpen(mifs, hrandom1_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << hrandom1_file << endl;
      return -1; // test failed
    }

  pattern.clear();
  while (getline(mifs, buffer))
    pattern += buffer;

  if (smarts.Init(pattern))
    cout << "not ok " << ++currentTest << " # random data #1 was parsed "
         << " but should have been rejected\n";
  else
    cout << "ok " << ++currentTest << " # random data #1 \n";
  cout << "# read " << pattern.size() << "\n";

  mifs.close();
  mifs.clear();

  // random2
  if (!SafeOpen(mifs, hrandom2_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << hrandom2_file << endl;
      return -1; // test failed
    }

  pattern.clear();
  while (getline(mifs, buffer))
    pattern += buffer;
  if (smarts.Init(pattern))
    cout << "not ok " << ++currentTest << " # random data #2 was parsed "
         << " but should have been rejected\n";
  else
    cout << "ok " << ++currentTest << " # random data #2 \n";
  cout << "# read " << pattern.size() << "\n";

  mifs.close();
  mifs.clear();

  // random3
  if (!SafeOpen(mifs, hrandom3_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << hrandom3_file << endl;
      return -1; // test failed
    }

  pattern.clear();
  while (getline(mifs, buffer))
    pattern += buffer;
  if (smarts.Init(pattern))
    cout << "not ok " << ++currentTest << " # random data #3 was parsed "
         << " but should have been rejected\n";
  else
    cout << "ok " << ++currentTest << " # random data #3 \n";
  cout << "# read " << pattern.size() << "\n";

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}
