 /**********************************************************************
 zipstream.cpp - Test compressed file reading of molecules

 This file is part of the Open Babel project.
 For more information, see <http://openbabel.org/>

 Copyright (C) 2001-2009 Geoffrey R. Hutchison

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
#include <openbabel/obconversion.h>
#include "../src/zipstream.h"

namespace OpenBabel
{
  bool SafeOpen(std::ifstream &fs, const char *filename);
  bool SafeOpen(std::ofstream &fs, const char *filename);
}

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string smilestypes_file = testdatadir + "ziptest.sdf.gz";
#else
   string smilestypes_file = "ziptest.sdf.gz";
#endif

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      cout << "Usage: zipstream\n";
      cout << "   Tests Open Babel compressed file reading." << endl;
      return 0;
    }
  
  std::istream* pIn;
  std::ifstream mifs;
  if (!SafeOpen(mifs, smilestypes_file.c_str()))
    {
      cout << "Bail out! Cannot read test data " << smilestypes_file << endl;
      return -1; // test failed
    }

  zlib_stream::zip_istream *zIn = new zlib_stream::zip_istream(mifs);
  pIn = zIn;
  //  pIn = &mifs;
  OBConversion conv(pIn, &cout);

  if (! conv.SetInFormat("SDF"))
    {
      cout << "Bail out! SDF format is not loaded" << endl;
      return -1;
    }

  unsigned int currentMol = 0;
  OBMol mol;
  std::vector<streampos> offsets; // positions of each molecule

  offsets.push_back(pIn->tellg());
  for (;pIn->good();)
    {
      mol.Clear();
      conv.Read(&mol);
      offsets.push_back(pIn->tellg());
      //      cout << " tellg: " << pIn->tellg() << endl;

      currentMol++;
    }

  // output the number of tests run
  cout << "1.." << (--currentMol * 2) << endl;
  pIn->clear();

  bool success;
  // First test seekoff()
  for (unsigned int i = 0; i < currentMol; ++i) {
    pIn->seekg(offsets[i],  std::ios_base::beg);
    success = conv.Read(&mol);
    if (!success)
      cout << "not ok " << i+1 << endl;
    else
      cout << "ok " << i+1 << " # " << mol.NumAtoms() << endl;
  }
  // Now test seekpos()
  for (unsigned int i = 0; i < currentMol; ++i) {
    pIn->seekg(offsets[i]);
    success = conv.Read(&mol);
    if (!success)
      cout << "not ok " << i+1+currentMol << endl;
    else
      cout << "ok " << i+1+currentMol << " # " << mol.NumAtoms() << endl;
  }

  // Passed Test
  return 0;
}
