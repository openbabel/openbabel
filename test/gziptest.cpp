/**********************************************************************
gziptest.cpp - Test handling of gzipped files.

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

//reads molecule(s) from file  and compares to results
//outputs error and terminates if not consistent
void checkResults(const string& file, const vector<string>& correctResults)
{
  if(correctResults.size() == 0)
    return;
  OBMol mol;
  OBConversion conv(file);
  conv.SetOutFormat(conv.GetInFormat(), conv.GetInGzipped());
  vector<string> results;
  while(conv.Read(&mol))
  {
    OBConversion strconv;
    string cansmi;
    if (conv.GetInFormat() == conv.FindFormat("mol2")) {
      // Don't test roundtripping for mol2 as the atom types are changed by the reader (TODO: Fix this)
      strconv.SetOutFormat("can");
      strconv.SetOptions("n", conv.OUTOPTIONS);
      cansmi = strconv.WriteString(&mol, true);
    }
    else {
      // Testing roundtrip through gz with Read/WriteString
      OBMol mol2;
      string gzippedstr = conv.WriteString(&mol);
      strconv.SetInAndOutFormats(conv.GetInFormat(), conv.FindFormat("CAN"), conv.GetInGzipped());
      strconv.SetOptions("n", OBConversion::OUTOPTIONS);
      strconv.ReadString(&mol2, gzippedstr);
      cansmi = strconv.WriteString(&mol2, true);
    }
    results.push_back(cansmi);
  }

  if(results.size() != correctResults.size())
  {
    cout << "Incorrect number of outputs (" << results.size() << " vs " << correctResults.size() << ") for " << file << "\n";
    exit(-1);
  }
  for(unsigned i = 0, n = results.size(); i < n; i ++)
  {
    if(results[i] != correctResults[i])
    {
      cout << file << " match failure: " << results[i] << " != " << correctResults[i] << "\n";
      exit(-1);
    }
  }
}

// 1
// reads gzip.infiles, which has a file name on each line,
// converts each molecule into a gzipped string, converts that
// string into a molecule, converts that molecule into canonical smiles
// and compares to gzip.out

int gziptest(int argc, char* argv[])
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

  cout << endl << "# Testing gzip handling...  " << endl;
 
  #ifdef TESTDATADIR
    string testdatadir(TESTDATADIR);
    string gzipin = testdatadir + "gzip.in";
  #else
    string gzipin = "files/gzip.in";
  #endif

  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  ifstream ifs(gzipin.c_str());
  if (!ifs)
    {
      cout << "Bail out! Cannot read gzip.in!" << endl;
      return(-1);
    }
  
  //gzip.in specifies what files to read and their canonical smiles
  string line;
  string filepath;
  vector<string> correctResults; //read from file

  OBConversion conv;
  conv.SetInAndOutFormats("smi", "can");
  conv.SetOptions("n", conv.OUTOPTIONS);
  OBMol mol;

  while(getline(ifs, line))
  {
    if(line.length() == 0)
      continue; //ignore blank lines
    if(line[0] == '#')
    {
      checkResults(filepath, correctResults);
      correctResults.clear();
      //next token should be file name
      stringstream str(line);
      string fname;
      str >> fname; //# character
      if(fname.length() == 1)
      {
        str >> fname; // actual file name
      }
      else
      {
        fname = fname.substr(1); //strip #
      }
#ifdef TESTDATADIR
    filepath = testdatadir + fname;
 #else
    filepath = "files/" + fname;
#endif
    }
    else
    {
      conv.ReadString(&mol, line);
      std::string can = conv.WriteString(&mol, true);
      correctResults.push_back(can);
    }
  }
  checkResults(filepath, correctResults);

  return(0);
}
