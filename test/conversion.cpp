/**********************************************************************
conversion.cpp - Unit tests for Open Babel OBConversion class

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <stdio.h>
#include <iostream>
#include <cstdlib>

using namespace std;
using namespace OpenBabel;

int conversion(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  cout << "# Unit tests for OBConversion \n";

  // the number of tests for "prove"
  cout << "1..9\n";

  cout << "ok 1\n"; // for loading tests

  OBMol obMol;
  OBConversion obConversion;
  obConversion.SetInAndOutFormats("smi", "mdl");
  cout << "ok 2\n";

  obConversion.ReadString(&obMol, "C1=CC=CS1");
  cout << "ok 3\n";

  if (obMol.NumAtoms() == 5) {
    cout << "ok 4\n";
  } else {
    cout << "not ok 4\n";
  }

  obMol.AddHydrogens();
  if (obMol.NumAtoms() == 9) {
    cout << "ok 5\n";
  } else {
    cout << "not ok 5\n";
  }

  if ( (obConversion.WriteString(&obMol)).length() > 0)
    cout << "ok 6\n";
  else
    cout << "not ok 6\n";

  // PR#1474265
  obConversion.WriteFile(&obMol, "test.mdl");
  ifstream ifs("test.mdl");
  if (ifs.good())
    cout << "ok 7\n";
  else
    cout << "not ok 7\n";

  // gzip input
  // gzip output

  // multi-molecule reading
  // PR#1465586
  // aromatics.smi
  // attype.00.smi

  //ReadFile()
  //Read()
  //WriteString()
  // GetOutputIndex()
  // IsLast

  //ReadString()
  //IsFirstInput
  //Read()

  // splitting
  
  // splitting using gzip-input
  // PR#1357705
  
  // size 0 input
  // PR#1250900
  
  // RegisterFormat
  // FindFormat
  // FormatFromExt
  // FormatFromMIME
  // GetNextFormat
  // GetDefaultFormat

  // BatchFileName
  // IncrementedFileName

  // option handling
  // AddOption
  // IsOption
  // RemoveOption
  // IsOption

  // SetOptions
  // IsOption

  // RegisterOptionParam
  // GetOptionParams

  // GetInStream
  // GetOutStream
  // SetInStream
  // SetOutStream

  // nasty tests
  obConversion.ReadString(&obMol, "");
  obConversion.Read(&obMol);
  cout << "ok 8\n";

  return(0);
}


