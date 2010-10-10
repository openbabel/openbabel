/**********************************************************************
logppsa.cpp - Unit tests for Open Babel OBLogP and OBPSA class

Copyright (C) 2007 Tim Vandermmeersch
 
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
#include <openbabel/groupcontrib.h>

#include <stdio.h>
#include <iostream>

using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  if (argc != 1)
    {
      cout << "Usage: logp_psa" << endl;
      cout << " Unit tests for OBLogP and OBPSA " << endl;
      return(-1);
    }

  cout << "# Unit tests for OBLogP and OBPSA \n";

  // the number of tests for "prove"
  cout << "1..7\n";

  OBConversion obConversion;
  obConversion.SetInAndOutFormats("smi", "mdl");
  OBMol obMol;
  
  OBDescriptor* obLogP = OBDescriptor::FindType("logP");
  OBDescriptor* obPSA  = OBDescriptor::FindType("TPSA");
  double logP, psa;

  cout << "ok 1\n"; // for loading tests
  
  obConversion.ReadString(&obMol, "Oc1ccccc1OC");
  obMol.AddHydrogens();
  
  logP = obLogP->Predict(&obMol);
  if (IsNear(logP , 1.4008)) { // value from JOELib2
    cout << "ok 2 # " << logP << '\n';
  } else {
    cout << "not ok 2 # " << logP << '\n';
  }
  
  psa = obPSA->Predict(&obMol);
  if (IsNear(psa , 29.46)) { // value from JOELib2
    cout << "ok 3 # " << psa << '\n';
  } else {
    cout << "not ok 3 # " << psa << '\n';
  }

  obConversion.ReadString(&obMol, "c1ccccc1CBr");
  obMol.AddHydrogens();
  
  logP = obLogP->Predict(&obMol);
  if (IsNear(logP, 2.5815)) { // Value from JOELib2
    cout << "ok 4 # " << logP << '\n';
  } else {
    cout << "not ok 4 # " << logP << '\n';
  }
  
  psa = obPSA->Predict(&obMol);
  if (IsNear(psa, 0.0)) { // Value from JOELib2
    cout << "ok 5 # " << psa << '\n';
  } else {
    cout << "not ok 5 # " << psa << '\n';
  }

  obConversion.ReadString(&obMol, "Cc1ccccc1NC(=O)C");
  obMol.AddHydrogens();
  
  logP = obLogP->Predict(&obMol);
  if (IsNear(logP, 2.0264)) { // JOELib2 = 1.9534, more H added on N
    cout << "ok 6 # " << logP << '\n';
  } else {
    cout << "not ok 6 # " << logP << '\n';
  }
  
  psa = obPSA->Predict(&obMol);
  if (IsNear(psa, 29.1)) { // Value from JOELib2
    cout << "ok 7 # " << psa << '\n';
  } else {
    cout << "not ok 7 # " << psa << '\n';
  }

  return(0);
}


