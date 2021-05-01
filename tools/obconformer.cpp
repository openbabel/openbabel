/**********************************************************************
obconformer.cpp - Run a Monte Carlo conformer search for a molecule

Copyright (C) 2005-2007 Geoffrey R. Hutchison

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
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obutil.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>

#include <fstream>
#include <iostream>

using namespace std;
using namespace OpenBabel;

int main(int argc, char* argv[]) {
  if (argc != 4 && argc != 5) {
    cout << "Usage: obconformer NSteps GeomSteps <file> [forcefield]" << endl;
    return (-1);
  }

  const int weightSteps = atoi(argv[1]);
  const int geomSteps = atoi(argv[2]);

  ifstream ifs{argv[3]};
  if (!ifs) {
    cerr << "Error! Cannot read input file!" << endl;
    return -1;
  }

  OBConversion conv{&ifs, &cout};
  OBFormat* pFormat = conv.FormatFromExt(argv[3]);

  if (!pFormat) {
    cerr << "Error! Cannot read file format!" << endl;
    return -1;
  }

  if (!conv.SetInAndOutFormats(pFormat, pFormat)) {
    cerr << "Error! File format isn't loaded" << endl;
    return -1;
  }

  // use this if a user doesn't specify forcefield
  const string default_forcefield = "MMFF94";
  // use this if a user doesn't specify forcefield and MMFF94 parameters are not found
  const string fallback_forcefield = "UFF";

  const bool is_forcefield_supplied = argc == 5;
  const string forcefield = is_forcefield_supplied ? argv[4] : default_forcefield;
  OBForceField* pFF = OBForceField::FindForceField(forcefield);
  if (!pFF) {
    cerr << "Error! Cannot find forcefield '" << forcefield << "'" << endl;
    return -1;
  }

  pFF->SetLogFile(&cerr);
  pFF->SetLogLevel(OBFF_LOGLVL_LOW);

  OBMol mol;
  while (ifs.peek() != EOF && ifs.good()) {
    mol.Clear();
    conv.Read(&mol);

    if (pFF->Setup(mol)) {
      pFF->WeightedRotorSearch(weightSteps, geomSteps);
      pFF->ConjugateGradients(geomSteps);  // final cleanup
      pFF->UpdateCoordinates(mol);
      conv.Write(&mol);
    } else {
      cerr << "Error! Cannot set up force field." << endl;
      if (is_forcefield_supplied) {
        return 1;
      }
      pFF = OBForceField::FindForceField(fallback_forcefield);
      assert(pFF);
      cerr << "Force field is switched to " << fallback_forcefield << '.' << endl;
      if (!pFF->Setup(mol)) {
        cerr << "Error! Cannot set up force field." << endl;
        return 1;
      }
      pFF->WeightedRotorSearch(weightSteps, geomSteps);
      pFF->ConjugateGradients(geomSteps);  // final cleanup
      pFF->UpdateCoordinates(mol);
      conv.Write(&mol);
      pFF = OBForceField::FindForceField(forcefield);  // switch back to MMFF94
    }
  }  // while reading molecules

  return 0;
}
