/**********************************************************************
obgen.cpp - test program for SMILES 3D coordinate generation
          - using systematic rotor search

Copyright (C) 2006 Tim Vandermeersch
Some portions Copyright (C) 2006 Geoffrey R. Hutchison

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

#include <list>
#include <algorithm>

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/builder.h>

using namespace std;
using namespace OpenBabel;

// PROTOTYPES /////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//! \brief  Generate rough 3D coordinates for SMILES (or other 0D files).
//
int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  string basename, filename = "", option, option2, ff = "MMFF94";

  list<string> argl(argv+1, argv+argc);

  list<string>::iterator optff = find(argl.begin(), argl.end(), "-ff");
  if (optff != argl.end()) {
    list<string>::iterator optffarg = optff;
    ++optffarg;

    if (optffarg != argl.end()) {
      ff = *optffarg;

      argl.erase(optff,++optffarg);
    } else {
      argl.erase(optff);
    }
  }

  if (argl.empty()) {
    cout << "Usage: obgen <filename> [options]" << endl;
    cout << endl;
    cout << "options:      description:" << endl;
    cout << endl;
    cout << "  -ff         select a forcefield" << endl;
    cout << endl;
    OBPlugin::List("forcefields", "verbose");
    exit(-1);
  }

  basename = filename = *argl.begin();
  size_t extPos = filename.rfind('.');

  if (extPos!= string::npos) {
    basename = filename.substr(0, extPos);
  }

  // Find Input filetype
  OBConversion conv;
  OBFormat *format_in = conv.FormatFromExt(filename.c_str());
  OBFormat *format_out = conv.FindFormat("sdf");

  if (!format_in || !format_out || !conv.SetInAndOutFormats(format_in, format_out)) {
    cerr << program_name << ": cannot read input/output format!" << endl;
    exit (-1);
  }

  ifstream ifs;
  ofstream ofs;

  // Read the file
  ifs.open(filename.c_str());
  if (!ifs) {
    cerr << program_name << ": cannot read input file!" << endl;
    exit (-1);
  }

  OBMol mol;

  for (c=1;;c++) {
      mol.Clear();
      if (!conv.Read(&mol, &ifs))
        break;
      if (mol.Empty())
        break;

      OBForceField* pFF = OBForceField::FindForceField(ff);
      if (!pFF) {
        cerr << program_name << ": could not find forcefield '" << ff << "'." <<endl;
        exit (-1);
      }

      //mol.AddHydrogens(false, true); // hydrogens must be added before Setup(mol) is called

      pFF->SetLogFile(&cerr);
      pFF->SetLogLevel(OBFF_LOGLVL_LOW);

      //pFF->GenerateCoordinates();
      OBBuilder builder;
      builder.Build(mol);

      mol.AddHydrogens(false, true); // hydrogens must be added before Setup(mol) is called
      if (!pFF->Setup(mol)) {
        cerr << program_name << ": could not setup force field." << endl;
        exit (-1);
      }

      pFF->SteepestDescent(500, 1.0e-4);
      pFF->WeightedRotorSearch(250, 50);
      pFF->SteepestDescent(500, 1.0e-6);

      pFF->UpdateCoordinates(mol);
      //pFF->ValidateGradients();
      //pFF->SetLogLevel(OBFF_LOGLVL_HIGH);
      //pFF->Energy();


      //char FileOut[32];
      //sprintf(FileOut, "%s_obgen.pdb", basename.c_str());
      //ofs.open(FileOut);
      //conv.Write(&mol, &ofs);
      //ofs.close();
      conv.Write(&mol, &cout);
  } // end for loop

  return(0);
}
