/**********************************************************************
obprobe.cpp - This program will create a grid around a molecule, place
              a probe atom with specified type and partial charge at
              each grid point and calculate the energy. This energy is
              stored in the grid. MMFF94 only for now.

Copyright (C) 2006 Tim Vandermeersch
 
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
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/griddata.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#ifndef _MSC_VER
  #include <unistd.h>
#endif
#include <cstdlib>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
  char *program_name = argv[0];
  char *type;
  int c;
  double step, padding, pchg;
  string basename, filename = "", option;

  step    = 0.5;
  padding = 5.0;

  if (argc < 4) {
    cout << "Usage: obprobe [options] <type> <pchg> <filename>" << endl;
    cout << endl;
    cout << "  <type>: the probe MMFF94 atom type" << endl;
    cout << endl;
    cout << "  <pchg>: the probe's partial charge" << endl;
    cout << endl;
    cout << endl;
    cout << "Options:          Description:" << endl;
    cout << endl;
    cout << "  -s <stepsize>   step size" << endl;
    cout << endl;
    cout << "  -p <padding>    padding" << endl;
    cout << endl;
    cout << "Example probes:" << endl;
    cout << endl;
    cout << "<type>  <pchg>    Description:" << endl;
    cout << "  7     -0.57     carbonyl oxygen (HBA)" << endl;
    cout << " 21      0.4      hydroxyl hydrogen (HBD)" << endl;
    cout << " 37      0.0      phenyl carbon (hydrophobic)" << endl;
    exit(-1);
  } else {
    int ifile = 1;
    for (int i = 1; i < argc; i++) {
      option = argv[i];
      
      if ((option == "-s") && (argc > (i+1))) {
        string stepstr = argv[i+1];
        step = atof(stepstr.c_str());
        ifile += 2;
      }

      if ((option == "-p") && (argc > (i+1))) {
        string paddingstr = argv[i+1];
        padding = atof(paddingstr.c_str());
        ifile += 2;
      }
    }
    
    type = argv[ifile];
    pchg = atof(argv[ifile+1]);
    basename = filename = argv[ifile+2];
    size_t extPos = filename.rfind('.');

    if (extPos!= string::npos) {
      basename = filename.substr(0, extPos);
    }


  }

  // Find Input filetype
  OBConversion conv;
  OBFormat *format_in = conv.FormatFromExt(filename.c_str());
  OBFormat *format_out = conv.FormatFromExt(".cube");
    
  if (!format_in || !format_out || !conv.SetInAndOutFormats(format_in, format_out)) {
    cerr << program_name << ": cannot read input/output format!" << endl;
    exit (-1);
  }

  ifstream ifs;

  // Read the file
  ifs.open(filename.c_str());
  if (!ifs) {
    cerr << program_name << ": cannot read input file!" << endl;
    exit (-1);
  }

  OBForceField* pFF = OBForceField::FindForceField("MMFF94");
  if (!pFF) {
    cerr << program_name << ": could not find forcefield 'MMFF94'." <<endl;
    exit (-1);
  }
  pFF->SetLogFile(&cerr);
  pFF->SetLogLevel(OBFF_LOGLVL_NONE);

  OBMol mol;
  char buffer[BUFF_SIZE];
  for (c=1;;c++) {
    mol.Clear();
    if (!conv.Read(&mol, &ifs))
      break;
    if (mol.Empty())
      break;

    if (!pFF->Setup(mol)) {
      cerr << program_name << ": could not setup force field." << endl;
      exit (-1);
    }
    
    OBGridData* gd = pFF->GetGrid(step, padding, type, pchg);
    mol.SetData(gd);

    ofstream ofs;
    snprintf(buffer, BUFF_SIZE, "%s_%s_%f.cube", basename.c_str(), type, pchg);
    ofs.open(buffer);
    if (!ofs) {
      cerr << program_name << ": cannot read input file!" << endl;
      exit (-1);
    }
    if (!conv.Write(&mol, &ofs)) {
      cerr << program_name << ": could not setup force field." << endl;
    }

    ofs.close();

  } // end for loop

  return(0);
}
