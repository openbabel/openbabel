/**********************************************************************
obgen.cpp - test program for SMILES 3D coordinate generation
          - using distance geometry

Copyright (C) 2006 Tim Vandermeersch
Some portions Copyright (C) 2006-2012 Geoffrey R. Hutchison

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
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/distgeom.h>

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
  string basename, filename = "";

  if (argc < 2) {
    cout << "Usage: obdistgen <filename>" << endl;
    cout << endl;
    exit(-1);
  } else {
    basename = filename = argv[1];
    size_t extPos = filename.rfind('.');

    if (extPos!= string::npos) {
      basename = filename.substr(0, extPos);
    }
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
  OBForceField* pFF = OBForceField::FindForceField("mmff94");

  for (c=1;;c++) {
      mol.Clear();
      if (!conv.Read(&mol, &ifs))
        break;
      if (mol.Empty())
        break;

      // hydrogens must be added before Setup(mol) is called
      //  to ensure stereochemistry
      mol.AddHydrogens();

      OBDistanceGeometry dg;
      dg.Setup(mol);

      for (unsigned int i = 1; i < 2; ++i) {
        //        cout << i << endl;
        dg.AddConformer();
      }
      dg.GetConformers(mol);
      //      cout << " Conformers: " << mol.NumConformers() << endl;
      // Check the energies
      pFF->Setup(mol);
      // load through all the conformers
      unsigned int minConf = 0;
      double e, minE = 1.0e10;
      for (unsigned int n = 0; n < mol.NumConformers();++n)
        {
          mol.SetConformer(n);
          pFF->SetCoordinates(mol);
          e = pFF->Energy(false);
          //          cout << " Conformer: " << n << " " << e << endl;
          if (e < minE) {
            minE = e;
            minConf = n;
          }
        }
      mol.SetConformer(minConf);
      conv.Write(&mol, &cout);
  } // end for loop

  return(0);
}
