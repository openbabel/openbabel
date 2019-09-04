/**********************************************************************
obenergy.cpp - calculate the energy for a molecule

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
  char *program_name= argv[0];
  int c;
  int verbose = 0;
  bool hydrogens = false;
  string basename, filename = "", option, option2, ff = "MMFF94";
  OBConversion conv;

  if (argc < 2) {
    cout << "Usage: obenergy [options] <filename>" << endl;
    cout << endl;
    cout << "options:      description:" << endl;
    cout << endl;
    cout << "  -v          verbose: print out indivual energy interactions" << endl;
    cout << endl;
    cout << "  -h          add hydrogens before calculating energy" << endl;
    cout << endl;
    cout << "  -ff ffid    select a forcefield" << endl;
    cout << endl;
    cout << "              available forcefields:" << endl;
    cout << endl;
    OBPlugin::List("forcefields", "verbose");
    exit(-1);
  } else {
    int ifile = 1;
    for (int i = 1; i < argc; i++) {
      option = argv[i];

      if (option == "-v") {
        verbose = 1;
        ifile++;
        break;
      }

      if (option == "-h") {
        hydrogens = true;
        ifile++;
      }

      if ((option == "-ff") && (argc > (i+1))) {
        ff = argv[i+1];
        ifile += 2;
      }
    }

    basename = filename = argv[ifile];
    size_t extPos = filename.rfind('.');

    if (extPos!= string::npos) {
      basename = filename.substr(0, extPos);
    }


  }

  // Find Input filetype
  OBFormat *format_in = conv.FormatFromExt(filename.c_str());

  if (!format_in || !conv.SetInFormat(format_in)) {
    cerr << program_name << ": cannot read input format!" << endl;
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

  OBForceField* pFF = OBForceField::FindForceField(ff);
  if (!pFF) {
    cerr << program_name << ": could not find forcefield '" << ff << "'." <<endl;
    exit (-1);
  }
  pFF->SetLogFile(&cout);
  if (verbose)
    pFF->SetLogLevel(OBFF_LOGLVL_HIGH);
  else
    pFF->SetLogLevel(OBFF_LOGLVL_MEDIUM);

  OBMol mol;
  double energy;
  for (c=1;;c++) {
    mol.Clear();
    if (!conv.Read(&mol, &ifs))
      break;
    if (mol.Empty())
      break;

    if (hydrogens)
      mol.AddHydrogens();

    if (!pFF->Setup(mol)) {
      cerr << program_name << ": could not setup force field." << endl;
      exit (-1);
    }

    energy = pFF->Energy(false);
    if (!isfinite(energy)) {
      cerr << " Title: " << mol.GetTitle() << endl;
      FOR_ATOMS_OF_MOL(atom, mol) {
        cerr << " x: " << atom->x() << " y: " << atom->y() << " z: " << atom->z() << endl;
      }
    }

  } // end for loop

  return(0);
}

/* obenergy man page*/
/** \page calculate the energy for a molecule
*
* \n
* \par SYNOPSIS
*
* \b obenergy [options] \<filename\>
*
* \par DESCRIPTION
*
* The obenergy tool can be used to calculate the energy for molecules
* inside (multi-)molecule files (e.g., MOL2, etc.)
*
* \par OPTIONS
*
* If no filename is given, obenergy will give all options including the
* available forcefields.
*
* \b -v:
*     Verbose: print out all individual energy interactions \n\n
* \b -ff \<forcefield\>:
*     Select the forcefield \n\n
*
* \par EXAMPLES
*  - View the possible options, including available forcefields:
*   obenergy
*  - Calculate the energy for the molecule(s) in file test.mol2:
*   obenergy test.mol2
*  - Calculate the energy for the molecule(s) in file test.mol2 using the Ghemical forcefield:
*   obenergy -ff Ghemical test.mol2
*  - Calculate the energy for the molecule(s) in file test.mol2 and print out all individual energy interactions:
*    obenergy -v test.mol2
*
* \par AUTHORS
*
* The obenergy program was contributed by \b Tim \b Vandermeersch.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.org/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 2007 by Tim Vandermeersch. \n \n
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation version 2 of the License.\n \n
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
* \par SEE ALSO
*   The web pages for Open Babel can be found at: http://openbabel.org/ \n
*   The web pages for Open Babel Molecular Mechanics can be found at:
*   http://openbabel.org/wiki/Molecular_mechanics \n
**/
