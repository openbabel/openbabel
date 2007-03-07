/**********************************************************************
obminimize.cpp - minimize the energy of a molecule

Copyright (C) 2006 Tim Vandermeersch
Some portions Copyright (C) 2006 Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

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
#include <unistd.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  int steps = 2500;
  string basename, filename = "", option, option2, ff = "Ghemical";

  if (argc < 2) {
    cout << "Usage: obminimize [options] <filename>" << endl;
    cout << endl;
    cout << "options:      description:" << endl;
    cout << endl;
    cout << "  -n          specify the maximum numer of steps" << endl;
    cout << endl;
    cout << "  -ff         select a forcefield:" << endl;
    cout << endl;
    FOR_EACH(OBForceField, iter) {
      //cout << "              " << iter.ID() << " - " << iter.Description() << endl;
      cout << "              " << iter.ID() << endl;
    }
    exit(-1);
  } else {
    int ifile = 1;
    for (int i = 1; i < argc; i++) {
      option = argv[i];
      
      if ((option == "-n") && (argc > (i+1))) {
        steps = atoi(argv[i+1]);
	ifile += 2;
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
  OBConversion conv;
  OBFormat *format_in = conv.FormatFromExt(filename.c_str());
  OBFormat *format_out = conv.FindFormat("pdb");
    
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
 
      pFF->SetLogFile(&cout);
      pFF->SetLogLevel(OBFF_LOGLVL_LOW);
      
      if (!pFF->Setup(mol)) {
        cerr << program_name << ": could not setup force field." << endl;
        exit (-1);
      }
      
      pFF->ConjugateGradients(steps);
      pFF->UpdateCoordinates(mol);

      //char FileOut[32];
      //sprintf(FileOut, "%s_obgen.pdb", basename.c_str());
      //ofs.open(FileOut);
      //conv.Write(&mol, &ofs);
      //ofs.close();
      conv.Write(&mol, &cout);
  } // end for loop

  return(1);
}

/* obminimize man page*/
/** \page minimize the energy for a molecule
*
* \n
* \par SYNOPSIS
*
* \b obminimize [options] \<filename\>
*
* \par DESCRIPTION
*
* The obminimize tool can be used to minimize the energy for molecules 
* inside (multi-)molecule files (e.g., MOL2, etc.)
*
* \par OPTIONS
*
* If no filename is given, obminimize will give all options including the
* available forcefields.
*
* \b -n \<steps\>:
*     Specify the maximum number of steps \n\n
* \b -ff \<forcefield\>:
*     Select the forcefield \n\n
*
* \par EXAMPLES
*  - View the possible options, including available forcefields: 
*   obminimize
*  - Minimize the energy for the molecule(s) in file test.mol2:
*   obminimize test.mol2
*  - Minimize the energy for the molecule(s) in file test.mol2 using the Ghemical forcefield:
*   obminimize -ff Ghemical test.mol2 
*  - Minimize the energy for the molecule(s) in file test.mol2 and set the maximum numer of steps to 300:
*    obenergy -n 300 test.mol2
*
* \par AUTHORS
*
* The obminimize program was contributed by \b Tim \b Vandermeersch.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.sourceforge.net/THANKS.shtml
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
*   The web pages for Open Babel can be found at: http://openbabel.sourceforge.net/ \n
*   The web pages for Open Babel Molecular Mechanics can be found at: 
*   http://openbabel.sourceforge.net/wiki/Molecular_mechanics \n
**/
