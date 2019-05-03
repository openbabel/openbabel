/**********************************************************************
obminimize.cpp - minimize the energy of a molecule

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
#include <cstdlib>
#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/obutil.h>
#ifndef _MSC_VER
  #include <unistd.h>
#endif

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  int steps = 2500;
  double crit = 1e-6;
  bool sd = false;
  bool cut = false;
  bool newton = false;
  bool hydrogens = false;
  double rvdw = 6.0;
  double rele = 10.0;
  int freq = 10;
  string basename, filename = "", option, option2, ff = "MMFF94";
  char *oext;
  OBConversion conv;
  OBFormat *format_out = conv.FindFormat("pdb"); // default output format

  if (argc < 2) {
    cout << "Usage: obminimize [options] <filename>" << endl;
    cout << endl;
    cout << "options:      description:" << endl;
    cout << endl;
    cout << "  -c crit     set convergence criteria (default=1e-6)" << endl;
    cout << endl;
    cout << "  -cg         use conjugate gradients algorithm (default)" << endl;
    cout << endl;
    cout << "  -sd         use steepest descent algorithm" << endl;
    cout << endl;
    cout << "  -newton     use Newton2Num linesearch (default=Simple)" << endl;
    cout << endl;
    cout << "  -ff ffid    select a forcefield:" << endl;
    cout << endl;
    cout << "  -h          add hydrogen atoms" << endl;
    cout << endl;
    cout << "  -n steps    specify the maximum numer of steps (default=2500)" << endl;
    cout << endl;
    cout << "  -cut        use cut-off (default=don't use cut-off)" << endl;
    cout << endl;
    cout << "  -rvdw rvdw  specify the VDW cut-off distance (default=6.0)" << endl;
    cout << endl;
    cout << "  -rele rele  specify the Electrostatic cut-off distance (default=10.0)" << endl;
    cout << endl;
    cout << "  -pf freq    specify the frequency to update the non-bonded pairs (default=10)" << endl;
    cout << endl;
    OBPlugin::List("forcefields", "verbose");
    exit(-1);
  } else {
    int ifile = 1;
    for (int i = 1; i < argc; i++) {
      option = argv[i];

      // steps
      if ((option == "-n") && (argc > (i+1))) {
        steps = atoi(argv[i+1]);
        ifile += 2;
      }
      // vdw cut-off
      if ((option == "-rvdw") && (argc > (i+1))) {
        rvdw = atof(argv[i+1]);
        ifile += 2;
      }
      // ele cut-off
      if ((option == "-rele") && (argc > (i+1))) {
        rele = atof(argv[i+1]);
        ifile += 2;
      }
      // pair update frequency
      if ((option == "-pf") && (argc > (i+1))) {
        freq = atoi(argv[i+1]);
        ifile += 2;
      }
      // steepest descent
      if (option == "-sd") {
        sd = true;
        ifile++;
      }
      // enable cut-off
      if (option == "-cut") {
        cut = true;
        ifile++;
      }
      // enable Newton2Num
      if (option == "-newton") {
        newton = true;
        ifile++;
      }

      if (strncmp(option.c_str(), "-o", 2) == 0) {
        oext = argv[i] + 2;
        if(!*oext) {
          oext = argv[++i]; //space left after -o: use next argument
          ifile++;
        }

        format_out = conv.FindFormat(oext);
        ifile++;
      }

      if (option == "-h") {
        hydrogens = true;
        ifile++;
      }

      if (option == "-cg") {
        sd = false;
        ifile++;
      }

      if ((option == "-c") && (argc > (i+1))) {
        crit = atof(argv[i+1]);
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
  OBFormat *format_in = conv.FormatFromExt(filename.c_str());

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

  OBForceField* pFF = OBForceField::FindForceField(ff);
  if (!pFF) {
    cerr << program_name << ": could not find forcefield '" << ff << "'." <<endl;
    exit (-1);
  }

  // set some force field variables
  pFF->SetLogFile(&cerr);
  pFF->SetLogLevel(OBFF_LOGLVL_LOW);
  pFF->SetVDWCutOff(rvdw);
  pFF->SetElectrostaticCutOff(rele);
  pFF->SetUpdateFrequency(freq);
  pFF->EnableCutOff(cut);
  if (newton)
    pFF->SetLineSearchType(LineSearchType::Newton2Num);

  OBMol mol;

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

    bool done = true;
    OBStopwatch timer;
    timer.Start();
    if (sd) {
      pFF->SteepestDescentInitialize(steps, crit);
    } else {
      pFF->ConjugateGradientsInitialize(steps, crit);
    }

    unsigned int totalSteps = 1;
    while (done) {
      if (sd)
        done = pFF->SteepestDescentTakeNSteps(1);
      else
        done = pFF->ConjugateGradientsTakeNSteps(1);
      totalSteps++;

      if (pFF->DetectExplosion()) {
        cerr << "explosion has occured!" << endl;
        conv.Write(&mol, &cout);
        return(1);
      } else
        pFF->GetCoordinates(mol);
    }
    double timeElapsed = timer.Elapsed();

    pFF->GetCoordinates(mol);

    conv.Write(&mol, &cout);
    cerr << "Time: " << timeElapsed << "seconds. Iterations per second: " <<  double(totalSteps) / timeElapsed << endl;
  } // end for loop

  return(0);
}

//
// THIS MAN PAGE IS NO LONGER UP TO DATE:
// SEE DOC/*.1 OR WIKI PAGES
//

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
