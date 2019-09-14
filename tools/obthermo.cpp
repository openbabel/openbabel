/**********************************************************************
obthermo.cpp - extract the thermochemistry for a molecule

Copyright (C) 2015 David van der Spoel

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
#include <openbabel/data_utilities.h>
#include <openbabel/pointgroup.h>
#include <cstdlib>
#ifndef _MSC_VER
  #include <unistd.h>
#endif

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
  char  *program_name = argv[0];
  int    Nsymm = 0, Nrot = 0;
  bool   bKJ   = false;
  double dBdT  = 0;
  string filename, option;
  OBConversion conv;
  double unit_factor = 1;
  string e_unit("kcal/mol");
  string s_unit("cal/mol K");
  
  if (argc < 2) {
    cout << "Usage: obthermo [options] <filename>" << endl;
    cout << endl;
    cout << "options:      description:" << endl;
    cout << endl;
    cout << "  --symm N    override symmetry number used in input file" << endl;
    cout << endl;
    cout << "  --nrot N    number of rotatable bonds for conformational entropy" << endl;
    cout << endl;
    cout << "  --dbdt x    temperature derivative of second virial coefficient for cp calculation" << endl;
    cout << endl;
    cout << "  --kj        output kJ/mol related units (default kcal/mol)" << endl;
    cout << endl;
    exit(-1);
  } else {
    int i;
    for (i = 1; i < argc; ) {
      option = argv[i];

      if ((option == "--symm") && (argc > (i+1))) {
        Nsymm = atoi(argv[i+1]);
        if (Nsymm < 1) {
            cerr << program_name << ": the symmetry number should be >= 1!" << endl;
            exit(-1);
        }
        i += 2;
      }
      else if ((option == "--nrot") && (argc > (i+1))) {
        Nrot = atoi(argv[i+1]);
        if (Nrot < 0) {
            cerr << program_name << ": the number of rotatable bonds should be >= 0!" << endl;
            exit(-1);
        }
        i += 2;
      }
      else if ((option == "--dbdt") && (argc > (i+1))) {
        dBdT = atof(argv[i+1]);
        if (dBdT < 0) {
            cerr << program_name << ": the derivative of the second virial coefficient with respect to temperature should be >= 0!" << endl;
            exit(-1);
        }
        i += 2;
      }
      else if (option == "--kj") {
        bKJ          = true;
        unit_factor  = 4.184;
        e_unit.assign("kJ/mol");
        s_unit.assign("J/mol K");
        i += 1;
      }
      else {
        filename.assign(argv[i]);
        i += 1;
      }
    }
  }
  if (filename.size() == 0) {
    cerr << program_name << ": no filename specified" << endl;
    exit (-1);
  }
  // Find Input filetype
  OBFormat *format_in = conv.FormatFromExt(filename.c_str());

  if (!format_in || !conv.SetInFormat(format_in)) {
    cerr << program_name << ": cannot read input format in file \"" << filename << "\"" << endl;
    exit (-1);
  }

  ifstream ifs;

  // Read the file
  ifs.open(filename.c_str());
  if (!ifs) {
    cerr << program_name << ": cannot read input file!" << endl;
    exit (-1);
  }
  OBMol mol;
  if ((conv.Read(&mol, &ifs)) && ! mol.Empty())
  {
      OBPointGroup obPG;
      double temperature, DeltaHf0, DeltaHfT, DeltaGfT, DeltaSfT, S0T, CVT, CPT, ZPVE;
      std::vector<double> Scomponents;
      
      obPG.Setup(&mol);
      printf("obthermo - extract thermochemistry data from quantum chemistry logfiles\n");
      printf("Number of rotatable bonds: %d\n", Nrot);
      if (dBdT == 0)
      {
          printf("Please supply --dbdt option to get reliable heat capacity at constant pressure.\n");
      }
      printf("Point group according to OpenBabel: %s\n", 
             obPG.IdentifyPointGroup());
      bool bVerbose = true;
      if (extract_thermochemistry(mol, 
                                  bVerbose,
                                  &Nsymm,
                                  Nrot,
                                  dBdT,
                                  &temperature,
                                  &DeltaHf0,
                                  &DeltaHfT,
                                  &DeltaGfT,
                                  &DeltaSfT,
                                  &S0T,
                                  &CVT,
                                  &CPT,
                                  Scomponents,
                                  &ZPVE))
      {
          double Rgas  = 1.9872041; // cal/mol K
          printf("DeltaHform(0K)  %10g  %s\n", DeltaHf0*unit_factor, e_unit.c_str());
          printf("Temperature     %10g  K\n", temperature);
          printf("DeltaHform(T)   %10g  %s\n", DeltaHfT*unit_factor, e_unit.c_str());
          printf("DeltaGform(T)   %10g  %s\n", DeltaGfT*unit_factor, e_unit.c_str());
          printf("DeltaSform(T)   %10g  %s\n", DeltaSfT*unit_factor, s_unit.c_str());
          printf("cv(T)           %10g  %s\n", CVT*unit_factor, s_unit.c_str());
          printf("cp(T)           %10g  %s\n", CPT*unit_factor, s_unit.c_str());
          printf("Strans(T)       %10g  %s\n", Scomponents[0]*unit_factor, s_unit.c_str());
          printf("Srot(T)         %10g  %s\n", Scomponents[1]*unit_factor, s_unit.c_str());
          printf("Svib(T)         %10g  %s\n", Scomponents[2]*unit_factor, s_unit.c_str());
          if (Scomponents[3] != 0) 
          {
              printf("Ssymm           %10g  %s\n", Scomponents[3]*unit_factor, s_unit.c_str());
          }
          if (Scomponents[4] != 0) 
          {
              printf("Sconf           %10g  %s\n", Scomponents[4]*unit_factor, s_unit.c_str());
          }
          printf("S0(T)           %10g  %s\n", S0T*unit_factor, s_unit.c_str());
      }
      else
      {
          printf("Could not find all necessary information to determine thermochemistry values.\n");
      }
  }
  ifs.close();
  
  return 0;
}

/* obthermo man page*/
/** \page extract thermochemistry data from a quantum chemistry calculation
*
* \n
* \par SYNOPSIS
*
* \b obthermo [options] \<filename\>
*
* \par DESCRIPTION
*
* The obthermo tool can be used to extract the thermochemistry, e.g. Delta H formation
* Delta G formation and the standard entropy from e.g. a Gaussian calculation.
*
* \par OPTIONS
*
* If no filename is given, obtherm will give all options
*
* \b --symm \<N\>:
*     Set the symmetry number manually and correct results if needed\n\n
*
* \b --nrot \<N\>:
*     Set the number of rotatable bond manually and correct results if needed\n\n
*
* \b --dbdt \<x\>:
*     Set the derivative of the second virial coefficient with temperature for cp calculation\n\n
*
* \b --kj:
*     Use Joules instead of Calories in the ouput\n\n
*
* \par EXAMPLES
*  - View the possible options:
*   obthermo
*  - Print thermochemistry for the molecule in file test.log:
*   obthermo test.log
*  - Id. but correct the symmetry number
*   obenergy --symm 12 test.log
*
* \par AUTHORS
*
* The obthermo program was contributed by \b David \b van \b der \b Spoel.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.org/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 2015 by David van der Spoel. \n \n
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
**/
