/**********************************************************************
obconformer.cpp - Run a Monte Carlo conformer search for a molecule

Copyright (C) 2005-2007 Geoffrey R. Hutchison
 
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
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/obutil.h>

#include <openbabel/forcefield.h>

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenBabel;

OBRotorList rl;
OBRandom generator;

int main(int argc,char *argv[])
{
  if (argc != 5)
    {
      cout << "Usage: obconformer NWeights NSteps GeomSteps <file>" << endl;
      return(-1);
    }

  int weightSteps, searchSteps, geomSteps;
  weightSteps = atoi(argv[1]);
  searchSteps = atoi(argv[2]);
  geomSteps = atoi(argv[3]);

  ifstream ifs(argv[4]);
  if (!ifs)
    {
      cerr << "Error! Cannot read input file!" << endl;
      return(-1);
    }
  
  OBConversion conv(&ifs, &cout);
  OBFormat* pFormat;
  
  pFormat = conv.FormatFromExt(argv[4]);
  if ( pFormat == NULL )
    {
      cerr << "Error! Cannot read file format!" << endl;
      return(-1);
    }
  
  // Finally, we can do some work!
  OBMol mol;
  if (! conv.SetInAndOutFormats(pFormat, pFormat))
    {
      cerr << "Error! File format isn't loaded" << endl;
      return (-1);
    }
  
  OBForceField *pFF = OBForceField::FindForceField("Ghemical");
  pFF->SetLogFile(&cerr);
  pFF->SetLogLevel(OBFF_LOGLVL_LOW);

  while(ifs.peek() != EOF && ifs.good())
    {
      mol.Clear();
      conv.Read(&mol);

      pFF->Setup(mol);
      pFF->RandomRotorSearch(weightSteps, geomSteps);
      // pFF->SystematicRotorSearch();
      pFF->UpdateCoordinates(mol);
      conv.Write(&mol);
    } // while reading molecules
  
  return(0);
}
