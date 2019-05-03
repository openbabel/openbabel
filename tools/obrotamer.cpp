/**********************************************************************
obrotamer.cpp - Generate a random rotamer for a given molecule

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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

#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>
#include <openbabel/obutil.h>
#include "../src/rand.h"
#include "../src/rand.cpp"

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
  if (argc != 2)
    {
      cout << "Usage: obrotamer <file>" << endl;
      cout << "  Outputs the number of rotable bonds and generate a random"
           << " rotamer" << endl;
      return(-1);
    }

  ifstream ifs(argv[1]);
  if (!ifs)
    {
      cerr << "Error! Cannot read input file!" << endl;
      return(-1);
    }
  
  OBConversion conv(&ifs, &cout);
  OBFormat* pFormat;
  
  pFormat = conv.FormatFromExt(argv[1]);
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
  
  OBRotorList rl;
  OBRotamerList rotamers;
  OBRandom rand;
  rand.TimeSeed();

  while(ifs.peek() != EOF && ifs.good())
    {
      mol.Clear();
      conv.Read(&mol);

      rl.Setup(mol);
      
      cerr << " Number of rotatable bonds: " << rl.Size() << endl;

      // indexed from 1, rotorKey[0] = 0
      std::vector<int> rotorKey(rl.Size() + 1, 0);

      // each entry represents the configuration of a rotor
      // e.g. indexes into OBRotor::GetResolution()
      //       (the different angles to sample from the OBRotorRules database)
      OBRotorIterator ri;
      OBRotor *rotor = rl.BeginRotor(ri);
      for (unsigned int i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri))
        rotorKey[i] = rand.NextInt() % rotor->GetResolution().size();

      rotamers.SetBaseCoordinateSets(mol);
      rotamers.Setup(mol, rl);
      rotamers.AddRotamer(rotorKey);

      // This is more useful when you have added many rotamers
      // rather than the one in this example
      rotamers.ExpandConformerList(mol, mol.GetConformers());
      
      // Change the molecule conformation -- from 0 .. rotamers.NumRotamers()
      mol.SetConformer(0);
      conv.Write(&mol);
    } // while reading molecules
  
  return(0);
}
