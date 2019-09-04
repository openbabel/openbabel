/**********************************************************************
roundtrip.cpp - Test "roundtrip" results for converting from one molec. format
                to another.

Copyright (C) 2003-2006 Geoffrey R. Hutchison
 
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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/data.h>
#include <openbabel/obconversion.h>

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
  OBConversion conv;
  OBFormat* pFormat1;
  OBFormat* pFormat2;

  if (argc != 3)
    {
      cout << "Usage: roundtrip <file1> <file2>" << endl;
      return(-1);
    }

  pFormat1 = conv.FormatFromExt(argv[1]);
  if ( pFormat1 == NULL )
    {
      cerr << argv[0] << ": Cannot read file #1 format!" << endl;
      return(-1);
    }
  pFormat2 = conv.FormatFromExt(argv[2]);
  if ( pFormat2 == NULL )
    {
      cerr << argv[0] << ": Cannot read file #2 format!" << endl;
      return(-1);
    }

  // Finally, we can do some work!
  OBMol mol, mol2;
  ifstream inFileStream1(argv[1]);
  ifstream inFileStream2(argv[2]);

  if (!inFileStream1)
    {
      cerr << argv[0] << ": Cannot read input file #1!" << endl;
      return(-1);
    }
  else if (!inFileStream2)
    {
      cerr << argv[0] << ": Cannot read input file #2!" << endl;
      return(-1);
    }

  OBAtom *atom1, *atom2;
  OBConversion conv1(&inFileStream1, &cout), conv2(&inFileStream2, &cout);

  if (! conv1.SetInAndOutFormats(pFormat1, pFormat2))
    {
      obErrorLog.ThrowError("", "File format #1 isn't loaded", obInfo);
      return (-1);
    }
  if (! conv2.SetInAndOutFormats(pFormat2, pFormat1))
    {
      obErrorLog.ThrowError("", "File format #2 isn't loaded", obInfo);
      return (-1);
    }

  int molCount = 0;
  while(inFileStream1.peek() != EOF && inFileStream1.good() && 
        inFileStream2.peek() != EOF && inFileStream2.good() )
    {
      mol.Clear();
      mol2.Clear();
      molCount++;
      //      cerr << " read " << molCount << " molecules " << endl;

      if(!conv1.Read(&mol) || !conv2.Read(&mol2))
        break; // make sure to check for failed reads

      const char* p1 = strrchr(argv[1],'.');
      const char* p2 = strrchr(argv[2],'.');

      if (p1 && strncasecmp(p1 + 1, "CML", 3) != 0
          && mol.NumAtoms() == 0)
        {
          cout << " ** ERROR ** molecule " << molCount 
               << " in file #1 has no atoms!" << endl;
          return(-1);
        }

      if (p1 && strncasecmp(p1 + 1, "CML", 3) != 0
          && mol2.NumAtoms() == 0)
        {
          cout << " ** ERROR ** molecule " << molCount 
               << " in file #2 has no atoms!" << endl;
          return(-1);
        }

      if (p1 && strncasecmp(p1 + 1, "BOX", 3) == 0)
        {
          if (mol.NumAtoms() != 8)
            {
              cout << " *** ERROR *** BOX file #1 without 8 atoms!" << endl;
              return(-1);
            }
          return(0);
        }

      if (p2 && strncasecmp(p2 + 1, "BOX", 3) == 0)
        {
          if (mol2.NumAtoms() != 8)
            {
              cout << " *** ERROR *** BOX file #2 without 8 atoms!" << endl;
              return(-1);
            }
          return(0);
        }

      if ( (p1 && strncasecmp(p1 + 1, "SMI", 3) == 0)
           || (p2 && strncasecmp(p2 + 1, "SMI", 3) == 0) )
        {
          if (mol.NumHvyAtoms() != mol2.NumHvyAtoms())
            {
              cout << " ** ERROR ** SMILES Number of heavy atoms differ: "
                   << mol.NumHvyAtoms() << " and " << mol2.NumHvyAtoms() << endl;
              return(-1);
            }
          return(0);
        }
      else
        {
          if (mol.NumAtoms() != mol2.NumAtoms())
            {
              cout << " ** ERROR ** Number of atoms differ: " << mol.NumAtoms()
                   << " and " << mol2.NumAtoms()
                   << " in molecule " << molCount << endl;
              return(-1);
            }
        }
    
      for(unsigned int i = 1;i <= mol.NumAtoms(); i++)
        {
          atom1 = mol.GetAtom(i);
          atom2 = mol2.GetAtom(i);
	
          if (atom1->GetAtomicNum() != atom2->GetAtomicNum())
            {
              cout << " ** ERROR ** Elements for atom " << i << " differ: " <<
                atom1->GetAtomicNum() << " and " << atom2->GetAtomicNum()
                   << " in molecule " << molCount << endl;
              return(-1);
            }
	
          if ( (p1 && strncasecmp(p1 + 1, "SMI", 3) == 0)
               || (p2 && strncasecmp(p2 + 1, "SMI", 3) == 0) )
            if ((atom1->GetX()-atom2->GetX()>1e-1) ||
                (atom1->GetY()-atom2->GetY()>1e-1) ||
                (atom1->GetZ()-atom2->GetZ()>1e-1))
              {
                cout << " ** ERROR ** Coordinates for atom " << i << " differ" 
                     << " in molecule " << molCount << endl;
                return(-1);
              }
        }

    } // while reading molecules

  char buffer[BUFF_SIZE];
  // try skipping any blank lines
  while(inFileStream1.peek() != EOF && inFileStream1.good() && 
        (inFileStream1.peek() == '\n' || inFileStream1.peek() == '\r'))
    inFileStream1.getline(buffer,BUFF_SIZE);

  while(inFileStream2.peek() != EOF && inFileStream2.good() && 
        (inFileStream2.peek() == '\n' || inFileStream2.peek() == '\r'))
    inFileStream2.getline(buffer,BUFF_SIZE);

  // now check to see if there's anything to read
  mol.Clear();
  mol2.Clear();
  if ( inFileStream1.good() && inFileStream1.peek() != EOF &&
       conv1.Read(&mol) )
    {
      if (mol.NumAtoms() > 0)
        {
          cout << " ** ERROR **  File 1 has more molecules! " << endl;
          cout << "   Already read " << molCount << " molecules from both" << endl;
          cout << "   New molecule " << mol.GetTitle()
               << " has " << mol.NumAtoms() << " atoms " << endl;
          exit(-1);
        }
    }
  else if ( inFileStream2.good() && inFileStream2.peek() != EOF &&
            conv2.Read(&mol2) )
    {
      if (mol2.NumAtoms() > 0)
        {
          cout << " ** ERROR **  File 2 has more molecules! " << endl;
          cout << "   Already read " << molCount << " molecules from both" << endl;
          cout << "   New molecule " << mol.GetTitle() 
               << " has " << mol2.NumAtoms() << " atoms " << endl;
          exit(-1);
        }
    }
    
  return(0);
}
