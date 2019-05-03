/**********************************************************************
obrotate = rotate a tortional bond matched by a SMART pattern
Copyright (C) 2003 Fabien Fontaine
Some portions Copyright (C) 2004-2005 Geoffrey R. Hutchison
Some portions Copyright (C) 2008 Tim Vandermeersch

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

/*
  Require a SMART pattern, a file containing molecule coordinates
  4 atoms of the SMART pattern to define the tortional, an angle value
  The angle value must be in degree
  the 2 atoms of the rotating bond must be bonded but the 2 others not
  the part of the molecule on the side of the second atom is kept fixed
  whereas the part on the side of the third atom is rotated.
  example of command line:
  obrotate "[nH]ccccc[O,C][C,O]" test.sdf 1 6 7 8 180.0
*/


// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/parsmart.h>
#include <openbabel/rotamer.h>
//#include <unistd.h>
#include <openbabel/obconversion.h>
#include <cstdlib>
using namespace std;
using namespace OpenBabel;

///////////////////////////////////////////////////////////////////////////////
//! \brief Set a tortional bond to a given angle
int main(int argc,char **argv)
{
  OBAtom *a1, *a2, *a3, *a4;
  unsigned int smartor[4]= {0,0,0,0};// atoms of the tortional in the SMART
  float angle =   0;      // tortional angle value to set in degree
  char *FileIn =NULL, *Pattern=NULL;
  unsigned int i, t, errflg = 0;
  int c;
  string err;
  bool changeAll = false; // default to only change the last matching torsion

  // parse the command line -- optional -a flag to change all matching torsions
  if (argc < 8 || argc > 9) {
    errflg++;
  } else {
    // Fetch the option and shift values after the option
    if (argc == 9) {
      int curArg = 0;
      while (curArg < 9) {
        if (strcmp(argv[curArg], "-a") == 0) {
          changeAll = true;
          break;
        }
        ++curArg;
      }
      // We expect -a and so changeAll should be true
      if (!changeAll)
        errflg++;

      // now let's shift values
      while (curArg < 8) {
        argv[curArg] = argv[curArg+1];
      }
    }
    FileIn = argv[2];
    Pattern = argv[1];
    // Read the atom position
    for(i=3, t=0; i<7; ++i, ++t) {
      c = sscanf(argv[i], "%u", &smartor[t]);
      if (c != 1) {
        errflg++;  // error in arguments, quit and warn user
        break;
      }
    }
    c = sscanf(argv[7], "%f", &angle);
    if (c != 1) {
      errflg++; // error in arguments, quit and warn user
    }
    if (argc == 9) {
      if (strcmp(argv[8], "-a") == 0)
        changeAll = true;
      else
        errflg++; // error in arguments, quit and warn user
    }
  }

  if (errflg) {
    cerr << "Usage: obrotate \"PATTERN\" <filename> <atom1> <atom2> <atom3> <atom4> <angle> [-a]" << endl;
    exit(-1);
  }

  // create pattern
  OBSmartsPattern sp;
  sp.Init(Pattern);
  if (sp.NumAtoms() < 4) {
    cerr << "obrotate: The number of atoms in the SMART pattern must be higher than 3." << endl;
    exit(-1);
  }

  for (i=0; i<4; ++i) {
    if ( smartor[i] < 1 || smartor[i] > sp.NumAtoms()) {
      cerr << "obrotate: The torsional atom values must be between 1 and "
           <<  sp.NumAtoms()
           << ", which is the number of atoms in the SMART pattern." << endl;
      exit(-1);
    }
  }

  OBConversion conv; //NF...
  OBFormat* format = conv.FormatFromExt(FileIn);
  if(!(format && conv.SetInAndOutFormats(format, format))) { //in and out formats same
    cerr << "obrotate: cannot read and/or write this file format!" << endl;
    exit (-1);
  } //...NF

  //Open the molecule file
  ifstream ifs;

  // Read the file
  ifs.open(FileIn);
  if (!ifs) {
    cerr << "obrotate: cannot read input file!" << endl;
    exit (-1);
  }

  OBMol mol;
  vector< vector <int> > maplist;      // list of matched atoms
  vector< vector <int> >::iterator m;  // and its iterators
  //   int tindex;

  // Set the angles
  for (;;) {
    mol.Clear();
    //NF      ifs >> mol;                   // Read molecule
    conv.Read(&mol,&ifs); //NF
    if (mol.Empty())
      break;

    if (sp.Match(mol)) {
      // if match perform rotation
      maplist = sp.GetUMapList(); // get unique matches

      if (maplist.size() > 1)
        cerr << "obrotate: Found " << maplist.size() << " matches. Only last one will be rotated." << endl;

      // look at all the mapping atom but save only the last one.
      for (m = maplist.begin(); m != maplist.end(); ++m) {
        a1 = mol.GetAtom( (*m)[ smartor[0] - 1] );
        a2 = mol.GetAtom( (*m)[ smartor[1] - 1] );
        a3 = mol.GetAtom( (*m)[ smartor[2] - 1] );
        a4 = mol.GetAtom( (*m)[ smartor[3] - 1] );
        if (changeAll)
          mol.SetTorsion(a1, a2, a3, a4, angle * DEG_TO_RAD);
      }

      if ( !a2->IsConnected(a3) ) {
        cerr << "obrotate: The atoms of the rotating bond must be bonded." << endl;
        exit(-1);
      }

      if (!changeAll)
        mol.SetTorsion(a1, a2, a3, a4, angle * DEG_TO_RAD);
    } else {
      cerr << "obrotate: Found 0 matches for the SMARTS pattern." << endl;
      exit(-1);
    }
    //NF      cout << mol;
    conv.Write(&mol,&cout); //NF
  }

  return(0);
}


/* obrotate man page*/
/** \page obrotate batch-rotate dihedral angles matching SMARTS patterns
*
* \n
* \par SYNOPSIS
*
* \b obrotate '<SMARTS-pattern>' \<filename\> \<atom1\> \<atom2\> \<atom3\> \<atom4\> \<angle\>
*
* \par DESCRIPTION
*
* The obrotate program rotates the torsional (dihedral) angle of a specified
* bond in molecules to that defined by the user. In other words, it does the
* same as a user setting an angle in a molecular modelling package, but much
* faster and in batch mode.
* \n\n
* The four atom IDs required are indexes into the SMARTS pattern, which starts
* at atom 1. The angle supplied is in degrees. The two atoms used to set
* the dihedral angle \<atom1\> and \<atom4\> do not need to be connected
* to the atoms of the bond \<atom2\> and \<atom3\> in any way.
*\n\n
* The order of the atoms matters -- the portion of the molecule attached to
* \<atom1\> and \<atom2\> remain fixed, but the portion bonded to \<atom3\> and
& \<atom4\> moves.
*
* \par EXAMPLES
*  - Let's say that you want to define the conformation of a large number of
*  molecules with a pyridyl scaffold and substituted with an aliphatic chain
*  at the 3-position, for example for docking or 3D-QSAR purposes.
* \n\n
*    To set the value of the first dihedral angle to 90 degrees:\n
*   obrotate "c1ccncc1CCC" pyridines.sdf 5 6 7 8 90
* \n
* Here 6 and 7 define the bond to rotate in the SMARTS patter, i.e., c1-C and
* atoms 5 and 8 define the particular dihedral angle to rotate.
*  - Since the atoms to define the dihedral do not need to be directly
*  connected, the nitrogen in the pyridine can be used:\n
*   obrotate "c1ccncc1CCC" pyridines.sdf 4 6 7 8 90
*
*  - Keep the pyridyl ring fixed and moves the aliphatic chain:\n
*   obrotate "c1ccncc1CCC" pyridines.sdf 5 6 7 8 90
*  - Keep the aliphatic chain fixed and move the pyridyl ring:\n
*   obrotate "c1ccncc1CCC" pyridines.sdf 8 7 6 5 90
*
* \par AUTHORS
*
* The obrotate program was contributed by \b Fabien \b Fontaine.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.org/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
*  Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison \n \n
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
*   A guide for constructing SMARTS patterns can be found at: http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
**/
