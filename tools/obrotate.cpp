/**********************************************************************
obrotate = rotate a tortional bond matched by a SMART pattern
Copyright (C) 2003 Fabien Fontaine

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

#include "mol.h"
#include "parsmart.h"
#include "binary.h"
#include <unistd.h>

using namespace std;
using namespace OpenBabel;

///////////////////////////////////////////////////////////////////////////////
//! \brief Set a tortional bond to a given angle
int main(int argc,char **argv)
{

  char *program_name=NULL;
  io_type inFileType = UNDEFINED, outFileType = UNDEFINED;
  float res=255.0f/360.0f; // constant to convert degree to unsigned char
  unsigned char tor[4]= {0,0,0,0};// atoms of the tortional in the molecule
  unsigned int smartor[4]= {0,0,0,0};// atoms of the tortional in the SMART
  float angle =   0;      // tortional angle value to set in degree
  int nrotamers = 1;
  int nrotors   = 1;
  unsigned char ang_array[2] = {0,0}; // coordinate set of reference and 
                                     //tortional angle
  char *FileIn =NULL, *Pattern=NULL;
  unsigned int i, t, errflg = 0;
  int c;
  string err;
  
  program_name= argv[0];
  // parse the command line
  if (argc!=8) {
    errflg++;
  }
  else {
    FileIn = argv[2];
    Pattern = argv[1]; 
    // Read the atom position
    for(i=3, t=0; i<7; i++, t++) {
      c = sscanf(argv[i], "%u", &smartor[t]);
      if (c != 1) {
	break;
	errflg++;
      }      
    }
    c = sscanf(argv[7], "%f", &angle);
  }
  
  if (errflg)
    {
      err = "Usage: ";
      err += program_name;
      err += " \"PATTERN\" <filename> <atom1> <atom2> <atom3> <atom4> <angle> \n";
      ThrowError(err);
      exit(-1);
    }

  
  // create pattern
  OBSmartsPattern sp;
  sp.Init(Pattern);
  if (sp.NumAtoms() < 4) {
    err = program_name;
    err += ": the number of atoms in the SMART pattern must be higher than 3\n";
    ThrowError(err);
    exit(-1);
  }

  // Check the tortional atoms
  // actually the atom can be bonded but not written consecutively in the SMART
  // but this is the simplest way to ensure bonding
  if( (smartor[1] +1) != smartor[2] ) {
    err = program_name;
    err += ": the atoms of the rotating bond must be bonded\n";
    ThrowError(err);
    exit(-1);
  }
  for (i=0; i<4; i++) {
    if ( smartor[i] < 1 || smartor[i] > sp.NumAtoms()) {
      cerr << program_name 
	   << ": The tortional atom values must be between 1 and "
	   <<  sp.NumAtoms()
	   << ", which is the number of atoms in the SMART pattern.\n";
      exit(-1);
    }
  }

  // set the angle array
  ang_array[0] = 0;
  while (angle < 0.0f)   angle += 360.0f;
  while (angle > 360.0f) angle -= 360.0f;
  ang_array[1] = (unsigned char)rint(angle*res);

  // Find Input filetype
  if (extab.CanReadExtension(FileIn))
    inFileType = extab.FilenameToType(FileIn);
  else
    {
      cerr << program_name << ": cannot read input format!" << endl;
      exit (-1);
    }
  if (extab.CanWriteExtension(FileIn))
    outFileType = extab.FilenameToType(FileIn);
  else
    {
      cerr << program_name << ": cannot write input format!" << endl;
      exit (-1);
    }

  //Open the molecule file
  ifstream ifs;
  
  // Read the file
  ifs.open(FileIn);
  if (!ifs)
    {
      cerr << program_name << ": cannot read input file!" << endl;
      exit (-1);
    }

  
  OBMol mol(inFileType,outFileType);
  OBRotamerList rlist;
  vector< vector <int> > maplist;      // list of matched atoms
  vector< vector <int> >::iterator m;  // and its iterators
  int tindex;

  // Set the angles
  for (;;)
    {
      mol.Clear();
      ifs >> mol;                   // Read molecule
      if (mol.Empty()) break;

      if (sp.Match(mol)) {          // if match perform rotation
	
	maplist = sp.GetUMapList(); // get unique matches
		
	// look at all the mapping atom but save only the last one. 	
	for (m = maplist.begin(); m != maplist.end(); m++)
	  {
	    for (i=0; i<4 ; i++) {
	      tindex = (int) smartor[i] - 1; // get the tortional atom number
	      tor[i] = (unsigned char) (*m)[tindex];   // save it
	    }
	  }
	
	// Set the coordinates of references for rotation
	rlist.SetBaseCoordinateSets(mol);
	
	// Set the atoms to rotate
	rlist.Setup(mol,tor,nrotors);

	// Set the tortional angle value
	rlist.AddRotamers(ang_array,nrotamers);
      
	// Rotate and save the conformers
	rlist.ExpandConformerList(mol,mol.GetConformers());
	
	//change the molecule conformation
	mol.SetConformer(0);
      }
      cout << mol;
    }

  return(0);
} 
