/**********************************************************************
obfit = Fit molecules according to a SMART pattern
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
  Require a fixed molecule, a set of molecules to move, 
  a SMART pattern to match the fixed and the moving molecules
  If the SMART is not found in the fixed molecule the program exits
  If the SMART is not found in a moving molecule, the molecule is not moved 
  example of command line:
  obfit "[nh]1c2c(=O)n(C)c(=O)n(C)c2cc1" testref.sdf testmv.sdf 
*/

#include "mol.h"
#include "parsmart.h"
#include <unistd.h>


using namespace std;
using namespace OpenBabel;


//find the center of mass of a list of atoms
vector3 mass_c( vector<int> &aindex, OBMol &mol);

///////////////////////////////////////////////////////////////////////////////
//! \brief superimpose a set of molecules on the atoms of a reference molecule
//! The atoms used for the overlay are defined by the SMART pattern 
int main(int argc,char **argv)
{

  int errflg=0;
  char *FileRef=NULL, *FileMove=NULL, *Pattern=NULL;
  string err;
  char *program_name=argv[0];
  io_type refFileType = UNDEFINED, mvFileType = UNDEFINED;

  // parse the command line
  if (argc!=4) {
    errflg++;
  }
  else {
    FileRef = argv[2];
    FileMove = argv[3];
    Pattern = argv[1]; 
  }

  if (errflg)
    {
      err = "Usage: ";
      err += program_name;
      err += " \"PATTERN\" <fixed_structure> <moving_structures>\n";
      ThrowError(err);
      exit(-1);
    }


  // create the pattern
  OBSmartsPattern sp;
  if (!sp.Init(Pattern)) {
    err = program_name;
    err += ": Unable to read the SMART: ";
    err += Pattern;
    ThrowError(err);
    exit(-1);
  }

  // Find Input filetypes
  if (extab.CanReadExtension(FileRef))
    refFileType = extab.FilenameToType(FileRef);
  else
    {
      cerr << program_name << ": cannot read fixed molecule format!" << endl;
      exit (-1);
    }
  
  if (extab.CanReadExtension(FileMove))
    mvFileType = extab.FilenameToType(FileMove);
  else
    {
      cerr << program_name << ": cannot read moving molecule(s) format!" << endl;
      exit (-1);
    }
  if (! extab.CanWriteExtension(FileMove))
    {
      cerr << program_name << ": cannot write moving molecule format!" << endl;
      exit (-1);
    }
  
  ifstream ifsref;
  OBMol molref(refFileType,mvFileType);
  vector< vector <int> > maplist;      // list of matched atoms
  vector< vector <int> >::iterator i;  // and its iterators
  vector< int >::iterator j;
  vector <int> refatoms;
  
  //Read the reference structure
  ifsref.open(FileRef);
  if (!ifsref)
    {
      cerr << program_name << ": cannot read fixed molecule file: " 
	   << FileRef << endl;
      exit (-1);
    }

  molref.Clear();
  ifsref >> molref; 

  // and check if the SMART match
  sp.Match(molref);
  maplist = sp.GetUMapList(); // get unique matches
  if (maplist.empty()) {
    err = program_name;
    err += ": Unable to map SMART: ";
    err += Pattern;
    err += " in reference molecule: ";
    err += FileRef;
    ThrowError(err);
    exit(-1);
  }

  // Find the matching atoms
  for (i = maplist.begin(); i != maplist.end(); i++) 
    {
      refatoms.clear(); // Save only the last set of atoms
      for(j= (*i).begin(); j != (*i).end(); j++)
	{
	  refatoms.push_back(*j);
	}
    }
  
  // set the translation vector
  vector3 tvref(0,0,0);
  OBAtom *atom;
  unsigned int c;

  tvref = mass_c(refatoms, molref);
  // center the molecule 
  molref.Translate(-tvref);
  
  //get the coordinates of the SMART atoms
  double *refcoor = new double[refatoms.size()*3];
  
  for(c=0; c<refatoms.size(); c++)
    {
      atom = molref.GetAtom(refatoms[c]);
      refcoor[c*3] = atom->x();
      refcoor[c*3+1] = atom->y();
      refcoor[c*3+2] = atom->z();
    }

  ifstream ifsmv;
  OBMol molmv(mvFileType,mvFileType);
  vector <int> mvatoms;
  vector3 tvmv;
  unsigned int size=0;
  double rmatrix[3][3];

  //Read the moving structures
  ifsmv.open(FileMove);
  if (!ifsmv)
    {
      cerr << program_name << ": cannot read file: " 
	   << FileMove << endl;
      exit (-1);
    }

  for (;;)
    {
      molmv.Clear();
      ifsmv >> molmv;                   // Read molecule
      if (molmv.Empty()) break;

      if (sp.Match(molmv)) {          // if match perform rotation
	
	maplist = sp.GetUMapList(); // get unique matches
	
	// Find the matching atoms
	for (i = maplist.begin(); i != maplist.end(); i++) 
	  {
	    mvatoms.clear(); // Save only the last set of atoms
	    for(j= (*i).begin(); j != (*i).end(); j++)
	      {
		mvatoms.push_back(*j);
	      }
	  }
  
	tvmv = mass_c(mvatoms, molmv);
	// center the molecule 
	molmv.Translate(-tvmv);
	
	//Find the rotation matrix
	size = mvatoms.size();
	if (size != refatoms.size()) {
	  err = program_name;
	  err += ": Error: not the same number of SMART atoms";
	  ThrowError(err);
	  exit(-1);
	}
	double *mvcoor = new double[size*3];

	for(c=0; c<size; c++)
	  {
	    atom = molmv.GetAtom(mvatoms[c]);
	    mvcoor[c*3] = atom->x();
	    mvcoor[c*3+1] = atom->y();
	    mvcoor[c*3+2] = atom->z();
	  }

	// quaternion fit
	qtrfit(refcoor, mvcoor, size, rmatrix);
		
	delete[] mvcoor;
	//rotate all the atoms
	molmv.Rotate(rmatrix);

	//translate the rotated molecule
	molmv.Translate(tvref);
      }
      cout << molmv;
    }

  delete[] refcoor;
  return(0);
} 

///////////////////////////////////////////////////////////////////////////////
//find the center of mass of a list of atoms
vector3 mass_c( vector<int> &aindex, OBMol &mol)
{
  vector3 center(0,0,0);
  vector< int >::iterator j;
  OBAtom *atom;

  for(j= aindex.begin(); j != aindex.end(); j++) 
    { 
      atom = mol.GetAtom(*j);
      center += atom->GetVector();
    }

  center /= (float) aindex.size();
  return (center);
}

