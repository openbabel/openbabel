/**********************************************************************
molvector.cpp - Vector to handle set of molecules like in multiple 
		files.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

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

#include "mol.h"
#include "obutil.h"

#include "molvector.h"

using namespace std;

namespace OpenBabel {

/** \class OBMolVector
    \brief Molecule Group Class

Functions for dealing with groups of molecules.  OBMolVector will read
either all molecules from a file or a set of conformers and can write
a multi-structure file.

Molecules are indexed as a normal C++ vector -- that is the first
molecule starts at index 0.
*/

OBMolVector::~OBMolVector()
{
  for (unsigned int i = 0; i < _molvec.size(); i++)
    {
      delete _molvec[i];
      _molvec[i] = NULL;
    }
}

//! Read molecules from a file into a OBMolVector. Input and output types
//! default to SDF. Can specify one particular molecule to read, or
//! -1 for all atoms.
void OBMolVector::Read(ifstream &ifs, const io_type in_type, const io_type out_type, int nToRead)
{	
  int nRead= 0;
  OBFileFormat ff;
  while (1)
    {
      if (nRead == nToRead) break;
      OBMol *mol;
      mol = new OBMol;
      (*mol).SetInputType(in_type);
      (*mol).SetOutputType(out_type);
      ff.ReadMolecule(ifs,*mol);
      nRead++;
      if (!(*mol).NumAtoms())
	{
	  delete mol;
	  mol = NULL;
	  break;
	}
      _molvec.push_back(mol);
    }
}

//! Write a OBMolVector to a file.  Output type defaults to SDF
void OBMolVector::Write(ostream &os, const char *options)
{
  vector<OBMol *>::iterator mol_i;
  OBFileFormat ff;
  char *dimension;

  //dimension and options
  for (mol_i = _molvec.begin(); mol_i != _molvec.end(); mol_i++)
    {
      if ((*mol_i)->Has3D())
	dimension = "3D";
      else
	dimension = "2D";

      ff.WriteMolecule(os,(**mol_i), dimension, options);
    }
}

//! Get a specific molecule from a OBMolVector.  Index starts at zero.
OBMol *OBMolVector::GetMol(unsigned int i)
{
  if (i >= 0 && i < _molvec.size())
    return(_molvec[i]);
  else
    {
      cerr << "Index " << i << " out of range in OBMolVector::GetMol " << endl;
      return(NULL);
    }
}

//! Read a set of conformers from an input file and put them into an OBMolVector.
//! This function read the first molecule and sets the current title (held
//! int the variable master) to be the current title.  It continues to read
//! molecules and push them into the vector util it reads a molecule with a
//! different name.  At this point it rewinds the file stream to the beginning
//! of the current molecule and returns
bool OBMolVector::ReadConfs(ifstream &ifs, const io_type in_type, const io_type out_type)
{
  OBMol *mol;
  OBFileFormat ff;
  string title,master;

  _molvec.resize(0);
  
  int i = 1;
  while (1)
    {
      mol = new OBMol;
      (*mol).SetInputType(in_type);
      (*mol).SetOutputType(out_type);
      streampos sp = ifs.tellg();
      ff.ReadMolecule(ifs,*mol);
      if (mol->NumAtoms() == 0)
        {
          delete mol; mol = NULL;
	  return(false);
        }
      
      title = mol->GetTitle();
      if (i == 1)
        {
          master = title;
          _molvec.push_back(mol);
        }
      else
        {
          if (title == master)
	    _molvec.push_back(mol);
          else
            {
              ifs.seekg(sp);
              delete mol; mol = NULL;
              break;
            }
        }
      i++;
    }
  return(true);
}


//! Push an OBMol onto the end of the OBMolVector
void OBMolVector::PushMol(OBMol *mol)
{
  _molvec.push_back(mol);
}

} // namespace OpenBabel
