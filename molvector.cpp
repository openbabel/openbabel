/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "oeutil.h"

#include "molvector.h"

namespace OpenEye {

//Functions for dealing with groups of molecules.  MolVec will read either all
//molecules from a file or a set of conformers.  

OEMolVector::~OEMolVector()
{
  for (unsigned int i = 0; i < _molvec.size(); i++)
    {
      delete _molvec[i];
    }
}

// Read all molecules from a file into a OEMolVector.  Input and output types
// default to SDF

void OEMolVector::Read(ifstream &ifs, const io_type in_type, const io_type out_type, int nToRead)
{	
  int nRead= 0;
  OEFileFormat ff;
  while (1)
    {
      if (nRead == nToRead) break;
      OEMol *mol;
      mol = new OEMol;
      (*mol).SetInputType(in_type);
      (*mol).SetOutputType(out_type);
      ff.ReadMolecule(ifs,*mol);
      nRead++;
      if (!(*mol).NumAtoms())
	{
	  delete mol;
	  break;
	}
      _molvec.push_back(mol);
    }
}

// Write a OEMolVector to a file.  Output type defaults to SDF

void OEMolVector::Write(ofstream &ofs)
{
  vector<OEMol *>::iterator mol_i;
  OEFileFormat ff;

  for (mol_i = _molvec.begin(); mol_i != _molvec.end(); mol_i++)
    {
      ff.WriteMolecule(ofs,(**mol_i));
    }
}

// Get a specific molecule from a OEMolVector.  Index starts at zero.
OEMol *OEMolVector::GetMol(int i)
{
  if (i >= 0 && i < (signed)_molvec.size())
    return(_molvec[i]);
  else
    {
      cerr << "Index " << i << " out of range in OEMolVector::GetMol " << endl;
      return(NULL);
    }
}

// Read a set of conformers from an input file and put them into a MolVec.
// This function read the first molecule and sets the current title (held
// int the variable master) to be the current title.  It continues to read
// molecules and push them into the vector util it reads a molecule with a
// different name.  At this point it rewinds the file stream to the beginning
// of the current molecule and returns
bool OEMolVector::ReadConfs(ifstream &ifs, const io_type in_type, const io_type out_type)
{
  OEMol *mol;
  OEFileFormat ff;
  string title,master;

  _molvec.resize(0);
  
  int i = 1;
  while (1)
    {
      mol = new OEMol;
      (*mol).SetInputType(in_type);
      (*mol).SetOutputType(out_type);
      streampos sp = ifs.tellg();
      ff.ReadMolecule(ifs,*mol);
      if (mol->NumAtoms() == 0)
        {
          delete mol;
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
              delete mol;
              break;
            }
        }
      i++;
    }
  return(true);
}


} // namespace OpenEye








