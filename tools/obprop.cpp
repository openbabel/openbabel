/**********************************************************************
obprop = Open Babel properties calculation
Copyright (C) 2003 Fabien Fontaine

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include <unistd.h>

using namespace std;
using namespace OpenBabel;

// PROTOTYPES /////////////////////////////////////////////////////////////////
int nrings(OBMol &mol);


///////////////////////////////////////////////////////////////////////////////
//! \brief Compute some properties easy to access from open babel
//


int main(int argc,char **argv)
{
  char *program_name= argv[0];
  io_type inFileType = UNDEFINED;
  int c;
  char *FileIn = NULL;

 
  if (argc != 2)
    {
      string err = "Usage: ";
      err += program_name;
      err += " <filename>\n";
      err += "Output format:\n";
      err += "name NAME\n";
      err += "mol_weight MOLECULAR_WEIGHT\n";
      err += "num_rings NUMBER_OF_RING_(SSSR)\n";
      err += "$$$$";
      ThrowError(err);
      exit(-1);
    }
  else {
    FileIn  = argv[1]; 
  }

  // Find Input filetype
  if (extab.CanReadExtension(FileIn))
    inFileType = extab.FilenameToType(FileIn);
  else
    {
      cerr << program_name << ": cannot read input format!" << endl;
      exit (-1);
    }

  ifstream ifs;
  
  // Read the file
  ifs.open(FileIn);
  if (!ifs)
    {
      cerr << program_name << ": cannot read input file!" << endl;
      exit (-1);
    }

  OBMol mol(inFileType, UNDEFINED);
  

  ////////////////////////////////////////////////////////////////////////////
  // List of properties
  // Name
  // Molecular weight (Standard molar mass given by IUPAC atomic masses)
  // Number of rings : the size of the smallest set of smallest rings (SSSR)

  //.....ADD YOURS HERE.....

  for (c=0;;)
    {
      mol.Clear();
      ifs >> mol;
      if (mol.Empty()) break;

      // Print the properties
      // The name should be enough self explaining

      cout << "name " << mol.GetTitle() << endl;
      cout << "mol_weight "<< mol.GetMolWt() << endl;
      cout << "num_rings " << nrings(mol) << endl;
      cout << "$$$$" << endl; // SDF like end of compound descriptor list

    } // end for loop

  return(1);
} 



///////////////////////////////////////////////////////////////////////////////
//! \brief Return the number of size of the set of smallest rings (SSSR)
int nrings(OBMol &mol)
{
  int nr;
  vector<OBRing*> vr;

  vr = mol.GetSSSR();
  nr = vr.size();

  return (nr);
}
