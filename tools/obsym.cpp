/**********************************************************************
obsym - 3D Point Group Symmetry

Copyright (C) 2007 Geoffrey R. Hutchison
Based on code (C) 1996,2003 by S. Patchkovskii
 
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
#include <openbabel/pointgroup.h>
#ifndef _MSC_VER
  #include <unistd.h>
#endif

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  char *FileIn = NULL;

  if (argc != 2) {
    cerr << " Usage: " << program_name << " <input file>\n";
    exit(-1);
  }
  else {
      FileIn  = argv[1];
  }

  // Find Input filetype
  OBConversion conv;
  OBFormat *inFormat = conv.FormatFromExt(FileIn);
    
  if (!inFormat || !conv.SetInFormat(inFormat)) {
    cerr << program_name << ": cannot read input format!" << endl;
    exit (-1);
  }

  ifstream ifs;

  // Read the file
  ifs.open(FileIn);
  if (!ifs) {
    cerr << program_name << ": cannot read input file!" << endl;
    exit (-1);
  }

  OBMol mol;
  OBPointGroup pg;

  for (c = 1;; ++c)
    {
      mol.Clear();
      conv.Read(&mol, &ifs);
      if (mol.Empty())
        break;
      
      // not needed by OBPointGroup, but useful for external programs
      mol.Center();
      mol.ToInertialFrame();

      pg.Setup(&mol);
      cout << "Point Group: " << pg.IdentifyPointGroup() << endl;

    } // end for loop
  
  return(1);
}



///////////////////////////////////////////////////////////////////////////////
//! \return the number of size of the set of smallest rings (SSSR)
int nrings(OBMol &mol)
{
  int nr;
  vector<OBRing*> vr;
  
  vr = mol.GetSSSR();
  nr = vr.size();
  return (nr);
}

//! \return the sequence of residues ordered by chain
string sequence(OBMol &mol)
{
  unsigned int currentChain = 0;
  string residueSequence;
  FOR_RESIDUES_OF_MOL(r, mol)
    {
      if (r->GetName().find("HOH") != string::npos)
        continue;
      
      if (r->GetChainNum() != currentChain) {
        if (residueSequence.size() != 0) { // remove the trailing "-"
          residueSequence.erase(residueSequence.size() - 1);
          residueSequence += ", "; // separate different chains
        }
        
        currentChain = r->GetChainNum();
      }
      residueSequence += r->GetName();
      residueSequence += "-";
    }
  if (residueSequence.size() != 0) // remove the trailing "-"
    residueSequence.erase(residueSequence.size() - 1);

  return residueSequence; 
}

/* obprop man page*/
/** \page obprop print standard molecular properties
*
* \n
* \par SYNOPSIS
*
* \b obprop \<filename\>
*
* \par DESCRIPTION
*
* The obprop program is a tool to print a set of standard molecular properties
* for all molecules in a file. It also serves as example code for using the
* Open Babel library (libopenbabel).
* 
* \par EXAMPLES
*
*   obprop pyridines.sdf
*
* \par AUTHORS
*
* The obprop program was contributed by \b Fabien \b Fontaine.
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
**/
