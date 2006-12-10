/**********************************************************************
obprop - Open Babel properties calculation

Copyright (C) 2003 Fabien Fontaine
Some portions Copyright (C) 2004-2005 Geoffrey R. Hutchison
 
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
#include <unistd.h>

using namespace std;
using namespace OpenBabel;
void some_tests(OBMol &mol);
// PROTOTYPES /////////////////////////////////////////////////////////////////
int nrings(OBMol &mol);


///////////////////////////////////////////////////////////////////////////////
//! \brief Compute some properties easy to access from open babel
//


int main(int argc,char **argv)
{
  char *program_name= argv[0];
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
  else
    {
      FileIn  = argv[1];
    }

  // Find Input filetype
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(FileIn);
    
  if (!format || !conv.SetInAndOutFormats(format, format))
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
  
  OBMol mol;
  
  
  ////////////////////////////////////////////////////////////////////////////
  // List of properties
  // Name
  // Molecular weight (Standard molar mass given by IUPAC atomic masses)
  // Number of rings : the size of the smallest set of smallest rings (SSSR)
  
  //.....ADD YOURS HERE.....
  
  for (c=0;;)
    {
      mol.Clear();
      conv.Read(&mol, &ifs);
      if (mol.Empty())
        break;
      // Print the properties
      // The name should be enough self explaining
      some_tests(mol);
        
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

void some_tests(OBMol &mol)
{
  for ( OBResidueIterator residue = mol.BeginResidues();
        residue != mol.EndResidues(); ++residue )
    {
      std::cout << "residue named: " 
                << (*residue)->GetName() << " has atoms\n";
      for ( OBAtomIterator atom = (*residue)->BeginAtoms();
                             atom != (*residue)->EndAtoms(); ++atom )
        {
          std::cout << (*atom)->GetIdx() << "\n";
        }
    }
  for ( OBBondIterator bond = mol.BeginBonds();
        bond != mol.EndBonds(); ++bond )
    {
        std::cout << "bond: ";
        std::cout << (*bond)->GetIdx() << " connects (";
        std::cout << (*bond)->GetBeginAtom()->GetIdx() << ";";
        std::cout << (*bond)->GetEndAtom()->GetIdx() << ")\n";
    }
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
* For more contributors to Open Babel, see http://openbabel.sourceforge.net/THANKS.shtml
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
*   The web pages for Open Babel can be found at: http://openbabel.sourceforge.net/ \n
**/
