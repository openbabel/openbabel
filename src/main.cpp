/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "molvector.h"
#include "obutil.h"
#include "parsmart.h"
#include "typer.h"
#include "rotor.h"
#include "binary.h"
#include "babelconfig.h"
#include "data.h"

#include <stdio.h>
#include <iostream>
#include <fstream>


#ifdef __BORLANDC__ 
 /* Borland c++ compiler does not have strncasecmp() */
 #include <string.h>
 #include <ctype.h>
 
 #ifdef DMALLOC
 #include "dmalloc.h"
 #endif
 
 int 
 strncasecmp(char *s1, char *s2, size_t n)
 {
         if (n == 0)
                 return 0;
 
         while ((n-- != 0)
             && (tolower(*(unsigned char *)s1) == tolower(*(unsigned char *)s2)))
         {
                 if (n == 0 || *s1 == '\0' || *s2 == '\0')
                         return 0;
                 s1++;
                 s2++;
         }
 
         return tolower(*(unsigned char *) s1) - tolower(*(unsigned char *) s2);
 }
 
#endif /* __BORLANDC__ */


using namespace std;
using namespace OpenBabel;

void usage();
void help();

// There isn't a great way to do this -- we need to save argv[0] for usage()
static char *program_name;

int main(int argc,char *argv[])
{
  io_type inFileType = UNDEFINED, outFileType = UNDEFINED;
  bool gotInType = false, gotOutType = false, removeHydrogens = false;
  bool addHydrogens = false, usePH = false, centerCoords = false;

  int arg, inFileArg, outFileArg;
  unsigned int firstMol = 1;
  unsigned int lastMol = INT_MAX;
  char *ext;
  char *formatOptions=NULL;
  OBFileFormat fileFormat;

  // Parse commandline
  program_name = argv[0];
  inFileArg = 0;
  outFileArg = 0;
  for (arg = 1; arg <= argc; arg++)
    {
      if (argv[arg])
	{
	  if (argv[arg][0] == '-')
	    {	      
	      switch (argv[arg][1])
		{
		case 'v':
		  {
		    cout << "Open Babel " << BABEL_VERSION << " -- " 
			 << __DATE__ << " -- " << __TIME__ << endl;
		    exit(0);
		  }
		case 'd':
		  removeHydrogens = true;
		  break;

		case 'f':
		  firstMol = atoi(argv[++arg]);
		  break;
		  
		case 'l':
		  lastMol = atoi(argv[++arg]);
		  break;

		case 'h':
		  addHydrogens = true;
		    if (strncmp(argv[arg],"-hpH",4) == 0) { usePH = true; }
		  break;

		case 'c':
		  centerCoords = true;
		  break;

		case 'i':
		  gotInType = true;
		  ext = argv[arg] + 2;
		  if (strncasecmp(ext, "MIME", 4) == 0)
		    outFileType = extab.MIMEToType(ext);
		  else if (extab.CanReadExtension(ext))
		    inFileType = extab.FilenameToType(ext);
		  else
		    {
		      cerr << program_name << ": cannot read input format!" << endl;
		      usage();
		    }
		  break;

		case 'o':
		  gotOutType = true;
		  ext = argv[arg] + 2;
		  if (strncasecmp(ext, "MIME", 4) == 0)
		    outFileType = extab.MIMEToType(ext);
		  else if (extab.CanWriteExtension(ext))
		    outFileType = extab.FilenameToType(ext);
		  else
		    {
		      cerr << program_name << ": cannot write output format!" << endl;
		      usage();
		    }

		  break;
		  
		case 'x':
		  formatOptions = argv[arg];
		  break;

		case 'H':
		  help();
		  exit(0);

		case '-':
		  if (inFileArg == 0)
		    inFileArg = -1;
		  else
		    outFileArg = -1;
		  break;
		  
		default:
		  cerr << program_name << ": unrecognized option `-" 
		       << argv[arg][1] << "'" << endl;
		  usage();
		  break;
		}
	    }
	  else if (inFileArg == 0)
	    inFileArg = arg;
	  else
	    outFileArg = arg;
	}
    }
  if (inFileArg == 0 || outFileArg == 0
      || (inFileArg < 0 && !gotInType)
      || (outFileArg < 0 && !gotOutType))
    usage();

  if (!gotInType)
    {
      if (extab.CanReadExtension(argv[inFileArg]))
	inFileType = extab.FilenameToType(argv[inFileArg]);
      else
	{
	  cerr << program_name << ": cannot read input format!" << endl;
	  usage();
	}
    }
  if (!gotOutType)
    {
      if (extab.CanWriteExtension(argv[outFileArg]))
	outFileType = extab.FilenameToType(argv[outFileArg]);
      else
	{
	  cerr << program_name << ": cannot write output format!" << endl;
	  usage();
	}
    }

  // Finally, we can do some work!
  OBMolVector moleculeList;
  ifstream inFileStream;
  bool usingStdin = false;
  bool canRead = true;
  int currentMol = 1;

  // read
  if (inFileArg > 0)
    {
      inFileStream.open(argv[inFileArg]);
      if (!inFileStream)
	{
	  cerr << program_name << ": cannot read input file!" << endl;
	  exit (-1);
	}
    }
  else
    usingStdin = true;

  while (canRead)
    {
      OBMol *mol = new OBMol(inFileType, outFileType);
      if (!usingStdin)
	fileFormat.ReadMolecule(inFileStream, *mol, argv[inFileArg]);
      else
	fileFormat.ReadMolecule(cin, *mol, "STDIN");

      if (mol->NumAtoms() != 0)
	{
	  // Perform any requested transformations
	  if (removeHydrogens)
	    mol->DeleteHydrogens();
	  if (addHydrogens)
	    mol->AddHydrogens(false, usePH);
	  if (centerCoords)
	    mol->Center();
      
	  if (currentMol >= firstMol && currentMol <= lastMol)
	    moleculeList.PushMol(mol);

	  if (!usingStdin && 
	      (inFileStream.peek() == EOF || !inFileStream.good()) )
	    canRead = false;
	  else if (usingStdin && (cin.peek() == EOF || !cin.good()) )
	    canRead = false;
	  else if (currentMol > lastMol)
	    canRead = false;
	}
      else // 0 atoms in this molecule! 
	{
	  //	  cerr << " error: has zero atoms! " << endl;
	  canRead = false;
	}

      currentMol++;
    }
 
  // write
  if (outFileArg > 0)
    {
      ofstream outFileStream(argv[outFileArg]);
      if (!outFileStream)
	{
	  cerr << program_name << ": cannot write to output file!" << endl;
	  exit (-1);
	}
      moleculeList.Write(outFileStream, formatOptions);
    }
  else
    moleculeList.Write(cout, formatOptions);

  return(0);
}

void usage()
{
  cout << "Open Babel " << BABEL_VERSION << " -- " << __DATE__ << " -- "
       << __TIME__ << endl;
  cout << "Usage: " << program_name
       << " [-i<input-type>] <name> [-o<output-type>] <name>" << endl;
  cout << "Try `babel -H' for more information." << endl;

  exit (0);
}

// Make sure to update the Doxygen-generated man documentation below! -GH
void help()
{
unsigned int i;

  cout << "Open Babel " << BABEL_VERSION << " -- " << __DATE__ << " -- "
       << __TIME__ << endl;
  cout << "Usage: " << program_name 
       << " [-i<input-type>] <name> [-o<output-type>] <name>" << endl;
  cout << endl << "Currently supported input types" << endl;
  for (i = 0; i < extab.Count(); i++)
    if (extab.IsReadable(i))
      cout << "\t" << extab.GetExtension(i) << " -- " 
		<< extab.GetDescription(i) << " file" << endl;
  cout << endl << "Currently supported output types" << endl;
  for (i = 0; i < extab.Count(); i++)
    if (extab.IsWritable(i))
      cout << "\t" << extab.GetExtension(i) << " -- " 
		<< extab.GetDescription(i) << " file" << endl;
  cout << "Additional options : " << endl;
  cout << " -f <#> Start import at molecule # specified " << endl;
  cout << " -l <#> End import at molecule # specified " << endl;
  cout << " -d Delete Hydrogens " << endl;
  cout << " -h Add Hydrogens " << endl;
  cout << " -hpH Add Hydrogens appropriate for pH (use transforms in phmodel.txt) " << endl; 
  cout << " -c Center Coordinates " << endl;
  cout << " -x[flags] XML.CML options (e.g. -x1ac)  " << endl;
  cout << "   1 output CML V1.0 (default)" << endl;
  cout << "   2 output CML V2.0 (Schema)" << endl;
  cout << "   a output array format for atoms and bonds (default <atom>)" << endl;
  cout << "   p prettyprint output (not implemented)" << endl;
  cout << "   n output namespace (default no namespace)" << endl;
  cout << "   c use 'cml' as output namespace prefix (else default) (forces n)" << endl;
  cout << "   d output DOCTYPE (default none)" << endl;
  cout << "   g debug output" << endl;
  cout << "   v add XML version (declaration)" << endl;
  cout << endl << "Report Bugs to <openbabel-discuss@lists.sourceforge.net>." << endl;

}

/* OpenBabel man page*/
/** \page openbabel a converter for molecular modeling data files
*
* \n
* \par SYNOPSIS
*
* \b openbabel [-d] [-h] [-hpH] [-f<#> -l<#>] [-c] [-x<xmlflags>] [-i<input-type>] <infile> [-o<output-type>] <outfile>  
*
* \par DESCRIPTION
*
* Open Babel is a program designed to interconvert a number of 
* file formats currently used in molecular modeling software. \n\n
* Note that OpenBabel can also be used as a library to import
* molecular formats in other software packages. See the OpenBabel
* web pages (http://openbabel.sourceforge.net) for more information.
*
* \par OPTIONS
*
* If only input and ouput files are given, Open Babel will guess 
* the file type from the filename extension. \n\n
* \b -i :
*     Specifies input format, see below for the available formats \n\n
* \b -o :
*     Specifies output format, see below for the available formats \n\n
* \b -d : 
*     Delete Hydrogens \n\n
* \b -h : 
*     Add Hydrogens \n\n
* \b -hpH : 
*     Add Hydrogens appropriate for pH (use transforms in phmodel.txt)  \n\n
* \b -f<#> : 
*     For multiple entries input, start import at molecule # \n\n
* \b -l<#> : 
*     For multiple entries input, stop import at molecule # \n\n
* \b -c : 
*     Center atomic coordinates at (0,0,0) \n\n
* \b -x<flags>:
*     XML.CML options (e.g. -x1ac) (see below) \n\n
*
*  \par INPUT FORMATS
*
*  Open Babel recognizes the following input formats:
*   \n     alc -- Alchemy file
*   \n     prep -- AMBER PREP file
*   \n     bs -- Ball and Stick file
*   \n     caccrt -- Cacao Cartesian file
*   \n     ccc -- CCC file
*   \n     box -- Dock 3.5 Box file
*   \n     feat -- Feature file
*   \n     gam -- GAMESS Output file
*   \n     gamout -- GAMESS Output file
*   \n     mm1gp -- Ghemical MM file
*   \n     qm1gp -- Ghemical QM file
*   \n     hin -- Hyperchem HIN file
*   \n     jout -- Jaguar Output file
*   \n     bin -- OpenEye Binary file
*   \n     mmd -- MacroModel file
*   \n     mmod -- MacroModel file
*   \n     out -- MacroModel file
*   \n     dat -- MacroModel file
*   \n     car -- MSI Biosym/Insight II CAR file
*   \n     sdf -- MDL Isis SDF file
*   \n     sd -- MDL Isis SDF file
*   \n     mdl -- MDL Molfile file
*   \n     mol -- MDL Molfile file
*   \n     mopcrt -- MOPAC Cartesian file
*   \n     mopout -- MOPAC Output file
*   \n     mpqc -- MPQC file
*   \n     bgf -- MSI BGF file
*   \n     nwo -- NWChem Output file
*   \n     pdb -- PDB file
*   \n     qcout -- QChem Output file
*   \n     smi -- SMILES file
*   \n     mol2 -- Sybyl Mol2 file
*   \n     unixyz -- UniChem XYZ file
*   \n     xyz -- XYZ file
*
* \par OUTPUT FORMATS
*
*  The following output formats may be written by Open Babel:
*   \n     alc -- Alchemy file
*   \n     bs -- Ball & Stick file
*   \n     caccrt -- Cacao Cartesian file
*   \n     cacint -- Cacao Internal file
*   \n     cache -- CAChe MolStruct file
*   \n     ct -- ChemDraw Connection Table file
*   \n     cssr -- CSD CSSR file
*   \n     box -- Dock 3.5 Box file
*   \n     dmol -- DMol3 Coordinates file
*   \n     feat -- Feature file
*   \n     fh -- Fenske-Hall Z-Matrix file
*   \n     gamin -- GAMESS Input file
*   \n     inp -- GAMESS Input file
*   \n     gcart -- Gaussian Cartesian file
*   \n     gau -- Gaussian Input file
*   \n     mm1gp -- Ghemical MM file
*   \n     gr96A -- GROMOS96 (A) file
*   \n     gr96N -- GROMOS96 (nm) file
*   \n     hin -- HyperChem HIN file
*   \n     jin -- Jaguar Input file
*   \n     bin -- OpenEye Binary file
*   \n     mmd -- MacroModel file
*   \n     mmod -- MacroModel file
*   \n     out -- MacroModel file
*   \n     dat -- MacroModel file
*   \n     sdf -- MDL Isis SDF file
*   \n     sd -- MDL Isis SDF file
*   \n     mdl -- MDL Molfile file
*   \n     mol -- MDL Molfile file
*   \n     mopcrt -- MOPAC Cartesian file
*   \n     bgf -- MSI BGF file
*   \n     csr -- MSI Quanta CSR file
*   \n     nw -- NWChem Input file
*   \n     pdb -- PDB file
*   \n     report -- Report file
*   \n     qcin -- QChem Input file
*   \n     smi -- SMILES file
*   \n     fix -- SMILES Fix file
*   \n     mol2 -- Sybyl Mol2 file
*   \n     txyz -- Tinker XYZ file
*   \n     txt -- Titles file
*   \n     unixyz -- UniChem XYZ file
*   \n     xed -- XED file
*   \n     xyz -- XYZ file
*
* \par XML.CML OPTIONS
*    Chemical Markup Language (http://www.xml-cml.org/) options \n \n
*     \n 1 output CML V1.0 (default)
*     \n 2 output CML V2.0 (Schema)
*     \n a output array format for atoms and bonds (default <atom>)
*     \n p prettyprint output (default no indent)
*     \n n output namespace (default no namespace)
*     \n c use 'cml' as output namespace prefix (else default) (forces n)
*     \n d output DOCTYPE (default none)
*     \n g debug output
*     
* \par EXAMPLE
*     openbabel -ixyz ethanol.xyz -opdb ethanol.pdb
* \par AUTHOR
*
* Open Babel is derived from OElib, written by \b Matt \b Stahl, which is a rewrite of the classic babel program. Open Babel is currently maintained by \b Geoff \b Hutchison and \b Michael \b Banck.
*
* \par COPYRIGHT
*  Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
*  Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison \n \n
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation version 2 of the License.\n \n
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
* \par SEE ALSO
*   The web pages for OpenBabel can be found at http://openbabel.sourceforge.net
**/
