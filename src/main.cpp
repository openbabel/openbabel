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

// There isn't a great way to do this -- we need to save argv[0] for usage()
static char *program_name;

int main(int argc,char *argv[])
{
  io_type inFileType = UNDEFINED, outFileType = UNDEFINED;
  bool gotInType = false, gotOutType = false, removeHydrogens = false;
  bool addHydrogens = false, usePH = false, centerCoords = false;

  int arg, inFileArg, outFileArg;
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
		      cerr << program_name << ": Cannot read input format!" << endl;
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
		      cerr << program_name << ": Cannot write output format!" << endl;
		      usage();
		    }

		  break;
		  
		case 'x':
		  formatOptions = argv[arg];
		  break;

		case '-':
		  if (inFileArg == 0)
		    inFileArg = -1;
		  else
		    outFileArg = -1;
		  break;
		  
		default:
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
	  cerr << program_name << ": Cannot read input format!" << endl;
	  usage();
	}
    }
  if (!gotOutType)
    {
      if (extab.CanWriteExtension(argv[outFileArg]))
	outFileType = extab.FilenameToType(argv[outFileArg]);
      else
	{
	  cerr << program_name << ": Cannot write output format!" << endl;
	  usage();
	}
    }

  // Finally, we can do some work!
  OBMol mol(inFileType, outFileType);

  // read
  if (inFileArg > 0)
    {
      ifstream inFileStream(argv[inFileArg]);
      if (!inFileStream)
	{
	  cerr << program_name << ": Cannot read input file!" << endl;
	  exit (-1);
	}
      fileFormat.ReadMolecule(inFileStream, mol, argv[inFileArg]);
    }
  else
    fileFormat.ReadMolecule(cin, mol, "STDIN");

  // Perform any requested transformations
  if (removeHydrogens)
    mol.DeleteHydrogens();
  if (addHydrogens)
    mol.AddHydrogens(false, usePH);
  if (centerCoords)
    mol.Center();
 
  // Dimensionality of file for output
  char *dimension;
  if (mol.Has3D())
    dimension = "3D";
  else
    dimension = "2D";

  // write
  if (outFileArg > 0)
    {
      ofstream outFileStream(argv[outFileArg]);
      if (!outFileStream)
	{
	  cerr << program_name << ": Cannot write to output file!" << endl;
	  exit (-1);
	}
      fileFormat.WriteMolecule(outFileStream, mol, dimension, formatOptions);
    }
  else
    fileFormat.WriteMolecule(cout, mol, dimension, formatOptions);

  return(0);
}

void usage()
{
unsigned int i;

  cout << "Open Babel " << BABEL_VERSION << " -- " << __DATE__ << " -- "
       << __TIME__ << endl;
  cout << "Usage is : " << endl << program_name 
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
  cout << " -d Delete Hydrogens " << endl;
  cout << " -h Add Hydrogens " << endl;
  cout << " -hpH Add Hydrogens appropriate for pH (use transforms in phmodel.txt) " << endl; 
  cout << " -c Center Coordinates " << endl;
  cout << " -x[flags] XML.CML options (e.g. -x1ac)  " << endl;
  cout << "   1 output CML V1.0 (default)" << endl;
  cout << "   2 output CML V2.0 (Schema)" << endl;
  cout << "   a output array format for atoms and bonds (default <atom>)" << endl;
  cout << "   p prettyprint output (default no indent)" << endl;
  cout << "   n output namespace (default no namespace)" << endl;
  cout << "   c use 'cml' as output namespace prefix (else default) (forces n)" << endl;
  cout << "   d output DOCTYPE (default none)" << endl;
  cout << "   g debug output" << endl;

  exit(0);
}
