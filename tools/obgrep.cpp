/**********************************************************************
obgrep = Open Babel molecule grep
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
#include "obutil.h"
#include "parsmart.h"
#include "obifstream.h"
#include <unistd.h>

using namespace std;
using namespace OpenBabel;

///////////////////////////////////////////////////////////////////////////////
//! \brief Find the molecule(s) with or without a given SMART pattern
int main(int argc,char **argv)
{
  char *program_name=NULL;
  io_type inFileType = UNDEFINED, outFileType = UNDEFINED;
  int c; 
  bool count=false, invert=false, full=false;
  char *FileIn = NULL, *Pattern = NULL;
  
  // Parse options
  while ((c = getopt(argc, argv, "vcf")) != -1)
    switch (c) {
    case 'c': // count the number of match
              count = true;
	      break;
    case 'v': // match only the molecules without the pattern 
              invert = true;
	      break;

    case 'f': full = true;
              break;

    case '?':
              if (isprint (optopt))
		fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	      else
		fprintf (stderr,
			 "Unknown option character `\\x%x'.\n",
			 optopt);
	      return 1;
    }

  int index = optind;
  program_name = argv[0];

  if (argc-index != 2)
    {
      string err = "Usage: ";
      err += program_name;
      err += " [options] \"PATTERN\" <filename>\n";
      err += "Options:\n";
      err += "   -v    Invert the matching, print non-matching molecules\n";
      err += "   -c    Print the number of matched molecules\n";
      err += "   -f    Full match, print matching-molecules when the number\n";
      err += "         of heavy atoms is equal to the number of PATTERN atoms\n";
      ThrowError(err);
      exit(0);
    }
  else {
    Pattern = argv[index++];
    FileIn  = argv[index]; 
  }

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


  // Match the SMART
  OBSmartsPattern sp;
  sp.Init(Pattern);
  ifstream ifs;
  
  // Read the file
  ifs.open(FileIn);
  if (!ifs)
    {
      cerr << program_name << ": cannot read input file!" << endl;
      exit (-1);
    }

  OBMol mol(inFileType, outFileType);
  
  bool impossible_match;

  // Search for pattern
  for (c=0;;)
    {
      mol.Clear();
      ifs >> mol;
      if (mol.Empty()) break;
      
      
      // impossible to make a full match if the number of atoms is
      // different
      if (full ) 
	impossible_match = (sp.NumAtoms() == mol.NumHvyAtoms()) ? false : true;
      else
	impossible_match = false;

      if (impossible_match) { // avoid useless SMART matching attempt 
	if (invert) {
	  if (!count) cout << mol;
	  c++;
	}
	continue;
      }
	
      if (sp.Match(mol) ) { // perform SMART matching
	if (!invert) {
	  if (!count) cout << mol;
	  c++;
	}
      }
      else {
	if (invert) {
	  if (!count) cout << mol;
	  c++;
	}
      }
    }

  if (count)
    cout << c << endl;

  return(1);
} 
