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
#include "parsmart.h"
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
  int ntimes=0;
  bool pattern_matched=false, ntimes_matched=true;
  bool count=false, invert=false, full=false, name_only=false;
  char *FileIn = NULL, *Pattern = NULL;
  program_name = argv[0];

  // Parse options
  while ((c = getopt(argc, argv, "t:nvcf")) != -1)
    switch (c) {
    case 't': // request ntimes unique matches
      
              c = sscanf(optarg, "%d", &ntimes);
	      if (c != 1 ) {
		cerr << program_name << ": unable to parse -t option" << endl;
		exit (-1);
	      }
	      break;
    
    case 'n': // print the molecule name only
              name_only = true;
	      break;
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
  

  if (argc-index != 2)
    {
      string err = "Usage: ";
      err += program_name;
      err += " [options] \"PATTERN\" <filename>\n";
      err += "Options:\n";
      err += "   -v      Invert the matching, print non-matching molecules\n";
      err += "   -c      Print the number of matched molecules\n";
      err += "   -f      Full match, print matching-molecules when the number\n";
      err += "           of heavy atoms is equal to the number of PATTERN atoms\n";
      err += "   -n      Only print the name of the molecules\n";
      err += "   -t NUM  Print a molecule only if the PATTERN occurs NUM times inside the molecule\n";
      ThrowError(err);
      exit(-1);
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
  vector< vector <int> > maplist;      // list of matched atoms
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
      
      
      ////////////////////////////////////////////////////////////////
      // Do not loose time trying to match the pattern if the matching
      // is impossible.
      // It is impossible to make a full match if the number of atoms is
      // different
      if (full ) 
	impossible_match = (sp.NumAtoms() == mol.NumHvyAtoms()) ? false : true;
      else
	impossible_match = false;

      if (impossible_match) { // -> avoid useless SMART matching attempt 
	if (invert) {
	  if (!count) {
	    if ( name_only )
	      cout << mol.GetTitle() << endl;
	    else
	      cout << mol;
	  }
	  c++;
	}
	continue;
      }


      ////////////////////////////////////////////////////////////////
      // perform SMART matching

      pattern_matched = sp.Match(mol);
      
      // the number of times the match occured may matter
      if ( ntimes ) { // ntimes is a positive integer of requested matches   
	// Here, a match mean a unique match (same set of atoms)
	// so we need to get the unique match list size

	maplist = sp.GetUMapList();                

	if( maplist.size() == ntimes )
	  ntimes_matched = true;
	else
	  ntimes_matched = false;
      } 
      else  {  // ntimes == 0, we don't care about the number of matches
	ntimes_matched = true; 
      }


      ////////////////////////////////////////////////////////////////
      // perform a set of tests to guess what to print out

      if ( pattern_matched == true && ntimes_matched == true) { 
	if (!invert) {      // do something only when invert flag is off
	  if (!count) {
	    if ( name_only )
	      cout << mol.GetTitle() << endl;
	    else
	      cout << mol;
	  }
	  c++;
	}
	
      }
    
      else { // The SMART pattern do not occur as many times as requested
	if (invert) {       // do something only if invert flag is on
	  if (!count) {
	    if ( name_only )
	      cout << mol.GetTitle() << endl;
	    else
	      cout << mol;
	  }
	  c++;
	}
      }
    } // end for loop


  ////////////////////////////////////////////////////////////////
  // Only print the number of matched molecules as requested
  if (count)
    cout << c << endl;

  return(1);
} 
