/**********************************************************************
obgrep - Open Babel molecule grep using SMARTS.

Copyright (C) 2003 Fabien Fontaine
Some portions Copyright (C) 2004-2005 Geoffrey R. Hutchison
 
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
#include <cstdlib>
#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/parsmart.h>

#ifdef _MSC_VER
	typedef char TCHAR;
	#include "getopt.h"
#else
	#include <unistd.h>
#endif

using namespace std;
using namespace OpenBabel;

///////////////////////////////////////////////////////////////////////////////
//! \brief Find the molecule(s) with or without a given SMART pattern
int main(int argc,char **argv)
{
  char c;
  unsigned int ntimes=0; // number of times SMARTS matches in a molecule
  unsigned int numMatching = 0; // number of matching molecules (for -c flag)
  bool pattern_matched=false, ntimes_matched=true;
  bool count=false, invert=false, full=false, name_only=false;
  char *FileIn = NULL, *Pattern = NULL;
  char *program_name = argv[0];
  char *iext;
  bool useInFile = true;

  OBConversion conv(&cin,&cout);
  OBFormat *pFormat = conv.FindFormat("smi"); // default format is SMILES
    
  // Parse options
  while ((c = getopt(argc, argv, "t:nvcfi:-")) != -1)
    {
#ifdef _WIN32
	    char optopt = c;
#endif
      switch (c)
        {
        case 't': // request ntimes unique matches

          c = sscanf(optarg, "%d", &ntimes);
          if (c != 1 )
            {
              cerr << program_name << ": unable to parse -t option" << endl;
              exit (-1);
            }
          break;

        case 'i':
          iext = optarg;

          //The ID provided by the OBFormat class is used as 
          // the identifying file extension. This is a slight
          // reduction in flexibility (which is not currently used)
          pFormat = conv.FindFormat(iext);
            
          if(pFormat==NULL)
            {
              cerr << program_name << ": cannot read input format!" << endl;
              exit(-1);
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

        case 'f':
          full = true;
          break;

        case '-':
          useInFile = false;
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
    }
  int index = optind;

  if (argc-index != 2 && argc-index != 1)
    {
      string err = "Usage: ";
      err += program_name;
      err += " [options] \"PATTERN\" <filename>\n";
      err += "If no filename is supplied, then obgrep will use stdin instead.\n";
      err += "Options:\n";
      err += "   -v      Invert the matching, print non-matching molecules\n";
      err += "   -c      Print the number of matched molecules\n";
      err += "   -i <format> Specify the input and output format\n";
      err += "   -f      Full match, print matching-molecules when the number\n";
      err += "           of heavy atoms is equal to the number of PATTERN atoms\n";
      err += "   -n      Only print the name of the molecules\n";
      err += "   -t NUM  Print a molecule only if the PATTERN occurs NUM times inside the molecule\n";
      cerr << err << ends;
      exit(-1);
    }
  else
    {
      Pattern = argv[index++];
      if (argc - index == 1)
        FileIn  = argv[index];
    }

  ifstream ifs;
  if (useInFile && FileIn != NULL)
    {
      // Read the file
      ifs.open(FileIn);
      if (!ifs)
        {
          cerr << program_name << ": cannot read input file!" << endl;
          exit (-1);
        }
      conv.SetInStream(&ifs);
	
	
      // Find Input filetype
      if (pFormat == NULL) {
          pFormat = conv.FormatFromExt(FileIn);
          if (pFormat == NULL)
            {
              cerr << program_name << ": cannot read input format!" << endl;
              return (-1);
            }
      }
    }

  if (! conv.SetInAndOutFormats(pFormat, pFormat))
    {
      cerr << program_name << ": cannot read or write to this file format" << endl;
      return (-1);
    }

  // Match the SMART
  OBSmartsPattern sp;
  vector< vector <int> > maplist;      // list of matched atoms
  sp.Init(Pattern);

  OBMol mol;

  bool impossible_match;

  // Search for pattern
  for (c=0;;)
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        break;


      ////////////////////////////////////////////////////////////////
      // Do not loose time trying to match the pattern if the matching
      // is impossible.
      // It is impossible to make a full match if the number of atoms is
      // different
      if (full )
        impossible_match = (sp.NumAtoms() == mol.NumHvyAtoms()) ? false : true;
      else
        impossible_match = false;

      if (impossible_match)
        { // -> avoid useless SMART matching attempt
          if (invert)
            {
              if (!count)
                {
                  if ( name_only )
                    cout << mol.GetTitle() << endl;
                  else
                    conv.Write(&mol, &cout);
                }
              numMatching++;
            }
          continue;
        }


      ////////////////////////////////////////////////////////////////
      // perform SMART matching

      pattern_matched = sp.Match(mol);

      // the number of times the match occured may matter
      if ( ntimes )
        { // ntimes is a positive integer of requested matches
          // Here, a match mean a unique match (same set of atoms)
          // so we need to get the unique match list size

          maplist = sp.GetUMapList();

          if( maplist.size() == ntimes )
            ntimes_matched = true;
          else
            ntimes_matched = false;
        }
      else
        {  // ntimes == 0, we don't care about the number of matches
          ntimes_matched = true;
        }


      ////////////////////////////////////////////////////////////////
      // perform a set of tests to guess what to print out

      if ( pattern_matched == true && ntimes_matched == true)
        {
          if (!invert)
            {      // do something only when invert flag is off
              if (!count)
                {
                  if ( name_only )
                    cout << mol.GetTitle() << endl;
                  else
                    conv.Write(&mol, &cout);
                }
              numMatching++;
            }

        }

      else
        { // The SMART pattern do not occur as many times as requested
          if (invert)
            {       // do something only if invert flag is on
              if (!count)
                {
                  if ( name_only )
                    cout << mol.GetTitle() << endl;
                  else
                    conv.Write(&mol, &cout);
                }
              numMatching++;
            }
        }
    } // end for loop


  ////////////////////////////////////////////////////////////////
  // Only print the number of matched molecules as requested
  if (count)
    {
      cout << numMatching << endl;
    }

  return(0);
}


/* obgrep man page*/
/** \page obgrep an advanced SMARTS grep program
*
* \n
* \par SYNOPSIS
*
* \b obgrep [options] '<SMARTS-pattern>' \<filename\>
*
* \par DESCRIPTION
*
* The obgrep tool can be used to search for molecules inside multi-molecule
* database files (e.g., SMILES, SDF, etc.).
*
* \par OPTIONS
*
* If only a filename is given, obgrep will attempt to guess
* the file type from the filename extension. \n\n
*
* \b -c:
*     Print the number of matches \n\n
* \b -f:
*     Full match, print matching-molecules only when the number
*     of heavy atoms is also equal to the number of atoms in the 
*     SMARTS pattern \n\n
* \b -i \<format\>:
*     Specifies input and output format, see "babel" for available formats \n\n
* \b -n:
*     Only print the name of the molecules\n\n
* \b -t \<NUM\>:
*     Print a molecule only if the pattern occurs NUM times inside the molecule\n\n
* \b -v:
*     Invert the matching, print non-matching molecules \n\n
*
* \par EXAMPLES
*  - Print all the molecules with a methylamine group: \n
*   obgrep "CN" database.smi
*  - Print all the molecules without a methylamine group: \n
*   obgrep -v "CN" database.smi
*  - Print the number of molecules without a methylamine group: \n
*   obgrep -v -c "CN" database.smi
*  - Print methylamine (if it exists in the file): \n
*   obgrep -f "CN" database.smi
*  - Print methylamine and/or methanol (if they exist): \n
*   obgrep -f "C[N,O]" database.smi
*
* \par AUTHORS
*
* The obgrep program was contributed by \b Fabien \b Fontaine.
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
*   A guide for constructing SMARTS patterns can be found at: http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
**/
