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
#include <cstdlib>
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
  OBConversion conv(&cin, &cout);
  OBFormat *inFormat = conv.FormatFromExt(FileIn);

  if (!inFormat || !conv.SetInFormat(inFormat)) {
    cerr << program_name << ": cannot read input format!" << endl;
    exit (-1);
  }
  // If we can't also use this for an output format, use XYZ
  if (!conv.SetOutFormat(inFormat))
    conv.SetOutFormat(conv.FindFormat("xyz"));

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
      pg.Setup(&mol);
      cerr << "Point Group: " << pg.IdentifyPointGroup() << endl;
      pg.Symmetrize(&mol);

      conv.Write(&mol, &cout);

    } // end for loop

  return(0);
}
