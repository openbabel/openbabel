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
#ifndef __BORLANDC__ 
#include "commandline.h"
#endif
#include "babelconfig.h"

#include "data.h"

#include <stdio.h>

#ifdef __sgi
#include <iostream.h>
#include <fstream.h>
#else
#include <iostream>
#include <fstream>
#endif

using namespace std;
using namespace OpenBabel;

__declspec(dllexport) int NewMain(char* inFile, char* inType, char* outFile, char* outType, int hydrogens, bool split)
{
  io_type inFileType = UNDEFINED, outFileType = UNDEFINED;
  bool gotInType = false, gotOutType = false, removeHydrogens = false;
  bool addHydrogens = false;
  //int arg, inFileArg, outFileArg;
  char *ext;
  OBFileFormat fileFormat;

  if (hydrogens == 1) {  // add hydrogens
	  addHydrogens = true;
  }
  if (hydrogens == 2) {  // remove hydrogens
	  removeHydrogens = true;
  }

  gotInType = true;
  ext = inType;
  if (extab.CanReadExtension(ext))
    inFileType = extab.FilenameToType(ext);
  else
	return 1;

  gotOutType = true;
  ext = outType;
  if (extab.CanWriteExtension(ext))
	outFileType = extab.FilenameToType(ext);
  else
	return 2;

  ifstream inFileStream(inFile);
  ofstream outFileStream(outFile);

  if (!inFileStream)
    {
      return 3;
    }
  if (!outFileStream)
    {
	  return 4;
    }

  // Finally, we can do some work!
  OBMol mol(inFileType, outFileType);

  fileFormat.ReadMolecule(inFileStream, mol, inFile);
  if (removeHydrogens)
    mol.DeleteHydrogens();
  if (addHydrogens)
    mol.AddHydrogens(false, false);
  fileFormat.WriteMolecule(outFileStream,mol);

  return 0;
}
