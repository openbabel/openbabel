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

void GenerateSmartsReference();
bool TestSmarts(void);
void GenerateRingReference();
bool TestRings(void);

using namespace std;
using namespace OpenBabel;

void usage();

int main(int argc,char *argv[])
{
  if (argc != 1) usage();

  //GenerateSmartsReference();
  //GenerateRingReference();

  bool alltests = true;


  if (!TestSmarts())
    {
      ThrowError("ERROR: ***SMARTS test failed***");
      alltests = false;
    }
  else ThrowError("SMARTS test passed");
    
  if (!TestRings())
    {
      ThrowError("ERROR: ***RING test failed***");
      alltests = false;
    }
  else ThrowError("RING test passed");
      

  if (alltests) ThrowError("\nAll tests passed");
  else          ThrowError("\nSome test(s) failed");  

  return(0);
}

void usage()
{
  ThrowError("Usage: obtest");
  exit(0);
}
