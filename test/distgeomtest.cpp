/**********************************************************************
distgeomtest.cpp - Unit tests for OBDistanceGeometry 3D generation and
                   stereo retention.

Copyright (C) 2024 by Geoffrey Hutchison

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

#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/stereo/stereo.h>

#ifdef HAVE_EIGEN3
#include <openbabel/distgeom.h>
#endif

#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace OpenBabel;

// Convert a 3D OBMol to canonical SMILES by doing a full SDF roundtrip so
// that StereoFrom3D is invoked and the SMILES reflects the 3D stereo.
static string canSmiFrom3D(OBMol& mol3D)
{
  OBConversion conv;
  conv.SetInAndOutFormats("sdf", "can");

  ostringstream sdfBuf;
  conv.SetOutFormat("sdf");
  conv.Write(&mol3D, &sdfBuf);

  OBMol mol2D;
  conv.SetInFormat("sdf");
  istringstream iss(sdfBuf.str());
  conv.Read(&mol2D, &iss);

  conv.SetOutFormat("can");
  string result = conv.WriteString(&mol2D, true);
  // trim trailing newlines/whitespace
  while (!result.empty() && (result.back() == '\n' || result.back() == '\r' || result.back() == '\t'))
    result.pop_back();
  return result;
}

// Read SMILES, get canonical form, generate 3D with distance geometry,
// do SDF roundtrip to force 3D stereo perception, compare canonical SMILES.
// Returns true if stereo is preserved.
#ifdef HAVE_EIGEN3
static bool doDistGeomStereoTest(const string& smiles)
{
  cout << " Testing: " << smiles << endl;

  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("can");

  // Get reference canonical SMILES from input
  OBMol mol;
  OB_REQUIRE(conv.ReadString(&mol, smiles));
  string refCan = conv.WriteString(&mol, true);
  while (!refCan.empty() && (refCan.back() == '\n' || refCan.back() == '\r'))
    refCan.pop_back();

  // Generate 3D geometry
  OBMol mol3D = mol;
  OBDistanceGeometry dg;
  bool ok = dg.GetGeometry(mol3D);
  if (!ok) {
    cout << "  FAILED: GetGeometry returned false" << endl;
    return false;
  }

  // Verify 3D coordinates were assigned
  OB_REQUIRE(mol3D.Has3D());
  OB_REQUIRE(mol3D.HasNonZeroCoords());

  // Check stereo via SDF roundtrip
  string can3D = canSmiFrom3D(mol3D);
  if (refCan != can3D) {
    cout << "  FAILED: stereo mismatch\n"
         << "  ref: " << refCan << "\n"
         << "  3D:  " << can3D << endl;
    return false;
  }
  cout << "  OK" << endl;
  return true;
}
#endif

int distgeomtest(int argc, char* argv[])
{
  int defaultchoice = 1;
  int choice = defaultchoice;
  if (argc > 1) {
    if (sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

#ifdef FORMATDIR
  char env[BUFF_SIZE];
  snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
  putenv(env);
#endif

#ifndef HAVE_EIGEN3
  cout << "HAVE_EIGEN3 not defined; distance geometry tests skipped." << endl;
  return 0;
#else

  switch (choice) {
  case 1:
    // Simple tetrahedral stereo -- both enantiomers
    OB_ASSERT( doDistGeomStereoTest("N[C@](Br)(O)C") );
    OB_ASSERT( doDistGeomStereoTest("N[C@@](Br)(O)C") );
    break;

  case 2:
    // Ring stereo: cis- and trans-decalin
    OB_ASSERT( doDistGeomStereoTest("C1CC[C@H]2[C@@H](C1)CCCC2") );   // cis
    OB_ASSERT( doDistGeomStereoTest("C1CC[C@@H]2[C@@H](C1)CCCC2") );  // trans
    break;

  case 3:
    // Bicyclic with multiple stereo centers
    OB_ASSERT( doDistGeomStereoTest("[C@H]1(NC[C@H]2[C@H]1N2)OC") );
    break;

  case 4:
    // Acyclic with multiple adjacent stereo centers
    OB_ASSERT( doDistGeomStereoTest("CCC[C@@H]([C@H](CC(C)C)C)C") );
    break;

  case 5:
    // Cis/trans double bond + tetrahedral stereo
    OB_ASSERT( doDistGeomStereoTest("C/C=C\\C") );
    OB_ASSERT( doDistGeomStereoTest("C/C=C/C") );
    break;

  case 6:
    // Moderately complex: tetracyclic terpenoid fragment
    OB_ASSERT( doDistGeomStereoTest("C[C@@H](CC(=O)OC)[C@@H]1CC[C@]2([C@@H]1CC[C@@H]3[C@@H]2C(=O)C=C4[C@@]3(C)CC[C@@H](C4)C(=O)O)C") );
    break;

  case 7:
    // Complex: glucuronide conjugate of a terpenoid with 19 stereo centers.
    // Slow (~10-60s); tests SetLowerBounds, 4D tunneling, and stereo checking.
    OB_ASSERT( doDistGeomStereoTest("O[C@@H]1[C@@H](O[C@@H]2O[C@H](C(=O)O)[C@H]([C@@H]([C@H]2O)O)O)"
                                    "[C@H](O[C@@H]([C@H]1O)C(=O)O)O[C@H]1CC[C@]2([C@H](C1(C)C)"
                                    "CC[C@@]1([C@@H]2C(=O)C=C2[C@@]1(C)CC[C@@]1([C@H]2C[C@](C)"
                                    "(CC1)C(=O)O)C)C)C") );
    break;

  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
#endif
}
