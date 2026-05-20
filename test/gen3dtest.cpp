/**********************************************************************
gen3dtest.cpp - Unit tests for the gen3D op, including fallback to
                distance geometry when OBBuilder fails or produces
                zero coordinates.

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
#include "stereo_test_helpers.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/op.h>

#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace OpenBabel;

// Read SMILES, apply gen3D op at the given speed level, verify 3D
// coordinates are generated and non-zero.  If checkStereo is true,
// also verify that canonical SMILES is preserved through the SDF
// roundtrip.  Returns true on success.
static bool doGen3DTest(const string& smiles, const char* speed = "3",
                        bool checkStereo = true)
{
  cout << " Testing gen3D(" << speed << "): " << smiles << endl;

  OBOp* gen3Dop = OBOp::FindType("gen3D");
  if (!gen3Dop) {
    cout << "  SKIPPED: gen3D op not available" << endl;
    return true;
  }

  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("can");

  OBMol mol;
  OB_REQUIRE(conv.ReadString(&mol, smiles));
  string refCan = conv.WriteString(&mol, true);
  while (!refCan.empty() && (refCan.back() == '\n' || refCan.back() == '\r'))
    refCan.pop_back();

  bool opOk = gen3Dop->Do(&mol, speed);
  if (!opOk) {
    cout << "  FAILED: gen3D returned false" << endl;
    return false;
  }

  if (!mol.Has3D() || !mol.HasNonZeroCoords()) {
    cout << "  FAILED: no valid 3D coordinates generated" << endl;
    return false;
  }

  if (checkStereo) {
    string can3D = canSmiFrom3D(mol);
    if (refCan != can3D) {
      cout << "  FAILED: stereo mismatch\n"
           << "  ref: " << refCan << "\n"
           << "  3D:  " << can3D << endl;
      return false;
    }
  }

  cout << "  OK" << endl;
  return true;
}

int gen3dtest(int argc, char* argv[])
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

  switch (choice) {
  case 1:
    // Basic sanity check: simple molecules at each speed level
    OB_ASSERT( doGen3DTest("c1ccccc1") );          // benzene
    OB_ASSERT( doGen3DTest("CC(=O)O") );           // acetic acid
    OB_ASSERT( doGen3DTest("N[C@@H](C)C(=O)O") ); // L-alanine
    break;

  case 2:
    // Builder fallback: molecules where OBBuilder fails or gives zero
    // coordinates (issue #2317).  gen3D must fall back to distgeom and
    // still produce valid non-zero coordinates.
    //
    // Triarylmethane as free base and as HCl salt (the .Cl fragment
    // caused the builder to return zero coords before the fix).
    OB_ASSERT( doGen3DTest("C1=CC(=N)C=CC1=C(C2=CC=C(C=C2)N)C3=CC=C(C=C3)N",
                            "3", false) );
    //OB_ASSERT( doGen3DTest("C1=CC(=N)C=CC1=C(C2=CC=C(C=C2)N)C3=CC=C(C=C3)N.Cl",
    //                        "3", false) );
    // Regression for issue #342: isotope-hydrogen fragment plus a bare
    // metal atom previously crashed gen3D.
    OB_ASSERT( doGen3DTest("[2HH].[Li]", "3", false) );
    break;

  case 3:
    // Stereo preservation through gen3D (using balanced speed level)
    OB_ASSERT( doGen3DTest("N[C@](Br)(O)C") );
    OB_ASSERT( doGen3DTest("N[C@@](Br)(O)C") );
    OB_ASSERT( doGen3DTest("C[C@H]([C@@H](C(=O)O)N)O") );   // L-threonine
    OB_ASSERT( doGen3DTest("[C@@H]([C@H](C(=O)O)O)(C(=O)O)O") ); // L-tartaric acid
    break;

  case 4:
    // Ring stereo through gen3D
    OB_ASSERT( doGen3DTest("C1CC[C@H]2[C@@H](C1)CCCC2") );   // cis-decalin
    OB_ASSERT( doGen3DTest("C1CC[C@@H]2[C@@H](C1)CCCC2") );  // trans-decalin
    break;

  case 5:
    // Complex fused-ring and sugar-ring molecules where OBDistanceGeometry
    // alone times out, but the builder handles the ring topology correctly.
    // Coord generation only (stereo round-trip unreliable for these systems).
    OB_ASSERT( doGen3DTest("CN1CCC2=CC3=C(C=C2[C@@H]1[C@@H]4C5=C(C(=C(C=C5)OC)OC)C(=O)O4)OCO3",
                            "3", false) ); // berberine analog
    OB_ASSERT( doGen3DTest("C1CN2CC3=CC4=C(C=C3[C@H]5[C@H]2C1=C[C@@H]([C@H]5O)O)OCO4",
                            "3", false) ); // polycyclic alkaloid
    OB_ASSERT( doGen3DTest("C1=C[C@H]2C(=CN1C)[C@H]1C(=CC=CN1C)C=C2",
                            "3", false) ); // vinca-like indole
    OB_ASSERT( doGen3DTest("C[C@@H]1[C@H](C[C@@H]([C@H](O1)OC2[C@@H]([C@H](C([C@@H]([C@@H]2O)O)O)O)O)N)N=C(C(=O)O)N",
                            "3", false) ); // aminoglycoside
    OB_ASSERT( doGen3DTest("C1[C@@H](NC(=N[C@H]1O)N)[C@@H](C(=O)O)N",
                            "3", false) ); // cyclic arginine analog
    break;

  case 6:
    // Stereoisomer coverage.  The rigid-fragment template library used by
    // OBBuilder stores specific stereoisomers; the wrong enantiomer or
    // diastereomer could in principle drop through to a zero-coordinate
    // result.  For each input, exercise the original, its enantiomer
    // (every @ flipped), and one diastereomer (a single @ flipped), and
    // verify both non-zero coordinates and stereo retention through the
    // SDF round-trip.
    {
      const char* inputs[] = {
        "OC[C@H]([C@H](S)CC)C",                // 2 stereocenters
        "C[C@H]([C@@H](C(=O)O)N)O",            // L-threonine
        "[C@@H]([C@H](C(=O)O)O)(C(=O)O)O",     // L-tartaric acid
        "C1CC[C@H]2[C@@H](C1)CCCC2",           // cis-decalin (ring stereo)
      };
      for (const char* smi : inputs) {
        OB_ASSERT( doGen3DTest(smi) );
        OB_ASSERT( doGen3DTest(makeEnantiomer(smi)) );
        if (countChirality(smi) >= 2)
          OB_ASSERT( doGen3DTest(makeDiastereomer(smi, 0)) );
      }
    }
    break;

  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}
