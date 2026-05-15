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
// Verify only that GetGeometry produces non-zero 3D coordinates.
// Use this when stereo cannot be verified via SMILES roundtrip (e.g. ring
// double bonds that StereoFrom3D doesn't perceive, or bridged bicyclics
// where the embedding is too slow to retry to convergence).
#ifdef HAVE_EIGEN3
static bool doDistGeomCoordsTest(const string& smiles)
{
  cout << " Testing coords: " << smiles << endl;

  OBConversion conv;
  conv.SetInFormat("smi");

  OBMol mol;
  OB_REQUIRE(conv.ReadString(&mol, smiles));

  OBDistanceGeometry dg;
  bool ok = dg.GetGeometry(mol);
  if (!ok) {
    cout << "  FAILED: GetGeometry returned false" << endl;
    return false;
  }

  OB_REQUIRE(mol.Has3D());
  OB_REQUIRE(mol.HasNonZeroCoords());

  cout << "  OK" << endl;
  return true;
}

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

  case 8:
    // Medium-sized rings with double bonds -- coord generation only.
    // E/Z notation on ring double bonds creates OBCisTransStereo constraints
    // that CheckStereoConstraints() can never satisfy, so we use plain SMILES
    // with the ring stereo stripped.
    OB_ASSERT( doDistGeomCoordsTest("C1CCCC=CCCC1") );               // cyclonon-4-ene
    OB_ASSERT( doDistGeomCoordsTest("C1CCCC=CCCCCCCCC(=O)CCC1") );   // 17-membered macrolide
    break;

  case 9:
    // Large rings and macrocyclic polyenes -- coord generation only.
    // Same reason as case 8: ring E/Z stereo stripped.
    OB_ASSERT( doDistGeomCoordsTest("C1=CC=CC=CC=CC=CC=CC=CC=C1") ); // [16]-annulene
    OB_ASSERT( doDistGeomCoordsTest("CC1=CCC(C=CCC(=CCC1)C)(C)C") ); // germacrene sesquiterpene
    break;

  case 10:
    // Open-chain monosaccharide stereochemistry
    OB_ASSERT( doDistGeomStereoTest("C([C@H]([C@@H]([C@@H]([C@H](CO)O)O)O)O)O") );     // galactose
    OB_ASSERT( doDistGeomStereoTest("C([C@H]([C@H]([C@@H]([C@H](C(=O)O)O)O)O)O)O") );  // glucuronic acid
    OB_ASSERT( doDistGeomStereoTest("C([C@H]([C@H]([C@@H]([C@H](CO)O)O)O)O)O") );       // glucose
    OB_ASSERT( doDistGeomStereoTest("C([C@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O)O") );      // mannose
    break;

  case 11:
    // Small chiral molecules: amino acids, hydroxy acids, diols
    OB_ASSERT( doDistGeomStereoTest("C[C@H]([C@@H](C)C(=O)O)C(=O)O") );           // dimethylsuccinic acid
    OB_ASSERT( doDistGeomStereoTest("[C@@H]([C@H](C(=O)O)O)(C(=O)O)O") );         // L-tartaric acid
    OB_ASSERT( doDistGeomStereoTest("C[C@H]([C@@H](C(=O)O)N)O") );                // L-threonine
    OB_ASSERT( doDistGeomStereoTest("C[C@H]([C@@H](C)O)O") );                     // butane-2,3-diol
    OB_ASSERT( doDistGeomStereoTest("[C@@H]([C@H](C(=O)N)O)(C(=O)N)O") );         // asparagine-diol
    break;

  case 12:
    // Halogenated stereocenters and chloramphenicol analog
    OB_ASSERT( doDistGeomStereoTest("[C@@H]([C@@H](C(=O)O)Br)(C(=O)O)Br") );      // dibromo succinic acid
    OB_ASSERT( doDistGeomStereoTest("C1=CC(=CC=C1[C@H]([C@@H](CO)NC(=O)C(Cl)Cl)O)[N+](=O)[O-]") ); // chloramphenicol analog
    break;

  case 13:
    // Quinine and quinidine (cinchona alkaloids).
    // NOT included in distgeom_parts: the bridged quinuclidine core causes
    // the L-BFGS inside each distgeom trial to exhaust its iteration budget
    // without converging, making all 10*N trials slow (~minutes total).
    // These molecules are tested via gen3dtest case 2, where OBBuilder
    // handles the ring topology and distgeom is only a fallback.
    OB_ASSERT( doDistGeomCoordsTest("C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@@H](C3=CC=NC4=CC=CC=C34)O") );  // quinine
    OB_ASSERT( doDistGeomCoordsTest("C=C[C@H]1CN2CC[C@H]1C[C@@H]2[C@H](C3=CC=NC4=CC=CC=C34)O") ); // quinidine
    break;

  case 14:
    // Steroid and terpenoid ring systems
    OB_ASSERT( doDistGeomStereoTest("C[C@]12CC[C@H]3[C@H]([C@@H]1C[C@H]([C@@H]2O)O)CCC4=C3C=CC(=C4)O") ); // estradiol-like
    OB_ASSERT( doDistGeomStereoTest("C[C@H]1C[C@@H](C(=O)[C@@H](C1)[C@@H](CC2CC(=O)NC(=O)C2)O)C") );      // terpenoid-lactam
    break;

  case 15:
    // Isoquinoline and indole polycyclic alkaloids.
    // NOT in distgeom_parts: fused ring distance constraints exhaust the
    // L-BFGS budget on every trial (same failure mode as quinine/case 13).
    // These are tested via gen3dtest case 5 using the builder path.
    OB_ASSERT( doDistGeomCoordsTest("CN1CCC2=CC3=C(C=C2[C@@H]1[C@@H]4C5=C(C(=C(C=C5)OC)OC)C(=O)O4)OCO3") ); // berberine analog
    OB_ASSERT( doDistGeomCoordsTest("C1CN2CC3=CC4=C(C=C3[C@H]5[C@H]2C1=C[C@@H]([C@H]5O)O)OCO4") );          // polycyclic alkaloid
    OB_ASSERT( doDistGeomCoordsTest("C1=C[C@H]2C(=CN1C)[C@H]1C(=CC=CN1C)C=C2") );                           // vinca-like indole
    break;

  case 16:
    // Aminoglycoside and cyclic guanidino stereo.
    // NOT in distgeom_parts: the pyranose ring system in the aminoglycoside
    // has the same L-BFGS convergence problem as case 15.
    // Tested via gen3dtest case 5 using the builder path.
    OB_ASSERT( doDistGeomCoordsTest("C[C@@H]1[C@H](C[C@@H]([C@H](O1)OC2[C@@H]([C@H](C([C@@H]([C@@H]2O)O)O)O)O)N)N=C(C(=O)O)N") ); // aminoglycoside
    OB_ASSERT( doDistGeomCoordsTest("C1[C@@H](NC(=N[C@H]1O)N)[C@@H](C(=O)O)N") );   // cyclic arginine analog
    break;

  case 17:
    // Amino acid derivatives and dipeptide fragments
    OB_ASSERT( doDistGeomStereoTest("CC(C)C[C@@H](C(=O)O)NC(=O)[C@H]([C@@H](CC1=CC=CC=C1)N)O") ); // dipeptide fragment
    OB_ASSERT( doDistGeomStereoTest("C[C@H]([C@@H](C(=O)O)N)OP(=O)(O)O") );                        // phosphoamino acid
    OB_ASSERT( doDistGeomStereoTest("C[C@H]([C@@H](C(=O)O)N)SC[C@@H](C(=O)O)N") );                 // cystine fragment
    break;

  case 18:
    // Complex multi-stereo-center molecules
    OB_ASSERT( doDistGeomStereoTest("Cc1nnc(CNC[C@@H]2CN(C(=O)[C@@]34CCCC[C@H]3C4)C[C@H]2C)n1C1CC1") );            // bicyclic proline-triazole
    OB_ASSERT( doDistGeomStereoTest("N1(C=C[C@@H](C=C1C)[C@H]1C=CN(C(=C1)C)CCCl)CCCl") );                          // bis-dihydropyridinium
    OB_ASSERT( doDistGeomStereoTest("Cc1ccc(-c2cccc([C@@H]3C[C@](C)(c4ccccc4)c4cc(C(=N)N)ccc4N3)c2)c(C(=O)O)c1") ); // biaryl amidine
    break;

  case 19:
    // Disaccharide-azo dye conjugate (many stereocenters).
    // NOT in distgeom_parts: 38+ heavy atoms → maxIter=380 trials; reliably
    // exceeds the 30 s wall-clock limit.  Coordinate generation is tested
    // via gen3dtest (builder path) instead.
    OB_ASSERT( doDistGeomCoordsTest("OC[C@H]1O[C@@H](Oc2ccc(N=Nc3ccccc3)cc2)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O") );
    break;

  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
#endif
}
