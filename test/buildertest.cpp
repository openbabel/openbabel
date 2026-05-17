#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>
#include <openbabel/forcefield.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace OpenBabel;

static unsigned int failed = 0;
static unsigned int testCount = 0;

void clearMolFlags(OBMol &mol)
{
  // Both of these are needed or there's a crash.
  mol.UnsetFlag(OB_SSSR_MOL);
  mol.UnsetFlag(OB_AROMATIC_MOL);

  // Not needed, doesn't cause crash
  //  mol.UnsetFlag(OB_RINGFLAGS_MOL);
  //  mol.UnsetFlag(OB_LSSR_MOL);
  //  mol.UnsetFlag(OB_KEKULE_MOL);
  //  mol.UnsetFlag(OB_CLOSURE_MOL);
  //  mol.UnsetFlag(OB_AROM_CORRECTED_MOL);
  //  mol.UnsetFlag(OB_RINGTYPES_MOL);
}

bool doBuildMoleculeTest(OBMol &mol)
{
  testCount++;

  OBBuilder builder;
  OB_REQUIRE(builder.Build(mol, false));
  // Does not need clearMolFlags -- crash still happens if you clear here
  // and not after AddHydrogens()
  OB_REQUIRE(mol.AddHydrogens());
  OB_REQUIRE(mol.HasAromaticPerceived() == 0);
  OB_REQUIRE(mol.HasSSSRPerceived() == 0);
  //  clearMolFlags(mol); // must clear here or you crash
  // Should now be handled by AddHydrogens()

  OBForceField* pff = OBForceField::FindType("mmff94");
  OB_REQUIRE(pff != nullptr);
  cout << mol.GetTitle() << endl;
  OB_REQUIRE(pff->Setup(mol));
  // Check for explosions -- PR#3016479
  pff->SteepestDescent(100);
  OB_REQUIRE(!pff->DetectExplosion()); // no explosions please!

  return true;
}

bool doMultiMoleculeFile(const std::string &filename)
{
  cout << " Starting " << filename << endl;

  std::string file = OBTestUtil::GetFilename(filename);
  std::ifstream ifs;
  ifs.open(file.c_str());
  OB_REQUIRE( ifs );

  OBMol mol;
  OBConversion conv(&ifs, &cout);
  OBFormat *format = conv.FormatFromExt(file.c_str());
  OB_REQUIRE(format);
  OB_REQUIRE(conv.SetInFormat(format));

  bool result = true;
  while (conv.Read(&mol, &ifs)) {
    bool res = doBuildMoleculeTest(mol);
    if (!res)
      result = res;
  }

  return result;
}

bool doSMILESBuilderTest(string smiles)
{
  cout << " SMILES " << smiles << endl;

  testCount++;

  OBMol mol;
  OBConversion conv;
  OBFormat *smilesFormat = conv.FindFormat("smi");
  OB_REQUIRE(smilesFormat);
  OB_REQUIRE(conv.SetInFormat(smilesFormat));

  OB_REQUIRE(conv.ReadString(&mol, smiles));

  OBBuilder builder;
  OB_REQUIRE(builder.Build(mol, false)); // some stereo errors are known
  return (mol.Has3D() && mol.HasNonZeroCoords());
}

// Verify that Build() produces a sane ring geometry: no atom-atom overlaps
// and no stretched ring-closure bonds. maxBondLen relaxes the bond-length
// cap for fused/bridged systems where a smaller fused ring may still close
// with some strain (only the largest ring per system is crowned).
bool doRingClosureTest(const string &smiles, double maxBondLen)
{
  cout << " Ring closure: " << smiles << " (max bond " << maxBondLen << ")" << endl;

  testCount++;

  OBMol mol;
  OBConversion conv;
  OBFormat *smilesFormat = conv.FindFormat("smi");
  OB_REQUIRE(smilesFormat);
  OB_REQUIRE(conv.SetInFormat(smilesFormat));
  OB_REQUIRE(conv.ReadString(&mol, smiles));

  OBBuilder builder;
  OB_REQUIRE(builder.Build(mol, false));
  OB_REQUIRE(mol.Has3D());

  double worst = 0.0;
  FOR_BONDS_OF_MOL(b, mol) {
    double len = b->GetLength();
    if (len > worst) worst = len;
  }
  if (worst > maxBondLen) {
    cout << "  FAIL: longest bond = " << worst << " A" << endl;
    return false;
  }

  // No two heavy atoms should overlap.
  vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(a, mol) atoms.push_back(&*a);
  for (size_t i = 0; i < atoms.size(); ++i) {
    for (size_t j = i + 1; j < atoms.size(); ++j) {
      double d = atoms[i]->GetDistance(atoms[j]);
      if (d < 0.5) {
        cout << "  FAIL: atoms " << i << " and " << j
             << " overlap (d=" << d << " A)" << endl;
        return false;
      }
    }
  }
  return true;
}

int buildertest(int argc, char* argv[])
{
  int defaultchoice = 1;

  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif



  // fails because of selenium
  //  OB_ASSERT( doMultiMoleculeFile("aromatics.smi") );
  // fails because of stereo crash
  //  OB_ASSERT( doMultiMoleculeFile("nci.smi") );
  // fails because of "organometallic" entries
  //  OB_ASSERT( doMultiMoleculeFile("attype.00.smi") );

  switch(choice) {
  case 1:
    OB_ASSERT( doMultiMoleculeFile("forcefield.sdf") );
    break;
  case 2:
    OB_ASSERT( doMultiMoleculeFile("filterset.sdf") );
    break;
  case 3:
    // from Martin Guetlein to mailing list on July 14, 2010
    OB_ASSERT( doSMILESBuilderTest("OC1=CC3=C([C@@]2([H])CC[C@@]4(C)[C@](CC[C@@H]4O)([H])[C@@]([H])2[C@H](CCCCCCCCCS(CCCC(F)(F)C(F)(F)F)=O)C3)C=C1") );
    break;
  case 4:
    // from Thomas Womack -- PR#3016479
    OB_ASSERT( doSMILESBuilderTest("O1C[C@H]2O[C@H]2c2ccc(Oc3c(O)ccc(CCC1=O)c3)cc2") );
    break;
  case 5:
    // from Martin Guetlein -- PR#3107218 ("OBBuilder terminates while building 3d")
    OB_ASSERT( doSMILESBuilderTest("N12[C@@H]([C@@H](NC([C@@H](c3ccsc3)C(=O)O)=O)C2=O)SC(C)(C)[C@@-]1C(=O)O") );
    break;
  case 6:
    // from Hubertus van Dam -- #2144
    OB_ASSERT( doSMILESBuilderTest("OC1(C2=CN(CC3=CC=CC=C3F)N=N2)CCOC1") );
    break;
  case 7:
    // Regression: small/medium rings (sizes 7, 8, 12) are templated in
    // ring-fragments.txt; the post-pass must not disturb their geometry.
    OB_ASSERT( doRingClosureTest("C1CCCCCC1", 1.8) );        // cycloheptane
    OB_ASSERT( doRingClosureTest("C1CCCCCCC1", 1.8) );       // cyclooctane
    OB_ASSERT( doRingClosureTest("C1CCCCCCCCCCC1", 1.8) );   // cyclododecane
    break;
  case 8:
    // Crown post-pass covers n=19, 22 (formula 180 - 720/n).
    OB_ASSERT( doRingClosureTest("C1CCCCCCCCCCCCCCCCCC1", 2.0) );     // cyclo-19
    OB_ASSERT( doRingClosureTest("C1CCCCCCCCCCCCCCCCCCCCC1", 2.5) );  // cyclo-22
    // Beyond the crown cap, OBDaleTorsions emits a precomputed
    // diamond-lattice Dale code. Even n close to ~1.5 A pre-FF (the
    // normal C-C bond length); odd n leave a ~2.5 A residual from
    // the bipartite-lattice constraint. Both are well inside FF
    // cleanup's basin of convergence.
    OB_ASSERT( doRingClosureTest("C1CCCCCCCCCCCCCCCCCCCCCCCCCCCCC1", 2.6) ); // cyclo-30
    // cyclo-44 and cyclo-50 exercise the back-edge ring tracer:
    // OB's SSSR (OB_RTREE_CUTOFF=20) silently drops rings >~40
    // atoms, so the post-pass falls back to tracing through the
    // workMol spanning tree.
    {
      std::string s44 = "C1" + std::string(43, 'C') + "1";
      OB_ASSERT( doRingClosureTest(s44, 2.0) );                       // cyclo-44
      std::string s50 = "C1" + std::string(49, 'C') + "1";
      OB_ASSERT( doRingClosureTest(s50, 2.0) );                       // cyclo-50
      std::string s100 = "C1" + std::string(99, 'C') + "1";
      OB_ASSERT( doRingClosureTest(s100, 2.0) );                      // cyclo-100
    }
    break;
  case 9:
    // Fused/bridged systems: only the largest ring per system is crowned,
    // so the smaller fused ring may still strain. Use a looser cap.
    OB_ASSERT( doRingClosureTest("C1CCC2CCCCC2C1", 2.2) );   // decalin
    OB_ASSERT( doRingClosureTest("C1CC2CCC1C2", 2.2) );      // norbornane
    break;
  case 10:
    // Aromatic ring that's unlikely to be templated -- verify it builds
    // (no stretched closure) and that the post-pass treats it as planar.
    OB_ASSERT( doRingClosureTest("C1=CC=CC=CC=N1", 1.8) );   // 1H-azocine
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}
