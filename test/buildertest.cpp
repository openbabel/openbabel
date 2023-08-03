#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>
#include <openbabel/forcefield.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

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
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}
