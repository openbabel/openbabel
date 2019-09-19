#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>

#include <openbabel/graphsym.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/canon.h>

#include <iostream>
#include <vector>
#include <algorithm>

//#define WITH_DEPICTION 1

using std::cout;
using std::endl;
using namespace OpenBabel;

static unsigned int N = 10;

void compareMolecules(OBMol *mol1, OBMol *mol2)
{
  std::vector<unsigned int> canlbls1, canlbls2;
  std::vector<unsigned int> symclasses1, symclasses2;

  CanonicalLabels(mol1, symclasses1, canlbls1);
  cout << "mol1.NumAtoms = " << mol1->NumAtoms() << endl;

  CanonicalLabels(mol2, symclasses2, canlbls2);
  cout << "mol2.NumAtoms = " << mol2->NumAtoms() << endl;

  std::vector<unsigned int> symclassesCopy1 = symclasses1;
  std::sort(symclassesCopy1.begin(), symclassesCopy1.end());
  std::vector<unsigned int>::iterator end1 = std::unique(symclassesCopy1.begin(), symclassesCopy1.end());
  unsigned int unique1 = end1 - symclassesCopy1.begin();

  std::vector<unsigned int> symclassesCopy2 = symclasses2;
  std::sort(symclassesCopy2.begin(), symclassesCopy2.end());
  std::vector<unsigned int>::iterator end2 = std::unique(symclassesCopy2.begin(), symclassesCopy2.end());
  unsigned int unique2 = end2 - symclassesCopy2.begin();

  OB_ASSERT( unique1 == unique2 );
  if (unique1 != unique2)
    cout << unique1 << " == " << unique2 << endl;

  FOR_ATOMS_OF_MOL (a1, mol1) {
    OBAtom *a2 = 0;
    unsigned int symClass1 = symclasses1.at(a1->GetIndex());
    for (unsigned int i = 0; i < symclasses2.size(); ++i)
      if (symclasses2.at(i) == symClass1) {
        a2 = mol2->GetAtom(i+1);
        break;
      }

    if (!a2)
      continue;

    OB_ASSERT( a1->GetAtomicNum() == a2->GetAtomicNum() );
    OB_ASSERT( a1->GetExplicitDegree() == a2->GetExplicitDegree() );
    OB_ASSERT( a1->GetHvyDegree() == a2->GetHvyDegree() );
    OB_ASSERT( a1->GetHeteroDegree() == a2->GetHeteroDegree() );
    OB_ASSERT( a1->GetImplicitHCount() == a2->GetImplicitHCount() );
  }

}




bool doRoundtrip(OBMol *mol)
{
  cout << "==================== start doRoundtrip =========================" << endl;

  OBConversion conv;
  conv.SetInAndOutFormats("can", "can");

  OBMol round2; // result of reading first output as canonical SMILES
  std::string output, roundtrip; // first output string, then the roundtrip

  mol->SetTitle("");
  output = conv.WriteString(mol, true); // trim whitespace
  OB_REQUIRE(conv.ReadString(&round2, output));

  bool result = true;
  round2.SetTitle("");
  roundtrip = conv.WriteString(&round2, true);
  if (roundtrip != output) {
    cout << "Failed roundtrip: \n  " << output << " -> \n  " << roundtrip << endl;
    result = false;
  }
  OB_ASSERT( roundtrip == output );

//  compareMolecules(mol, &round2);

  cout << "==================== end doRoundtrip =========================" << endl;

  return result;
}








bool doShuffleTest(const std::string &smiles)
{
  cout << " True Shuffling: " << smiles << endl;
  // read a smiles string
  OBMol mol;
  OBConversion canConv, smiConv;
  OB_REQUIRE( canConv.SetInFormat("smi") );
  OB_REQUIRE( canConv.SetOutFormat("can") );
  OB_REQUIRE( smiConv.SetOutFormat("smi") );
  // read a smiles string
  OB_REQUIRE( canConv.ReadString(&mol, smiles) );


  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);


  std::string ref = canConv.WriteString(&mol, true); // FIXME
  cout << "ref = " << ref << endl;

  bool result = true;
  for (unsigned int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string cansmi = canConv.WriteString(&mol, true);
    OB_ASSERT( cansmi == ref );

    doRoundtrip(&mol);
  }

  bool roundtrip = doRoundtrip(&mol);
  OB_ASSERT( result == roundtrip );

  return result;
}

bool doShuffleTestFile(const std::string &filename)
{
  cout << " True Shuffling: " << filename << endl;
  std::string file = OBTestUtil::GetFilename(filename);
  // read a smiles string
  OBMol mol;
  OBConversion canConv, smiConv;
  OBFormat *format = canConv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_REQUIRE( canConv.SetInFormat(format) );
  OB_REQUIRE( canConv.ReadFile(&mol, file) );
  OB_REQUIRE( canConv.SetOutFormat("can") );
  OB_REQUIRE( smiConv.SetOutFormat("smi") );


  std::string smiles = canConv.WriteString(&mol, true);

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);

  std::string ref = canConv.WriteString(&mol, true); // FIXME
  cout << "ref = " << ref << endl;

  bool result = true;
  for (unsigned int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);

    // get can smiles
    std::string cansmi = canConv.WriteString(&mol, true); // FIXME
    OB_ASSERT( cansmi == ref );
  }

  bool roundtrip = doRoundtrip(&mol);
  OB_ASSERT( result == roundtrip );

  return result;
}

bool doShuffleTestOnMultiFile(const std::string &filename)
{
  cout << " True Shuffling: " << filename << endl;
  std::string file = OBTestUtil::GetFilename(filename);
  // read a smiles string
  OBMol mol;
  OBConversion canConv;
  OBFormat *format = canConv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_REQUIRE( canConv.SetInFormat(format) );
  OB_REQUIRE( canConv.SetOutFormat("can") );


  std::ifstream ifs;
  ifs.open(file.c_str());
  OB_REQUIRE( ifs );

  bool result = true;
  while (canConv.Read(&mol, &ifs)) {

    std::vector<OBAtom*> atoms;
    FOR_ATOMS_OF_MOL(atom, mol)
      atoms.push_back(&*atom);

    mol.SetTitle("");
    std::string ref = canConv.WriteString(&mol, true);
    cout << "ref = " << ref << endl;

    bool subresult = true;
    for (unsigned int i = 0; i < N; ++i) {
      // shuffle the atoms
      std::random_shuffle(atoms.begin(), atoms.end());
      mol.RenumberAtoms(atoms);

      // get can smiles
      mol.SetTitle("");
      std::string cansmi = canConv.WriteString(&mol, true);
      // comapare with ref
      OB_ASSERT( cansmi == ref );
    }

    bool roundtrip = doRoundtrip(&mol);
    OB_ASSERT( subresult == roundtrip );
  }

  return result;
}

int shuffletest(int argc, char* argv[])
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

  /*
  if (argc == 3) {
    OB_ASSERT( doShuffleTestOnMultiFile(argv[2]) );
    cout << "PASSED TESTS: " << testCount - failed << "/" << testCount << endl;
    return 0;
  }*/

  switch(choice) {
  case 1:
    OB_ASSERT( doShuffleTest("COc1cc2c(cc1O)CN1C[C@@H]2c2cc(OC)c(cc2C1)O") );
    OB_ASSERT( doShuffleTest("C(C1C=CCCC1)(C(=O)N/N=C/c1ccc(cc1)C#N)C(=O)N/N=C\\c1ccc(cc1)C#N") );
    OB_ASSERT( doShuffleTest("C(C1C=CCCC1)(C(=O)N/N=C\\c1ccc(cc1)C#N)C(=O)N/N=C/c1ccc(cc1)C#N") );


    OB_ASSERT( doShuffleTest("CC(=O)O.N=C\\1/C=C/C(=C(/c2ccc(N)cc2)\\c2ccc(N)c(C)c2)/C=C1") );
    OB_ASSERT( doShuffleTest("CCNc1ccc(cc1)/C(=C/1\\C=C/C(=N)/C=C1)/c1ccc(N(CC)CC)c(C)c1") );
    OB_ASSERT( doShuffleTest("C[C@H]1CCCCN1/C=C/1\\C(=O)O[C@@](C)(CC)OC1=O") );
    OB_ASSERT( doShuffleTest("N=C\\1/C=C/C(=C(\\c2ccc(N)cc2)/c2ccc(N)c(C)c2)/C=C1") );
    OB_ASSERT( doShuffleTest("N=C\\1/C=C/C(=N\\Cc2ccc(N)cc2)/C=C1") );

    OB_ASSERT( doShuffleTest("[O+]#C[Fe+]1234(C#[O+])(/C=C/CCC/C=C/[Fe+]5678(C#[O+])(C#[O+])[C@H]9C7=C6C5=C89)[C@H]5C3=C2C1=C45.F[B-](F)(F)F") );
    OB_ASSERT( doShuffleTest("[O+]#C[Fe+]1234(C#[O+])(/C=C/CCCC/C=C/[Fe+]5678(C#[O+])(C#[O+])C9=C6C7=C5[C@@H]89)C5=C3C2=C1[C@H]45.F[B-](F)(F)F") );


    OB_ASSERT( doShuffleTest("c1[14cH]cccc1") );
    OB_ASSERT( doShuffleTest("[14cH]1[14cH]cccc1") );
    OB_ASSERT( doShuffleTest("[14cH]1[14cH]ccc[14cH]1") );
    OB_ASSERT( doShuffleTest("[14cH]1[14cH]cc[14cH][14cH]1") );
    OB_ASSERT( doShuffleTest("[14cH]1[14cH]c[14cH][14cH][14cH]1") );
    OB_ASSERT( doShuffleTest("[14cH]1[14cH][14cH][14cH][14cH][14cH]1") );


    OB_ASSERT( doShuffleTestOnMultiFile("stereo/error2.smi") ); // FYI, kekulization error somewhere in here

    OB_ASSERT( doShuffleTest("F[Po@SP1](Cl)(Br)I") );
    OB_ASSERT( doShuffleTest("F[Po@SP2](Br)(Cl)I") );
    OB_ASSERT( doShuffleTest("F[Po@SP3](Cl)(I)Br") );
    break;
  case 2:
    //OB_ASSERT( doShuffleTest("C12=C3C4=C5[C@@H]1[Fe]16782345C2=C7[C@@]8(C6=C12)[Ge]([C@]12C3=C4C5=C1[Fe]16782345C2=C7C6=C1[C@H]82)([C@@]12C3=C4C5=C1[Fe]16782345C2=C7C6=C1[C@H]82)[C@]12C3=C4C5=C1[Fe]16782345C2=C7C6=C1[C@H]82") );
    OB_ASSERT( doShuffleTest("Br.CC(=O)Nc1ccc(cc1)S(=O)(=O)c1ccc(cc1)NC(=O)c1ccccc1SC(=O)CCCCn1ccccc1") );
    OB_ASSERT( doShuffleTest("Br.CC(=O)Nc1cccc(c1)S(=O)(=O)c1cccc(NC(=O)c2ccccc2SC(=O)CCCCn2ccccc2)c1") );
    OB_ASSERT( doShuffleTest("Br.COc1ccc(cc1)Cn1cccc(O)c1") );
    OB_ASSERT( doShuffleTest("Br.COc1ccc2c(c1)[C@]13CCCC[C@@H]3[C@H](N(C)CC1)[C@@]12CC1") );
    OB_ASSERT( doShuffleTest("Br.NCCCCNCCCC[C@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82") );
    OB_ASSERT( doShuffleTest("Br.NCCCCNCCCC[C@@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82") );
    OB_ASSERT( doShuffleTest("Br.NCCCNCCCCNCCCC[C@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82") );
    OB_ASSERT( doShuffleTest("Br.NCCCNCCCCNCCCC[C@@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82") );
    OB_ASSERT( doShuffleTest("Br.NCCCNCCCCNCCCNCCCC[C@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82") );
    OB_ASSERT( doShuffleTest("Br.NCCCNCCCCNCCCNCCCC[C@@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82") );
    OB_ASSERT( doShuffleTest("Br.NCCCNCCCC[C@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82") );
    OB_ASSERT( doShuffleTest("Br.NCCCNCCCC[C@@-]12C3=C4C5=C1[Fe+2]16782345[C-]2C7=C6C1=C82") );
    OB_ASSERT( doShuffleTest("CC(=O)O.Oc1nc(N)c2nc3c(nc2n1)c1cccc2cccc3c12.OS(=O)(=O)O") );
    OB_ASSERT( doShuffleTest("C12=C3C4=C5[C@H]1[Fe]16782345C2=C6C7=C1C82") );
    OB_ASSERT( doShuffleTest("C12([C@@H]3CCC[C@@H]1CCC[C@@H]2CCC3)OC(=O)C") );
    OB_ASSERT( doShuffleTest("[Fe]12345678(C9(C1=C4C3=C29)C=O)C1=C6C7=C5C81") );
    OB_ASSERT( doShuffleTest("[Fe+2]12345678([C-]9(C1=C7C6=C29)CCCCNCCCCN)C1=C4C5=C3[C-]81.Br") );
    OB_ASSERT( doShuffleTest("O[C@H]1CC[C@@H](O)CC1") );
    OB_ASSERT( doShuffleTest("c1csc(c1)c1nnc(o1)C12C[C@@H]3C[C@@H](C[C@@H](C3)C2)C1") );
    OB_ASSERT( doShuffleTest("CC(=O)OC12[C@@H]3CCC[C@@H]1CCC[C@@H]2CCC3") );
    OB_ASSERT( doShuffleTest("CC(=O)OC12[C@@H]3CCC[C@@H]1CCC[C@@H]2CCC3") );
    OB_ASSERT( doShuffleTest("C/C(=N/c1ccncc1)/C(=N/c1ccncc1)/C") );
    OB_ASSERT( doShuffleTest("CC(=CCN[C@]12C[C@@H]3C[C@@H](C[C@@H](C3)C2)C1)C") );
    OB_ASSERT( doShuffleTest("Brc1ccc(nc1)NC(=S)N[C@]12C[C@@H]3C[C@H](C[C@@H](C3)C2)C1") );
    OB_ASSERT( doShuffleTest("c1ccc2c(c1)[C@H]1[C@@H]3O[C@@H]3[C@@H]2c2ccccc12") );
    OB_ASSERT( doShuffleTest("C(C/N=C\\c1ccccc1)/N=C/c1ccccc1") );
    OB_ASSERT( doShuffleTest("O=C(CCCC(=O)N/N=C(/C)\\c1ccccc1)N/N=C(\\C)/c1ccccc1") );
    OB_ASSERT( doShuffleTest("C/C(=N/c1ccccc1)/C(=N/c1ccccc1)/C") );
    OB_ASSERT( doShuffleTest("OC[C@H](O)COP1(=O)O[Pt@]2(N[C@H]3CCCC[C@@H]3N2)O1") );
    OB_ASSERT( doShuffleTest("n1ccc(cc1)C12C[C@H]3C[C@H](C[C@H](C3)C2)C1") );
    OB_ASSERT( doShuffleTest("C12=C3C4=C5[C@H]1[Fe]16782345C2=C6C7=C1C82") );
    OB_ASSERT( doShuffleTest("CC(=O)O.Oc1nc2nc3c(nc2c(N)n1)c1cccc2cccc3c12") );
    OB_ASSERT( doShuffleTestFile("stereo/incorrect_canon1.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/incorrect_canon2.mol") );
    break;
  case 3:
    OB_ASSERT( doShuffleTest("O[C@H]1CC[C@@H](O)CC1") );
    OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@H](O)C1") );
    OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@@H](O)C1") );

    OB_ASSERT( doShuffleTest("[C@@H]1([C@H]([C@H]([C@H]1C)C)C)C") );
    OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D1.smi") );

    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_5_spec.mol") );

    OB_ASSERT( doShuffleTestFile("stereo/cyclohexanediol_D1.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/cyclohexanediol_D2.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/cyclohexanetriol_D1.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/cyclohexanetriol_D2.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D1.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D2.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D3.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D4.mol") );

    OB_ASSERT( doShuffleTest("[C@@H]1([C@H]([C@H]([C@H]1C)C)C)C") );
    OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@H]([C@H]1C)C)C)C") );	
    OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@H]([C@@H]1C)C)C)C") );	
    OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@@H]([C@H]1C)C)C)C") );	

    OB_ASSERT( doShuffleTestFile("stereo/inositol_cis.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/inositol_epi.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/inositol_allo.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/inositol_myo.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/inositol_muco.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/inositol_neo.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/inositol_scyllo.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/inositol_chiroD.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/inositol_chiroL.mol") );
    OB_ASSERT( doShuffleTest("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O") );
    break;
  case 4:
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec1.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec2.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec3.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec4.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec5.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_20_spec1.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_26_spec1.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_59_spec1.mol") );

    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_30.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_34.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_35.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_36.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_37.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_38.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_39.mol") );
    break;
  case 5:
    OB_ASSERT( doShuffleTestFile("stereo/canon1.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon2.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon3.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon4.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon5.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon6.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon7.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon8.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon9.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon10.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon11.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon12.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon13.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon14.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon15.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon16.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon17.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon18.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon19.mol") );
    OB_REQUIRE( doShuffleTestFile("stereo/canon19.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon20.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon21.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon22.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon23.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon24.mol") );

    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans1.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans2.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans3.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans4.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans5.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans6.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans7.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans8.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans9.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans15.mol") );
    OB_ASSERT( doShuffleTestFile("stereo/canon_cistrans16.mol") );
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}

