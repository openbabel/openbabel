#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/graphsym.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

#include <openbabel/canon.h>

using namespace std;
using namespace OpenBabel;

/*
 * Stereo classes have their own tests. This file tests if the smiles
 * format uses them correctly.
 */

void genericGraphSymTest(const std::string &smiles)
{
  cout << "Testing generic smiles <-> canonical smiles" << endl;
  // read a smiles string
  OBMol mol1, mol2;
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("smi") );
  OB_REQUIRE( conv.SetOutFormat("can") );
  cout << "smiles: " << smiles << endl;
  // read a smiles string
  OB_REQUIRE( conv.ReadString(&mol1, smiles) );


  std::vector<unsigned int> canlbls1, canlbls2;
  std::vector<unsigned int> symclasses1, symclasses2;

  OBGraphSym gs1(&mol1);
  gs1.GetSymmetry(symclasses1);
  CanonicalLabels(&mol1, symclasses1, canlbls1);
  cout << "mol1.NumAtoms = " << mol1.NumAtoms() << endl;

  // write to can smiles
  std::string canSmiles = conv.WriteString(&mol1);
  cout << "canSmiles: " << canSmiles;
  // read can smiles in again
  OB_REQUIRE( conv.ReadString(&mol2, canSmiles) );

  OBGraphSym gs2(&mol2);
  gs2.GetSymmetry(symclasses2);
  CanonicalLabels(&mol2, symclasses2, canlbls2);
  cout << "mol2.NumAtoms = " << mol2.NumAtoms() << endl;

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
        a2 = mol2.GetAtom(i+1);
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

  cout << "." << endl << endl;
}

void countGraphSymClassesTest(const std::string &filename, int numberOfClasses)
{
  cout << filename << endl;
  std::string file = OBTestUtil::GetFilename(filename);
  OBMol mol;
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_REQUIRE( conv.SetInFormat(format) );
  OB_REQUIRE( conv.ReadFile(&mol, file) );

  OBGraphSym graphSym(&mol);
  std::vector<unsigned int> symmetry_classes;
  graphSym.GetSymmetry(symmetry_classes);
  for (unsigned int i = 0; i < symmetry_classes.size(); ++i)
    cout << i+1 << ": " << symmetry_classes[i] << endl;

  std::sort(symmetry_classes.begin(), symmetry_classes.end());
  std::vector<unsigned int>::iterator end = std::unique(symmetry_classes.begin(), symmetry_classes.end());
  unsigned int n = end - symmetry_classes.begin();

  OB_ASSERT( n == numberOfClasses);
  if (n != numberOfClasses) {
    cout << n << " == " << numberOfClasses << endl;
  }
}

int graphsymtest(int argc, char* argv[])
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

  switch(choice) {
  case 1:
    genericGraphSymTest("C[C@H](O)N");
    genericGraphSymTest("Cl[C@@](CCl)(I)Br");
    genericGraphSymTest("Cl/C=C/F");
    genericGraphSymTest("CCC[C@@H](O)CC\\C=C\\C=C\\C#CC#C\\C=C\\CO");
    genericGraphSymTest("O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5");
    genericGraphSymTest("OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1");
    genericGraphSymTest("OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2");
    genericGraphSymTest("CC(=O)OCCC(/C)=C\\C[C@H](C(C)=C)CCC=C");
    genericGraphSymTest("CC[C@H](O1)CC[C@@]12CCCO2");
    genericGraphSymTest("CN1CCC[C@H]1c2cccnc2");
    genericGraphSymTest("C(CS[14CH2][14C@@H]1[14C@H]([14C@H]([14CH](O1)O)O)O)[C@@H](C(=O)O)N");
    genericGraphSymTest("CCC[C@@H]1C[C@H](N(C1)C)C(=O)NC([C@@H]2[C@@H]([C@@H]([C@H]([C@H](O2)SC)OP(=O)(O)O)O)O)C(C)Cl");
    genericGraphSymTest("CC(C)[C@H]1CC[C@]([C@@H]2[C@@H]1C=C(COC2=O)C(=O)O)(CCl)O");
    genericGraphSymTest("CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2");
    break;

  case 2:
    // ring gets converted to aromatic ring, adding H on n (i.e. N -> [nH])
    //genericGraphSymTest("CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)CNCC3=CC=CC=C3)O");
    // This is the aromatic form -- GRH, it passes
    genericGraphSymTest("Cc1cn(c(=O)[nH]c1=O)[C@H]1C[C@@H]([C@H](O1)CNCc1ccccc1)O");
    break;

  case 3:
    countGraphSymClassesTest("stereo/razinger_fig3.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_1.mol", 5);
    countGraphSymClassesTest("stereo/razinger_fig7_2.mol", 7);
    countGraphSymClassesTest("stereo/razinger_fig7_3.mol", 7);
    countGraphSymClassesTest("stereo/razinger_fig7_4.mol", 7);
    countGraphSymClassesTest("stereo/razinger_fig7_5.mol", 5);
    countGraphSymClassesTest("stereo/razinger_fig7_6.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_7.mol", 7);
    countGraphSymClassesTest("stereo/razinger_fig7_8.mol", 7);
    countGraphSymClassesTest("stereo/razinger_fig7_9.mol", 10);
    countGraphSymClassesTest("stereo/razinger_fig7_10.mol", 13);
    countGraphSymClassesTest("stereo/razinger_fig7_11.mol", 9);
    countGraphSymClassesTest("stereo/razinger_fig7_12.mol", 17);
    countGraphSymClassesTest("stereo/razinger_fig7_13.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_14.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_15.mol", 5);
    countGraphSymClassesTest("stereo/razinger_fig7_16.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_17.mol", 5);
    countGraphSymClassesTest("stereo/razinger_fig7_18.mol", 3);
    countGraphSymClassesTest("stereo/razinger_fig7_19.mol", 6);
    break;
  case 4:
    countGraphSymClassesTest("stereo/razinger_fig7_20.mol", 5);
    countGraphSymClassesTest("stereo/razinger_fig7_21.mol", 7);
    countGraphSymClassesTest("stereo/razinger_fig7_22.mol", 2);
    countGraphSymClassesTest("stereo/razinger_fig7_23.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_24.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_25.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_26.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_27.mol", 8);
    countGraphSymClassesTest("stereo/razinger_fig7_28.mol", 3);
    countGraphSymClassesTest("stereo/razinger_fig7_29.mol", 6);
    break;
    /*
    countGraphSymClassesTest("stereo/razinger_fig7_30.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_31.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_32.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_33.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_34.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_35.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_36.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_37.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_38.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_39.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_40.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_41.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_42.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_43.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_44.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_45.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_46.mol", );
    countGraphSymClassesTest("stereo/razinger_fig7_47.mol", );
    */
  case 5:
    countGraphSymClassesTest("stereo/razinger_fig7_48.mol", 3);
    countGraphSymClassesTest("stereo/razinger_fig7_49.mol", 11);
    countGraphSymClassesTest("stereo/razinger_fig7_50.mol", 5);
    countGraphSymClassesTest("stereo/razinger_fig7_51.mol", 7);
    countGraphSymClassesTest("stereo/razinger_fig7_52.mol", 4);
    //countGraphSymClassesTest("stereo/razinger_fig7_53.mol", 3); // missing
    //countGraphSymClassesTest("stereo/razinger_fig7_54.mol", 5); // missing
    //countGraphSymClassesTest("stereo/razinger_fig7_55.mol", 9); // missing
    //countGraphSymClassesTest("stereo/razinger_fig7_56.mol", 9); // missing
    //countGraphSymClassesTest("stereo/razinger_fig7_57.mol", ); // missing
    countGraphSymClassesTest("stereo/razinger_fig7_58.mol", 5);
    countGraphSymClassesTest("stereo/razinger_fig7_59.mol", 4);
    countGraphSymClassesTest("stereo/razinger_fig7_60.mol", 7);
    //countGraphSymClassesTest("stereo/razinger_fig7_61.mol", ); // missing
    //countGraphSymClassesTest("stereo/razinger_fig7_62.mol", ); // missing
    //countGraphSymClassesTest("stereo/razinger_fig7_63.mol", ); // missing
    countGraphSymClassesTest("stereo/razinger_fig7_64.mol", 8);
    //countGraphSymClassesTest("stereo/razinger_fig7_65.mol", 2); // missing
    //countGraphSymClassesTest("stereo/razinger_fig7_66.mol", 2); // missing
    //countGraphSymClassesTest("stereo/razinger_fig7_67.mol", 3); // missing
    //countGraphSymClassesTest("stereo/razinger_fig7_68.mol", 3); // missing
    countGraphSymClassesTest("stereo/razinger_fig7_69.mol", 2);
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }
  return 0;
}


