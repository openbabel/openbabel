#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>

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
  OBBitVec allbits(mol1.NumAtoms());
  FOR_ATOMS_OF_MOL(a, mol1) {
    allbits.SetBitOn(a->GetIdx());
  }

  CanonicalLabels(&mol1, allbits, symclasses1, canlbls1); 
  cout << "mol1.NumAtoms = " << mol1.NumAtoms() << endl;
    
  // write to can smiles
  std::string canSmiles = conv.WriteString(&mol1);
  cout << "canSmiles: " << canSmiles;
  // read can smiles in again
  OB_REQUIRE( conv.ReadString(&mol2, canSmiles) );

  CanonicalLabels(&mol2, allbits, symclasses2, canlbls2); 
  cout << "mol2.NumAtoms = " << mol2.NumAtoms() << endl;
 
  std::vector<unsigned int> symclassesCopy1 = symclasses1;
  std::vector<unsigned int>::iterator end1 = std::unique(symclassesCopy1.begin(), symclassesCopy1.end());
  unsigned int unique1 = end1 - symclassesCopy1.begin();
  
  std::vector<unsigned int> symclassesCopy2 = symclasses2;
  std::vector<unsigned int>::iterator end2 = std::unique(symclassesCopy2.begin(), symclassesCopy2.end());
  unsigned int unique2 = end2 - symclassesCopy2.begin();

  std::cout << "unique1 = " << unique1 << std::endl;
  std::cout << "unique2 = " << unique2 << std::endl;
  OB_ASSERT( unique1 == unique2 );

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
    OB_ASSERT( a1->GetValence() == a2->GetValence() );
    OB_ASSERT( a1->GetHvyValence() == a2->GetHvyValence() );
    OB_ASSERT( a1->GetHeteroValence() == a2->GetHeteroValence() );
    OB_ASSERT( a1->GetImplicitValence() == a2->GetImplicitValence() );
  }

  cout << "." << endl << endl; 
}


int main() 
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

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
  
  // FAILING:

  // ring gets converted to aromatic ring, adding H on n (i.e. N -> [nH])
  //genericGraphSymTest("CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)CNCC3=CC=CC=C3)O");

  //genericGraphSymTest("CC(C)[C@H]1CC[C@]([C@@H]2[C@@H]1C=C(COC2=O)C(=O)O)(CCl)O");
  
  //genericGraphSymTest("CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2");

  cout << "end" << endl;

  return 0;
}

                
