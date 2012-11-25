#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;

/*
void dumpJSON(OBMol &mol)
{
  mol.GetAtom(1)->GetHyb();
  mol.GetAtom(1)->GetType();
  FOR_ATOMS_OF_MOL (atom, mol)
    std::cout << "    " << atom->JSON() << std::endl;
}

void compareJSON(OBMol &mol1, OBMol &mol2)
{
  for (std::size_t i = 0; i < mol1.NumAtoms(); ++i) {
    std::string json1 = mol1.GetAtom(i+1)->JSON();
    std::string json2 = mol2.GetAtom(i+1)->JSON();
    OB_COMPARE(json1.substr(0, json1.find(", flags")) + " }", json2.substr(0, json2.find(", flags")) + " }");
  }
}
*/


void test_parser(const std::string &smiles)
{
  std::cout << "Testing: " << smiles << std::endl;
  OBConversion obConv, smileyConv;

  OB_REQUIRE(obConv.SetInFormat("smi"));
  OB_REQUIRE(obConv.SetOutFormat("smi"));
  OB_REQUIRE(smileyConv.SetInFormat("smy"));

  OBMol obMol, smileyMol;

  OB_REQUIRE(obConv.ReadString(&obMol, smiles));
  OB_REQUIRE(smileyConv.ReadString(&smileyMol, smiles));

  std::string obCanSmiles = obConv.WriteString(&obMol, true);
  std::string smileyCanSmiles = obConv.WriteString(&smileyMol, true);

  /*
  std::cout << "obMol:" << std::endl;
  dumpJSON(obMol);
  std::cout << "smileyMol:" << std::endl;
  dumpJSON(smileyMol);

  compareJSON(obMol, smileyMol);
  */

  OB_COMPARE(obCanSmiles, smileyCanSmiles);
}

int main()
{
  // simple test...
  test_parser("CCC");

  // test if aromaticity is working
  test_parser("c1ccccc1");
  test_parser("[nH]1cccc1");
  test_parser("c1ncncc1");

  // tetrahedral stereochemistry
  test_parser("C[C@](F)(Cl)Br");
  test_parser("C1CCCC[C@H]1[C@@H]1CCCCC1");
  test_parser("O=C1OC(=O)[C@@H]2[C@H]1C1C=CC2C1");
  test_parser("CO[C@@H](C)[C@@H](N)C(=O)O");

  // cis/trans stereochemistry
  test_parser("F/C=C/F");
  test_parser("F/C=C\\F");
  test_parser("F\\C=C/F");
  test_parser("F\\C=C\\F");
  test_parser("O=C1CCCC/C/1=N\\Nc1ccccc1");

  // mixed stereochemistry
  test_parser("CC(C)CCC[C@@H](C)[C@H]1CCC2C3CC=C4CC(CC[C@@]4(C)C3CC[C@@]12C)OC(=O)/C=C/c1ccccc1");

  // OB's smiles parser assigns hybridization 2 to 2nd O???
  test_parser("O=Nc1ccc(O)cc1");


  test_parser("Cn1nccc1");

  test_parser("NC1CC(C)(C)N([O])C(C)(C)C1");

  return 0;

  std::ifstream ifs("emolecules1M.smi");
  std::string line;
  while (std::getline(ifs, line))
    test_parser(line);

  return 0;
}
