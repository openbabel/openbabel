/**********************************************************************
aliastest.cpp - unit testing code for alias.h and alias.cpp
Copyright (C) 2019 by Alex Ustinov

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
#include <openbabel/alias.h>

using namespace std;
using namespace OpenBabel;

// a class for a test case should contain smiles of the molecule, number of aliases that
// are expected to be produced, and number of irreducible atoms that are not parts of aliases
// note this test uses standard superatom.txt file shipped with OpenBabel
/*
CO2Et    EtO2C    C(=O)OCC
COOEt    EtOOC    C(=O)OCC
OiBu     iBuO     OCC(C)C
tBu      tBu      C(C)(C)C
nBu      nBu      CCCC
iPr      iPr      C(C)C
nPr      nPr      CCC
Et       Et       CC
NCF3     F3CN     NC(F)(F)F
CF3      F3C      C(F)(F)F
CCl3     Cl3C     C(Cl)(Cl)Cl
CN       NC       C#N
NC       CN       [N+]#[C-]
N(OH)CH3 CH3(OH)N N(O)C
NO2      O2N      N(=O)=O
NO2    O2N      [N+]([O-])=O
NO       ON       N=O
SO3H     HO3S     S(=O)(=O)O
COOH     HOOC     C(=O)O                blue
OEt      EtO      OCC
OAc      AcO      OC(=O)C
NHAc     AcNH     NC(=O)C
Ac       Ac       C(=O)C
CHO      OHC      C=O
NMe      MeN      NC
SMe      MeS      SC
OMe      MeO      OC
COO-     -OOC     C(=O)[O-]
*/

class AliasTestExample {
  string _smiles;
  unsigned int _num_aliases;
  unsigned int _num_irreducibles;
public:
  AliasTestExample(const string smiles, const unsigned int num_aliases, const unsigned int num_irreducibles):
    _smiles(smiles), _num_aliases(num_aliases), _num_irreducibles(num_irreducibles) {};
  void test() {
    OBMol mol;
    AliasData ad;
    OBConversion conv;
    cout << _smiles << endl;
    OB_REQUIRE( conv.SetInFormat("smi") );
    OB_REQUIRE( conv.ReadString(&mol, _smiles) );
    cout << mol.NumAtoms() << endl;
    OB_ASSERT( ad.AddAliases(&mol) == (_num_aliases != 0) );
    AliasData::RevertToAliasForm(mol);
    OB_ASSERT( mol.NumAtoms() == _num_irreducibles );
    cout << mol.NumAtoms() << endl;
    OB_REQUIRE( conv.SetOutFormat("smi") );
    cout << conv.WriteString(&mol) << endl;
  }
};

void testAliases()
{
  AliasTestExample test_set[] = {
    // methyl tert-butyl ether -> COtBu
    AliasTestExample("COC(C)(C)C", 1, 2),
    // diethylmalonate -> EtO2CCCO2Et
    AliasTestExample("CCOC(=O)CC(=O)OCC", 2, 1),
    // thioanisole -> PhSMe
    AliasTestExample("c1ccccc1SC", 1, 2)
  };

  for (auto i : test_set) {
    i.test();
  }
  cout << endl;
}

int aliastest(int argc, char* argv[])
{
  int defaultchoice = 1;

  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }
  switch(choice) {
  case 1:
    testAliases();
    break;
  case 2:
    break;
  case 3:
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}
