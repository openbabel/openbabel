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
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include <openbabel/alias.h>

using namespace std;
using namespace OpenBabel;

// a class for a test case should contain smiles of the molecule, number of aliases that
// are expected to be produced, and number of unaliasable atoms that are not parts of aliases
// note this test uses standard superatom.txt file included in OpenBabel

class AliasTestExample {
  string _smiles;
  unsigned int _num_aliases;
  unsigned int _num_nonaliased;
public:
  AliasTestExample(const string smiles, const unsigned int num_aliases, const unsigned int num_nonaliased):
    _smiles(smiles), _num_aliases(num_aliases), _num_nonaliased(num_nonaliased) {};

  // this will create aliases in the molecule, then remove all atoms present in aliases,
  // count aliases and atoms that are not present in aliases, and compare with the expected result
  void test() {
    OBMol mol;
    AliasData ad;
    OBConversion conv;
    OB_REQUIRE( conv.SetInFormat("smi") );
    OB_REQUIRE( conv.ReadString(&mol, _smiles) );
    ad.AddAliases(&mol);
    AliasData::RevertToAliasForm(mol);
    unsigned int alias_count = 0;
    unsigned int nonaliased_count = 0;

    FOR_ATOMS_OF_MOL(a, mol)
    {
      AliasData* ad = nullptr;
      if ( static_cast<AliasData*>(a->GetData(AliasDataType)) )
        alias_count++;
      else
        nonaliased_count++;
    }
    //cout << "Testing smiles " << _smiles << endl;
    //cout << "number of aliases " << alias_count << ", number of nonaliased atoms " << nonaliased_count << endl;

    OB_ASSERT( nonaliased_count == _num_nonaliased );
    OB_ASSERT( alias_count == _num_aliases );
  }
};

void testAliases()
{
  AliasTestExample test_set[] = {
    // methyl tert-butyl ether -> MeOtBu
    AliasTestExample("COC(C)(C)C", 2, 0),
    // diethylmalonate -> EtO2CCCO2Et
    AliasTestExample("CCOC(=O)CC(=O)OCC", 2, 1),
    // thioanisole -> PhSMe
    AliasTestExample("c1ccccc1SC", 1, 6),
    // oxalate dianion -> -O2CCO2-
    AliasTestExample("[O-]C(=O)C(=O)[O-]", 2, 0),
    // diphenyl - nothing to alias
    AliasTestExample("c1ccccc1c2ccccc2", 0, 12)
  };

  for (auto i : test_set) {
    i.test();
  }
  cout << endl;
}

int aliastest(int argc, char* argv[])
{
  // Define location of file formats for testing
#ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
#endif

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
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}
