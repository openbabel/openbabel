#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/generic.h>

#include <iostream>
#include <vector>
#include <algorithm>

using std::cout;
using std::endl;
using namespace OpenBabel;

/**
 * Check if writing the molecule from orig_smiles using frag_atoms results in ref_smiles.
 */
void test_smiles_fragment(const std::string &orig_smiles, const std::string &frag_atoms,
    const std::string &ref_smiles)
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("can");

  mol.Clear();
  OB_REQUIRE( conv.ReadString(&mol, orig_smiles) );
  OBPairData *pd = new OBPairData;
  pd->SetAttribute("SMILES_Fragment");
  pd->SetValue(frag_atoms);
  mol.SetData(pd);
  std::string smiles = conv.WriteString(&mol, true);
  OB_COMPARE(smiles, ref_smiles);
}

int canonfragmenttest(int argc, char *argv[])
{

  // Define location of file formats for testing
#ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
#endif

  test_smiles_fragment("c1ccccc1CC", "1 2 3 4 5 6", "c1ccccc1");
  test_smiles_fragment("CCc1cc(C)ccc1CC", "3 4 5 7 8 9", "c1ccccc1");
  test_smiles_fragment("CC.c1ccccc1.CC", "1 2 3 4 5 6 7 8", "c1ccccc1.CC");
  test_smiles_fragment("CC.c1ccccc1.CC", "3 4 5 6 7 8 9 10", "c1ccccc1.CC");
  test_smiles_fragment("[nH]1cccc1.c1ccccc1", "1 2 3 4 5", "c1ccc[nH]1");
  test_smiles_fragment("[nH]1cccc1.c1ccccc1", "6 7 8 9 10 11", "c1ccccc1");

  return 0;
}

