#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/atomclass.h>
#include <openbabel/obconversion.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace OpenBabel;


// A segfault was occuring when a Universal SMILES was output after an InChIfied SMILES.
// This was due to short-circuit caching of InChIs on reading. The fix was to limit
// the situations when the cached value was used, but also to delete the cached value
// in this particular instance.
void test_Issue135_UniversalSmiles()
{
  // Test writing U smiles after I smiles
  OBConversion conv;
  conv.SetInFormat("smi");
  OBMol mol;
  conv.ReadString(&mol, "C(=O)([O-])C(=O)O");
  conv.SetOutFormat("smi");
  conv.SetOptions("I", OBConversion::OUTOPTIONS);
  std::string res = conv.WriteString(&mol, true);
  OB_COMPARE(res, "C(=O)(C(=O)O)[O-]");
  conv.SetOptions("U", OBConversion::OUTOPTIONS);
  res = conv.WriteString(&mol, true);
  OB_COMPARE(res, "C(=O)(C(=O)[O-])O");
}

// Reading an InChI and then adding hydrogens messed up the structure
void test_Issue134_InChI_addH()
{
  OBConversion conv;
  conv.SetInFormat("inchi");
  OBMol mol;
  conv.ReadString(&mol, "InChI=1S/C2H7NO/c1-2(3)4/h2,4H,3H2,1H3/t2-/m0/s1");
  OB_ASSERT(!mol.HasData(OBGenericDataType::VirtualBondData));
  mol.AddHydrogens();
  conv.SetOutFormat("smi");
  std::string res = conv.WriteString(&mol, true);
  OB_COMPARE(res, "C[C@@H](N)O");
}

// Delete hydrogens should not remove charged or isotopic hydrogens or [H][H] or [Cu][H][Cu]
// or hydrogens with assigned atom classes
void test_Issue178_DeleteHydrogens()
{
  OBConversion conv;
  conv.SetInFormat("smi");
  OBMol mol;
  // Test DeleteHydrogens() and DeleteNonPolarHydrogens()
  static const char *smi[] = { "C[H]", "[H][H]", "C[1H]", "C[H]C", "C[H+]" };
  int numHs[] = { 0, 2, 1, 1, 1 };
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 2; ++j) {
      conv.ReadString(&mol, smi[i]);
      if (j == 0)
        mol.DeleteHydrogens();
      else
        mol.DeleteNonPolarHydrogens();
      int myNumHs = 0;
      FOR_ATOMS_OF_MOL(atom, mol)
        if (atom->IsHydrogen())
          myNumHs++;
      OB_COMPARE(myNumHs, numHs[i]);
    }
  }
  // Test DeletePolarHydrogens()
  static const char *smiB[] = { "N[H]", "[H][H]", "N[1H]", "N[H]C", "N[H+]" };
  int numHsB[] = { 0, 2, 1, 1, 1 };
  for (int i = 0; i < 5; ++i) {
    conv.ReadString(&mol, smiB[i]);
    mol.DeletePolarHydrogens();
    int myNumHs = 0;
    FOR_ATOMS_OF_MOL(atom, mol)
      if (atom->IsHydrogen())
        myNumHs++;
    OB_COMPARE(myNumHs, numHsB[i]);
  }
  // Test atom class
  // Currently, the SMILES parser does not retain atom classes for hydrogens on reading so...
  conv.ReadString(&mol, "C[H]");
  OBAtomClassData *ac = new OBAtomClassData;
  ac->Add(2, 99); // Assign the hydrogen (atom 2) a class of 99
  mol.SetData(ac);
  mol.DeleteHydrogens();
  int myNumHs = 0;
  FOR_ATOMS_OF_MOL(atom, mol)
    if (atom->IsHydrogen())
      myNumHs++;
  OB_COMPARE(myNumHs, 1);
}

void test_Issue305_NumRotors()
{
  OBMolPtr mol = OBTestUtil::ReadFile("regressiontest_numrotors.mol");
  OB_COMPARE(mol->NumRotors(), 9); // was returning 4
}

void test_PR329_Molfile_RGroups()
{
  // There are several way to get an R Group into a mol file
  // 1. The existing use of atom maps on dummy atoms in SMILES
  OBConversion conv;
  OBMol mol;
  conv.SetInAndOutFormats("smi", "mol");
  conv.ReadString(&mol, "C([*:1]CO[*:2]");
  obErrorLog.SetOutputLevel(obError); // avoid warning about no 2D or 3D coords
  std::string molfileWithRGP = conv.WriteString(&mol);
  obErrorLog.SetOutputLevel(obWarning);
  OB_ASSERT( molfileWithRGP.find("R#") != std::string::npos );
  OB_ASSERT( molfileWithRGP.find("M  RGP  2   2   1   5   2") != std::string::npos); // i.e. atom 2 is labelled R1, atom 5 is labelled R2
  // Check negative case
  conv.ReadString(&mol, "C([*]CO[*]");
  std::string molfileb = conv.WriteString(&mol);
  OB_ASSERT( molfileb.find("R#") == std::string::npos );
  OB_ASSERT( molfileb.find("M  RGP") == std::string::npos);

  // 2. By reading a molfile that use the R#, RGP notation
  conv.SetInAndOutFormats("mol", "mol");
  conv.ReadString(&mol, molfileWithRGP);
  molfileb = conv.WriteString(&mol);
  OB_ASSERT( molfileb.find("R#") != std::string::npos );
  OB_ASSERT( molfileb.find("M  RGP  2   2   1   5   2") != std::string::npos); // i.e. atom 2 is labelled R1, atom 5 is labelled R2

  // 3. By reading a molfile that specifies the atom alias as Rn, where n is an integer
  std::string molfileWithAlias = "\n"
" OpenBabel07211621152D\n"
"\n"
"  2  1  0  0  0  0  0  0  0  0999 V2000\n"
"    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
"    0.0000    1.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n"
"  1  2  1  0  0  0  0\n"
"A    2\n"
"R1\n"
"M  END";
  conv.SetInAndOutFormats("mol", "mol");
  conv.ReadString(&mol, molfileWithAlias);
  std::string molfile = conv.WriteString(&mol);
  OB_ASSERT( molfile.find("R#") != std::string::npos );
  OB_ASSERT( molfile.find("M  RGP  1   2   1") != std::string::npos); // i.e. atom 2 is labelled R1
  // Check negative case
  molfileWithAlias = "\n"
" OpenBabel07211621152D\n"
"\n"
"  2  1  0  0  0  0  0  0  0  0999 V2000\n"
"    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
"    0.0000    1.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n"
"  1  2  1  0  0  0  0\n"
"A    2\n"
"R\n"
"M  END";
  conv.SetInAndOutFormats("mol", "mol");
  obErrorLog.SetOutputLevel(obError); // avoid warning "Alias R was not chemically interpreted"
  conv.ReadString(&mol, molfileWithAlias);
  obErrorLog.SetOutputLevel(obWarning);
  molfile = conv.WriteString(&mol);
  OB_ASSERT( molfile.find("R#") == std::string::npos );
  OB_ASSERT( molfile.find("M  RGP") == std::string::npos);

  // 4. By reading a molfile that specifies the element name as R1, etc.
  std::string molfileWithRGroupElementName = "\n"
" OpenBabel07211621152D\n"
"\n"
"  2  1  0  0  0  0  0  0  0  0999 V2000\n"
"    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
"    0.0000    1.0000    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0\n"
"  1  2  1  0  0  0  0\n"
"M  END";
  conv.SetInAndOutFormats("mol", "mol");
  conv.ReadString(&mol, molfileWithRGroupElementName);
  molfile = conv.WriteString(&mol);
  OB_ASSERT( molfile.find("R#") != std::string::npos );
  OB_ASSERT( molfile.find("M  RGP  1   2   1") != std::string::npos); // i.e. atom 2 is labelled R1
}

int regressionstest(int argc, char* argv[])
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
    test_Issue135_UniversalSmiles();
    break;
  case 221:
    test_Issue134_InChI_addH();
    break;
  case 222:
    test_Issue178_DeleteHydrogens();
    break;
  case 223:
    test_Issue305_NumRotors();
    break;
  case 224:
    test_PR329_Molfile_RGroups();
    break;
    //case N:
  //  YOUR_TEST_HERE();
  //  Remember to update CMakeLists.txt with the number of your test
  //  break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}

