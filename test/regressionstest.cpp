#include "obtest.h"
#include <openbabel/mol.h>
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

