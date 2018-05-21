#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/phmodel.h>
#include <openbabel/elements.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace OpenBabel;

// PR #1831
void testMolToCdxmlConversion()
{
  OBConversion conv;
  OBMol mol;
  conv.SetInFormat("mol");
  conv.SetOutFormat("cdxml");
  conv.SetOutputIndex(1);

  conv.ReadFile(&mol, OBTestUtil::GetFilename("alanine.mol"));
  std::string cdxmlFromMol = conv.WriteString(&mol, true);

  std::string cdxmlTarget = OBTestUtil::ReadFileContent("alanine.cdxml");

  OB_COMPARE(cdxmlFromMol, cdxmlTarget);
}

int conversiontest(int argc, char* argv[])
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
    testMolToCdxmlConversion();
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
