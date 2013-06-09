#include "obtest.h"
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;

typedef vector<vector3> vv3;

void testLossOfHydrogen(string filename)
{
  string testfile = OBTestUtil::GetFilename(filename);
  ifstream ifs(testfile.c_str());

  OBConversion conv(&ifs);
  OB_REQUIRE(conv.SetInFormat("sdf"));
  OBMol mol;
  OB_REQUIRE(conv.Read(&mol));
  bool success = true;
  int i = 0;
  while (success) {
    unsigned int Natoms = mol.NumAtoms();
    mol.DeleteHydrogens();
    mol.AddHydrogens();
    unsigned int newNatoms = mol.NumAtoms();
    cout << "Mol#" << i << ", Title " << mol.GetTitle() << ", Original atoms vs New atoms: ";
    cout << Natoms << " vs " << newNatoms << "\n";
    OB_ASSERT( Natoms == newNatoms);
    cout << "\n";
    success = conv.Read(&mol);
    i += 1;
  }

  
}

int implicitHtest(int argc, char* argv[])
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif  

  testLossOfHydrogen("implicitH.sdf");

  return 0;
}
