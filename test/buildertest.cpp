#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>
#include <openbabel/forcefield.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace OpenBabel;

std::string GetFilename(const std::string &filename)
{
  std::string path = TESTDATADIR + filename;
  return path;
}

static unsigned int failed = 0;
static unsigned int testCount = 0;

void clearMolFlags(OBMol &mol)
{
  // Both of these are needed or there's a crash.
  mol.UnsetFlag(OB_SSSR_MOL);
  mol.UnsetFlag(OB_AROMATIC_MOL);

  // Not needed, doesn't cause crash
  //  mol.UnsetFlag(OB_RINGFLAGS_MOL);
  //  mol.UnsetFlag(OB_LSSR_MOL);
  //  mol.UnsetFlag(OB_KEKULE_MOL);
  //  mol.UnsetFlag(OB_CLOSURE_MOL);
  //  mol.UnsetFlag(OB_AROM_CORRECTED_MOL);
  //  mol.UnsetFlag(OB_RINGTYPES_MOL);
}

bool doBuildMoleculeTest(OBMol &mol)
{
  testCount++;

  OBBuilder builder;
  OB_REQUIRE(builder.Build(mol));
  // Does not need clearMolFlags -- crash still happens if you clear here
  // and not after AddHydrogens()
  OB_REQUIRE(mol.AddHydrogens());
  clearMolFlags(mol); // must clear here or you crash

  OBForceField* pff = OBForceField::FindType("mmff94");
  cout << mol.GetTitle() << endl;
  OB_REQUIRE(pff->Setup(mol));

  return true;
}

bool doMultiMoleculeFile(const std::string &filename)
{
  cout << " Starting " << filename << endl;

  std::string file = GetFilename(filename);
  std::ifstream ifs;
  ifs.open(file.c_str());
  OB_REQUIRE( ifs );

  OBMol mol;
  OBConversion conv(&ifs, &cout);
  OBFormat *format = conv.FormatFromExt(file.c_str());
  OB_REQUIRE(format);
  OB_REQUIRE(conv.SetInFormat(format));
 
  bool result = true;
  while (conv.Read(&mol, &ifs)) {
    bool res = doBuildMoleculeTest(mol);
    if (!res)
      result = res;
  }

  return result;
}


int main(int argc, char **argv)
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif  

  OB_ASSERT( doMultiMoleculeFile("forcefield.sdf") );
  OB_ASSERT( doMultiMoleculeFile("filterset.sdf") );

  // fails because of selenium
  //  OB_ASSERT( doMultiMoleculeFile("aromatics.smi") );
  // fails because of stereo crash
  //  OB_ASSERT( doMultiMoleculeFile("nci.smi") );
  // fails because of "organometallic" entries
  //  OB_ASSERT( doMultiMoleculeFile("attype.00.smi") );

  cout << "PASSED TESTS: " << testCount - failed << "/" << testCount << endl;

  return 0;
}

