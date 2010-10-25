#include "obbench.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

std::string GetFilename(const std::string &filename)
{
  std::string path = TESTDATADIR + filename;
  return path;
}

using namespace OpenBabel;

void benchmarkOBMol1()
{
  OB_NAMED_BENCHMARK("OBMol 1: create/destroy OBMol with 10000 atoms") {
    OBMol mol;
    for (unsigned int i = 0; i < 10000; ++i)
      mol.NewAtom();  
  }
}

void benchmarkOBMol2()
{
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("pdb") );
  OB_NAMED_BENCHMARK("OBMol 2: reading pdb file with 1788 atoms") {
    OBMol mol;
    conv.ReadFile(&mol, GetFilename("1DRF.pdb"));
  }
}

void benchmarkOBMol3()
{
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("pdb") );
  OB_NAMED_BENCHMARK("OBMol 2: reading pdb file with 18448 atoms") {
    OBMol mol;
    conv.ReadFile(&mol, GetFilename("3G61.pdb"));
  }
}

int main()
{
  benchmarkOBMol1();
  benchmarkOBMol2();
  benchmarkOBMol3();
}
