#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

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


std::vector< std::vector<unsigned long> > getIdRingPaths(OBMol &mol)
{
  mol.UnsetFlag(OB_LSSR_MOL);
  mol.DeleteData("LSSR");
  std::vector<OBRing*> lssr = mol.GetLSSR();

  std::vector< std::vector<unsigned long> > idPaths;

  for (unsigned int i = 0; i < lssr.size(); ++i) {
    OBRing *ring = lssr[i];
    std::vector<unsigned long> idPath;
    for (unsigned int j = 0; j < ring->_path.size(); ++j) {
      idPath.push_back(mol.GetAtom(ring->_path[j])->GetId());
    }

    std::sort(idPath.begin(), idPath.end());  
    idPaths.push_back(idPath);
  }

  /*
  cout << "# idPaths = " << idPaths.size() << endl;
  for (unsigned int i = 0; i < idPaths.size(); ++i) {
    cout << "    ring: ";
    for (unsigned int j = 0; j < idPaths[i].size(); ++j) {
      cout << idPaths[i][j] << " ";    
    }
    cout << endl;
  }
  */

  return idPaths;
}


bool doShuffleTestMolecule(OBMol &mol)
{
  //cout << "-------------------------------------------------------------------------" << endl;
  int N = 100;
  testCount++;

  std::vector< std::vector<unsigned long> > ref = getIdRingPaths(mol);

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);
  
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get rings
    std::vector< std::vector<unsigned long> > rings = getIdRingPaths(mol);
    OB_ASSERT( rings.size() == ref.size() );
    if (rings.size() == ref.size()) {
      for (unsigned int j = 0; j < rings.size(); ++j) {
        bool found = false;
        for (unsigned int k = 0; k < ref.size(); ++k) {
          if (rings[j] == ref[k]) {
            found = true;
            break;
          }
        }
        OB_ASSERT( found );
        if (!found) {
          failed++;
          return false;
        }
      }
    }
  }

  return true;
}

bool doShuffleTestMultiFile(const std::string &filename)
{
  cout << "Shuffling: " << filename << endl;
  std::string file = GetFilename(filename);
  // read a smiles string
  OBMol mol;
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_REQUIRE( conv.SetInFormat(format) );

  testCount++;

  std::ifstream ifs;
  ifs.open(file.c_str());
  OB_REQUIRE( ifs );
 
  bool result = true;
  while (conv.Read(&mol, &ifs)) {
    bool res = doShuffleTestMolecule(mol);
    if (!res)
      result = res;
  }

  return result;
}

int main(int argc, char **argv)
{

  OB_ASSERT( doShuffleTestMultiFile("aromatics.smi") );
  OB_ASSERT( doShuffleTestMultiFile("nci.smi") );
  OB_ASSERT( doShuffleTestMultiFile("attype.00.smi") );

  //OB_ASSERT( doShuffleTest("") );

  cout << "PASSED TESTS: " << testCount - failed << "/" << testCount << endl;

  return 0;
}

