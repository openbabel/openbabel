#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/ring.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace OpenBabel;


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
      }
    }
  }

  return true;
}

bool doShuffleTestMultiFile(const std::string &filename)
{
  cout << "Shuffling: " << filename << endl;
  std::string file = OBTestUtil::GetFilename(filename);
  // read a smiles string
  OBMol mol;
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_REQUIRE( conv.SetInFormat(format) );

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

class LSSR 
{
  public:
    struct Size_Count 
    {
      Size_Count(int _ringSize, int _ringCount) : ringSize(_ringSize), ringCount(_ringCount) {}
      int ringSize;
      int ringCount;
    };

    std::vector<Size_Count> size_count;
  
    LSSR(const Size_Count &g1) 
    {
      size_count.push_back(g1);
    }
    LSSR(const Size_Count &g1, const Size_Count g2) 
    {
      size_count.push_back(g1);
      size_count.push_back(g2);
    }
    LSSR(const Size_Count &g1, const Size_Count g2, const Size_Count &g3)
    {
      size_count.push_back(g1);
      size_count.push_back(g2);
      size_count.push_back(g3);
    }


};

bool verifyLSSR(const std::string &filename, const LSSR &ref)
{
  cout << "Verify LSSR: " << filename << endl;
  std::string file = OBTestUtil::GetFilename(filename);
  // read a smiles string
  OBMol mol;
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_REQUIRE( conv.SetInFormat(format) );

  std::ifstream ifs;
  ifs.open(file.c_str());
  OB_REQUIRE( ifs );
  OB_REQUIRE( conv.Read(&mol, &ifs) );

  std::vector<int> ringSizeCount(20, 0); 
  std::vector<OBRing*> lssr = mol.GetLSSR();

  for (unsigned int i = 0; i < lssr.size(); ++i) {
    ringSizeCount[lssr[i]->_path.size()]++;
  }

  /*
  cout << "ringSize: ringCount" << endl;
  cout << "3: " << ringSizeCount[3] << endl;
  cout << "4: " << ringSizeCount[4] << endl;
  cout << "5: " << ringSizeCount[5] << endl;
  cout << "6: " << ringSizeCount[6] << endl;
  */

  bool fail = false;
  for (unsigned int i = 0; i < ref.size_count.size(); ++i) {
    const LSSR::Size_Count &size_count = ref.size_count[i];
    OB_ASSERT( ringSizeCount[size_count.ringSize] == size_count.ringCount );
  }

  return true;
}

int lssrtest(int argc, char* argv[])
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
    OB_ASSERT( doShuffleTestMultiFile("aromatics.smi") );
    break;
  case 2:
    OB_ASSERT( doShuffleTestMultiFile("nci.smi") );
    break;
  case 3:
    OB_ASSERT( doShuffleTestMultiFile("attype.00.smi") );
    break;
  case 4:    
    OB_ASSERT( doShuffleTestMultiFile("rings/tetrahedron.mdl") );
    OB_ASSERT( doShuffleTestMultiFile("rings/cubane.mdl") );
    OB_ASSERT( doShuffleTestMultiFile("rings/cubane2.mdl") );
    OB_ASSERT( doShuffleTestMultiFile("rings/octahedron.mdl") );
    OB_ASSERT( doShuffleTestMultiFile("rings/bridged1.mdl") );
    OB_ASSERT( doShuffleTestMultiFile("rings/fullerene20.mdl") );
    OB_ASSERT( doShuffleTestMultiFile("rings/fullerene60.mdl") );
    break;
  case 5:
    // 4x 3-ring, 1x 5-ring
    OB_ASSERT( verifyLSSR("rings/tetrahedron.mdl", LSSR(LSSR::Size_Count(3, 4), LSSR::Size_Count(5, 1))) );
    // 6x 4-ring, 1x 6-ring
    OB_ASSERT( verifyLSSR("rings/cubane.mdl", LSSR(LSSR::Size_Count(4, 6), LSSR::Size_Count(6, 1))) );
    OB_ASSERT( verifyLSSR("rings/cubane2.mdl", LSSR(LSSR::Size_Count(4, 6), LSSR::Size_Count(6, 1))) );
    // 8x 3-ring
    OB_ASSERT( verifyLSSR("rings/octahedron.mdl", LSSR(LSSR::Size_Count(3, 8))) );
    // 3x 6-ring
    OB_ASSERT( verifyLSSR("rings/bridged1.mdl", LSSR(LSSR::Size_Count(6, 3))) );
    // 12x 5-ring
    OB_ASSERT( verifyLSSR("rings/fullerene20.mdl", LSSR(LSSR::Size_Count(5, 12))) );
    // 12x 5-ring, 20x 6-ring
    OB_ASSERT( verifyLSSR("rings/fullerene60.mdl", LSSR(LSSR::Size_Count(5, 12), LSSR::Size_Count(6, 20))) );
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}

