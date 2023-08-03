#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>
#include <openbabel/elements.h>
#include <openbabel/obiter.h>

using namespace std;
using namespace OpenBabel;

void testIsomorphism1()
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "CC1CCC(C)CC1");

  OBQuery *query = CompileMoleculeQuery(&mol);
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
  OBIsomorphismMapper::Mappings maps;
  mapper->MapAll(&mol, maps);

  OB_ASSERT( maps.size() == 4 );

  delete query;
  delete mapper;

  query = CompileSmilesQuery("C1(C)CCC(C)CC1");
  mapper = OBIsomorphismMapper::GetInstance(query);

  OBIsomorphismMapper::Mapping map;
  mapper->MapFirst(&mol, map);
  OB_ASSERT( map.size() == 8 );

  mapper->MapUnique(&mol, maps);
  OB_ASSERT( maps.size() == 1 );

  mapper->MapAll(&mol, maps);
  OB_ASSERT( maps.size() == 4 );

  delete query;
  delete mapper;
}

void testIsomorphism2()
{
  cout << "testIsomorphism2" << endl;
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "Cc1ccc(C)cc1");

  OBQuery *query = CompileSmilesQuery("c1ccccc1");
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
  OBIsomorphismMapper::Mappings maps;
  mapper->MapUnique(&mol, maps);

  cout << maps.size() << endl;

  OB_ASSERT( maps.size() == 1 );

  delete query;
  delete mapper;
}

void testIsomorphism3()
{
  cout << "testIsomorphism3" << endl;
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "C1CC1CC1CC1");

  OBQuery *query = CompileSmilesQuery("C1CC1");
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
  OBIsomorphismMapper::Mappings maps;
  mapper->MapUnique(&mol, maps);

  cout << maps.size() << endl;

  OB_ASSERT( maps.size() == 2 );

  delete query;
  delete mapper;
}

void testIsomorphism4()
{
  cout << "testIsomorphism4" << endl;
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "C12C(C2)C1");

  OBQuery *query = CompileSmilesQuery("C1CC1");
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
  OBIsomorphismMapper::Mappings maps;
  mapper->MapUnique(&mol, maps);

  cout << maps.size() << endl;

  OB_ASSERT( maps.size() == 2 );

  delete query;
  delete mapper;
}

void testIsomorphismMask()
{
  // read file: 3 6-rings
  //
  //     /\ /\ /\
  //    |  |  |  |
  //     \/ \/ \/
  //
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("cml");
  std::ifstream ifs(OBTestUtil::GetFilename("isomorphism1.cml").c_str());
  OB_REQUIRE( ifs );
  conv.Read(&mol, &ifs);

  OBQuery *query = CompileSmilesQuery("C1CCCCC1");
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);

  // no mask
  OBIsomorphismMapper::Mappings maps;
  mapper->MapUnique(&mol, maps);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 3 );

  // mask first ring
  OBBitVec mask;
  for (int i = 0; i < 6; ++i)
    mask.SetBitOn(i+1);
  mapper->MapUnique(&mol, maps, mask);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 1 );

  // mask second ring also
  for (int i = 6; i < 10; ++i)
    mask.SetBitOn(i+1);
  mapper->MapUnique(&mol, maps, mask);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 2 );

  // just mask last ring (atomIds 7-8, 10-13)
  mask.Clear();
  for (int i = 10; i < 14; ++i)
    mask.SetBitOn(i+1);
  mask.SetBitOn(7 + 1); mask.SetBitOn(8 + 1);
  mapper->MapUnique(&mol, maps, mask);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 1 ); // Should be same result as masking just the first ring

  delete query;
  delete mapper;
}

void testAutomorphismMask() {
  // read file: 3 6-rings
  //
  //     /\ /\ /\
  //    |  |  |  |
  //     \/ \/ \/
  //
  cout <<  "testAutomorphismMask" << endl;
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("cml");
  std::ifstream ifs(OBTestUtil::GetFilename("isomorphism1.cml").c_str());
  OB_REQUIRE( ifs );
  conv.Read(&mol, &ifs);

  OBIsomorphismMapper::Mappings maps;

  // First of all, how many automorphisms are there without any mask?
  // This takes about 20 seconds, so you may want to comment this out while debugging
  FindAutomorphisms(&mol, maps);
  cout << maps.size() << endl;
  OB_ASSERT( maps.size() == 4 );

  // Now, let's remove the bridge (atomId 6) of the central ring.
  //
  //     /\ /\ /\
  //    |  |  |  |
  //     \/    \/
  // both rings can be flipped around exocyclic bond, the whole molecule can be mirrored
  // horizontally, this results in 2 x 2 x 2 = 8 automorphisms
  OBBitVec mask;
  mask.SetRangeOn(1, mol.NumAtoms());
  mask.SetBitOff(6+1);
  FindAutomorphisms(&mol, maps, mask);
  cout << maps.size() << endl;
  for (unsigned int i = 0; i < maps.size(); ++i) {
    OBIsomorphismMapper::Mapping::const_iterator j;
    for (j = maps[i].begin(); j != maps[i].end(); ++j)
      cout << j->second << " ";
    cout << endl;
  }
  OB_ASSERT( maps.size() == 8 );

  // Verify that atom Id 6 does not occur anywhere in the mappings
  OBIsomorphismMapper::Mappings::const_iterator a;
  OBIsomorphismMapper::Mapping::const_iterator b;
  for (a = maps.begin(); a != maps.end(); ++a)
    for (b = a->begin(); b!= a->end(); ++b) {
      OB_ASSERT( b->first != 6 );
      OB_ASSERT( b->second != 6 );
    }
}

void testAutomorphismMask2()
{
  // The test molecule is progesterone, a steroid (four fused non-planar rings)
  cout <<  "testAutomorphismMask2" << endl;
  OBMol mol;
  OBConversion conv;

  conv.SetInFormat("sdf");
  std::ifstream ifs(OBTestUtil::GetFilename("progesterone.sdf").c_str());
  OB_REQUIRE( ifs );
  OB_REQUIRE( conv.Read(&mol, &ifs) );

  Automorphisms _aut;
  OBBitVec _frag_atoms;
  FOR_ATOMS_OF_MOL(a, mol) {
    if (a->GetAtomicNum() != OBElements::Hydrogen)
      _frag_atoms.SetBitOn(a->GetIdx());
  }
  FindAutomorphisms((OBMol*)&mol, _aut, _frag_atoms);
  OB_ASSERT( _aut.size() == 1 );

}

void testAutomorphismPreMapping()
{
  cout <<  "testAutomorphismPreMapping" << endl;
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "c1(C)c(C)c(C)c(C)c(C)c1");

  Automorphisms aut;
  FindAutomorphisms((OBMol*)&mol, aut);
  cout << aut.size() << endl;
  OB_ASSERT( aut.size() == 2 );
}

// https://github.com/openbabel/openbabel/issues/1929
void testIsomorphism9()
{
  cout << "testIsomorphism9" << endl;
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "Cc1onc(C)c1");

  OBQuery *query = CompileMoleculeQuery(&mol);
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
  OBIsomorphismMapper::Mappings maps;
  mapper->MapAll(&mol, maps);

  OB_ASSERT(maps.size() == 1);

  OBIsomorphismMapper::Mapping map;
  OBIsomorphismMapper::Mapping::const_iterator iter;
  map = maps[0];
  for (iter=map.begin(); iter!=map.end(); ++iter)
    OB_ASSERT( iter->first == iter->second);

  delete query;
  delete mapper;
}

int isomorphismtest(int argc, char* argv[])
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
    testIsomorphism1();
    break;
  case 2:
    testIsomorphism2();
    break;
  case 3:
    testIsomorphism3();
    break;
  case 4:
    testIsomorphism4();
    break;
  case 5:
    testIsomorphismMask();
    break;
  case 6:
    testAutomorphismMask();
    break;
  case 7:
    testAutomorphismMask2();
    break;
  case 8:
    testAutomorphismPreMapping();
    break;
  case 9:
    testIsomorphism9();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }
  return 0;
}


