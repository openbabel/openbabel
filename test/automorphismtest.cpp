#include "obtest.h"
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;

bool doAutomorphismTest(OBMol &mol, int numAutomorphisms)
{
  OBIsomorphismMapper::Mappings G;
  FindAutomorphisms(&mol, G);

  return (G.size() == numAutomorphisms);
}

void testAutomorphisms()
{
  cout <<  "testAutomorphisms" << endl;
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "C1C(CC2CC2)C1");

  Automorphisms aut;
  FindAutomorphisms((OBMol*)&mol, aut);
  cout << aut.size() << endl;
  OB_ASSERT( aut.size() == 8 );
}

/**
 * Test detection of stereoisomers
 */

int automorphismtest(int argc, char* argv[])
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

  OBMol mol;
  OBConversion conv;
  OB_ASSERT( conv.SetInFormat("mol") );

  switch(choice) {
  case 1:
    testAutomorphisms();
    break;
  case 2:
    /*
     * Computers & Chemistry 26 (2002) 119-123
     *
     * Figure 2. Test graphs
     */
    cout << "Hao, Xu paper, fig. 2: structure 1" << endl;
    OB_ASSERT( conv.ReadFile(&mol, OBTestUtil::GetFilename("hao_xu_1.mol")) );
    OB_ASSERT( doAutomorphismTest(mol, 8) );
    break;
  case 3:
    cout << "Hao, Xu paper, fig. 2: structure 2" << endl;
    OB_ASSERT( conv.ReadFile(&mol, OBTestUtil::GetFilename("hao_xu_2.mol")) );
    OB_ASSERT( doAutomorphismTest(mol, 2) );
    break;
  case 4:
    cout << "Hao, Xu paper, fig. 2: structure 3" << endl;
    OB_ASSERT( conv.ReadFile(&mol, OBTestUtil::GetFilename("hao_xu_3.mol")) );
    OB_ASSERT( doAutomorphismTest(mol, 48) );
    break;
  case 5:
    cout << "Hao, Xu paper, fig. 2: structure 4" << endl;
    OB_ASSERT( conv.ReadFile(&mol, OBTestUtil::GetFilename("hao_xu_4.mol")) );
    OB_ASSERT( doAutomorphismTest(mol, 2) );
    break;
  case 6:
    cout << "Hao, Xu paper, fig. 2: structure 5" << endl;
    OB_ASSERT( conv.ReadFile(&mol, OBTestUtil::GetFilename("hao_xu_5.mol")) );
    OB_ASSERT( doAutomorphismTest(mol, 2) );
    break;
  case 7:
    cout << "Hao, Xu paper, fig. 2: structure 6" << endl;
    OB_ASSERT( conv.ReadFile(&mol, OBTestUtil::GetFilename("hao_xu_6.mol")) );
    OB_ASSERT( doAutomorphismTest(mol, 6) );
    break;
  case 8:
    cout << "Hao, Xu paper, fig. 2: structure 7" << endl;
    OB_ASSERT( conv.ReadFile(&mol, OBTestUtil::GetFilename("hao_xu_7.mol")) );
    OB_ASSERT( doAutomorphismTest(mol, 1) );
    break;
  case 9:
    cout << "Hao, Xu paper, fig. 2: structure 8" << endl;
    OB_ASSERT( conv.ReadFile(&mol, OBTestUtil::GetFilename("hao_xu_8.mol")) );
    OB_ASSERT( doAutomorphismTest(mol, 1) );
    break;
  case 10:
    cout << "Hao, Xu paper, fig. 2: structure 9" << endl;
    OB_ASSERT( conv.ReadFile(&mol, OBTestUtil::GetFilename("hao_xu_9.mol")) );
    OB_ASSERT( doAutomorphismTest(mol, 20) );
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }
  return 0;
}



