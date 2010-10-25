#include "obtest.h"
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;

std::string GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;
  return path;
}

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
int main(int argc, char **argv)
{
  // Define location of file formats for testing
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif  

  testAutomorphisms();

  OBMol mol;
  OBConversion conv;
  OB_ASSERT( conv.SetInFormat("mol") );

  /*
   * Computers & Chemistry 26 (2002) 119-123
   *
   * Figure 2. Test graphs
   */
  cout << "Hao, Xu paper, fig. 2: structure 1" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("hao_xu_1.mol")) );
  OB_ASSERT( doAutomorphismTest(mol, 8) );

  cout << "Hao, Xu paper, fig. 2: structure 2" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("hao_xu_2.mol")) );
  OB_ASSERT( doAutomorphismTest(mol, 2) );

  cout << "Hao, Xu paper, fig. 2: structure 3" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("hao_xu_3.mol")) );
  OB_ASSERT( doAutomorphismTest(mol, 48) );

  cout << "Hao, Xu paper, fig. 2: structure 4" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("hao_xu_4.mol")) );
  OB_ASSERT( doAutomorphismTest(mol, 2) );

  cout << "Hao, Xu paper, fig. 2: structure 5" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("hao_xu_5.mol")) );
  OB_ASSERT( doAutomorphismTest(mol, 2) );

  cout << "Hao, Xu paper, fig. 2: structure 6" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("hao_xu_6.mol")) );
  OB_ASSERT( doAutomorphismTest(mol, 6) );

  cout << "Hao, Xu paper, fig. 2: structure 7" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("hao_xu_7.mol")) );
  OB_ASSERT( doAutomorphismTest(mol, 1) );

  cout << "Hao, Xu paper, fig. 2: structure 8" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("hao_xu_8.mol")) );
  OB_ASSERT( doAutomorphismTest(mol, 1) );

  cout << "Hao, Xu paper, fig. 2: structure 9" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("hao_xu_9.mol")) );
  OB_ASSERT( doAutomorphismTest(mol, 20) );



}



