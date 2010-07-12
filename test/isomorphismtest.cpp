#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>

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
  OBIsomorphismMapper::Mappings maps = mapper->MapAll(&mol);

  OB_ASSERT( maps.size() == 4 );

  delete query;
  delete mapper;

  query = CompileSmilesQuery("C1(C)CCC(C)CC1");
  mapper = OBIsomorphismMapper::GetInstance(query);
  
  OB_ASSERT( mapper->MapFirst(&mol).size() == 8 );
  OB_ASSERT( mapper->MapUnique(&mol).size() == 1 );
  OB_ASSERT( mapper->MapAll(&mol).size() == 4 );

  delete query;
  delete mapper;
}

void testIsomorphism2()
{
  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, "Cc1ccc(C)cc1");

  OBQuery *query = CompileSmilesQuery("C1=CC=CC=C1");
  OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
  OBIsomorphismMapper::Mappings maps = mapper->MapUnique(&mol);

  cout << maps.size() << endl;

  OB_ASSERT( maps.size() == 1 );

  delete query;
  delete mapper;
}

int main() 
{
//  testIsomorphism1();
  testIsomorphism2();

  return 0;
}

                
