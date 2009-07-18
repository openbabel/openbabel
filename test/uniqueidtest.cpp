#include "obtest.h"

#include <openbabel/mol.h>

using namespace std;
using namespace OpenBabel;

// OBMol::NewAtom()
void testIdsNewAtom1()
{
  OBMol mol;
  for (int i = 0; i < 10; ++i) {
    OBAtom *atom = mol.NewAtom();
    OB_REQUIRE(atom->GetId() == i);
  }

  OB_REQUIRE( mol.GetAtomById(0) );
  OB_REQUIRE( mol.GetAtomById(4) );
  OB_REQUIRE( mol.GetAtomById(9) );
  OB_REQUIRE( !mol.GetAtomById(10) );

  OB_REQUIRE( mol.GetAtomById(0)->GetId() == 0 );
  OB_REQUIRE( mol.GetAtomById(4)->GetId() == 4 );
  OB_REQUIRE( mol.GetAtomById(9)->GetId() == 9 );
}

// OBMol::NewAtom(unsigned long id)
void testIdsNewAtom2()
{
  OBMol mol;
  for (int i = 0; i < 10; ++i) {
    OBAtom *atom = mol.NewAtom(i*2);
    OB_REQUIRE(atom->GetId() == i*2);
  }

  OB_REQUIRE( mol.GetAtomById(0) );
  OB_REQUIRE( !mol.GetAtomById(7) );
  OB_REQUIRE( mol.GetAtomById(8) );
  OB_REQUIRE( !mol.GetAtomById(9) );
  OB_REQUIRE( mol.GetAtomById(18) );
  OB_REQUIRE( !mol.GetAtomById(19) );

  OB_REQUIRE( mol.GetAtomById(0)->GetId() == 0 );
  OB_REQUIRE( mol.GetAtomById(8)->GetId() == 8 );
  OB_REQUIRE( mol.GetAtomById(18)->GetId() == 18 );
}

void testIdsDeleteAtom()
{
  OBMol mol;
  for (int i = 0; i < 10; ++i)
    mol.NewAtom();

  OB_REQUIRE( mol.GetAtomById(3) );
  OB_REQUIRE( mol.GetAtomById(4) );
  OB_REQUIRE( mol.GetAtomById(5) );

  mol.DeleteAtom(mol.GetAtomById(4));

  OB_REQUIRE( mol.GetAtomById(3) );
  OB_REQUIRE( mol.GetAtomById(3)->GetId() == 3 );
  OB_REQUIRE( !mol.GetAtomById(4) );
  OB_REQUIRE( mol.GetAtomById(5) );
  OB_REQUIRE( mol.GetAtomById(5)->GetId() == 5 );
}

void testIdsAddAtom()
{
  OBMol mol;
  // add 5 atoms
  for (int i = 0; i < 5; ++i)
    mol.NewAtom();

  OBAtom a;
  a.SetAtomicNum(6);
  // add a sixth atom
  mol.AddAtom(a);

  OB_REQUIRE( mol.NumAtoms() == 6 );
  OB_REQUIRE( mol.GetAtomById(5) );
  OB_REQUIRE( mol.GetAtomById(5)->GetId() == 5 );
}


int main() 
{
  testIdsNewAtom1();
  testIdsNewAtom2();
  testIdsDeleteAtom();
  testIdsAddAtom();
  
  return 0;
}

                
