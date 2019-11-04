#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>

using namespace std;
using namespace OpenBabel;

//////////////////////////////////////////////////////////////////////
//
//  Atom
//
//////////////////////////////////////////////////////////////////////

// OBMol::NewAtom()
void testIdsNewAtom1()
{
  OBMol mol;
  for (unsigned long i = 0; i < 10; ++i) {
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
  for (unsigned long i = 0; i < 10; ++i) {
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

//////////////////////////////////////////////////////////////////////
//
//  Bond
//
//////////////////////////////////////////////////////////////////////

// OBMol::NewBond()
void testIdsNewBond1()
{
  OBMol mol;
  for (unsigned long i = 0; i < 10; ++i) {
    OBBond *bond = mol.NewBond();
    OB_REQUIRE(bond->GetId() == i);
  }

  OB_REQUIRE( mol.GetBondById(0) );
  OB_REQUIRE( mol.GetBondById(4) );
  OB_REQUIRE( mol.GetBondById(9) );
  OB_REQUIRE( !mol.GetBondById(10) );

  OB_REQUIRE( mol.GetBondById(0)->GetId() == 0 );
  OB_REQUIRE( mol.GetBondById(4)->GetId() == 4 );
  OB_REQUIRE( mol.GetBondById(9)->GetId() == 9 );
}

// OBMol::NewBond(unsigned long id)
void testIdsNewBond2()
{
  OBMol mol;
  for (unsigned long i = 0; i < 10; ++i) {
    OBBond *bond = mol.NewBond(i*2);
    OB_REQUIRE(bond->GetId() == i*2);
  }

  OB_REQUIRE( mol.GetBondById(0) );
  OB_REQUIRE( !mol.GetBondById(7) );
  OB_REQUIRE( mol.GetBondById(8) );
  OB_REQUIRE( !mol.GetBondById(9) );
  OB_REQUIRE( mol.GetBondById(18) );
  OB_REQUIRE( !mol.GetBondById(19) );

  OB_REQUIRE( mol.GetBondById(0)->GetId() == 0 );
  OB_REQUIRE( mol.GetBondById(8)->GetId() == 8 );
  OB_REQUIRE( mol.GetBondById(18)->GetId() == 18 );
}

void testIdsDeleteBond()
{
  OBMol mol;
  for (int i = 0; i < 10; ++i)
    mol.NewBond();

  OB_REQUIRE( mol.GetBondById(3) );
  OB_REQUIRE( mol.GetBondById(4) );
  OB_REQUIRE( mol.GetBondById(5) );

  OBBond *bond4 = mol.GetBondById(4);
  OBAtom *a = mol.NewAtom();
  OBAtom *b = mol.NewAtom();
  bond4->SetBegin(a);
  bond4->SetEnd(b);
 
  mol.DeleteBond(mol.GetBondById(4));

  OB_REQUIRE( mol.GetBondById(3) );
  OB_REQUIRE( mol.GetBondById(3)->GetId() == 3 );
  OB_REQUIRE( !mol.GetBondById(4) );
  OB_REQUIRE( mol.GetBondById(5) );
  OB_REQUIRE( mol.GetBondById(5)->GetId() == 5 );
}

void testIdsAddBond()
{
  OBMol mol;
  // add 5 bonds
  for (int i = 0; i < 5; ++i)
    mol.NewBond();

  OBBond bond;
  OBAtom *a = mol.NewAtom();
  OBAtom *b = mol.NewAtom();
  bond.SetBegin(a);
  bond.SetEnd(b);
  // add a sixth bond
  mol.AddBond(bond);

  OB_REQUIRE( mol.NumBonds() == 6 );
  OB_REQUIRE( mol.GetBondById(5) );
  OB_REQUIRE( mol.GetBondById(5)->GetId() == 5 );
}

int uniqueidtest(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  switch(choice) {
  case 1:
    // atoms
    testIdsNewAtom1();
    testIdsNewAtom2();
    testIdsDeleteAtom();
    testIdsAddAtom();
    break;

  case 2:
    // bonds
    testIdsNewBond1();
    testIdsNewBond2();
    testIdsDeleteBond();
    testIdsAddBond();
    break;
  
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }
  
  return 0;
}

                
