#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/rotor.h>
#include <openbabel/bond.h>
#include <openbabel/obutil.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace OpenBabel;

//! Comparison for doubles with a modulus: returns mod(a - b,m) < epsilon
bool IsNear_mod(const double &a, const double &b, const double &m, const double epsilon)
{
  double arg=a-b;
  while(arg<-m/2)
    arg+=m;
  while(arg>=m/2)
    arg-=m;

  return (fabs(arg) < epsilon);
}

void testOBRotorGetSet()
{
  OBBond bond;
  OBRotor rotor;
  
  // SetBond/GetBond
  rotor.SetBond(&bond);
  OB_ASSERT(rotor.GetBond()->GetIdx() == bond.GetIdx());
  // SetIdx/GetIdx
  rotor.SetIdx(45);
  OB_ASSERT(rotor.GetIdx() == 45);
  // SetDihedralAtoms/GetDihedralAtoms
  int ref_ptr[4] = {1, 2, 3, 4};
  rotor.SetDihedralAtoms(ref_ptr);
  rotor.GetDihedralAtoms(ref_ptr);
  OB_ASSERT(ref_ptr[0] == 1);
  OB_ASSERT(ref_ptr[1] == 2);
  OB_ASSERT(ref_ptr[2] == 3);
  OB_ASSERT(ref_ptr[3] == 4);
  std::vector<int> ref_vector = rotor.GetDihedralAtoms();
  OB_ASSERT(ref_vector[0] == 1);
  OB_ASSERT(ref_vector[1] == 2);
  OB_ASSERT(ref_vector[2] == 3);
  OB_ASSERT(ref_vector[3] == 4);
  rotor.SetDihedralAtoms(ref_vector);
  ref_vector = rotor.GetDihedralAtoms();
  OB_ASSERT(ref_vector[0] == 1);
  OB_ASSERT(ref_vector[1] == 2);
  OB_ASSERT(ref_vector[2] == 3);
  OB_ASSERT(ref_vector[3] == 4);
  // SetTorsionValues/GetTorsionValues/Size
  std::vector<double> angles;
  angles.push_back(0.0);
  angles.push_back(3.1415);
  rotor.SetTorsionValues(angles);
  OB_ASSERT(rotor.GetTorsionValues().size() == 2);
  OB_ASSERT(rotor.Size() == 2);
}

void testOBRotorSetToAngle()
{
  // 1 2 3 4 5 6 7 8
  // C-C-C-C-C-C-C-C
  //  0 1 2 3 4 5 6
  OBMolPtr mol = OBTestUtil::ReadFile("octane.cml");
  
  OBBond *bond = mol->GetBond(3);

  OBRotor rotor;
  // set the bond to rotate
  rotor.SetBond(bond);
  // set the dihedral indexes (these start from 1)
  int ref[4] = {3, 4, 5, 6};
  rotor.SetDihedralAtoms(ref);

  // find atoms to rotate
  std::vector<int> atoms;
  mol->FindChildren(atoms, 4, 5);
  // convert to coordinate indexes (start from 0, multiplied by 3)
  for (unsigned int i = 0; i < atoms.size(); ++i)
    atoms[i] = (atoms[i] - 1) * 3;
  rotor.SetRotAtoms(atoms);

  OB_ASSERT(IsNear_mod(fabs(RAD_TO_DEG * rotor.CalcTorsion(mol->GetCoordinates())), 180.0, 360.0, 1.0));

  // rotate
  rotor.SetToAngle(mol->GetCoordinates(), 60.0 * DEG_TO_RAD);

  OB_ASSERT(IsNear(RAD_TO_DEG * rotor.CalcTorsion(mol->GetCoordinates()), 60.0, 1.0));
}

void testOBRotorSetRotor()
{
  // 1 2 3 4 5 6 7 8
  // C-C-C-C-C-C-C-C
  //  0 1 2 3 4 5 6
  OBMolPtr mol = OBTestUtil::ReadFile("octane.cml");
  
  OBBond *bond = mol->GetBond(3);

  OBRotor rotor;
  // set the bond to rotate
  rotor.SetBond(bond);
  // set the dihedral indexes (these start from 1)
  int ref[4] = {3, 4, 5, 6};
  rotor.SetDihedralAtoms(ref);

  // find atoms to rotate
  std::vector<int> atoms;
  mol->FindChildren(atoms, 4, 5);
  // convert to coordinate indexes (start from 0, multiplied by 3)
  for (unsigned int i = 0; i < atoms.size(); ++i)
    atoms[i] = (atoms[i] - 1) * 3;
  rotor.SetRotAtoms(atoms);

  OB_ASSERT(IsNear_mod(fabs(RAD_TO_DEG * rotor.CalcTorsion(mol->GetCoordinates())), 180.0, 360.0, 1.0));
  rotor.SetToAngle(mol->GetCoordinates(), 60.0 * DEG_TO_RAD);

  // set torsion values
  std::vector<double> angles;
  angles.push_back(0.0);
  angles.push_back(3.1415);
  rotor.SetTorsionValues(angles);

  // rotate to 0.0 radians
  rotor.SetRotor(mol->GetCoordinates(), 0);
  OB_ASSERT(IsNear(RAD_TO_DEG * rotor.CalcTorsion(mol->GetCoordinates()), 0.0, 1.0));
  // rotate to 3.1415 radians
  rotor.SetRotor(mol->GetCoordinates(), 1);
  OB_ASSERT(IsNear_mod(RAD_TO_DEG * rotor.CalcTorsion(mol->GetCoordinates()), 180.0, 360.0, 1.0));
   // rotate to 0.0 radians
  rotor.SetRotor(mol->GetCoordinates(), 0, 1);
  OB_ASSERT(IsNear(RAD_TO_DEG * rotor.CalcTorsion(mol->GetCoordinates()), 0.0, 1.0)); 
}


void testOBRotorListFixedBonds()
{
  // 1 2 3 4 5 6 7 8
  // C-C-C-C-C-C-C-C
  //  0 1 2 3 4 5 6
  OBMolPtr mol = OBTestUtil::ReadFile("octane.cml");

  // test with no bonds fixed
  OBRotorList rlist1;
  rlist1.Setup(*mol);
  OB_ASSERT(rlist1.Size() == 5);

  // test with bond 3 fixed
  OBBitVec fixedBonds;
  fixedBonds.SetBitOn(3);
  rlist1.SetFixedBonds(fixedBonds);
  rlist1.Setup(*mol);
  OB_ASSERT(rlist1.Size() == 4);

  // test with bond 1, 3, 5 fixed
  fixedBonds.SetBitOn(1);
  fixedBonds.SetBitOn(5);
  rlist1.SetFixedBonds(fixedBonds);
  rlist1.Setup(*mol);
  OB_ASSERT(rlist1.Size() == 2);

  // test with bond 1, 2, 3, 5 fixed
  fixedBonds.SetBitOn(2);
  rlist1.SetFixedBonds(fixedBonds);
  rlist1.Setup(*mol);
  OB_ASSERT(rlist1.Size() == 1);




}


int rotortest(int argc, char* argv[])
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
  // OBRotor
  case 1:
    testOBRotorGetSet();
    break;
  case 2:
    testOBRotorSetToAngle();
    break;
  case 3:
    testOBRotorSetRotor();
    break;
  // OBRotorList
  case 4:
    testOBRotorListFixedBonds();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}

