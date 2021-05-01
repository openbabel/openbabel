/**********************************************************************
periodictest.cpp - Unit tests to check implementation of periodic boundary
conditions via the minimum image convention.

Copyright (C) 2018 by Ben Bucior

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "obtest.h"
#include <openbabel/babelconfig.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/generic.h>
#include <openbabel/obconversion.h>

#include <string>
#include <algorithm>

using namespace std;
using namespace OpenBabel;


class PeriodicTester {
public:
  PeriodicTester();
  void TestLengths(double a, double b, double c);
  void TestAngles(double a, double b);
  void TestTorsion(double a);
  void MakePeriodic(double a = 10.0);
  OBMol* GetMol() { return &tmol; }
  bool near(double a, double b, double tol = 0.001) {
    return (fabs(a - b) < tol);
  }

private:
  OBMol tmol;
  std::vector<OBAtom*> atom_list;  // atoms in the consistent order
};


PeriodicTester::PeriodicTester() {
  // Builds a made-up test molecule that straddles periodic boundaries
  OBAtom* a;
  OBBond* b;

  tmol.BeginModify();
  a = tmol.NewAtom();
  a->SetVector(3.0, 3.0, 1.0);
  a->SetAtomicNum(8);
  atom_list.push_back(a);
  a = tmol.NewAtom();
  a->SetVector(1.0, 1.0, 1.0);
  a->SetAtomicNum(7);
  atom_list.push_back(a);
  a = tmol.NewAtom();
  a->SetVector(-1.0, 1.0, 1.0);
  a->SetAtomicNum(6);
  atom_list.push_back(a);
  a = tmol.NewAtom();
  a->SetVector(-1.0, 1.0, -1.0);
  a->SetAtomicNum(35);
  atom_list.push_back(a);

  for (int i=0; i<3; ++i) {
    OBAtom *a1, *a2;
    a1 = atom_list[i];
    a2 = atom_list[i+1];
    b = tmol.NewBond();
    b->SetBegin(a1);
    b->SetEnd(a2);
    b->SetBondOrder(1);
    a1->AddBond(b);
    a2->AddBond(b);
  }
  tmol.GetBond(atom_list[0], atom_list[1])->SetBondOrder(2);

  tmol.EndModify();
  OB_COMPARE( tmol.NumAtoms(), 4);
  OB_COMPARE( tmol.NumBonds(), 3);
}


void PeriodicTester::TestLengths(double a, double b, double c) {
  std::vector<double> expected;
  expected.push_back(a);
  expected.push_back(b);
  expected.push_back(c);
  for (int i=0; i<3; ++i) {
    OBAtom* a1 = atom_list[i];
    OBAtom* a2 = atom_list[i+1];
    OB_ASSERT( near( a1->GetDistance(a2), expected[i] ) );
  }
}


void PeriodicTester::TestAngles(double a, double b) {
  std::vector<double> expected;
  expected.push_back(a);
  expected.push_back(b);
  for (int i=0; i<2; ++i) {
    OBAtom* a1 = atom_list[i];
    OBAtom* a2 = atom_list[i+1];
    OBAtom* a3 = atom_list[i+2];
    OB_ASSERT( near( a1->GetAngle(a2, a3), expected[i] ) );
  }
}


void PeriodicTester::TestTorsion(double a) {
  double torsion = tmol.GetTorsion(atom_list[0], atom_list[1], atom_list[2], atom_list[3]);
  OB_ASSERT( near( torsion, a ) );
}


void PeriodicTester::MakePeriodic(double a) {
  OBUnitCell *uc = new OBUnitCell;
  uc->SetData(a, a, a, 90, 90, 90);
  tmol.SetData(uc);
  tmol.SetPeriodicMol();

  // Wrap coordinates into the UC (<3,3,1>, <1,1,1>, <9,1,1,>, and <9,1,9>)
  tmol.BeginModify();
  FOR_ATOMS_OF_MOL(a, tmol) {
    a->SetVector(uc->WrapCartesianCoordinate(a->GetVector()));
  }
  tmol.EndModify();
}



void testNonperiodicNegative() {
  // Base case to check that the base of the test is working, sans periodicity
  PeriodicTester testobj;
  testobj.TestLengths(2.0*sqrt(2), 2.0, 2.0);  // diagonal hypotenuse, straight, straight
  testobj.TestAngles(135.0, 90.0);  // diagonal within xy plane, then down in the z direction
  testobj.TestTorsion(90.0);  // orthogonal planes
}


void testPeriodicFlag() {
  // Check that periodic code isn't activated unless specified
  PeriodicTester testobj;
  testobj.MakePeriodic();

  // With periodicity, the code should return the same values as above in
  // testNonPeriodicNegative(), even though the coordinates are now wrapped.
  testobj.TestLengths(2.0*sqrt(2), 2.0, 2.0);
  testobj.TestAngles(135.0, 90.0);
  testobj.TestTorsion(90.0);

  // If the flag is not activated, the wrapped coordinates behave as-is.
  testobj.GetMol()->UnsetFlag(OB_PERIODIC_MOL);
  testobj.TestLengths(2.0*sqrt(2), 8.0, 8.0);
  testobj.TestAngles(45.0, 90.0);
  testobj.TestTorsion(90.0);
}


void testPeriodicCIFWrite() {
  PeriodicTester testobj;
  testobj.MakePeriodic();
  OBMol *mol = testobj.GetMol();

  mol->SetTitle("Test for periodic CIFs");
  OBConversion conv;
  conv.SetOutFormat("cif");
  conv.AddOption("g");
  std::istringstream full_test_cif(conv.WriteString(mol));
  std::string line;
  std::string bond_section;
  bool found_bonds = false;
  while (std::getline(full_test_cif, line)) {
    if (found_bonds || line.find("bond") != std::string::npos) {
      bond_section.append(line);
      bond_section.append("\n");  // This also adds a newline to the end
      found_bonds = true;
    }
  }

  const std::string expected_cif_bonds =
"    _geom_bond_atom_site_label_1\n\
    _geom_bond_atom_site_label_2\n\
    _geom_bond_distance\n\
    _geom_bond_site_symmetry_2\n\
    _ccdc_geom_bond_type\n\
    O0     N1        2.82843      .   D\n\
    N1     C2        2.00000  1_455   S\n\
    C2     Br3       2.00000  1_554   S\n\
";

  OB_COMPARE(expected_cif_bonds, bond_section);
}


void testPeriodicNoncubic() {
  PeriodicTester testobj;
  OBMol *mol = testobj.GetMol();
  OBUnitCell *uc = new OBUnitCell;
  uc->SetData(15, 20, 11, 60.0, 78.8, 128.2);  // non-special triclinic parameters
  mol->SetData(uc);
  mol->SetPeriodicMol();

  mol->BeginModify();
  FOR_ATOMS_OF_MOL(a, *mol) {
    a->SetVector(uc->WrapCartesianCoordinate(a->GetVector()));
  }
  mol->EndModify();

  // When properly wrapped, the original coordinates should be invariant to
  // the selected unit cell, as long as it's large enough, etc.
  testobj.TestLengths(2.0*sqrt(2), 2.0, 2.0);
  testobj.TestAngles(135.0, 90.0);
  testobj.TestTorsion(90.0);
}



int periodictest(int argc, char* argv[])
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
    testNonperiodicNegative();
    break;
  case 2:
    testPeriodicFlag();
    break;
  case 3:
    testPeriodicCIFWrite();
    break;
  case 4:
    testPeriodicNoncubic();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return(0);
}
