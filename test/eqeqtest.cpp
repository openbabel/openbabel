/**********************************************************************
eqeqtest.cpp - Unit tests for the EQEq periodic charge model.

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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include <openbabel/generic.h>
#include <openbabel/chargemodel.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace std;
using namespace OpenBabel;

// The model assembles and solves its linear system in single precision,
// so charges are only reproducible to ~1e-6 e. The defect being guarded
// against moves charges by 1e-2 e or more.
static const double CHARGE_TOL = 1.0e-4;

// Rock-salt Na4Cl4 in fractional coordinates.
static const int    kNumAtoms = 8;
static const int    kElem[kNumAtoms] = { 11, 11, 11, 11, 17, 17, 17, 17 };
static const double kFrac[kNumAtoms][3] = {
  {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
  {0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}, {0.5, 0.5, 0.5}
};

// Build the crystal on lattice vectors a, b, c (Cartesian, Angstrom),
// with every atom and every lattice vector rotated by `rot`.
static void BuildCrystal(OBMol &mol, const vector3 &a, const vector3 &b,
                         const vector3 &c, const matrix3x3 &rot)
{
  mol.Clear();
  mol.BeginModify();
  for (int i = 0; i < kNumAtoms; ++i) {
    OBAtom *atom = mol.NewAtom();
    atom->SetAtomicNum(kElem[i]);
    vector3 r = kFrac[i][0] * a + kFrac[i][1] * b + kFrac[i][2] * c;
    atom->SetVector(rot * r);
  }
  mol.EndModify();

  OBUnitCell *cell = new OBUnitCell;
  cell->SetData(rot * a, rot * b, rot * c);
  mol.SetData(cell);
}

static vector<double> ComputeEQEq(OBMol &mol)
{
  OBChargeModel *model = OBChargeModel::FindType("eqeq");
  OB_REQUIRE(model != nullptr);
  OB_REQUIRE(model->ComputeCharges(mol));
  vector<double> q;
  FOR_ATOMS_OF_MOL(atom, mol)
    q.push_back(atom->GetPartialCharge());
  return q;
}

static double MaxAbsDiff(const vector<double> &x, const vector<double> &y)
{
  OB_REQUIRE(x.size() == y.size());
  double m = 0.0;
  for (size_t i = 0; i < x.size(); ++i)
    m = max(m, fabs(x[i] - y[i]));
  return m;
}

// Rigidly rotating a crystal must not change its charges.
static void CheckRotationInvariance(const vector3 &a, const vector3 &b,
                                    const vector3 &c, const matrix3x3 &rot,
                                    const char *label)
{
  matrix3x3 identity(1.0);
  OBMol ref, rotated;
  BuildCrystal(ref, a, b, c, identity);
  BuildCrystal(rotated, a, b, c, rot);

  vector<double> q0 = ComputeEQEq(ref);
  vector<double> q1 = ComputeEQEq(rotated);
  double diff = MaxAbsDiff(q0, q1);
  printf("%s: max |q(rotated) - q(original)| = %.3e e\n", label, diff);
  OB_ASSERT(diff < CHARGE_TOL);
}

// Orthorhombic cell given in a rotated frame (all angles still 90 degrees),
// as a POSCAR or other vector-based format would supply it.
void testRotatedOrthorhombic()
{
  matrix3x3 rot;
  rot.RotAboutAxisByAngle(vector3(0.0, 0.0, 1.0), 37.0);
  CheckRotationInvariance(vector3(6.0, 0.0, 0.0), vector3(0.0, 7.0, 0.0),
                          vector3(0.0, 0.0, 8.0), rot, "rotated orthorhombic");
}

// Triclinic cell, rotated about a general axis.
void testRotatedTriclinic()
{
  matrix3x3 rot;
  rot.RotAboutAxisByAngle(vector3(1.0, 2.0, 3.0), 50.0);
  CheckRotationInvariance(vector3(6.0, 0.0, 0.0), vector3(1.2, 7.0, 0.0),
                          vector3(0.9, 1.7, 8.0), rot, "rotated triclinic");
}

// The same triclinic crystal entered once as lattice vectors and once as
// cell parameters (a, b, c, alpha, beta, gamma), the path CIF input takes.
void testTriclinicVectorsVsParameters()
{
  vector3 a(6.0, 0.0, 0.0), b(1.2, 7.0, 0.0), c(0.9, 1.7, 8.0);
  matrix3x3 rot;
  rot.RotAboutAxisByAngle(vector3(-2.0, 1.0, 0.5), 71.0);

  OBMol fromVectors;
  BuildCrystal(fromVectors, a, b, c, rot);

  OBMol fromParams;
  BuildCrystal(fromParams, a, b, c, matrix3x3(1.0));
  OBUnitCell *cell = static_cast<OBUnitCell *>(
      fromParams.GetData(OBGenericDataType::UnitCell));
  cell->SetData(a.length(), b.length(), c.length(),
                vectorAngle(b, c), vectorAngle(a, c), vectorAngle(a, b));
  // Re-place the atoms in the frame the parameter form defines.
  for (int i = 0; i < kNumAtoms; ++i)
    fromParams.GetAtom(i + 1)->SetVector(
        cell->FractionalToCartesian(vector3(kFrac[i][0], kFrac[i][1], kFrac[i][2])));

  double diff = MaxAbsDiff(ComputeEQEq(fromVectors), ComputeEQEq(fromParams));
  printf("vectors vs parameters: max |dq| = %.3e e\n", diff);
  OB_ASSERT(diff < CHARGE_TOL);
}

// Charges must sum to zero (neutrality is imposed as a constraint row).
void testNeutrality()
{
  matrix3x3 rot;
  rot.RotAboutAxisByAngle(vector3(1.0, 2.0, 3.0), 50.0);
  OBMol mol;
  BuildCrystal(mol, vector3(6.0, 0.0, 0.0), vector3(1.2, 7.0, 0.0),
               vector3(0.9, 1.7, 8.0), rot);
  vector<double> q = ComputeEQEq(mol);
  double total = 0.0;
  for (size_t i = 0; i < q.size(); ++i)
    total += q[i];
  printf("total charge = %.3e e\n", total);
  OB_ASSERT(fabs(total) < CHARGE_TOL);
}

// ---- Non-periodic path (issue #2771) ----

// Water, no unit cell.
static void BuildWater(OBMol &mol)
{
  mol.Clear();
  mol.BeginModify();
  const int    z[3]    = { 8, 1, 1 };
  const double r[3][3] = { {0.0, 0.0, 0.0}, {0.757, 0.586, 0.0}, {-0.757, 0.586, 0.0} };
  for (int i = 0; i < 3; ++i) {
    OBAtom *atom = mol.NewAtom();
    atom->SetAtomicNum(z[i]);
    atom->SetVector(r[i][0], r[i][1], r[i][2]);
  }
  mol.EndModify();
}

// An isolated molecule must get a successful, neutral, chemically sensible result.
void testNonperiodicWater()
{
  OBMol mol;
  BuildWater(mol);
  vector<double> q = ComputeEQEq(mol);
  printf("water: q(O) = %.4f, q(H) = %.4f, %.4f\n", q[0], q[1], q[2]);
  OB_ASSERT(fabs(q[0] + q[1] + q[2]) < CHARGE_TOL);
  OB_ASSERT(q[0] < 0.0);
  OB_ASSERT(q[1] > 0.0);
  OB_ASSERT(fabs(q[1] - q[2]) < CHARGE_TOL);  // the two H are equivalent
}

// An isolated molecule and the same molecule alone in a very large periodic
// box should get nearly the same charges.
void testNonperiodicMatchesLargeBox()
{
  OBMol isolated, boxed;
  BuildWater(isolated);
  BuildWater(boxed);
  OBUnitCell *cell = new OBUnitCell;
  cell->SetData(vector3(100.0, 0.0, 0.0), vector3(0.0, 100.0, 0.0),
                vector3(0.0, 0.0, 100.0));
  boxed.SetData(cell);

  double diff = MaxAbsDiff(ComputeEQEq(isolated), ComputeEQEq(boxed));
  printf("isolated vs 100 A box: max |dq| = %.3e e\n", diff);
  OB_ASSERT(diff < 1.0e-3);
}

int eqeqtest(int argc, char *argv[])
{
  int choice = 1;
  if (argc > 1) {
    if (sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  // Charge models are loadable plugins in shared builds.
  #ifdef FORMATDIR
    char env[BUFF_SIZE];
    snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
    putenv(env);
  #endif

  switch (choice) {
  case 1:
    testRotatedOrthorhombic();
    break;
  case 2:
    testRotatedTriclinic();
    break;
  case 3:
    testTriclinicVectorsVsParameters();
    break;
  case 4:
    testNeutrality();
    break;
  case 5:
    testNonperiodicWater();
    break;
  case 6:
    testNonperiodicMatchesLargeBox();
    break;
  default:
    cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }
  return 0;
}
