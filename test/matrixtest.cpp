/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "math/matrix3x3.h"
#include "obutil.h"

using namespace std;
using namespace OpenBabel;

OBRandom randomizer;

matrix3x3 randomMatrix(void)
{
  matrix3x3 A;

  do{
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
	A.Set(i, j, randomizer.NextFloat());
  }while(A.determinant() < 0.1);

  return A;
}

// Test the inversion of matrices and the matrix product. Get an
// invertible random matrix, invert it, multiply and check if the
// resulting matrix is the unit matrix.
bool testInversion() 
{
  matrix3x3 rnd = randomMatrix();
  matrix3x3 result = rnd * rnd.inverse();

  if (!result.isUnitMatrix()) {
    cout << "Matrix inversion test failed." << endl;
    return false;
  }
  return true;
}


// Test the eigenvalue finder. Set up a diagonal matrix and conjugate
// by a rotation. That way we obtain a symmetric matrix that can be
// diagonalized. Check if the eigenvalue finder reveals the original
// diagonal elements.
bool testEigenvalues()
{
  matrix3x3 Diagonal;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      Diagonal.Set(i, j, 0.0);
  Diagonal.Set(0, 0, randomizer.NextFloat());
  Diagonal.Set(1, 1, Diagonal.Get(0,0)+fabs(randomizer.NextFloat()));
  Diagonal.Set(2, 2, Diagonal.Get(1,1)+fabs(randomizer.NextFloat()));

  matrix3x3 rndRotation;
  rndRotation.randomRotation(randomizer);

  matrix3x3 toDiagonalize = rndRotation * Diagonal * rndRotation.inverse();
  if (!toDiagonalize.isSymmetric()) {
    cout << "Matrix eigenvalue test failed, conjugation of a diagonal matrix"
	 << "with a rotation is not symmetric." << endl;
    return false;
  }
  
  vector3 eigenvals;
  toDiagonalize.findEigenvectorsIfSymmetric(eigenvals);

  for(unsigned int j=0; j<3; j++)
    if ( fabs(eigenvals[j]-Diagonal.Get(j,j)) > 2e-6 )
    {
      cout << "Matrix eigenvalue test(" << j << 
	") failed, wrong eigenvalues computed." << endl;
      cout << "Expected: " << eigenvals[j] << " and got instead: " 
	   << Diagonal.Get(j,j) << endl;
      return false;
    }

  if ( (eigenvals[0] >= eigenvals[1]) || (eigenvals[1] >= eigenvals[2]) ) {
    cout << "Matrix eigenvalue test failed, eigenvalues not not ordered." << endl;
    return false;
  }
  return true;
}


// Test the eigenvector finder. Set up a symmetric diagonal matrix and
// diagonalize it. The that conjugation with the computed eigenvectors
// really gives a diagonal matrix. Calculate the matrix-vector
// products directly to see if the vectors found are eigenvectors
// indeed, and if the computed eigenvalues are correct.
bool testEigenvectors() 
{
  matrix3x3 rnd = randomMatrix();
  rnd.Set(0,1, rnd.Get(1,0));
  rnd.Set(0,2, rnd.Get(2,0));
  rnd.Set(1,2, rnd.Get(2,1));
  vector3 eigenvals;
  matrix3x3 eigenvects = rnd.findEigenvectorsIfSymmetric(eigenvals);

  // Check if the eigenvects are normalized and mutually orthogonal
  if (!eigenvects.isOrthogonal()) {
    cout << "Matrix eigenvector test failed, findEigenvectorsIfSymmetric()"
	 << " returned a matrix that is not orthogonal." << endl;
    return false;
  }

  matrix3x3 shouldBeDiagonal = eigenvects.inverse() * rnd * eigenvects;
  
  if (!shouldBeDiagonal.isDiagonal()) {
    cout << "Matrix eigenvector test failed, matrix is not diagonalized." << endl;
    return false;
  }
  
  for(unsigned int j=0; j<3; j++)
    if ( fabs(shouldBeDiagonal.Get(j,j) - eigenvals[j]) > 2e-6 )
    {
      cout << "Matrix eigenvector test(" << j << 
	") failed, wrong eigenvalues computed." << endl;
      cout << "Expected: " << eigenvals[j] << " and got instead: " 
	   << shouldBeDiagonal.Get(j,j) << endl;
      return false;
    }

  for(unsigned int i=0; i<3; i++ ) {
    vector3 EV = eigenvects.GetColumn(i);
    EV = rnd*EV-eigenvals[i]*EV;
    if (EV.length() > 1e-4) {
      cout << EV.length() << " Matrix eigenvector test failed,"
	   << " wrong eigenvector computed for column " << i << "." << endl;
      return false;
    }
  }
  return true;
}


bool TestMatrixAlgebra(void)
{
  randomizer.Seed(346534);

  for (int i=0; i < 50000 ; i++) {
    if (!testInversion())
      return false;
    if (!testEigenvalues())
      return false;
    if (!testEigenvectors())
      return false;
  }

  return true;
}
