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

#include <math.h>

#include "mol.h"
#include "math/matrix3x3.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

using namespace std;

namespace OpenBabel {

  /*! \class matrix3x3
     \brief Represents a real 3x3 matrix.

Rotating points in space can be performed by a vector-matrix
multiplication. The matrix3x3 class is designed as a helper to the
vector3 class for rotating points in space. The rotation matrix may be
initialised by passing in the array of doubleing point values, by
passing euler angles, or a rotation vector and angle of rotation about
that vector. Once set, the matrix3x3 class can be used to rotate
vectors by the overloaded multiplication operator. The following
demonstrates the usage of the matrix3x3 class:
/code
  matrix3x3 mat;
  mat.SetupRotMat(0.0,180.0,0.0); //rotate theta by 180 degrees
  vector3 v = VX;
  v *= mat; //apply the rotation
/endcode
*/

  /*! the axis of the rotation will be uniformly distributed on
    the unit sphere, the angle will be uniformly distributed in
    the interval 0..360 degrees. */
void matrix3x3::randomRotation(OBRandom &rnd)
{ 
  double rotAngle;
  vector3 v1;
 
  v1.randomUnitVector(&rnd);
  rotAngle = 360.0 * rnd.NextFloat();
  this->RotAboutAxisByAngle(v1,rotAngle);
}

void matrix3x3::SetupRotMat(double phi,double theta,double psi)
{
  double p  = phi * DEG_TO_RAD;
  double h  = theta * DEG_TO_RAD;
  double b  = psi * DEG_TO_RAD;

  double cx = cos(p);  double sx = sin(p);
  double cy = cos(h);  double sy = sin(h);
  double cz = cos(b);  double sz = sin(b);

  ele[0][0] = cy*cz;
  ele[0][1] = cy*sz;
  ele[0][2] = -sy;
  
  ele[1][0] = sx*sy*cz-cx*sz;
  ele[1][1] = sx*sy*sz+cx*cz;
  ele[1][2] = sx*cy;

  ele[2][0] = cx*sy*cz+sx*sz;
  ele[2][1] = cx*sy*sz-sx*cz;
  ele[2][2] = cx*cy;
}

/*! Replaces *this with a matrix that represents reflection on
  the plane through 0 which is given by the normal vector norm.
	
  \warning If the vector norm has length zero, this method will
  generate the 0-matrix. If the length of the axis is close to
  zero, but not == 0.0, this method may behave in unexpected
  ways and return almost random results; details may depend on
  your particular doubleing point implementation. The use of this
  method is therefore highly discouraged, unless you are certain
  that the length is in a reasonable range, away from 0.0
  (Stefan Kebekus)
	
  \deprecated This method will probably replaced by a safer
  algorithm in the future.
	
  \todo Replace this method with a more fool-proof version.
	
  @param norm specifies the normal to the plane
*/
void matrix3x3::PlaneReflection(const vector3 &norm)
{
  //@@@ add a safety net

  vector3 normtmp = norm;
  normtmp.normalize();

  SetColumn(0, vector3(1,0,0) - 2*normtmp.x()*normtmp);
  SetColumn(1, vector3(0,1,0) - 2*normtmp.y()*normtmp);
  SetColumn(2, vector3(0,0,1) - 2*normtmp.z()*normtmp);
}

#define x vtmp.x()
#define y vtmp.y()
#define z vtmp.z()

/*! Replaces *this with a matrix that represents rotation about the
  axis by a an angle. 
	
  \warning If the vector axis has length zero, this method will
  generate the 0-matrix. If the length of the axis is close to
  zero, but not == 0.0, this method may behave in unexpected ways
  and return almost random results; details may depend on your
  particular doubleing point implementation. The use of this method
  is therefore highly discouraged, unless you are certain that the
  length is in a reasonable range, away from 0.0 (Stefan
  Kebekus)
	
  \deprecated This method will probably replaced by a safer
  algorithm in the future.
	
  \todo Replace this method with a more fool-proof version.
	
  @param axis specifies the axis of the rotation
  @param angle angle in degrees (0..360)
*/
void matrix3x3::RotAboutAxisByAngle(const vector3 &v,const double angle)
{
  double theta = angle*DEG_TO_RAD;
  double s = sin(theta);
  double c = cos(theta);
  double t = 1 - c;
  
  vector3 vtmp = v;
  vtmp.normalize();
  
  ele[0][0] = t*x*x + c;
  ele[0][1] = t*x*y + s*z;
  ele[0][2] = t*x*z - s*y;
  
  ele[1][0] = t*x*y - s*z;
  ele[1][1] = t*y*y + c;
  ele[1][2] = t*y*z + s*x;

  ele[2][0] = t*x*z + s*y;
  ele[2][1] = t*y*z - s*x;
  ele[2][2] = t*z*z + c;
}

#undef x
#undef y
#undef z

void matrix3x3::SetColumn(int col, const vector3 &v) throw(OBError)
{
  if (col > 2) {
    OBError er("matrix3x3::SetColumn(int col, const vector3 &v)",
               "The method was called with col > 2.",
               "This is a programming error in your application.");
    throw er;
  }

  ele[0][col] = v.x();
  ele[1][col] = v.y();
  ele[2][col] = v.z();
}

void matrix3x3::SetRow(int row, const vector3 &v) throw(OBError)
{
  if (row > 2) {
    OBError er("matrix3x3::SetRow(int row, const vector3 &v)",
               "The method was called with row > 2.",
               "This is a programming error in your application.");
    throw er;
  }

  ele[row][0] = v.x();
  ele[row][1] = v.y();
  ele[row][2] = v.z();
}

vector3 matrix3x3::GetColumn(unsigned int col) const throw(OBError)
{
  if (col > 2) {
    OBError er("matrix3x3::GetColumn(unsigned int col) const",
               "The method was called with col > 2.",
               "This is a programming error in your application.");
    throw er;
  }

  return vector3(ele[0][col], ele[1][col], ele[2][col]);
}

vector3 matrix3x3::GetRow(unsigned int row) const throw(OBError)
{
  if (row > 2) {
    OBError er("matrix3x3::GetRow(unsigned int row) const",
               "The method was called with row > 2.",
               "This is a programming error in your application.");
    throw er;
  }

  return vector3(ele[row][0], ele[row][1], ele[row][2]);
}

/*! calculates the product m*v of the matrix m and the column
  vector represented by v
*/
vector3 operator *(const matrix3x3 &m,const vector3 &v)
{
  vector3 vv;

  vv._vx = v._vx*m.ele[0][0] + v._vy*m.ele[0][1] + v._vz*m.ele[0][2];
  vv._vy = v._vx*m.ele[1][0] + v._vy*m.ele[1][1] + v._vz*m.ele[1][2];
  vv._vz = v._vx*m.ele[2][0] + v._vy*m.ele[2][1] + v._vz*m.ele[2][2];

  return(vv);
}

matrix3x3 operator *(const matrix3x3 &A,const matrix3x3 &B)
{
  matrix3x3 result;

  result.ele[0][0] = A.ele[0][0]*B.ele[0][0] + A.ele[0][1]*B.ele[1][0] + A.ele[0][2]*B.ele[2][0];
  result.ele[0][1] = A.ele[0][0]*B.ele[0][1] + A.ele[0][1]*B.ele[1][1] + A.ele[0][2]*B.ele[2][1];
  result.ele[0][2] = A.ele[0][0]*B.ele[0][2] + A.ele[0][1]*B.ele[1][2] + A.ele[0][2]*B.ele[2][2];

  result.ele[1][0] = A.ele[1][0]*B.ele[0][0] + A.ele[1][1]*B.ele[1][0] + A.ele[1][2]*B.ele[2][0];
  result.ele[1][1] = A.ele[1][0]*B.ele[0][1] + A.ele[1][1]*B.ele[1][1] + A.ele[1][2]*B.ele[2][1];
  result.ele[1][2] = A.ele[1][0]*B.ele[0][2] + A.ele[1][1]*B.ele[1][2] + A.ele[1][2]*B.ele[2][2];

  result.ele[2][0] = A.ele[2][0]*B.ele[0][0] + A.ele[2][1]*B.ele[1][0] + A.ele[2][2]*B.ele[2][0];
  result.ele[2][1] = A.ele[2][0]*B.ele[0][1] + A.ele[2][1]*B.ele[1][1] + A.ele[2][2]*B.ele[2][1];
  result.ele[2][2] = A.ele[2][0]*B.ele[0][2] + A.ele[2][1]*B.ele[1][2] + A.ele[2][2]*B.ele[2][2];

  return(result);
}

/*! calculates the product m*(*this) of the matrix m and the
  column vector represented by *this
*/
vector3 &vector3::operator *= (const matrix3x3 &m)
{
  vector3 vv;
  
  vv.SetX(_vx*m.Get(0,0) + _vy*m.Get(0,1) + _vz*m.Get(0,2));
  vv.SetY(_vx*m.Get(1,0) + _vy*m.Get(1,1) + _vz*m.Get(1,2));
  vv.SetZ(_vx*m.Get(2,0) + _vy*m.Get(2,1) + _vz*m.Get(2,2));
  _vx = vv.x();
  _vy = vv.y();
  _vz = vv.z();
  
  return(*this);
}

/*! This method checks if the absolute value of the determinant is smaller than 1e-6. If
  so, nothing is done and an exception is thrown. Otherwise, the
  inverse matrix is calculated and returned. *this is not changed.
	
  \warning If the determinant is close to zero, but not == 0.0,
  this method may behave in unexpected ways and return almost
  random results; details may depend on your particular doubleing
  point implementation. The use of this method is therefore highly
  discouraged, unless you are certain that the determinant is in a
  reasonable range, away from 0.0 (Stefan Kebekus)
*/
matrix3x3 matrix3x3::inverse(void) const throw(OBError)
{
  double det = determinant();
  if (fabs(det) <= 1e-6) {
    OBError er("matrix3x3::invert(void)",
               "The method was called on a matrix with |determinant| <= 1e-6.",
               "This is a runtime or a programming error in your application.");
    throw er;
  }

  matrix3x3 inverse;
  inverse.ele[0][0] = ele[1][1]*ele[2][2] - ele[1][2]*ele[2][1];
  inverse.ele[1][0] = ele[1][2]*ele[2][0] - ele[1][0]*ele[2][2];
  inverse.ele[2][0] = ele[1][0]*ele[2][1] - ele[1][1]*ele[2][0];
  inverse.ele[0][1] = ele[2][1]*ele[0][2] - ele[2][2]*ele[0][1];
  inverse.ele[1][1] = ele[2][2]*ele[0][0] - ele[2][0]*ele[0][2];
  inverse.ele[2][1] = ele[2][0]*ele[0][1] - ele[2][1]*ele[0][0];
  inverse.ele[0][2] = ele[0][1]*ele[1][2] - ele[0][2]*ele[1][1];
  inverse.ele[1][2] = ele[0][2]*ele[1][0] - ele[0][0]*ele[1][2];
  inverse.ele[2][2] = ele[0][0]*ele[1][1] - ele[0][1]*ele[1][0];
  
  inverse /= det;
  
  return(inverse);
}

/* This method returns the transpose of a matrix. The original
   matrix remains unchanged. */
matrix3x3 matrix3x3::transpose(void) const 
{
  matrix3x3 transpose;

  for(unsigned int i=0; i<3; i++)
    for(unsigned int j=0; j<3; j++)
      transpose.ele[i][j] = ele[j][i];
  
  return(transpose);
}

double matrix3x3::determinant(void) const
{
  double x,y,z;

  x = ele[0][0] * (ele[1][1] * ele[2][2] - ele[1][2] * ele[2][1]);
  y = ele[0][1] * (ele[1][2] * ele[2][0] - ele[1][0] * ele[2][2]);
  z = ele[0][2] * (ele[1][0] * ele[2][1] - ele[1][1] * ele[2][0]);

  return(x + y + z);
}

/*! This method returns false if there are indices i,j such that
  fabs(*this[i][j]-*this[j][i]) > 1e-6. Otherwise, it returns
  true. */
bool matrix3x3::isSymmetric(void) const
{
  if (fabs(ele[0][1] - ele[1][0]) > 1e-6)
    return false;
  if (fabs(ele[0][2] - ele[2][0]) > 1e-6)
    return false;
  if (fabs(ele[1][2] - ele[2][1]) > 1e-6)
    return false;
  return true;
}

/*! This method returns false if there are indices i != j such
  that fabs(*this[i][j]) > 1e-6. Otherwise, it returns true. */
bool matrix3x3::isDiagonal(void) const
{
  if (fabs(ele[0][1]) > 1e-6)
    return false;
  if (fabs(ele[0][2]) > 1e-6)
    return false;
  if (fabs(ele[1][2]) > 1e-6)
    return false;

  if (fabs(ele[1][0]) > 1e-6)
    return false;
  if (fabs(ele[2][0]) > 1e-6)
    return false;
  if (fabs(ele[2][1]) > 1e-6)
    return false;

  return true;
}

/*! This method returns false if there are indices i != j such
  that fabs(*this[i][j]) > 1e-6, or if there is an index i such
  that fabs(*this[i][j]-1) > 1e-6. Otherwise, it returns
  true. */
bool matrix3x3::isUnitMatrix(void) const
{
  if (!isDiagonal())
    return false;

  if (fabs(ele[0][0]-1) > 1e-6)
    return false;
  if (fabs(ele[1][1]-1) > 1e-6)
    return false;
  if (fabs(ele[2][2]-1) > 1e-6)
    return false;

  return true;
}

/*! This method employs the static method matrix3x3::jacobi(...)
  to find the eigenvalues and eigenvectors of a symmetric
  matrix. On entry it is checked if the matrix really is
  symmetric: if isSymmetric() returns 'false', an OBError is
  thrown.

  \note The jacobi algorithm is should work great for all
  symmetric 3x3 matrices. If you need to find the eigenvectors
  of a non-symmetric matrix, you might want to resort to the
  sophisticated routines of LAPACK.

  @param eigenvals a reference to a vector3 where the
  eigenvalues will be stored. The eigenvalues are ordered so
  that eigenvals[0] <= eigenvals[1] <= eigenvals[2].

  @return an orthogonal matrix whose ith column is an
  eigenvector for the eigenvalue eigenvals[i]. Here 'orthogonal'
  means that all eigenvectors have length one and are mutually
  orthogonal. The ith eigenvector can thus be conveniently
  accessed by the GetColumn() method, as in the following
  example.
  \code
  // Calculate eigenvectors and -values
  vector3 eigenvals;
  matrix3x3 eigenmatrix = somematrix.findEigenvectorsIfSymmetric(eigenvals);
  
  // Print the 2nd eigenvector
  cout << eigenmatrix.GetColumn(1) << endl;
  \endcode
  With these conventions, a matrix is diagonalized in the following way:
  \code
  // Diagonalize the matrix
  matrix3x3 diagonalMatrix = eigenmatrix.inverse() * somematrix * eigenmatrix;
  \endcode
  
*/
matrix3x3 matrix3x3::findEigenvectorsIfSymmetric(vector3 &eigenvals) const throw(OBError)
{
  matrix3x3 result;
  
  if (!isSymmetric()) {
    OBError er("matrix3x3::findEigenvectorsIfSymmetric(vector3 &eigenvals) const throw(OBError)",
               "The method was called on a matrix that was not symmetric, i.e. where isSymetric() == false.",
               "This is a runtime or a programming error in your application.");
    throw er;
  }

  double d[3];
  matrix3x3 copyOfThis = *this;

  jacobi(3, copyOfThis.ele[0], d, result.ele[0]);
  eigenvals.Set(d);

  return result;
}

matrix3x3 &matrix3x3::operator/=(const double &c)
{
  for (int row = 0;row < 3; row++)
    for (int col = 0;col < 3; col++)
      ele[row][col] /= c;

  return(*this);
}

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

void matrix3x3::FillOrth(double Alpha,double Beta, double Gamma, 
			 double A, double B, double C)
{
  double V;
  
  Alpha *= DEG_TO_RAD; Beta  *= DEG_TO_RAD; Gamma *= DEG_TO_RAD;
  
  V= 1.0 - SQUARE(cos(Alpha)) - SQUARE(cos(Beta)) - SQUARE(cos(Gamma)) 
    + 2.0 * cos(Alpha) * cos(Beta) *  cos(Gamma);
  V = sqrt(fabs(V))/sin(Gamma);
  
  ele[0][0] = A;
  ele[0][1] = B*cos(Gamma);
  ele[0][2] = C*cos(Beta);

  ele[1][0] = 0.0;
  ele[1][1] = B*sin(Gamma);
  ele[1][2] = C*(cos(Alpha)-cos(Beta)*cos(Gamma))/sin(Gamma);

  ele[2][0] = 0.0; 
  ele[2][1] = 0.0; 
  ele[2][2] = C*V;
}

ostream& operator<< ( ostream& co, const matrix3x3& m )

{
  co << "[ "
     << m.ele[0][0]
     << ", "
     << m.ele[0][1]
     << ", "
     << m.ele[0][2]
     << " ]" << endl;

  co << "[ "
     << m.ele[1][0]
     << ", "
     << m.ele[1][1]
     << ", "
     << m.ele[1][2]
     << " ]" << endl;

  co << "[ "
     << m.ele[2][0]
     << ", "
     << m.ele[2][1]
     << ", "
     << m.ele[2][2]
     << " ]" << endl;

  return co ;
}

/*! This static function computes the eigenvalues and
  eigenvectors of a SYMMETRIC nxn matrix. This method is used
  internally by OpenBabel, but may be useful as a general
  eigenvalue finder.
  
  The algorithm uses Jacobi transformations. It is described
  e.g. in Wilkinson, Reinsch "Handbook for automatic computation,
  Volume II: Linear Algebra", part II, contribution II/1. The
  implementation is also similar to the implementation in this
  book. This method is adequate to solve the eigenproblem for
  small matrices, of size perhaps up to 10x10. For bigger
  problems, you might want to resort to the sophisticated routines
  of LAPACK.
  
  \note If you plan to find the eigenvalues of a symmetric 3x3
  matrix, you will probably prefer to use the more convenient
  method findEigenvectorsIfSymmetric()
  
  @param n the size of the matrix that should be diagonalized
  
  @param a array of size n^2 which holds the symmetric matrix
  whose eigenvectors are to be computed. The convention is that
  the entry in row r and column c is addressed as a[n*r+c] where,
  of course, 0 <= r < n and 0 <= c < n. There is no check that the
  matrix is actually symmetric. If it is not, the behaviour of
  this function is undefined.  On return, the matrix is
  overwritten with junk.
  
  @param d pointer to a field of at least n doubles which will be
  overwritten. On return of this function, the entries d[0]..d[n-1]
  will contain the eigenvalues of the matrix.
  
  @param v an array of size n^2 where the eigenvectors will be
  stored. On return, the columns of this matrix will contain the
  eigenvectors. The eigenvectors are normalized and mutually
  orthogonal.
*/
#define MAX_SWEEPS 50
void matrix3x3::jacobi(unsigned int n, double *a, double *d, double *v)
{
  double onorm, dnorm;
  double b, dma, q, t, c, s;
  double  atemp, vtemp, dtemp;
  register int i, j, k, l;
  int nrot;
  

  // Set v to the identity matrix, set the vector d to contain the
  // diagonal elements of the matrix a
  for (j = 0; j < static_cast<int>(n); j++) {
    for (i = 0; i < static_cast<int>(n); i++) 
      v[n*i+j] = 0.0;
    v[n*j+j] = 1.0;
    d[j] = a[n*j+j];
  }
  
  nrot = MAX_SWEEPS;
  for (l = 1; l <= nrot; l++) {
    // Set dnorm to be the maximum norm of the diagonal elements, set
    // onorm to the maximum norm of the off-diagonal elements
    dnorm = 0.0;
    onorm = 0.0;
    for (j = 0; j < static_cast<int>(n); j++) {
      dnorm += (double)fabs(d[j]);
      for (i = 0; i < j; i++)
	onorm += (double)fabs(a[n*i+j]);
    }
    // Normal end point of this algorithm.
    if((onorm/dnorm) <= 1.0e-12)
      goto Exit_now;
    
    for (j = 1; j < static_cast<int>(n); j++) {
      for (i = 0; i <= j - 1; i++) {

	b = a[n*i+j];
	if(fabs(b) > 0.0) {
	  dma = d[j] - d[i];
	  if((fabs(dma) + fabs(b)) <=  fabs(dma))
	    t = b / dma;
	  else {
	    q = 0.5 * dma / b;
	    t = 1.0/((double)fabs(q) + (double)sqrt(1.0+q*q));
	    if (q < 0.0)
	      t = -t;
	  }
	  
	  c = 1.0/(double)sqrt(t*t + 1.0);
	  s = t * c;
	  a[n*i+j] = 0.0;
	  
	  for (k = 0; k <= i-1; k++) {
	    atemp = c * a[n*k+i] - s * a[n*k+j];
	    a[n*k+j] = s * a[n*k+i] + c * a[n*k+j];
	    a[n*k+i] = atemp;
	  }
	  
	  for (k = i+1; k <= j-1; k++) {
	    atemp = c * a[n*i+k] - s * a[n*k+j];
	    a[n*k+j] = s * a[n*i+k] + c * a[n*k+j];
	    a[n*i+k] = atemp;
	  }
	  
	  for (k = j+1; k < static_cast<int>(n); k++) {
	    atemp = c * a[n*i+k] - s * a[n*j+k];
	    a[n*j+k] = s * a[n*i+k] + c * a[n*j+k];
	    a[n*i+k] = atemp;
	  }
	  
	  for (k = 0; k < static_cast<int>(n); k++) {
	    vtemp = c * v[n*k+i] - s * v[n*k+j];
	    v[n*k+j] = s * v[n*k+i] + c * v[n*k+j];
	    v[n*k+i] = vtemp;
	  }
	  
	  dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
	  d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
	  d[i] = dtemp;
	} /* end if */
      } /* end for i */
    } /* end for j */
  } /* end for l */
  
 Exit_now:
  
  // Now sort the eigenvalues (and the eigenvectors) so that the
  // smallest eigenvalues come first.
  nrot = l;
  
  for (j = 0; j < static_cast<int>(n)-1; j++) {
    k = j;
    dtemp = d[k];
    for (i = j+1; i < static_cast<int>(n); i++)
      if(d[i] < dtemp) {
	k = i;
	dtemp = d[k];
      }

    if(k > j) {
      d[k] = d[j];
      d[j] = dtemp;
      for (i = 0; i < static_cast<int>(n); i++) {
	dtemp = v[n*i+k];
	v[n*i+k] = v[n*i+j];
	v[n*i+j] = dtemp;
      }
    }
  }
}


}

