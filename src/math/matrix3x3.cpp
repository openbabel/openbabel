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


void matrix3x3::randomRotation(OBRandom &rnd)
{ 
  float rotAngle;
  vector3 v1;
 
  v1.randomUnitVector(&rnd);
  rotAngle = 360.0f * rnd.NextFloat();
  this->RotAboutAxisByAngle(v1,rotAngle);
}

void matrix3x3::SetupRotMat(float phi,float theta,float psi)
{
  float p  = phi * DEG_TO_RAD;
  float h  = theta * DEG_TO_RAD;
  float b  = psi * DEG_TO_RAD;

  float cx = cos(p);  float sx = sin(p);
  float cy = cos(h);  float sy = sin(h);
  float cz = cos(b);  float sz = sin(b);

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

void matrix3x3::RotAboutAxisByAngle(const vector3 &v,const float angle)
{
  float theta = angle*DEG_TO_RAD;
  float s = sin(theta);
  float c = cos(theta);
  float t = 1 - c;
  
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

matrix3x3 matrix3x3::inverse(void) const throw(OBError)
{
  float det = determinant();
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

matrix3x3 matrix3x3::transpose(void) const 
{
  matrix3x3 transpose;

  for(unsigned int i=0; i<3; i++)
    for(unsigned int j=0; j<3; j++)
      transpose.ele[i][j] = ele[j][i];
  
  return(transpose);
}

float matrix3x3::determinant(void) const
{
  float x,y,z;

  x = ele[0][0] * (ele[1][1] * ele[2][2] - ele[1][2] * ele[2][1]);
  y = ele[0][1] * (ele[1][2] * ele[2][0] - ele[1][0] * ele[2][2]);
  z = ele[0][2] * (ele[1][0] * ele[2][1] - ele[1][1] * ele[2][0]);

  return(x + y + z);
}

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

matrix3x3 matrix3x3::findEigenvectorsIfSymmetric(vector3 &eigenvals) const throw(OBError)
{
  matrix3x3 result;

  if (!isSymmetric()) {
    OBError er("matrix3x3::findEigenvectorsIfSymmetric(vector3 &eigenvals) const throw(OBError)",
               "The method was called on a matrix that was not symmetric, i.e. where isSymetric() == false.",
               "This is a runtime or a programming error in your application.");
    throw er;
  }

  float d[3];
  matrix3x3 copyOfThis = *this;

  jacobi(3, copyOfThis.ele[0], d, result.ele[0]);
  eigenvals.Set(d);

  return result;
}

matrix3x3 &matrix3x3::operator/=(const float &c)
{
  for (int row = 0;row < 3; row++)
    for (int col = 0;col < 3; col++)
      ele[row][col] /= c;

  return(*this);
}

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

void matrix3x3::FillOrth(float Alpha,float Beta, float Gamma, 
			 float A, float B, float C)
{
  float V;
  
  Alpha *= DEG_TO_RAD; Beta  *= DEG_TO_RAD; Gamma *= DEG_TO_RAD;
  
  V= 1.0f - SQUARE(cos(Alpha)) - SQUARE(cos(Beta)) - SQUARE(cos(Gamma)) 
    + 2.0f * cos(Alpha) * cos(Beta) *  cos(Gamma);
  V = sqrt(fabs(V))/sin(Gamma);
  
  ele[0][0] = A;
  ele[0][1] = B*cos(Gamma);
  ele[0][2] = C*cos(Beta);

  ele[1][0] = 0.0f;
  ele[1][1] = B*sin(Gamma);
  ele[1][2] = C*(cos(Alpha)-cos(Beta)*cos(Gamma))/sin(Gamma);

  ele[2][0] = 0.0f; 
  ele[2][1] = 0.0f; 
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


#define MAX_SWEEPS 50
void matrix3x3::jacobi(unsigned int n, float *a, float *d, float *v)
{
  float onorm, dnorm;
  float b, dma, q, t, c, s;
  float  atemp, vtemp, dtemp;
  register int i, j, k, l;
  int nrot;
  

  // Set v to the identity matrix, set the vector d to contain the
  // diagonal elements of the matrix a
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) 
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
    for (j = 0; j < n; j++) {
      dnorm += (float)fabs(d[j]);
      for (i = 0; i < j; i++)
	onorm += (float)fabs(a[n*i+j]);
    }
    // Normal end point of this algorithm.
    if((onorm/dnorm) <= 1.0e-12)
      goto Exit_now;
    
    for (j = 1; j < n; j++) {
      for (i = 0; i <= j - 1; i++) {

	b = a[n*i+j];
	if(fabs(b) > 0.0) {
	  dma = d[j] - d[i];
	  if((fabs(dma) + fabs(b)) <=  fabs(dma))
	    t = b / dma;
	  else {
	    q = 0.5f * dma / b;
	    t = 1.0f/((float)fabs(q) + (float)sqrt(1.0f+q*q));
	    if (q < 0.0f)
	      t = -t;
	  }
	  
	  c = 1.0f/(float)sqrt(t*t + 1.0f);
	  s = t * c;
	  a[n*i+j] = 0.0f;
	  
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
	  
	  for (k = j+1; k < n; k++) {
	    atemp = c * a[n*i+k] - s * a[n*j+k];
	    a[n*j+k] = s * a[n*i+k] + c * a[n*j+k];
	    a[n*i+k] = atemp;
	  }
	  
	  for (k = 0; k < n; k++) {
	    vtemp = c * v[n*k+i] - s * v[n*k+j];
	    v[n*k+j] = s * v[n*k+i] + c * v[n*k+j];
	    v[n*k+i] = vtemp;
	  }
	  
	  dtemp = c*c*d[i] + s*s*d[j] - 2.0f*c*s*b;
	  d[j] = s*s*d[i] + c*c*d[j] +  2.0f*c*s*b;
	  d[i] = dtemp;
	} /* end if */
      } /* end for i */
    } /* end for j */
  } /* end for l */
  
 Exit_now:
  
  // Now sort the eigenvalues (and the eigenvectors) so that the
  // smallest eigenvalues come first.
  nrot = l;
  
  for (j = 0; j < n-1; j++) {
    k = j;
    dtemp = d[k];
    for (i = j+1; i < n; i++)
      if(d[i] < dtemp) {
	k = i;
	dtemp = d[k];
      }

    if(k > j) {
      d[k] = d[j];
      d[j] = dtemp;
      for (i = 0; i < n; i++) {
	dtemp = v[n*i+k];
	v[n*i+k] = v[n*i+j];
	v[n*i+j] = dtemp;
      }
    }
  }
}


}

