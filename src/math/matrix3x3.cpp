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

void matrix3x3::SetColumn(int column, const vector3 v)
{
  ele[0][column] = v.x();
  ele[1][column] = v.y();
  ele[2][column] = v.z();
}

vector3 operator *(const matrix3x3 &m,const vector3 &v)
{
  vector3 vv;

  vv._vx = v._vx*m.ele[0][0] + v._vy*m.ele[0][1] + v._vz*m.ele[0][2];
  vv._vy = v._vx*m.ele[1][0] + v._vy*m.ele[1][1] + v._vz*m.ele[1][2];
  vv._vz = v._vx*m.ele[2][0] + v._vy*m.ele[2][1] + v._vz*m.ele[2][2];

  return(vv);
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

matrix3x3 matrix3x3::invert()
{
  float t[3][3];
  float det;
  
  det = determinant();
  
  if (det != 0.0) {
    t[0][0] = ele[1][1]*ele[2][2] - ele[1][2]*ele[2][1];
    t[1][0] = ele[1][2]*ele[2][0] - ele[1][0]*ele[2][2];
    t[2][0] = ele[1][0]*ele[2][1] - ele[1][1]*ele[2][0];
    t[0][1] = ele[2][1]*ele[0][2] - ele[2][2]*ele[0][1];
    t[1][1] = ele[2][2]*ele[0][0] - ele[2][0]*ele[0][2];
    t[2][1] = ele[2][0]*ele[0][1] - ele[2][1]*ele[0][0];
    t[0][2] = ele[0][1]*ele[1][2] - ele[0][2]*ele[1][1];
    t[1][2] = ele[0][2]*ele[1][0] - ele[0][0]*ele[1][2];
    t[2][2] = ele[0][0]*ele[1][1] - ele[0][1]*ele[1][0];
    
    register int i,j;
    for (i = 0;i < 3;i++)
      for (j = 0;j < 3;j++)
	ele[i][j] = t[i][j];
    
    *this /= det;
  }
  
  return(*this);
}

float matrix3x3::determinant()
{
  float x,y,z;

  x = ele[0][0] * (ele[1][1] * ele[2][2] - ele[1][2] * ele[2][1]);
  y = ele[0][1] * (ele[1][2] * ele[2][0] - ele[1][0] * ele[2][2]);
  z = ele[0][2] * (ele[1][0] * ele[2][1] - ele[1][1] * ele[2][0]);

  return(x + y + z);
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

}

