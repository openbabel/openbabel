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
#include <time.h>
#include "mol.h"
#include "Vector.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif


namespace OpenEye {

// create a random unit vector
// if seed is nonNegative then use this as the seed,
// otherwise do not seed
void Vector::randomUnitVector(OERandom *oeRandP)
{
   bool doFree= false;
   if (oeRandP == NULL){
      doFree= true;
      oeRandP= new OERandom(true);
      oeRandP->TimeSeed();
   }

   // make sure to sample in the unit sphere
   float f1= 0.0f, f2= 0.0f, f3= 0.0f;
   bool b= false;
   while(!b){
       f1= oeRandP->NextFloat();
       f2= oeRandP->NextFloat();
       f3= oeRandP->NextFloat();
       if(b=f1*f1 + f2*f2 + f3*f3 <= 1.0){
				 if (oeRandP->NextInt() % 2 == 0) f1 *= -1.0f;
				 if (oeRandP->NextInt() % 2 == 0) f2 *= -1.0f;
				 if (oeRandP->NextInt() % 2 == 0) f3 *= -1.0f;
			 }
   }
   this->Set(f1,f2,f3);
   this->normalize();
   if (doFree) { delete oeRandP; oeRandP= NULL; }
}

Vector :: Vector ( const float x, const float y, const float z )
{
  _vx = x ;
  _vy = y ;
  _vz = z ;
}

Vector :: Vector ( const Vector& v ) 
{
  _vx = v._vx ;
  _vy = v._vy ;
  _vz = v._vz ;
}

Vector :: ~Vector () { }

Vector& Vector :: operator= ( const Vector& v ) 
{
  if ( this == &v ) return ( *this ) ;
  
  _vx = v._vx ;
  _vy = v._vy ;
  _vz = v._vz ;
  
  return ( *this ) ;
}

ostream& operator<< ( ostream& co, const Vector& v )
{
  co << "< " << v._vx << ", " << v._vy << ", " << v._vz << " >" ;
  return co ;
}

int operator== ( const Vector& v1, const Vector& v2 ) 
{
  if ( ( v1._vx == v2._vx ) &&
       ( v1._vy == v2._vy ) &&
       ( v1._vz == v2._vz ) )
     return ( true ) ;
  else
     return ( false ) ;
}

int operator!= ( const Vector& v1, const Vector& v2 ) 
{
  if ( ( v1._vx != v2._vx ) ||
       ( v1._vy != v2._vy ) ||
       ( v1._vz != v2._vz ) )
     return ( true ) ;
  else
     return ( false ) ;
}

Vector operator+ ( const Vector& v1, const Vector& v2 ) 
{
  Vector vv ;

  vv._vx = v1._vx + v2._vx ;
  vv._vy = v1._vy + v2._vy ;
  vv._vz = v1._vz + v2._vz ;

  return ( vv ) ;
}

Vector operator- ( const Vector& v1, const Vector& v2 ) 
{
  Vector vv ;

  vv._vx = v1._vx - v2._vx ;
  vv._vy = v1._vy - v2._vy ;
  vv._vz = v1._vz - v2._vz ;

  return ( vv ) ;
}

Vector operator- ( const Vector& v ) 
{
  Vector vv ;

  vv._vx = - v._vx ;
  vv._vy = - v._vy ;
  vv._vz = - v._vz ;

  return ( vv ) ;
}

Vector operator* ( const float& c, const Vector& v ) 
{
  Vector vv ;

  vv._vx = c * v._vx ;
  vv._vy = c * v._vy ;
  vv._vz = c * v._vz ;

  return ( vv ) ;
}

Vector operator* ( const Vector& v, const float& c ) 
{
  Vector vv ;

  vv._vx = c * v._vx ;
  vv._vy = c * v._vy ;
  vv._vz = c * v._vz ;

  return ( vv ) ;
}

Vector operator*(const Vector &v,const Vector &v1)
{
  Vector vv ;

  vv._vx = v1._vx * v._vx ;
  vv._vy = v1._vy * v._vy ;
  vv._vz = v1._vz * v._vz ;

  return ( vv ) ;
}

Vector operator/ ( const Vector& v, const float& c ) 
{
  Vector vv ;

  vv._vx = v._vx / c ;
  vv._vy = v._vy / c ;
  vv._vz = v._vz / c ;

  return ( vv ) ;
}

Vector& Vector :: operator+= ( const Vector& v ) 
{
  _vx += v._vx ;
  _vy += v._vy ;
  _vz += v._vz ;

  return *this ;
}

Vector& Vector :: operator-= ( const Vector& v ) 
{
  _vx -= v._vx ;
  _vy -= v._vy ;
  _vz -= v._vz ;

  return *this ;
}

Vector& Vector :: operator+= ( const float *f ) 
{
  _vx += f[0];
  _vy += f[1];
  _vz += f[2];

  return *this ;
}

Vector& Vector :: operator-= ( const float *f ) 
{
  _vx -= f[0];
  _vy -= f[1];
  _vz -= f[2];

  return *this ;
}

Vector& Vector :: operator*= ( const float& c ) 
{
  _vx *= c ;
  _vy *= c ;
  _vz *= c ;

  return(*this);
}

Vector& Vector :: operator/= ( const float& c ) 
{
  _vx /= c ;
  _vy /= c ;
  _vz /= c ;

  return(*this);
}

Vector& Vector :: normalize ()  
{
  float l =  length () ;

  if (l == 0) return(*this);

  _vx = _vx / l ;
  _vy = _vy / l ;
  _vz = _vz / l ;

  return(*this);
}

float Vector :: length ()  const
{
  float l;

  l =  sqrt ( _vx * _vx + _vy * _vy + _vz * _vz ) ;
  return ( l ) ;
}

float Vector :: length_2 ()  const
{
  float l;

  l =  ( _vx * _vx + _vy * _vy + _vz * _vz ) ;
  return ( l ) ;
}

float dot ( const Vector& v1, const Vector& v2 ) 
{
  float d;

  d =  v1._vx * v2._vx + v1._vy * v2._vy + v1._vz * v2._vz ;
  return ( d ) ;
}

Vector cross ( const Vector& v1, const Vector& v2 ) 
{
  Vector vv ;

  vv._vx = v1._vy * v2._vz - v1._vz * v2._vy ;
  vv._vy = - v1._vx * v2._vz + v1._vz * v2._vx ;
  vv._vz = v1._vx * v2._vy - v1._vy * v2._vx ;

  return ( vv ) ;
}


// ***angle***

float VectorAngle ( const Vector& v1, const Vector& v2 ) 
{
  float mag;
  float dp;

  mag = v1.length() * v2.length();
  dp = dot(v1,v2)/mag;

  if (dp < -0.999999)
     dp = -0.9999999f;

  if (dp > 0.9999999)
     dp = 0.9999999f;

  if (dp > 1.0)
     dp = 1.0f;

  return((RAD_TO_DEG * acos(dp)));
}

float CalcTorsionAngle(Vector &a,Vector &b,Vector &c,Vector &d)
{
  float torsion;
  Vector b1,b2,b3,c1,c2,c3;
  
  b1 = a - b;
  b2 = b - c;
  b3 = c - d;
  
  c1 = cross(b1,b2);
  c2 = cross(b2,b3);
  c3 = cross(c1,c2);
  
  if (c1.length() * c2.length() < 0.001)
     torsion = 0.0;
  else
  {
    torsion = VectorAngle(c1,c2);
    if (dot(b2,c3) > 0.0)
       torsion *= -1.0;
  }

  return(torsion);
}

//MATRIX ROUTINES


void Matrix3x3::randomRotation(OERandom &rnd)
{ 
    Vector v1; 
    v1.randomUnitVector(&rnd);
    float rotAngle= float(rnd.NextInt() % 36000)/100.0;
    if (rnd.NextInt() % 2 == 0) rotAngle *= -1.0f;
    this->RotAboutAxisByAngle(v1,rotAngle);
}

void Matrix3x3::SetupRotMat(float phi,float theta,float psi)
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

#define x vtmp.x()
#define y vtmp.y()
#define z vtmp.z()

void Matrix3x3::RotAboutAxisByAngle(const Vector &v,const float angle)
{
  float theta = angle*DEG_TO_RAD;
  float s = sin(theta);
  float c = cos(theta);
  float t = 1 - c;
  
  Vector vtmp = v;
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

void Matrix3x3::RotateCoords(float *c,int noatoms)
{
  int i,idx;
  float x,y,z;
  for (i = 0;i < noatoms;i++)
    {
			idx = i*3;
      x = c[idx]*ele[0][0] + c[idx+1]*ele[0][1] + c[idx+2]*ele[0][2];
      y = c[idx]*ele[1][0] + c[idx+1]*ele[1][1] + c[idx+2]*ele[1][2];
      z = c[idx]*ele[2][0] + c[idx+1]*ele[2][1] + c[idx+2]*ele[2][2];
      c[idx] = x;c[idx+1] = y;c[idx+2] = z;
    }
}

Vector operator *(const Vector &v,const Matrix3x3 &m)
{
  Vector vv;

  vv._vx = v._vx*m.ele[0][0] + v._vy*m.ele[0][1] + v._vz*m.ele[0][2];
  vv._vy = v._vx*m.ele[1][0] + v._vy*m.ele[1][1] + v._vz*m.ele[1][2];
  vv._vz = v._vx*m.ele[2][0] + v._vy*m.ele[2][1] + v._vz*m.ele[2][2];

  return(vv);
}

Vector operator *(const Matrix3x3 &m,const Vector &v)
{
  Vector vv;

  vv._vx = v._vx*m.ele[0][0] + v._vy*m.ele[0][1] + v._vz*m.ele[0][2];
  vv._vy = v._vx*m.ele[1][0] + v._vy*m.ele[1][1] + v._vz*m.ele[1][2];
  vv._vz = v._vx*m.ele[2][0] + v._vy*m.ele[2][1] + v._vz*m.ele[2][2];

  return(vv);
}

Vector &Vector::operator *= (const Matrix3x3 &m)
{
  Vector vv;
  
  vv.SetX(_vx*m.Get(0,0) + _vy*m.Get(0,1) + _vz*m.Get(0,2));
  vv.SetY(_vx*m.Get(1,0) + _vy*m.Get(1,1) + _vz*m.Get(1,2));
  vv.SetZ(_vx*m.Get(2,0) + _vy*m.Get(2,1) + _vz*m.Get(2,2));
  _vx = vv.x();
  _vy = vv.y();
  _vz = vv.z();
  
  return(*this);
}

Matrix3x3 Matrix3x3::invert()
{
  float t[3][3];
  float det;

  det = determinant();

  if (det != 0.0)
  {
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

float Matrix3x3::determinant()
{
  float x,y,z;

  x = ele[0][0] * (ele[1][1] * ele[2][2] - ele[1][2] * ele[2][1]);
  y = ele[0][1] * (ele[1][2] * ele[2][0] - ele[1][0] * ele[2][2]);
  z = ele[0][2] * (ele[1][0] * ele[2][1] - ele[1][1] * ele[2][0]);

  return(x + y + z);
}

Matrix3x3 &Matrix3x3::operator/=(const float &c)
{
  register int i,j;

  for (i = 0;i < 3;i++)
     for (j = 0;j < 3;j++)
	ele[i][j] /= c;

  return(*this);
}

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

void Matrix3x3::FillOrth(float Alpha,float Beta, float Gamma, 
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

ostream& operator<< ( ostream& co, const Matrix3x3& m )

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

void Vector::createOrthoVector(Vector &res) const
{
    Vector cO;

    if ((this->x() == 0.0)&&(this->y() == 0.0))
    {
        if (this->z() == 0.0){
            cerr << "makeorthovec zero vector" << endl;
            exit(0);
        }
        cO.SetX(1.0);
    }
    else
    {
        cO.SetZ(1.0);
    }
    res= cross(cO,*this);
    res.normalize(); 
}

}          
