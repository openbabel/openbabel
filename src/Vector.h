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

#ifndef VECTOR_H
#define VECTOR_H


#ifdef __sgi
#include <iostream.h>
#include <fstream.h>
#else
#include <iostream>
#include <fstream>
#endif

#include <math.h>
#include "obutil.h"

#ifndef PI
#define PI 3.1415926535897932384626433f
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG 180.0f/PI
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD PI/180.0f
#endif 

namespace OpenBabel {

class Matrix3x3;

class	Vector {

   private :
	
       float		_vx, _vy, _vz ;

   public :

			//  Constructors

     //Vector () {_vx = 0.0;_vy = 0.0;_vz = 0.0;}

      Vector ( const float=0.0f,const float=0.0f,const float=0.0f) ;
      Vector ( const Vector& ) ;//  Copy Constructor
			//  Destructor
      virtual ~Vector () ;

			//  Assignment

      void Set(  const float x, const float y, const float z ) 
			  {_vx = x ;_vy = y ;_vz = z ;}
      void Set(const float *c) {_vx = c[0];_vy = c[1];_vz = c[2];}
      void SetX(const float x) {_vx = x;};
      void SetY(const float y) {_vy = y;};
      void SetZ(const float z) {_vz = z;};
      void Get(float *c) {c[0]=_vx;c[1]=_vy;c[2]=_vz;}

      Vector& operator= ( const Vector& ) ;

			//  Output

      friend std::ostream& operator<< ( std::ostream&, const Vector& ) ;

			//  Comparison

      friend int operator== ( const Vector&, const Vector& ) ;
      friend int operator!= ( const Vector&, const Vector& ) ;

			//  Sum, Difference, Scalar Product

      friend Vector operator+ ( const Vector&, const Vector& ) ;
      friend Vector operator- ( const Vector&, const Vector& ) ;
      friend Vector operator- ( const Vector& ) ;
      friend Vector operator* ( const float&, const Vector& ) ;
      friend Vector operator* ( const Vector&, const float& ) ;
      friend Vector operator/ ( const Vector&, const float& ) ;
      friend Vector operator*(const Vector &,const Vector &);

      //vector and matrix ops

      friend Vector operator *(const Vector &v,const Matrix3x3 &m);
      friend Vector operator *(const Matrix3x3 &m,const Vector &v);


			//  Immediate Sum, Difference, Scalar Product

      Vector& operator+= ( const Vector& ) ;
      Vector& operator-= ( const Vector& ) ;
      Vector& operator*= ( const float& ) ;
      Vector& operator*= (const Matrix3x3 &);
      Vector& operator/= ( const float& ) ;
      Vector& operator+= ( const float* ) ;
      Vector& operator-= ( const float* ) ;


      // create a random unit vector in R3
      void randomUnitVector(OBRandom *oeRand= NULL);

			//  Member Functions
			   //  Dot Product


      friend float dot ( const Vector&, const Vector& ) ;

			   //  Cross Product

      friend Vector cross ( const Vector&, const Vector& ) ;

			   //  Normalization, Make it a unit Vector

      friend float VectorAngle ( const Vector& v1, const Vector& v2 );
			   
			   // calculate angle between vectors

      friend float CalcTorsionAngle(Vector &a,Vector &b,Vector &c,Vector &d);

      Vector& normalize () ;

			   //  Vector Length

      float length () const ;

			   //  Vector Length Squared

      float length_2 () const ;

			//  Access Functions to get 
			//    x-coordinate, y-coordinate or
			//    z-coordinate of the vector

      float x () const { return _vx ; } ;
      float y () const { return _vy ; } ;
      float z () const { return _vz ; } ;

      inline float distSq(const Vector &vv) const { return( (_vx - vv.x() )*(_vx - vv.x() ) + 
                                                            (_vy - vv.y() )*(_vy - vv.y() ) + 
                                                            (_vz - vv.z() )*(_vz - vv.z() ) ); 
      }

      // create a vector orthogonal to me
      void createOrthoVector(Vector &v) const;

   } ;

			//  The global constant Vectors

const Vector VZero ( 0.0, 0.0, 0.0 ) ;
const Vector VX    ( 1.0, 0.0, 0.0 ) ;
const Vector VY    ( 0.0, 1.0, 0.0 ) ;
const Vector VZ    ( 0.0, 0.0, 1.0 ) ;

class Matrix3x3
{
  float ele[3][3];
  public:
  Matrix3x3(void) 
    {
      ele[0][0] = 0.0;ele[0][1] = 0.0;ele[0][2] = 0.0;
      ele[1][0] = 0.0;ele[1][1] = 0.0;ele[1][2] = 0.0;
      ele[2][0] = 0.0;ele[2][1] = 0.0;ele[2][2] = 0.0;
    }

  Matrix3x3(Vector a,Vector b,Vector c)
    {
      ele[0][0] = a.x();ele[0][1] = a.y();ele[0][2] = a.z();
      ele[1][0] = b.x();ele[1][1] = b.y();ele[1][2] = b.z();
      ele[2][0] = c.x();ele[2][1] = c.y();ele[2][2] = c.z();
    }

  Matrix3x3(float d[3][3])
    {
      ele[0][0] = d[0][0];ele[0][1] = d[0][1];ele[0][2] = d[0][2];
      ele[1][0] = d[1][0];ele[1][1] = d[1][1];ele[1][2] = d[1][2];
      ele[2][0] = d[2][0];ele[2][1] = d[2][1];ele[2][2] = d[2][2];
    }

  void GetArray(float *m) 
    {
      m[0] = ele[0][0];m[1] = ele[0][1];m[2] = ele[0][2];
      m[3] = ele[1][0];m[4] = ele[1][1];m[5] = ele[1][2];
      m[6] = ele[2][0];m[7] = ele[2][1];m[8] = ele[2][2];
    }
  Matrix3x3 invert();
  void randomRotation(OBRandom &rnd);
  float determinant();
  float Get(int i,int j) const {return(ele[i][j]);}
  void  Set(int i,int j, float v) {ele[i][j]= v;}
  Matrix3x3 &operator/=(const float &c);
  void SetupRotMat(float,float,float);
  void RotAboutAxisByAngle(const Vector &,const float);
  void RotateCoords(float *,int);
  void FillOrth(float,float,float,float,float,float);
  friend Vector operator *(const Vector &,const Matrix3x3 &);
  friend Vector operator *(const Matrix3x3 &,const Vector &);
  friend std::ostream& operator<< ( std::ostream&, const Matrix3x3 & ) ;
};

Vector center_coords(float*,int);
}

#endif //VECTOR_H
