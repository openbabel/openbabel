/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2002 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_VECTOR_H
#define OB_VECTOR_H


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

/*!
 * \brief Represents a vector in the 3-dimensional real space.
 */

 class	Vector {
   private :
     float		_vx, _vy, _vz ;

   public :
     //! Constructor
     Vector (const float x=0.0f, const float y=0.0f, const float z=0.0f) {_vx = x; _vy = y; _vz = z;};
     //! Copy Constructor
     Vector (const Vector& v) { _vx = v._vx; _vy = v._vy; _vz = v._vz; }; 

     //! set x,y and z-component of a vector
     void Set(const float x, const float y, const float z) {_vx = x ;_vy = y ;_vz = z ;};
     //! set x,y and z-component of a vector from c[0]..c[2]
     void Set(const float *c) {_vx = c[0];_vy = c[1];_vz = c[2];}
     //! access function to get the x-coordinate of the vector
     void SetX(const float x) {_vx = x;};
     //! access function to get the y-coordinate of the vector
     void SetY(const float y) {_vy = y;};
     //! access function to get the z-coordinate of the vector
     void SetZ(const float z) {_vz = z;};
     //! set c[0]..c[2] to the components of the vector
     void Get(float *c) {c[0]=_vx; c[1]=_vy; c[2]=_vz;};

     //! assignment
     Vector& operator= ( const Vector& v) {_vx = v._vx; _vy = v._vy; _vz = v._vz; return *this;};

     //! prints a representation of the vector as a row vector of the form "<0.1,1,2>"
     friend std::ostream& operator<< ( std::ostream&, const Vector& ) ;

     //  Comparison
     friend int operator== ( const Vector&, const Vector& ) ;
     friend int operator!= ( const Vector&, const Vector& ) ;

     //  Sum, Difference, Scalar Product
     //! vector addition
     friend Vector operator+ ( const Vector& v1, const Vector& v2) { return Vector(v1._vx+v2._vx, v1._vy+v2._vy, v1._vz+v2._vz); };
     //! vector subtraction
     friend Vector operator- ( const Vector& v1, const Vector& v2) { return Vector(v1._vx-v2._vx, v1._vy-v2._vy, v1._vz-v2._vz); };
     //! unary minus
     friend Vector operator- ( const Vector& v) { return Vector(-v._vx, -v._vy, -v._vz); };
     //! multiplication with a scalar
     friend Vector operator* ( const float& c, const Vector& v) { return Vector( c*v._vx, c*v._vy, c*v._vz); };
     //! multiplication with a scalar
     friend Vector operator* ( const Vector& v, const float& c) { return Vector( c*v._vx, c*v._vy, c*v._vz); };
     //! division by a scalar
     friend Vector operator/ ( const Vector& v, const float& c) { return Vector( v._vx/c, v._vy/c, v._vz/c); };
     // @removed@ misleading operation
     // friend Vector operator* ( const Vector &,const Vector &);

     //vector and matrix ops
     // @removed@ misleading operation; matrix multiplication is not commutitative
     //     friend Vector operator *(const Vector &v,const Matrix3x3 &m);

     //! multiplication of matrix and vector
     /*! calculates the product m*v of the matrix m and the column
         vector represented by v
      */
     friend Vector operator *(const Matrix3x3 &m,const Vector &v);
     
     //  Immediate Sum, Difference, Scalar Product
     Vector& operator+= ( const Vector& v) {_vx += v._vx; _vy += v._vy; _vz += v._vz; return *this;};
     Vector& operator-= ( const Vector& v) {_vx -= v._vx; _vy -= v._vy; _vz -= v._vz; return *this;};
     Vector& operator+= ( const float* f)  {_vx += f[0]; _vy += f[1]; _vz += f[2]; return *this;};
     Vector& operator-= ( const float* f)  {_vx -= f[0]; _vy -= f[1]; _vz -= f[2]; return *this;};
     Vector& operator*= ( const float& c)  {_vx *= c; _vy *= c; _vz *= c; return *this; };
     Vector& operator/= ( const float& c)  {_vx /= c; _vy /= c; _vz /= c; return *this; };
     //! multiplication of matrix and vector
     /*! calculates the product m*(*this) of the matrix m and the
         column vector represented by *this
      */
     Vector& operator*= ( const Matrix3x3 &);

     //! create a random unit vector
     /*! replaces *this with a random unit vector, which is (supposed
       to be) uniformly distributed over the unit sphere. Uses the
       random number generator oeRand, or uses the system number
       generator with a time seed if oeRand == NULL.
       
       @param oeRand random number generator to use, or 0L, if the
       system random number generator (with time seed) should be used
     */
     void randomUnitVector(OBRandom *oeRand= 0L);
     
     //  Member Functions
      
     //! dot product of two vectors
     friend float dot ( const Vector&, const Vector& ) ;
     
     //! cross product of two vectors
     friend Vector cross ( const Vector&, const Vector& ) ;
     
     //! calculate angle between vectors
     /*! This method calculates the angle between two vectors
       
         \warning If length() of any of the two vectors is == 0.0f,
         this method will divide by zero. If the product of the
         length() of the two vectors is very close to 0,0f, but not ==
         0.0f, this method may behave in unexpected ways and return
         almost random results; details may depend on your particular
         floating point implementation. The use of this method is
         therefore highly discouraged, unless you are certain that the
         length()es are in a reasonable range, away from 0.0f (Stefan
         Kebekus)

	 \deprecated This method will probably replaced by a safer
	 algorithm in the future.

	 \todo Replace this method with a more fool-proof version.

         @returns the angle in degrees (0-360)
     */
     friend float VectorAngle ( const Vector& v1, const Vector& v2 );
			   
     //! calculate the torsion angle between vectors
     friend float CalcTorsionAngle(const Vector &a, const Vector &b,
				    const Vector &c, const Vector &d);

      //! scales a vector to give it length one.
      /*! This method checks if the current vector has length() ==
        0.0f.  If so, *this remains unchanged. Otherwise, *this is
        scaled by 1.0/length().

	\warning If length() is very close to zero, but not == 0.0f,
	this method may behave in unexpected ways and return almost
	random results; details may depend on your particular floating
	point implementation. The use of this method is therefore
	highly discouraged, unless you are certain that length() is in
	a reasonable range, away from 0.0f (Stefan Kebekus)

	\deprecated This method will probably replaced by a safer
	algorithm in the future.

	\todo Replace this method with a more fool-proof version.

        @returns a reference to *this
      */
      Vector& normalize () ;

      //! vector length
      float length () const { return sqrt(_vx*_vx + _vy*_vy + _vz*_vz);};
      //! vector length squared
      float length_2 () const { return _vx*_vx + _vy*_vy + _vz*_vz; };
      //! access function to get the x-coordinate of the vector
      float x () const { return _vx ; } ;
      //! access function to get the y-coordinate of the vector
      float y () const { return _vy ; } ;
      //! access function to get the z-coordinate of the vector
      float z () const { return _vz ; } ;

      //! square to the distance between *this and vv
      /*! equivalent to length_2(*this-vv)
       */
      inline float distSq(const Vector &vv) const 
	{ return( (_vx - vv.x() )*(_vx - vv.x() ) + 
		  (_vy - vv.y() )*(_vy - vv.y() ) + 
		  (_vz - vv.z() )*(_vz - vv.z() ) ); 
	}
      
      //! creates a vector of length one, orthogonal to *this.
      /*! This method checks if the current vector *this is zero
        (i.e. if all entries == 0.0f). If so, a warning message is
        printed, and the whole program is aborted with exit(0).
        Otherwise, a vector of length one is generated, which is
        orthogonal to *this, and stored in v. The resulting vector is
        not random.

	\warning If the entries of the *this (in particular the
	z-component) are very close to zero, but not == 0.0f, this
	method may behave in unexpected ways and return almost random
	results; details may depend on your particular floating point
	implementation. The use of this method is therefore highly
	discouraged, unless you are certain that all components of
	*this are in a reasonable range, away from 0.0f (Stefan
	Kebekus)

	\deprecated This method will probably replaced by a safer
	algorithm in the future.

	\todo Replace this method with a more fool-proof version that
	does not call exit()


       
        @param v a reference to a vector where the result will be
        stored
      */
      void createOrthoVector(Vector &v) const;

   } ;

			//  The global constant Vectors

const Vector VZero ( 0.0f, 0.0f, 0.0f ) ;
const Vector VX    ( 1.0f, 0.0f, 0.0f ) ;
const Vector VY    ( 0.0f, 1.0f, 0.0f ) ;
const Vector VZ    ( 0.0f, 0.0f, 1.0f ) ;

class Matrix3x3
{
  float ele[3][3];
  public:
  Matrix3x3(void) 
    {
      ele[0][0] = 0.0f; ele[0][1] = 0.0f; ele[0][2] = 0.0f; 
      ele[1][0] = 0.0f; ele[1][1] = 0.0f; ele[1][2] = 0.0f; 
      ele[2][0] = 0.0f; ele[2][1] = 0.0f; ele[2][2] = 0.0f; 
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

  
  //! Replaces *this with it's inverse.
  /*! This method checks if the determinant of *this is == 0.0f. If
      so, nothing is done and no warning is issued. Otherwise, the
      inverse matrix is calculated, and *this is overwritten with the
      result.

      \warning If the determinant is close to zero, but not == 0.0f,
      this method may behave in unexpected ways and return almost
      random results; details may depend on your particular floating
      point implementation. The use of this method is therefore highly
      discouraged, unless you are certain that the determinant is in a
      reasonable range, away from 0.0f (Stefan Kebekus)

      \deprecated This method will probably replaced by a safer
      algorithm in the future.

      \todo Replace this method with a more fool-proof version.
  */
  Matrix3x3 invert();

  void randomRotation(OBRandom &rnd);
  //! returns the determinant of the matrix
  float determinant();
  //! access function
  /*! \warning i or j are not in the range 0..2, random results are
      returned, and your program may even segfault. (Stefan Kebekus)

      \todo Replace this method with a more fool-proof version.
  */
  float Get(int i,int j) const {return(ele[i][j]);}
  //! access function
  /*! \warning i or j are not in the range 0..2, random variable are
      overwritten, and your program may segfault. (Stefan Kebekus)

      \todo Replace this method with a more fool-proof version.
  */
  void  Set(int i,int j, float v) {ele[i][j]= v;}
  Matrix3x3 &operator/=(const float &c);
  void SetupRotMat(float,float,float);
  //! calculates rotation matrix
  /*! Replaces *this with a matrix that represents rotation about the
      axis by a an angle. 
   
      \warning If the vector axis has length zero, this method will
      generate the 0-matrix. If the length of the axis is close to
      zero, but not == 0.0f, this method may behave in unexpected ways
      and return almost random results; details may depend on your
      particular floating point implementation. The use of this method
      is therefore highly discouraged, unless you are certain that the
      length is in a reasonable range, away from 0.0f (Stefan
      Kebekus)

      \deprecated This method will probably replaced by a safer
      algorithm in the future.

      \todo Replace this method with a more fool-proof version.

      @param axis specifies the axis of the rotation
      @angle angle in degrees (0..360)
  */
  void RotAboutAxisByAngle(const Vector &axis, const float angle);
  void RotateCoords(float *,int);
  void FillOrth(float,float,float,float,float,float);
  friend Vector operator *(const Vector &,const Matrix3x3 &);
  friend Vector operator *(const Matrix3x3 &,const Vector &);
  friend std::ostream& operator<< ( std::ostream&, const Matrix3x3 & ) ;
};

Vector center_coords(float*,int);
}

#endif // OB_VECTOR_H
