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

class matrix3x3;

/*!
 * \brief Represents a vector in the 3-dimensional real space.
 */

 class	vector3 {
   private :
     float		_vx, _vy, _vz ;

   public :
     //! Constructor
     vector3 (const float x=0.0f, const float y=0.0f, const float z=0.0f) {_vx = x; _vy = y; _vz = z;};
     //! Copy Constructor
     vector3 (const vector3& v) { _vx = v._vx; _vy = v._vy; _vz = v._vz; }; 

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
     //! access function
     /*! This (slow) method allows to access the elements of the
       vector as if it were an array of floats. If the index is > 2,
       then a warning is printed, and the program is terminated via
       exit(-1). Otherwise, if i is 0, 1 or 2, then a reference to x,
       y or z is returned, respectively.
       
       \warning This method is primarily designed to facilitate the
       integration ('OpenBabelization') of code that uses arrays of
       floats rather than the vector class. Due to the error checks
       the method is of course very slow and should therefore be
       avoided in production code.
     */
     float& operator[] ( unsigned int i);

     //! assignment
     vector3& operator= ( const vector3& v) {_vx = v._vx; _vy = v._vy; _vz = v._vz; return *this;};

     //! prints a representation of the vector as a row vector of the form "<0.1,1,2>"
     friend std::ostream& operator<< ( std::ostream&, const vector3& ) ;

     //  Comparison
     friend int operator== ( const vector3&, const vector3& ) ;
     friend int operator!= ( const vector3&, const vector3& ) ;

     //  Sum, Difference, Scalar Product
     //! vector addition
     friend vector3 operator+ ( const vector3& v1, const vector3& v2) { return vector3(v1._vx+v2._vx, v1._vy+v2._vy, v1._vz+v2._vz); };
     //! vector subtraction
     friend vector3 operator- ( const vector3& v1, const vector3& v2) { return vector3(v1._vx-v2._vx, v1._vy-v2._vy, v1._vz-v2._vz); };
     //! unary minus
     friend vector3 operator- ( const vector3& v) { return vector3(-v._vx, -v._vy, -v._vz); };
     //! multiplication with a scalar
     friend vector3 operator* ( const float& c, const vector3& v) { return vector3( c*v._vx, c*v._vy, c*v._vz); };
     //! multiplication with a scalar
     friend vector3 operator* ( const vector3& v, const float& c) { return vector3( c*v._vx, c*v._vy, c*v._vz); };
     //! division by a scalar
     friend vector3 operator/ ( const vector3& v, const float& c) { return vector3( v._vx/c, v._vy/c, v._vz/c); };
     // @removed@ misleading operation
     // friend vector3 operator* ( const vector3 &,const vector3 &);

     //vector and matrix ops
     // @removed@ misleading operation; matrix multiplication is not commutitative
     //     friend vector3 operator *(const vector3 &v,const matrix3x3 &m);

     //! multiplication of matrix and vector
     /*! calculates the product m*v of the matrix m and the column
         vector represented by v
      */
     friend vector3 operator *(const matrix3x3 &m,const vector3 &v);
     
     //  Immediate Sum, Difference, Scalar Product
     vector3& operator+= ( const vector3& v) {_vx += v._vx; _vy += v._vy; _vz += v._vz; return *this;};
     vector3& operator-= ( const vector3& v) {_vx -= v._vx; _vy -= v._vy; _vz -= v._vz; return *this;};
     vector3& operator+= ( const float* f)  {_vx += f[0]; _vy += f[1]; _vz += f[2]; return *this;};
     vector3& operator-= ( const float* f)  {_vx -= f[0]; _vy -= f[1]; _vz -= f[2]; return *this;};
     vector3& operator*= ( const float& c)  {_vx *= c; _vy *= c; _vz *= c; return *this; };
     vector3& operator/= ( const float& c)  {_vx /= c; _vy /= c; _vz /= c; return *this; };
     //! multiplication of matrix and vector
     /*! calculates the product m*(*this) of the matrix m and the
         column vector represented by *this
      */
     vector3& operator*= ( const matrix3x3 &);

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
     friend float dot ( const vector3&, const vector3& ) ;
     
     //! cross product of two vectors
     friend vector3 cross ( const vector3&, const vector3& ) ;
     
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
     friend float vectorAngle ( const vector3& v1, const vector3& v2 );
			   
     //! calculate the torsion angle between vectors
     friend float CalcTorsionAngle(const vector3 &a, const vector3 &b,
				    const vector3 &c, const vector3 &d);

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
      vector3& normalize () ;

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
      inline float distSq(const vector3 &vv) const 
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
      void createOrthoVector(vector3 &v) const;

   } ;

			//  The global constant vector3s

const vector3 VZero ( 0.0f, 0.0f, 0.0f ) ;
const vector3 VX    ( 1.0f, 0.0f, 0.0f ) ;
const vector3 VY    ( 0.0f, 1.0f, 0.0f ) ;
const vector3 VZ    ( 0.0f, 0.0f, 1.0f ) ;


vector3 center_coords(float*,int);
}

#endif // OB_VECTOR_H
