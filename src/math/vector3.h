/**********************************************************************
vector3.h - Handle 3D coordinates.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#include <ostream>
#include <math.h>

#ifndef PI
#define PI 3.1415926535897932384626433
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG 180.0/PI
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD PI/180.0
#endif

namespace OpenBabel
{

  class matrix3x3; // declared in math/matrix3x3.h
  class OBRandom; // declared in obutil.h

  // class introduction in vector3.cpp
  class	OBAPI vector3
    {
      private :
	double		_vx, _vy, _vz ;

      public :
	//! Constructor
	vector3 (const double x=0.0, const double y=0.0, const double z=0.0)
	{
	  _vx = x;
	  _vy = y;
	  _vz = z;
	};
      //! Copy Constructor
      vector3 (const vector3& v)
	{
	  _vx = v._vx;
	  _vy = v._vy;
	  _vz = v._vz;
	};

      //! set x,y and z-component of a vector
      void Set(const double x, const double y, const double z)
	{
	  _vx = x ;
	  _vy = y ;
	  _vz = z ;
	};
      //! set x,y and z-component of a vector from c[0]..c[2]
      void Set(const double *c)
	{
	  _vx = c[0];
	  _vy = c[1];
	  _vz = c[2];
	}
      //! access function to get the x-coordinate of the vector
      void SetX(const double x)
	{
	  _vx = x;
	};
      //! access function to get the y-coordinate of the vector
      void SetY(const double y)
	{
	  _vy = y;
	};
      //! access function to get the z-coordinate of the vector
      void SetZ(const double z)
	{
	  _vz = z;
	};
      //! set c[0]..c[2] to the components of the vector
      void Get(double *c)
	{
	  c[0]=_vx;
	  c[1]=_vy;
	  c[2]=_vz;
	};
      //! access function
      double& operator[] ( unsigned int i);

      //! assignment
      vector3& operator= ( const vector3& v)
	{
	  _vx = v._vx;
	  _vy = v._vy;
	  _vz = v._vz;
	  return *this;
	};

      //! prints a representation of the vector as a row vector of the form "<0.1,1,2>"
      friend OBAPI std::ostream& operator<< ( std::ostream&, const vector3& ) ;

      //  Comparison
      friend OBAPI int operator== ( const vector3&, const vector3& ) ;
      friend OBAPI int operator!= ( const vector3&, const vector3& ) ;

      //  Sum, Difference, Scalar Product
      //! vector addition
      friend OBAPI vector3 operator+ ( const vector3& v1, const vector3& v2)
	{
	  return vector3(v1._vx+v2._vx, v1._vy+v2._vy, v1._vz+v2._vz);
	};
      //! vector subtraction
      friend OBAPI vector3 operator- ( const vector3& v1, const vector3& v2)
	{
	  return vector3(v1._vx-v2._vx, v1._vy-v2._vy, v1._vz-v2._vz);
	};
      //! unary minus
      friend OBAPI vector3 operator- ( const vector3& v)
	{
	  return vector3(-v._vx, -v._vy, -v._vz);
	};
      //! multiplication with a scalar
      friend OBAPI vector3 operator* ( const double& c, const vector3& v)
	{
	  return vector3( c*v._vx, c*v._vy, c*v._vz);
	};
      //! multiplication with a scalar
      friend OBAPI vector3 operator* ( const vector3& v, const double& c)
	{
	  return vector3( c*v._vx, c*v._vy, c*v._vz);
	};
      //! division by a scalar
      friend OBAPI vector3 operator/ ( const vector3& v, const double& c)
	{
	  return vector3( v._vx/c, v._vy/c, v._vz/c);
	};
      // @removed@ misleading operation
      // friend vector3 operator* ( const vector3 &,const vector3 &);

      //vector and matrix ops
      // @removed@ misleading operation; matrix multiplication is not commutitative
      //     friend vector3 operator *(const vector3 &v,const matrix3x3 &m);

      //! multiplication of matrix and vector
      friend OBAPI vector3 operator *(const matrix3x3 &m,const vector3 &v);

      //  Immediate Sum, Difference, Scalar Product
      vector3& operator+= ( const vector3& v)
	{
	  _vx += v._vx;
	  _vy += v._vy;
	  _vz += v._vz;
	  return *this;
	};
      vector3& operator-= ( const vector3& v)
	{
	  _vx -= v._vx;
	  _vy -= v._vy;
	  _vz -= v._vz;
	  return *this;
	};
      vector3& operator+= ( const double* f)
	{
	  _vx += f[0];
	  _vy += f[1];
	  _vz += f[2];
	  return *this;
	};
      vector3& operator-= ( const double* f)
	{
	  _vx -= f[0];
	  _vy -= f[1];
	  _vz -= f[2];
	  return *this;
	};
      vector3& operator*= ( const double& c)
	{
	  _vx *= c;
	  _vy *= c;
	  _vz *= c;
	  return *this;
	};
      vector3& operator/= ( const double& c)
	{
	  _vx /= c;
	  _vy /= c;
	  _vz /= c;
	  return *this;
	};
      //! multiplication of matrix and vector
      vector3& operator*= ( const matrix3x3 &);

      //! create a random unit vector
      void randomUnitVector(OBRandom *oeRand= 0L);

      //  Member Functions

      //! dot product of two vectors
      friend OBAPI double dot ( const vector3&, const vector3& ) ;

      //! cross product of two vectors
      friend OBAPI vector3 cross ( const vector3&, const vector3& ) ;

      //! calculate angle between vectors
      friend OBAPI double vectorAngle ( const vector3& v1, const vector3& v2 );

      //! calculate the torsion angle between vectors
      friend OBAPI double CalcTorsionAngle(const vector3 &a, const vector3 &b,
					   const vector3 &c, const vector3 &d);

      //! scales a vector to give it length one.
      vector3& normalize () ;

      //! vector length
      double length () const
	{
	  return sqrt(_vx*_vx + _vy*_vy + _vz*_vz);
	};
      //! vector length squared
      double length_2 () const
	{
	  return _vx*_vx + _vy*_vy + _vz*_vz;
	};
      //! access function to get the x-coordinate of the vector
      double x () const
	{
	  return _vx ;
	} ;
      //! access function to get the y-coordinate of the vector
      double y () const
	{
	  return _vy ;
	} ;
      //! access function to get the z-coordinate of the vector
      double z () const
	{
	  return _vz ;
	} ;

      //! square to the distance between *this and vv
      /*! equivalent to length_2(*this-vv)
       */
      inline double distSq(const vector3 &vv) const
	{
	  return( (_vx - vv.x() )*(_vx - vv.x() ) +
		  (_vy - vv.y() )*(_vy - vv.y() ) +
		  (_vz - vv.z() )*(_vz - vv.z() ) );
	}

      //! creates a vector of length one, orthogonal to *this.
      void createOrthoVector(vector3 &v) const;

    } ;

  //! \brief Calculate the distance of point a to the plane determined by b,c,d
  OBAPI double Point2Plane(vector3 a, vector3 b, vector3 c, vector3 d);

  //  The global constant vector3s
  extern OBAPI const vector3 VZero;
  extern OBAPI const vector3 VX;
  extern OBAPI const vector3 VY;
  extern OBAPI const vector3 VZ;

#ifndef SWIG
  OBAPI vector3 center_coords(double*,int);
#endif
}

#endif // OB_VECTOR_H

//! \file
//! \brief Handle 3D coordinates.
