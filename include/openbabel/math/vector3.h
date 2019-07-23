/**********************************************************************
vector3.h - Handle 3D coordinates.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2006 by Benoit Jacob

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

#ifndef OB_VECTOR_H
#define OB_VECTOR_H

#include <ostream>
#include <math.h>
#include <iostream>

#include <openbabel/babelconfig.h>

#ifndef RAD_TO_DEG
#define RAD_TO_DEG (180.0/M_PI)
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD (M_PI/180.0)
#endif

namespace OpenBabel
{

  class matrix3x3; // declared in math/matrix3x3.h

  // class introduction in vector3.cpp
  class OBAPI vector3
  {
  private :
    double _vx, _vy, _vz ;

  public :
    //! Constructor
    vector3 (const double inX=0.0, const double inY=0.0, const double inZ=0.0):
      _vx(inX), _vy(inY), _vz(inZ)
      {}
    vector3 (double inV[3]):
      _vx(inV[0]), _vy(inV[1]), _vz(inV[2])
      {}
    //! Copy Constructor
    vector3 (const vector3& v):
      _vx(v._vx), _vy(v._vy), _vz(v._vz)
        { }

    //! Destructor
    ~vector3() { }

    //! A random access iterator over x, y, z
    typedef double* iterator;

    //! A random access iterator over const x, y, z
    typedef const double* const_iterator;

    //! A signed integral type for differences between two iterators
    typedef std::ptrdiff_t difference_type;

    //! \return iterator to beginning
    iterator begin() { return &_vx; }

    //! \return iterator to end
    iterator end() { return &_vx + 3; }

    //! /return const_iterator to beginning
    const_iterator begin() const { return &_vx; }

    //! /return const_iterator to end
    const_iterator end() const { return &_vx + 3; }

    //! Set x,y and z-component of a vector
    void Set(const double inX, const double inY, const double inZ)
    {
      _vx = inX;
      _vy = inY;
      _vz = inZ;
    }
    //! Set x,y and z-component of a vector from c[0]..c[2]
    void Set(const double *c)
    {
      _vx = c[0];
      _vy = c[1];
      _vz = c[2];
    }
    //! Access function to set the x-coordinate of the vector
    void SetX(const double inX)
    {
      _vx = inX;
    }
    //! Access function to set the y-coordinate of the vector
    void SetY(const double inY)
    {
      _vy = inY;
    }
    //! Access function to set the z-coordinate of the vector
    void SetZ(const double inZ)
    {
      _vz = inZ;
    }
    //! Access function to get the x-coordinate of the vector
    double GetX() const
    {
      return _vx;
    }
    //! Access function to get the y-coordinate of the vector
    double GetY() const
    {
      return _vy;
    }
    //! Access function to get the z-coordinate of the vector
    double GetZ() const
    {
      return _vz;
    }
    //! \brief Set c[0]..c[2] to the components of the vector
    //! \warning No error checking is performed
    void Get(double *c)
    {
      c[0]=_vx;
      c[1]=_vy;
      c[2]=_vz;
    }
    //! Access function to x: [0], y: [1], and z[2]
    double operator[] ( unsigned int i) const;

    //! Assignment
    vector3& operator= ( const vector3& v)
      {
        _vx = v._vx;
        _vy = v._vy;
        _vz = v._vz;
        return *this;
      }

    //! \return the vector as a const double *
    const double *AsArray() const
    {
      return &_vx;
    }

    //! \brief Vector addition (add @p v to *this)
    //! \return *this + v
    vector3& operator+= ( const vector3& v)
      {
        _vx += v._vx;
        _vy += v._vy;
        _vz += v._vz;
        return *this;
      };
    //! \brief Vector subtraction (subtract @p v from *this)
    //! \return *this - v
    vector3& operator-= ( const vector3& v)
      {
        _vx -= v._vx;
        _vy -= v._vy;
        _vz -= v._vz;
        return *this;
      };
    //! \brief Scalar addition (add @p f to *this)
    //! \return *this + f
    vector3& operator+= ( const double* f)
      {
        _vx += f[0];
        _vy += f[1];
        _vz += f[2];
        return *this;
      };
    //! \brief Scalar subtraction (subtract @p f from *this)
    //! \return *this - f
    vector3& operator-= ( const double* f)
      {
        _vx -= f[0];
        _vy -= f[1];
        _vz -= f[2];
        return *this;
      };
    //! \brief Scalar multiplication (multiply *this by @p c)
    //! \return *this * c
    vector3& operator*= ( const double& c)
      {
        _vx *= c;
        _vy *= c;
        _vz *= c;
        return *this;
      };

    //! \brief Scalar division (divide *this by @p c)
    //! \return *this divided by c
    vector3& operator/= ( const double& c)
      {
        double inv = 1.0 / c;
        return( (*this) *= inv );
      };
    //! Multiplication of matrix and vector
    //! \return the result (i.e., the updated vector)
    //! \todo Currently unimplemented
    vector3& operator*= ( const matrix3x3 &);

    //! Create a random unit vector
    void randomUnitVector();

    //  Member Functions

    //! Scales a vector to give it length one.
    //! \return the result (i.e., the normalized vector)
    vector3& normalize () ;

    //! \return Whether a vector can be normalized
    bool CanBeNormalized () const;

    //! \return The length of the vector squared
    inline double length_2 () const
    {
      return _vx*_vx + _vy*_vy + _vz*_vz;
    };
    //! \return The vector length
    double length () const
    {
      return sqrt( length_2() );
    };
    //! Access function to get the x-coordinate of the vector
    const double & x () const
    {
      return _vx ;
    } ;
    //! Access function to get the y-coordinate of the vector
    const double & y () const
    {
      return _vy ;
    } ;
    //! Access function to get the z-coordinate of the vector
    const double & z () const
    {
      return _vz ;
    } ;
    //! Access function to set the x-coordinate of the vector
    double & x ()
    {
      return _vx ;
    } ;
    //! Access function to set the y-coordinate of the vector
    double & y ()
    {
      return _vy ;
    } ;
    //! Access function to set the z-coordinate of the vector
    double & z ()
    {
      return _vz ;
    } ;

    //! Comparison Methods
    // @{
    //! \brief Equivalence of vectors
    //! \deprecated This method uses unreliable floating point == comparisons
    //!    Use vector3::IsApprox() instead.
    //! \return true if every component is equal
    int operator== ( const vector3& ) const;
    //! \deprecated This method uses unreliable floating point == comparisons
    //!    Use vector3::IsApprox() instead.
    //! \return true if at least one component of the two vectors is !=
    int operator!= ( const vector3& other ) const
    {
      return ! ( (*this) == other );
    }
    //! \brief Safe comparison for floating-point vector3
    //! \return true if the vector *this is approximately equal to the vector
    //!         @p other, to the precision @p precision. More specifically,
    //!         this method works exactly like the OpenBabel::IsApprox()
    //!         function, replacing the absolute value for doubles by the norm
    //!         for vectors.
    //! \param other The vector for comparison
    //! \param precision This parameter plays the same role as in
    //!        OpenBabel::IsApprox().
    bool IsApprox( const vector3 & other, const double & precision ) const;
    //! }@

    //! \return square of the distance between *this and vv
    /*! equivalent to length_2(*this-vv)
     */
    double distSq(const vector3 &vv) const
    {
      double dx = x() - vv.x();
      double dy = y() - vv.y();
      double dz = z() - vv.z();
      return( dx*dx + dy*dy + dz*dz );
    }

    //! Creates a vector of length one, orthogonal to *this.
    //! \return Whether the method was successful
    bool createOrthoVector(vector3 &v) const;

  };

  //! Prints a representation of the vector as a row vector of the form "<0.1,1,2>"
  OBAPI std::ostream& operator<< ( std::ostream&, const vector3& );

  //  Sum, Difference, Scalar Product
  //! Vector addition
  inline OBAPI vector3 operator+ ( const vector3& v1, const vector3& v2)
  {
    return vector3(v1.x()+v2.x(), v1.y()+v2.y(), v1.z()+v2.z());
  }
  //! Vector subtraction
  inline OBAPI vector3 operator- ( const vector3& v1, const vector3& v2)
  {
    return vector3(v1.x()-v2.x(), v1.y()-v2.y(), v1.z()-v2.z());
  }
  //! Unary minus
  inline OBAPI vector3 operator- ( const vector3& v)
  {
    return vector3(-v.x(), -v.y(), -v.z());
  }
  //! Multiplication with a scalar
  inline OBAPI vector3 operator* ( const double& c, const vector3& v)
    {
      return vector3( c*v.x(), c*v.y(), c*v.z());
    }
  //! Multiplication with a scalar
  inline OBAPI vector3 operator* ( const vector3& v, const double& c)
    {
      return vector3( c*v.x(), c*v.y(), c*v.z());
    }
  //! Division by a scalar
  inline OBAPI vector3 operator/ ( const vector3& v, const double& c)
  {
    return vector3( v.x()/c, v.y()/c, v.z()/c);
  }
  // @removed@ misleading operation
  // friend vector3 operator* ( const vector3 &,const vector3 &);

  //vector and matrix ops
  // @removed@ misleading operation; matrix multiplication is not commutitative
  //     friend vector3 operator *(const vector3 &v,const matrix3x3 &m);

  //! Multiplication of matrix and vector
  OBAPI vector3 operator *(const matrix3x3 &m, const vector3 &v);

  //! Dot product of two vectors
  inline OBAPI double dot ( const vector3& v1, const vector3& v2 )
  {
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z() ;
  }
  //! Cross product of two vectors
  OBAPI vector3 cross ( const vector3&, const vector3& );

  //! Calculate the angle between vectors (in degrees)
  OBAPI double vectorAngle ( const vector3& v1, const vector3& v2 );

  //! Calculate the torsion angle between vectors (in degrees)
  OBAPI double CalcTorsionAngle(const vector3 &a, const vector3 &b,
                                        const vector3 &c, const vector3 &d);

  //! Calculate the signed distance of point a to the plane determined by b,c,d
  OBAPI double Point2PlaneSigned(vector3 a, vector3 b, vector3 c, vector3 d);
  //! Calculate the distance of point a to the plane determined by b,c,d
  OBAPI double Point2Plane(vector3 a, vector3 b, vector3 c, vector3 d);
  //! Calculate the angle between point a and the plane determined by b,c,d
  OBAPI double Point2PlaneAngle(const vector3 a, const vector3 b, const vector3 c, const vector3 d);

  //! Calculate the distance of a point a to a line determined by b and c
  OBAPI double Point2Line(const vector3& a, const vector3& b, const vector3& c);

  //  The global constant vector3 objects
  //! The zero vector: <0.0, 0.0, 0.0>
  extern OBAPI const vector3 VZero;
  //! The x unit vector: <1.0, 0.0, 0.0>
  extern OBAPI const vector3 VX;
  //! The y unit vector: <0.0, 1.0, 0.0>
  extern OBAPI const vector3 VY;
  //! The z unit vector: <0.0, 0.0, 1.0>
  extern OBAPI const vector3 VZ;

}

#endif // OB_VECTOR_H

//! \file
//! \brief Handle 3D coordinates.
