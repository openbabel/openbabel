/**********************************************************************
vector3.cpp - Handle 3D coordinates.

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

#include <openbabel/babelconfig.h>

#include <iostream>
#include <float.h>

#include <openbabel/math/vector3.h>
#include "../rand.h"
#include <openbabel/obutil.h>

using namespace std;

namespace OpenBabel
{

  /*! \class vector3 vector3.h <openbabel/math/vector3.h>
    \brief Represents a vector in 3-dimensional real space.

    The vector3 class was designed to simplify operations with floating
    point coordinates. To this end many of the common operations have been
    overloaded for simplicity. Vector addition, subtraction, scalar
    multiplication, dot product, cross product, magnitude and a number of
    other utility functions are built in to the vector class. For a full
    description of the class member functions please consult the
    documentation. The following code demonstrates several of the
    functions of the vector class:
    \code
    vector3 v1,v2,v3;
    v1 = VX;
    v2 = VY;
    v3 = cross(v1,v2);
    v3 *= 2.5;
    v3.normalize();
    \endcode
  */

  /*! This (slow) method allows one to access the elements of the
    vector as if it were an array of doubles. If the index is > 2,
    then a warning is printed and 0.0 is returned.
    Otherwise, if i is 0, 1 or 2, then a reference to x,
    y or z is returned, respectively.

    \warning This method is primarily designed to facilitate the
    integration ('Open Babelization') of code that uses arrays of
    doubles rather than the vector class. Due to the error checks
    the method is of course very slow and should therefore be
    avoided in production code.
  */
  double vector3::operator[] ( unsigned int i) const
  {
#ifdef OB_OLD_MATH_CHECKS
    if (i > 2)
      {
        cerr << "ERROR in OpenBabel::vector3::operator[]" << endl
             << "The method has been called with an illegal index i=" << i << "." << endl
             << "Please contact the author of the offending program immediately." << endl;
        return 0.0;
      }
#endif
    if (i == 0)
      return _vx;
    if (i == 1)
      return _vy;
    else return _vz;
  }

  /*! Replaces *this with a random unit vector, which is (supposed
    to be) uniformly distributed over the unit sphere. Uses the
    system number generator with a time seed.

  */
  void vector3::randomUnitVector()
  {
    OBRandom *ptr;
    static OBRandom singleRand(true);
    ptr = &singleRand;

    // obtain a random vector with 0.001 <= length^2 <= 1.0, normalize
    // the vector to obtain a random vector of length 1.0.
    double l;
    do
      {
        this->Set(ptr->NextFloat()-0.5, ptr->NextFloat()-0.5, ptr->NextFloat()-0.5);
        l = length_2();
      }
    while ( (l > 1.0) || (l < 1e-4) );
    this->normalize();
  }

  OBAPI ostream& operator<< ( ostream& co, const vector3& v )
  {
    co << "< " << v.x() << ", " << v.y() << ", " << v.z() << " >" ;
    return co ;
  }

  OBAPI int vector3::operator== ( const vector3& other ) const
  {
    return ( ( x() == other.x() ) &&
             ( y() == other.y() ) &&
             ( z() == other.z() ) );
  }

  bool vector3::IsApprox(const vector3 & other, const double & precision) const
  {
    return( distSq( other )
            <= precision * precision
               * std::min( length_2(), other.length_2() ) );
  }

  /*! This method returns true if *this can be safely normalized.
   * Vectors that can't be safely normalized are:
   * - the zero vector (0,0,0)
   * - vectors having coords that can't be squared without triggering
   * an overflow or underflow. This means doubles having absolute
   * value greater than 1e150 or smaller than 1e-150.
   */
  bool vector3::CanBeNormalized () const
  {
    if( _vx == 0.0 && _vy == 0.0 && _vz == 0.0 ) return false;
    return( CanBeSquared(_vx)
         && CanBeSquared(_vy)
         && CanBeSquared(_vz) );
  }

  /*! This method normalizes *this. In other words, it divides
   * the x,y,z coords of *this by this->length().
   * If *this can't be safely normalized, it remains unchanged.
   * See CanBeNormalized().

   @returns a reference to *this

   */
  vector3& vector3 :: normalize ()
  {
#ifdef OB_OLD_MATH_CHECKS
    if( CanBeNormalized() )
      (*this) /= length();
#else
    (*this) /= length();
#endif
    return(*this);
  }

  OBAPI vector3 cross ( const vector3& v1, const vector3& v2 )
  {
    vector3 vv ;

    vv.x() =   v1.y()*v2.z() - v1.z()*v2.y() ;
    vv.y() = - v1.x()*v2.z() + v1.z()*v2.x() ;
    vv.z() =   v1.x()*v2.y() - v1.y()*v2.x() ;

    return ( vv ) ;
  }


  /*! This method calculates the angle between two vectors

  \warning If length() of any of the two vectors is == 0.0,
  this method will divide by zero. If the product of the
  length() of the two vectors is very close to 0.0, but not ==
  0.0, this method may behave in unexpected ways and return
  almost random results; details may depend on your particular
  floating point implementation. The use of this method is
  therefore highly discouraged, unless you are certain that the
  length()es are in a reasonable range, away from 0.0 (Stefan
  Kebekus)

  \deprecated This method will probably replaced by a safer
  algorithm in the future.

  \todo Replace this method with a more fool-proof version.

  @returns the angle in degrees (0-360)
  */
  OBAPI double vectorAngle ( const vector3& v1, const vector3& v2 )
  {
    double dp;

    dp = dot(v1,v2)/ ( v1.length() * v2.length() );


    if (dp < -0.999999)
      dp = -0.9999999;

    if (dp > 0.9999999)
      dp = 0.9999999;


    return((RAD_TO_DEG * acos(dp)));
  }

  /*!  This function calculates the torsion angle of three vectors, represented
    by four points A--B--C--D, i.e. B and C are vertexes, but none of A--B,
    B--C, and C--D are colinear.  A "torsion angle" is the amount of "twist"
    or torsion needed around the B--C axis to bring A--B into the same plane
    as B--C--D.  The torsion is measured by "looking down" the vector B--C so
    that B is superimposed on C, then noting how far you'd have to rotate
    A--B to superimpose A over D.  Angles are + in theanticlockwise
    direction.  The operation is symmetrical in that if you reverse the image
    (look from C to B and rotate D over A), you get the same answer.
  */

  OBAPI double CalcTorsionAngle(const vector3 &a, const vector3 &b,
                                const vector3 &c, const vector3 &d)
  {

    double torsion;
    vector3 b1,b2,b3,c1,c2,c3;

    b1 = a - b;
    b2 = b - c;
    b3 = c - d;

#ifdef OB_OLD_MATH_CHECKS
    c1 = cross(b1,b2);
    c2 = cross(b2,b3);
    c3 = cross(c1,c2);


    if (c1.length() * c2.length() < 0.001)
    {
      torsion = 0.0;
      return torsion;
    }
#endif

    double rb2 = sqrt(dot(b2, b2));

    vector3 b2xb3 = cross(b2, b3);
    vector3 b1xb2 = cross(b1, b2);
    torsion = - atan2(dot(rb2 * b1, b2xb3), dot(b1xb2, b2xb3));

    return(torsion * RAD_TO_DEG);
  }

  /*! \brief Construct a unit vector orthogonal to *this.

   It requires that *this is normalizable; otherwise it just
   returns false. See CanBeNormalized()

   @param res reference by which to pass the result.

   @returns always true. (Return value kept for compatibility,
            as old versions of OpenBabel used to check for
            normalizability).
  */
  bool vector3::createOrthoVector(vector3 &res) const
  {
#ifdef OB_OLD_MATH_CHECKS
    // sanity check
    if( ! CanBeNormalized() ) return false;
#endif

    /* Let us compute the crossed product of *this with a vector
       that is not too close to being colinear to *this.
    */

    /* unless the _vx and _vy coords are both close to zero, we can
     * simply take ( -_vy, _vx, 0 ) and normalize it.
     */
    if( ! IsNegligible( _vx, _vz ) || ! IsNegligible( _vy, _vz ) )
    {
      double norm = sqrt( _vx*_vx + _vy*_vy );
      res._vx = -_vy / norm;
      res._vy = _vx / norm;
      res._vz = 0.0;
    }
    /* if both _vx and _vy are close to zero, then the vector is close
     * to the z-axis, so it's far from colinear to the x-axis for instance.
     * So we take the crossed product with (1,0,0) and normalize it.
     */
    else
    {
      double norm = sqrt( _vy*_vy + _vz*_vz );
      res._vx = 0.0;
      res._vy = -_vz / norm;
      res._vz = _vy / norm;
    }

    return true;
  }

  const vector3 VZero ( 0.0, 0.0, 0.0 ) ;
  const vector3 VX    ( 1.0, 0.0, 0.0 ) ;
  const vector3 VY    ( 0.0, 1.0, 0.0 ) ;
  const vector3 VZ    ( 0.0, 0.0, 1.0 ) ;

  /* Calculate the signed distance of point a to the plane determined by b,c,d */
  double Point2PlaneSigned(vector3 a, vector3 b, vector3 c, vector3 d)
  {
    vector3 v_ba = a-b;
    vector3 v_normal = cross(c-b, d-b);
    return dot( v_normal, v_ba ) / v_normal.length();
  }

  /* Calculate the distance of point a to the plane determined by b,c,d */
  double Point2Plane(vector3 a, vector3 b, vector3 c, vector3 d)
  {
    return fabs( Point2PlaneSigned(a, b, c, d) );
  }

  /* Calculate the angle between point a and the plane determined by b,c,d */
  double Point2PlaneAngle(const vector3 a, const vector3 b, const vector3 c, const vector3 d)
  {
    vector3 ac, bc, cd, normal;
    double angle;

    ac = a - c;
    bc = b - c;
    cd = c - d;

    normal = cross(bc, cd);
    angle = 90.0 - vectorAngle(normal, ac);

    return angle;
  }

  double Point2Line(const vector3& a, const vector3& b, const vector3& c)
  {
    //http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    vector3 v_cb = c-b;
    vector3 v_normal = cross(a-b, a-c);
    return fabs(v_normal.length() / v_cb.length() );
  }

} // namespace OpenBabel

//! \file vector3.cpp
//! \brief Handle 3D coordinates.
