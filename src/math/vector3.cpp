/**********************************************************************
vector3.cpp - Handle 3D coordinates.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

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

#include <math.h>

#include "mol.h"
#include "math/vector3.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

using namespace std;

namespace OpenBabel {

  /*! \class vector3
     \brief Represents a vector in the 3-dimensional real space.

The vector3 class was designed to simplify operations with doubleing
point coordinates. To this end many of the common operations have been
overloaded for simplicity. Vector addition, subtraction, scalar
multiplication, dot product, cross product, magnitude and a number of
other utility functions are built in to the vector class. For a full
description of the class member functions please consult the header
file vector3.h. The following code demonstrates several of the
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

     /*! This (slow) method allows to access the elements of the
       vector as if it were an array of doubles. If the index is > 2,
       then a warning is printed, and the program is terminated via
       exit(-1). Otherwise, if i is 0, 1 or 2, then a reference to x,
       y or z is returned, respectively.
       
       \warning This method is primarily designed to facilitate the
       integration ('Open Babelization') of code that uses arrays of
       doubles rather than the vector class. Due to the error checks
       the method is of course very slow and should therefore be
       avoided in production code.
     */
  double& vector3::operator[] ( unsigned int i)
  {
    if (i > 2) {
      cerr << "ERROR in OpenBabel::vector3::operator[]" << endl
	   << "The method has been called with an illegal index i=" << i << "." << endl
	   << "Please contact the author of the offending program immediately." << endl;
      exit(-1);
    }
    if (i == 0)
      return _vx;
    if (i == 1)
      return _vy;
    return _vz;
  }

  /*! replaces *this with a random unit vector, which is (supposed
    to be) uniformly distributed over the unit sphere. Uses the
    random number generator obRand, or uses the system number
    generator with a time seed if obRand == NULL.
       
    @param obRandP random number generator to use, or 0L, if the
    system random number generator (with time seed) should be used
  */
  void vector3::randomUnitVector(OBRandom *obRandP)
  {
    OBRandom *ptr;
    if (!obRandP) {
      ptr = new OBRandom(true);
      ptr->TimeSeed();
    } else
      ptr = obRandP;
    
    // obtain a random vector with 0.001 <= length^2 <= 1.0, normalize
    // the vector to obtain a random vector of length 1.0.
    double l;
    do {
      this->Set(ptr->NextFloat()-0.5, ptr->NextFloat()-0.5, ptr->NextFloat()-0.5);
      l = length_2();
    } while ( (l > 1.0) || (l < 1e-4) );
    this->normalize();
    
    if (!obRandP) 
      delete ptr;
  }

  ostream& operator<< ( ostream& co, const vector3& v )
  {
    co << "< " << v._vx << ", " << v._vy << ", " << v._vz << " >" ;
    return co ;
  }
  
  int operator== ( const vector3& v1, const vector3& v2 ) 
  {
    if ( ( v1._vx == v2._vx ) &&
	 ( v1._vy == v2._vy ) &&
	 ( v1._vz == v2._vz ) )
      return ( true ) ;
    else
      return ( false ) ;
  }

  int operator!= ( const vector3& v1, const vector3& v2 ) 
  {
    if ( ( v1._vx != v2._vx ) ||
	 ( v1._vy != v2._vy ) ||
	 ( v1._vz != v2._vz ) )
      return ( true ) ;
    else
      return ( false ) ;
  }

  /*! This method checks if the current vector has length() ==
    0.0.  If so, *this remains unchanged. Otherwise, *this is
    scaled by 1.0/length().

    \warning If length() is very close to zero, but not == 0.0,
    this method may behave in unexpected ways and return almost
    random results; details may depend on your particular doubleing
    point implementation. The use of this method is therefore
    highly discouraged, unless you are certain that length() is in
    a reasonable range, away from 0.0 (Stefan Kebekus)

    \deprecated This method will probably replaced by a safer
    algorithm in the future.

    \todo Replace this method with a more fool-proof version.

    @returns a reference to *this
  */
  vector3& vector3 :: normalize ()  
  {
    double l = length ();
    
    if (IsNearZero(l)) 
      return(*this);
    
    _vx = _vx / l ;
    _vy = _vy / l ;
    _vz = _vz / l ;
    
    return(*this);
  }
  
  double dot ( const vector3& v1, const vector3& v2 ) 
  {
    return v1._vx*v2._vx + v1._vy*v2._vy + v1._vz*v2._vz ;
  }
  
  vector3 cross ( const vector3& v1, const vector3& v2 ) 
  {
    vector3 vv ;
    
    vv._vx =   v1._vy*v2._vz - v1._vz*v2._vy ;
    vv._vy = - v1._vx*v2._vz + v1._vz*v2._vx ;
    vv._vz =   v1._vx*v2._vy - v1._vy*v2._vx ;
    
    return ( vv ) ;
  }


  /*! This method calculates the angle between two vectors
       
    \warning If length() of any of the two vectors is == 0.0,
    this method will divide by zero. If the product of the
    length() of the two vectors is very close to 0.0, but not ==
    0.0, this method may behave in unexpected ways and return
    almost random results; details may depend on your particular
    doubleing point implementation. The use of this method is
    therefore highly discouraged, unless you are certain that the
    length()es are in a reasonable range, away from 0.0 (Stefan
    Kebekus)

    \deprecated This method will probably replaced by a safer
    algorithm in the future.

    \todo Replace this method with a more fool-proof version.

    @returns the angle in degrees (0-360)
  */
  double vectorAngle ( const vector3& v1, const vector3& v2 ) 
  {
    double mag;
    double dp;

    mag = v1.length() * v2.length();
    dp = dot(v1,v2)/mag;

    if (dp < -0.999999)
      dp = -0.9999999;

    if (dp > 0.9999999)
      dp = 0.9999999;

    if (dp > 1.0)
      dp = 1.0;

    return((RAD_TO_DEG * acos(dp)));
  }

  double CalcTorsionAngle(const vector3 &a, const vector3 &b,
			 const vector3 &c, const vector3 &d)
  {
    double torsion;
    vector3 b1,b2,b3,c1,c2,c3;
  
    b1 = a - b;
    b2 = b - c;
    b3 = c - d;
    
    c1 = cross(b1,b2);
    c2 = cross(b2,b3);
    c3 = cross(c1,c2);
  
    if (c1.length() * c2.length() < 0.001)
      torsion = 0.0;
    else {
      torsion = vectorAngle(c1,c2);
      if (dot(b2,c3) > 0.0)
	torsion *= -1.0;
    }
    
    return(torsion);
  }
  
  /*! This method checks if the current vector *this is zero
    (i.e. if all entries == 0.0). If so, a warning message is
    printed, and the whole program is aborted with exit(0).
    Otherwise, a vector of length one is generated, which is
    orthogonal to *this, and stored in v. The resulting vector is
    not random.

    \warning If the entries of the *this (in particular the
    z-component) are very close to zero, but not == 0.0, this
    method may behave in unexpected ways and return almost random
    results; details may depend on your particular doubleing point
    implementation. The use of this method is therefore highly
    discouraged, unless you are certain that all components of
    *this are in a reasonable range, away from 0.0 (Stefan
    Kebekus)

    \deprecated This method will probably replaced by a safer
    algorithm in the future.

    \todo Replace this method with a more fool-proof version that
    does not call exit()

    @param res a reference to a vector where the result will be
    stored
  */
  void vector3::createOrthoVector(vector3 &res) const
  {
    vector3 cO;
    
    if ( ( IsNearZero(this->x())) && (IsNearZero(this->y())) ) {
      if ( IsNearZero(this->z()) ) {
	cerr << "makeorthovec zero vector" << endl;
	exit(0);
      }
      cO.SetX(1.0);
    } else {
      cO.SetZ(1.0);
    }
    res= cross(cO,*this);
    res.normalize(); 
  }

}
