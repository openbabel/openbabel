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
#include "math/vector3.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

using namespace std;

namespace OpenBabel {

  float& vector3::operator[] ( unsigned int i)
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
    float l;
    do {
      this->Set(ptr->NextFloat()-0.5f, ptr->NextFloat()-0.5f, ptr->NextFloat()-0.5f);
      l = length_2();
    } while ( (l > 1.0) || (l < 0.0001) );
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

  vector3& vector3 :: normalize ()  
  {
    float l = length ();
    
    if (l == 0) 
      return(*this);
    
    _vx = _vx / l ;
    _vy = _vy / l ;
    _vz = _vz / l ;
    
    return(*this);
  }
  
  float dot ( const vector3& v1, const vector3& v2 ) 
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


  // ***angle***

  float vectorAngle ( const vector3& v1, const vector3& v2 ) 
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

  float CalcTorsionAngle(const vector3 &a, const vector3 &b,
			 const vector3 &c, const vector3 &d)
  {
    float torsion;
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
  
  void vector3::createOrthoVector(vector3 &res) const
  {
    vector3 cO;
    
    if ((this->x() == 0.0)&&(this->y() == 0.0)) {
      if (this->z() == 0.0){
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
