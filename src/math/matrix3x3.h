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

#ifndef OB_MATRIX3x3_H
#define OB_MATRIX3x3_H


#ifdef __sgi
#include <iostream.h>
#include <fstream.h>
#else
#include <iostream>
#include <fstream>
#endif

#include <math.h>

#include "obutil.h"
#include "math/vector3.h"

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

  class matrix3x3
    {
      float ele[3][3];
    public:
      matrix3x3(void) 
	{
	  ele[0][0] = 0.0f; ele[0][1] = 0.0f; ele[0][2] = 0.0f; 
	  ele[1][0] = 0.0f; ele[1][1] = 0.0f; ele[1][2] = 0.0f; 
	  ele[2][0] = 0.0f; ele[2][1] = 0.0f; ele[2][2] = 0.0f; 
	}
      
      matrix3x3(vector3 a,vector3 b,vector3 c)
	{
	  ele[0][0] = a.x();ele[0][1] = a.y();ele[0][2] = a.z();
	  ele[1][0] = b.x();ele[1][1] = b.y();ele[1][2] = b.z();
	  ele[2][0] = c.x();ele[2][1] = c.y();ele[2][2] = c.z();
	}

      matrix3x3(float d[3][3])
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
      matrix3x3 invert();
      
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
      matrix3x3 &operator/=(const float &c);
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
      void RotAboutAxisByAngle(const vector3 &axis, const float angle);
      void RotateCoords(float *,int);
      void FillOrth(float,float,float,float,float,float);
      friend vector3 operator *(const vector3 &,const matrix3x3 &);
      friend vector3 operator *(const matrix3x3 &,const vector3 &);
      friend std::ostream& operator<< ( std::ostream&, const matrix3x3 & ) ;
    };
  
  vector3 center_coords(float*,int);
}

#endif // OB_MATRIX3x3_H
