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

/*!
 * \brief Represents a real 3x3 matrix.
 */

  class matrix3x3
    {
      //! Elements of the matrix
      /*! This array holds the matrix. The first index refers to the
	row, the second the column. */
      float ele[3][3];
      
    public:
      //! constructs the zero-matrix
      matrix3x3(void) 
	{
	  ele[0][0] = 0.0f; ele[0][1] = 0.0f; ele[0][2] = 0.0f; 
	  ele[1][0] = 0.0f; ele[1][1] = 0.0f; ele[1][2] = 0.0f; 
	  ele[2][0] = 0.0f; ele[2][1] = 0.0f; ele[2][2] = 0.0f; 
	}

      //! constructs a matrix from row vectors
      matrix3x3(vector3 row1,vector3 row2,vector3 row3)
	{
	  ele[0][0] = row1.x();ele[0][1] = row1.y();ele[0][2] = row1.z();
	  ele[1][0] = row2.x();ele[1][1] = row2.y();ele[1][2] = row2.z();
	  ele[2][0] = row3.x();ele[2][1] = row3.y();ele[2][2] = row3.z();
	}

      //! constructs a matrix from a 3x3-array of floats
      /*! constructs a matrix from a 3x3-array of floats. The first
	index represents the row, the second index the column */
      matrix3x3(float d[3][3])
	{
	  ele[0][0] = d[0][0];ele[0][1] = d[0][1];ele[0][2] = d[0][2];
	  ele[1][0] = d[1][0];ele[1][1] = d[1][1];ele[1][2] = d[1][2];
	  ele[2][0] = d[2][0];ele[2][1] = d[2][1];ele[2][2] = d[2][2];
	}
      
      //! access function
      /*! writes the matrix into the 1-dimensional array m, row by
	row. The array must be able to hold 9 floats, otherwise your
	prgram will segfault. */ 
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

      //! generates a matrix for a random rotation
      /*! the axis of the rotation will be uniformly distributed on
	the unit sphere, the angle will be uniformly distributed in
	the interval 0..360 degrees. */
      void randomRotation(OBRandom &rnd);
      
      //! returns the determinant of the matrix
      float determinant();
      
      //! access function
      /*! \warning row or column are not in the range 0..2, random
	results are returned, and your program may even
	segfault. (Stefan Kebekus)
	
	\todo Replace this method with a more fool-proof version.
      */
      float Get(int row,int column) const {return(ele[row][column]);}
      
      //! access function
      /*! \warning row or column are not in the range 0..2, random
	variables are overwritten, and your program may
	segfault. (Stefan Kebekus)
	
	\todo Replace this method with a more fool-proof version.
      */
      void  Set(int row,int column, float v) {ele[row][column]= v;}
      
      //! divides all entries of the matrix by a scalar c
      matrix3x3 &operator/=(const float &c);

      void SetupRotMat(float,float,float);

      //! calculates a rotation matrix
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
	@param angle angle in degrees (0..360)
      */
      void RotAboutAxisByAngle(const vector3 &axis, const float angle);

      void FillOrth(float,float,float,float,float,float);

      //! matrix-vector multiplication
      friend vector3 operator *(const matrix3x3 &,const vector3 &);

      friend std::ostream& operator<< ( std::ostream&, const matrix3x3 & ) ;
    };
  
  vector3 center_coords(float*,int);
}

#endif // OB_MATRIX3x3_H
