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

#include "oberror.h"

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

      //! constructs s times the unit matrix
      matrix3x3(float s) 
	{
	  ele[0][0] = s;    ele[0][1] = 0.0f; ele[0][2] = 0.0f; 
	  ele[1][0] = 0.0f; ele[1][1] = s;    ele[1][2] = 0.0f; 
	  ele[2][0] = 0.0f; ele[2][1] = 0.0f; ele[2][2] = s; 
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
      
      //! Calculates the inverse of a matrix.
      /*! This method checks if the absolute value of the determinant is smaller than 1e-6. If
	so, nothing is done and an exception is thrown. Otherwise, the
	inverse matrix is calculated and returned. *this is not changed.
	
	\warning If the determinant is close to zero, but not == 0.0f,
	this method may behave in unexpected ways and return almost
	random results; details may depend on your particular floating
	point implementation. The use of this method is therefore highly
	discouraged, unless you are certain that the determinant is in a
	reasonable range, away from 0.0f (Stefan Kebekus)
      */
      matrix3x3 inverse(void) const throw(OBError);

      //! Calculates the transpose of a matrix.
      /* This method returns the transpose of a matrix. The original
	 matrix remains unchanged. */
      matrix3x3 transpose(void) const;

      //! generates a matrix for a random rotation
      /*! the axis of the rotation will be uniformly distributed on
	the unit sphere, the angle will be uniformly distributed in
	the interval 0..360 degrees. */
      void randomRotation(OBRandom &rnd);
      
      //! returns the determinant of the matrix
      float determinant() const;
      
      //! Checks if a matrix is symmetric
      /*! This method returns false if there are indices i,j such that
	fabs(*this[i][j]-*this[j][i]) > 1e-6. Otherwise, it returns
	true. */
      bool isSymmetric(void) const;

      //! Checks if a matrix is orthogonal
      /*! This method checks if a matrix describes an orthogonal
	transformation, i.e. if all column vectors are normalized and
	are mutually orthogonal. An orthogonal transformation is a
	transformation the preserves length and angle. 

	The check is performed using the method isUnitMatrix() to
	check if
	\code
	*this * transpose()
	\endcode
	is a unit matrix. The criterion is therefore numerically quite
	tight. */
      bool isOrthogonal(void) const {return (*this * transpose()).isUnitMatrix();};

      //! Checks if a matrix is diagonal
      /*! This method returns false if there are indices i != j such
	that fabs(*this[i][j]) > 1e-6. Otherwise, it returns true. */
      bool isDiagonal(void) const;

      //! Checks if a matrix is the unit matrix
      /*! This method returns false if there are indices i != j such
	that fabs(*this[i][j]) > 1e-6, or if there is an index i such
	that fabs(*this[i][j]-1) > 1e-6. Otherwise, it returns
	true. */
      bool isUnitMatrix(void) const;

      //! access function
      /*! \warning row or column are not in the range 0..2, random
	results are returned, and your program may even
	segfault. (Stefan Kebekus)
	
	\todo Replace this method with a more fool-proof version.
      */
      float Get(int row,int column) const {return(ele[row][column]);}
      
      //! access function
      /*! \warning if row or column are not in the range 0..2, random
	variables are overwritten, and your program may
	segfault. (Stefan Kebekus)
	
	\todo Replace this method with a more fool-proof version.
      */
      void Set(int row,int column, float v) {ele[row][column]= v;}

      //! access function
      /*! \warning If column is not in the range 0..2, the vector
	remains unchanged and an exception is thrown. */
      void SetColumn(int column, const vector3 &v) throw(OBError);

      //! access function
      /*! \warning If col is not in the range 0..2, an exception is
	thrown. */
      vector3 GetColumn(unsigned int col) const throw(OBError);

      //! access function
      /*! \warning If row is not in the range 0..2, an exception is
	thrown. */
      vector3 GetRow(unsigned int row) const throw(OBError);

      
      //! divides all entries of the matrix by a scalar c
      matrix3x3 &operator/=(const float &c);

      void SetupRotMat(float,float,float);

      //! calculates a matrix that represents reflection on a plane
      /*! Replaces *this with a matrix that represents reflection on
	the plane through 0 which is given by the normal vector norm.
	
	\warning If the vector norm has length zero, this method will
	generate the 0-matrix. If the length of the axis is close to
	zero, but not == 0.0f, this method may behave in unexpected
	ways and return almost random results; details may depend on
	your particular floating point implementation. The use of this
	method is therefore highly discouraged, unless you are certain
	that the length is in a reasonable range, away from 0.0f
	(Stefan Kebekus)
	
	\deprecated This method will probably replaced by a safer
	algorithm in the future.
	
	\todo Replace this method with a more fool-proof version.
	
	@param norm specifies the normal to the plane
      */
      void PlaneReflection(const vector3 &norm);

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

      //! find the eigenvalues and -vectors of a symmetric matrix
      /*! This method employs the static method matrix3x3::jacobi(...)
	to find the eigenvalues and eigenvectors of a symmetric
	matrix. On entry it is checked if the matrix really is
	symmetric: if isSymmetric() returns 'false', an OBError is
	thrown.

	\note The jacobi algorithm is should work great for all
	symmetric 3x3 matrices. If you need to find the eigenvectors
	of a non-symmetric matrix, you might want to resort to the
	sophisticated routines of LAPACK.

	@param eigenvals a reference to a vector3 where the
	eigenvalues will be stored. The eigenvalues are ordered so
	that eigenvals[0] <= eigenvals[1] <= eigenvals[2].

	@return an orthogonal matrix whose ith column is an
	eigenvector for the eigenvalue eigenvals[i]. Here 'orthogonal'
	means that all eigenvectors have length one and are mutually
	orthogonal. The ith eigenvector can thus be conveniently
	accessed by the GetColumn() method, as in the following
	example.
	\code
	// Calculate eigenvectors and -values
	vector3 eigenvals;
	matrix3x3 eigenmatrix = somematrix.findEigenvectorsIfSymmetric(eigenvals);
	
	// Print the 2nd eigenvector
	cout << eigenmatrix.GetColumn(1) << endl;
	\endcode
	With these conventions, a matrix is diagonalized in the following way:
	\code
	// Diagonalize the matrix
	matrix3x3 diagonalMatrix = eigenmatrix.inverse() * somematrix * eigenmatrix;
	\endcode

       */
      matrix3x3 findEigenvectorsIfSymmetric(vector3 &eigenvals) const throw(OBError);

      //! matrix-vector multiplication
      friend vector3 operator *(const matrix3x3 &,const vector3 &);

      //! matrix-matrix multiplication
      friend matrix3x3 operator *(const matrix3x3 &,const matrix3x3 &);

      friend std::ostream& operator<< ( std::ostream&, const matrix3x3 & ) ;

      //! eigenvalue calculation
      /*! This static function computes the eigenvalues and
      eigenvectors of a SYMMETRIC nxn matrix. This method is used
      internally by OpenBabel, but may be useful as a general
      eigenvalue finder.

      The algorithm uses Jacobi transformations. It is described
      e.g. in Wilkinson, Reinsch "Handbook for automatic computation,
      Volume II: Linear Algebra", part II, contribution II/1. The
      implementation is also similar to the implementation in this
      book. This method is adequate to solve the eigenproblem for
      small matrices, of size perhaps up to 10x10. For bigger
      problems, you might want to resort to the sophisticated routines
      of LAPACK.

      \note If you plan to find the eigenvalues of a symmetric 3x3
      matrix, you will probably prefer to use the more convenient
      method findEigenvectorsIfSymmetric()

      @param n the size of the matrix that should be diagonalized

      @param a array of size n^2 which holds the symmetric matrix
      whose eigenvectors are to be computed. The convention is that
      the entry in row r and column c is addressed as a[n*r+c] where,
      of course, 0 <= r < n and 0 <= c < n. There is no check that the
      matrix is actually symmetric. If it is not, the behaviour of
      this function is undefined.  On return, the matrix is
      overwritten with junk.

      @param d pointer to a field of at least n floats which will be
      overwritten. On return of this function, the entries d[0]..d[n-1]
      will contain the eigenvalues of the matrix.

      @param v an array of size n^2 where the eigenvectors will be
      stored. On return, the columns of this matrix will contain the
      eigenvectors. The eigenvectors are normalized and mutually
      orthogonal.
      */
      static void jacobi(unsigned int n, float *a, float *d, float *v);
    };

  vector3 center_coords(float*,int);
}

#endif // OB_MATRIX3x3_H
