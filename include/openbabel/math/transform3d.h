/**********************************************************************
transform3d.h - Handle 3D transformations in space groups.

Copyright (C) 2007 by Jean Bréfort

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_TRANSFORM_3D_H
#define OB_TRANSFORM_3D_H

#include <openbabel/math/matrix3x3.h>
#include <list>
#include <string>

namespace OpenBabel
{

  /** \class transform3d transform3d.h <openbabel/math/transform3d.h>
      \brief Handle 3D transformations, such as space group definitions
      \since version 2.2
      \sa SpaceGroup
  */
  class OBAPI transform3d: private matrix3x3, private vector3
    {
    public:
      transform3d(void): matrix3x3(), vector3()
        {
        }

      transform3d(const matrix3x3 &m, const vector3 &v): matrix3x3(m), vector3(v)
        {
          Normalize();
        }

      transform3d(double s): matrix3x3(s), vector3()
        {
        }

      //! Constructs a matrix from row vectors
      transform3d(vector3 row1,vector3 row2,vector3 row3,vector3 translation):
        matrix3x3(row1, row2, row3), vector3(translation)
        {
          Normalize();
        }

      //! \brief Constructs a matrix from a 3x3-array of doubles
      /*! The first index represents the row, the second index the column */
      transform3d(double d[3][3], double t[3]): matrix3x3(d), vector3(t)
       {
         Normalize();
       }

      vector3 operator *(const vector3 &) const;

      transform3d operator *(const transform3d &) const;

      std::string DescribeAsString() const;
      std::string DescribeAsValues() const;

      void Normalize();
	};

  typedef std::list<transform3d*>::const_iterator transform3dIterator;

}

#endif // OB_TRANSFORM_3D_H

//! \file transform3d.h
//! \brief Handle 3D transformations in space groups.
