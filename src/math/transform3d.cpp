/**********************************************************************
transform3d.cpp - Handle 3D transformations in space groups.

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
#include <openbabel/babelconfig.h>

#include <openbabel/math/transform3d.h>

using namespace std;

namespace OpenBabel
{

  /*!
  */
  vector3 transform3d::operator *(const vector3 &v) const
    {
      return *static_cast <const matrix3x3 *> (this) * v + *static_cast <const vector3 *> (this);
    }

  /*!
  */
  transform3d transform3d::operator *(const transform3d &t) const
    {
      return transform3d(
                    *static_cast <const matrix3x3 *> (this) * *static_cast <const matrix3x3 *> (&t),
                    *this * *static_cast <const vector3 *> (&t)
                  );
    }

  /*!
  */
  string transform3d::DescribeAsString() const
    {
      ostringstream r;
      int n, i, j;
      const matrix3x3 *m = static_cast <const matrix3x3 *> (this);
      const vector3 *v = static_cast <const vector3 *> (this);
      bool neg, first;
      for (i = 0; i < 3; i++)
        {
          if (i)
            r << ",";
          n = static_cast<int> (floor((*v)[i] * 12.0 + 0.1));
          j = 0;
          while ((*m)(i, j) == 0.)
            j++;
          neg = (*m)(i, j) < 0.;
          switch (n)
            {
              case 2:
                r << ((neg)? "1/6": "1/6+");
                break;
              case 3:
                r << ((neg)? "1/4": "1/4+");
                break;
              case 4:
                r << ((neg)? "1/3": "1/3+");
                break;
              case 6:
                r << ((neg)? "1/2": "1/2+");
                break;
              case 8:
                r << ((neg)? "2/3": "2/3+");
                break;
              case 9:
                r << ((neg)? "3/4": "3/4+");
                break;
              case 10:
                r << ((neg)? "5/6": "5/6+");
                break;
            }
            first = true;
            while (j < 3) {
              if ((*m)(i, j) != 0.)
                {
                  neg = (*m)(i, j) < 0.;
                  switch (j)
                    {
                      case 0:
                        r << ((neg)? "-x": (first? "x": "+x"));
                        break;
                     case 1:
                        r << ((neg)? "-y": (first? "y": "+y"));
                        break;
                     case 2:
                        r << ((neg)? "-z": (first? "z": "+z"));
                        break;
                    }
                  first = false;
                }
              j++;
            }
        }
      return r.str();
    }

  /*!
  */
  string transform3d::DescribeAsValues() const
    {
      ostringstream oss;
      const matrix3x3 *m = static_cast <const matrix3x3 *> (this);
      const vector3 *v = static_cast <const vector3 *> (this);
      oss << (*m)(0,0) << " " << (*m)(0,1) << " " << (*m)(0,2) << " " << v->x() << " ";
      oss << (*m)(1,0) << " " << (*m)(1,1) << " " << (*m)(1,2) << " " << v->y() << " ";
      oss << (*m)(2,0) << " " << (*m)(2,1) << " " << (*m)(2,2) << " " << v->z();
      return oss.str();
    }

  /*!
  */
  void transform3d::Normalize()
	{
      vector3 *vv = static_cast<vector3*>(this);
      vv->x() -= floor (vv->x() + .01); /* .01 should work in all cases in
                                        this context */
      vv->y() -= floor (vv->y() + .01);
      vv->z() -= floor (vv->z() + .01);
	}

}

//! \file transform3d.cpp
//! \brief Handle 3D transformations in space groups.
