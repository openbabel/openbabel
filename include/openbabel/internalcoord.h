/**********************************************************************
internalcoord.h - Handle OBInternalCoord class, Cartesian <=> Z-matrix

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck

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

#ifndef OB_INTERNALCOORD_H
#define OB_INTERNALCOORD_H

#include <openbabel/babelconfig.h>
#include <stddef.h>

#ifndef OB_EXTERN
#  define OB_EXTERN extern
#endif

namespace OpenBabel
{

  class OBAtom;

  /** \class OBInternalCoord internalcoord.h <openbabel/internalcoord.h>
      \brief Used to transform from z-matrix to cartesian coordinates.

      Used with OpenBabel::InternalToCartesian and
      OpenBabel::CartesianToInternal methods. Does not perform any actions
      itself. You must create or free OBAtom pointers yourself.

      The z-matrix representation uses coordinates relative to up to three
      atoms, which need not be bonded in any fashion. A rough sketch of the
      a, b, and c atoms would be:

      \code
          '*'
         /
        /
       a----b
           /
          /
         c
      \endcode

      where the OBInternalCoord record reflects the '*' atom.

      \warning Does not detect if NULL pointers are used. You should be careful.
   **/
  class OBAPI OBInternalCoord
  {
  public:
    //class members
    OBAtom *_a;   //!< First connection for this atom (i.e., distance)
    OBAtom *_b;   //!< Second reference atom (i.e., angle)
    OBAtom *_c;   //!< Third reference atom (i.e., dihedral / torsion angle)
    double  _dst; //!< Distance between this atom and _a
    double  _ang; //!< Angle between this, _a, and _b (i.e., _a is the vertex)
    double  _tor; //!< Torsional/dihedral angle between this, _a, _b, and _c

    //! Constructor
  OBInternalCoord(OBAtom *a= nullptr, OBAtom *b= nullptr, OBAtom *c= nullptr,
                  double dst = 0.0, double ang = 0.0, double tor = 0.0) :
    _a(a), _b(b), _c(c), _dst(dst), _ang(ang), _tor(tor)
      {}
  };

} // end namespace

#endif // OB_INTERNALCOORD_H

//! \file internalcoord.h
//! \brief Declaration of OBInternalCoord class, conversion between Cartesian
//!        and Z-matrix form
