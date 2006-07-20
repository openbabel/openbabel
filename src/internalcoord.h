/**********************************************************************
internalcoord.h - Handle OBInternalCoord class, Cartesian <=> Z-matrix

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck
 
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

#ifndef OB_INTERNALCOORD_H
#define OB_INTERNALCOORD_H

#include "babelconfig.h"

#ifndef EXTERN
#  define EXTERN extern
#endif

namespace OpenBabel
{
  
  class OBAtom;
  
  //! \brief Used to transform from z-matrix to cartesian coordinates.
  class OBAPI OBInternalCoord
  {
  public:
    //class members
    OBAtom *_a,*_b,*_c;
    double   _dst,_ang,_tor;
    //! Constructor
  OBInternalCoord(OBAtom *a= NULL,
                  OBAtom *b= NULL,
                  OBAtom *c= NULL) :
    _a(a), _b(b), _c(c), _dst(0.0), _ang(0.0), _tor(0.0)
      {}
  };
  
} // end namespace

#endif // OB_INTERNALCOORD_H

//! \file internalcoord.h
//! \brief Declaration of OBInternalCoord class, conversion between Cartesian
//!        and Z-matrix form
