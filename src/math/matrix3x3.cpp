/**********************************************************************
matrix3x3.cpp - Handle 3D rotation matrix.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2006 by Benoit Jacob

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
#include <openbabel/babelconfig.h>

#include <openbabel/math/matrix3x3.h>
#include <openbabel/obutil.h>

using namespace std;

namespace OpenBabel
{
  static inline double SQUARE( double x ) { return x*x; }
} // end namespace OpenBabel

//! \file matrix3x3.cpp
//! \brief Handle 3D rotation matrix.

