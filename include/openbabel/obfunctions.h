/**********************************************************************
oberror.h - Handle error messages, warnings, notices, etc.

Copyright (C) 2002 by Stefan Kebekus
Some portions Copyright (C) 2003-2006 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
***********************************************************************/

#ifndef OB_FUNCTIONS_H
#define OB_FUNCTIONS_H

#include <openbabel/babelconfig.h>

#include <openbabel/atom.h>

namespace OpenBabel
{

  OBAPI void OBAtomAssignTypicalImplicitHydrogens(OBAtom* atom);

} // end namespace OpenBabel

#endif

//! \file obfunctions.h
//! \brief A collection of functions that operate on *OBAtom
