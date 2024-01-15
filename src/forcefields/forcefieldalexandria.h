/**********************************************************************
forcefieldalexandriaf.h - Alexandria force field.

Copyright (C) 2009 by Frank Peters <e.a.j.f.peters@tue.nl>
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
Copyright (C) 2021 by David van der Spoel <david.vanderspoel@icm.uu.se>

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

#include <vector>
#include <string>
#include <map>

#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

#include "forcefieldgaff.h"

namespace OpenBabel
{

  // Class OBForceFieldAlexandria
  // class introduction in forcefieldalexandria.cpp
  class OBForceFieldAlexandria: public OBForceFieldGaff
  {
    public:
      //! Constructor
      explicit OBForceFieldAlexandria(const char* ID, bool IsDefault=true) : OBForceFieldGaff(ID, IsDefault)
      {
          SetPrmFile("alexandria.prm");
      }

  }; // class OBForceFieldAlexandria

}// namespace OpenBabel

//! \file forcefieldalexandria.h
//! \brief Alexandria force field
