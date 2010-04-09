/**********************************************************************
partialcharge.h - Base class for partial charge models
 
Copyright (C) 2010 by Geoffrey Hutchison
 
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

#ifndef OB_PARTIALCHARGE_H
#define OB_PARTIALCHARGE_H

#include <string>
#include <sstream>
#include <limits>

#include <openbabel/babelconfig.h>
#include <openbabel/plugin.h>

namespace OpenBabel
{
class OBBase; //Forward declaration; used only as pointer.
class OBMol;

class OBAPI OBPartialCharge : public OBPlugin
{
  MAKE_PLUGIN(OBPartialCharge)

  public:
    virtual const char* TypeID(){return "charges";};

    /// \return whether partial charges were successfully assigned to this molecule
    virtual bool AssignPartialCharges(OBMol &mol) { return false; }

protected:
};

}//namespace
#endif

//! \file partialcharge.h
//! \brief Base class for molecular partial charge models
