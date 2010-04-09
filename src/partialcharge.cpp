/**********************************************************************
partialcharge.cpp - Base class for partial charge models
 
Copyright (C) 2010 Geoffrey Hutchison
 
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
#include <openbabel/partialcharge.h>

using namespace std;
namespace OpenBabel
{
#if defined(__CYGWIN__) || defined(__MINGW32__)
  // macro to implement static OBPlugin::PluginMapType& Map()
  PLUGIN_CPP_FILE(OBPartialCharge)
#endif
}


//! \file partialcharge.cpp
//! \brief Base class for molecular partial charge models
