/**********************************************************************
op.cpp -

Copyright (C) 2009 Tim Vandermeersch

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

#include <openbabel/babelconfig.h>
#include <openbabel/op.h>

using namespace std;
namespace OpenBabel
{
#if defined(__CYGWIN__) || defined(__MINGW32__)
  // macro to implement static OBPlugin::PluginMapType& Map()
  PLUGIN_CPP_FILE(OBOp)
#endif
}


//! \file op.cpp
//! \brief Base plugin class for operations on molecules
