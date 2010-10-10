/**********************************************************************
locale.h - Handle internal numeric locale issues -- parse data in "C"

Copyright (C) 2008 by Geoffrey R. Hutchison

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

#ifndef OB_LOCALE_H
#define OB_LOCALE_H

#include <locale>
#include <openbabel/babelconfig.h>

#ifndef OBERROR
#define OBERROR
#endif

namespace OpenBabel
{
  class OBLocalePrivate;

  // more detailed descriptions and documentation in locale.cpp
  //! \brief Handle the locale for numeric data parsing
  class OBERROR OBLocale {
  public:

    OBLocale();
    ~OBLocale();

    void SetLocale();
    void RestoreLocale();

  protected:
    OBLocalePrivate* d;
  };

  //global definitions
  //! Global OBLocale for setting and restoring locale information
  OBERROR extern  OBLocale   obLocale;

} // namespace OpenBabel
#endif // OB_LOCALE_H

//! \file locale.h
//! \brief Handle internal numeric locale issues -- parse data in "C"
