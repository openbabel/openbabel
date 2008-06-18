/**********************************************************************
locale.cpp - Handle internal numeric locale issues -- parse data in "C"
 
Copyright (C) 2008 by Geoffrey R. Hutchison
 
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

#include <openbabel/locale.h>

#if HAVE_XLOCALE_H
#include <xlocale.h>
#endif
#if HAVE_LOCALE_H
#include <locale.h>
#endif

namespace OpenBabel
{
  class OBLocalePrivate {
  public:
    char *old_locale_string;
#if HAVE_USELOCALE
    locale_t new_c_num_locale;
    locale_t old_locale;
#endif
    unsigned int counter; // Reference counter -- ensures balance in SetLocale/RestoreLocale calls

    OBLocalePrivate(): counter(0)
    {
#if HAVE_USELOCALE
      new_c_num_locale = newlocale(LC_NUMERIC_MASK, NULL, NULL);
#endif
    }
    
    ~OBLocalePrivate()
    {    }
  }; // class definition for OBLocalePrivate

  OBLocale::OBLocale()
  {
    d = new OBLocalePrivate;
  }
  
  OBLocale::~OBLocale()
  {
    if (d) {
      delete d;
      d = NULL;
    }
  }

  void OBLocale::SetLocale()
  {
    if (d->counter == 0) {
      // Set the locale for number parsing to avoid locale issues: PR#1785463
#if HAVE_USELOCALE
      // Extended per-thread interface
      d->old_locale = uselocale(d->new_c_num_locale);
#else
      // Original global POSIX interface
      d->old_locale_string = strdup (setlocale (LC_NUMERIC, NULL));
  	  setlocale(LC_NUMERIC, "C");
#endif
    }
    
    ++d->counter;
  }
  
  void OBLocale::RestoreLocale()
  {
    --d->counter;
    if(d->counter == 0) {
      // return the locale to the original one
#ifdef HAVE_USELOCALE
      uselocale(d->old_locale);
#else
      setlocale(LC_NUMERIC, d->old_locale_string);
      free (d->old_locale_string);
#endif
    }
  }

  //global definitions
  // Global OBLocale for setting and restoring locale information
  OBLocale   obLocale;

} // namespace OpenBabel

//! \file locale.cpp
//! \brief Handle internal numeric locale issues -- parse data in "C"
