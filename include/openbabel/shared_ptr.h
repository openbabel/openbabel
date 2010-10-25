/**********************************************************************
shared_ptr.h - shared_ptr class.

Copyright (C) Copyright (C) 2007 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_SHARED_PTR_H
#define OB_SHARED_PTR_H

#ifdef USE_BOOST
  #include <boost/shared_ptr.hpp>
  #define shared_ptr boost::shared_ptr
#else
  #include <memory>
  #if __GNUC__ == 4  //&& __GNUC_MINOR__ < 3  removed at the suggestion of Konstantin Tokarev
    #include <tr1/memory>
  #endif
  using std::tr1::shared_ptr;
#endif

#endif // OB_SHARED_PTR_H

/// @file shared_ptr.h
/// @brief Shared pointer.
