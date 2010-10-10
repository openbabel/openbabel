/**********************************************************************
tokenst.h - Tokenize and trim strings; Open data files

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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

#ifndef OB_TOKENST_H
#define OB_TOKENST_H

#include <openbabel/babelconfig.h>

#include <vector>
#include <sstream>
#include <string>
#include <fstream>

namespace OpenBabel
{
  // Utility function prototypes
  OBERROR bool tokenize(std::vector<std::string>&, const char *buf, const char *delimstr=" \t\n\r");
  OBERROR bool tokenize(std::vector<std::string>&, std::string&, const char *delimstr=" \t\n\r", int limit=-1);
  //! Remove leading and trailing whitespace from a string (docs in tokenst.cpp)
  OBERROR std::string& Trim(std::string& txt);

  // \return a string representation of a variable
  template<typename T>
  std::string toString(T val)
  {
    std::ostringstream s;
    s << val;
    return s.str();
  }

  //!Functor for deleting containers of pointers. Scott Meyers, Effective STL, Item 7
  struct DeleteObject
  {
    template<typename T>
    void operator()(const T* ptr)const { delete ptr; }
  };

  //!Read and discard all characters from input stream up to, and including, a string
  OBERROR std::istream& ignore(std::istream& ifs, const std::string& txt);

  //! Opens a datafile in a directory where OpenBabel expects to find it.
  // full documentation in tokenst.cpp
  OBERROR std::string OpenDatafile(std::ifstream& fs,
                                 const std::string& filename,
                                 const std::string& envvar = "BABEL_DATADIR");

  // Used by other code for reading files
#ifdef WIN32
  #define FILE_SEP_CHAR "\\"
#else
#define FILE_SEP_CHAR "/"
#endif


}//namespace
#endif

//! \file tokenst.h
//! \brief Tokenize strings, open data files
