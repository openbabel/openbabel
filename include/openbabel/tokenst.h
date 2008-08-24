/**********************************************************************
tokenst.h - Tokenize and trim strings; Open data files
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
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

#ifndef OB_TOKENST_H
#define OB_TOKENST_H

#include <openbabel/babelconfig.h>
#include <vector>
#include <string>
#include <fstream>

namespace OpenBabel
{
  // Utility function prototypes
  /** 
   * @brief Break a string (supplied as the second argument) into tokens, 
   * returned in @p vcr.  Tokens are determined by the delimiters supplied
   * (defaults to whitespace (i.e., spaces, tabs, newlines)
   * 
   * @param[out] vcr 
   * @param[in] buf The buffer to tokenize
   * @param[in] delimstr The delimiters.
   */
  OBERROR bool tokenize(std::vector<std::string>&vcr, const char *buf, const char *delimstr=" \t\n\r");
  /** 
   * @brief Break a string (supplied as the second argument) into tokens, 
   * returned in @p vcr. Tokens are determined by the delimiters supplied
   * (defaults to whitespace (i.e., spaces, tabs, newlines)
   * Only breaks at most 'limit' tokens and the last item in the vector may
   * include un-parsed tokens.
   * 
   * @param[out] vcr
   * @param[in] s The string to tokenize.
   * @param[in] delimstr The delimiters.
   * @param[in] limit Maximum number of tokens.
   * 
   * @return 
   */
  OBERROR bool tokenize(std::vector<std::string>&vcr, std::string&s, const char *delimstr=" \t\n\r", int limit=-1);
  /** 
   * @brief Trim any trailing spaces at the end of the supplied string.
   * 
   * @param txt The string to trim.
   * 
   * @return The trimed string.
   */
  OBERROR std::string& Trim(std::string& txt);

  //! Opens a datafile in a directory where OpenBabel expects to find it.
  // full documentation in tokenst.cpp
  /** 
   * Opens the filestream with the first file called @p filename
   * found by looking successively in the following directories:
   * - the current directory
   * - in a subdirectory (of the directory below) with the version of 
   * OpenBabel as its name
   * - the parent directory specified by the environment variable 
   * named @p envvar 
   * or "BABEL_DATADIR" if @p envvar is not specified, or the compiled-in
   * macro BABEL_DATADIR if the environment variable is not set
   *
   * @param ifs        Stream to load
   * @param filename   Name of the data file to load
   * @param envvar     Name of the environment variable 
   *
   * @return the name of the file that was opened. This includes the path
   *  unless it is in current directory
   */
  OBERROR std::string OpenDatafile(std::ifstream& fs, 
                                 const std::string& filename,
                                 const std::string& envvar = "BABEL_DATADIR");

  // Used by other code for reading files
  #ifdef WIN32
    #define FILE_SEP_CHAR "\\"
  #else
    #define FILE_SEP_CHAR "/"
  #endif

} //namespace

#endif

//! @file tokenst.h
//! @brief Tokenize strings, open data files
