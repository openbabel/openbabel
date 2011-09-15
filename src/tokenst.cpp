/**********************************************************************
tokenst.cpp - Tokenize a string.

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
#include <openbabel/babelconfig.h>

#include <algorithm>
#include <vector>
#include <string>
#include <limits>
#include <cstring>
#include <openbabel/tokenst.h>

using namespace std;

namespace OpenBabel
{

  //! Break a string (supplied as the second argument) into tokens, returned
  //! in the first argument. Tokens are determined by the delimiters supplied
  //! (defaults to whitespace (i.e., spaces, tabs, newlines)
  bool tokenize(std::vector<std::string> &vcr, const char *buf,
                      const char *delimstr)
  {
    vcr.clear();
    if (!buf || !delimstr)
      return false;

    string s(buf);
    s += delimstr[0]; //forces last token to be parsed
    size_t startpos=0,endpos=0;

    for (;;)
      {
        startpos = s.find_first_not_of(delimstr,startpos);
        endpos = s.find_first_of(delimstr,startpos);

        if (endpos <= s.size() && startpos <= s.size())
          vcr.push_back(s.substr(startpos,endpos-startpos));
        else
          break;

        startpos = endpos+1;
      }

    return(true);
  }

  //! Trim any trailing spaces at the end of the supplied string.
  char *trim_spaces(char *string)
  {
    if (!string)
      return(NULL);

    size_t length;

    length = strlen(string);
    if (length == 0)
      return string;

    while ((length > 0) && (string[0] == ' '))
      {
        string++;
        --length;
      }

    if (length > 0)
      {
        while ((length > 0) && (string[length-1] == ' '))
          {
            string[length-1] = '\0';
            --length;
          }
      }

    return(string);
  }

  //! Break a string (supplied as the second argument) into tokens, returned
  //! in the first argument. Tokens are determined by the delimiters supplied
  //! (defaults to whitespace (i.e., spaces, tabs, newlines)
  //! Only breaks at most 'limit' tokens and the last item in the vector may
  //! include un-parsed tokens.
  bool tokenize(std::vector<std::string> &vcr, std::string &s,
                      const char *delimstr, int limit)
  {
    vcr.clear();
    size_t startpos=0,endpos=0;

    int matched=0;
    size_t s_size = s.size();
    for (;;)
      {
        startpos = s.find_first_not_of(delimstr,startpos);
        endpos = s.find_first_of(delimstr,startpos);
        if (endpos <= s_size && startpos <= s_size)
          {
            vcr.push_back(s.substr(startpos,endpos-startpos));

            matched++;
            if (matched == limit)
              {
                startpos = endpos+1;
                vcr.push_back(s.substr(startpos,s_size));
                break;
              }
          }
        else
          {
            if (startpos < s_size)
              vcr.push_back(s.substr(startpos,s_size-startpos));
            break;
          }

        startpos = endpos+1;
      }
    return(true);
  }

  //! Removes white space from front and back of string
  std::string& Trim(std::string& txt)
  {
    string::size_type pos = txt.find_last_not_of(" \t\n\r");
    if(pos!=string::npos)
      txt.erase(pos+1);
    else
      txt.erase();

    pos = txt.find_first_not_of(" \t\n\r");
    if(pos!=string::npos)
      txt.erase(0, pos);
    else
      txt.erase();
    return txt;
  }

  //!Read and discard all characters from input stream upto the occurence of a string
  //! \param ifs The input file stream.
  //! \param txt (which is also discarded), or the end of the stream.
  //! \returns the stream.
  std::istream& ignore(std::istream& ifs, const std::string& txt)
  {
    while(ifs)
    {
      ifs.ignore(numeric_limits<std::streamsize>::max(),txt[0]);
      size_t i=1;
      while(ifs && ifs.get()==txt[i++])
      {
        if(i==txt.size())
          return ifs;
      }
      ifs.unget();
    }
    return ifs; //at eof
  }


    /** Opens the filestream with the first file called @p filename
     found by looking
     successively in the following directories:
     - the current directory
     - in a subdirectory (of the directory below) with the version of
     OpenBabel as its name
     - the parent directory specified by the environment variable
     named @p envvar
     or "BABEL_DATADIR" if @p envvar is not specified, or the compiled-in
     macro BABEL_DATADIR if the environment variable is not set

     \param ifs        Stream to load
     \param filename   Name of the data file to load
     \param envvar     Name of the environment variable

     \return the name of the file that was opened. This includes the path
     unless it is in current directory

  **/
  std::string OpenDatafile(std::ifstream& ifs, const std::string& filename,
                           const std::string& envvar)
  {
    ios_base::openmode imode = ios_base::in;
    #ifdef ALL_READS_BINARY //Makes unix files compatible with VC++6
      imode = ios_base::in|ios_base::binary;
    #endif

    // check the current directory
    ifs.close();
    ifs.clear();
    ifs.open(filename.c_str(),imode);
    if(ifs)
      return filename;

    string file;
    const char* datadir = getenv(envvar.c_str());
    if(!datadir)
      datadir = BABEL_DATADIR;

    // check the subdirectory for this version number
    file = datadir;
    file += FILE_SEP_CHAR;
    file += BABEL_VERSION;
    file += FILE_SEP_CHAR + filename;

    ifs.clear();
    ifs.open(file.c_str(),imode);
    if(ifs)
      return file;

    // couldn't find it with the version built in, so try the parent
    file = datadir;
    file += FILE_SEP_CHAR;
    file += filename;

    ifs.clear();
    ifs.open(file.c_str(),imode);

    if (ifs)
      return file;

    ifs.clear();
    ifs.close();
    return(""); // error
  }


} // end namespace OpenBabel

//! \file tokenst.cpp
//! \brief Tokenize a string.
