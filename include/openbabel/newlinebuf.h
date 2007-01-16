/**********************************************************************
newlinebuf.h - Filter line endings, converting \r or \r\n -> \n
 
Copyright (C) 2005-2006 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
***********************************************************************/

#ifndef OB_NEWLINEBUF_H
#define OB_NEWLINEBUF_H

#include <openbabel/babelconfig.h>

#include <streambuf>
#include <istream>

#ifndef OBCONV
#define OBCONV
#endif

namespace OpenBabel
{

  //! \brief A minimal streambuf derivative to convert line endings
 class OBCONV newlinebuf : public std::streambuf
  {
  public:
    newlinebuf(std::streambuf *);
    virtual ~newlinebuf();
    
  protected:
    virtual int underflow();
  private:
    std::streambuf * const _internalBuf; //!< the internal streambuf to filter
    bool _returnChar;            //!< whether we've just seen a '\r' character
    int _chcount;
    char _buffer[8192];
    std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way, std::ios_base::openmode which);
  };

  //! \brief A convenience istream wrapper which calls an internal newlinebuf
  //!   to filter/convert line endings
 class OBCONV NewlineInput : public std::istream
    {
    private:
      newlinebuf _nBuf; //!< the internal newline filter buffer
    public:
    NewlineInput(std::istream &i) : std::istream(&_nBuf), _nBuf(i.rdbuf()) {}
    };

} // end namespace OpenBabel

#endif

//! \file newlinebuf.h
//! \brief Translate line endings automatically (UNIX, Classic Mac, DOS)
