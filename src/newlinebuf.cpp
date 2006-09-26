/**********************************************************************
newlinebuf.cpp - Filter line endings, converting \r or \r\n -> \n
 
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

#include "babelconfig.h"
#include "newlinebuf.h"

#include <streambuf>

namespace OpenBabel
{
  newlinebuf::newlinebuf(std::streambuf *sb) :
    _internalBuf(sb), _returnChar(false)
  {
    setg(_buffer, _buffer, _buffer);
  }

  newlinebuf::~newlinebuf()
  {
  }
    
  int newlinebuf::underflow()
  {
    char *ptr = _buffer;
    int ch;
    while ( (ptr != _buffer + sizeof(_buffer) - 2) 
            && ((ch = _internalBuf->sgetc()) != EOF) )
      {
        _internalBuf->sbumpc();

        switch (ch)
          {
          case 10:
            if (_returnChar) // so this is a \r\n
              _returnChar = false; // reset for normal output
            else
              *ptr++ = '\n';              
            break;
          case 13:
            _returnChar = true;
            *ptr++ = '\n';
            break;
          default:
            _returnChar = false;
            *ptr++ = ch;
          }
      }

    if (ptr == _buffer)
      return EOF;

    setg(_buffer, _buffer, ptr);
    return *_buffer;
  }

} // end namespace OpenBabel

//! \file newlinebuf.cpp
//! \brief Filter line endings -> '\n'
