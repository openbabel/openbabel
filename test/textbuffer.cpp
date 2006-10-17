/**********************************************************************
textbuffer.cpp - Unit tests for Open Babel buffer/stream classes

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include "babelconfig.h"
#include "zipstream.h"
#include "newlinebuf.h"

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenBabel;

#ifndef BUFF_SIZE
#define BUFF_SIZE 8095
#endif

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      cout << "Usage: textbuffer" << endl;
      cout << " Unit tests for buffer/stream reading/writing " << endl;
      return(-1);
    }

  cout << "# Unit tests for text buffer classes \n";

  istream* pInStream;
  newlinebuf filter(cin.rdbuf());
  cin.rdbuf(&filter);
  
  zlib_stream::zip_istream zIn(cin);
  if (zIn.is_gzip())
    pInStream = &zIn;
  else
    pInStream = &cin;

  //  char buffer[BUFF_SIZE];
  string buffer;
  unsigned int linecount = 0;
  while(getline(*pInStream, buffer))
    {
      ++linecount;
    }
  cout << linecount << endl;

  return(0);
}
