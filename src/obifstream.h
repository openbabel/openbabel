/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2002 by Geoffrey R. Hutchison

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

#ifndef OB_IFSTREAM_H
#define OB_IFSTREAM_H

#ifdef __sgi
#include <iostream.h>
#include <fstream.h>
#else
#include <iostream>
#include <fstream>
#endif

class obifstream : public std::ifstream
{
public:
  bool SafeGetline(char *buffer,int size)
  {
    if (!getline(buffer,size)) return(false);

      char badchar1 = '\015';
      char badchar2 = '\032';

      char *p1,*p2;

      for (p1 = p2 = buffer;*p2 != '\0';)
	{
	  if (*p2 == badchar1 || *p2 == badchar2)
	    {
	      *p1 = '\0';
	      p2++; continue;
	    }
	  *p1++ = *p2++;
	}

      return(true);
    }
};

#endif //OB_IFSTREAM_H
