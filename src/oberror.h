/**********************************************************************
Copyright (C) 2002 by Stefan Kebekus

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


#ifndef OB_error_H
#define OB_error_H

#include <iostream>
#include <string>

class OBError 
{
 public:
  OBError( const string &method,
	   const string &errorMsg,
	   const string &explanation = "", 
	   const string &possibleCause = "", 
	   const string &suggestedRemedy ="" );

  string message(void) const;

  friend std::ostream& operator<< ( std::ostream &os, const OBError &er ){ os << er.message(); };

  string _method;
  string _errorMsg;
  string _explanation;
  string _possibleCause;
  string _suggestedRemedy;
};

#endif
