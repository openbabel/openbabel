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

namespace OpenBabel {

class OBError 
{
 public:
  OBError( const std::string &method,
	   const std::string &errorMsg,
	   const std::string &explanation = "", 
	   const std::string &possibleCause = "", 
	   const std::string &suggestedRemedy ="" );

  std::string message(void) const;

  friend std::ostream& operator<< ( std::ostream &os, const OBError &er ){ return os << er.message(); };

  std::string _method;
  std::string _errorMsg;
  std::string _explanation;
  std::string _possibleCause;
  std::string _suggestedRemedy;
};

} // end namespace OpenBabel

#endif
