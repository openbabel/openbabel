/**********************************************************************
oberror.cpp - Handle error messages.

Copyright (C) 2002 by Stefan Kebekus
Some portions Copyright (c) 2002-2003 by Geoffrey R. Hutchison

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

#include <iostream>
#include <string>

#include "oberror.h"

using namespace std;

namespace OpenBabel {

OBError::OBError( const string &method, const string &errorMsg, const string &explanation, 
	   const string &possibleCause, const string &suggestedRemedy)
{
  _method          = method;
  _errorMsg        = errorMsg;
  _explanation     = explanation;
  _possibleCause   = possibleCause;
  _suggestedRemedy = suggestedRemedy;
  
  cerr << message();
}

string OBError::message() const
{
  string tmp = "==============================\n";

  tmp += "OPENBABEL-ERROR in " + _method + string("\n  ") + _errorMsg + "\n";
  if (_explanation.size() != 0)
    tmp += "  " + _explanation + "\n";
  if (_possibleCause.size() != 0)
    tmp += "  Possible reason: " + _possibleCause + "\n";
  if (_suggestedRemedy.size() != 0)
    tmp += "  Suggestion: " + _suggestedRemedy + "\n";
  tmp += "==============================\n";
 return tmp;
}

} // end namespace OpenBabel
