/**********************************************************************
oberror.cpp - Handle error messages.
 
Copyright (C) 2002 by Stefan Kebekus
Some portions Copyright (c) 2002-2005 by Geoffrey R. Hutchison
 
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

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#endif

#include <string>

#include "oberror.h"

using namespace std;

namespace OpenBabel
{

  OBMessageHandler obErrorLog;

OBError::OBError( const string &method, const string &errorMsg, const string &explanation,
                  const string &possibleCause, const string &suggestedRemedy) :
  _method(method), _errorMsg(errorMsg), _explanation(explanation),
  _possibleCause(possibleCause), _suggestedRemedy(suggestedRemedy)
{

}

string OBError::message() const
{
    string tmp = "==============================\n";

    tmp += "*** Open Babel Error in " + _method + string("\n  ") + _errorMsg + "\n";
    if (_explanation.size() != 0)
        tmp += "  " + _explanation + "\n";
    if (_possibleCause.size() != 0)
        tmp += "  Possible reason: " + _possibleCause + "\n";
    if (_suggestedRemedy.size() != 0)
        tmp += "  Suggestion: " + _suggestedRemedy + "\n";
    tmp += "==============================\n";
    return tmp;
}


OBMessageHandler::OBMessageHandler() :
  _outputStream(&cerr), _outputLevel(obWarning)
{ }

OBMessageHandler::~OBMessageHandler()
{ }
    
void OBMessageHandler::ThrowError(OBError err, obMessageLevel level)
{
  pair <OBError, obMessageLevel> p(err, level);
  _messageList.push_back(p);

  if (level <= _outputLevel)
    *_outputStream << err;
}

void OBMessageHandler::ThrowError(const std::string &method, 
				  const std::string &errorMsg,
				  obMessageLevel level)
{
  OBError err(method, errorMsg);

  ThrowError(err, level);
}

std::vector<std::string> OBMessageHandler::GetMessagesOfLevel(const obMessageLevel level)
{
  vector<string> results;
  vector<pair<OBError, obMessageLevel> >::iterator i;
  pair<OBError, obMessageLevel> message;

  for (i = _messageList.begin(); i != _messageList.end(); i++)
    {
      message = (*i);
      if (message.second == level)
	results.push_back( (message.first).message() );
    }

  return results;
}


} // end namespace OpenBabel
