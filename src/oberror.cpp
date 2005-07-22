/**********************************************************************
oberror.cpp - Handle error messages.
 
Copyright (C) 2002 by Stefan Kebekus
Some portions Copyright (C) 2002-2005 by Geoffrey R. Hutchison
 
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

OBError::OBError( const string &method, const string &errorMsg, 
		  const string &explanation,
                  const string &possibleCause, const string &suggestedRemedy) :
  _method(method), _errorMsg(errorMsg), _explanation(explanation),
  _possibleCause(possibleCause), _suggestedRemedy(suggestedRemedy)
{

}

string OBError::message() const
{
    string tmp = "==============================\n";

    tmp += "*** Open Babel Error ";
    if (_method.length() != 0)
      {
	tmp += " in " + _method + string("\n  ");
	}
    tmp += _errorMsg + "\n";
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
  _outputLevel(obWarning), _outputStream(&clog), _logging(true)
{
  //  StartErrorWrap();
}

OBMessageHandler::~OBMessageHandler()
{ 
  StopErrorWrap();

  // free the internal filter streambuf
  if (_filterStreamBuf != NULL)
    delete _filterStreamBuf;
}
    
void OBMessageHandler::ThrowError(OBError err, obMessageLevel level)
{
  pair <OBError, obMessageLevel> p(err, level);

  _messageList.push_back(p);
  if (_maxEntries != 0 && _messageList.size() > _maxEntries)
    _messageList.pop_front();
  
  if (_logging && level <= _outputLevel)
    *_outputStream << err;
}

void OBMessageHandler::ThrowError(const std::string &method, 
				  const std::string &errorMsg,
				  obMessageLevel level)
{
  if (errorMsg.length() > 1)
    {
      OBError err(method, errorMsg);
      ThrowError(err, level);
    }
}

std::vector<std::string> OBMessageHandler::GetMessagesOfLevel(const obMessageLevel level)
{
  vector<string> results;
  deque<pair<OBError, obMessageLevel> >::iterator i;
  pair<OBError, obMessageLevel> message;

  for (i = _messageList.begin(); i != _messageList.end(); i++)
    {
      message = (*i);
      if (message.second == level)
	results.push_back( (message.first).message() );
    }

  return results;
}

bool OBMessageHandler::StartErrorWrap()
{
  if (_inWrapStreamBuf != NULL)
    return true; // already wrapped cerr  -- don't go into loops!

  _inWrapStreamBuf = cerr.rdbuf();
 
  if (_filterStreamBuf == NULL)
    {
      _filterStreamBuf = new(obLogBuf);
    }

  cerr.rdbuf(_filterStreamBuf);
  return true;
}

bool OBMessageHandler::StopErrorWrap()
{
  if (_inWrapStreamBuf == NULL)
    return true; // never wrapped cerr

  // don't free the filter streambuf yet -- we might start wrapping later
  // it's freed in the dtor

  cerr.rdbuf(_inWrapStreamBuf);
  return true;
}

} // end namespace OpenBabel

//! \file oberror.cpp
//! \brief Handle error messages, warnings, notices, etc.
//!  Implements OBMessageHandler class.
