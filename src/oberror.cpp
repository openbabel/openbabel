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

OBError::OBError( const string &method, 
		  const string &errorMsg, 
		  const string &explanation,
                  const string &possibleCause,
		  const string &suggestedRemedy,
		  const obMessageLevel level) :
  _method(method), _errorMsg(errorMsg), _explanation(explanation),
  _possibleCause(possibleCause), _suggestedRemedy(suggestedRemedy),
  _level(level)
{ }

string OBError::message() const
{
    string tmp = "==============================\n";

    if (_level == obError)
      tmp += "*** Open Babel Error ";
    else if (_level == obWarning)
      tmp += "*** Open Babel Warning ";
    else if (_level == obInfo)
      tmp += "*** Open Babel Information ";
    else if (_level == obAuditMsg)
      tmp += "*** Open Babel Audit Log ";
    else
      tmp += "*** Open Babel Debugging Message ";

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
    return tmp;
}


/** \class OBMessageHandler

 OBMessageHandler represents a configurable error system for Open Babel.

A global error log is defined by the Open Babel library for use of
internal code as well as code built on top of Open Babel. This class
allows flexible filtering based on urgency (defined by the
obMessageLevel type), an "audit log" of molecular changes, including
recall using the GetMessagesOfLevel method, etc.

The default is to only log and output errors of priority
obMessageLevel::obError or obMessageLevel::obWarning.

Long-running code may wish to set the size of the in-memory error log
using the StartLogging / StopLogging methods and
SetMaxLogEntries. Otherwise, the error log may easily fill up,
requiring large amounts of memory.

If you wish to divert error output to a different std::ostream (i.e.,
for graphical display, or a file log), use the SetOutputStream method
-- the default goes to the std::clog stream. Furthermore, some older
code uses std::cerr for direct error output, rather than the
ThrowError() methods in this class. To prevent this, you can turn on
"error wrapping" using the StartErrorWrap method -- this behavior is
turned off by default.

To make it easy to use the OBMessageHandler class and error logging
facilities, a global log is defined:

\code
EXTERN OBMessageHandler obErrorLog;
\endcode

Therefore, it is very easy to log errors:

\code
   if (atomIndex < 1 || atomIndex > mol.NumAtoms() )
    obErrorLog.ThrowError(__FUNCTION__, "Requested Atom Out of Range", obDebug);
\endcode

or

\code
     stringstream errorMsg;
     errorMsg << " Could not parse line in type translation table types.txt -- incorect number of columns";
     errorMsg << " found " << vc.size() << " expected " << _ncols << ".";
     obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
\endcode

The __FUNCTION__ builtin is defined by many compilers (e.g., <a
href="http://gcc.gnu.org/">GCC</a>) but can be defined to an empty
string on some platforms without this compiler extension.

Output from the error log typically looks like:
\code
==============================
*** Open Babel Audit Log  in ReadChemObject
  OpenBabel::Read molecule Protein Data Bank format
==============================
*** Open Babel Information  in ParseConectRecord
  WARNING: Problems reading a PDB file
  Problems reading a CONECT record.
  According to the PDB specification,
  the record should have 70 columns, but OpenBabel found 61 columns.
\endcode

**/

OBMessageHandler::OBMessageHandler() :
  _outputLevel(obWarning), _outputStream(&clog), _logging(true), _maxEntries(100)
{
  //  StartErrorWrap(); // (don't turn on error wrapping by default)
}

OBMessageHandler::~OBMessageHandler()
{ 
  StopErrorWrap();

  // free the internal filter streambuf
  if (_filterStreamBuf != NULL)
    delete _filterStreamBuf;
}
    
void OBMessageHandler::ThrowError(OBError err)
{
  _messageList.push_back(err);
  if (_maxEntries != 0 && _messageList.size() > _maxEntries)
    _messageList.pop_front();
  
  if (_logging && err.GetLevel() <= _outputLevel)
    {
      *_outputStream << err;
    }
}

void OBMessageHandler::ThrowError(const std::string &method, 
				  const std::string &errorMsg,
				  obMessageLevel level)
{
  if (errorMsg.length() > 1)
    {
      OBError err(method, errorMsg, "", "", "", level);
      ThrowError(err);
    }
}

std::vector<std::string> OBMessageHandler::GetMessagesOfLevel(const obMessageLevel level)
{
  vector<string> results;
  deque<OBError>::iterator i;
  OBError error;

  for (i = _messageList.begin(); i != _messageList.end(); i++)
    {
      error = (*i);
      if (error.GetLevel() == level)
	results.push_back( error.message() );
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
