/**********************************************************************
oberror.cpp - Handle error messages.

Copyright (C) 2002 by Stefan Kebekus
Some portions Copyright (C) 2002-2006 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>

#include <iostream>
#include <string>
#include <algorithm>

#include <openbabel/oberror.h>

using namespace std;

namespace OpenBabel
{
  // Initialize the global obErrorLog declared in oberror.h
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
    if (!_explanation.empty())
      tmp += "  " + _explanation + "\n";
    if (!_possibleCause.empty())
      tmp += "  Possible reason: " + _possibleCause + "\n";
    if (!_suggestedRemedy.empty())
      tmp += "  Suggestion: " + _suggestedRemedy + "\n";
    return tmp;
  }

  bool OBError::operator== (const OBError& other)const {return GetError()==other.GetError();}

  /** \class OBMessageHandler oberror.h <openbabel/oberror.h>

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
  OBERROR extern OBMessageHandler obErrorLog;
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
    _messageCount[0] = _messageCount[1] = _messageCount[2] = 0;
    _messageCount[3] = _messageCount[4] = 0;
    _filterStreamBuf = _inWrapStreamBuf = NULL;
    //  StartErrorWrap(); // (don't turn on error wrapping by default)
  }

  OBMessageHandler::~OBMessageHandler()
  {
    StopErrorWrap();

    // free the internal filter streambuf
    if (_filterStreamBuf != NULL)
      delete _filterStreamBuf;
  }

  void OBMessageHandler::ThrowError(OBError err, errorQualifier qualifier)
  {
    if (!_logging)
      return;

    //Output error message if level sufficiently high and, if onceOnly set, it has not been logged before
    if (err.GetLevel() <= _outputLevel &&
      (qualifier!=onceOnly || find(_messageList.begin(), _messageList.end(), err)==_messageList.end()))
    {
      *_outputStream << err;
    }

    _messageList.push_back(err);
    _messageCount[err.GetLevel()]++;
    if (_maxEntries != 0 && _messageList.size() > _maxEntries)
      _messageList.pop_front();
  }

  void OBMessageHandler::ThrowError(const std::string &method,
                                    const std::string &errorMsg,
                                    obMessageLevel level, errorQualifier qualifier)
  {
    if (errorMsg.length() > 1)
      {
        OBError err(method, errorMsg, "", "", "", level);
        ThrowError(err, qualifier);
      }
  }

  std::vector<std::string> OBMessageHandler::GetMessagesOfLevel(const obMessageLevel level)
  {
    vector<string> results;
    deque<OBError>::iterator i;
    OBError error;

    for (i = _messageList.begin(); i != _messageList.end(); ++i)
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

    cerr.rdbuf(_inWrapStreamBuf);
    _inWrapStreamBuf=NULL;//shows not wrapped

    // don't delete the filter streambuf yet -- we might start wrapping later
    // it's freed in the dtor

    return true;
  }

  string OBMessageHandler::GetMessageSummary()
  {
    stringstream summary;
    if (_messageCount[obError] > 0)
      summary << _messageCount[obError] << " errors ";
    if (_messageCount[obWarning] > 0)
      summary << _messageCount[obWarning] << " warnings ";
    if (_messageCount[obInfo] > 0)
      summary << _messageCount[obInfo] << " info messages ";
    if (_messageCount[obAuditMsg] > 0)
      summary << _messageCount[obAuditMsg] << " audit log messages ";
    if (_messageCount[obDebug] > 0)
      summary << _messageCount[obDebug] << " debugging messages ";

    return summary.str();
  }

} // end namespace OpenBabel

//! \file oberror.cpp
//! \brief Handle error messages, warnings, notices, etc.
//!  Implements OBMessageHandler class.
