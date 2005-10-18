/**********************************************************************
oberror.h - Handle error messages, warnings, notices, etc.
 
Copyright (C) 2002 by Stefan Kebekus
Some portions Copyright (C) 2003-2005 by Geoffrey R. Hutchison
 
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

#ifndef OB_ERROR_H
#define OB_ERROR_H

#include "babelconfig.h"

#ifndef EXTERN
#  define EXTERN extern
#endif

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#endif
#if HAVE_SSTREAM
        #include <sstream>
#elif
        #include <sstream.h>
#endif
#include <string>
#include <vector>
#include <deque>

namespace OpenBabel
{

//! \brief Levels of error and audit messages to allow filtering
enum obMessageLevel {
   obError,     //!< for critical errors (e.g., cannot read a file)
   obWarning,   //!< for non-critical problems (e.g., molecule appears empty)
   obInfo,      //!< for informative messages (e.g., file is a non-standard format)
   obAuditMsg,  //!< for messages auditing methods which destroy or perceive molecular data (e.g., kekulization, atom typing, etc.)
   obDebug      //!< for messages only useful for debugging purposes
};

//! \brief Customizable error handling and logging -- store a message,
//!        including the method yielding the error, causes, etc.
class OBAPI OBError
{
public:

  //! Constructor for an error message e.g. OBError(__FUNCTION__, " message ")
  OBError( const std::string &method = "",
	   const std::string &errorMsg = "",
	   const std::string &explanation = "",
	   const std::string &possibleCause = "",
	   const std::string &suggestedRemedy = "",
	   const obMessageLevel = obDebug );

  //! \return a formatted message string, including optional explanations, etc.
  std::string message(void) const;
  
  //! Output a formatted message string
  friend std::ostream& operator<< ( std::ostream &os, const OBError &er )
    { return os << er.message(); };

  std::string    GetMethod()            { return _method;          }
  std::string    GetError()             { return _errorMsg;        }
  std::string    GetExplanation()       { return _explanation;     }
  std::string    GetPossibleCause()     { return _possibleCause;   }
  std::string    GetSuggestedRemedy()   { return _suggestedRemedy; }
  obMessageLevel GetLevel()             { return _level;           }

 protected:

    std::string _method;
    std::string _errorMsg;
    std::string _explanation;
    std::string _possibleCause;
    std::string _suggestedRemedy;

    obMessageLevel _level;
};

 //! \brief Handle error messages, warnings, debugging information and the like
class OBAPI OBMessageHandler
  {
  public:
    OBMessageHandler();
    ~OBMessageHandler();
    
    //! Throw an error with an already-formatted OBError object
    void ThrowError(OBError err);
    //! Throw an error in the specified method with an appropriate level
    void ThrowError(const std::string &method, const std::string &errorMsg, 
		    obMessageLevel level = obDebug);

    //! \return all messages matching a specified level
    std::vector<std::string> GetMessagesOfLevel(const obMessageLevel);

    //! Start logging messages (default)
    void StartLogging() { _logging = true; }
    //! Stop logging messages completely
    void StopLogging()  { _logging = false; }

    //! Set the maximum number of entries (or 0 for no limit)
    void SetMaxLogEntries(unsigned int max) { _maxEntries = max; }
    //! \return the current maximum number of entries (default = 0 for no limit)
    unsigned int GetMaxLogEntries() { return _maxEntries; }

    //! Clear the current message log entirely
    void ClearLog() { _messageList.clear(); }

    //! \brief Set the level of messages to output
    //! (i.e., messages with at least this priority will be output)
    void SetOutputLevel(const obMessageLevel level) { _outputLevel = level; }
    //! \return the current output level
    obMessageLevel GetOutputLevel() { return _outputLevel; }

    void SetOutputStream(std::ostream *os) { _outputStream = os; }
    std::ostream* GetOutputStream() { return _outputStream; }

    //! Start "wrapping" messages to cerr into ThrowError calls
    bool StartErrorWrap();
    //! Turn off "wrapping" messages, restoring normal cerr use (default)
    bool StopErrorWrap();

  protected:
    //! Log of messages for later retrieval via GetMessagesOfLevel()
    std::deque<OBError>    _messageList;

    //! Filtering level for messages and logging (messages of lower priority will be ignored
    obMessageLevel         _outputLevel;

    // self-explanatory
    std::ostream          *_outputStream;

    //! Whether messages will be logged into _messageList
    bool                   _logging;
    //! The maximum size of _messageList log
    unsigned int           _maxEntries;

    //! The default stream buffer for the output stream (saved if wrapping is ued)
    std::streambuf        *_inWrapStreamBuf;
    //! The filtered obLogBuf stream buffer to wrap error messages
    std::streambuf        *_filterStreamBuf;
  }; 

EXTERN OBMessageHandler obErrorLog;

//! \brief A minimal streambuf derivative to wrap calls to cerr into calls to OBMessageHandler as needed
 class OBAPI obLogBuf : public std::stringbuf
  {
    public:
      virtual ~obLogBuf() { sync(); }
    
    protected:
    //! Call OBMessageHandler::ThrowError() and flush the buffer
    int sync()
    {
      obErrorLog.ThrowError("", str(), obInfo);
      str(std::string()); // clear the buffer
			return 0;
    }
  };

} // end namespace OpenBabel

#endif

//! \file oberror.h
//! \brief Handle error messages, warnings, notices, etc. 
