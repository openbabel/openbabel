/**********************************************************************
oberror.h - Handle error messages, warnings, notices, etc.

Copyright (C) 2002 by Stefan Kebekus
Some portions Copyright (C) 2003-2006 by Geoffrey R. Hutchison

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

#ifndef OB_ERROR_H
#define OB_ERROR_H

#include <openbabel/babelconfig.h>

#include <iosfwd>
#include <sstream>
#include <string>
#include <vector>
#include <deque>

#ifndef OBERROR
#define OBERROR
#endif

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

  enum errorQualifier {always, onceOnly};

  /** \class OBError oberror.h <openbabel/oberror.h>
      \brief Customizable error handling and logging -- store a message,
      including the method yielding the error, causes, etc. **/
  class OBERROR OBError
    {
    public:

      //! Constructor for an error message e.g. OBError(__FUNCTION__, " message ")
      OBError( const std::string &method = "",
               const std::string &errorMsg = "",
               const std::string &explanation = "",
               const std::string &possibleCause = "",
               const std::string &suggestedRemedy = "",
               const obMessageLevel = obDebug );

      //! \return A formatted message string, including optional explanations, etc.
      std::string message(void) const;

      //! Output a formatted message string
      friend std::ostream& operator<< ( std::ostream &os, const OBError &er )
        { return os << er.message(); };

      /** \return The method which caused this error
          (typically supplied via the __FUNCTION__ compiler macro **/
      std::string    GetMethod() const           { return _method;          }
      //! \return The main error message
      std::string    GetError() const            { return _errorMsg;        }
      //! \return A more detailed explanation of the error (optional)
      std::string    GetExplanation() const      { return _explanation;     }
      //! \return A possible cause for the error (optional)
      std::string    GetPossibleCause() const    { return _possibleCause;   }
      //! \return The suggested workaround or remedy for the error (optional)
      std::string    GetSuggestedRemedy() const  { return _suggestedRemedy; }
      //! \return The severity level of this error
      obMessageLevel GetLevel() const            { return _level;           }

      bool operator== (const OBError&)const;

    protected:

      //! The method causing the error (typically from the compiler macro __FUNCTION__)
      std::string _method;
      //! The actual error message
      std::string _errorMsg;
      //! Optional explanation message: more detailed than the brief error
      std::string _explanation;
      //! Optional cause message
      std::string _possibleCause;
      //! Optional workaround or remedy
      std::string _suggestedRemedy;

      //! The severity level: used for filtering via OBMessageHandler
      obMessageLevel _level;
    };

  //! \brief Handle error messages, warnings, debugging information and the like
  // More documentation in oberror.cpp
  class OBERROR OBMessageHandler
    {
    protected:
      //! Count of messages at each message level
      unsigned int           _messageCount[5];

    public:
      OBMessageHandler();
      ~OBMessageHandler();

      //! Throw an error with an already-formatted OBError object
      void ThrowError(OBError err, errorQualifier qqualifier = always);
      //! Throw an error in the specified method with an appropriate level
      void ThrowError(const std::string &method, const std::string &errorMsg,
                      obMessageLevel level = obDebug, errorQualifier qualifier = always);

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

      //! \return Count of messages received at the obError level
      unsigned int GetErrorMessageCount() { return _messageCount[obError];}
      //! \return Count of messages received at the obWarning level
      unsigned int GetWarningMessageCount() { return _messageCount[obWarning];}
      //! \return Count of messages received at the obInfo level
      unsigned int GetInfoMessageCount() { return _messageCount[obInfo];}
      //! \return Count of messages received at the obAuditMsg level
      unsigned int GetAuditMessageCount() { return _messageCount[obAuditMsg];}
      //! \return Count of messages received at the obDebug level
      unsigned int GetDebugMessageCount() { return _messageCount[obDebug];}
      //! \return Summary of messages received at all levels
      std::string GetMessageSummary();

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

  //! Global OBMessageHandler error handler
  OBERROR extern  OBMessageHandler obErrorLog;

  //! \class obLogBuf oberror.h <openbabel/oberror.h>
  //! \brief A minimal streambuf derivative to wrap calls to cerr into calls to OBMessageHandler as needed
  /** This class is used for internal use, via OBMessageHandler::StartErrorWrap()
      To use this class, use the global OBMessageHandler object @p obErrorLog:
      \code
      obErrorLog.StartErrorWrap(); // All output to cerr will become OBErrors
      cerr << " This is error 1" << endl; // flush output, create a new error
      cerr << " Error 2" << endl;
      cerr << " Error 3: Done with output wrapping." << endl;
      obErrorLog.StopErrorWrap(); // return to default behavior
      \endcode
  **/
  class OBERROR obLogBuf : public std::stringbuf
    {
    public:
      //! Close the output buffer, flush, and call OBMessageHandler::ThrowError()
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
