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

#include <string>
#include <vector>

namespace OpenBabel
{

//! \brief Customizable error handling and logging -- store a message, including the method yielding the error.
class OBAPI OBError
{
public:

  //! Constructor for an error message e.g. OBError(__FUNCTION__, " message ")
  OBError( const std::string &method = "",
	   const std::string &errorMsg = "",
	   const std::string &explanation = "",
	   const std::string &possibleCause = "",
	   const std::string &suggestedRemedy = "" );

    std::string message(void) const;

    friend std::ostream& operator<< ( std::ostream &os, const OBError &er )
    {
        return os << er.message();
    };

 protected:

    std::string _method;
    std::string _errorMsg;
    std::string _explanation;
    std::string _possibleCause;
    std::string _suggestedRemedy;

};

//! \brief Levels of error and audit messages to allow filtering
 enum obMessageLevel { 
   obError,    //!< for critical errors (e.g., cannot read a file)
   obWarning,  //!< for non-critical problems (e.g., molecule appears empty)
   obInfo,     //!< for informative messages (e.g., file is a non-standard format)
   obAuditMsg, //!< for messages auditing methods which destroy or perceive molecular data (e.g., kekulization)
   obDebug     //!< for messages only useful for debugging purposes
 };

 //! \brief Handle error messages, warnings, debugging information and the like
class OBAPI OBMessageHandler
  {
  public:
    OBMessageHandler();
    ~OBMessageHandler();
    
    void ThrowError(OBError err, obMessageLevel level = obDebug);
    void ThrowError(const std::string &method, const std::string &errorMsg, 
		    obMessageLevel level = obDebug);
    std::vector<std::string> GetMessagesOfLevel(const obMessageLevel);

    void SetOutputLevel(const obMessageLevel level) { _outputLevel = level; }
    obMessageLevel GetOutputLevel() { return _outputLevel; }

    void SetOutputStream(std::ostream *os) { _outputStream = os; }
    std::ostream* GetOutputStream() { return _outputStream; }

  protected:
    std::vector<std::pair<OBError, obMessageLevel> >  _messageList;
    obMessageLevel         _outputLevel;
    std::ostream           *_outputStream;
  }; 

EXTERN OBMessageHandler obErrorLog;

} // end namespace OpenBabel

#endif

//! \file oberror.h
//! \brief Handle error messages, warnings, notices, etc. 
