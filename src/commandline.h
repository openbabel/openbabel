/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2002 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_COMMANDLINE_H
#define OB_COMMANDLINE_H

#include <stdlib.h>

#include <string>
#include <vector>
#include <map>

class ArgumentInfo
{
 public:
  int      argtyp;
  int      numargs;
  void    *argvar;
  bool     found;
  bool     list;
  bool     hasArgument;
  std::string   helpString;
  std::vector<std::string> arguments;

  ArgumentInfo() {
    numargs = 1;
    argvar  = NULL;
    argtyp  = -2;
    found   = false;
    list    = false;
    hasArgument = false;
    helpString= "";
  }

  // constructor with Help string 
  ArgumentInfo( const std::string & help) : argtyp(-2), numargs(1), argvar(NULL), 
                                found(false), list(false), hasArgument(false),
                                helpString(help) {}
};

typedef std::map<std::string, ArgumentInfo, std::less<std::string> > maptype;

class CommandLine
{
 protected:
  
  int argc, switches;
  char **argv;
  std::string arglist;

  maptype cmdline;
  void (*usage)(void);
  
  void SetVariable(ArgumentInfo &, char[]);
  bool _exitOnError;
 public:
  
  CommandLine() { switches = 0; usage = NULL; arglist = ""; _exitOnError= false; }
  CommandLine(bool exitOnError) { switches = 0; usage = NULL; arglist = ""; _exitOnError= exitOnError; }
  virtual ~CommandLine() {}

  // Switch/Flag with no Arguments
  
  void AddFlag( const char *, bool &var, bool list = false );
  void AddFlag( const char *, char *hlp, bool list = false );

  // Switch/Flag followed by Multiple Arguments

  void AddSwitch( const char *, int numargs,          bool list = false );
  void AddSwitch( const char *, std::vector<int>    &,     bool list = false );
  void AddSwitch( const char *, std::vector<float>  &,     bool list = false );
  void AddSwitch( const char *, std::vector<std::string> &,     bool list = false );
  void AddSwitch( const char *, void (*)(void *),     bool list = false );

  // Switch/Flag followed by Single Argument (with default value)

  void AddSwitch( const char *, bool   &, bool   def, bool list = false );
  void AddSwitch( const char *, int    &, int    def, bool list = false );
  void AddSwitch( const char *, float  &, float  def, bool list = false );
  void AddSwitch( const char *, char   *, char  *def, bool list = false );
  void AddSwitch( const char *, std::string &, std::string def, bool list = false );

  int ProcessArguments( char *argv[], int argc );
  int ProcessArguments( int argc, char *argv[] ) { return ProcessArguments(argv, argc); }

  bool WasCalledWith(const char *);
  bool HasArgument(const char *arg);

  char *GetArgument(const char *,int i = 0);
  char *GetArgList(void) { return (char *)arglist.c_str(); }

  void SetUsageFunction(void (*fxn)(void));
  virtual void Usage();

  // Add switched with help string
  void AddSwitch(const char *arg, std::string &var, std::string def, const char *help, bool list= false);
  void AddSwitch(const char *arg, int &var, int def, const char *help, bool list= false);
  void AddSwitch(const char *arg, float &var, float def, const char *help, bool list= false);

  void printHelp();

  bool getValue(const char *s, int &f);     
  bool getValue(const char *s, float &f);   
  bool getValue(const char *s, std::string &ret);
};

#endif // OB_COMMANDLINE_H

