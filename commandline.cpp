/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include <stdio.h>
#include <stdlib.h>

#include "mol.h"
#include "commandline.h"

using namespace OpenBabel;
using namespace std;

enum type { t_bool = -1, t_int, t_float, t_char, t_string, t_vector_int, t_vector_float, 
	    t_vector_string, t_func };

void CommandLine::AddFlag(const char *arg, bool &var,bool list)
{
  ArgumentInfo args;
  args.numargs = 0;
  args.argvar = (void *) &var;
  args.argtyp = t_bool;
  var = false;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}

void CommandLine::AddFlag( const char *arg, char *hlp, bool list)
{
  ArgumentInfo args(hlp);
  args.list    = list;
  args.numargs = 0;
  cmdline.insert(make_pair(string(arg), args));
  switches++;

}

void CommandLine::AddSwitch(const char *arg, int numargs, bool list )
{
  ArgumentInfo args;
  args.numargs = numargs;
  args.list    = list;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}


void CommandLine::AddSwitch(const char *arg, bool &var, bool def, bool list)
{
  ArgumentInfo args;
  args.argvar = (void *) &var;
  args.argtyp = t_bool;
  args.list   = list;
  var = def;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}


void CommandLine::AddSwitch(const char *arg, int &var, int def, bool list)
{
  ArgumentInfo args;
  args.argvar = (void *) &var;
  args.argtyp = t_int;
  args.list   = list;
  var = def;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}


void CommandLine::AddSwitch(const char *arg, float &var, float def, bool list)
{
  ArgumentInfo args;
  args.argvar = (void *) &var;
  args.argtyp = t_float;
  args.list   = list;
  var = def;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}


void CommandLine::AddSwitch(const char *arg, char *var, char *def, bool list)
{
  ArgumentInfo args;
  args.argvar = (void *) var;
  args.argtyp = t_char;
  args.list   = list;
  var = strdup(def);
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}


void CommandLine::AddSwitch(const char *arg, string &var, string def, bool list)
{
  ArgumentInfo args;
  args.argvar = (void *) &var;
  args.argtyp = t_string;
  args.list   = list;
  var = def;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}


void CommandLine::AddSwitch(const char *arg, vector<int> &var, bool list)
{
  ArgumentInfo args;
  args.argvar = (void *) &var;
  args.argtyp = t_vector_int;
  args.list   = list;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}


void CommandLine::AddSwitch(const char *arg, vector<float> &var, bool list)
{
  ArgumentInfo args;
  args.argvar = (void *) &var;
  args.argtyp = t_vector_float;
  args.list   = list;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}


void CommandLine::AddSwitch(const char *arg, vector<string> &var, bool list)
{
  ArgumentInfo args;
  args.argvar = (void *) &var;
  args.argtyp = t_vector_string;
  args.list   = list;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}

void CommandLine::AddSwitch(const char *arg, void (*fxn)(void *), bool list)
{
  ArgumentInfo args;
  args.argvar  = (void *) fxn;
  args.argtyp  = t_func;
  args.list    = list;
  args.numargs = 0;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}
 

void CommandLine::SetVariable(ArgumentInfo &argInfo, char *arg)
{
  switch (argInfo.argtyp)
    {
    case t_bool:
      *((bool *)argInfo.argvar) = true;
      break;
    case t_int:
      *((int *)argInfo.argvar) = atoi(arg);      
      break;
    case t_float:
      *((float *)argInfo.argvar) = atof(arg);
      break;
    case t_char:
      argInfo.argvar = strdup(arg);
      break;
    case t_string:
      *((string *)argInfo.argvar) = arg;
      break;
    case t_vector_int:
      ((vector<int> *)argInfo.argvar)->push_back(atoi(arg));
      break;
    case t_vector_float:
      ((vector<float> *)argInfo.argvar)->push_back(atof(arg));
      break;
    case t_vector_string:
      ((vector<string> *)argInfo.argvar)->push_back(string(arg));
      break;
    }
  if ( argInfo.list && arg && strcmp(arg, "") != 0 )
    arglist += " " + string(arg);
} 


int CommandLine::ProcessArguments( char *argv[], int argc )
{
  int i, j;
  unsigned int ii;
  maptype::iterator iter;
  char **argvP= NULL;

  vector<void (*)(void *)> functions;

  // check to see if we are reading from a file
  if (argc > 1 ){
      // if we have at at sign then 
      int j= 1;
      if (argv[1][0] == '@') {
          // we
          if (strlen(argv[1]) > 1) {
             ifstream inp;
             inp.open(&(argv[1][1]));
             if (!inp){
                 cerr << "Unable to open input file " << &argv[1][1] << endl;
                 exit(0);
             }
             char s[1000];
             vector <string> argvS;
             while(inp.getline(s,1000)){
                vector <string> res;
                strcat(s," ");
                tokenize(res,s);
                for(ii=0; ii != res.size(); ii++){
                    argvS.push_back(res[ii]);
                }
             }
             argvP= (char **) calloc(argc+argvS.size(),sizeof(char*));
             for(j= 1, i= 0; (unsigned)i < argvS.size(); i++)
	       {
                argvP[j]=  (char*) calloc(argvS[i].length()+1,sizeof(char));
                strcpy(argvP[j],argvS[i].c_str());
                j++;
             }
             inp.close();
          }
          argv= argvP;
          argc= j;
      }
  }

  for ( i = 1 ; i < argc ; i++ )
    {
      iter = cmdline.find(string(argv[i]));
      if ( iter == cmdline.end() )
	{
	  string errstr = "Unknown flag: "; errstr += argv[i];
	  ThrowError(errstr);
	  Usage();
	  if (_exitOnError) { exit(-1); }
	  continue;
	}
      (*iter).second.found = true;

      if ( (*iter).second.list )
	arglist += " " + string(argv[i]);

      if ( (*iter).second.argtyp >= t_vector_int && 
	   (*iter).second.argtyp <= t_vector_string )
	{
	  for ( j = i ; j + 1 < argc ; j++ )
	    {
	      if (argv[j+1][0] == '-' && (argv[j+1][1] < '0' || argv[j+1][1] > '9'))
		break;
	      SetVariable((*iter).second, argv[j+1]);
	    }
	  i = j;
	}
      else if ( (*iter).second.argtyp == t_func )
	{
	  void (*fxn)(void *) = (void (*)(void *))(*iter).second.argvar;
	  functions.push_back(fxn);
	}
      else if ( (*iter).second.numargs == 0 )
	{
	  if ( (*iter).second.argtyp == t_bool )
	    SetVariable((*iter).second,"");
	}
      else if ( (*iter).second.numargs >= 1 )
	{
	  for ( j = 0 ; j < (*iter).second.numargs ; j++ )
	    {
	      if ( (*iter).second.argtyp == t_bool )
		SetVariable((*iter).second, argv[i]);
	      else if (i + 1 < argc && !(argv[i+1][0] == '-' && (argv[i+1][1] < '0' || argv[i+1][1] > '9')))
		{
		  (*iter).second.hasArgument = true;
		  (*iter).second.arguments.push_back(argv[i+1]);
		  SetVariable((*iter).second, argv[i+1]);
		  i++;
		}
	      else
		{
		  cerr << (*iter).first << ": expected argument" << endl;
		  Usage();
		  return -1;
		}
	    }
	}
    }

  for ( i = 0 ; (unsigned)i < functions.size() ; i++ )
    functions[i]((void *)this);

  return 1;
}


bool CommandLine::WasCalledWith(const char *arg)
{
  maptype::iterator iter = cmdline.find(string(arg));
  if ( iter != cmdline.end() )
    return (*iter).second.found;
  else
    return false;
}

bool CommandLine::HasArgument(const char *arg)
{
  maptype::iterator iter = cmdline.find(string(arg));
  if ( iter != cmdline.end() )
    return (*iter).second.hasArgument;
  else
    return false;
}

char *CommandLine::GetArgument(const char *arg, int i)
{
  maptype::iterator iter = cmdline.find(string(arg));
  if ( iter != cmdline.end() )
    return (char *)(*iter).second.arguments[i].c_str();
  return NULL;
}

bool CommandLine::getValue(const char *arg, float &v)
{
  bool ret= false;
  maptype::iterator iter = cmdline.find(string(arg));
  if ( iter != cmdline.end() ){
      if ((*iter).second.argtyp == t_float) {
          void *vp= (*iter).second.argvar;
          v= *((float*)(vp));
          ret= true;
      }
  }
  return ret;
}

bool CommandLine::getValue(const char *arg, int &v)
{
  bool ret= false;
  maptype::iterator iter = cmdline.find(string(arg));
  if ( iter != cmdline.end() ){
      if ((*iter).second.argtyp == t_int) {
          void *vp= (*iter).second.argvar;
          v= *((int*)(vp));
          ret= true;
      }
  }
  return ret;
}

bool CommandLine::getValue(const char *arg, string &s)
{
  bool ret= false;
  maptype::iterator iter = cmdline.find(string(arg));
  if ( iter != cmdline.end() ){
      if ((*iter).second.argtyp == t_string) {
          void *vp= (*iter).second.argvar;
          s= *((string*)(vp));
          ret= true;
      }
  }
  return ret;
}

void CommandLine::SetUsageFunction(void (*usageFxn)(void))
{
  usage = usageFxn;
}

void CommandLine::Usage()
{
  if ( usage ) usage();
  else         printHelp();
}

void CommandLine::AddSwitch(const char *arg, string &var, string def, const char *help, bool list)
{
  const string hlp= string(help);
  ArgumentInfo args(hlp);
  args.argvar = (void *) &var;
  args.argtyp = t_string;
  args.list   = list;
  var = def;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}

void CommandLine::AddSwitch(const char *arg, int &var, int def, const char *help, bool list)
{
  const string hlp= string(help);
  ArgumentInfo args(hlp);
  args.argvar = (void *) &var;
  args.argtyp = t_int;
  args.list   = list;
  var = def;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}

void CommandLine::AddSwitch(const char *arg, float &var, float def, const char *help, bool list)
{
  const string hlp= string(help);
  ArgumentInfo args(hlp);
  args.argvar = (void *) &var;
  args.argtyp = t_float;
  args.list   = list;
  var = def;
  cmdline.insert(make_pair(string(arg), args));
  switches++;
}

void CommandLine::printHelp()
{
  maptype::iterator iter;
  cerr << " Available Commands " << endl;

  for(iter= cmdline.begin(); iter != cmdline.end(); iter++){
      cerr << (*iter).first << ": " << (*iter).second.helpString << endl;
  }

  cerr << "@fname:  read input from file fname" << endl;
}


