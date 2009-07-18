/**********************************************************************
dlhandler_win32.cpp - Dynamic loader for Win32 (handles file format DDLs)

Copyright (C) 2004 by Chris Morley
 
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
#include <openbabel/babelconfig.h>

#include <vector>
#include <cstdarg>
#include <iostream>
//# define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
//#include <winbase.h>
//#include <direct.h>
#include <windows.h>
#include <openbabel/dlhandler.h>
using namespace std;

namespace OpenBabel {
  OBAPI bool tokenize(vector<string> &, const char *, const char *);
}

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

bool DLHandler::getConvDirectory(string& convPath)
{
    char path[MAX_PATH+1];
    if (!GetModuleFileName(0, path, MAX_PATH)) //gets exe name NEEDS REVISITING 
      return false;

    // strip of appname.exe
    convPath = path;
    std::string::size_type p = convPath.rfind('\\');
    if (p == std::string::npos)
      return false;

    convPath = convPath.substr(0, p+1);
#if defined(__CYGWIN__) || defined(__MINGW32__)
    convPath += "..\\lib\\openbabel\\";
    convPath += BABEL_VERSION;
    convPath += "\\";
#endif

    return true;
}


int DLHandler :: findFiles (std::vector<std::string>& file_list,const std::string &filename)
{
    if(filename.find_first_of("*?")==string::npos)
    {
        //no wildcard in filename
        file_list.push_back(filename);
        return -1;
    }
    size_t pos = filename.find_last_of("\\/");
    if(pos!=string::npos)
        return findFiles(file_list,filename.substr(pos+1), filename.substr(0,pos+1));
    else
        return findFiles(file_list,filename, "");
}

int DLHandler :: findFiles (std::vector<std::string>& file_list,const std::string &pattern,const std::string &path)
{
  vector<string> paths, vs;
  char buffer[BUFF_SIZE];
  
  if (!path.empty())
    paths.push_back(path);
 
  /*
  if (getenv("BABEL_LIBDIR") != NULL) 
  {
    // environment variable should override built-in path
    paths.clear();
      
    strncpy(buffer,getenv("BABEL_LIBDIR"), BUFF_SIZE - 1);
    // add a trailing NULL just in case
    buffer[BUFF_SIZE - 1] = '\0';
      
    OpenBabel::tokenize(vs, buffer, "\r\n\t :");
      
    if (!vs.empty())
    {
      for (unsigned int i = 0; i < vs.size(); ++i) {
        paths.push_back(vs[i]);
	  cout << "path[" << i << "] = " << paths[i] << endl;
      }
    }
  }
  */
  
  if (paths.empty())
    paths.push_back("./"); // defaults to current directory
  
  string currentPath;
  for (unsigned int i = 0; i < paths.size(); ++i)
  {
    currentPath = paths.at(i);
    WIN32_FIND_DATA file_data;
    HANDLE handle;
    handle = FindFirstFile ((currentPath + pattern).c_str(), &file_data);
    while(handle!=INVALID_HANDLE_VALUE)
    {
      ULONG value = (file_data.dwFileAttributes) & FILE_ATTRIBUTE_DIRECTORY;
      if (value != FILE_ATTRIBUTE_DIRECTORY)
      {
        string s = file_data.cFileName;
        file_list.push_back(currentPath + s);
      }
      if(!FindNextFile (handle, &file_data))
        break;
    }
    if (! (FindClose (handle)))
      return 0;// couldn't close search handle
  }

  return file_list.size();
}

const char* DLHandler::getFormatFilePattern()
{
    return "*.obf";
}


bool DLHandler :: openLib(const string& lib_name)
{

    if(LoadLibrary(lib_name.c_str()))
        return true;

    unsigned long err = GetLastError();
    return false;
}

char DLHandler::getSeparator()
{
    return '\\';
}

void DLHandler::Sleep(int n)
{
	::Sleep(n);
}

//! \file dlhandler_win32.cpp
//! \brief Dynamic loader for Win32 (handles file format DDLs)
