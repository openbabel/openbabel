/**********************************************************************
dlhandler_unix.cpp - Dynamic loader for UNIX (handles file format shared obj.)

Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2004-2006 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifdef __MINGW32__
 #include <windows.h>
#else
 #include <dlfcn.h>
#endif

#include <openbabel/dlhandler.h>
#include <openbabel/babelconfig.h>
#include <openbabel/oberror.h>

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include <iostream>

using namespace std;

namespace OpenBabel {
  OBAPI bool tokenize(vector<string> &, const char *, const char *);
}

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

//Globals for scandir()

int matchFiles (SCANDIR_CONST struct dirent *entry_p)
{
	string filename(entry_p->d_name);
	string::size_type extPos = filename.rfind(DLHandler::getFormatFilePattern());

	if(extPos!=string::npos && filename.substr(extPos) == DLHandler::getFormatFilePattern())
		return true;

	return false;
}

bool DLHandler::getConvDirectory(string& convPath)
{
  //Need to provide the directory from which this shared library was loaded.
  //This is the default directory for format shared library files.

  string testPath;
  testPath += OB_MODULE_PATH; // defined in src/config.h.cmake -> babelconfig.h
  convPath = testPath;

  return true;
}

int DLHandler::findFiles (std::vector <std::string>& file_list,
                          const std::string& pattern,
                          const std::string& path)
{
  vector<string> paths, vs;
  char buffer[BUFF_SIZE];

  if (!path.empty())
    paths.push_back(path);

  if (getenv("BABEL_LIBDIR") != NULL)
    {
      // environment variable should override built-in path
      paths.clear();

      strncpy(buffer,getenv("BABEL_LIBDIR"), BUFF_SIZE - 1);
      // add a trailing NULL just in case
      buffer[BUFF_SIZE - 1] = '\0';

      OpenBabel::tokenize(vs, buffer, "\r\n:");

      if (!vs.empty())
        {
          for (unsigned int i = 0; i < vs.size(); ++i) {
            paths.push_back(vs[i]);
          }
        }
    }

  if (paths.empty())
    paths.push_back("./"); // defaults to current directory

  /* Our old code used scandir. Replaced with readdir (below) as for example
   * Solaris pre 10 doesn't implement scandir.
   */

  DIR *dp;
  struct dirent *entry;

  string currentPath;
  for (unsigned int i = 0; i < paths.size(); ++i)
    {
      currentPath=paths[i];

      if ((dp = opendir(currentPath.c_str())) == NULL)
        continue; // no big deal, this path causes an error
      else
        {
          while((entry = readdir(dp)) != NULL)
            {
              if (matchFiles(entry) != 0)
                file_list.push_back(currentPath + getSeparator() + (entry)->d_name);
            }
          closedir(dp); // calls free(dp) -- no memory leak
        }
    }

  if (file_list.empty())
    return(-1); // error, didn't find any files at all
  return file_list.size();
}

int DLHandler::findFiles (std::vector<std::string>& file_list,
                          const std::string &filename)
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

#ifdef __MINGW32__
bool DLHandler :: openLib(const string& lib_name)
{

    if(LoadLibrary(lib_name.c_str()))
        return true;

    unsigned long err = GetLastError();
    return false;
}
#else
bool DLHandler::openLib(const string& lib_name)
{
  bool success = (dlopen(lib_name.c_str(), RTLD_LAZY | RTLD_GLOBAL) != 0);
  if (!success) {
    char buffer[BUFF_SIZE];
    sprintf(buffer, "%s did not load properly.\n Error: %s",
              lib_name.c_str(), dlerror());
    OpenBabel::obErrorLog.ThrowError(__FUNCTION__, buffer, OpenBabel::obError);
  }
  return success;
}
#endif

const char* DLHandler::getFormatFilePattern()
{
  return MODULE_EXTENSION;
}

char DLHandler::getSeparator()
{
  return '/';
}

//! \file dlhandler_unix.cpp
//! \brief Dynamic loader for UNIX (handles file format shared obj.)
