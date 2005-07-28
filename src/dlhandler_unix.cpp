/**********************************************************************
dlhandler_unix.cpp - Dynamic loader for UNIX (handles file format shared obj.)

Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2004-2005 by Geoffrey R. Hutchison
 
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

#include "dlhandler.h"
#include "babelconfig.h"

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <dlfcn.h>

using namespace std;

//Globals for scandir()
string targetPattern;

int matchFiles (SCANDIR_CONST struct dirent *entry_p)
{
  return (strstr(entry_p->d_name, targetPattern.c_str()) != 0);
}

bool DLHandler::getConvDirectory(string& convPath)
{
    //Need to provide the directory from which this shared library was loaded.
    //This is the default directory for format shared library files.
  
  string testPath;
  testPath += OB_MODULE_PATH;
  convPath = testPath;

  return true;
}

int DLHandler::findFiles (vector <string>& file_list,const string& pattern, const string& path)
{
  string internalPath = path;

  if (path.empty())
    internalPath="./"; //use current directory if path is empty

  targetPattern = pattern; //make accessible to global function
  struct dirent **entries_pp;
  int count = scandir (internalPath.c_str(), &entries_pp, SCANDIR_T matchFiles, NULL);

  for(int i=0; i<count; i++)
    file_list.push_back(path + getSeparator() + (entries_pp[i])->d_name);

  return count;
}

bool DLHandler::openLib(const string& lib_name)
{
    return dlopen(lib_name.c_str(), RTLD_LAZY) != 0;
}

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
