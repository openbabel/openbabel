/**********************************************************************
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

#include <cstdarg>
#include <iostream>
# define WIN32_LEAN_AND_MEAN
//#include "winbase.h"
#include <direct.h>
#include "windows.h"
#include "dlhandler.h"
using namespace std;



bool DLHandler::getConvDirectory(string& convPath)
{
	char path[MAX_PATH+1];

	HMODULE hConv = GetModuleHandle("OBConv");
	if(GetModuleFileName(hConv,path,MAX_PATH)==0) return false;
	char* p = strrchr(path,'\\');
	if(p && *p) *(p+1)='\0'; //chop off name after last '\'
	convPath=path;
	return true;
}


int DLHandler :: findFiles (vector<string>& file_list,const string &filename)
{
	if(filename.find_first_of("*?")==string::npos)
	{
		//no wildcard in filename
		file_list.push_back(filename);
		return -1;
	}
	int pos = filename.find_last_of("\\/");
	if(pos!=string::npos)
		return findFiles(file_list,filename.substr(pos+1), filename.substr(0,pos+1));
	else
		return findFiles(file_list,filename, "");
}

int DLHandler :: findFiles (vector<string>& file_list,const string &pattern,const string &path)
{
	WIN32_FIND_DATA file_data;
	HANDLE handle;
	if ((handle = FindFirstFile ((path + pattern).c_str(), &file_data)) == INVALID_HANDLE_VALUE)
		return 0;

	ULONG value = (file_data.dwFileAttributes) & FILE_ATTRIBUTE_DIRECTORY;
	if (value != FILE_ATTRIBUTE_DIRECTORY)
	{
		string s = file_data.cFileName;
		file_list.push_back(path + s);

		while (FindNextFile (handle, &file_data))
		{
			value = (file_data.dwFileAttributes) & FILE_ATTRIBUTE_DIRECTORY;
			if (value != FILE_ATTRIBUTE_DIRECTORY)
			{
				s = file_data.cFileName;
				file_list.push_back(path + s);
			}
		}
	}

	if (! (FindClose (handle)))
		return 0;// couldn't close search handle
	return file_list.size();
}

const char* DLHandler::getFormatFilePattern()
{
	return "*.obf";
}


bool DLHandler :: openLib(const string& lib_name)
{
	
	if(LoadLibrary(lib_name.c_str())) return true;

	unsigned long err = GetLastError();
	return false;
}

char DLHandler::getSeparator()
{ return '\\';}
