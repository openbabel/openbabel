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

#include <algorithm>
#include <vector>
#include <string>

using namespace std;

bool tokenize(vector<string> &, char *,char *);
char *trim_spaces(char *string);
bool tokenize(vector<string> &vcr, string &s,char *delimstr,int limit=-1);

namespace OpenEye {

bool tokenize(vector<string> &vcr, char *buf,char *delimstr)
{
    vcr.clear();
    string s = buf;
    s += "\n";
    size_t startpos=0,endpos=0;
    
    for (;;)
    {
	startpos = s.find_first_not_of(delimstr,startpos);
	endpos = s.find_first_of(delimstr,startpos);

	if (endpos <= s.size() && startpos <= s.size())
	    vcr.push_back(s.substr(startpos,endpos-startpos));
	else
	    break;
	
	startpos = endpos+1;
    }

    return(true);
}


char *trim_spaces(char *string)
{
  int length;

  length = strlen(string);
  if (length == 0)
    return string;

  while ((length > 0) && (string[0] == ' '))
  {
    string++;
    --length;
  }

  if (length > 0)
  {
    while ((length > 0) && (string[length-1] == ' '))
    {
      string[length-1] = '\0';
      --length;
    }
  }

  return(string);
}

bool tokenize(vector<string> &vcr, string &s,char *delimstr,int limit)
{
  vcr.clear();
  size_t startpos=0,endpos=0;

  int matched=0;
  unsigned int s_size = s.size();
  for (;;)
    {
        startpos = s.find_first_not_of(delimstr,startpos);
        endpos = s.find_first_of(delimstr,startpos);
        if (endpos <= s_size && startpos <= s_size)
          {
            vcr.push_back(s.substr(startpos,endpos-startpos));

            matched++;
            if (matched == limit)
              {
                startpos = endpos+1;
                vcr.push_back(s.substr(startpos,s_size));
                break;
              }
          }
        else
          {
            if (startpos < s_size)
              vcr.push_back(s.substr(startpos,s_size-startpos));
            break;
          }

        startpos = endpos+1;
    }
    return(true);
}

}



