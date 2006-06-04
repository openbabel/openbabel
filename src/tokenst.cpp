/**********************************************************************
tokenst.cpp - Tokenize a string.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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
#include "babelconfig.h"

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include <algorithm>
#include <vector>
#include <string>


using namespace std;
/*
OBAPI bool tokenize(vector<string> &, const char *, const char *);
OBAPI char *trim_spaces(char *string);
OBAPI bool tokenize(vector<string> &vcr, string &s, const char *delimstr,int limit=-1);
*/
namespace OpenBabel
{

  //! Break a string (supplied as the second argument) into tokens, returned 
  //! in the first argument. Tokens are determined by the delimiters supplied
  //! (defaults to whitespace (i.e., spaces, tabs, newlines)
OBAPI bool tokenize(vector<string> &vcr, const char *buf, const char *delimstr)
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

  //! Trim any trailing spaces at the end of the supplied string.
OBAPI char *trim_spaces(char *string)
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

//! Break a string (supplied as the second argument) into tokens, returned 
//! in the first argument. Tokens are determined by the delimiters supplied
//! (defaults to whitespace (i.e., spaces, tabs, newlines)
//! Only breaks at most 'limit' tokens and the last item in the vector may
//! include un-parsed tokens.
OBAPI bool tokenize(vector<string> &vcr, string &s, const char *delimstr, int limit)
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

/// Removes white space from front and back of string
OBAPI string& Trim(string& txt)
{
	string::size_type pos = txt.find_last_not_of(" \t\n\r");
	if(pos!=string::npos)
		txt.erase(pos+1);
	else
		txt.erase();

	pos = txt.find_first_not_of(" \t\n\r");
	if(pos!=string::npos)
		txt.erase(0, pos);
	else
		txt.erase();
	return txt;
}

} // end namespace OpenBabel

//! \file tokenst.cpp
//! \brief Tokenize a string.
