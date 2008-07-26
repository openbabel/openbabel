/**********************************************************************
plugin.cpp - facilitates construction of plugin classes
 
Copyright (C) 2007 by Chris Morley
 
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

#include <openbabel/babelconfig.h>
#include <openbabel/plugin.h>

#include <iterator>

using namespace std;
namespace OpenBabel
{

  OBPlugin::PluginMapType& OBPlugin::GetTypeMap(const char* PluginID)
{
  PluginMapType::iterator itr;
  itr = PluginMap().find(PluginID);
  if(itr!=PluginMap().end())
    return itr->second->GetMap();
  return PluginMap();//error: type not found; return plugins map
}

bool OBPlugin::ListAsVector(const char* PluginID, const char* param, vector<string>& vlist)
{
  PluginMapType::iterator itr;
  bool ret=true;
  if(PluginID)
  {
    if(*PluginID!=0 && strcmp(PluginID, "plugins"))
    {
      //List the sub classes of the specified type
      itr = PluginMap().find(PluginID);
      if(itr!=PluginMap().end())
      {
        bool onlyIDs = param!=NULL && strstr(param,"ids")!=NULL;
        //Get map of plugin type (like OBFingerprint) and output its contents
        PluginMapType Map = itr->second->GetMap();
        for(itr=Map.begin(); itr!=Map.end(); ++itr)
        {
          if(*(itr->first)=='_')//no listing when ID starts with '_'
            continue;
          if(onlyIDs)
            vlist.push_back(itr->first);
          else
          {
            string txt;
            if((itr->second)->Display(txt, param, itr->first))
              vlist.push_back(txt);
          }
        }
        return true;
      }
      ret=false; //asked for a type not available; provide plugin types instead
    }
  }
  //List the plugin types
  for(itr=PluginMap().begin();itr!= PluginMap().end();++itr)
    vlist.push_back(itr->first);
  return ret;
}

void OBPlugin::List(const char* PluginID, const char* param, ostream* os)
{
  vector<string> vlist;
  if(!ListAsVector(PluginID,param, vlist))
    *os << PluginID << " is not a recognized plugin type. Those with instances of sub-types loaded are:" << endl;
  copy(vlist.begin(), vlist.end(), std::ostream_iterator<string>(*os, "\n"));
}

string OBPlugin::ListAsString(const char* PluginID, const char* param)
{
  stringstream ss;
  List(PluginID, param, &ss);
  return ss.str();
}

string OBPlugin::FirstLine(const char* txt)
{
  string stxt(txt);
  string::size_type pos = stxt.find('\n');
  if(pos==string::npos)
    return stxt;
  else
    return stxt.substr(0,pos);
}

//Default version
bool OBPlugin::Display(string& txt, const char* param, const char* ID)
{
  //Use the provided ID if possible.
  if(ID)
    txt = ID;
  else
    txt = GetID();
  txt += "    ";// was '\t'; but caused problems in GUI menu
  if(param && !strcasecmp(param, "verbose"))
    txt += Description();
  else
    txt += FirstLine(Description());
  return true;
}

}//namespace

//! \file plugin.cpp
//! \brief Simplify 'plugin' classes to be discovered and/or loaded at runtime
