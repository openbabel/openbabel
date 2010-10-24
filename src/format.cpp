/**********************************************************************
obconversion.cpp -  Definitions for OBFormat

Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2005-2007 by Geoffrey Hutchison

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
#include <openbabel/babelconfig.h>
#include <openbabel/format.h>
#include <typeinfo>

using namespace std;
namespace OpenBabel
{
#if defined(__CYGWIN__) || defined(__MINGW32__)
  // macro to implement static OBPlugin::PluginMapType& Map()
  PLUGIN_CPP_FILE(OBFormat)
#endif

int OBFormat::RegisterFormat(const char* ID, const char* MIME)
{
  GetMap()[ID] = this;
  if (MIME)
    FormatsMIMEMap()[MIME] = this;
  if(Flags() & DEFAULTFORMAT)
    Default() = this;

  //ensure "formats" is registered as a plugin
  PluginMap()[TypeID()] =this;
  _id=ID;
  return GetMap().size();
}

//////////////////////////////////////////////////////////
const char* OBFormat::TargetClassDescription()
{
  //Provides class of default format unless overridden
  if(OBFormat::FindType(NULL))
    return OBFormat::FindType(NULL)->TargetClassDescription();
  else
    return "";
}

//////////////////////////////////////////////////////////
const type_info& OBFormat::GetType()
{
  //Provides info on class of default format unless overridden
  if(OBFormat::FindType(NULL))
    return OBFormat::FindType(NULL)->GetType();
  else
    return typeid(this); //rubbish return if DefaultFormat not set
}

//////////////////////////////////////////////////////////
OBFormat* OBFormat::FormatFromMIME(const char* MIME)
{
  if(FormatsMIMEMap().find(MIME) == FormatsMIMEMap().end())
    return NULL;
  else
    return static_cast<OBFormat*>(FormatsMIMEMap()[MIME]);
}

//////////////////////////////////////////////////////////
bool OBFormat::Display(std::string& txt, const char* param, const char* ID)
{
  //No output for formats which can't be written or read
  if((Flags() & NOTREADABLE) && (Flags() & NOTWRITABLE))
    return false;

  bool justread=false, justwrite=false;
  //No output if formats is not readable or writable if this was requested
  if(param)
  {
    if((!strncasecmp(param, "in", 2) || !strncasecmp(param, "read",4)))
    {
      justread=true;
      if(Flags() & NOTREADABLE)
        return false;
    }
    if((!strncasecmp(param, "out",3) || !strncasecmp(param, "write",5)))
    {
      justwrite=true;
      if(Flags() & NOTWRITABLE)
        return false;
    }
  }

  //Use the provided ID if possible. If more than one ID has been registed
  //for the format, e.g. "smiles" and "smi", the contents of the member
  //variable _id, returned by GetID() is the last one.
  if(ID)
    txt = ID;
  else
    txt = GetID();
  txt += " -- ";
  txt += FirstLine(Description());
  if(!justread && (Flags() & NOTWRITABLE))
    txt += " [Read-only]";
  if(!justwrite && (Flags() & NOTREADABLE))
    txt += " [Write-only]";

  if(param && strstr(param, "verbose"))
  {
    const char* nl = strchr(Description(), '\n');
    if(nl)
    {
      txt += '\n';
      txt += ++nl; // add rest of description
      if(strlen(SpecificationURL()))
      {
        txt += "\nSpecification at: ";
        txt += SpecificationURL();
      }
      txt += "\n";
    }
  }
  return true;//means txt has been updated
}
}//namespace

//! \file format.cpp
//! \brief Base class OBFormat for file formats
