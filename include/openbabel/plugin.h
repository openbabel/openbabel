/**********************************************************************
plugin.h - facilitates construction of plugin classes
 
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

#ifndef OB_PLUGIN_H
#define OB_PLUGIN_H

#include <openbabel/babelconfig.h>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <cstring>

#ifndef OBERROR
 #define OBERROR
#endif

namespace OpenBabel
{

/// @brief Case insensitive string comparison for PluginMapType key.
struct OBERROR CharPtrLess : public std::binary_function<const char*,const char*, bool>
{
  bool operator()(const char* p1,const char* p2) const
  { return strcasecmp(p1,p2)<0; }
};

/** \class OBPlugin plugin.h <openbabel/plugin.h>
    \brief Base class for all types of dynamic classes discovered at runtime
    \since version 2.2
    \sa plugin.h
 */
class OBERROR OBPlugin
{
public:

  //Maps of thistype are used to store 
  // (a)a list of the plugin types in OBPlugin, and
  // (b)a list of the sub-types in each type class derived from OBPlugin.
  typedef std::map<const char*, OBPlugin*, CharPtrLess> PluginMapType;
  typedef PluginMapType::const_iterator PluginIterator;

  ///Virtual destructor necessary for classes with virtual functions
  virtual ~OBPlugin(){};

  ///Required description of a sub-type
  virtual const char* Description() { return NULL;} ;

  ///Write information on a plugin class to the string txt.
  ///Return false if not written.
  ///The default implementation outputs:
  /// the ID, a tab character, and the first line of the Description.
  ///The param string can be used in derived types to provide different outputs.
  virtual bool Display(std::string&txt, const char* param, const char* ID=NULL);

  ///Make a new instance of the class.
  ///See OpTransform, OBGroupContrib, SmartsDescriptor classes for derived versions.
  ///Usually, the first parameter is the classname, the next three are
  ///parameters(ID, filename, description) for a constructor, and the rest data.
  virtual OBPlugin* MakeInstance(const std::vector<std::string>&){return NULL;}

  ///Get a pointer to a plugin from its type and ID. Return NULL if not found. Not cast to Type*
  static OBPlugin* GetPlugin(const char* Type, const char* ID)
  { return BaseFindType(GetTypeMap(Type), ID); }

  ///Return the ID of the sub-type instance.
  const char* GetID()const{return _id;};

  ///Output a list of sub-type classes of the the type PluginID,
  ///or, if PluginID is "plugins" or empty, a list of the base types.
  ///If PluginID is not recognized or is NULL, the base types are output and the return is false.
  static bool ListAsVector(const char* PluginID, const char* param, std::vector<std::string>& vlist);

  ///As ListAsVector but sent to an ostream with a default of cout if not specified
  static void List(const char* PluginID, const char* param=NULL, std::ostream* os=&std::cout);

  ///As ListAsVector but returns a string containing the list
  static std::string ListAsString(const char* PluginID, const char* param=NULL);

  ///Utility function to return only the first line of a string
  static std::string FirstLine(const char* txt);

  ///Return an iterator at the start of the map of the plugin types PluginID or,
  /// if there is no such map, the end of the top level plugin map.
  static PluginIterator Begin(const char* PluginID)
  {
    if( !strcmp(PluginID, "plugins") || GetTypeMap(PluginID)!=PluginMap())
      return GetTypeMap(PluginID).begin();
    else
      return PluginMap().end();
  }

  static PluginIterator End(const char* PluginID)
  {
    return GetTypeMap(PluginID).end();
  }

  ///Returns the map of the subtypes
  virtual PluginMapType& GetMap() const =0;

protected:
  ///\brief Returns a reference to the map of the plugin types.
  /// Is a function rather than a static member variable to avoid initialization problems. 
  static PluginMapType& PluginMap()
  {
    static PluginMapType m;
    return m;
  }

  ///Returns the map of a particular plugin type, e.g. GetMapType("fingerprints")
  static PluginMapType& GetTypeMap(const char* PluginID);

  ///\brief Returns the type with the specified ID, or NULL if not found.
  ///Will be cast to the appropriate class in the calling routine.
  static OBPlugin* BaseFindType(PluginMapType& Map, const char* ID)
  {
    PluginMapType::iterator itr = Map.find(ID);
    if(itr==Map.end())
      return NULL;
    else
      return itr->second;
  }

protected:
  const char* _id;
};

#if defined(__CYGWIN__) || defined(__MINGW32__)

//Macro to be added to definition of the base class
#define MAKE_PLUGIN(BaseClass)\
protected:\
static PluginMapType& Map();\
virtual PluginMapType& GetMap()const{return Map();}\
public:\
static BaseClass*& Default(){static BaseClass* d;return d;}\
  BaseClass(const char* ID, bool IsDefault=false)\
 {_id=ID;if(ID&&*ID){if(IsDefault || Map().empty()) Default() = this;\
 Map()[ID]=this;PluginMap()[TypeID()] =this;}}\
static BaseClass* FindType(const char* ID)\
 {if(!ID || *ID==0) return Default();\
 return static_cast<BaseClass*>(BaseFindType(Map(),ID));}

#define PLUGIN_CPP_FILE(BaseClass)\
OBPlugin::PluginMapType& BaseClass::Map()\
{ static OBPlugin::PluginMapType map; return map; }

#else // __CYGWIN__ || __MINGW32__

//Macro to be added to definition of the base class
#define MAKE_PLUGIN(BaseClass)\
protected:\
static PluginMapType& Map(){static PluginMapType m;return m;}\
virtual PluginMapType& GetMap()const{return Map();}\
public:\
static BaseClass*& Default(){static BaseClass* d;return d;}\
  BaseClass(const char* ID, bool IsDefault=false)\
 {_id=ID;if(ID&&*ID){if(IsDefault || Map().empty()) Default() = this;\
 Map()[ID]=this;PluginMap()[TypeID()] =this;}}\
static BaseClass* FindType(const char* ID)\
 {if(!ID || *ID==0) return Default();\
 return static_cast<BaseClass*>(BaseFindType(Map(),ID));}

#endif // __CYGWIN__ || __MINGW32__

/** \file plugin.h
   \brief Simplify 'plugin' classes to be discovered and/or loaded at runtime.

The code in this file makes it easy to make 'plugin' classes. These classes are
derived from a base class, like OBFingerprint. The derived classes 
('sub-types' like fingerprint2) usually have a single instance. Plugin classes 
are only discovered at runtime, so no existing code needs to be changed when 
adding a new derived class. In some builds the new code can be added or removed 
by just moving a DLL or so file. The plugin classes derived from any base class
(including new ones) type can be listed from the commandline.

<h2>Step-by-Step Instructions</h2>

1) In the header file for YourBaseClass (which handles whatsits).
Make sure to include the plugin.h header ,
derive the class from OBPlugin
and in its definition add 
  the MAKE_PLUGIN macro 
  and a function TypeID() containing a simple descriptor of the type
\code
#include <openbabel/plugin.h>
class YourBaseClass : public OBPlugin
{
  MAKE_PLUGIN(YourBaseClass)
  
  const char* TypeID()
  { return "whatsits"; };

  ...rest of implementation, probably involving virtual functions redefined
  in the sub-type classes
};
\endcode
See below for what the macro contains.

2) Declare each sub-type in a class derived from the base class
and give it a constructor which calls OBPlugin constructor as shown:
\code
class YourSubType1 : public YourBaseClass
{
public:
  YourSubtype1(const char* ID, bool IsDefault=false) 
    : YourBaseClass(ID, IsDefault){}

  virtual string Description()
  { return "A description with one or more lines";};

  ...rest of implementation
};
\endcode
Only the first line of the description is used when the subclasses are listed.

3) Declare a global instance of the sub-type class which specifies its ID.
and, optionally, whether it is to be regarded as the default type of YourBaseClass.
\code
YourSubType1 theType1("whatsit2",true);
\endcode

4) The following functions are available:

YourBaseClass* YourBaseClass::FindType(const char* ID);
This returns the default type when ID is NULL or empty.

To list the sub-types of any plugin class use the List which sends to cout 
by default (or any other ostream if specified).
\code
  OBPlugin::List("whatsits")
\endcode
The ListAsString and ListAsVector functions are alternatives, usable with scripting.

It is also possible to iterate through each sub-type by the following code:
\code
  OBPlugin::PluginIterator itr;
  for(itr=OBPlugin::Begin("whatsits");itr!=OBPlugin::End("whatsits");++itr)
  {
    itr is a std::map::const_iterator
    itr->first is the ID of the subtype;
    itr->second is The OBPlugin* which you will have to cast to your type
  }
\endcode
Since this is not the most beautiful code, it is probably better to use the
List methods if possible.

YourBaseClass* MakeNewInstance();



<h2>How it works</h2>

MAKE_PLUGIN(YourBaseClass) inserts the following code into YourBaseClass:
\code
protected:
  
  //The collection of sub-types is in a local static variable to avoid
  //any difficulties with the order of initialization of static objects. 
  static PluginMapType& Map()
  {
    static PluginMapType m;
    return m;
  }

  //Making the map accessible to the base class (Cannot be used during construction)
  virtual PluginMapType& GetMap()const
  {
   return Map();
  }
   
public:
  static YourBaseClass*& Default()
  {
    static YourBaseClass* d;
    return d;
  }

  //Constructor registers the sub-type 
  YourBaseClass(const char* ID, bool IsDefault=false)
  {
    _id = ID;
    if(ID && *ID) //do not register if ID is empty
    {
      if(IsDefault || Map().empty())
        Default() = this;
      Map()[ID]=this;
      //Ensure YourBaseClass is registered in OBPlugin so it can be accessed from the commandline
      PluginMap()[TypeID()] =this;
    }
  }
   
  ///Returns the sub-type associated with the ID, or the default subtype if ID NULL or empty.
  static YourBaseClass* FindType(const char* ID)
  {
    if(!ID || *ID==0)
      return Default();
    return static_cast<YourBaseClass*>(BaseFindType(Map(),ID));
  }
\endcode
*/

} // end namespce

#endif
