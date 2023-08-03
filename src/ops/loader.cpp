/**********************************************************************
loader.cpp Plugin type "loaders"
  and derived class OBDefine to make plugin instances from a text file
Copyright (C) 2008 Chris Morley

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
#include <iosfwd>
#include <string>
#include <vector>
#include <cstdlib>
#include <openbabel/oberror.h>
#include <openbabel/tokenst.h>
#include <openbabel/plugin.h>
#include <openbabel/locale.h>

using namespace std;
namespace OpenBabel
{
//A top-level plugin (with only basic capabilities)
class OBLoader : public OBPlugin
{
  MAKE_PLUGIN(OBLoader)
public:
  const char* TypeID(){return "loaders";};
};

#if defined(__CYGWIN__) || defined(__MINGW32__)
// macro to implement static OBPlugin::PluginMapType& Map()
PLUGIN_CPP_FILE(OBLoader)
#endif

//*********************************************************
///Class which makes instances of plugin classes from information in text file.
///This allows the commandline and GUI interfaces to be extended without recompiling.
///The class is itself a plugin but needs a short piece of code as a hook in
///OBPlugin::LoadAllPlugins(). This does nothing if the plugin is not loaded.
class OBDefine : public OBLoader
{
public:
  //Constructor for placeholder
  OBDefine() : OBLoader("define", false), _descr("*"){}

  //Main constructor
  OBDefine(const char* ID, const char* filename)
    : OBLoader(ID, false), _filename(filename)
  {
    ifstream ifs;
    bool filefound = !OpenDatafile(ifs, filename).empty();
    if(!ifs)
    {
      if(filefound)
        obErrorLog.ThrowError(__FUNCTION__,string(filename) + " found but could not be opened", obError);
      return;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    string ln;
    while(ifs) //read entries for multiple objects
    {
      //More general version
      vector<string> textlines(1);
      while(ifs && !ReadLine(ifs, textlines[0], true)); //get first non-blank line
      if(!ifs)
        break;
      while(ReadLine(ifs, ln, true )) //get subsequent lines with comments removed
      {
        //Concatenate with next line if ends with "\n". Otherwise each in separate string
        if(ln.size()>=2 && ln.substr(ln.size()-2)=="\\n")
          ln.replace(ln.size()-2, 2,"\n");//replace "\n" by '\n'
        if(textlines.back()[textlines.back().size()-1]=='\n')
          textlines.back() += ln; //concatenate
        else
          textlines.push_back(ln);
      }
      //look up class name in map maintained in OBPlugin
      if(!textlines.empty() && FindDef(textlines[0].c_str()) != nullptr)
      {
        //Save a copy of textlines so that const char* pointers in new instance point
        //to something which has not been deleted
        _text.push_back(textlines);
        // Previously, instances of the plugins were made here:
        //_instances.push_back(pdef->MakeInstance(_text.back()));
        // However, subsequent _text.push_back calls can cause vector reallocations, 
        // which invalidates the const char* pointers in the instances. So instead we
        // populate _text fully before making the plugin instances.
      }
      else
        obErrorLog.ThrowError(__FUNCTION__, "Failed to make an instance " + textlines[0], obError);
      textlines.clear();
    }

    // Iterate through _text and make instances of the plugins.
    // They will be deleted in the destructor.
    vector<vector<string> >::iterator iter;
    for(iter=_text.begin();iter!=_text.end();++iter) {
      OBPlugin* pdef = FindDef((*iter)[0].c_str());
      _instances.push_back(pdef->MakeInstance(*iter));
    }

    // return the locale to the original one
    obLocale.RestoreLocale();
  }


  virtual ~OBDefine()
  {
    std::vector<OBPlugin*>::iterator iter;
    for(iter=_instances.begin();iter!=_instances.end();++iter)
      delete *iter;
  }

  virtual const char* Description(){ return "Makes plugin classes from a datafile"; }

  virtual OBDefine* MakeInstance(const std::vector<std::string>& textlines)
  {
    OBDefine* pdef = new OBDefine(textlines[1].c_str(),textlines[2].c_str());

    //The following line links the new instance to its parent, and eventually
    //the dummy global instance. When the global instance is deleted at the end,
    //all the other instances are too, avoiding (pseudo) memory leaks. Versions
    //for other OBPlugin classes don't do this.
    _instances.push_back(pdef);
    return pdef;
  }

private:
  ///Returns the address of an instance of a plugin class that
  /// has been registered as user definable
  static OBPlugin* FindDef(const char* ID);

  //return true if non-empty line read, false if empty line read or eof or error
  static bool ReadLine(std::istream& ifs, std::string& ln, bool removeComments);

private:
  const char* _filename;
  const char* _descr;
  std::vector<OBPlugin*> _instances;
  std::vector<std::vector<std::string> > _text;
};

//***********************************************************
/* This dummy global instance has no content but is necessary to inform the
system of the existence of OBDefine. It cannot do the work of the plugin and
parse the datafile, because the plugins referred to there may not have been
loaded yet. Another instance with the same ID is made using MakeInstance() in
OBPlugin::LoadAllPlugins() after all the plugins are present.*/

OBDefine placeholderOBDefine;
//************************************************************

//Read line and remove comments
  //Return true if non-empty line read, false if empty line read or eof or error
  bool OBDefine::ReadLine(istream& ifs, string& ln, bool removeComments)
  {
    if(getline(ifs, ln))
    {
      if(removeComments)
      {
        //Remove rest of line after # in first column or # followed by whitespace
        string::size_type pos = ln.find('#');
        if(pos!=string::npos && (pos==0 || isspace(ln[pos+1])))
          ln.erase(pos);
      }
      Trim(ln);
      return !ln.empty();
    }
    return false;
  }

OBPlugin* OBDefine::FindDef(const char* ID)
{
  PluginIterator typeiter, iter;
  for(typeiter=PluginMap().begin(); typeiter!=PluginMap().end();++typeiter)
  {
    PluginMapType map = typeiter->second->GetMap();
    for(iter=map.begin();iter!=map.end();++iter)
    {
      const char* pdescr = iter->second->Description();
      if(!pdescr)
        continue;
      string descr(pdescr);
      //Matches if the ID is before the last occurrence of "definable"
      //in the description on the same line.
      string::size_type pos, pos2;
      pos= descr.rfind("definable");
      if(pos==string::npos)
        continue;
      pos2 = descr.rfind('\n', pos);
      if(pos2!=string::npos && descr.substr(pos2, pos-pos2).find(ID)!=string::npos)
        return iter->second;
    }
  }
  return nullptr;
}

}//namespace

/*
Derived class OBDDefine
global variable with minimal constructor - no other action

In OBConversion get pointer to this class OBPlugin::GetPlugin(), which
has ASCII parameters, and call the vitual OBPlugin function MakeInstance
with plugindefines.txt as parameter to make a working instance of OBDefine.
This removes the dependence of OBConversion on the new OBDefine code.

Definable plugin classes have a line in their descripions
<classname> is definable
So in main constructor of OBDefine the datafile is parsed and the
class names looked for by scanning the descriptions of all plugin
classes for "is definable" and checking the name. This removes any
mention of OBDefine in any of the other plugin classes.
*/
