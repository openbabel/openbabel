/**********************************************************************
transform.cpp - Perform command-line requested transformations

Copyright (C) 2004-2005 by Chris Morley

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
#include <sstream>
#include <openbabel/mol.h>
#include <openbabel/generic.h>
#include <openbabel/oberror.h>
#include <openbabel/descriptor.h>
#include <openbabel/op.h>
#include <openbabel/parsmart.h>

#include <cstdlib>

using namespace std;
namespace OpenBabel
{
  class OBConversion; //used only as a pointer

  OBBase* OBMol::DoTransformations(const std::map<std::string, std::string>* pOptions, OBConversion* pConv)
  {
    // Perform any requested transformations
    // on a OBMol
    //The input map has option letters or name as the key and
    //any associated text as the value.
    //For normal(non-filter) transforms:
    // returns a pointer to the OBMol (this) if ok or NULL if not.
    //For filters returns a pointer to the OBMol (this) if there is a  match,
    //and NULL when not and in addition the OBMol object is deleted NULL.

    //This is now a virtual function. The OBBase version just returns the OBMol pointer.
    //This is declared in mol.h

    //The filter options, s and v allow a obgrep facility.
    //Used together they must both be true to allow a molecule through.

    //Parse GeneralOptions
    if(pOptions->empty())
      return this;

    // DoOps calls Do() for each of the plugin options in the map
    // It normally returns true, even if there are no options but
    // can return false if one of the options decides that the
    // molecule should not be output. If it is a filtering op, it
    // should delete the molecule itself (unlike the -s, --filter options,
    // which delete it in this function).
    if(!OBOp::DoOps(this, pOptions, pConv))
      return (OBBase *)NULL;

    bool ret=true;

    map<string,string>::const_iterator itr, itr2;

    if(pOptions->find("b")!=pOptions->end())
      ConvertDativeBonds();
    if(pOptions->find("B")!=pOptions->end())
      MakeDativeBonds();

    if(pOptions->find("d")!=pOptions->end())
      if(!DeleteHydrogens())
        ret=false;

    if(pOptions->find("h")!=pOptions->end())
      if(!AddHydrogens(false, false))
        ret=false;

    if(pOptions->find("r")!=pOptions->end()) {
      StripSalts();
      ret = true;
    }

    itr = pOptions->find("p");
    if(itr!=pOptions->end()) {
      if(pOptions->find("h")!=pOptions->end()){
        stringstream errorMsg;
        errorMsg << "Both -p and -h options are set. "
                 << "All implicit hydrogens (-h) will be added without considering pH."
                 << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      }
      else {
        double pH = strtod(itr->second.c_str(), 0);
        if(!AddHydrogens(false, true, pH))
          ret=false;
      }
    }

    if(pOptions->find("c")!=pOptions->end())
      Center();

    itr = pOptions->find("title"); //Replaces title
    if(itr!=pOptions->end())
      SetTitle(itr->second.c_str());

    itr = pOptions->find("addtotitle"); //Appends text to title
    if(itr!=pOptions->end())
      {
        string title(GetTitle());
        title += itr->second;
        SetTitle(title.c_str());
      }

    //Add an extra property to the molecule.
    //--property attribute value   (no spaces in attribute)
    //or
    // --property attribute1=value1; multi word attribute2 = multi word value2;
    //For both forms, no ';' in either attribute or value.
    //For second form, no '=' in attribute or value
    itr = pOptions->find("property");
    if(itr!=pOptions->end()) {
        string txt(itr->second);

        vector<string> vec;
        vector<string>::iterator it;
        tokenize(vec, txt,";");
        for(it=vec.begin();it!=vec.end();++it) {
          string attr, val;
          string::size_type pos = it->find('=');
          if(pos==string::npos) {
            //form with space
            pos = it->find(' ');
            if(pos==string::npos) {
                obErrorLog.ThrowError(__FUNCTION__, "Missing property value", obError);
                ret=false;
                break;
            }
          }
          val  = it->substr(pos+1);
          Trim(val);
          attr = it->erase(pos);
          Trim(attr);

          //Update value if it already exists
          OBPairData* dp = dynamic_cast<OBPairData*>(GetData(attr));
          if(dp) {
            dp->SetValue(val);
            dp->SetOrigin(userInput);
          }
          else {
            // Pair did not exist; make new one
            dp = new OBPairData;
            dp->SetAttribute(attr);
            dp->SetValue(val);
            dp->SetOrigin(userInput);
            SetData(dp);
          }
        }
      }

    itr = pOptions->find("add");  //adds new properties from descriptors in list
    if(itr!=pOptions->end())
      OBDescriptor::AddProperties(this, itr->second);

    itr = pOptions->find("delete"); //deletes the specified properties
    if(itr!=pOptions->end())
      OBDescriptor::DeleteProperties(this, itr->second);

    itr = pOptions->find("append"); //Appends values of descriptors or properties to title
    if(itr!=pOptions->end())
      {
        string title(GetTitle());
        title += OBDescriptor::GetValues(this, itr->second);
        if(ispunct(title[0]))
          title[0]=' ';//a leading punct char is used only as a separator, not at start
        SetTitle(Trim(title).c_str());
      }



      //Filter using OBDescriptor comparison and (older) SMARTS tests
    //Continue only if previous test was true.
    bool fmatch = true;
    itr = pOptions->find("filter");
    if(itr!=pOptions->end())
      {
        std::istringstream optionText(itr->second);
        fmatch = OBDescriptor::FilterCompare(this, optionText, false);
      }

    if(fmatch)
      {
        itr = pOptions->find("v");
        if(itr!=pOptions->end() && !itr->second.empty())
          {
            //inverse match quoted SMARTS string which follows
            OBSmartsPattern sp;
            sp.Init(itr->second);
            fmatch = !sp.Match(*this); //(*pmol) ;
          }
      }
    if(fmatch)
    {
      itr = pOptions->find("s");
      if(itr!=pOptions->end() && !itr->second.empty())
        {
          //SMARTS filter
          //If exactmatch option set (probably in fastsearchformat) the
          //number of atoms in the pattern (passed as a string in the option text)
          //has to be the same as in the molecule.
          itr2 = pOptions->find("exactmatch");
          if(itr2!=pOptions->end() && (int)NumHvyAtoms()!=atoi(itr2->second.c_str()))
            fmatch=false;
          else
            {
              //match quoted SMARTS string which follows
              OBSmartsPattern sp;
                sp.Init(itr->second.c_str());
                fmatch = sp.Match(*this);
            }
        }
    }

    if(!fmatch)
      {
        //filter failed: delete OBMol and return NULL
        delete this;
        return NULL;
      }
    else
      {
        if(ret==false)
          {
            obErrorLog.ThrowError(__FUNCTION__, "Error executing an option", obError);
            delete this; //added 9March2006
            return NULL;
          }
        else
          return this;
      }
  }

  ///////////////////////////////////////////////////
  const char* OBMol::ClassDescription()
  {
    static string ret;
    ret = "For conversions of molecules\n"
"Additional options :\n"
"-d Delete hydrogens (make implicit)\n"
"-h Add hydrogens (make explicit)\n"
"-p <pH> Add hydrogens appropriate for this pH\n"
"-b Convert dative bonds e.g.-[N+]([O-])=O to -N(=O)=O\n"
"-B Make dative bonds e.g.-[N+]([O-])=O from -N(=O)=O\n"
"-r Remove all but the largest contiguous fragment\n"
"-c Center Coordinates\n"
"-C Combine mols in first file with others by name\n"
"--filter <filterstring> Filter: convert only when tests are true:\n"
"--add <list> Add properties from descriptors\n"
"--delete <list> Delete properties in list\n"
"--append <list> Append properties or descriptors in list to title:\n"
"-s\"smarts\" Convert only if match SMARTS or mols in file:\n"
"-v\"smarts\" Convert only if NO match to SMARTS or mols in file(not displayed in GUI)\n"
"--join Join all input molecules into a single output molecule\n"
"--separate Output disconnected fragments separately\n"
"--property <attrib> <value> add or replace a property (SDF)\n"
"--title <title> Add or replace molecule title\n"
"--addtotitle <text> Append text to title\n"
"--writeconformers Output multiple conformers separately\n"
"--addoutindex Append output index to title\n" ;

    //Append lines from OBOp plugins that work with OBMol
    OBMol dummymol; //just needed to carry class type information; messy!
    ret += OBOp::OpOptions(&dummymol);

      return ret.c_str();
  }

} //namespace OpenBabel

//! \file transform.cpp
//! \brief Perform command-line requested transformations for OBMol
//!  and SMARTS filtering
