/**********************************************************************
transform.cpp - Perform command-line requested transformations

Copyright (C) 2004-2005 by Chris Morley
 
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
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>

using namespace std;
namespace OpenBabel
{

  OBBase* OBMol::DoTransformations(const std::map<std::string, std::string>* pOptions)
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

    bool ret=true;
    bool smatch=true, vmatch=true;

    map<string,string>::const_iterator itr;

    if(pOptions->find("b")!=pOptions->end())
      if(!ConvertDativeBonds())
        ret=false;

    if(pOptions->find("d")!=pOptions->end())
      if(!DeleteHydrogens())
        ret=false;

    if(pOptions->find("h")!=pOptions->end())
      if(!AddHydrogens(false, false))
        ret=false;

    if(pOptions->find("p")!=pOptions->end())
      if(!AddHydrogens(false, true))
        ret=false;

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

    itr = pOptions->find("addformula"); //Appends tab + formula to title
    if(itr!=pOptions->end())
      {
        string title(GetTitle());
        title += '\t' + GetSpacedFormula(1,"");//actually unspaced
        SetTitle(title.c_str());
      }

    //Add an extra property to the molecule.
    //Parameter has atrribute and value separated by a space
    itr = pOptions->find("property");
    if(itr!=pOptions->end())
      {
        string txt(itr->second);
        string::size_type pos = txt.find(' ');
        if(pos==string::npos)
          {
            obErrorLog.ThrowError(__FUNCTION__, "Missing property value", obError);
            ret=false;
          }
        else
          {
            string attr(txt.substr(0,pos)), val(txt.substr(pos+1));
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

    itr = pOptions->find("v");
    if(itr!=pOptions->end())
      {
        //inverse match quoted SMARTS string which follows
        OBSmartsPattern sp;
        sp.Init(itr->second);
        vmatch = !sp.Match(*this); //(*pmol) ;
      }

    itr = pOptions->find("s");
    if(itr!=pOptions->end())
      {
        //match quoted SMARTS string which follows
        OBSmartsPattern sp;
        sp.Init(itr->second.c_str());
        smatch = sp.Match(*this); //(*pmol) ;
      }

    if(!smatch || !vmatch)
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
    return "For conversions of molecules\n \
Additional options :\n \
-d Delete hydrogens (make implicit)\n \
-h Add hydrogens (make explicit)\n \
-p Add Hydrogens appropriate for pH model\n \
-b Convert dative bonds e.g.[N+]([O-])=O to N(=O)=O\n \
-c Center Coordinates\n \
--join Join all input molecules into a single output molecule\n \
--separate Output disconnected fragments separately\n \
-C Combine mols in first file with others having same name\n \
-s\"smarts\" Convert only molecules matching SMARTS:\n \
-v\"smarts\" Convert only molecules NOT matching SMARTS:\n \
--property <attrib> <value> add or replace a property (SDF)\n \
--title <title> Add or replace molecule title\n \
--addtotitle <text> Append to title\n \
--addformula Append formula to title\n\n" ;
  }

} //namespace OpenBabel

//! \file transform.cpp
//! \brief Perform command-line requested transformations for OBMol
//!  and SMARTS filtering
