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

#include "mol.h"
#include <iostream>

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

		int ret=1;
    bool smatch=true, vmatch=true;

		map<string,string>::const_iterator itr;

		if(pOptions->find("b")!=pOptions->end())
			ret=ConvertDativeBonds();

		if(pOptions->find("d")!=pOptions->end())
			ret=DeleteHydrogens();

		if(pOptions->find("h")!=pOptions->end())
      ret=AddHydrogens(false, false);

		if(pOptions->find("p")!=pOptions->end())
      ret=AddHydrogens(false, true);

		if(pOptions->find("c")!=pOptions->end())
		{
			Center();
			ret=1;
		}

		itr = pOptions->find("addtotitle"); //Appends text to title
		if(itr!=pOptions->end())
    {
			string title(GetTitle());
			title += itr->second;
			SetTitle(title.c_str());
			ret=1;
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
        return ret ? this : NULL;
}

///////////////////////////////////////////////////
const char* OBMol::ClassDescription()
{
    return "For conversions of molecules\n \
  Additional options :\n \
  -d Delete Hydrogens\n \
  -h Add Hydrogens\n \
  -p Add Hydrogens appropriate for pH\n \
  -b Convert dative bonds e.g.[N+]([O-])=O to N(=O)=O\n \
  -c Center Coordinates\n \
  -j Join all input molecules into a single output molecule\n \
  -s\"smarts\" Convert only molecules matching SMARTS:\n \
  -v\"smarts\" Convert only molecules NOT matching SMARTS:\n\n" ;
}

} //namespace OpenBabel

//! \file transform.cpp
//! \brief Perform command-line requested transformations for OBMol
//!  and SMARTS filtering
