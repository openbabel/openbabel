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

#include "mol.h"
#include <iostream>

using namespace std;
namespace OpenBabel {

OBBase* OBMol::DoTransformations(const char* Opts)
{
	// Perform any requested transformations
	// on a OBMol
	//The input string has option letters, some of which may be
	//followed by a quoted string, e.g. dcs"CN"v"CO"
	//For normal(non-filter) transforms:
	// returns a pointer to the OBMol (this) if ok or NULL if not.
	//For filters returns a pointer to the OBMol (this) if there is a  match,
	//and NULL when not and in addition the OBMol object is deleted NULL.

	//This is now a virtual function. The OBBase version just returns the OBMol pointer.
 	//This is declared in mol.h 

	//The filter options, s and v allow a obgrep facility. They can be used together.
	
	//Parse GeneralOptions
	int ret=1;
	bool invert=false, smatch=true, vmatch=true;

	const char* p = Opts; //GetGeneralOptions();
	while(*p)
	{
		switch(*(p++))
		{
		case '\"' : //ignore quoted charcters
			p=strchr(p++,'\"');
			if(p++==NULL)
			{
				cerr << "Missing \" in options" <<endl;
				return NULL;
			}
			break;
		case 'd':
			ret=DeleteHydrogens();
			break;
		case 'h':
			ret=AddHydrogens(false, false);
			break;
		case 'c':
			Center(); //has void return
			break;
		case 'p':
			ret=AddHydrogens(false, true);
			break;
		case 'v':
			invert=true;
		case 's':
			//match quoted SMARTS string which follows
			OBSmartsPattern sp;
			string smarts(p+1); //after "
			sp.Init(smarts.substr(0,smarts.find('\"')));
			bool match = sp.Match(*this); //(*pmol) ;
			if(invert)vmatch=!match;
			else smatch=match;
			invert=false;
		}
	}
	if(!smatch || !vmatch)
	{
		//filter failed delete OBMol and return NULL
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
 -p Add Hydrogens appropriate for pH (use transforms in phmodel.txt)\n \
 -c Center Coordinates\n \
 -s\"smarts\" Convert only molecules matching SMARTS:\n \
 -v\"smarts\" Convert only molecules NOT matching SMARTS:\n\n" ;
}

} //namespace OpenBabel
