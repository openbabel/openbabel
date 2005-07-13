/**********************************************************************
APIInterface.cpp - Pseudo-format to transfer data from user interface to API

Copyright (C) 2005 by Chris Morley

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
#include "oberror.h"
#include "obconversion.h"

using namespace std;
namespace OpenBabel
{

class OBAPIInterface : public OBFormat
{
public:
  OBAPIInterface()
	{		OBConversion::RegisterFormat("obapi",this); }

	const char* Description(){return "Interface to OBAPI internals";}

	unsigned int Flags(){ return (NOTWRITABLE | NOTREADABLE);}

	bool WriteMolecule(OBBase* , OBConversion* pConv)
	{
		string txt = pConv->IsOption("w");
		if(!txt.empty())
		{
			stringstream ss(txt);
			int ilevel=-1;
			ss >> ilevel;	
			if(ilevel>=0)
				obErrorLog.SetOutputLevel((obMessageLevel)ilevel);
		}
		return true;
	}
};

OBAPIInterface theOBAPIInterface;

} //namespace