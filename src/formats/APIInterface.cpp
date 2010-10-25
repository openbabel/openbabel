/**********************************************************************
APIInterface.cpp - Pseudo-format to transfer data from user interface to API

Copyright (C) 2005 by Chris Morley

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
#include <openbabel/obconversion.h>

using namespace std;
namespace OpenBabel
{

class OBAPIInterface : public OBFormat
{
public:
  OBAPIInterface()
	{
		OBConversion::RegisterFormat("obapi",this);
		OBConversion::RegisterOptionParam("-errorlevel", this, 1, OBConversion::GENOPTIONS);
	}

	const char* Description(){
    return
    "Interface to OBAPI internals\n"
    "API options, e.g. ---errorlevel 2\n"
    " errorlevel # min warning level displayed\n\n";
  }

	unsigned int Flags(){ return (NOTWRITABLE | NOTREADABLE);}

	bool WriteMolecule(OBBase* , OBConversion* pConv)
	{
		const char* txt = pConv->IsOption("errorlevel",OBConversion::GENOPTIONS);
		if(txt)
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
