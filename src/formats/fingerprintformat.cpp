/**********************************************************************
Copyright (C) 2005 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include "mol.h"
#include "obconversion.h"
#include "obmolecformat.h"
#include "fingerprint.h"

using namespace std;
namespace OpenBabel {

class FingerprintFormat : public OBMoleculeFormat
{
public:
	//Register this format type ID
	FingerprintFormat() {OBConversion::RegisterFormat("fpt",this);}

	virtual const char* Description() //required
	{ return
"Fingerprint format\n \
See Fabien Fontaine's source code\n \
";
	};

	virtual unsigned int Flags(){return NOTREADABLE;};

	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
};

////////////////////////////////////////////////////
//Make an instance of the format class
FingerprintFormat theFingerprintFormat;
////////////////////////////////////////////////////////////////

bool FingerprintFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;

	fingerprint fpt(mol.GetTitle());
  fpt.HashMol(mol);
  fpt.printFingerprint(ofs);

	return true;
}

} //namespace OpenBabel
