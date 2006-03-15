/**********************************************************************
obmolecformat.h - Subclass of OBFormat for conversion of OBMol.

Copyright (C) 2005 Chris Morley

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

#ifndef OB_MOLECULEFORMAT_H
#define OB_MOLECULEFORMAT_H

#include "mol.h"
#include "obconversion.h"

namespace OpenBabel {

//! \brief An OBFormat convenience subclass for conversion to/from OBMol data
//!
//! An OBFormat which converts to and/or from OBMol can derive from this class
//! to save duplicating the ReadChemObject() and/or WriteChemObject() methods.
//! Derive directly from OBFormat if the object converted is not OBMol or 
//! if interaction with the framework is required during the execution 
//! of ReadMolecule() or WriteMolecule(), as for example in CMLFormat
class OBMoleculeFormat : public OBFormat
{
private:
	static std::map<std::string, OBMol*> IMols;
	static OBMol* _jmol; //Accumulates molecules with the -j option

public:

	OBMoleculeFormat()
	{
		OBConversion::RegisterOptionParam("b", this, 0, OBConversion::INOPTIONS);
		OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);
		OBConversion::RegisterOptionParam("title", this, 1,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("addtotitle", this, 1,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("property", this, 2, OBConversion::GENOPTIONS);
		//The follow are OBMol options, which should not be in OBConversion.
		//But here isn't entirely appropriate either, since could have
		//OBMol formats loaded but none of them derived from this class.
		//However, this possibility is remote.
		OBConversion::RegisterOptionParam("s", NULL, 1,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("v", NULL, 1,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("h", NULL, 0,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("d", NULL, 0,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("b", NULL, 0,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("c", NULL, 0,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("p", NULL, 0,OBConversion::GENOPTIONS); 
		OBConversion::RegisterOptionParam("t", NULL, 0,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("j", NULL, 0,OBConversion::GENOPTIONS);
		OBConversion::RegisterOptionParam("C", NULL, 0,OBConversion::GENOPTIONS);
	};

	//Static routines which can be called from elsewhere
	static bool ReadChemObjectImpl(OBConversion* pConv, OBFormat*);
	static bool WriteChemObjectImpl(OBConversion* pConv, OBFormat*);

	/// The "Convert" interface functions
	virtual bool ReadChemObject(OBConversion* pConv)
	{ return ReadChemObjectImpl(pConv, this);}
		
	virtual bool WriteChemObject(OBConversion* pConv)
	{ return WriteChemObjectImpl(pConv, this);}
	
	/// Routines to handle the -C option for combining data from several OBMols
	static bool DeferMolOutput(OBMol* pmol, OBConversion* pConv, OBFormat* pF);
	static bool OutputDeferredMols(OBConversion* pConv);
	static bool DeleteDeferredMols();
	static OBMol* MakeCombinedMolecule(OBMol* pFirst, OBMol* pSecond);

	const std::type_info& GetType()
	{
		return typeid(OBMol*);
	}
//////////////////////////////////////////////////////////////

};

}
#endif //OB_MOLECULEFORMAT_H

//! \file obmolecformat.h
//! \brief Subclass of OBFormat for conversion of OBMol.
