/**********************************************************************
Copyright (C) 2002-2003 Peter Murray-Rust.
Some portions Copyright (c) 2003 by Geoffrey R. Hutchison
Some Portions Copyright (C) 2004 by Chris Morley

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
#include "obconversion.h"

using namespace std;
namespace OpenBabel {
//The original routines
extern bool ReadCML(istream &ifs,OBMol &mol, const char *title);
extern bool WriteCML(ostream &ofs,OBMol &mol,const char *dim,const char* xmlOptions);

class CMLFormat : public OBFormat
{
public:
	//Register this format type ID
	CMLFormat() {OBConversion::RegisterFormat("CML",this);}

	virtual const char* Description()
	{ return
"Chemical Markup Language\n \
XML format\n \
Additional command line option for CML: -x[flags] (e.g. -x1ac)\n \
 1  output CML V1.0 (default) or\n \
 2  output CML V2.0 (Schema)\n \
 a  output array format for atoms and bonds (default <atom>)\n \
 p  prettyprint output (default no indent)\n \
 n  output namespace (default no namespace)\n \
 c  use 'cml' as output namespace prefix (default no namspace)\n \
 d  output DOCTYPE (default none)\n \
 g  debug output\n \
";
	};
	virtual const char* SpecificationURL(){return
		"http://wwmm.ch.cam.ac.uk/moin/ChemicalMarkupLanguage";};

	virtual unsigned int Flags(){return READONEONLY;};

	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

////////////////////////////////////////////////////
	/// The "Convert" interface functions
	virtual bool ReadChemObject(OBConversion* pConv)
	{
		OBMol* pmol = new OBMol;
		bool ret=ReadMolecule(pmol,pConv);
		if(ret) //Do transformation and return molecule
			pConv->AddChemObject(pmol->DoTransformations(pConv->GetGeneralOptions()));
		else
			pConv->AddChemObject(NULL);
		return ret;
	};
	
	virtual bool WriteChemObject(OBConversion* pConv)
	{
		//Retrieve the target OBMol
		OBBase* pOb = pConv->GetChemObject();
		OBMol* pmol = dynamic_cast<OBMol*> (pOb);
		bool ret=false;
		if(pmol)
			ret=WriteMolecule(pmol,pConv);
		delete pOb; 
		return ret;
	};
};

//Make an instance of the format class
CMLFormat theCMLFormat;

/////////////////////////////////////////////////////////////////
bool CMLFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;

	//Just call the old routine
	 return ReadCML(ifs, mol, pConv->GetTitle());
}

////////////////////////////////////////////////////////////////

bool CMLFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();
//****

	//If more than one molecule to be output, write <cml> at start and </cml> at end.
	if((pConv->GetOutputIndex()==1) && !pConv->IsLast())
		ofs << "<cml>" << endl;
	
	//Just call the old routine
	bool ret = WriteCML(ofs, mol, dimension, pConv->GetOptions());

	if((pConv->GetOutputIndex()>1) && pConv->IsLast())
		ofs << "</cml>" << endl;

	return ret;
}

}
