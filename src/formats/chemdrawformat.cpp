/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

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

class ChemdrawFormat : public OBFormat
{
public:
	//Register this format type ID
	ChemdrawFormat() {OBConversion::RegisterFormat("CHEMDRAW",this);}

	virtual const char* Description() //required
	{ return
"ChemDraw format \n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://??";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return NOTREADABLE;};

//*** This section identical for most OBMol conversions ***
	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

////////////////////////////////////////////////////
	/// The "Convert" interface functions

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
//***

//Make an instance of the format class
ChemdrawFormat theChemdrawFormat;

////////////////////////////////////////////////////////////////

bool ChemdrawFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  char buffer[BUFF_SIZE];

  ofs << mol.GetTitle() << endl;
  sprintf(buffer," %d %d",mol.NumAtoms(),mol.NumBonds());
  ofs << buffer << endl;
  
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;

  for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    sprintf(buffer," %9.4f %9.4f    0.0000 %-1s",
	    atom->x(),
	    atom->y(),
	    etab.GetSymbol(atom->GetAtomicNum()));
    ofs << buffer << endl;
  }

  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;

  for(bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
  {
    sprintf(buffer,"%3d%3d%3d%3d",
	    bond->GetBeginAtomIdx(),
	    bond->GetEndAtomIdx(),
	    bond->GetBO(),1);
    ofs << buffer << endl;
  }
  return(true);
}

} //namespace OpenBabel
