/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2004 by Geoffrey R. Hutchison
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

class GAUSSFormat : public OBFormat
{
public:
	//Register this format type ID
	GAUSSFormat() {OBConversion::RegisterFormat("GAUSS",this);}

	virtual const char* Description() //required
	{ return
"Gaussian98 Cartesian\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://www.gaussian.com/g_ur/m_input.htm";};

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
GAUSSFormat theGAUSSFormat;


////////////////////////////////////////////////////////////////

bool GAUSSFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  int charge = 0;
  unsigned int multiplicity = 0;
  char buffer[BUFF_SIZE];
  
  ofs << "%cmem=20000000" << endl << '\045';
  ofs << "#Put Keywords Here, check Charge and Multiplicity" << endl << endl;
  ofs << "XX " << mol.GetTitle() << endl << endl;

  OBAtom *atom;
  string str,str1;
  sprintf(buffer,"  %d  %d", mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity());
  ofs << buffer << endl;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%-3s      %10.5f      %10.5f      %10.5f ",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX(), atom->GetY(), atom->GetZ());
    ofs << buffer << endl;
  }
  // file should end with a blank line
  ofs << endl;
  return(true);
}

} //namespace OpenBabel
