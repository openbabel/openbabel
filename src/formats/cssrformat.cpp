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

class CSSRFormat : public OBFormat
{
public:
	//Register this format type ID
	CSSRFormat() {OBConversion::RegisterFormat("CSSR",this);}

	virtual const char* Description() //required
	{ return
"CSSR Format\n \
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
CSSRFormat theCSSRFormat;

////////////////////////////////////////////////////////////////

// \todo Needs update for unit cell information
bool CSSRFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  char buffer[BUFF_SIZE];

  sprintf(buffer,
	  " REFERENCE STRUCTURE = 00000   A,B,C =  %6.3f  %6.3f  %6.3f",
	  1.0,1.0,1.0); 
  ofs << buffer << endl;
  sprintf(buffer,
	  "   ALPHA,BETA,GAMMA =  90.000  90.000  90.000    SPGR =    P1");
  ofs << buffer << endl;
  sprintf(buffer,"%4d\n",mol.NumAtoms());
  ofs << buffer << endl;

  OBAtom *atom,*nbr;
  vector<OBNodeBase*>::iterator i;
  vector<OBEdgeBase*>::iterator j;
  vector<int> vtmp(106,0);

  for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    //assign_pdb_number(pdb_types,atom->GetIdx());
    vtmp[atom->GetAtomicNum()]++;
    sprintf(buffer," %3d%2s%-3d  %8.4f  %8.4f  %8.4f ",
	    atom->GetIdx(),
	    etab.GetSymbol(atom->GetAtomicNum()),
	    vtmp[atom->GetAtomicNum()],
	    atom->x(),
	    atom->y(),
	    atom->z());
    ofs << buffer;
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      {
	sprintf(buffer,"%4d",nbr->GetIdx());
	ofs << buffer;
      }
    ofs << endl;
  }

  return(true);
}

} //namespace OpenBabel
