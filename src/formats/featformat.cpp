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

class FEATFormat : public OBFormat
{
public:
	//Register this format type ID
	FEATFormat() {OBConversion::RegisterFormat("FEAT",this);}

	virtual const char* Description() //required
	{ return
"FEAT format\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://??";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return READONEONLY;};

//*** This section identical for most OBMol conversions ***
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
//***

//Make an instance of the format class
FEATFormat theFEATFormat;

/////////////////////////////////////////////////////////////////
bool FEATFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

  char buffer[BUFF_SIZE];
  int i,natoms;

  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%d",&natoms);

  mol.ReserveAtoms(natoms);
  mol.BeginModify();

  if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
  mol.SetTitle(buffer);

  double x,y,z;
  char type[20];
  OBAtom *atom;
  for (i = 0; i < natoms;i++)
  {
    if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
    sscanf(buffer,"%s %lf %lf %lf",
	   type,
	   &x,
	   &y,
	   &z);
    CleanAtomType(type);
    atom = mol.NewAtom();
    atom->SetVector(x,y,z);
    atom->SetAtomicNum(etab.GetAtomicNum(type));
  }

  mol.EndModify();
  return(true);
}

////////////////////////////////////////////////////////////////

bool FEATFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  char buffer[BUFF_SIZE];
  
  ofs << mol.NumAtoms() << endl;
  ofs << mol.GetTitle() << endl;

  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  {
    sprintf(buffer,"%-3s %8.5f  %8.5f  %8.5f ",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->x(),
	    atom->y(),
	    atom->z());
    ofs << buffer << endl;
  }

  return(true);
}

} //namespace OpenBabel
