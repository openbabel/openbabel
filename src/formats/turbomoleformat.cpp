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
#include "obconversion.h"

using namespace std;
namespace OpenBabel {

class TurbomoleFormat : public OBFormat
{
public:
	//Register this format type ID
	TurbomoleFormat() {OBConversion::RegisterFormat("TURBOMOLE",this);}

	virtual const char* Description() //required
	{ return
"TurboMole Coordinate file\n \
Incomplete - no bonds\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://www.turbomole.com/";}; //optional


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
TurbomoleFormat theTurbomoleFormat;

/////////////////////////////////////////////////////////////////
bool TurbomoleFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;

  char buff[BUFF_SIZE];
	do
	{
		ifs.getline(buff,BUFF_SIZE);
	}while(strncmp(buff,"$coord",6));

  mol.BeginModify();
	OBAtom atom;
	while(!(!ifs))
	{
		ifs.getline(buff,BUFF_SIZE);
		if(*buff=='$') break;
		if(*buff=='#') continue;
		float x,y,z;
		char atomtype[5];
		if(sscanf(buff,"%f %f %f %s",&x,&y,&z,atomtype)!=4)
			return false;

		atom.SetVector(x, y, z);
    atom.SetAtomicNum(etab.GetAtomicNum(atomtype));
    atom.SetType(atomtype);

		if(!mol.AddAtom(atom))
			return false;
		atom.Clear();
	}
	while(!(!ifs) && strncmp(buff,"$end",4)) ifs.getline(buff,BUFF_SIZE);
		
  mol.EndModify();
	return true;
}

////////////////////////////////////////////////////////////////

bool TurbomoleFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;

  ofs << "$coord" <<endl;

	char buff[BUFF_SIZE];
	OBAtom *atom;
	vector<OBNodeBase*>::iterator i;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
	{
		sprintf(buff,"%20.14f  %20.14f  %20.14f %6s",
        atom->GetX(),
        atom->GetY(),
        atom->GetZ(),
        strlwr(etab.GetSymbol(atom->GetAtomicNum())));
		ofs << buff << endl;
	}
	ofs << "$end" << endl;

	return true;
}

} //namespace OpenBabel
