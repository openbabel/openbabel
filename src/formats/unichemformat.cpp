/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
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

class UniChemFormat : public OBFormat
{
public:
	//Register this format type ID
	UniChemFormat() {OBConversion::RegisterFormat("UNICHEM",this);}

	virtual const char* Description() //required
	{ return
"UniChem format\n \
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
UniChemFormat theUniChemFormat;

/////////////////////////////////////////////////////////////////
bool UniChemFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

  int i;
  int natoms;
  char buffer[BUFF_SIZE];

  ifs.getline(buffer,BUFF_SIZE);
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%d", &natoms);
  if (!natoms) return(false);

  mol.ReserveAtoms(natoms);

  string str;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;

  for (i = 1; i <= natoms; i ++)
  {
    if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
    tokenize(vs,buffer);
    if (vs.size() != 4) return(false);
    atom = mol.NewAtom();
    x = atof((char*)vs[1].c_str());
    y = atof((char*)vs[2].c_str());
    z = atof((char*)vs[3].c_str());
    atom->SetVector(x,y,z); //set coordinates

    //set atomic number
    atom->SetAtomicNum(atoi((char*)vs[0].c_str()));
  }
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  mol.SetTitle(title);
  return(true);
}

////////////////////////////////////////////////////////////////

bool UniChemFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  char buffer[BUFF_SIZE];

  ofs << mol.GetTitle() << endl;
  ofs << mol.NumAtoms() << endl;

  OBAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%3d%15.5f%15.5f%15.5f",
	    atom->GetAtomicNum(),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer << endl;
  }

  return(true);
}

} //namespace OpenBabel
