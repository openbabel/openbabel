/**********************************************************************
Copyright (C) 2000-2003 by Geoffrey Hutchison
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

#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989


class DMolFormat : public OBFormat
{
public:
	//Register this format type ID
	DMolFormat() {OBConversion::RegisterFormat("DMOL",this);}

	virtual const char* Description() //required
	{ return
"DMol format\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://wwmm.ch.cam.ac.uk/moin/ChemicalMarkupLanguage";}; //optional

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
DMolFormat theDMolFormat;

/////////////////////////////////////////////////////////////////
bool DMolFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

  char buffer[BUFF_SIZE];
  string str;
  double x,y,z;
  OBAtom *atom;
  vector3 v1,v2,v3;
  vector<string> vs;

  ifs.getline(buffer,BUFF_SIZE);
  while (strstr(buffer,"$coordinates") == NULL &&
	 strstr(buffer,"$cell vectors") == NULL)
    ifs.getline(buffer,BUFF_SIZE);

  if (strstr(buffer,"$cell vectors") != NULL)
    {
      ifs.getline(buffer,BUFF_SIZE);
      tokenize(vs,buffer); // we really need to check that it's 3 entries only
      x = atof((char*)vs[0].c_str()) * BOHR_TO_ANGSTROM;
      y = atof((char*)vs[1].c_str()) * BOHR_TO_ANGSTROM;
      z = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
      v1.Set(x,y,z);
      ifs.getline(buffer,BUFF_SIZE);
      tokenize(vs,buffer);
      x = atof((char*)vs[0].c_str()) * BOHR_TO_ANGSTROM;
      y = atof((char*)vs[1].c_str()) * BOHR_TO_ANGSTROM;
      z = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
      v2.Set(x,y,z);
      ifs.getline(buffer,BUFF_SIZE);
      tokenize(vs,buffer);
      x = atof((char*)vs[0].c_str()) * BOHR_TO_ANGSTROM;
      y = atof((char*)vs[1].c_str()) * BOHR_TO_ANGSTROM;
      z = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
      v3.Set(x,y,z);

      OBUnitCell *uc = new OBUnitCell;
      uc->SetData(v1,v2,v3);
      mol.SetData(uc);

      ifs.getline(buffer,BUFF_SIZE); // next line
    }

  while (strstr(buffer,"$end") == NULL)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) break;
      tokenize(vs,buffer);
      if (vs.size() != 4) break;
      atom = mol.NewAtom();
      //set atomic number
      atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
      x = atof((char*)vs[1].c_str()) * BOHR_TO_ANGSTROM;
      y = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
      z = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
      atom->SetVector(x,y,z); //set coordinates
    }
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  mol.SetTitle(title);
  return(true);
}

////////////////////////////////////////////////////////////////

bool DMolFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  char buffer[BUFF_SIZE];

  if (mol.HasData(obUnitCell))
    {
      OBUnitCell *uc = (OBUnitCell*)mol.GetData(obUnitCell);
      vector<vector3> v = uc->GetCellVectors();
      vector3 v1;

      ofs << "$cell vectors" << endl;
      v1 = v[0] * ANGSTROM_TO_BOHR;
      sprintf(buffer,"%-3s% 27.14f% 20.14f% 20.14f","", v1.x(), v1.y(), v1.z());
      ofs << buffer << endl;
      v1 = v[1] * ANGSTROM_TO_BOHR;
      sprintf(buffer,"%-3s% 27.14f% 20.14f% 20.14f","", v1.x(), v1.y(), v1.z());
      ofs << buffer << endl;
      v1 = v[2] * ANGSTROM_TO_BOHR;
      sprintf(buffer,"%-3s% 27.14f% 20.14f% 20.14f","", v1.x(), v1.y(), v1.z());
      ofs << buffer << endl;
    }
  
  ofs << "$coordinates" << endl;

  OBAtom *atom;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%-3s% 27.14f% 20.14f% 20.14f",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX() * ANGSTROM_TO_BOHR,
	    atom->GetY() * ANGSTROM_TO_BOHR,
	    atom->GetZ() * ANGSTROM_TO_BOHR);
    ofs << buffer << endl;
  }

  ofs << "$end" << endl;

  return(true);
}

} //namespace OpenBabel
