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

class MOPACFormat : public OBFormat
{
public:
	//Register this format type ID
	MOPACFormat() {OBConversion::RegisterFormat("MOPAC",this);}

	virtual const char* Description() //required
	{ return
"Mopac file\n \
No comments yet";
	};

	virtual unsigned int Flags(){return NOTWRITABLE;};

//*** This section identical for most OBMol conversions ***
	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
//	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv); Is Read Only

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
};
//***

//Make an instance of the format class
MOPACFormat theMOPACFormat;

/////////////////////////////////////////////////////////////////
bool MOPACFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title= pConv->GetTitle();

	char buffer[BUFF_SIZE];
  string str,str1;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;
  vector<double> charges;
  bool hasPartialCharges = false;
  double energy;

  mol.BeginModify();
  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"CARTESIAN COORDINATES") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  mol.BeginModify();
	  ifs.getline(buffer,BUFF_SIZE);	// blank
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);	// blank
	  ifs.getline(buffer,BUFF_SIZE);	  
	  tokenize(vs,buffer);
	  while (vs.size() == 5)
	    {
	      atom = mol.NewAtom();
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));
	      x = atof((char*)vs[2].c_str());
	      y = atof((char*)vs[3].c_str());
	      z = atof((char*)vs[4].c_str());
	      atom->SetVector(x,y,z);

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      else if(strstr(buffer,"FINAL HEAT") != NULL)
	{
	  sscanf(buffer,"%*s%*s%*s%*s%*s%lf",&energy);
	  mol.SetEnergy(energy);
	}
      else if(strstr(buffer,"NET ATOMIC CHARGES") != NULL)
	{
	  hasPartialCharges = true;
	  charges.clear();
	  ifs.getline(buffer,BUFF_SIZE);	// blank
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 4)
	    {
	      atom = mol.GetAtom(atoi(vs[0].c_str()));
	      atom->SetPartialCharge(atof(vs[2].c_str()));
	      charges.push_back(atof(vs[2].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
    }
  mol.EndModify();
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  if (hasPartialCharges)
    {
      mol.SetPartialChargesPerceived();
      for (unsigned int i = 1; i <= mol.NumAtoms(); i++)
	{
	  atom = mol.GetAtom(i);
	  atom->SetPartialCharge(charges[i-1]);
	}
    }
  mol.SetTitle(title);

  return(true);
}

//************************************************************
class MOPACCARTFormat : public OBFormat
{
public:
	//Register this format type ID
	MOPACCARTFormat() {OBConversion::RegisterFormat("MOPACCART",this);}

	virtual const char* Description() //required
	{ return
"MOPACCART file\n \
No comments yet";
	};

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
MOPACCARTFormat theMOPACCARTFormat;

/////////////////////////////////////////////////////////////////
bool MOPACCARTFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title= pConv->GetTitle();

  char buffer[BUFF_SIZE];
  string str;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;

  ifs.getline(buffer,BUFF_SIZE); // keywords
  ifs.getline(buffer,BUFF_SIZE); // filename
  ifs.getline(buffer,BUFF_SIZE); // title (currently ignored)

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      tokenize(vs,buffer);
      if (vs.size() == 0) break;
      else if (vs.size() < 7) return false;
      atom = mol.NewAtom();
      x = atof((char*)vs[1].c_str());
      y = atof((char*)vs[3].c_str());
      z = atof((char*)vs[5].c_str());
      atom->SetVector(x,y,z); //set coordinates

      //set atomic number
      atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
    }

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  mol.SetTitle(title);

 return(true);
}

////////////////////////////////////////////////////////////////

bool MOPACCARTFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  char buffer[BUFF_SIZE];
  
  ofs << "PUT KEYWORDS HERE" << endl;
  ofs << endl;
  ofs << mol.GetTitle() << endl;

  OBAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%-3s%8.5f 1 %8.5f 1 %8.5f 1",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer << endl;
  }
	return(true);
}

} //namespace OpenBabel
