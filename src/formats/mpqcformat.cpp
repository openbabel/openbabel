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

class MPQCFormat : public OBFormat
{
public:
	//Register this format type ID
	MPQCFormat() {OBConversion::RegisterFormat("MPQC",this);}

	virtual const char* Description() //required
	{ return
"MPQC format\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://??";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return NOTWRITABLE;};

	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

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
MPQCFormat theMPQCFormat;

/////////////////////////////////////////////////////////////////
bool MPQCFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

  char buffer[BUFF_SIZE];
  string str,str1;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;
  bool bohr = true;

  mol.BeginModify();
  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"<Molecule>:") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  while	(strstr(buffer,"geometry") == NULL)
	    {
	      if (strstr(buffer,"angstrom") != NULL)
		bohr = false;
	      if (!ifs.getline(buffer,BUFF_SIZE))
		return(false);
	    }
	  ifs.getline(buffer,BUFF_SIZE); // Now we're on the atoms
	  tokenize(vs,buffer);
	  while (vs.size() == 6)
	    {
	      if (bohr)
		{
		  x = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
		  y = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
		  z = atof((char*)vs[5].c_str()) * BOHR_TO_ANGSTROM;
		}
	      else {
		  x = atof((char*)vs[3].c_str());
		  y = atof((char*)vs[4].c_str());
		  z = atof((char*)vs[5].c_str());
		}
	      atom = mol.NewAtom();
	      atom->SetVector(x,y,z);
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
    }

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  mol.SetTitle(title);
  
  return(true);
}


} //namespace OpenBabel
