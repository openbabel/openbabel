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

class CARFormat : public OBFormat
{
public:
	//Register this format type ID
	CARFormat() {OBConversion::RegisterFormat("CAR",this);}

	virtual const char* Description() //required
	{ return
"Biosym chemical modeller input\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://??";}; //optional

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

//Make an instance of the format class
CARFormat theCARFormat;

/////////////////////////////////////////////////////////////////
bool CARFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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
  vector<string> vs;

  mol.BeginModify();

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"PBC") != NULL)
	{
	  if(strstr(buffer,"ON") != NULL)
	    {
	      ifs.getline(buffer,BUFF_SIZE);
	      ifs.getline(buffer,BUFF_SIZE);
	      ifs.getline(buffer,BUFF_SIZE);
	    }
	  else
	    {
	      ifs.getline(buffer,BUFF_SIZE);
	      ifs.getline(buffer,BUFF_SIZE);
	    }
	  break;
	}
    }

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"end") != NULL)
	break;

      atom = mol.NewAtom();

      tokenize(vs,buffer);
      atom->SetAtomicNum(etab.GetAtomicNum(vs[7].c_str()));
      x = atof((char*)vs[1].c_str());
      y = atof((char*)vs[2].c_str());
      z = atof((char*)vs[3].c_str());
      atom->SetVector(x,y,z);
    }
  mol.EndModify();

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  mol.SetTitle(title);
  return(true);
}

} //namespace OpenBabel
