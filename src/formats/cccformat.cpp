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

class CCCFormat : public OBFormat
{
public:
	//Register this format type ID
	CCCFormat() {OBConversion::RegisterFormat("CCC",this);}

	virtual const char* Description() //required
	{ return
"CCC format\n \
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

//Make an instance of the format class
CCCFormat theCCCFormat;

/////////////////////////////////////////////////////////////////
bool CCCFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

  char buffer[BUFF_SIZE];
  ifs.getline(buffer,BUFF_SIZE);

  if (strlen(buffer) > 5) mol.SetTitle(&buffer[5]);
  mol.SetEnergy(0.0);

  int natoms;
  ifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%*s%d",&natoms);
  mol.ReserveAtoms(natoms);

  int end,order;
  double x,y,z;
  OBAtom atom;
  vector3 v;
  vector<string> vs;
  char element[3];
  element[2] = '\0';
  
  for (int i = 1;i <= natoms;i++)
    {
      if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
      atom.Clear();
      element[0] = buffer[0];
      element[1] = (buffer[1] != ' ') ? buffer[1]:'\0';
      atom.SetAtomicNum(etab.GetAtomicNum(element));
      sscanf(&buffer[15],"%lf%lf%lf",&x,&y,&z);
      v.Set(x,y,z);
      atom.SetVector(v);

      if (!mol.AddAtom(atom)) return(false);
      tokenize(vs,&buffer[60]);
      vector<string>::iterator j;

      for (j = vs.begin();j != vs.end();j++)
	if (!j->empty())
	{
	  //get the bond order
	  switch((char)(*j)[j->size()-1])
	    {
	    case 'S':  order = 1; break;
	    case 'D':  order = 2; break;
	    case 'T':  order = 3; break;
	    default: order = 1;
	    }
	  (*j)[j->size()-1] = ' ';
	  end = atoi(j->c_str());
	  if (i>end)
	    mol.AddBond(i,end,order);
	}
    }

  return(true);
}

} //namespace OpenBabel
