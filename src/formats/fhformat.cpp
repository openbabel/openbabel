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

class FenskeZmatFormat : public OBFormat
{
public:
	//Register this format type ID
	FenskeZmatFormat() {OBConversion::RegisterFormat("FH",this);}

	virtual const char* Description() //required
	{ return
"FenskeZmat format\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://??";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return NOTREADABLE;};

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

//Make an instance of the format class
FenskeZmatFormat theFenskeZmatFormat;

/////////////////////////////////////////////////////////////////

bool FenskeZmatFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  OBAtom *atom,*a,*b,*c;
  char type[10],buffer[BUFF_SIZE];
  vector<OBNodeBase*>::iterator i;
  
  vector<OBInternalCoord*> vic;
  vic.push_back((OBInternalCoord*)NULL);
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    vic.push_back(new OBInternalCoord);

  CartesianToInternal(vic,mol);

  ofs << endl << mol.NumAtoms() << endl;
    
  double r,w,t;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
  { 
    a = vic[atom->GetIdx()]->_a;
    b = vic[atom->GetIdx()]->_b;
    c = vic[atom->GetIdx()]->_c;
    r = vic[atom->GetIdx()]->_dst;
    w = vic[atom->GetIdx()]->_ang;
    t = vic[atom->GetIdx()]->_tor;
    strcpy(type,etab.GetSymbol(atom->GetAtomicNum()));

    if (atom->GetIdx() == 1)
    {
      sprintf(buffer,"%-2s  1",type); 
      ofs << buffer << endl;
      continue;
    }

    if (atom->GetIdx() == 2)
    {
      sprintf(buffer,"%-2s%3d%6.3f",
	      type,
	      a->GetIdx(),r);
      ofs << buffer << endl;
      continue;
    }

    if (atom->GetIdx() == 3)
    {
      sprintf(buffer,"%-2s%3d%6.3f%3d%8.3f",
	      type,
	      a->GetIdx(),r, b->GetIdx(),w);
      ofs << buffer << endl;
      continue;
    }

    if (t < 0) t += 360;
    sprintf(buffer,"%-2s%3d%6.3f%3d%8.3f%3d%6.1f",
	    type,
	    a->GetIdx(),r,b->GetIdx(),w,c->GetIdx(),t);
    ofs << buffer << endl;
  }

  return(true);
}

} //namespace OpenBabel
