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
#include <ctype.h>
#include "obconversion.h"

using namespace std;
namespace OpenBabel {

class JaguarFormat : public OBFormat
{
public:
	//Register this format type ID
	JaguarFormat() {OBConversion::RegisterFormat("JAGUAR",this);}

	virtual const char* Description() //required
	{ return
"Jaguar format\n \
In and Out formats?\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http:/??";}; //optional

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
JaguarFormat theJaguarFormat;

/////////////////////////////////////////////////////////////////
bool JaguarFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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
  unsigned int i;
  OBAtom *atom;
  vector<string> vs;

  mol.BeginModify();
  while (ifs.getline(buffer,BUFF_SIZE))
  {
    if (strstr(buffer,"Input geometry:") != NULL
        || strstr(buffer,"symmetrized geometry:") != NULL
        || strstr(buffer,"new geometry:") != NULL
        || strstr(buffer,"final geometry:") != NULL)
    {
      // mol.EndModify();
      mol.Clear();
      mol.BeginModify();
      ifs.getline(buffer,BUFF_SIZE);  //     angstroms
      ifs.getline(buffer,BUFF_SIZE);  // column headings
      ifs.getline(buffer,BUFF_SIZE);
      tokenize(vs,buffer);
      while (vs.size() == 4)
      {
        atom = mol.NewAtom();
        str = vs[0]; // Separate out the Symbol# into just Symbol ...
        for (i = 0;i < str.size();i++)
          if (isdigit(str[i]))
            str[i] = '\0';

        atom->SetAtomicNum(etab.GetAtomicNum(str.c_str()));
        x = atof((char*)vs[1].c_str());
        y = atof((char*)vs[2].c_str());
        z = atof((char*)vs[3].c_str());
        atom->SetVector(x,y,z);

        if (!ifs.getline(buffer,BUFF_SIZE)) break;
        tokenize(vs,buffer);
      }
    }
    if (strstr(buffer, "Atomic charges from electrostatic potential") != NULL)
    { 
      mol.SetAutomaticPartialCharge(false);
      unsigned int chgcount=0;
      while (chgcount<mol.NumAtoms())
      {
        ifs.getline(buffer,BUFF_SIZE);  // blank line
        ifs.getline(buffer,BUFF_SIZE);  // header line
        ifs.getline(buffer,BUFF_SIZE);  // data line
        tokenize(vs,buffer);
        for (vector<string>::size_type icount=1;icount<vs.size();++icount)
        {
          chgcount=chgcount+1;
          mol.GetAtom(chgcount)->SetPartialCharge(atof((char*)vs[icount].c_str()));
        }
      }
    }
  }
  mol.EndModify();

  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  mol.SetTitle(title);
  return(true);
}

////////////////////////////////////////////////////////////////

bool JaguarFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  char buffer[BUFF_SIZE];
  OBAtom *atom;

  ofs << mol.GetTitle() << endl << endl;
  ofs << "&gen" << endl;
  ofs << "&" << endl;
  ofs << "&zmat" << endl;

  for (i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"  %s%d   %12.7f  %12.7f  %12.7f",
            etab.GetSymbol(atom->GetAtomicNum()), i,
            atom->GetX(),
            atom->GetY(),
            atom->GetZ());
    ofs << buffer << endl;
  }

  ofs << "&" << endl;
  return(true);
}

} //namespace OpenBabel
