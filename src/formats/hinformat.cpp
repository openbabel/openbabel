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

class HINFormat : public OBFormat
{
public:
	//Register this format type ID
	HINFormat() {OBConversion::RegisterFormat("HIN",this);}

	virtual const char* Description() //required
	{ return
"HyperChem INput format\n \
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
HINFormat theHINFormat;

/////////////////////////////////////////////////////////////////
bool HINFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

  // Right now only read in the first molecule
  int i;
  int max, bo;
  char buffer[BUFF_SIZE];
  string str,str1;
  double x,y,z;
  OBAtom *atom;
  vector<string> vs;

  ifs.getline(buffer, BUFF_SIZE);
  while (strstr(buffer,"mol") == NULL)
    ifs.getline(buffer, BUFF_SIZE);
  ifs.getline(buffer, BUFF_SIZE);
  
  mol.BeginModify();
  while (strstr(buffer,"endmol") == NULL)
    {
      tokenize(vs,buffer); // Don't really know how long it'll be
      if (vs.size() < 11) break;
      atom = mol.NewAtom();
      atom->SetAtomicNum(etab.GetAtomicNum(vs[3].c_str()));
      x = atof((char*)vs[7].c_str());
      y = atof((char*)vs[8].c_str());
      z = atof((char*)vs[9].c_str());
      atom->SetVector(x,y,z);
      
      max = 11 + 2 * atoi((char *)vs[10].c_str());
      for (i = 11; i < max; i+=2)
	{
	  switch(((char*)vs[i+1].c_str())[0]) // First char in next token
	    {
	    case 's': bo = 1; break;
	    case 'd': bo = 2; break;
	    case 't': bo = 3; break;
	    case 'a': bo = 5; break;
	    default : bo = 1; break;
	    }
	  mol.AddBond(mol.NumAtoms(), atoi((char *)vs[i].c_str()), bo);
	}
      ifs.getline(buffer, BUFF_SIZE);
    }
  mol.EndModify();

  mol.SetTitle(title);
  return(true);
}

////////////////////////////////////////////////////////////////

bool HINFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i, file_num = 1;
  string str,str1;
  char buffer[BUFF_SIZE];
  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;
  char bond_char;

  ofs << "mol " << file_num << " " << mol.GetTitle() << endl;;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"atom %d - %-3s **  - %8.5f %8.5f  %8.5f  %8.5f %d ",
	    i,
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetPartialCharge(),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ(),
	    atom->GetValence());
    ofs << buffer;
    for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j))
    {
      switch(bond->GetBO())
      {
      case 1 : bond_char = 's'; break;
      case 2 : bond_char = 'd'; break;
      case 3 : bond_char = 't'; break;
      case 5 : bond_char = 'a'; break;
      default: bond_char = 's'; break;
      }
      sprintf(buffer,"%d %c ", (bond->GetNbrAtom(atom))->GetIdx(), bond_char);
      ofs << buffer;
    }
    ofs << endl;
  }
  ofs << "endmol " << file_num << endl;
  return(true);
}

} //namespace OpenBabel
