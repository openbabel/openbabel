/**********************************************************************
Copyright (C) 2002-2003 by Geoffrey Hutchison
Based on code Copyright (C) 1999 Joerg-Ruediger Hill
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

class ViewMolFormat : public OBFormat
{
public:
	//Register this format type ID
	ViewMolFormat() {OBConversion::RegisterFormat("VIEWMOL",this);}

	virtual const char* Description() //required
	{ return
"ViewMol format\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://viewmol.sourceforge.net/";}; //optional

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
ViewMolFormat theViewMolFormat;

/////////////////////////////////////////////////////////////////
bool ViewMolFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	istream &ifs = *pConv->GetInStream();
	OBMol &mol = *pmol;
	const char* title = pConv->GetTitle();

  char buffer[BUFF_SIZE];
  OBAtom *atom;
  double x,y,z, border;
  double factor = 1.0;
  int bgn, end, order;
  vector<string> vs;
  bool foundTitle = false;
  bool foundBonds = false;

  mol.BeginModify();

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if (strstr(buffer,"$title") != NULL)
	{
	  if (!ifs.getline(buffer,BUFF_SIZE)) return (false);
	  mol.SetTitle(buffer);
	  foundTitle = true;
	}
      else if (strstr(buffer,"$coord") != NULL)
	{
	  tokenize(vs,buffer);
	  if (vs.size() == 2)
	    factor = atof((char*)vs[1].c_str()); // conversion to angstrom
	  while (ifs.getline(buffer,BUFF_SIZE))
	    {
	      if (buffer[0] == '$') break;
	      tokenize(vs,buffer);
	      if (vs.size() != 4) return(false);
	      atom = mol.NewAtom();
	      x = atof((char*)vs[0].c_str()) * factor;
	      y = atof((char*)vs[1].c_str()) * factor;
	      z = atof((char*)vs[2].c_str()) * factor;
	      atom->SetVector(x,y,z); //set coordinates
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[3].c_str()));
	    }
	}
      else if (strstr(buffer,"$bonds") != NULL) 
	{
	  foundBonds = true;
	  while (ifs.getline(buffer,BUFF_SIZE))
	    {
	      if (buffer[0] == '$') break;
	      sscanf(buffer,"%d %d %lf",&bgn,&end, &border);
	      if (border > 1.0)
		order = int(border);
	      else
		order = 1;
	      mol.AddBond(bgn+1,end+1,order);
	    }
	}
      else if (strstr(buffer,"$end") != NULL)
	break;
    } // while
  
  mol.EndModify();

  if (!foundTitle)
    mol.SetTitle(title);
  if (!foundBonds)
    {
      mol.ConnectTheDots();
      mol.PerceiveBondOrders();
    }
  return(true);
}

////////////////////////////////////////////////////////////////

bool ViewMolFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  char buffer[BUFF_SIZE];
  
  if (strlen(mol.GetTitle()) > 0)
    ofs << "$title" << endl << mol.GetTitle() << endl;

  ofs << "$coord 1.0" << endl;

  OBAtom *atom;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%22.14f%22.14f%22.14f %s",
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ(),
	    etab.GetSymbol(atom->GetAtomicNum()));
    ofs << buffer << endl;
  }

  ofs << "$end" << endl;

  return(true);
}

} //namespace OpenBabel
