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

#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989

class GAMESSFormat : public OBFormat
{
public:
	//Register this format type ID
	GAMESSFormat() {OBConversion::RegisterFormat("GAMESS",this);}

	virtual const char* Description() //required
	{ return
"GAMESS\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://www.msg.ameslab.gov/GAMESS/doc.menu.html";}; //optional

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
GAMESSFormat theGAMESSFormat;

/////////////////////////////////////////////////////////////////
bool GAMESSFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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
  bool hasPartialCharges = false;

  mol.BeginModify();
  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"ATOMIC                      COORDINATES (BOHR)") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  mol.BeginModify();
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 5)
	    {
	      atom = mol.NewAtom();
	      atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
	      x = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
	      y = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
	      z = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
	      atom->SetVector(x,y,z);
	      vs[1].erase(vs[1].size() - 2, 2);

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      else if(strstr(buffer,"COORDINATES OF ALL ATOMS ARE (ANGS)") != NULL)
	{
	  // mol.EndModify();
	  mol.Clear();
	  mol.BeginModify();
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);	// ---------------
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 5)
	    {
	      atom = mol.NewAtom();
	      atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
	      x = atof((char*)vs[2].c_str());
	      y = atof((char*)vs[3].c_str());
	      z = atof((char*)vs[4].c_str());
	      atom->SetVector(x,y,z);
	      vs[1].erase(vs[1].size() - 2, 2);

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      else if(strstr(buffer,"MOPAC CHARGES") != NULL)
	{
	  hasPartialCharges = true;
	  ifs.getline(buffer,BUFF_SIZE);	// ---------------
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() == 4)
	    {
	      atom = mol.GetAtom(atoi(vs[0].c_str()));
	      atom->SetPartialCharge(atof(vs[2].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
    }
  mol.EndModify();
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();

  if (hasPartialCharges)
    mol.SetPartialChargesPerceived();
  mol.SetTitle(title);
  return(true);
}

////////////////////////////////////////////////////////////////

bool GAMESSFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  char buffer[BUFF_SIZE];

  ofs << " $CONTRL COORD=CART UNITS=ANGS $END" << endl;
  ofs << " $DATA" << endl;
  ofs << mol.GetTitle() << endl;
  ofs << "Put symmetry info here" << endl << endl;
  
  OBAtom *atom;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    sprintf(buffer,"%-3s %4d.0    %8.5f  %8.5f  %8.5f ",
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetAtomicNum(),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ());
    ofs << buffer << endl;
  }

  ofs << " $END" << endl << endl << endl;
  return(true);
}

} //namespace OpenBabel
