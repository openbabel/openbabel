/**********************************************************************
Copyright (C) 2000-2004 by Geoffrey Hutchison
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

class QChemFormat : public OBFormat
{
public:
	//Register this format type ID
	QChemFormat() {OBConversion::RegisterFormat("QCHEM",this);}

	virtual const char* Description() //required
	{ return
"QChem\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://www.q-chem.com/";}; //optional

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
QChemFormat theQChemFormat;

/////////////////////////////////////////////////////////////////
bool QChemFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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
//  OBInternalCoord *coord; CM
  vector<string> vs;
  vector<OBInternalCoord *> internals; // If we get a z-matrix
  int charge = 0;
  unsigned int spin = 1;
  bool hasPartialCharges = false;

  mol.BeginModify();

  while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"Standard Nuclear Orientation") != NULL)
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
	      atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));
	      x = atof((char*)vs[2].c_str());
	      y = atof((char*)vs[3].c_str());
	      z = atof((char*)vs[4].c_str());
	      atom->SetVector(x,y,z);

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      else if(strstr(buffer,"Mulliken Net Atomic Charges") != NULL)
	{
	  hasPartialCharges = true;
	  ifs.getline(buffer,BUFF_SIZE);	// (blank)
	  ifs.getline(buffer,BUFF_SIZE);	// column headings
	  ifs.getline(buffer,BUFF_SIZE);	// -----------------
	  ifs.getline(buffer,BUFF_SIZE);
	  tokenize(vs,buffer);
	  while (vs.size() >= 3)
	    {
	      atom = mol.GetAtom(atoi(vs[0].c_str()));
	      atom->SetPartialCharge(atof(vs[2].c_str()));

	      if (!ifs.getline(buffer,BUFF_SIZE)) break;
	      tokenize(vs,buffer);
	    }
	}
      // In principle, the final geometry in cartesians is exactly the same
      // as this geometry, so we shouldn't lose any information
      // but we grab the charge and spin from this $molecule block
      else if(strstr(buffer,"OPTIMIZATION CONVERGED") != NULL)
 	{
  	  // Unfortunately this will come in a Z-matrix, which would
	  // change our frame of reference -- we'll ignore this geometry
  	  ifs.getline(buffer,BUFF_SIZE); // *****
  	  ifs.getline(buffer,BUFF_SIZE); // blank
 	  ifs.getline(buffer,BUFF_SIZE); // Z-matrix Print:
 	  ifs.getline(buffer,BUFF_SIZE); // $molecule
 	  ifs.getline(buffer,BUFF_SIZE); // Charge,Spin
	  tokenize(vs, buffer, ", \t\n");
	  if (vs.size() == 2)
	    {
	      charge = atoi(vs[0].c_str());
	      spin = atoi(vs[1].c_str());
	    }

 	  ifs.getline(buffer,BUFF_SIZE);
	  
	  // Uncomment the following to use the z-matrix coordinates
	  // 	  index = 1;
	  // 	  while (strstr(buffer, "$end") == NULL)
	  // 	    {
	  // 	      tokenize(vs,buffer);
	  // 	      atom = mol.NewAtom();
	  // 	      atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));
	  //
	  // 	      tokenize(vs,buffer);
	  // 	      coord = new OBInternalCoord();
	  // 	      if (index > 1)
	  // 		{
	  // 		  coord->_a = mol.GetAtom(atoi(vs[2].c_str()));
	  // 		  coord->_dst = atof(vs[3].c_str());
	  // 		}
	  // 	      if (index > 2)
	  // 		{
	  // 		  coord->_b = mol.GetAtom(atoi(vs[4].c_str()));
	  // 		  coord->_ang = atof(vs[5].c_str()) * DEG_TO_RAD;
	  // 		}
	  // 	      if (index > 3)
	  // 		{
	  // 		  coord->_c = mol.GetAtom(atoi(vs[6].c_str()));
	  // 		  coord->_tor = atof(vs[7].c_str()) * DEG_TO_RAD;
	  // 		}
	  // 	      index++;
	  // 	      ifs.getline(buffer,BUFF_SIZE);
	  // 	      internals.push_back(coord);
	  //   }  // end while
	  //	  InternalToCartesian(internals, mol);
	} // end (OPTIMIZATION CONVERGED)
    } // end while
  mol.EndModify();
  mol.ConnectTheDots();
  mol.PerceiveBondOrders();
  if (hasPartialCharges)
    mol.SetPartialChargesPerceived();
  mol.SetTotalCharge(charge);
  mol.SetTotalSpinMultiplicity(spin);

  mol.SetTitle(title);
  return(true);
}

////////////////////////////////////////////////////////////////

bool QChemFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  OBAtom *atom;
  
  ofs << "$comment" << endl;
  ofs << mol.GetTitle() << endl;
  ofs << "$end" << endl;
  ofs << endl << "$molecule" << endl;
  ofs << mol.GetTotalCharge() << " " << mol.GetTotalSpinMultiplicity() << endl;

  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    ofs << atom->GetAtomicNum() << " "
	<< atom->GetX() << " " << atom->GetY() << " " << atom->GetZ() << endl;
  }
  ofs << "$end" << endl;
  ofs << endl << "$rem" << endl << "$end" << endl;

  return(true);
}

} //namespace OpenBabel
