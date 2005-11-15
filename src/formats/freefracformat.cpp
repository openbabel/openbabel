/**********************************************************************
Copyright (C) 2005 by Geoffrey R. Hutchison
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
#include "obmolecformat.h"
#include "math/matrix3x3.h"

using namespace std;
namespace OpenBabel
{

class FreeFormFractionalFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    FreeFormFractionalFormat()
    {
        OBConversion::RegisterFormat("fract",this);
    }

  virtual const char* Description() //required
  {
    return
      "Free Form Fractional format\n \
       Read Options e.g. -as\n\
        s  Output single bonds only\n\
        b  Disable bonding entirely\n\n";
  };

  virtual const char* SpecificationURL()
  {return "http://openbabel.sourceforge.net/formats/fract.shtml";}; //optional

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
};
//***

//Make an instance of the format class
FreeFormFractionalFormat theFreeFormFractionalFormat;

const char * TrimErrors(const std::string data)
{
  string temp = data;
  size_t stdErr = temp.rfind("(");
  
  if(stdErr!=string::npos)
    temp.erase(stdErr);
  
  return temp.c_str();
}

/////////////////////////////////////////////////////////////////
bool FreeFormFractionalFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];

    if (!ifs.getline(buffer,BUFF_SIZE))
    {
      obErrorLog.ThrowError(__FUNCTION__,
			    "Problems reading a free form fractional file: Could not read the first line (title/comments).", obWarning);
        return(false);
    }
    if (strlen(buffer) != 0)
        mol.SetTitle(buffer);
    else
        mol.SetTitle(title);

    if (!ifs.getline(buffer,BUFF_SIZE))
    {
      obErrorLog.ThrowError(__FUNCTION__,
			    "Problems reading a free form fractional file: Could not read the second line (unit cell parameters a b c alpha beta gamma).",
			    obWarning);
        return(false);
    }
    vector<string> vs;
    tokenize(vs,buffer," \n\t,");
    if (vs.size() != 6)
        return(false);

    //parse cell values
    double A, B, C, Alpha, Beta, Gamma;
    string temp; // used to trim ending (xx) data from strings
    
    A = atof(TrimErrors(vs[0]));
    B = atof(TrimErrors(vs[1]));
    C = atof(TrimErrors(vs[2]));
    Alpha = atof(TrimErrors(vs[3]));
    Beta  = atof(TrimErrors(vs[4]));
    Gamma = atof(TrimErrors(vs[5]));
    OBUnitCell *uc = new OBUnitCell;
    uc->SetData(A, B, C, Alpha, Beta, Gamma);
    mol.SetData(uc);
    matrix3x3 m = uc->GetOrthoMatrix();

    mol.BeginModify();

    string str;
    double x,y,z;
    vector3 v;
    int atomicNum;
    OBAtom *atom;

    while(ifs.getline(buffer,BUFF_SIZE))
      {
	if (strlen(buffer) == 0 || *buffer == 0x0D) //incl Windows kludge
		// blank line -- consider it the end of this molecule
	  break;

	tokenize(vs,buffer);
        if (vs.size() != 4)
            return(false);

        atom = mol.NewAtom();

	// check to see if first column is number or element symbol
	// (PCModel has files of the form X Y Z symbol)
        atomicNum = etab.GetAtomicNum(vs[0].c_str());
	if (atomicNum == 0 && (isdigit(vs[0][0]) || ispunct(vs[0][0])))
	  {
	    x = atof(TrimErrors(vs[0]));
	    y = atof(TrimErrors(vs[1]));
	    z = atof(TrimErrors(vs[2]));
	    atomicNum = etab.GetAtomicNum(vs[3].c_str());
	  }
	else
	  {
	    x = atof(TrimErrors(vs[1]));
	    y = atof(TrimErrors(vs[2]));
	    z = atof(TrimErrors(vs[3]));
	  }
	v.Set(x, y, z);
	v *= m;	// get cartesian coordinates -- multiply by orthogonalization matrix
        atom->SetVector(v);

        atom->SetAtomicNum(atomicNum);
    }

    // clean out any remaining blank lines
    while(ifs.peek() != EOF && ifs.good() && 
	  (ifs.peek() == '\n' || ifs.peek() == '\r'))
      ifs.getline(buffer,BUFF_SIZE);
    
    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) 
	&& !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    return(true);
}

////////////////////////////////////////////////////////////////

bool FreeFormFractionalFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    OBUnitCell *uc = NULL;

    ofs << mol.GetTitle() << endl;

    if (!mol.HasData(OBGenericDataType::UnitCell))
      {
	sprintf(buffer,
		"%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f",
		1.0f, 1.0f, 1.0f, 90.0f, 90.0f, 90.0f);
	ofs << buffer << endl;
      }
    else
      {
	uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
	sprintf(buffer,
		"%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f",
		uc->GetA(), uc->GetB(), uc->GetC(),
		uc->GetAlpha() , uc->GetBeta(), uc->GetGamma());
	ofs << buffer << endl;
      }

    vector3 v;
    FOR_ATOMS_OF_MOL(atom, mol)
    {
      v = atom->GetVector();
      if (uc != NULL)
	v *= uc->GetFractionalMatrix();

      sprintf(buffer,"%s %10.5f%10.5f%10.5f",
	      etab.GetSymbol(atom->GetAtomicNum()),
	      v.x(),
	      v.y(),
	      v.z());
        ofs << buffer << endl;
    }
    ofs << endl; // add a blank line between molecules
    return(true);
}

} //namespace OpenBabel
