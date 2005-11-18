/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
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

using namespace std;
namespace OpenBabel
{

class MOPACFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    MOPACFormat()
    {
        OBConversion::RegisterFormat("mopout",this);
    }

  virtual const char* Description() //required
  {
    return
      "MOPAC Output format\n \
       Read Options e.g. -as\n\
        s  Output single bonds only\n\
        b  Disable bonding entirely\n\n";
  };

    virtual unsigned int Flags()
    {
        return NOTWRITABLE;
    };

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    //	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv); Is Read Only
};

//Make an instance of the format class
MOPACFormat theMOPACFormat;

/////////////////////////////////////////////////////////////////
bool MOPACFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title= pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    vector<double> charges;
    bool hasPartialCharges = false;
    double energy;

    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
    {
        if(strstr(buffer,"CARTESIAN COORDINATES") != NULL)
        {
            // mol.EndModify();
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);	// blank
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

                if (!ifs.getline(buffer,BUFF_SIZE))
                    break;
                tokenize(vs,buffer);
            }
        }
        else if(strstr(buffer,"FINAL HEAT") != NULL)
        {
            sscanf(buffer,"%*s%*s%*s%*s%*s%lf",&energy);
            mol.SetEnergy(energy);
        }
        else if(strstr(buffer,"NET ATOMIC CHARGES") != NULL)
        {
            hasPartialCharges = true;
            charges.clear();
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 4)
            {
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                atom->SetPartialCharge(atof(vs[2].c_str()));
                charges.push_back(atof(vs[2].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE))
                    break;
                tokenize(vs,buffer);
            }
        }
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();

    if (hasPartialCharges)
    {
        mol.SetPartialChargesPerceived();
        for (unsigned int i = 1; i <= mol.NumAtoms(); i++)
        {
            atom = mol.GetAtom(i);
            atom->SetPartialCharge(charges[i-1]);
        }
    }
    mol.SetTitle(title);

    return(true);
}

//************************************************************
class MOPACCARTFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    MOPACCARTFormat()
    {
        OBConversion::RegisterFormat("mopcrt",this, "chemical/x-mopac-input");
    }

  virtual const char* Description() //required
  {
    return
      "MOPAC Cartesian format\n \
       Options e.g. -xs\n\
        s  Output single bonds only\n\
        b  Disable bonding entirely\n";
  };

     /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ////////////////////////////////////////////////////
};

//Make an instance of the format class
MOPACCARTFormat theMOPACCARTFormat;

/////////////////////////////////////////////////////////////////
bool MOPACCARTFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title= pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    ifs.getline(buffer,BUFF_SIZE); // keywords
    ifs.getline(buffer,BUFF_SIZE); // filename
    ifs.getline(buffer,BUFF_SIZE); // title (currently ignored)

    while	(ifs.getline(buffer,BUFF_SIZE))
    {
        tokenize(vs,buffer);
        if (vs.size() == 0)
            break;
        else if (vs.size() < 7)
            return false;
        atom = mol.NewAtom();
        x = atof((char*)vs[1].c_str());
        y = atof((char*)vs[3].c_str());
        z = atof((char*)vs[5].c_str());
        atom->SetVector(x,y,z); //set coordinates

        //set atomic number
        atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();
    mol.SetTitle(title);

    return(true);
}

////////////////////////////////////////////////////////////////

bool MOPACCARTFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];

    ofs << "PUT KEYWORDS HERE" << endl;
    ofs << endl;
    ofs << mol.GetTitle() << endl;

    OBAtom *atom;
    string str,str1;
    for(i = 1;i <= mol.NumAtoms(); i++)
    {
        atom = mol.GetAtom(i);
        sprintf(buffer,"%-3s%8.5f 1 %8.5f 1 %8.5f 1",
                etab.GetSymbol(atom->GetAtomicNum()),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs << buffer << endl;
    }
    return(true);
}

} //namespace OpenBabel
