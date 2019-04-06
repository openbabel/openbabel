/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/internalcoord.h>

#include <cstdlib>

using namespace std;
namespace OpenBabel
{

class AmberPrepFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    AmberPrepFormat()
    {
        OBConversion::RegisterFormat("prep",this);
    }

  virtual const char* Description() //required
  {
    return
      "Amber Prep format\n"
      " Read Options e.g. -as\n"
      " s  Output single bonds only\n"
      " b  Disable bonding entirely\n\n";
  };

  virtual const char* SpecificationURL()
  {return "http://amber.scripps.edu/doc/prep.html";};

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return NOTWRITABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
};

//Make an instance of the format class
AmberPrepFormat theAmberPrepFormat;

/////////////////////////////////////////////////////////////////
bool AmberPrepFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    OBAtom *atom;
    OBInternalCoord *coord;
    vector<string> vs;
    vector<OBInternalCoord*> internals;

    mol.BeginModify();

    while	(ifs.getline(buffer,BUFF_SIZE))
    {
        tokenize(vs,buffer);
        if (vs.size() == 10)
        {
            atom = mol.NewAtom();
            coord = new OBInternalCoord();
            if (mol.NumAtoms() > 1)
                coord->_a = mol.GetAtom(atoi(vs[4].c_str()));
            if (mol.NumAtoms() > 2)
                coord->_b = mol.GetAtom(atoi(vs[5].c_str()));
            if (mol.NumAtoms() > 3)
                coord->_c = mol.GetAtom(atoi(vs[6].c_str()));
            coord->_dst = atof(vs[7].c_str());
            coord->_ang = atof(vs[8].c_str());
            coord->_tor = atof(vs[9].c_str());
            internals.push_back(coord);

            atom->SetAtomicNum(OBElements::GetAtomicNum(vs[1].c_str()));

            if (!ifs.getline(buffer,BUFF_SIZE))
                break;
            tokenize(vs,buffer);
        }
    }

    if (internals.size() > 0)
      InternalToCartesian(internals,mol);

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    mol.SetTitle(title);
    return(true);
}

} //namespace OpenBabel
