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
#include <openbabel/bond.h>
#include <openbabel/data.h>

#include <cstdlib>

using namespace std;
namespace OpenBabel
{

class XEDFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    XEDFormat()
    {
        OBConversion::RegisterFormat("xed",this);
    }

    virtual const char* Description() //required
    {
        return
          "XED format\n"
          "No comments yet\n";
    };

  virtual const char* SpecificationURL()
  {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return NOTREADABLE | WRITEONEONLY;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

};

//Make an instance of the format class
XEDFormat theXEDFormat;

////////////////////////////////////////////////////////////////

bool XEDFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    int type_name, mass;
    OBAtom *atom;
    OBBond *bond;
    string str,str1;

    ttab.SetFromType("INT");
    ttab.SetToType("XED");
    snprintf(buffer, BUFF_SIZE, "%10.3f%10i%10i",
            mol.GetEnergy(),mol.NumAtoms(),mol.NumBonds());
    ofs << buffer << endl;
    ofs << "File conversion by Open Babel" << endl;

    for (i = 0; i < mol.NumBonds(); i++)
    {
        bond = mol.GetBond(i);
        snprintf(buffer, BUFF_SIZE, "%8i%8i",
                bond->GetBeginAtomIdx(),
                bond->GetEndAtomIdx());
        ofs << buffer;
        if ( !((i+1) % 5) )
            ofs << endl;
    }
    if (mol.NumBonds()%5)
        ofs << endl;

    for(i = 1;i <= mol.NumAtoms(); i++)
    {
        atom = mol.GetAtom(i);
        str = atom->GetType();
        ttab.Translate(str1,str);

        type_name = atoi((char*) str1.c_str());
        switch (type_name)
        {
        case 1:
        case 2:
        case 3:
        case 4:
            mass=6;
            break;
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 23:
        case 25:
            mass=7;
            break;
        case 10:
        case 11:
        case 22:
        case 24:
        case 26:
            mass=8;
            break;
        case 12:
        case 13:
            mass=16;
            break;
        case 14:
            mass=15;
            break;
        case 15:
            mass=1;
            break;
        case 16:
            mass=9;
            break;
        case 17:
            mass=17;
            break;
        case 18:
            mass=35;
            break;
        case 19:
            mass=53;
            break;
        default:
            mass=0;
        }

        snprintf(buffer, BUFF_SIZE, "%6i%15.6f%15.6f%15.6f%6i%12.4f",
                mass, atom->GetX(),atom->GetY(),atom->GetZ(), type_name, 0.0);
        ofs << buffer << endl;
    }
    ofs << "    1         0.0000    0         0.0000" << endl;

    return(true);
}

} //namespace OpenBabel
