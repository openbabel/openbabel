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

using namespace std;
namespace OpenBabel
{

class TinkerFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    TinkerFormat()
    {
        OBConversion::RegisterFormat("txyz",this);
    }

    virtual const char* Description() //required
    {
        return
            "Tinker MM2 format\n \
            No comments yet\n";
    };

  virtual const char* SpecificationURL()
  {return "http://dasher.wustl.edu/tinker/";}; //optional

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
TinkerFormat theTinkerFormat;

/////////////////////////////////////////////////////////////////

bool TinkerFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    OBBond *bond;
    vector<OBBond*>::iterator j;

    snprintf(buffer, BUFF_SIZE, "%6d %-20s   MM2 parameters\n",mol.NumAtoms(),mol.GetTitle());
    ofs << buffer;

    ttab.SetFromType("INT");

    OBAtom *atom;
    string str,str1;
    for(i = 1;i <= mol.NumAtoms(); i++)
    {
        atom = mol.GetAtom(i);
        str = atom->GetType();
        ttab.SetToType("MM2");
        ttab.Translate(str1,str);
        snprintf(buffer, BUFF_SIZE, "%6d %2s  %12.6f%12.6f%12.6f %5d",
                i,
                etab.GetSymbol(atom->GetAtomicNum()),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ(),
                atoi((char*)str1.c_str()));
        ofs << buffer;

        for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j))
        {
            snprintf(buffer, BUFF_SIZE, "%6d", (bond->GetNbrAtom(atom))->GetIdx());
            ofs << buffer;
        }

        ofs << endl;
    }

    return(true);
}

} //namespace OpenBabel
