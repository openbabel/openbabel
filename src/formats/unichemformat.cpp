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

#include <cstdlib>

using namespace std;
namespace OpenBabel
{

class UniChemFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    UniChemFormat()
    {
        OBConversion::RegisterFormat("unixyz",this);
    }

  virtual const char* Description() //required
  {
    return
      "UniChem XYZ format\n"
      "Read Options e.g. -as\n"
      " s  Output single bonds only\n"
      " b  Disable bonding entirely\n\n";
  };

  virtual const char* SpecificationURL()
  {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | WRITEONEONLY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
};
//***

//Make an instance of the format class
UniChemFormat theUniChemFormat;

/////////////////////////////////////////////////////////////////
bool UniChemFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int i;
    int natoms;
    char buffer[BUFF_SIZE];

    ifs.getline(buffer,BUFF_SIZE);
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%d", &natoms);
    if (!natoms)
        return(false);

    mol.ReserveAtoms(natoms);
    mol.BeginModify();

    string str;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    for (i = 1; i <= natoms; i ++)
    {
        if (!ifs.getline(buffer,BUFF_SIZE))
            return(false);
        tokenize(vs,buffer);
        if (vs.size() != 4)
            return(false);
        atom = mol.NewAtom();
        x = atof((char*)vs[1].c_str());
        y = atof((char*)vs[2].c_str());
        z = atof((char*)vs[3].c_str());
        atom->SetVector(x,y,z); //set coordinates

        //set atomic number
        atom->SetAtomicNum(atoi((char*)vs[0].c_str()));
    }
    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    // clean out remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    mol.EndModify();
    mol.SetTitle(title);
    return(true);
}

////////////////////////////////////////////////////////////////

bool UniChemFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];

    ofs << mol.GetTitle() << endl;
    ofs << mol.NumAtoms() << endl;

    OBAtom *atom;
    string str,str1;
    for(i = 1;i <= mol.NumAtoms(); i++)
    {
        atom = mol.GetAtom(i);
        snprintf(buffer, BUFF_SIZE,
		"%3d%15.5f%15.5f%15.5f",
                atom->GetAtomicNum(),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs << buffer << endl;
    }

    return(true);
}

} //namespace OpenBabel
