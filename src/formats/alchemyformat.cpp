/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
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
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include <openbabel/data.h>
#include <cstdlib>


using namespace std;
namespace OpenBabel
{

  class AlchemyFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    AlchemyFormat()
    {
      OBConversion::RegisterFormat("alc",this, "chemical/x-alchemy");
    }

    virtual const char* Description() //required
    {
      return
        "Alchemy format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    { return "";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-alchemy"; };

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
  AlchemyFormat theAlchemyFormat;

  /////////////////////////////////////////////////////////////////
  bool AlchemyFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int i;
    int natoms = 0, nbonds = 0;
    char buffer[BUFF_SIZE];

    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer," %d %*s %d", &natoms, &nbonds);
    if (!natoms)
      {
        ifs.getline(buffer,BUFF_SIZE);
        sscanf(buffer," %d %*s %d", &natoms, &nbonds);
        if (!natoms)
          return(false);
      }

    mol.ReserveAtoms(natoms);
    mol.BeginModify();
    ttab.SetFromType("ALC");

    string str;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    for (i = 1; i <= natoms; i ++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        tokenize(vs,buffer);
        if (vs.size() < 5)
          return(false);
        atom = mol.NewAtom();
        x = atof((char*)vs[2].c_str());
        y = atof((char*)vs[3].c_str());
        z = atof((char*)vs[4].c_str());
        atom->SetVector(x,y,z); //set coordinates

        //set atomic number
        ttab.SetToType("ATN");
        ttab.Translate(str,vs[1]);
        atom->SetAtomicNum(atoi(str.c_str()));
        //set type
        ttab.SetToType("INT");
        ttab.Translate(str,vs[1]);
        atom->SetType(str);
      }

    char bobuf[100];
    string bostr;
    int bgn,end,order;

    for (i = 0; i < nbonds; i++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        sscanf(buffer," %*d%d%d%99s",&bgn,&end,bobuf);
        bostr = bobuf;
        order = 1;
        if      (bostr == "DOUBLE")
          order = 2;
        else if (bostr == "TRIPLE")
          order = 3;
        else if (bostr == "AROMATIC")
          order = 5;
        mol.AddBond(bgn,end,order);
      }

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

  bool AlchemyFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    char bond_string[10];

    snprintf(buffer, BUFF_SIZE, "%5d ATOMS, %5d BONDS,     0 CHARGES",
             mol.NumAtoms(),
             mol.NumBonds());
    ofs << buffer << endl;

    OBAtom *atom;
    string str,str1;
    for(i = 1;i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        str = atom->GetType();
        ttab.SetFromType("INT");
        ttab.SetToType("ALC");
        ttab.Translate(str1,str);
        snprintf(buffer, BUFF_SIZE, "%5d %-6s%8.4f %8.4f %8.4f     0.0000",
                 i,
                 (char*)str1.c_str(),
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ());
        ofs << buffer << endl;
      }

    OBBond *bond;
    vector<OBBond*>::iterator j;

    for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
      {
        switch(bond->GetBondOrder())
          {
          case 1 :
            strcpy(bond_string,"SINGLE");
            break;
          case 2 :
            strcpy(bond_string,"DOUBLE");
            break;
          case 3 :
            strcpy(bond_string,"TRIPLE");
            break;
          case 5 :
            strcpy(bond_string,"AROMATIC");
            break;
          default :
            strcpy(bond_string,"SINGLE");
          }
        snprintf(buffer, BUFF_SIZE, "%5d  %4d  %4d  %s",
                 bond->GetIdx()+1,
                 bond->GetBeginAtomIdx(),
                 bond->GetEndAtomIdx(),
                 bond_string);
        ofs << buffer << endl;
      }
    return(true);
  }

} //namespace OpenBabel
