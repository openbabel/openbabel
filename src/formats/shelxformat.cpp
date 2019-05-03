/**********************************************************************
Copyright (C) 1998-2003 by OpenEye Scientific Software, Inc.
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

#include <openbabel/math/matrix3x3.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class ShelXFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ShelXFormat()
    {
      OBConversion::RegisterFormat("res",this, "chemical/x-shelx");
      OBConversion::RegisterFormat("ins",this, "chemical/x-shelx");
    }

    virtual const char* Description() //required
    {
      return
        "ShelX format\n"
        "Read Options e.g. -as\n"
        " s  Output single bonds only\n"
        " b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://shelx.uni-ac.gwdg.de/SHELX/";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-shelx"; };

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
  ShelXFormat theShelXFormat;

  /////////////////////////////////////////////////////////////////
  bool ShelXFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename

    char buffer[BUFF_SIZE];
    //  int natoms; CM
    double A,B,C,Alpha,Beta,Gamma;
    matrix3x3 m;

    ifs.getline(buffer,BUFF_SIZE);
    mol.SetTitle(buffer);
    while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"CELL",4))
      ;

    if (!EQn(buffer,"CELL",4))
      return(false);
    vector<string> vs;
    tokenize(vs,buffer," \n\t,");
    if (vs.size() != 8)
      return(false);

    //parse cell values
    A = atof((char*)vs[2].c_str());
    B = atof((char*)vs[3].c_str());
    C = atof((char*)vs[4].c_str());
    Alpha = atof((char*)vs[5].c_str());
    Beta  = atof((char*)vs[6].c_str());
    Gamma = atof((char*)vs[7].c_str());
    OBUnitCell *uc = new OBUnitCell;
    uc->SetOrigin(fileformatInput);
    uc->SetData(A, B, C, Alpha, Beta, Gamma);
    mol.SetData(uc);

    //  int i; CM
    double x,y,z;
    char type[16], *j;
    OBAtom *atom;
    vector3 v;

    //skip until FVAR found
    while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"FVAR",4))
      ;

    mol.BeginModify();

    //read atom records
    while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"HKLF",4))
      {
        tokenize(vs,buffer," \n\t,");

        //skip AFIX and PART instructions
        if (vs.size() < 7)
          continue;
        atom = mol.NewAtom();

        x = atof((char*)vs[2].c_str());
        y = atof((char*)vs[3].c_str());
        z = atof((char*)vs[4].c_str());
        v.Set(x,y,z);
        v = uc->FractionalToCartesian(v);

        strncpy(type,vs[0].c_str(), sizeof(type));
        type[sizeof(type) - 1] = '\0';
        j = strpbrk(type, "0123456789");
        j[0] = '\0';
        atom->SetAtomicNum(OBElements::GetAtomicNum(type));
        atom->SetVector(v);

        //skip next line if anisotropic atoms.
        if (vs.size() == 9)
          ifs.getline(buffer,BUFF_SIZE);
      } //while

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    return(true);
  }

} //namespace OpenBabel
