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
#include <openbabel/elements.h>
#include <openbabel/obutil.h>

using namespace std;
namespace OpenBabel
{

  class FEATFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    FEATFormat()
    {
      OBConversion::RegisterFormat("feat",this);
    }

    virtual const char* Description() //required
    {
      return
        "Feature format\n"
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
  FEATFormat theFEATFormat;

  /////////////////////////////////////////////////////////////////
  bool FEATFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename

    char buffer[BUFF_SIZE];
    int i,natoms;

    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%d",&natoms);

    mol.ReserveAtoms(natoms);
    mol.BeginModify();

    if (!ifs.getline(buffer,BUFF_SIZE))
      return(false);
    mol.SetTitle(buffer);

    double x,y,z;
    char type[32];
    OBAtom *atom;
    for (i = 0; i < natoms;i++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        sscanf(buffer,"%30s %lf %lf %lf",
               type,
               &x,
               &y,
               &z);
        CleanAtomType(type);
        atom = mol.NewAtom();
        atom->SetVector(x,y,z);
        atom->SetAtomicNum(OBElements::GetAtomicNum(type));
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

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool FEATFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];

    ofs << mol.NumAtoms() << endl;
    ofs << mol.GetTitle() << endl;

    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        snprintf(buffer, BUFF_SIZE, "%-3s %8.5f  %8.5f  %8.5f ",
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 atom->x(),
                 atom->y(),
                 atom->z());
        ofs << buffer << endl;
      }

    return(true);
  }

} //namespace OpenBabel
