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


using namespace std;
namespace OpenBabel
{

  class CacheFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CacheFormat()
    {
      OBConversion::RegisterFormat("cac",this);
      OBConversion::RegisterFormat("cache",this);
    }

    virtual const char* Description() //required
    {
      return
        "CAChe MolStruct format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-cache"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ////////////////////////////////////////////////////
  };

  //Make an instance of the format class
  CacheFormat theCacheFormat;

  ////////////////////////////////////////////////////////////////

  bool CacheFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char type_name[16];
    char buffer[BUFF_SIZE];

    ofs << "molstruct88_Apr_30_1993_11:02:29 <molecule> 0x1d00\n";
    ofs << "Written by Molecular Editor on <date>\n";
    ofs << "Using data dictionary         9/9/93  4:47 AM\n";
    ofs << "Version 6\n";
    ofs << "local_transform\n";
    ofs << "0.100000 0.000000 0.000000 0.000000\n";
    ofs << "0.000000 0.100000 0.000000 0.000000\n";
    ofs << "0.000000 0.000000 0.100000 0.000000\n";
    ofs << "0.000000 0.000000 0.000000 1.000000\n";
    ofs << "object_class atom\n";
    ofs << "property xyz_coordinates MoleculeEditor angstrom 6 3 FLOAT\n";
    ofs << "property anum MoleculeEditor unit 0 1 INTEGER\n";
    ofs << "property sym MoleculeEditor noUnit 0 2 STRING\n";
    ofs << "property chrg MoleculeEditor charge_au 0 1 INTEGER\n";
    ofs << "property rflag MoleculeEditor noUnit 0 1 HEX\n";
    ofs << "ID xyz_coordinates             anum sym	chrg rflag\n";

    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        // 16 = sizeof(type_name)
        strncpy(type_name,OBElements::GetSymbol(atom->GetAtomicNum()), 16);
        // sizeof(type_name) - 1)
        type_name[15] = '\0';

        snprintf(buffer, BUFF_SIZE, "%3d %10.6f %10.6f %10.6f %2d %2s %2d 0x7052",
                atom->GetIdx(),
                atom->x(),
                atom->y(),
                atom->z(),
                atom->GetAtomicNum(),
                type_name,
                atom->GetFormalCharge());
        ofs << buffer << endl;
      }

    ofs << "property_flags:\n";
    ofs << "object_class bond\n";
    ofs << "property rflag MoleculeEditor noUnit 0 1 HEX\n";
    ofs << "property type MoleculeEditor noUnit 0 1 NAME\n";
    ofs << "property bond_order MoleculeEditor noUnit 4 1 FLOAT\n";
    ofs << "ID rflag type bond_order\n";

    char bstr[16];
    OBBond *bond;
    vector<OBBond*>::iterator j;
    for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
      {
        switch (bond->GetBondOrder())
          {
          case 1:
            strcpy(bstr,"single");
            break;
          case 2:
            strcpy(bstr,"double");
            break;
          case 3:
            strcpy(bstr,"triple");
            break;
          default:
            strcpy(bstr,"weak");
          }

        snprintf(buffer, BUFF_SIZE, "%3d 0x7005 %s\n", bond->GetIdx()+1,bstr);
        ofs << buffer;
      }

    ofs << "property_flags:\n";
    ofs << "object_class connector\n";
    ofs << "property dflag MoleculeEditor noUnit 0 1 HEX\n";
    ofs << "property objCls1 MoleculeEditor noUnit 0 1 NAME\n";
    ofs << "property objCls2 MoleculeEditor noUnit 0 1 NAME\n";
    ofs << "property objID1 MoleculeEditor noUnit 0 1 INTEGER\n";
    ofs << "property objID2 MoleculeEditor noUnit 0 1 INTEGER\n";
    ofs << "ID dflag objCls1 objCls2 objID1 objID2\n";


    int k;
    for (bond = mol.BeginBond(j),k=1;bond;bond = mol.NextBond(j))
      {
        snprintf(buffer, BUFF_SIZE, "%3d 0xa1 atom bond %d %d\n",
                k++,bond->GetBeginAtomIdx(),bond->GetIdx()+1);
        ofs << buffer;
        snprintf(buffer, BUFF_SIZE, "%3d 0xa1 atom bond %d %d\n",
                k++,bond->GetEndAtomIdx(),bond->GetIdx()+1);
        ofs << buffer;
      }

    ofs << "property_flags:\n";
    return(true);
  }

} //namespace OpenBabel
