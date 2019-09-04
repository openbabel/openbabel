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

#include <openbabel/math/matrix3x3.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/data.h>
#include <cstdlib>


using namespace std;
namespace OpenBabel
{

  class CHEM3D1Format : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CHEM3D1Format()
    {
      OBConversion::RegisterFormat("c3d1",this);
    }

    virtual const char* Description() //required
    {
      return
        "Chem3D Cartesian 1 format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    static bool ReadChem3d(istream &ifs,OBMol &mol,bool mmads,const char *type_key);
    static bool WriteChem3d(ostream &ofs,OBMol &mol, const char *mol_typ);

  };
  //***

  //Make an instance of the format class
  CHEM3D1Format theCHEM3D1Format;

  /////////////////////////////////////////////////////////////////
  bool CHEM3D1Format::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename

    return(ReadChem3d(ifs,mol,false,"MM2"));
  }

  ////////////////////////////////////////////////////////////////

  bool CHEM3D1Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    return(WriteChem3d(ofs,mol,"MM2"));
  }

  //***********************************************************
  class CHEM3D2Format : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CHEM3D2Format()
    {
      OBConversion::RegisterFormat("c3d2",this);
    }

    virtual const char* Description() //required
    {
      return
        "Chem3D Cartesian 2 format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  CHEM3D2Format theCHEM3D2Format;

  /////////////////////////////////////////////////////////////////
  bool CHEM3D2Format::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename

    return(CHEM3D1Format::ReadChem3d(ifs,mol,false,"C3D"));
  }

  ////////////////////////////////////////////////////////////////

  bool CHEM3D2Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    return(CHEM3D1Format::WriteChem3d(ofs,mol,"C3D"));
  }

  //****************************************************************

  bool CHEM3D1Format::ReadChem3d(istream &ifs,OBMol &mol,bool mmads,const char *type_key)
  {
    char buffer[BUFF_SIZE];
    int natoms,i;
    char tmp[16],tmp1[16];
    char atomic_type[16];
    double exponent = 0.0;
    double divisor = 1.0;
    double Alpha,Beta,Gamma,A,B,C;
    bool has_fractional = false, has_divisor = false;
    matrix3x3 m;

    vector<string> vs;
    ifs.getline(buffer,BUFF_SIZE);
    tokenize(vs,buffer);

    if (mmads)
      {
        if (vs.empty())
          return(false);
        natoms = atoi((char*)vs[0].c_str());
        if (vs.size() == 2)
          mol.SetTitle(vs[1]);
      }
    else
      {
        switch(vs.size())
          {
          case 7 :
            sscanf(buffer,"%d%lf%lf%lf%lf%lf%lf",
                   &natoms,&Alpha,&Beta,&Gamma,&A,&B,&C);
            m.FillOrth(Alpha,Beta,Gamma,A,B,C);
            has_fractional = true;
            break;
          case 8 :
            sscanf(buffer,"%d%lf%lf%lf%lf%lf%lf%lf",
                   &natoms,&Alpha,&Beta,&Gamma,&A,&B,&C,&exponent);
            m.FillOrth(Alpha,Beta,Gamma,A,B,C);
            has_fractional = true;
            has_divisor = true;
            break;
          default :
            sscanf(buffer,"%d",&natoms);
            break;
          }
      }

    if (!natoms)
      return(false);
    divisor = pow(10.0,exponent);
    mol.ReserveAtoms(natoms);

    ttab.SetToType("INT");
    ttab.SetFromType(type_key);

    OBAtom *atom;
    double x,y,z;
    vector3 v;

    unsigned int k;
    for (i = 1; i <= natoms; i++)
      {
        ifs.getline(buffer,BUFF_SIZE);
        sscanf(buffer,"%15s%*d%lf%lf%lf%15s",
               atomic_type,
               &x,
               &y,
               &z,
               tmp);
        v.Set(x,y,z);
        if (has_fractional)
          v *= m;
        if (has_divisor)
          v/= divisor;

        tokenize(vs,buffer);
        if (vs.empty())
          return(false);

        atom = mol.NewAtom();
        ttab.Translate(tmp1,tmp);
        atom->SetType(tmp1);
        atom->SetVector(v);
        atom->SetAtomicNum(OBElements::GetAtomicNum(atomic_type));

        for (k = 6;k < vs.size(); k++)
          mol.AddBond(atom->GetIdx(),atoi((char*)vs[k].c_str()),1);
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

    mol.PerceiveBondOrders();

    return(true);
  }

  bool CHEM3D1Format::WriteChem3d(ostream &ofs,OBMol &mol, const char *mol_typ)
  {
    int atnum;
    int type_num;
    char buffer[BUFF_SIZE],type_name[16],ele_type[16];

    ofs << mol.NumAtoms();
    if (EQ(mol_typ,"MMADS"))
      {
        ofs << " " << mol.GetTitle();
        ttab.SetToType("MM2");
      }
    else
      ttab.SetToType(mol_typ);
    ofs << endl;

    ttab.SetFromType("INT");

    OBAtom *atom,*nbr;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        if (!ttab.Translate(type_name,atom->GetType()))
          {
            snprintf(buffer, BUFF_SIZE,
                     "Unable to assign %s type to atom %d type = %s\n",
                     mol_typ,atom->GetIdx(),atom->GetType());
            obErrorLog.ThrowError(__FUNCTION__, buffer, obInfo);
            atnum = atom->GetAtomicNum();
            type_num = atnum * 10 + atom->GetExplicitDegree();
            snprintf(type_name, sizeof(type_num), "%d",type_num);
          }
        strncpy(ele_type, OBElements::GetSymbol(atom->GetAtomicNum()), sizeof(ele_type));
        ele_type[sizeof(ele_type) - 1] = '\0';
        snprintf(buffer, BUFF_SIZE, "%-3s %-5d %8.4f  %8.4f  %8.4f %5s",
                ele_type,
                atom->GetIdx(),
                atom->x(),
                atom->y(),
                atom->z(),
                type_name);
        ofs << buffer;

        for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
          {
            snprintf(buffer, BUFF_SIZE, "%6d",nbr->GetIdx());
            ofs << buffer;
          }
        ofs << endl;
      }
    return(true);
  }

} //namespace OpenBabel
