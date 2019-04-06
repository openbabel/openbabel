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
#include <openbabel/generic.h>
#include <openbabel/internalcoord.h>
#include <openbabel/math/matrix3x3.h>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class CacaoFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CacaoFormat()
    {
      OBConversion::RegisterFormat("caccrt",this);
    }

    virtual const char* Description() //required
    {
      return
        "Cacao Cartesian format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.chembio.uoguelph.ca/oakley/310/cacao/cacao.htm";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | WRITEONEONLY;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    static void SetHilderbrandt(OBMol&,vector<OBInternalCoord*>&);
  };

  //Make an instance of the format class
  CacaoFormat theCacaoFormat;

  /////////////////////////////////////////////////////////////////
  bool CacaoFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename

    char buffer[BUFF_SIZE];
    int natoms;
    double A,B,C,Alpha,Beta,Gamma;
    matrix3x3 m;

    ifs.getline(buffer,BUFF_SIZE);
    mol.SetTitle(buffer);
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%d",&natoms);

    while (ifs.getline(buffer,BUFF_SIZE) &&!EQn(buffer,"CELL",4))
      ;

    if (!EQn(buffer,"CELL",4))
      return(false);
    vector<string> vs;
    tokenize(vs,buffer," \n\t,");
    if (vs.size() != 7)
      return(false);

    //parse cell values
    A = atof((char*)vs[1].c_str());
    B = atof((char*)vs[2].c_str());
    C = atof((char*)vs[3].c_str());
    Alpha = atof((char*)vs[4].c_str());
    Beta  = atof((char*)vs[5].c_str());
    Gamma = atof((char*)vs[6].c_str());

    OBUnitCell *uc = new OBUnitCell;
    uc->SetData(A, B, C, Alpha, Beta, Gamma);
    uc->SetOrigin(fileformatInput);
    mol.SetData(uc);

    int i;
    double x,y,z;
    OBAtom *atom;
    vector3 v;

    mol.BeginModify();

    for (i = 1; i <= natoms;i++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        tokenize(vs,buffer," \n\t,");
        if (vs.size() < 4)
          return(false);
        atom = mol.NewAtom();

        x = atof((char*)vs[1].c_str());
        y = atof((char*)vs[2].c_str());
        z = atof((char*)vs[3].c_str());
        v.Set(x,y,z);
        v = uc->FractionalToCartesian(v);

        atom->SetAtomicNum(OBElements::GetAtomicNum(vs[0].c_str()));
        atom->SetVector(v);
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
    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool CacaoFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    OBAtom *atom;
    char buffer[BUFF_SIZE];
    vector<OBAtom*>::iterator i;

    snprintf(buffer, BUFF_SIZE, "%s\n",mol.GetTitle());
    ofs << buffer;
    snprintf(buffer, BUFF_SIZE, "%3d   DIST  0  0  0\n",mol.NumAtoms());
    ofs << buffer;

    if (!mol.HasData(OBGenericDataType::UnitCell))
      ofs << "CELL 1.,1.,1.,90.,90.,90.\n";
    else
      {
        OBUnitCell *uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
        snprintf(buffer, BUFF_SIZE, "CELL %f,%f,%f,%f,%f,%f\n",
                 uc->GetA(), uc->GetB(), uc->GetC(),
                 uc->GetAlpha(), uc->GetBeta(), uc->GetGamma());
        ofs << buffer;
      }

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        snprintf(buffer,BUFF_SIZE,"%2s %7.4f, %7.4f, %7.4f\n",
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 atom->x(),
                 atom->y(),
                 atom->z());
        ofs << buffer;
      }

    return(true);
  }

  //! \todo Make this method bulletproof. Currently it causes segfaults sometimes
  void CacaoFormat::SetHilderbrandt(OBMol &mol,vector<OBInternalCoord*> &vit)
  {
    // Roundtrip testing shows that some atoms are NULL
    //  which causes segfaults when dereferencing later
    //   (e.g. in the last "segment" of this routine)
    double sum,r;

    OBAtom dummy1,dummy2;
    dummy1.SetVector(0.0,0.0,1.0);
    dummy2.SetVector(1.0,0.0,0.0);

    OBAtom *atom,*a1,*a2,*ref;
    vector<OBAtom*>::iterator ai;

    vit.push_back((OBInternalCoord*)NULL);
    for (atom = mol.BeginAtom(ai);atom;atom = mol.NextAtom(ai))
      vit.push_back(new OBInternalCoord (atom));

    vit[1]->_a = &dummy1;
    vit[1]->_b = &dummy2;
    if (vit.size() > 2) {
      vit[2]->_b = &dummy1;
      vit[2]->_c = &dummy2;
      if  (vit.size() > 3) {
        vit[3]->_c = &dummy1;
      }
    }

    unsigned int i,j;
    for (i = 2;i <= mol.NumAtoms();i++)
      {
        ref = (OBAtom*)NULL;
        a1 = mol.GetAtom(i);
        sum = 100.0;
        for (j = 1;j < i;j++)
          {
            a2 = mol.GetAtom(j);
            r = (a1->GetVector()-a2->GetVector()).length_2();
            if ((r < sum) && (vit[j]->_a != a2) && (vit[j]->_b != a2))
              {
                sum = r;
                ref = a2;
              }
          }
        vit[i]->_a = ref;
      }

    for (i = 3;i <= mol.NumAtoms();i++)
      vit[i]->_b = vit[vit[i]->_a->GetIdx()]->_a;

    for (i = 4;i <= mol.NumAtoms();i++)
      {
        if (vit[i]->_b && vit[i]->_b->GetIdx())
          vit[i]->_c = vit[vit[i]->_b->GetIdx()]->_b;
        else
          vit[i]->_c = &dummy1;
      }

    OBAtom *a,*b,*c;
    vector3 v1,v2;
    for (i = 2;i <= mol.NumAtoms();i++)
      {
        atom = mol.GetAtom(i);
        a = vit[i]->_a;
        b = vit[i]->_b;
        c = vit[i]->_c;
        v1 = atom->GetVector() - a->GetVector();
        v2 = b->GetVector() - a->GetVector();
        vit[i]->_ang = vectorAngle(v1,v2);
        vit[i]->_tor = CalcTorsionAngle(atom->GetVector(),
                                        a->GetVector(),
                                        b->GetVector(),
                                        c->GetVector());
        vit[i]->_dst = (vit[i]->_a->GetVector() - atom->GetVector()).length();
      }

  }
  //***************************************************************
  class CacaoInternalFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CacaoInternalFormat()
    {
      OBConversion::RegisterFormat("cacint",this);
    }

    virtual const char* Description() //required
    {
      return
        "Cacao Internal format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL(){return
        "http://www.chembio.uoguelph.ca/oakley/310/cacao/cacao.htm";}; //optional

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
  CacaoInternalFormat theCacaoInternalFormat;

  bool CacaoInternalFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    vector3 v;
    char tmptype[16],buffer[BUFF_SIZE];

    if (mol.Empty())
      return(false);

    //translate first atom to origin
    v = mol.GetAtom(1)->GetVector();
    v *= -1.0;
    mol.Translate(v);

    vector<OBInternalCoord*> vit;
    CacaoFormat::SetHilderbrandt(mol,vit);
    strncpy(tmptype,OBElements::GetSymbol(mol.GetAtom(1)->GetAtomicNum()), sizeof(tmptype));
    tmptype[sizeof(tmptype) - 1] = '\0';

    ofs << " # TITLE\n";
    snprintf(buffer, BUFF_SIZE, "%3d  0DIST  0  0  0\n",mol.NumAtoms());
    ofs << "  EL\n";
    snprintf(buffer, BUFF_SIZE, "0.,0.,0., %s\n",tmptype);
    ofs << buffer;
    for (i = 2; i <= mol.NumAtoms(); i++)
      {
        strncpy(tmptype,OBElements::GetSymbol(mol.GetAtom(i)->GetAtomicNum()), sizeof(tmptype));
        tmptype[sizeof(tmptype) - 1] = '\0';

        if (vit[i]->_tor < 0.0)
          vit[i]->_tor += 360.0;
        snprintf(buffer, BUFF_SIZE, "%2d,%d,%2s%7.3f,%7.3f,%7.3f",
                 vit[i]->_a->GetIdx(),i,tmptype,
                 vit[i]->_dst,
                 vit[i]->_ang,
                 vit[i]->_tor);
        ofs << buffer << endl;
      }

    vector<OBInternalCoord*>::iterator j;
    for (j = vit.begin();j != vit.end();++j)
      if (*j)
        {
          delete *j;
          *j = NULL;
        }

    return(true);
  }

} //namespace OpenBabel
