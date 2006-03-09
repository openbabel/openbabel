/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include "babelconfig.h"

#include "mol.h"
#include "obconversion.h"
#include "obmolecformat.h"

using namespace std;
namespace OpenBabel
{

class BGFFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    BGFFormat()
    {
        OBConversion::RegisterFormat("bgf",this);
    }

    virtual const char* Description() //required
    {
        return
            "MSI BGF format\n \
            No comments yet\n";
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
BGFFormat theBGFFormat;

/////////////////////////////////////////////////////////////////
bool BGFFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename
    mol.BeginModify();

    char buffer[BUFF_SIZE];
    char tmp[15],tmptyp[15];

    while (ifs.getline(buffer,BUFF_SIZE))
        if (EQn(buffer,"FORMAT",6))
            break;

    ttab.SetFromType("DRE");
    ttab.SetToType("INT");
    OBAtom *atom;
    double x,y,z,chrg;
    for (;;)
    {
        if (!ifs.getline(buffer,BUFF_SIZE))
            break;
        if (EQn(buffer,"FORMAT",6))
            break;

        sscanf(buffer,"%*s %*s %*s %*s %*s %*s %lf %lf %lf %s %*s %*s %lf",
               &x,&y,&z,
               tmptyp,
               &chrg);
        atom = mol.NewAtom();

        ttab.Translate(tmp,tmptyp);
        atom->SetType(tmp);

	CleanAtomType(tmptyp);
        atom->SetAtomicNum(etab.GetAtomicNum(tmptyp));

        atom->SetVector(x,y,z);
    }
    unsigned int i;
    vector<int> vtmp;
    vector<vector<int> > vcon;
    vector<vector<int> > vord;

    for (i = 0; i < mol.NumAtoms();i++)
    {
        vcon.push_back(vtmp);
        vord.push_back(vtmp);
    }

    unsigned int bgn;
    vector<string> vs;
    for (;;)
    {
        if (!ifs.getline(buffer,BUFF_SIZE) || EQn(buffer,"END",3))
            break;

        tokenize(vs,buffer);
        if (vs.empty() || vs.size() < 3 || vs.size() > 10)
            continue;

        if (EQn(buffer,"CONECT",6))
        {
            bgn = atoi((char*)vs[1].c_str()) - 1;
            if (bgn < 1 || bgn > mol.NumAtoms())
                continue;
            for (i = 2;i < vs.size();i++)
            {
                vcon[bgn].push_back(atoi((char*)vs[i].c_str()));
                vord[bgn].push_back(1);
            }
        }
        else
            if (EQn(buffer,"ORDER",5))
            {
                bgn = atoi((char*)vs[1].c_str()) - 1;
                if (bgn < 1 || bgn > mol.NumAtoms())
                    continue;
                if (vs.size() > vord[bgn].size()+2)
                    continue;
                for (i = 2;i < vs.size();i++)
                    vord[bgn][i-2] = atoi((char*)vs[i].c_str());
            }
    }

    unsigned int j;
    for (i = 1;i <= mol.NumAtoms();i++)
        if (!vcon[i - 1].empty())
            for (j = 0;j < vcon[i - 1].size();j++)
            {
                mol.AddBond(i,vcon[i - 1][j],vord[i - 1][j]);
            }

    //load up the next line after the END marker
    ifs.getline(buffer,BUFF_SIZE);

    mol.EndModify();
    return(true);
}

////////////////////////////////////////////////////////////////

bool BGFFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    vector<OBNodeBase*>::iterator i;
    int max_val;
    OBAtom *atom;
    char buffer[BUFF_SIZE];
    char elmnt_typ[5], dreid_typ[5], atm_sym[10], max_val_str[5];

    mol.Kekulize();

    ofs << "BIOGRF 200" << endl;
    sprintf(buffer,"DESCRP %s",mol.GetTitle());
    ofs << buffer << endl;
    sprintf(buffer,"REMARK BGF file created by Open Babel %s",BABEL_VERSION);
    ofs << buffer << endl;
    ofs << "FORCEFIELD DREIDING  " << endl;
    ofs << "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)" << endl;

    ttab.SetFromType("INT");

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
        strcpy(elmnt_typ,etab.GetSymbol(atom->GetAtomicNum()));
        ToUpper(elmnt_typ);

        ttab.SetToType("DRE");
        ttab.Translate(dreid_typ,atom->GetType());
        ttab.SetToType("HAD");
        ttab.Translate(max_val_str,atom->GetType());
        max_val = atoi(max_val_str);
        if (max_val == 0)
            max_val = 1;
        sprintf(atm_sym,"%s%d",elmnt_typ,atom->GetIdx());
        sprintf(buffer,"%6s %5d %-5s %3s %1s %5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f",
                "HETATM",
                atom->GetIdx(),
                atm_sym,
                "RES",
                "A",
                "444",
                atom->GetX(),
                atom->GetY(),
                atom->GetZ(),
                dreid_typ,
                max_val,
                0,
                atom->GetPartialCharge());
        ofs << buffer << endl;
    }
    sprintf(buffer,"FORMAT CONECT (a6,12i6)\n");
    ofs << buffer << endl;

    OBAtom *nbr;
    vector<OBEdgeBase*>::iterator j;
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
        if (atom->GetValence())
        {
            sprintf(buffer,"CONECT%6d",atom->GetIdx());
            ofs << buffer;
            for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
            {
                sprintf(buffer,"%6d",nbr->GetIdx());
                ofs << buffer;
            }
            ofs << endl;

            sprintf(buffer,"ORDER %6d",atom->GetIdx());
            ofs << buffer;
            for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
            {
                sprintf(buffer,"%6d",(*j)->GetBO());
                ofs << buffer;
            }
            ofs << endl;
        }

    sprintf(buffer,"END");
    ofs << buffer << endl;
    return(true);
}

} //namespace OpenBabel
