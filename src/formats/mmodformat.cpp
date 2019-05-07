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
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/generic.h>
#include <cstdlib>


using namespace std;
namespace OpenBabel
{

  class MacroModFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MacroModFormat()
    {
      OBConversion::RegisterFormat("mmd",this, "chemical/x-macromodel-input");
      OBConversion::RegisterFormat("mmod",this, "chemical/x-macromodel-input");
    }

    virtual const char* Description() //required
    {
      return
        "MacroModel format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-macromodel-input"; };

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
  MacroModFormat theMacroModFormat;

  /////////////////////////////////////////////////////////////////
  bool MacroModFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* defaultTitle = pConv->GetTitle();

    // Get Title
    char buffer[BUFF_SIZE];
    int natoms;
    vector<vector<pair<int,int> > > connections;

    if (ifs.getline(buffer,BUFF_SIZE))
      {
        vector<string> vs;
        tokenize(vs,buffer," \n");

        if ( !vs.empty() && vs.size() > 0)
          sscanf(buffer,"%i%*s",&natoms);

        if (natoms == 0)
          return false;

        if ( !vs.empty() && vs.size() > 1)
          mol.SetTitle(vs[1]);
        else
          {
            string s = defaultTitle;
            mol.SetTitle(defaultTitle);
          }
      }
    else
      return(false);

    mol.BeginModify();
    mol.ReserveAtoms(natoms);
    connections.resize(natoms+1);

    /***********************************************************************/

    // Get Type Bonds, BondOrder, X, Y, Z

    double x,y,z;
    vector3 v;
    char temp_type[10];
    int i,j;
    double charge;
    OBAtom atom;

    ttab.SetFromType("MMD");
    for (i = 1; i <= natoms; i++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          break;

        int end[6], order[6];

        sscanf(buffer,"%9s%d%d%d%d%d%d%d%d%d%d%d%d%lf%lf%lf",
               temp_type,&end[0],&order[0],&end[1],&order[1],&end[2],&order[2],
               &end[3], &order[3], &end[4], &order[4], &end[5], &order[5],
               &x, &y, &z);

        pair<int,int> tmp;
        for ( j = 0 ; j <=5 ; j++ )
          {
            if ( end[j] > 0  && end[j] > i)
              {
                tmp.first = end[j];
                tmp.second = order[j];
                connections[i].push_back(tmp);
              }
          }

        v.SetX(x);
        v.SetY(y);
        v.SetZ(z);
        atom.SetVector(v);

        string str = temp_type,str1;
        ttab.SetToType("ATN");
        ttab.Translate(str1,str);
        atom.SetAtomicNum(atoi(str1.c_str()));
        ttab.SetToType("INT");
        ttab.Translate(str1,str);
        atom.SetType(str1);

        // stuff for optional fields

        buffer[109]='\0';
        sscanf(&buffer[101],"%lf", &charge);
        atom.SetPartialCharge(charge);
        mol.AddAtom(atom);
      }

    for (i = 1; i <= natoms; i++)
      for (j = 0; j < (signed)connections[i].size(); j++)
        mol.AddBond(i, connections[i][j].first, connections[i][j].second);

    mol.EndModify();
    mol.SetPartialChargesPerceived();
    // Annotate origin of partial charges
    OBPairData *dp = new OBPairData;
    dp->SetAttribute("PartialCharges");
    dp->SetValue("MACROMODEL");
    dp->SetOrigin(fileformatInput);
    mol.SetData(dp);

    OBBond *bond;
    vector<OBBond*>::iterator bi;
    for (bond = mol.BeginBond(bi);bond;bond = mol.NextBond(bi))
      if (bond->GetBondOrder() == 5 && !bond->IsInRing())
        bond->SetBondOrder(1);

    if ( natoms != (signed)mol.NumAtoms() )
      return(false);

    // clean out remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool MacroModFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    snprintf(buffer, BUFF_SIZE, " %5d %6s      E = %7.3f KJ/mol",
            mol.NumAtoms(),mol.GetTitle(),4.184*mol.GetEnergy());
    ofs << buffer << endl;

    int type,k;
    OBAtom *atom,*nbr;
    string from,to;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;
    ttab.SetFromType("INT");
    ttab.SetToType("MMD");

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i)) {
        if (atom->GetAtomicNum() == OBElements::Hydrogen) {
          type = 41;
          if ((nbr = atom->BeginNbrAtom(j))) {
            if (nbr->GetAtomicNum() == OBElements::Oxygen) {
              type = 42;
            }
            else if (nbr->GetAtomicNum() == OBElements::Nitrogen) {
              type = 43;
            }
          }
        } else {
            from = atom->GetType();
            ttab.Translate(to,from);
            type = atoi((char*)to.c_str());
        }
        snprintf(buffer, BUFF_SIZE, "%4d",type);
        ofs << buffer;
        for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
          {
            snprintf(buffer, BUFF_SIZE, " %5d %1d",nbr->GetIdx(),(*j)->GetBondOrder());
            ofs << buffer;
          }
        for (k=atom->GetExplicitDegree();k < 6;k++)
          {
            snprintf(buffer, BUFF_SIZE," %5d %1d",0,0);
            ofs << buffer;
          }

        snprintf(buffer, BUFF_SIZE, " %11.6f %11.6f %11.6f %5d %5d %8.5f \n",
                atom->x(), atom->y(),atom->z(),0,0,
                atom->GetPartialCharge());
        ofs << buffer;
      }

    return(true);
  }

} //namespace OpenBabel
