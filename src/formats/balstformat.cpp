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
#include <cstdlib>


using namespace std;
namespace OpenBabel
{

  class BallStickFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    BallStickFormat()
    {
      OBConversion::RegisterFormat("bs",this);
    }

    virtual const char* Description() //required
    {
      return
        "Ball and Stick format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.orc.uni-linz.ac.at/mueller/ball_and_stick.shtml"; };

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
  BallStickFormat theBallStickFormat;

  /////////////////////////////////////////////////////////////////
  bool BallStickFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int i,natoms;
    char buffer[BUFF_SIZE];

    if (!ifs.getline(buffer,BUFF_SIZE))
      return(false);
    if (!ifs.getline(buffer,BUFF_SIZE))
      return(false);
    sscanf(buffer,"%d",&natoms);
    mol.ReserveAtoms(natoms);
    mol.BeginModify();

    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    vector<string>::iterator j;

    for (i = 1; i <= natoms;i ++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        tokenize(vs,buffer);
        if (vs.size() < 4)
          return(false);
        if (vs[0].size() > 1)
          vs[0][1] = tolower(vs[0][1]);
        atom = mol.NewAtom();
        x = atof((char*)vs[1].c_str());
        y = atof((char*)vs[2].c_str());
        z = atof((char*)vs[3].c_str());
        atom->SetVector(x,y,z); //set coordinates
        atom->SetAtomicNum(OBElements::GetAtomicNum(vs[0].c_str()));

        for (j = vs.begin()+4;j != vs.end();++j)
          mol.AddBond(atom->GetIdx(),atoi((char*)j->c_str()),1);
      }

    // clean out any remaining blank lines
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

  bool BallStickFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char tmptype[16];
    char buffer[BUFF_SIZE];

    if (strlen(mol.GetTitle()) > 0)
      ofs << mol.GetTitle() << endl;
    else
      ofs << "Untitled" << endl;

    snprintf(buffer,BUFF_SIZE,"%d",mol.NumAtoms());
    ofs << buffer << endl;

    OBAtom *atom,*nbr;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;

    for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        strncpy(tmptype,OBElements::GetSymbol(atom->GetAtomicNum()), sizeof(tmptype));
        tmptype[15] = '\0';
        if (strlen(tmptype) > 1)
          tmptype[1] = toupper(tmptype[1]);
        snprintf(buffer,BUFF_SIZE,"%-3s %8.4f  %8.4f  %8.4f",
                 tmptype,
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ());
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
