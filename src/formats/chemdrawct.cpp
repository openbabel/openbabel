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

using namespace std;
namespace OpenBabel
{

  class ChemDrawFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ChemDrawFormat()
    {
      OBConversion::RegisterFormat("ct",this);
    }

    virtual const char* Description() //required
    {
      return
        "ChemDraw Connection Table format\n"
        "No comments yet\n";
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
  ChemDrawFormat theChemDrawFormat;

  ////////////////////////////////////////////////////////////////

  bool ChemDrawFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];

    ofs << mol.GetTitle() << endl;
    ofs << " " << mol.NumAtoms() << " " << mol.NumBonds() << endl;

    OBAtom *atom;
    vector<OBAtom*>::iterator i;

    for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        snprintf(buffer, BUFF_SIZE, " %9.4f %9.4f    0.0000 %-1s",
                 atom->x(),
                 atom->y(),
                 etab.GetSymbol(atom->GetAtomicNum()));
        ofs << buffer << endl;
      }

    OBBond *bond;
    vector<OBBond*>::iterator j;

    for(bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
      {
        snprintf(buffer, BUFF_SIZE, "%3d%3d%3d%3d",
                 bond->GetBeginAtomIdx(),
                 bond->GetEndAtomIdx(),
                 bond->GetBO(), bond->GetBO());
        ofs << buffer << endl;
      }
    return(true);
  }

  bool ChemDrawFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    unsigned int natoms, nbonds;
    unsigned int bgn, end, order;
    vector<string> vs;
    OBAtom *atom;
    double x, y, z;

    mol.SetDimension(2);
    mol.BeginModify();

    ifs.getline(buffer,BUFF_SIZE);
    if (strlen(buffer) == 0)
      mol.SetTitle(buffer);
    else
      mol.SetTitle(title);

    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer," %d %d", &natoms, &nbonds);

    for (unsigned int i = 1; i <= natoms; i ++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
        tokenize(vs,buffer);
        if (vs.size() != 4) return(false);
        atom = mol.NewAtom();

        x = atof((char*)vs[0].c_str());
        y = atof((char*)vs[1].c_str());
        z = atof((char*)vs[2].c_str());

        atom->SetVector(x,y,z); //set coordinates
        atom->SetAtomicNum(etab.GetAtomicNum(vs[3].c_str()));
      }

    if (nbonds != 0)
      for (unsigned int i = 0; i < nbonds; i++)
        {
          if (!ifs.getline(buffer,BUFF_SIZE)) return(false);
          tokenize(vs,buffer);
          if (vs.size() != 4) return(false);
          if (!sscanf(buffer,"%d%d%d%*d",&bgn,&end,&order)) return (false);
          mol.AddBond(bgn,end,order);
        }

    // clean out remaining blank lines
    while(ifs.peek() != EOF && ifs.good() &&
          (ifs.peek() == '\n' || ifs.peek() == '\r'))
      ifs.getline(buffer,BUFF_SIZE);

    mol.EndModify();
    return(true);
  }


} //namespace OpenBabel
