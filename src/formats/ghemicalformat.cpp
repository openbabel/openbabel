/**********************************************************************
  Copyright (C) 2000-2006 by Geoffrey Hutchison
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
#include <openbabel/obiter.h>
#include <openbabel/generic.h>
#include <openbabel/bond.h>
#include <cstdlib>


using namespace std;
namespace OpenBabel
{

  class GhemicalFormat : public OBMoleculeFormat
  {
    public:
      //Register this format type ID
      GhemicalFormat()
      {
        //        OBConversion::RegisterFormat("mm1gp",this);
        //        OBConversion::RegisterFormat("qm1gp",this);
        OBConversion::RegisterFormat("gpr",this);
      }

      virtual const char* Description() //required
      {
        return
          "Ghemical format\n"
          "Open source molecular modelling\n";
      };

      virtual const char* SpecificationURL()
      { return "http://www.uku.fi/~thassine/ghemical/"; }; //optional

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
  GhemicalFormat theGhemicalFormat;

  /////////////////////////////////////////////////////////////////
  bool GhemicalFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int i;
    int natoms, nbonds;
    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    char bobuf[100];
    string bostr;
    int bgn,end,order;
    bool hasPartialCharges = false;

    mol.BeginModify();

    // Get !Header line with version number
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%*s %*s %d", &i);
    if (!i)
      return false;

    // Get !Info line with number of coord sets
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%*s %d", &i);
    if (!i)
      return false;

    // Get !Atoms line with number
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%*s %d", &natoms);
    if (!natoms)
      return(false);

    for (i = 1; i <= natoms; i ++)
    {
      if (!ifs.getline(buffer,BUFF_SIZE))
        return(false);
      tokenize(vs,buffer);
      if (vs.size() < 2)
        return(false);
      atom = mol.NewAtom();
      atom->SetAtomicNum(atoi(vs[1].c_str()));
    }

    // Get !Bonds line with number
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%*s %d", &nbonds);
    if (nbonds != 0)
      for (i = 0; i < nbonds; i++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        if (!sscanf(buffer,"%d%d%2s",&bgn,&end,bobuf))
          return (false);
        bostr = bobuf;
        order = 1;
        if      (bostr == "D")
          order = 2;
        else if (bostr == "T")
          order = 3;
        else if (bostr == "C")
          order = 5; // Conjugated ~= Aromatic
        mol.AddBond(bgn+1,end+1,order);
      }

    // Get !Coord line
    ifs.getline(buffer,BUFF_SIZE);
    for (i = 1; i <= natoms; i ++)
    {
      if (!ifs.getline(buffer,BUFF_SIZE))
        return(false);
      tokenize(vs,buffer);
      if (vs.size() != 4)
        return(false);
      atom = mol.GetAtom(i);
      x = 10.0*atof((char*)vs[1].c_str());
      y = 10.0*atof((char*)vs[2].c_str());
      z = 10.0*atof((char*)vs[3].c_str());
      atom->SetVector(x,y,z); //set coordinates
    }

    if (ifs.getline(buffer,BUFF_SIZE) && (strstr(buffer, "!Charges") != NULL
          || strstr(buffer, "!PartialCharges") != NULL))
    {
      hasPartialCharges = true;
      for (i = 1; i <= natoms; i ++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        tokenize(vs,buffer);
        if (vs.size() != 2)
          return(false);
        atom = mol.GetAtom(i);
        atom->SetPartialCharge(atof((char*)vs[1].c_str()));
      }
    }

    // look for the !End block if it exists
    while (ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"!End") != NULL)
        break;
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
    if (hasPartialCharges)
      mol.SetPartialChargesPerceived();
    mol.SetTitle(title);
    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool GhemicalFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    // delete dummy atoms
    FOR_ATOMS_OF_MOL(atom, pmol)
      if (atom->GetAtomicNum() == 0)
        mol.DeleteAtom(&*atom);

    // Ghemical header -- here "version 1.0" format
    ofs << "!Header gpr 100\n";

    // Number of coordinate sets
    ofs << "!Info 1\n";

    // Atom definitions
    ofs << "!Atoms " << mol.NumAtoms() << '\n';

    FOR_ATOMS_OF_MOL(atom, mol)
      ofs << (atom->GetIdx() - 1) << " " << atom->GetAtomicNum() << '\n';

    // Bond definitions
    ofs << "!Bonds " << mol.NumBonds() << '\n';

    char bond_char;
    FOR_BONDS_OF_MOL(bond, mol)
    {
      switch(bond->GetBondOrder())
      {
        case 1 :
          bond_char = 'S';
          break;
        case 2 :
          bond_char = 'D';
          break;
        case 3 :
          bond_char = 'T';
          break;
        case 5 :
          bond_char = 'C';
          break;
        default :
          bond_char = 'S';
      }
      if (bond->IsAromatic())
        bond_char = 'C';
      ofs << bond->GetBeginAtomIdx()-1 << ' '
        << bond->GetEndAtomIdx()-1 << ' '
        <<  bond_char << '\n';
    }

    // Coordinate sets (here only 1)
    ofs << "!Coord\n";

    FOR_ATOMS_OF_MOL(atom, mol)
    {
      ofs << atom->GetIdx()-1 << ' '
        << atom->GetX()/10.0 << ' '
        << atom->GetY()/10.0 << ' '
        << atom->GetZ()/10.0 << '\n';
    }

    // Calculated charges
    ofs << "!Charges\n";

    FOR_ATOMS_OF_MOL(atom, mol)
    {
      ofs << atom->GetIdx()-1 << ' '
        << atom->GetPartialCharge() << '\n';
    }

    OBSetData *gmsset = (OBSetData *)pmol->GetData("gamess");
    if(gmsset)
    {
      ofs << "!GAMESS" << endl;
      std::vector<OBGenericData*>::iterator i,j;

      for(i = gmsset->GetBegin(); i != gmsset->GetEnd(); ++i)
      {
        OBSetData *cset = (OBSetData *)(*i);
        if(cset)
        {
          string section = cset->GetAttribute();
          for(j = cset->GetBegin(); j != cset->GetEnd(); ++j)
          {
            OBPairData *pd = (OBPairData *) (*j);
            if(pd)
            {
              ofs << section << " " << pd->GetAttribute() << " " << pd->GetValue() << endl;
            }
          }
        }
      }
    }

    ofs << "!End\n";

    return(true);
  }

} //namespace OpenBabel
