/**********************************************************************
  Copyright (C) 2013 by Casper Steinmann

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

#include <algorithm>

using namespace std;
namespace Dalton
{
}

namespace OpenBabel
{
#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989

  class DALTONInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    DALTONInputFormat()
    {
      OBConversion::RegisterFormat("dalmol",this, "chemical/x-dalton-input");
    }


    virtual const char* Description() //required
    {
      return
        "DALTON Input\n";
    };

    virtual const char* SpecificationURL()
    {return "http://daltonprogram.org/www/resources/dalton2011manual.pdf";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-dalton-input"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return WRITEONEONLY; // | NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };

  //Make an instance of the format class
  DALTONInputFormat theDALTONInputFormat;

  bool DALTONInputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    const char bb = '=';
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    int atomtypes = 0;
    int atomcount = 0; // for each atom type
    int atomcharge = 0; // for each atom type
    double factor = 1.0f;

    mol.BeginModify();
    while(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"ATOMBASIS") != NULL)
      {
        cout << "Cannot read ATOMBASIS format" << endl;
        return(false);
      }
      if(strstr(buffer,"BASIS") != NULL)
      {
        ifs.getline(buffer,BUFF_SIZE); // basis
        ifs.getline(buffer,BUFF_SIZE); // title line 1
        mol.SetTitle(buffer);
        ifs.getline(buffer,BUFF_SIZE); // title line 2

        // reading the first real option
        ifs.getline(buffer,BUFF_SIZE); // options

        // first check if there is a NoSymmetry line. otherwise bail out
        if(strstr(buffer,"NoSymmetry") == NULL)
        {
          cout << "Only molecules with NoSymmetry specified can be read" << endl;
          return(false);
        }

        // if input is in bohr, convert to angstrom
        if(strstr(buffer,"Angstrom") == NULL)
          factor = BOHR_TO_ANGSTROM;

        if(strstr(buffer,"AtomTypes") != NULL)
        {
           tokenize(vs,(strstr(buffer,"AtomTypes=")), " \t\n=");
           atomtypes = atoi(vs[1].c_str());
        }
        while(atomtypes >= 0 && ifs.getline(buffer, BUFF_SIZE))
        {
          if(strstr(buffer, "Atoms") != NULL && strstr(buffer, "Charge") != NULL)
          {
             tokenize(vs,(strstr(buffer,"Atoms=")), " \t\n=");
             atomcount = atoi(vs[1].c_str());
             tokenize(vs,(strstr(buffer,"Charge=")), " \t\n=");
             atomcharge = atoi(vs[1].c_str());
             atomtypes--;
             continue;
          }
          tokenize(vs,buffer);
          if(vs.size() == 4)
          {
            atom = mol.NewAtom();
            atom->SetAtomicNum(atomcharge);
            x = atof((char*)vs[1].c_str()) * factor;
            y = atof((char*)vs[2].c_str()) * factor;
            z = atof((char*)vs[3].c_str()) * factor;
            atom->SetVector(x,y,z);
          }
        }
      }
    }
    // always find bonds and figure out bond orders
    mol.EndModify();
    mol.ConnectTheDots();
    mol.PerceiveBondOrders();
    return(true);
  }

  bool DALTONInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    //   unsigned int i;
    char buffer[BUFF_SIZE];

    // containers for atomtypes and charges
    std::vector<int> groupcounts;
    std::vector<int> groupcharges;

    ofs << "BASIS" << endl;
    ofs << "6-31G*" << endl;
    ofs << mol.GetTitle() << endl;
    ofs << "Generated by Open Babel" << endl;

    // dalton needs some additional information
    //   AtomTypes: the number of atom types you want to include.
    //              it is the number of groups of atoms you want
    //              to read.
    // furthermore, each group has an atom count, so we might as
    // well collect that information here along with the charge.
    int atomtypes = 0;
    int atomtype = -1;
    FOR_ATOMS_OF_MOL(atom, mol)
    {
      if(atom->GetAtomicNum() != atomtype)
      {
        atomtype = atom->GetAtomicNum();
        atomtypes++;
        groupcounts.push_back(0);
        groupcharges.push_back(atom->GetAtomicNum());
      }
      groupcounts[atomtypes-1] += 1;
    }
    ofs << "AtomTypes=" << atomtypes << " NoSymmetry Angstrom" << endl;

    // now let us dump those atoms to the .mol file
    atomtypes = 0;
    atomtype = -1;
    FOR_ATOMS_OF_MOL(atom, mol)
    {
      if(atom->GetAtomicNum() != atomtype)
      {
        atomtype = atom->GetAtomicNum();
        atomtypes++;
        snprintf(buffer, BUFF_SIZE, "Charge=%d.0 Atoms=%i",
                 groupcharges[atomtypes-1],
                 groupcounts[atomtypes-1]);
        ofs << buffer << endl;
      }
      snprintf(buffer, BUFF_SIZE, "%-3s %22.10f  %14.10f  %14.10f ",
               etab.GetSymbol(atom->GetAtomicNum()),
               atom->GetX(),
               atom->GetY(),
               atom->GetZ());
      ofs << buffer << endl;
    }
    return(true);
  }

} //namespace OpenBabel
