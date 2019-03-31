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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>


#include <algorithm>

using namespace std;
namespace Dalton
{
}

namespace OpenBabel
{
#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989
  class DALTONOutputFormat : public OBMoleculeFormat
  {
  public:

    DALTONOutputFormat()
    {
      OBConversion::RegisterFormat("dallog", this, "chemical/x-dalton-output");
    }

    virtual const char* Description()
    {
      return
        "DALTON output format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n";
    };

    virtual const char* SpecificationURL()
    {return "http://daltonprogram.org/www/resources/dalton2016manual.pdf";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-dalton-output"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

  private:
  };

  DALTONOutputFormat theDALTONOutputFormat;

  class DALTONInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    DALTONInputFormat()
    {
      OBConversion::RegisterFormat("dalmol",this, "chemical/x-dalton-input");
      OBConversion::RegisterOptionParam("a", NULL, 0, OBConversion::OUTOPTIONS); // write atomic units
      OBConversion::RegisterOptionParam("b", NULL, 0, OBConversion::OUTOPTIONS); // write atombasis format
      OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS); // specify basis set in .mol file
    }


    virtual const char* Description() //required
    {
      return
        "DALTON input format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n"

        "Write Options e.g. -xa\n"
        "  a                write input in atomic units instead of Angstrom\n"
        "  b                write input using the ATOMBASIS format\n"
        "  k <basis>        specify basis set to use\n"
        "                     e.g. ``-xk STO-3G``\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://daltonprogram.org/www/resources/dalton2016manual.pdf";}; //optional

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
  private:
    enum BasisFormat_t {BASIS, ATOMBASIS, INTGRL};
    BasisFormat_t basisformat;
  };

  //Make an instance of the format class
  DALTONInputFormat theDALTONInputFormat;


  // ---
  // INPUT READ
  // ---
  bool DALTONInputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    int molcharge = 0; // overall molecular charge.
    int atomtypes = 0;
    int atomcount = 0; // for each atom type
    int atomcharge = 0; // for each atom type. nuclear charge
    double factor = 1.0f;

    basisformat = BASIS; // Always assume BASIS format. We change it below otherwise

    mol.BeginModify();
    while(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"INTGRL") != NULL)
      {
        basisformat = INTGRL;
        cout << "Cannot read INTGRL format" << endl;
        return(false);
      }
      if(strstr(buffer,"ATOMBASIS") != NULL)
      {
        basisformat = ATOMBASIS;
      }

      // parse the file
      if (basisformat == BASIS) {
        ifs.getline(buffer,BUFF_SIZE); // read basis set line
      }
      // There are always two title lines in the DALTON .mol file format
      ifs.getline(buffer,BUFF_SIZE); // title line 1
      mol.SetTitle(buffer);
      ifs.getline(buffer,BUFF_SIZE); // title line 2. Ignore.

      // reading the first real option
      ifs.getline(buffer,BUFF_SIZE); // options

      // first check if there are any atoms specified. otherwise bail out
      if(strstr(buffer,"AtomTypes") != NULL)
      {
        tokenize(vs,(strstr(buffer,"AtomTypes=")), " \t\n=");
        atomtypes = atoi(vs[1].c_str());
      }
      else
      {
        cout << "AtomTypes not specified in file." << endl;
        return(false);
      }

      // then check if there is a NoSymmetry line. otherwise bail out
      if(strstr(buffer,"NoSymmetry") == NULL)
      {
        cout << "Only molecules with NoSymmetry can be read" << endl;
        return(false);
      }

      if(strstr(buffer,"Charge") != NULL)
      {
        tokenize(vs,(strstr(buffer,"Charge=")), " \t\n=");
        molcharge = atoi(vs[1].c_str());
      }

      // if input is in bohr, convert to angstrom
      if(strstr(buffer,"Angstrom") == NULL)
        factor = BOHR_TO_ANGSTROM;

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
        if(strstr(buffer, "ZMAT") != NULL)
        {
           cout << "ZMAT format not supported" << endl;
           return(false);
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

    mol.EndModify();
    if (!pConv->IsOption("b",OBConversion::INOPTIONS)) {
      mol.ConnectTheDots();
    }

    if (!pConv->IsOption("s",OBConversion::INOPTIONS)
        && !pConv->IsOption("b",OBConversion::INOPTIONS)) {
      mol.PerceiveBondOrders();
    }

    // Set properties that would otherwise be overwritten
    // by EndModify
    mol.SetTotalCharge(molcharge);

    return(true);
  }

  // ---
  // INPUT WRITE
  // ---
  bool DALTONInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    basisformat = BASIS; // Always assume BASIS format. We change it below otherwise

    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    double factor = 1.0f;
    bool writeatomicunit = pConv->IsOption("a", OBConversion::OUTOPTIONS) != NULL;
    const char *keywords = pConv->IsOption("k",OBConversion::OUTOPTIONS);
    string atombasis_str = "";
    string basisset = "6-31G*";

    if (pConv->IsOption("b",OBConversion::OUTOPTIONS)) {
        basisformat = ATOMBASIS;
    }

    if (keywords) {
      basisset = keywords;
    }

    if (writeatomicunit)
      factor *= ANGSTROM_TO_BOHR;

    // containers for atomtypes and charges
    std::vector<int> groupcounts;
    std::vector<int> groupcharges;

    if (basisformat == ATOMBASIS) {
      ofs << "ATOMBASIS" << endl;
      atombasis_str = " Basis=" + basisset;
    } else {
      ofs << "BASIS" << endl;
      ofs << basisset << endl;
    }
    ofs << mol.GetTitle() << endl;
    ofs << "Generated by Open Babel. Check overall charge below." << endl;

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
    ofs << "AtomTypes=" << atomtypes;
    ofs << " Charge=" << mol.GetTotalCharge();
    ofs << " NoSymmetry";
    if (!writeatomicunit) ofs << " Angstrom";
    ofs << endl;

    // now let us dump those atoms to the .mol file
    atomtypes = 0;
    atomtype = -1;
    FOR_ATOMS_OF_MOL(atom, mol)
    {
      if(atom->GetAtomicNum() != atomtype)
      {
        atomtype = atom->GetAtomicNum();
        atomtypes++;
        snprintf(buffer, BUFF_SIZE, "Charge=%d.0 Atoms=%i%s",
                 groupcharges[atomtypes-1],
                 groupcounts[atomtypes-1],
                 atombasis_str.c_str());
        ofs << buffer << endl;
      }
      snprintf(buffer, BUFF_SIZE, "%-3s %22.10f  %14.10f  %14.10f ",
               OBElements::GetSymbol(atom->GetAtomicNum()),
               atom->GetX()*factor,
               atom->GetY()*factor,
               atom->GetZ()*factor);
      ofs << buffer << endl;
    }
    return(true);
  }

  // ---
  // OUTPUT READ
  // ---
  bool DALTONOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    int atomcount = 0;

    mol.BeginModify();
    while(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"Cartesian Coordinates (a.u.)") != NULL)
      {
        cout << "Reading coordinates." << endl;
        ifs.getline(buffer,BUFF_SIZE); // ---
        ifs.getline(buffer,BUFF_SIZE); // whitespace
        ifs.getline(buffer,BUFF_SIZE); // number of coordinates
        tokenize(vs,buffer);
        atomcount = atoi(vs[4].c_str()) / 3; // number of atoms to read
        while(atomcount > 0)
        {
          atomcount --;
          ifs.getline(buffer,BUFF_SIZE); // line with data
          tokenize(vs,buffer);
          cout << vs.size() << endl;
          if(vs.size() == 11)
          {
            atom = mol.NewAtom();
            atom->SetAtomicNum(OBElements::GetAtomicNum(vs[0].c_str()));
            x = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
            y = atof((char*)vs[7].c_str()) * BOHR_TO_ANGSTROM;
            z = atof((char*)vs[10].c_str()) * BOHR_TO_ANGSTROM;
            atom->SetVector(x,y,z);
          }
        }
      }
    }

    mol.EndModify();
    if (!pConv->IsOption("b",OBConversion::INOPTIONS)) {
      mol.ConnectTheDots();
    }

    if (!pConv->IsOption("s",OBConversion::INOPTIONS)
        && !pConv->IsOption("b",OBConversion::INOPTIONS)) {
      mol.PerceiveBondOrders();
    }

    return(true);
  }
} //namespace OpenBabel
