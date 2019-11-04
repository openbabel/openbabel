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


using namespace std;
namespace OpenBabel
{

#define BOHR_TO_ANGSTROM 0.529177249

  class MPQCFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MPQCFormat()
    {
      OBConversion::RegisterFormat("mpqc",this);
    }

    virtual const char* Description() //required
    {
      return
        "MPQC output format\n"
        "Read Options e.g. -as\n"
        " s  Output single bonds only\n"
        " b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.mpqc.org/mpqc-html/mpqcinp.html";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTWRITABLE | READONEONLY;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };

  //Make an instance of the format class
  MPQCFormat theMPQCFormat;



  class MPQCInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MPQCInputFormat()
    {
      OBConversion::RegisterFormat("mpqcin",this);
    }

    virtual const char* Description() //required
    {
      return
        "MPQC simplified input format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.mpqc.org/mpqc-html/mpqcinp.html";}; //optional

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
  MPQCInputFormat theMPQCInputFormat;

  /////////////////////////////////////////////////////////////////
  bool MPQCFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    bool bohr = true;

    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"<Molecule>:") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            while	(strstr(buffer,"geometry") == NULL)
              {
                if (strstr(buffer,"angstrom") != NULL)
                  bohr = false;
                if (!ifs.getline(buffer,BUFF_SIZE))
                  return(false);
              }
            ifs.getline(buffer,BUFF_SIZE); // Now we're on the atoms
            tokenize(vs,buffer);
            while (vs.size() == 6)
              {
                if (bohr)
                  {
                    x = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
                    y = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
                    z = atof((char*)vs[5].c_str()) * BOHR_TO_ANGSTROM;
                  }
                else
                  {
                    x = atof((char*)vs[3].c_str());
                    y = atof((char*)vs[4].c_str());
                    z = atof((char*)vs[5].c_str());
                  }
                atom = mol.NewAtom();
                atom->SetVector(x,y,z);
                atom->SetAtomicNum(OBElements::GetAtomicNum(vs[1].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
      }

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();

    mol.SetTitle(title);

    return(true);
  }

  bool MPQCInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    //   unsigned int i;
    char buffer[BUFF_SIZE];

    ofs << "% " << mol.GetTitle() << "\n";
    ofs << "\n"; // keywords/direction lines here
    ofs << "molecule:\n";

    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "%4s  %8.5f  %8.5f  %8.5f \n",
                OBElements::GetSymbol(atom->GetAtomicNum()),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs << buffer;
      }

    ofs << "\n\n\n";
    return(true);
  }

} //namespace OpenBabel
