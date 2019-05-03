/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
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
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <ctype.h>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class JaguarOutputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    JaguarOutputFormat()
    {
      OBConversion::RegisterFormat("jout",this);
    }

    virtual const char* Description() //required
    {
      return
        "Jaguar output format\n"
        "Read Options e.g. -as\n"
        " s  Output single bonds only\n"
        " b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.schrodinger.com/"; }; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags(){return READONEONLY | NOTWRITABLE;};

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  JaguarOutputFormat theJaguarOutputFormat;



  class JaguarInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    JaguarInputFormat()
    {
      OBConversion::RegisterFormat("jin",this);
    }

    virtual const char* Description() //required
    {
      return
        "Jaguar input format\n"
        "Read Options e.g. -as\n"
        " s  Output single bonds only\n"
        " b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.schrodinger.com/"; }; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags() {return WRITEONEONLY;};

    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  JaguarInputFormat theJaguarInputFormat;



  /////////////////////////////////////////////////////////////////
  bool JaguarOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL) return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    unsigned int i;
    OBAtom *atom;
    vector<string> vs;

    mol.BeginModify();
    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if (strstr(buffer,"Input geometry:") != NULL
            || strstr(buffer,"symmetrized geometry:") != NULL
            || strstr(buffer,"new geometry:") != NULL
            || strstr(buffer,"final geometry:") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);  //     angstroms
            ifs.getline(buffer,BUFF_SIZE);  // column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 4)
              {
                atom = mol.NewAtom();
                str = vs[0]; // Separate out the Symbol# into just Symbol ...
                for (i = 0;i < str.size();i++)
                  if (isdigit(str[i]))
                    str[i] = '\0';

                atom->SetAtomicNum(OBElements::GetAtomicNum(str.c_str()));
                x = atof((char*)vs[1].c_str());
                y = atof((char*)vs[2].c_str());
                z = atof((char*)vs[3].c_str());
                atom->SetVector(x,y,z);

                if (!ifs.getline(buffer,BUFF_SIZE)) break;
                tokenize(vs,buffer);
              }
          }
        if (strstr(buffer, "Atomic charges from electrostatic potential") != NULL)
          {
            mol.SetAutomaticPartialCharge(false);
            unsigned int chgcount=0;
            while (chgcount<mol.NumAtoms())
              {
                ifs.getline(buffer,BUFF_SIZE);  // blank line
                ifs.getline(buffer,BUFF_SIZE);  // header line
                ifs.getline(buffer,BUFF_SIZE);  // data line
                tokenize(vs,buffer);
                for (vector<string>::size_type icount=1;icount<vs.size();++icount)
                  {
                    chgcount=chgcount+1;
                    mol.GetAtom(chgcount)->SetPartialCharge(atof((char*)vs[icount].c_str()));
                  }
              }
          }
        else if(strstr(buffer,"Dipole Moments (Debye)") != NULL)
          {
            ifs.getline(buffer,BUFF_SIZE); // actual components   X ###  Y #### Z ###
            tokenize(vs,buffer);
            if (vs.size() >= 8)
              {
                OBVectorData *dipoleMoment = new OBVectorData;
                dipoleMoment->SetAttribute("Dipole Moment");
                double x, y, z;
                x = atof(vs[1].c_str());
                y = atof(vs[3].c_str());
                z = atof(vs[5].c_str());
                dipoleMoment->SetData(x, y, z);
                dipoleMoment->SetOrigin(fileformatInput);
                mol.SetData(dipoleMoment);
              }
            if (!ifs.getline(buffer,BUFF_SIZE)) break;
          }
      }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS)
        && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    mol.SetTitle(title);
    return(true);
  }

  ////////////////////////////////////////////////////////////////
  bool JaguarInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL) return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    OBAtom *atom;

    ofs << mol.GetTitle() << endl << endl;
    ofs << "&gen" << endl;
    ofs << "&" << endl;
    ofs << "&zmat" << endl;

    for (i = 1;i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        snprintf(buffer, BUFF_SIZE, "  %s%d   %12.7f  %12.7f  %12.7f",
                 OBElements::GetSymbol(atom->GetAtomicNum()), i,
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ());
        ofs << buffer << endl;
      }

    ofs << "&" << endl;
    return(true);
  }

  /////////////////////////////////////////////////////////////////
  bool JaguarInputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL) return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    unsigned int i;
    OBAtom *atom;
    vector<string> vs;

    mol.BeginModify();
    mol.Clear();
    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if (strstr(buffer,"&zmat") != NULL)
          {
            ifs.getline(buffer,BUFF_SIZE);  //   AtomLabel, X, Y, Z
            tokenize(vs,buffer);
            while (vs.size() == 4)
              {
                atom = mol.NewAtom();
                str = vs[0]; // Separate out the Symbol# into just Symbol ...
                for (i = 0;i < str.size();i++)
                  if (isdigit(str[i]))
                    str[i] = '\0';

                atom->SetAtomicNum(OBElements::GetAtomicNum(str.c_str()));
                x = atof((char*)vs[1].c_str());
                y = atof((char*)vs[2].c_str());
                z = atof((char*)vs[3].c_str());
                atom->SetVector(x,y,z);

                if (!ifs.getline(buffer,BUFF_SIZE)) break;
                tokenize(vs,buffer);
              }
          }
      }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS)
        && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    mol.SetTitle(title);
    return(true);
  }

} //namespace OpenBabel
