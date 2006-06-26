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
#include "babelconfig.h"

#include "obmolecformat.h"

using namespace std;
namespace OpenBabel
{

#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989

  class GAMESSOutputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    GAMESSOutputFormat()
    {
      OBConversion::RegisterFormat("gam",this);
      OBConversion::RegisterFormat("gamout",this);
    }

    virtual const char* Description() //required
    {
      return
        "GAMESS Output\n \
       Read Options e.g. -as\n\
        s  Output single bonds only\n\
        b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.msg.ameslab.gov/GAMESS/doc.menu.html";}; //optional

    virtual const char* GetMIMEType() 
    { return "chemical/x-gamess-output"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

  };
  //***

  //Make an instance of the format class
  GAMESSOutputFormat theGAMESSOutputFormat;


  class GAMESSInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    GAMESSInputFormat()
    {
      OBConversion::RegisterFormat("inp",this);
      OBConversion::RegisterFormat("gamin",this);
    }

    virtual const char* Description() //required
    {
      return
        "GAMESS Input\n \
            No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.msg.ameslab.gov/GAMESS/doc.menu.html";}; //optional

    virtual const char* GetMIMEType() 
    { return "chemical/x-gamess-input"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return WRITEONEONLY | NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  GAMESSInputFormat theGAMESSInputFormat;

  /////////////////////////////////////////////////////////////////
  bool GAMESSOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
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
    bool hasPartialCharges = false;

    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"ATOMIC                      COORDINATES (BOHR)") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 5)
              {
                atom = mol.NewAtom();
                atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
                x = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
                y = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
                z = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
                atom->SetVector(x,y,z);
                vs[1].erase(vs[1].size() - 2, 2);

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"COORDINATES OF ALL ATOMS ARE (ANGS)") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);	// ---------------
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 5)
              {
                atom = mol.NewAtom();
                atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
                x = atof((char*)vs[2].c_str());
                y = atof((char*)vs[3].c_str());
                z = atof((char*)vs[4].c_str());
                atom->SetVector(x,y,z);
                vs[1].erase(vs[1].size() - 2, 2);

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"MOPAC CHARGES") != NULL)
          {
            hasPartialCharges = true;
            ifs.getline(buffer,BUFF_SIZE);	// ---------------
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 4)
              {
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                atom->SetPartialCharge(atof(vs[2].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
      }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    if (hasPartialCharges)
      mol.SetPartialChargesPerceived();
    mol.SetTitle(title);
    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool GAMESSInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];

    ofs << " $CONTRL COORD=CART UNITS=ANGS $END" << endl;
    ofs << " $DATA" << endl;
    ofs << mol.GetTitle() << endl;
    ofs << "Put symmetry info here" << endl << endl;

    OBAtom *atom;
    for(i = 1;i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        snprintf(buffer, BUFF_SIZE, "%-3s %4d.0    %8.5f  %8.5f  %8.5f ",
                etab.GetSymbol(atom->GetAtomicNum()),
                atom->GetAtomicNum(),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs << buffer << endl;
      }

    ofs << " $END" << endl << endl << endl;
    return(true);
  }

} //namespace OpenBabel
