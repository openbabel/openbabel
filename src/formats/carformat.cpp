/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2013 by Schrodinger Inc.

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
#include <openbabel/generic.h>
#include <cstdlib>


using namespace std;
namespace OpenBabel
{

  class CARFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CARFormat()
    {
      OBConversion::RegisterFormat("car",this, "chemical/x-msi-car");
      OBConversion::RegisterFormat("arc",this, "chemical/x-msi-car");
    }

    virtual const char* Description() //required
    {
      return
        "Accelrys/MSI Biosym/Insight II CAR format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.centrcn.umontreal.ca/accelrys/life/insight2000.1/formats980/Files980TOC.doc.html";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-msi-car"; };

    virtual unsigned int Flags()
    {
      return NOTWRITABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };

  //Make an instance of the format class
  CARFormat theCARFormat;

  /////////////////////////////////////////////////////////////////
  bool CARFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    bool hasPartialCharges = false;
    char buffer[BUFF_SIZE];
    string str;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    mol.BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"end") != NULL)
          {
            if (mol.NumAtoms() > 0) // we've already read in a molecule, so exit
              break;
            // else, we hit the end of the previous molecular system
            // (in a multimolecule file)
            ifs.getline(buffer,BUFF_SIZE); // title
            ifs.getline(buffer,BUFF_SIZE); // DATE
          }

        if (strncmp(buffer, "!BIOSYM", 7) == 0)
          {
            continue;
          }

        if(strstr(buffer,"PBC") != NULL)
          {
            if(strstr(buffer,"ON") != NULL)
              {
                ifs.getline(buffer,BUFF_SIZE); // title
                ifs.getline(buffer,BUFF_SIZE); // DATE
                ifs.getline(buffer,BUFF_SIZE); // PBC a b c alpha beta gamma SG

                string str = buffer;
                // parse cell parameters
                tokenize(vs,str," \t\r\n", 7);
                if (vs.size() >= 7)
                  {
                    //parse cell values
                    double A,B,C,Alpha,Beta,Gamma;
                    A = atof((char*)vs[1].c_str());
                    B = atof((char*)vs[2].c_str());
                    C = atof((char*)vs[3].c_str());
                    Alpha = atof((char*)vs[4].c_str());
                    Beta  = atof((char*)vs[5].c_str());
                    Gamma = atof((char*)vs[6].c_str());
                    OBUnitCell *uc = new OBUnitCell;
                    uc->SetOrigin(fileformatInput);
                    uc->SetData(A, B, C, Alpha, Beta, Gamma);
                    if(vs.size() > 7) 
                      {
                        string& space_group = vs[7];

                        // Remove parentheses enclosing the space
                        // group and remove white space from front
                        // and back of string.
                        Trim(space_group);
                        if(space_group[0] == '(')
                          {
                            space_group.erase(0, 1);
                            space_group.erase(space_group.size()-1); 
                          }
                        Trim(space_group);

                        uc->SetSpaceGroup(space_group);
                      }
                    mol.SetData(uc);
                  }
              }
            else // PBC=OFF
              {
                ifs.getline(buffer,BUFF_SIZE); // title
                ifs.getline(buffer,BUFF_SIZE); // !DATE
              }
            continue;
          } // PBC

        // reading real data!
        tokenize(vs,buffer);
        if (vs.size() < 8) {
          break;
        }

        atom = mol.NewAtom();

        atom->SetAtomicNum(OBElements::GetAtomicNum(vs[7].c_str()));
        x = atof((char*)vs[1].c_str());
        y = atof((char*)vs[2].c_str());
        z = atof((char*)vs[3].c_str());
        atom->SetVector(x,y,z);

        // vs[0] contains atom label
        // vs[4] contains "type of residue containing atom"
        // vs[5] contains "residue sequence name"
        // vs[6] contains "potential type of atom"

        if (vs.size() == 9)
          {
            atom->SetPartialCharge(atof((char*)vs[8].c_str()));
            hasPartialCharges = true;
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

} //namespace OpenBabel
