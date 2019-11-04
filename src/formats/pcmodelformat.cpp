/**********************************************************************
Copyright (C) 2005-2006 by Geoffrey R. Hutchison

Thanks to Kevin Gilbert from Serena Software for documentation and examples!

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/data.h>
#include <openbabel/obiter.h>
#include <openbabel/bond.h>
#include <cstdlib>


using namespace std;
namespace OpenBabel
{

  class PCModelFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    PCModelFormat()
    {
      OBConversion::RegisterFormat("pcm", this);
    }

    virtual const char* Description() //required
    {
      return
        "PCModel Format\n"
        "No comments yet\n";
    };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    //  virtual unsigned int Flags()
    //  {
    //    return NOTREADABLE;
    //  };

    virtual const char* SpecificationURL()
    {return "http://www.serenasoft.com/";}; //optional

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  PCModelFormat thePCModelFormat;

  /////////////////////////////////////////////////////////////////
  bool PCModelFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();
    char buffer[BUFF_SIZE];
    string temp, temp2;

    OBAtom *atom;
    vector<string> vs;
    double x, y, z;
    unsigned int token;
    int bondNbr, bondOrder;
    bool parsingBonds, readingMol = false;
    bool hasPartialCharges = false;

    ttab.SetFromType("PCM");

    mol.BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if(strncmp(buffer,"{PCM", 4) == 0)
          {
            temp = buffer;
            temp = temp.substr(4, temp.length());
            mol.SetTitle(temp);
            readingMol = true;
          }
        else if (readingMol && strncmp(buffer,"}", 1) == 0)
          {
            readingMol = false;
            break;
          }
        else if (readingMol && strncmp(buffer,"AT ",3) == 0)
          {
            tokenize(vs,buffer, "\n\r\t ,:");
            if (vs.size() < 3)
              return false;

            atom = mol.NewAtom();
            temp = vs[2].c_str();
            ttab.SetToType("INT");
            ttab.Translate(temp2, temp);
            atom->SetType(temp2);

            ttab.SetToType("ATN");
            ttab.Translate(temp2, temp);
            atom->SetAtomicNum(atoi(temp2.c_str()));
            x = atof(vs[3].c_str());
            y = atof(vs[4].c_str());
            z = atof(vs[5].c_str());
            atom->SetVector(x,y,z); //set coordinates

            token = 6;
            parsingBonds = false;
            while(token < vs.size())
              {
                if (vs[token] == "B")
                  parsingBonds = true;
                else if (vs[token][0] == 'C')
                  {
                    parsingBonds = false;
                    hasPartialCharges = true;
                    if (vs[token].size() > 1)
                      temp = vs[token].substr(1,vs[token].size());
                    else
                      {
                        token++;
                        temp = vs[token];
                      }
                    atom->SetPartialCharge(atof(temp.c_str()));
                  }

                else if (parsingBonds && token < vs.size() - 1 &&
                         isdigit(vs[token][0]))
                  {
                    bondNbr = atoi(vs[token++].c_str()); // advance to bond order
                    bondOrder = atoi(vs[token].c_str());
                    if (bondOrder == 9)
                      bondOrder = 1;
                    mol.AddBond(atom->GetIdx(), bondNbr, bondOrder, 0);
                  }
                else
                  parsingBonds = false; // any other token

                token++;
              } // end atom fields
          } // end AT line
      } // end reading

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

  bool PCModelFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;
    OBAtom *nbr;
    vector<OBBond*>::iterator j;
    string type, temp;
    int nbrIdx, atomIdx;

    temp = mol.GetTitle();
    ofs << "{PCM " << temp.substr(0,60) << endl;
    ofs << "NA " << mol.NumAtoms() << endl;
    ofs << "ATOMTYPES 1" << endl; // MMX atom types

    ttab.SetFromType("INT");
    ttab.SetToType("PCM");

    string str,str1;
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        ttab.Translate(type,atom->GetType());
        atomIdx = atom->GetIdx();

        ofs << "AT " << atom->GetIdx() << "," << type << ":";
        ofs << atom->GetX() << "," << atom->GetY() << "," << atom->GetZ();

        if (atom->GetExplicitDegree() > 0)
          {
            ofs << " B";
            for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
              {
                nbrIdx = nbr->GetIdx();
                ofs << " " << nbrIdx << ","
                    << (mol.GetBond(nbrIdx, atomIdx))->GetBondOrder();
              }
          }

        ofs << " C " << atom->GetPartialCharge();

        ofs << endl;
      }

    ofs << "}" << endl;

    return(true);
  }

} //namespace OpenBabel
