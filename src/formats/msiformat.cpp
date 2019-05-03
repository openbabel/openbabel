/**********************************************************************
Copyright (C) 2008 by Geoffrey Hutchison

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


#include <sstream>

using namespace std;
namespace OpenBabel
{

  class MSIFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MSIFormat()
    {
      OBConversion::RegisterFormat("msi", this, "chemical/x-msi-msi");
    }

    virtual const char* Description() //required
    {
      return
        "Accelrys/MSI Cerius II MSI format\n";
    };

    virtual const char* SpecificationURL()
    {return "http://openbabel.org/wiki/MSI_format";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-msi-msi"; };

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
  MSIFormat theMSIFormat;

  /////////////////////////////////////////////////////////////////
  bool MSIFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();
    char buffer[BUFF_SIZE];

    stringstream errorMsg;

    if (!ifs)
      return false; // we're attempting to read past the end of the file

    if (!ifs.getline(buffer,BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading an MSI file: Cannot read the first line.", obWarning);
        return(false);
      }

    if (!EQn(buffer, "# MSI CERIUS2 DataModel File", 28))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problems reading an MSI file: The first line must contain the MSI header.", obWarning);
        return(false);
      }

    // "records" start with
    // (1 Model
    // ....
    //   and end with
    // ....
    // )
    unsigned int openParens = 0; // the count of "open parentheses" tags
    unsigned int startBondAtom, endBondAtom, bondOrder;
    bool atomRecord = false;
    bool bondRecord = false;
    OBAtom *atom;
//    OBBond *bond;
    vector<string> vs;
    const SpaceGroup *sg;
    bool setSpaceGroup = false;
    double x,y,z;
    vector3 translationVectors[3];
    int numTranslationVectors = 0;

    mol.BeginModify();
    while (ifs.getline(buffer,BUFF_SIZE))
      {
        // model record
        if (strstr(buffer, "Model") != NULL) {
          openParens++;
          continue;
        }

        // atom record
        if (!bondRecord && strstr(buffer, "Atom") != NULL) {
          atomRecord = true;
          openParens++;
          continue;
        }

        if (strstr(buffer, "Bond") != NULL) {
          bondRecord = true;
          startBondAtom = endBondAtom = 0;
          bondOrder = 1;
          openParens++;
          continue;
        }

        /* (A I PeriodicType 100)
           (A D A3 (6.2380000000000004 0 0))
           (A D B3 (0 6.9909999999999997 0))
           (A D C3 (0 0 6.9960000000000004))
           (A C SpaceGroup "63 5")
        */
        if (strstr(buffer, "PeriodicType") != NULL) {
          ifs.getline(buffer,BUFF_SIZE); // next line should be translation vector
          tokenize(vs,buffer);
            while (vs.size() == 6) {
              x = atof((char*)vs[3].erase(0,1).c_str());
              y = atof((char*)vs[4].c_str());
              z = atof((char*)vs[5].c_str());

              translationVectors[numTranslationVectors++].Set(x, y, z);
              if (!ifs.getline(buffer,BUFF_SIZE))
                break;
              tokenize(vs,buffer);
            }
        }

        if (strstr(buffer, "SpaceGroup") != NULL) {
          tokenize(vs, buffer);
          if (vs.size() != 5)
            continue; // invalid space group
          setSpaceGroup = true;
          sg = SpaceGroup::GetSpaceGroup(vs[4]); // remove the initial " character
        }

        // atom information
        if (atomRecord) {
          if (strstr(buffer, "ACL") != NULL) {
            tokenize(vs, buffer);
            // size should be 5 -- need a test here
            if (vs.size() != 5) return false; // timvdm 18/06/2008
            vs[3].erase(0,1); // "6 => remove the first " character
            unsigned int atomicNum = atoi(vs[3].c_str());
            if (atomicNum == 0)
              atomicNum = 1; // hydrogen ?

            // valid element, so create the atom
            atom = mol.NewAtom();
            atom->SetAtomicNum(atomicNum);
            continue;
          }
          else if (strstr(buffer, "XYZ") != NULL) {
            tokenize(vs, buffer);
            // size should be 6 -- need a test here
            if (vs.size() != 6) return false; // timvdm 18/06/2008
            vs[3].erase(0,1); // remove ( character
            vs[5].erase(vs[5].length()-2, 2); // remove trailing )) characters
            atom->SetVector(atof(vs[3].c_str()),
                            atof(vs[4].c_str()),
                            atof(vs[5].c_str()));
            continue;
          }
        } // end of atom records

        // bond information
        if (bondRecord) {
          if (strstr(buffer, "Atom1") != NULL) {
            tokenize(vs, buffer);
            if (vs.size() < 4) return false; // timvdm 18/06/2008
            vs[3].erase(vs[3].length()-1,1);
            startBondAtom = atoi(vs[3].c_str());
            continue;
          }
          else if (strstr(buffer, "Atom2") != NULL) {
            tokenize(vs, buffer);
            if (vs.size() < 4) return false; // timvdm 18/06/2008
            vs[3].erase(vs[3].length()-1,1);
            endBondAtom = atoi(vs[3].c_str());
            continue;
          }
          else if (strstr(buffer, "Type") != NULL) {
            tokenize(vs, buffer);
            if (vs.size() < 4) return false; // timvdm 18/06/2008
            vs[3].erase(vs[3].length()-1,1);
            bondOrder = atoi(vs[3].c_str());
            if (bondOrder == 4) // triple bond?
              bondOrder = 3;
            else if (bondOrder == 8) // aromatic?
              bondOrder = 5;
            else if (bondOrder != 2) // 1 OK, 2 OK, others unknown
              bondOrder = 1;
            continue;
          }
        }

        // ending a "tag" -- a lone ")" on a line
        if (strstr(buffer,")") != NULL && strstr(buffer, "(") == NULL) {
          openParens--;
          if (atomRecord) {
            atomRecord = false;
          }
          if (bondRecord) {
            // Bond records appear to be questionable
            mol.AddBond(startBondAtom - 1, endBondAtom - 1, bondOrder);
            bondRecord = false;
          }

          if (openParens == 0) {
            ifs.getline(buffer, BUFF_SIZE);
            break; // closed this molecule
          }
        }
      }

    mol.EndModify();

    // clean out any remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    /*
    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();
    */

    if (numTranslationVectors > 0) {
      OBUnitCell* uc = new OBUnitCell;
      uc->SetData(translationVectors[0], translationVectors[1], translationVectors[2]);
      uc->SetOrigin(fileformatInput);
      if (setSpaceGroup) {
        uc->SetSpaceGroup(sg);
      }
      mol.SetData(uc);
    }

    return(true);
  }

} //namespace OpenBabel
