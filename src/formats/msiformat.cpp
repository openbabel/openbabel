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
        "MSI text format\n";
    };

    virtual const char* SpecificationURL()
    {return "http://openbabel.sourceforge.net/wiki/MSI_format";}; //optional

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
    bool atomRecord;
    bool bondRecord;
    OBAtom *atom;
    OBBond *bond;
    vector<string> vs;

    cerr << " starting read " << endl;

    mol.BeginModify();
    while (ifs.getline(buffer,BUFF_SIZE))
      {
        cerr << "buffer: " << buffer << endl;

        // model record
        if (strstr(buffer, " Model") != NULL) {
          openParens++;
          continue;
        }

        // atom record
        if (strstr(buffer, " Atom") != NULL) {
          atomRecord = true;
          openParens++;
          cerr << " starting atom " << endl;
          continue;
        }

        if (strstr(buffer, " Bond") != NULL) {
          bondRecord = true;
          startBondAtom = endBondAtom = 0;
          bondOrder = 1;
          openParens++;
          cerr << " stating bond " << endl;
          continue;
        }
        
        // atom information
        if (atomRecord) {
          if (strstr(buffer, "ACL") != NULL) {
            atom = mol.NewAtom();
            tokenize(vs, buffer);
            // size should be 5 -- need a test here
            vs[3].erase(0,1); // "6 => remove the first " character
            int atomicNum = atoi(vs[3].c_str());
            if (atomicNum == 0)
              atomicNum = 1; // hydrogen atom
            atom->SetAtomicNum(atoi(vs[3].c_str()));
            continue;
          }
          else if (strstr(buffer, "XYZ") != NULL) {
            tokenize(vs, buffer);
            // size should be 6 -- need a test here
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
            vs[3].erase(vs[3].length()-1,1);
            startBondAtom = atoi(vs[3].c_str());
            continue;
          }
          else if (strstr(buffer, "Atom2") != NULL) {
            tokenize(vs, buffer);
            vs[3].erase(vs[3].length()-1,1);
            endBondAtom = atoi(vs[3].c_str());
            continue;
          }
          else if (strstr(buffer, "Type") != NULL) {
            tokenize(vs, buffer);
            vs[3].erase(vs[3].length()-1,1);
            bondOrder = atoi(vs[3].c_str());
            continue;
          }
        }

        // ending a "tag" -- a lone ")" on a line
        if (strstr(buffer,")") != NULL && strstr(buffer, ")") == NULL) {
          cerr << " closing tag " << openParens << endl;
          openParens--;
          if (atomRecord) {
            atomRecord = false;
          }
          if (bondRecord) {
            mol.AddBond(startBondAtom, endBondAtom, bondOrder);
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
    while(ifs.peek() != EOF && ifs.good() && 
          (ifs.peek() == '\n' || ifs.peek() == '\r'))
      ifs.getline(buffer,BUFF_SIZE);

    return(true);
  }

} //namespace OpenBabel
