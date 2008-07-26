/**********************************************************************
Copyright (C) 2008 Geoffrey R. Hutchison
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

#include <vector>
#include <map>

#include <sstream>

using namespace std;
namespace OpenBabel
{

  class PQRFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    PQRFormat()
    {
      OBConversion::RegisterFormat("pqr",this, "chemical/x-pqr");
    }

    virtual const char* Description() //required
    {
      return
        "PQR format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "";};

    virtual const char* GetMIMEType() 
    { return "chemical/x-pqr"; };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
  	virtual int SkipObjects(int n, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    //    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };
  //***

  //Make an instance of the format class
  PQRFormat thePQRFormat;

  /////////////////////////////////////////////////////////////////

  // In atomrecord.cpp (shared with PDB format)
  bool ParseAtomRecord(char *, OBMol &, int);
  double ParseAtomCharge(char *, OBMol &);

  /////////////////////////////////////////////////////////////////
 	int PQRFormat::SkipObjects(int n, OBConversion* pConv)
  {
    if (n == 0)
      ++ n;
    istream &ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    while (n && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6))
          -- n;
      }
      
    return ifs.good() ? 1 : -1;       
  }
  /////////////////////////////////////////////////////////////////
  bool PQRFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int chainNum = 1;
    char buffer[BUFF_SIZE];
    OBBitVec bs;
    vector<double> charges;
    string line, key, value;

    mol.SetTitle(title);
    mol.SetChainsPerceived(); // It's a PDB-like file, we read all chain/res info.

    mol.BeginModify();
    while (ifs.good() && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6))
          break;
        if (EQn(buffer,"END",3)) {
          // eat anything until the next ENDMDL
          while (ifs.getline(buffer,BUFF_SIZE) && !EQn(buffer,"ENDMDL",6));
          break;
        }
        if (EQn(buffer,"TER",3)) {
          chainNum++;
          continue;
        }
        if (EQn(buffer,"ATOM",4) || EQn(buffer,"HETATM",6))
          {
            ParseAtomRecord(buffer,mol,chainNum);
            if (EQn(buffer,"ATOM",4))
              bs.SetBitOn(mol.NumAtoms());

            // Read in the partial charge too
            charges.push_back( ParseAtomCharge(buffer, mol) );
            continue;
          }
        }

    if (!mol.NumAtoms()) { // skip the rest of this processing
      mol.EndModify();
      return(false);
    }

    // Use residue definitions to assign bond orders
    resdat.AssignBonds(mol,bs);

    mol.EndModify();

    /*Now assign hetatm bonds based on distance*/
    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();

    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    FOR_ATOMS_OF_MOL(a, mol) {
      // WARNING: Atom index issue here
      a->SetPartialCharge(charges[a->GetIdx() - 1]);

      cerr << " charge : " << charges[a->GetIdx() - 1] << endl;
    }
    mol.SetPartialChargesPerceived();

    // clean out remaining blank lines
    while(ifs.peek() != EOF && ifs.good() && 
          (ifs.peek() == '\n' || ifs.peek() == '\r'))
      ifs.getline(buffer,BUFF_SIZE);

    return(true);
  }

  double ParseAtomCharge(char *buffer, OBMol &mol)
  // In PQR format, either:
  // Field name, atom number, atom name, residue name, residue number
  //    x y z charge radius
  // OR
  // Field, atom number, atom name, chain id, residue number, X, Y, Z, chg, rad
  {
    vector<string> vs;
    tokenize(vs,buffer);

    OBAtom *atom = mol.GetAtom(mol.NumAtoms());

    if (vs.size() == 10)
      return atof(vs[8].c_str());
    else if (vs.size() == 11)
      return atof(vs[9].c_str());

    return 0.0;
  }


} //namespace OpenBabel
