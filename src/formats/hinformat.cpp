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

#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class HINFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    HINFormat()
    {
      OBConversion::RegisterFormat("hin",this, "chemical/x-hin");
    }

    virtual const char* Description() //required
    {
      return "HyperChem HIN format\n"
             "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    { return "";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-hin";}; //optional

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  HINFormat theHINFormat;

  /////////////////////////////////////////////////////////////////
  bool HINFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    // Right now only read in the first molecule
    int i;
    int max, bo;
    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    ifs.getline(buffer, BUFF_SIZE);
    while (ifs.good() && (strstr(buffer,"mol") == NULL || buffer[0]==';') ) //The "mol" in comment line should be ignored.
      {
        ifs.getline(buffer, BUFF_SIZE);
        if (ifs.peek() == EOF || !ifs.good())
          return false;
      }
    ifs.getline(buffer, BUFF_SIZE);
    if (!ifs.good())
      return false; // ended early

    mol.BeginModify();
    while (ifs.good() && strstr(buffer,"endmol") == NULL)
      {
	if(buffer[0]==';'){
		 ifs.getline(buffer, BUFF_SIZE);
		 continue; //The comment Line in HIN should be ignored.
	}

        tokenize(vs,buffer); // Don't really know how long it'll be
        if (vs.size() < 11)
          {
            ifs.getline(buffer, BUFF_SIZE);
            continue;
          }

        atom = mol.NewAtom();
        atom->SetAtomicNum(OBElements::GetAtomicNum(vs[3].c_str()));
        atom->SetPartialCharge(atof(vs[6].c_str()));
        x = atof((char*)vs[7].c_str());
        y = atof((char*)vs[8].c_str());
        z = atof((char*)vs[9].c_str());
        atom->SetVector(x,y,z);

        max = 11 + 2 * atoi((char *)vs[10].c_str());
        for (i = 11; i < max; i+=2)
          {
            switch(((char*)vs[i+1].c_str())[0]) // First char in next token
              {
              case 's':
                bo = 1;
                break;
              case 'd':
                bo = 2;
                break;
              case 't':
                bo = 3;
                break;
              case 'a':
                bo = 5;
                break;
              default :
                bo = 1;
                break;
              }
            mol.AddBond(mol.NumAtoms(), atoi((char *)vs[i].c_str()), bo);
          }
        ifs.getline(buffer, BUFF_SIZE);
      }

    // clean out remaining blank lines
    // blank lines cleaning codes rewritten for avoiding peek() and tellg() bugs
    // https://github.com/openbabel/openbabel/issues/1569
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);


    mol.EndModify();

    mol.SetTitle(title);
    mol.SetPartialChargesPerceived();

    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool HINFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i, file_num = 1;
    string str,str1;
    char buffer[BUFF_SIZE];
    OBAtom *atom;
    OBBond *bond;
    vector<OBBond*>::iterator j;
    char bond_char;

    // make sure to escape titles in double quotes
    // PR#1501694
    ofs << "mol " << file_num << " \"" << mol.GetTitle() << "\"\n";

    for(i = 1;i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        snprintf(buffer, BUFF_SIZE, "atom %d - %-3s **  - %8.5f %8.5f  %8.5f  %8.5f %d ",
                i,
                OBElements::GetSymbol(atom->GetAtomicNum()),
                atom->GetPartialCharge(),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ(),
                atom->GetExplicitDegree());
        ofs << buffer;
        for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j))
          {
            switch(bond->GetBondOrder())
              {
              case 1 :
                bond_char = 's';
                break;
              case 2 :
                bond_char = 'd';
                break;
              case 3 :
                bond_char = 't';
                break;
              case 5 :
                bond_char = 'a';
                break;
              default:
                bond_char = 's';
                break;
              }
            if (bond->IsAromatic())
              bond_char = 'a';

            snprintf(buffer,BUFF_SIZE, "%d %c ", (bond->GetNbrAtom(atom))->GetIdx(), bond_char);
            ofs << buffer;
          }
        ofs << endl;
      }
    ofs << "endmol " << file_num << endl;
    return(true);
  }

} //namespace OpenBabel
