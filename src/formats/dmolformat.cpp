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
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989


  class DMolFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    DMolFormat()
    {
      OBConversion::RegisterFormat("dmol",this);
      OBConversion::RegisterFormat("outmol",this, "chemical/x-dmol");
    }

    virtual const char* Description() //required
    {
      return
        "DMol3 coordinates format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "" ;}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | WRITEONEONLY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  DMolFormat theDMolFormat;

  /////////////////////////////////////////////////////////////////
  bool DMolFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str;
    double x,y,z;
    OBAtom *atom;
    vector3 v1,v2,v3;
    vector<string> vs;

    ifs.getline(buffer,BUFF_SIZE);
    while (strstr(buffer,"$coordinates") == NULL &&
           strstr(buffer,"$cell vectors") == NULL)
      {
        if (ifs.peek() == EOF || !ifs.good())
          return false;
        ifs.getline(buffer,BUFF_SIZE);
      }

    if (strstr(buffer,"$cell vectors") != NULL)
      {
        ifs.getline(buffer,BUFF_SIZE);
        tokenize(vs,buffer); // we really need to check that it's 3 entries only
	if (vs.size() < 3) return false; // timvdm 18/06/2008
        x = atof((char*)vs[0].c_str()) * BOHR_TO_ANGSTROM;
        y = atof((char*)vs[1].c_str()) * BOHR_TO_ANGSTROM;
        z = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
        v1.Set(x,y,z);
        ifs.getline(buffer,BUFF_SIZE);
        tokenize(vs,buffer);
	if (vs.size() < 3) return false; // timvdm 18/06/2008
        x = atof((char*)vs[0].c_str()) * BOHR_TO_ANGSTROM;
        y = atof((char*)vs[1].c_str()) * BOHR_TO_ANGSTROM;
        z = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
        v2.Set(x,y,z);
        ifs.getline(buffer,BUFF_SIZE);
        tokenize(vs,buffer);
	if (vs.size() < 3) return false; // timvdm 18/06/2008
        x = atof((char*)vs[0].c_str()) * BOHR_TO_ANGSTROM;
        y = atof((char*)vs[1].c_str()) * BOHR_TO_ANGSTROM;
        z = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
        v3.Set(x,y,z);

        OBUnitCell *uc = new OBUnitCell;
        uc->SetOrigin(fileformatInput);
        uc->SetData(v1,v2,v3);
        mol.SetData(uc);

        ifs.getline(buffer,BUFF_SIZE); // next line
      }

    mol.BeginModify();

    while (strstr(buffer,"$end") == NULL)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          break;
        tokenize(vs,buffer);
        if (vs.size() != 4)
          break;
        atom = mol.NewAtom();
        //set atomic number
        atom->SetAtomicNum(OBElements::GetAtomicNum(vs[0].c_str()));
        x = atof((char*)vs[1].c_str()) * BOHR_TO_ANGSTROM;
        y = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
        z = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
        atom->SetVector(x,y,z); //set coordinates
      }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    // clean out any remaining blank lines
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
    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool DMolFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];

    if (mol.HasData(OBGenericDataType::UnitCell))
      {
        OBUnitCell *uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
        vector<vector3> v = uc->GetCellVectors();
        vector3 v1;

        ofs << "$cell vectors" << endl;
        v1 = v[0] * ANGSTROM_TO_BOHR;
        snprintf(buffer, BUFF_SIZE,
                 "%-3s% 27.14f% 20.14f% 20.14f","", v1.x(), v1.y(), v1.z());
        ofs << buffer << endl;
        v1 = v[1] * ANGSTROM_TO_BOHR;
        snprintf(buffer, BUFF_SIZE,
                 "%-3s% 27.14f% 20.14f% 20.14f","", v1.x(), v1.y(), v1.z());
        ofs << buffer << endl;
        v1 = v[2] * ANGSTROM_TO_BOHR;
        snprintf(buffer, BUFF_SIZE,
                 "%-3s% 27.14f% 20.14f% 20.14f","", v1.x(), v1.y(), v1.z());
        ofs << buffer << endl;
      }

    ofs << "$coordinates" << endl;

    OBAtom *atom;
    for(i = 1;i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        snprintf(buffer, BUFF_SIZE, "%-3s% 27.14f% 20.14f% 20.14f",
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 atom->GetX() * ANGSTROM_TO_BOHR,
                 atom->GetY() * ANGSTROM_TO_BOHR,
                 atom->GetZ() * ANGSTROM_TO_BOHR);
        ofs << buffer << endl;
      }

    ofs << "$end" << endl;

    return(true);
  }

} //namespace OpenBabel
