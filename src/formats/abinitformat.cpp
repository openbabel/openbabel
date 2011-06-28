/**********************************************************************
Copyright (C) 2011 by Geoffrey Hutchison

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

using namespace std;
namespace OpenBabel
{

#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989


  class ABINITFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ABINITFormat()
    {
      OBConversion::RegisterFormat("abinit",this);
    }

    virtual const char* Description() //required
    {
      return
        "ABINIT Output Format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://abinit.org/" ;}; //optional

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
  ABINITFormat theABINITFormat;

  /////////////////////////////////////////////////////////////////
  bool ABINITFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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
        atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
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
    while(ifs.peek() != EOF && ifs.good() &&
          (ifs.peek() == '\n' || ifs.peek() == '\r'))
      ifs.getline(buffer,BUFF_SIZE);

    mol.EndModify();
    mol.SetTitle(title);
    return(true);
  }

} //namespace OpenBabel
