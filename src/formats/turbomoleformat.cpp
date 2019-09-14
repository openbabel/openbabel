/**********************************************************************
Copyright (C) 2004 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

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


const double AAU = 0.5291772108;  // �ngstr�m per bohr (CODATA 2002)

using namespace std;
namespace OpenBabel
{

class TurbomoleFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    TurbomoleFormat()
    {
      OBConversion::RegisterFormat("tmol",this);
      OBConversion::RegisterOptionParam("a", this, OBConversion::INOPTIONS);
    }

  virtual const char* Description() //required
  {
    return
      "TurboMole Coordinate format\n"
      "Write options e.g.-xa\n"
      " a  Output Angstroms\n\n"
      "Read Options e.g. -as\n"
      " s  Output single bonds only\n"
      " b  Disable bonding entirely\n"
      " a  Input in Angstroms\n\n" ;
  };

  virtual const char* SpecificationURL()
  {return "http://www.cosmologic.de/QuantumChemistry/main_qChemistry.html";};

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
TurbomoleFormat theTurbomoleFormat;

/////////////////////////////////////////////////////////////////
bool TurbomoleFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    double UnitConv=AAU;
    if(pConv->IsOption("a", OBConversion::INOPTIONS))
      UnitConv=1;


    char buffer[BUFF_SIZE];
    do
    {
        ifs.getline(buffer,BUFF_SIZE);
	if (ifs.peek() == EOF || !ifs.good())
	  return false;
    }
    while(strncmp(buffer,"$coord",6));

    mol.BeginModify();
    OBAtom atom;
    while(!(!ifs))
    {
        ifs.getline(buffer,BUFF_SIZE);
        if(*buffer=='$')
            break;
        if(*buffer=='#')
            continue;
        float x,y,z;
        char atomtype[8];
        if(sscanf(buffer,"%f %f %f %7s",&x,&y,&z,atomtype)!=4)
            return false;

        atom.SetVector(x*UnitConv, y*UnitConv, z*UnitConv);
        atom.SetAtomicNum(OBElements::GetAtomicNum(atomtype));
        atom.SetType(atomtype);

        if(!mol.AddAtom(atom))
            return false;
        atom.Clear();
    }
    while(!(!ifs) && strncmp(buffer,"$end",4))
        ifs.getline(buffer,BUFF_SIZE);

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

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
    return true;
}

////////////////////////////////////////////////////////////////

char *strlwr(char *s)
{
    if (s != NULL)
      {
        char *p;
        for (p = s; *p; ++p)
            *p = tolower(*p);
      }
    return s;
}


bool TurbomoleFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;
    double UnitConv=AAU;
    if(pConv->IsOption("a"))
      UnitConv=1;

    ofs << "$coord" <<endl;

    char buffer[BUFF_SIZE];
    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      char symb[8];
      strcpy(symb,OBElements::GetSymbol(atom->GetAtomicNum()));
        snprintf(buffer, BUFF_SIZE, "%20.14f  %20.14f  %20.14f      %s",
                atom->GetX()/UnitConv,
                atom->GetY()/UnitConv,
                atom->GetZ()/UnitConv,
                strlwr(symb) );
        ofs << buffer << endl;
    }
    ofs << "$end" << endl;

    return true;
}

} //namespace OpenBabel
