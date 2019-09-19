/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
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
#include <openbabel/elements.h>
#include <openbabel/generic.h>

using namespace std;
namespace OpenBabel
{

  class CSSRFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CSSRFormat()
    {
      OBConversion::RegisterFormat("cssr",this);
    }

    virtual const char* Description() //required
    {
      return
        "CSD CSSR format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  CSSRFormat theCSSRFormat;

  ////////////////////////////////////////////////////////////////

  bool CSSRFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];

    if (!mol.HasData(OBGenericDataType::UnitCell))
      {
        snprintf(buffer, BUFF_SIZE,
                 " REFERENCE STRUCTURE = 00000   A,B,C =%8.3f%8.3f%8.3f",
                 1.0,1.0,1.0);
        ofs << buffer << endl;
        snprintf(buffer, BUFF_SIZE,
                 "   ALPHA,BETA,GAMMA =%8.3f%8.3f%8.3f    SPGR =    P1"
                 , 90.0f, 90.0f, 90.0f);
        ofs << buffer << endl;
      }
    else
      {
        OBUnitCell *uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
        snprintf(buffer, BUFF_SIZE,
                 " REFERENCE STRUCTURE = 00000   A,B,C =%8.3f%8.3f%8.3f",
                 uc->GetA(), uc->GetB(), uc->GetC());
        ofs << buffer << endl;
        snprintf(buffer, BUFF_SIZE,
                 "   ALPHA,BETA,GAMMA =%8.3f%8.3f%8.3f    SPGR =    P1",
                 uc->GetAlpha() , uc->GetBeta(), uc->GetGamma());
        ofs << buffer << endl;
      }

    snprintf(buffer, BUFF_SIZE, "%4d   1 %s\n",mol.NumAtoms(), mol.GetTitle());
    ofs << buffer << endl << endl;

    OBAtom *atom,*nbr;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;
    vector<int> vtmp(106,0);
    int bonds;

    for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        //assign_pdb_number(pdb_types,atom->GetIdx());
        vtmp[atom->GetAtomicNum()]++;
        snprintf(buffer, BUFF_SIZE, "%4d%2s%-3d  %9.5f %9.5f %9.5f ",
                 atom->GetIdx(),
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 vtmp[atom->GetAtomicNum()],
                 atom->x(),
                 atom->y(),
                 atom->z());
        ofs << buffer;
        bonds = 0;
        for (nbr = atom->BeginNbrAtom(j); nbr; nbr = atom->NextNbrAtom(j))
          {
            if (bonds > 8) break;
            snprintf(buffer, BUFF_SIZE, "%4d",nbr->GetIdx());
            ofs << buffer;
            bonds++;
          }
        for (; bonds < 8; bonds ++)
          {
            snprintf(buffer, BUFF_SIZE, "%4d",0);
            ofs << buffer;
          }
        snprintf(buffer, BUFF_SIZE, " %7.3f%4d", atom->GetPartialCharge(), 1);
        ofs << buffer << endl;
      }

    return(true);
  }

} //namespace OpenBabel
