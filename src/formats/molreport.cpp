/**********************************************************************
molreport.cpp - Report information about the molecule:
     atom list, bond list, etc.

Copyright (C) 2006 by Geoffrey R. Hutchison

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

using namespace std;
namespace OpenBabel
{

  class MolReportFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MolReportFormat()
    {
      OBConversion::RegisterFormat("molreport",this);
    }

    virtual const char* Description() //required
    {
      return
        "Open Babel molecule report\n"
        "Generates a summary of the atoms and bonds in a molecule\n"
        "Example output::\n\n"
        " TITLE: Ethanol.mopout\n"
" FORMULA: C2H6O\n"
" MASS: 46.0684\n"
" ATOM:         1   C TYPE: C3     HYB:  3 CHARGE:  -0.2151\n"
" ATOM:         2   C TYPE: C3     HYB:  3 CHARGE:  -0.0192\n"
" ATOM:         3   O TYPE: O3     HYB:  3 CHARGE:  -0.3295\n"
" ATOM:         4   H TYPE: HC     HYB:  0 CHARGE:   0.0771\n"
" ATOM:         5   H TYPE: HC     HYB:  0 CHARGE:   0.0873\n"
" ATOM:         6   H TYPE: HC     HYB:  0 CHARGE:   0.0874\n"
" ATOM:         7   H TYPE: HC     HYB:  0 CHARGE:   0.0577\n"
" ATOM:         8   H TYPE: HC     HYB:  0 CHARGE:   0.0577\n"
" ATOM:         9   H TYPE: HC     HYB:  0 CHARGE:   0.1966\n"
" BOND:         0 START:         8 END:         2 ORDER:   1\n"
" BOND:         1 START:         6 END:         1 ORDER:   1\n"
" BOND:         2 START:         1 END:         2 ORDER:   1\n"
" BOND:         3 START:         1 END:         4 ORDER:   1\n"
" BOND:         4 START:         1 END:         5 ORDER:   1\n"
" BOND:         5 START:         2 END:         3 ORDER:   1\n"
" BOND:         6 START:         2 END:         7 ORDER:   1\n"
" BOND:         7 START:         3 END:         9 ORDER:   1\n\n"

".. seealso::\n\n"

" :ref:`Open_Babel_report_format`\n\n"

;
    };

    virtual const char* SpecificationURL()
    {return "http://openbabel.org/wiki/MolReport";}; //optional

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
  MolReportFormat theMolReportFormat;

  ////////////////////////////////////////////////////////////////

  bool MolReportFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    ofs << "TITLE: " << mol.GetTitle() << "\n";
    ofs << "FORMULA: " << mol.GetFormula() << "\n";
    ofs << "MASS: ";
    snprintf(buffer, BUFF_SIZE, "%5.4f\n", mol.GetMolWt());
    ofs << buffer;

    if (mol.GetTotalCharge() != 0)
      {
        ofs << "TOTAL CHARGE: ";
        snprintf(buffer, BUFF_SIZE, "%d", mol.GetTotalCharge());
        ofs << buffer << "\n";
      }
    if (mol.GetTotalSpinMultiplicity() != 1)
      {
        ofs << "TOTAL SPIN: ";
        snprintf(buffer, BUFF_SIZE, "%d", mol.GetTotalSpinMultiplicity());
        ofs << buffer << "\n";
      }

    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "ATOM: %9d %3s TYPE: %-6s HYB: %2d CHARGE: %8.4f",
                 atom->GetIdx(),
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 atom->GetType(),
                 atom->GetHyb(),
                 atom->GetPartialCharge());
        ofs << buffer << "\n";
      }

    FOR_BONDS_OF_MOL(bond, mol)
      {
        snprintf(buffer, BUFF_SIZE, "BOND: %9d START: %9d END: %9d ORDER: %3d",
                 bond->GetIdx(),
                 bond->GetBeginAtomIdx(),
                 bond->GetEndAtomIdx(),
                 bond->GetBondOrder());
        ofs << buffer << "\n";
      }

    return(true);
  }

} //namespace OpenBabel
