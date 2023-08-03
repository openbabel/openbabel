/**********************************************************************
neutralize.cpp - The option --neutralize neutralizes charged atoms

Copyright (C) 2020 by Noel O'Boyle

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
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>

namespace OpenBabel
{

class OpNeutralize : public OBOp
{
public:
  OpNeutralize(const char* ID) : OBOp(ID, false){};
  const char* Description() {
    return "Neutralize +1 and -1 charges\n\n"

      "Neutralize uses a simple procedure to generate the neutral form of a\n"
      "molecule. It does not attempt to balance charges but simply to convert\n"
      "all atoms with a +1 or -1 charge to neutral by addition or subtraction\n"
      "of H+.\n\n"

      "To a first approximation the procedure is simply to identify all atoms with\n"
      "either a +1 or -1 charge, set their charges to zero and adjust their\n"
      "hydrogen counts by -1 or +1 (i.e. we are adding/removing H+). The first\n"
      "minor issue is that +1 charged atoms must have a hydrogen, or otherwise\n"
      "we can't remove H+. The second issue is that we must avoid altering\n"
      "charge-separated representations of multiple bonds, such as nitro which\n"
      "is often represented as [N+](=O)[O-]. It does this by checking whether a\n"
      "neighbor atom has a negative charge (for the +1 atoms) or a positive\n"
      "charge (for the -1 atoms).\n\n"

      "If specified, the optional argument 'changed' causes the method to return\n"
      "True if the molecule was changed by the neutralize operation, and False\n"
      "otherwise. This is mostly useful if using the API, but for command-line\n"
      "usage (e.g. via obabel) this filters out molecules that are unchanged\n"
      "by the operation and only retains those that are changed.";
  }

  virtual bool WorksWith(OBBase* pOb) const { return dynamic_cast<OBMol*>(pOb) != nullptr; }
  virtual bool Do(OBBase* pOb, const char* OptionText=nullptr, OpMap* pOptions=nullptr, OBConversion* pConv=nullptr);
  bool NoNegativelyChargedNbr(OBAtom *atm);
  bool NoPositivelyChargedNbr(OBAtom *atm);
};

/////////////////////////////////////////////////////////////////
OpNeutralize theOpNeutralize("neutralize"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpNeutralize::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  // The algorithm assumes that hydrogens are suppressed
  pmol->DeleteHydrogens();

  const char *p = OptionText;
  // If "changed" is given as an option, then return true if the molecule was altered
  // by neutralizing some charges, and false if left unchanged.
  bool report_changes = (p && p[0] == 'c' && p[1] == 'h' && p[2] == 'a' && p[3] == 'n' &&
                              p[4] == 'g' && p[5] == 'e' && p[6] == 'd' && p[7] == '\0');
  
  bool changed = false;
  FOR_ATOMS_OF_MOL(atmit, pmol) {
    OBAtom* atm = &*atmit;
    int chg = atm->GetFormalCharge();
    unsigned char hcount;
    switch(chg) {
    case 1:
      hcount = atm->GetImplicitHCount();
      if (hcount >= 1 && NoNegativelyChargedNbr(atm)) {
        atm->SetFormalCharge(0);
        atm->SetImplicitHCount(hcount - 1);
        changed = true;
      }
      break;
    case -1:
      hcount = atm->GetImplicitHCount();
      if (NoPositivelyChargedNbr(atm)) {
        atm->SetFormalCharge(0);
        atm->SetImplicitHCount(hcount + 1);
        changed = true;
      }
      break;
    }
  }

  return report_changes ? changed : true;
}

bool OpNeutralize::NoNegativelyChargedNbr(OBAtom *atm)
{
  FOR_NBORS_OF_ATOM(nbr, atm) {
    if (nbr->GetFormalCharge() < 0)
      return false;
  }
  return true;
}

bool OpNeutralize::NoPositivelyChargedNbr(OBAtom *atm)
{
  FOR_NBORS_OF_ATOM(nbr, atm) {
    if (nbr->GetFormalCharge() > 0)
      return false;
  }
  return true;
}

}//namespace
