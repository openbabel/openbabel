#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace OpenBabel;


int main()
{
  OBMol mol;
  OBAtom *begin = mol.NewAtom();
  OBAtom *end = mol.NewAtom();
  mol.AddBond(begin->GetIdx(), end->GetIdx(), 1, OB_WEDGE_BOND);
  mol.GetBond(0)->SetWedge();


  OBConversion conv;
  conv.SetInAndOutFormats("mol", "mol");
  
  // write the molecule
  std::ofstream ofs("stereobonds.mol");
  conv.Write(&mol, &ofs);
  ofs.close();

  std::ifstream ifs("stereobonds.mol");
  conv.Read(&mol, &ifs);

  OB_COMPARE( mol.GetBond(0)->IsWedge(), true );
}



