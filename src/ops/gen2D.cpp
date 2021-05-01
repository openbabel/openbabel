/**********************************************************************
gen2d.cpp - A OBOp for generation of 2D coordinates

Copyright (C) 2007,2008 by Sergei V. Trepalin sergey_trepalin@chemical-block.com
Copyright (C) 2007,2008 by Andrei Gakh andrei.gakh@nnsa.doe.gov
          (C) 2007 by Chris Morley

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
#include <iostream>
#include<openbabel/op.h>
#include<openbabel/mol.h>
#include <openbabel/mcdlutil.h>
#include <openbabel/stereo/stereo.h>

namespace OpenBabel
{

class OpGen2D : public OBOp
{
public:
  OpGen2D(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return
    "Generate 2D coordinates\n"
    "Trepalin, S. V.; Yarkov, A. V.; Pletnev, I. V.; Gakh, A.A."
    "A Java Chemical Structure Editor Supporting the"
    "Modular Chemical Descriptor Language (MCDL)."
    "Molecules, 2006, 11, 219-231"; }

  virtual bool WorksWith(OBBase* pOb) const { return dynamic_cast<OBMol*>(pOb) != nullptr; }
  virtual bool Do(OBBase* pOb, const char* OptionText=nullptr, OpMap* pOptions=nullptr, OBConversion* pConv=nullptr);
};

/////////////////////////////////////////////////////////////////
OpGen2D theOpGen2D("gen2D"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpGen2D::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  // If we are coming from a 0D structure, then we need to perceive the cis/trans
  // bonds *now*, before adding the coordinates, to mark unspecified cis/trans
  // as such. Otherwise, when writing a MOL file it will be missing the '3', or
  // similarly when depicting it would be presented as specified.
  // Really, all we need to do is handle the cis/trans bond cases.
  // However, the current API requires us to also reperceive the tet stereo,
  // and to remove any stereo that is not real.

  if (pmol->GetDimension() == 0) {
    pmol->UnsetFlag(OB_CHIRALITY_MOL);
    StereoFrom0D(pmol);
  }

  generateDiagram(pmol);
  pmol->SetDimension(2);

  return true;
}
}//namespace
