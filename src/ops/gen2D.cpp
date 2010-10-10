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

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
};

/////////////////////////////////////////////////////////////////
OpGen2D theOpGen2D("gen2D"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpGen2D::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  generateDiagram(pmol);
  pmol->SetDimension(2);

  return true;
}
}//namespace
