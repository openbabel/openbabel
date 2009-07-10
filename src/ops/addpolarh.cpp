/**********************************************************************
AddPolarH.cpp - The option --AddPolarHAdds  adds hydrogen to polar atoms only.

Copyright(C) 2007 by Chris Morley
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
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

namespace OpenBabel
{

class OpAddPolarH : public OBOp
{
public:
  OpAddPolarH(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return "Adds hydrogen to polar atoms only"; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, OpMap* pmap, const char* OptionText);
};

/////////////////////////////////////////////////////////////////
OpAddPolarH theOpAddPolarH("AddPolarH"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpAddPolarH::Do(OBBase* pOb, OpMap* pmap, const char* OptionText)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  pmol->AddHydrogens(true, false);

  return true;
}
}//namespace
