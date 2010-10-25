/**********************************************************************
partialcharges.cpp - The option --partialcharges (type) adds charges using different models

Copyright(C) 2010 by Geoffrey R. Hutchison

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

#include<openbabel/op.h>
#include<openbabel/mol.h>
#include<openbabel/chargemodel.h>
#include <openbabel/obconversion.h>

namespace OpenBabel
{

class OpPartialCharge : public OBOp
{
public:
  OpPartialCharge(const char* ID) : OBOp(ID, false) {
    OBConversion::RegisterOptionParam(ID, NULL, 1, OBConversion::GENOPTIONS);
  }

  const char* Description(){ return "<method> Calculate partial charges by specified method"; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);

  OBChargeModel *_pChargeModel;
};

/////////////////////////////////////////////////////////////////
OpPartialCharge theOpPartialCharge("partialcharge"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpPartialCharge::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  _pChargeModel = OBChargeModel::FindType(OptionText);
  if(!_pChargeModel)
    {
      obErrorLog.ThrowError(__FUNCTION__,
                            std::string("Unknown charge model ") + OptionText, obError, onceOnly);
      return false;
    }

  return _pChargeModel->ComputeCharges(*pmol);
}
}//namespace
