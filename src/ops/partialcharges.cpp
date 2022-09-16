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
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include<openbabel/chargemodel.h>
#include <openbabel/obconversion.h>

#include <cstring>
namespace OpenBabel
{

class OpPartialCharge : public OBOp
{
public:
  OpPartialCharge(const char* ID) : OBOp(ID, false) {
    OBConversion::RegisterOptionParam(ID, nullptr, 1, OBConversion::GENOPTIONS);
  }

  const char* Description() override { return "<method> Calculate partial charges by specified method"; }

  bool WorksWith(OBBase* pOb) const override { return dynamic_cast<OBMol*>(pOb) != nullptr; }
  bool Do(OBBase* pOb, const char* OptionText=nullptr, OpMap* pOptions=nullptr,
      OBConversion* pConv=nullptr) override;

  OBChargeModel *_pChargeModel;
};

/////////////////////////////////////////////////////////////////
OpPartialCharge theOpPartialCharge("partialcharge"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpPartialCharge::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv)
{
  char *arg = nullptr;
  const char *tok1= nullptr;
  const char *tok2= nullptr;
  OpMap::const_iterator iter;
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  bool print = false;

  if(!pmol)
    return false;

  iter = pmap->find("print");
  if(iter!=pmap->end())
    print=true;

	if( OptionText ) {
		arg = strdup( OptionText );
		tok1 = strtok( arg, ":" );
		tok2 = strtok( nullptr, "\0" );
	}
	else {
		tok1 = OptionText;
	}

  _pChargeModel = OBChargeModel::FindType(tok1);


  if(!_pChargeModel)
    {
      obErrorLog.ThrowError(__FUNCTION__,
                            std::string("Unknown charge model ") + tok1, obError, onceOnly);
			return false;
    }

  bool success = _pChargeModel->ComputeCharges(*pmol, tok2);

  if (print) {
    // print them on stdout
    FOR_ATOMS_OF_MOL(atom, pmol) {
      std::cout << atom->GetPartialCharge() << '\n';
    }
    std::cout << std::endl; // extra blank and flush the stream
  }

  return success;
}
}//namespace
