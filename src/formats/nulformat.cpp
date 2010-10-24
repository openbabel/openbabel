/**********************************************************************
Copyright (C) 2010 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#include <openbabel/base.h>

namespace OpenBabel
{

class NulFormat : public OBFormat
{
public:
  NulFormat() { OBConversion::RegisterFormat("nul",this); }

  virtual const char* Description() { return "Outputs nothing"; };

  virtual unsigned int Flags() { return NOTREADABLE; };

  virtual bool WriteChemObject(OBConversion* pConv)
  {
    delete pConv->GetChemObject();
    return true;
  }
  bool WriteMolecule(OBBase* /*pOb*/, OBConversion* /*pConv*/) { return false; }
};

NulFormat theNulFormat;

} //namespace OpenBabel

