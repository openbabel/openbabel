/**********************************************************************
cansmi.cpp - OBDescriptor class accessing Canonical SMILES

Copyright (C) 2009 by Chris Morley

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
#include <openbabel/oberror.h>
#include <openbabel/tokenst.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>

using namespace std;
namespace OpenBabel
{
class CanSmiles : public OBDescriptor
{
public:
  CanSmiles(const char* ID, bool noStereo) : OBDescriptor(ID), _noStereo(noStereo){};
  virtual const char* Description()
  {
    return _noStereo ? "Canonical SMILES without isotopes or stereo" : "Canonical SMILES";
  };
  virtual bool Compare(OBBase* pOb, istream& optionText, bool noEval);
  virtual double GetStringValue(OBBase* pOb, std::string& svalue, std::string* = NULL);
private:
  bool _noStereo;
};

bool CanSmiles::Compare(OBBase* pOb, istream& optionText, bool noEval)
{
  string can;
  GetStringValue(pOb, can);
  return CompareStringWithFilter(optionText, can, noEval);
}

double CanSmiles::GetStringValue(OBBase* pOb, std::string& svalue, std::string*)
{
  OBConversion conv;
  conv.AddOption("n"); //no name
  if(_noStereo)
    conv.AddOption("i");
  if(conv.SetOutFormat("can"))
    svalue = conv.WriteString(pOb);
  else
    obErrorLog.ThrowError(__FUNCTION__, "SmilesFormat is not loaded" , obError);
  Trim(svalue);

  return std::numeric_limits<double>::quiet_NaN();
}

//Make a global instance
CanSmiles theCanSmiles("cansmi", false);
CanSmiles theCanSmilesNS("cansmiNS", true);

//**************************************************************

} //namespace
