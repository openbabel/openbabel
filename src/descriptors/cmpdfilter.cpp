/**********************************************************************
cmpdfilter.cpp - OBDescriptor class to

Copyright (C) 2008 by Chris Morley

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
#include <openbabel/descriptor.h>

using namespace std;
namespace OpenBabel
{
/**
Compound filter provides a macro capability:
a descriptor returning a boolean value can be defined in terms of conditional
statements involvingother descriptors or properties(OBPairData probably from SDF).
For example Lipinsky's Rule of Five would have a definition (macrotext in the constructor)
 "HBD<5 HBA1<10 MW<500 logP<5"

An instance can be made in code by defining a global variable:
 CompoundFilter theL5("L5", "HBD<5 HBA1<10 MW<500 logP<5", Lipinski Rule of Five);
passing an ID, the macro string and a description.

Alternatively, if the pluginclass OBDefine has been loaded, CompoundFilters
are definable via an entry in plugindefines.txt:

CompoundFilter
L5                           # ID
HBD<5 HBA1<10 MW<500 logP<5  # definition in terms of other descriptors or properties
Lipinski Rule of Five        # description

For this non-code method to work, at least one instance of a CompoundFilter must have
made using the global variable method. See below for a suitable dummy instance, which
would not be necessary if a more useful instance was defined instead.
**/

class CompoundFilter : public OBDescriptor
{
public:
///Constructor defining ID, the filter string that the descriptor is short for, and a description.
  CompoundFilter(const char* ID, const char* macrotext, const char* descr)
    : OBDescriptor(ID), _descr(descr), _macroText(macrotext){}

  virtual const char* Description()
  {
    static string txt;
    txt = _descr;
    txt += '\n';
    txt += _macroText;
    txt += "\nCompoundFilter is definable";//Entries in plugindefines.txt can start "CompoundFilter"
    return txt.c_str();
  }

///
  virtual CompoundFilter* MakeInstance(const std::vector<std::string>& textlines)
  {
    return new CompoundFilter(textlines[1].c_str(),textlines[2].c_str(),textlines[3].c_str());
  }

///Returns the result of evaluating the conditional expressions in the macrotext
  virtual bool Compare(OBBase* pOb, istream&, bool noEval, string* param)
  {
    stringstream ss;
    ss.str(_macroText);
    return FilterCompare(pOb, ss, noEval);
  }

private:
  const char* _descr;
  const string _macroText;
};

//*********************************************************************
//This global instance is needed to make system aware of the existence
// of CompoundFilter. ID starts with a _ so it is not listed
CompoundFilter dummyCmpFilter("_", "", "dummyCompoundFilter");
//*********************************************************************

}//namespace
