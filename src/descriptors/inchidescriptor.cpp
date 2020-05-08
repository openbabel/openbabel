/**********************************************************************
inchifilter.cpp - A descriptor giving an InChI string

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
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
#include <openbabel/inchiformat.h>

using namespace std;
namespace OpenBabel
{

class InChIFilter : public OBDescriptor
{
public:
  InChIFilter(const char* ID, bool useKey=false) : OBDescriptor(ID), bKey(useKey) {};
  virtual const char* Description()
  {
    return bKey ? "InChIKey" : "IUPAC InChI identifier";
  }
  virtual bool Compare(OBBase* pOb, istream& optionText, bool noEval, std::string* param=nullptr);
  virtual double GetStringValue(OBBase* pOb, std::string& svalue, std::string* param=nullptr);
  virtual bool Order(std::string s1, std::string s2)
  {
    InChIFormat::InchiLess f;
    return f(s1, s2);
  }
private:
  bool bKey;
};

double InChIFilter::GetStringValue(OBBase* pOb, std::string& svalue, std::string*)
{
  OBConversion conv;
  conv.AddOption("w");//suppress trivial warnings
  if(bKey)
    conv.AddOption("K");
  if(conv.SetOutFormat("inchi"))
    svalue = conv.WriteString(pOb);
  else
    obErrorLog.ThrowError(__FUNCTION__, "InChIFormat is not loaded" , obError);
  Trim(svalue);

  return std::numeric_limits<double>::quiet_NaN();
}

bool InChIFilter::Compare(OBBase* pOb, istream& optionText, bool noEval, std::string*)
{
  string InchiFilterString, inchi;
  string::size_type filterpos=0, inchipos, len;
  bool ret;
  bool matchornegate = ReadStringFromFilter(optionText, InchiFilterString);
  if(noEval)
    return false;
  GetStringValue(pOb, inchi);
  if(!bKey)
  {
    inchipos = inchi.find('/');

    //See if filterstring starts with "InChI=1/"
    if(InchiFilterString.find(inchi.substr(0,inchipos))==0)
      filterpos=inchipos+1;
    //If filterstring starts with a number, set filterpos after the next '/'(for pasted InChIs)
    if(isdigit(InchiFilterString[0]))
      filterpos=InchiFilterString.find('/')+1;

    //Considering only the significant parts,
    //compare InChI and filter string, only to length of filter string
    len = InchiFilterString.size() - filterpos;
    ret = inchi.compare(inchipos+1, len, InchiFilterString, filterpos, len)==0;
  }
  else
    // compare up to length of the provided filter string,
    // so can match ignoring stereo by providing only the first part of the key.
    ret = (inchi.compare(0,InchiFilterString.size(),InchiFilterString)==0);

  if(!matchornegate)
    ret = !ret;
  return ret;
}

//**********************************************
//Make global instances
InChIFilter theInChIFilter("InChI");
InChIFilter keyInChIFilter("InChIKey",true);
//**********************************************

}//namespace

