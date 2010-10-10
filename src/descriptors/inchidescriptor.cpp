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
  InChIFilter(const char* ID) : OBDescriptor(ID){};
  virtual const char* Description(){return "IUPAC InChI identifier";};
  virtual bool Compare(OBBase* pOb, istream& optionText, bool noEval, std::string* param=NULL);
  virtual double GetStringValue(OBBase* pOb, std::string& svalue, std::string* param=NULL);
//  virtual bool LessThan(OBBase* pOb1, OBBase* pOb2);
//  virtual void Init(){cache.clear();}
//  void GetCachedValue(OBBase* pOb, string& s);
  virtual bool Order(std::string s1, std::string s2)
  {
    InChIFormat::InchiLess f;
    return f(s1, s2);
  }
private:
  map<OBBase*, string> cache;
};

double InChIFilter::GetStringValue(OBBase* pOb, std::string& svalue, std::string*)
{
  OBConversion conv;
  conv.AddOption("w");//suppress trivial warnings

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
  bool matchornegate = ReadStringFromFilter(optionText, InchiFilterString);
  if(noEval)
    return false;
  GetStringValue(pOb, inchi);
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
  bool ret = inchi.compare(inchipos+1, len, InchiFilterString, filterpos, len)==0;

  if(!matchornegate)
    ret = !ret;
  return ret;
}

/*
///Compare InChI strings in a chemically sensible way
// See InChIFormat::InchiLess()
bool InChIFilter::LessThan(OBBase* pOb1, OBBase* pOb2)
{
  string s1, s2;
  //There is a cache, indexed by the molecule pointer, to reduce
  //the number of InChI evaluations. This is potentially large when sorting.
  GetCachedValue(pOb1, s1);
  GetCachedValue(pOb2, s2);
  InChIFormat::InchiLess f;
  return f(s1, s2);
}
void InChIFilter::GetCachedValue(OBBase* pOb, string& s)
{
  map<OBBase*, string>::iterator itr;
  itr = cache.find(pOb);
  if(itr!=cache.end())
    s = itr->second;
  else
  {
    GetStringValue(pOb, s);
    cache[pOb] =s;
  }
}
*/
//**********************************************
//Make a global instance
InChIFilter theInChIFilter("InChI");

}//namespace

