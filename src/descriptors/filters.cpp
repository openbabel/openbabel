/**********************************************************************
filters.cpp - Some classes derived from OBDescriptor
 
Copyright (C) 2007 by Chris Morley
 
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
#include <openbabel/oberror.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>

using namespace std;
namespace OpenBabel
{
//**************************************************************

  class MWFilter : public OBDescriptor
{ 
public:
  MWFilter(const char* ID) : OBDescriptor(ID){};
  virtual const char* Description(){return "Molecular Weight filter";};
  double Predict(OBBase* pOb)
  {
    OBMol* pmol = dynamic_cast<OBMol*> (pOb);
    if(!pmol)
      return 0;
    return pmol->GetMolWt();
  }
};
//Make a global instance
MWFilter theMWFilter("MW");

//**************************************************************

/// Class to wrap SMARTS comparisons for use with --filter option
class SmartsFilter : public OBDescriptor
{
public:
  SmartsFilter(const char* ID) : OBDescriptor(ID){};
  virtual const char* Description(){return "SMARTS filter";};
  virtual bool Compare(OBBase* pOb, istream& optionText, bool noEval);

}; 

/// For interpreting conditions like s!=c1ccccc1CN
/** The descriptor name can be s or smarts and is case independent
    The operator to return true for a match can be:
    one or more spaces, =, ==,  or nothing if the SMARTS string
    starts with a letter.
    To return true for a mismatch the operator is !=
    A space or tab should follow the SMARTS string.
 **/      
bool SmartsFilter::Compare(OBBase* pOb, istream& optionText, bool noEval)
{
  OBMol* pmol = dynamic_cast<OBMol*> (pOb);
  if(!pmol)
    return false;
  
  string smarts;
  bool matchornegate = ReadStringFromFilter(optionText, smarts);
  if(noEval)
    return false;
  OBSmartsPattern sp;
  sp.Init(smarts);
  bool ret = sp.Match(*pmol,true);//single match  
  if(!matchornegate)
    ret = !ret;
  return ret;
}

//Make a global instances with alternative IDs
SmartsFilter firstSmartsFilter("smarts");
SmartsFilter secondSmartsFilter("s");

//**************************************************************

/// Class to filter on molecule title
class TitleFilter : public OBDescriptor
{
public:
  TitleFilter(const char* ID) : OBDescriptor(ID){};
  virtual const char* Description(){return "For comparing a molecule's title";};
  virtual bool Compare(OBBase* pOb, istream& optionText, bool noEval);
  virtual double GetStringValue(OBBase* pOb, std::string& svalue);
};

bool TitleFilter::Compare(OBBase* pOb, istream& optionText, bool noEval)
{
  OBMol* pmol = dynamic_cast<OBMol*> (pOb);
  if(!pmol)
    return false;

  string title(pmol->GetTitle());
  return CompareStringWithFilter(optionText, title, noEval);
}

double TitleFilter::GetStringValue(OBBase* pOb, std::string& svalue)
{
  OBMol* pmol = dynamic_cast<OBMol*> (pOb);
  if(pmol)
    svalue = pmol->GetTitle();
  return std::numeric_limits<double>::quiet_NaN();
}

//Make a global instance
TitleFilter theTitleFilter("title");

//**************************************************************

/// Temporary class to test spin multiplicity
class SpinFilter : public OBDescriptor
{
public:
  SpinFilter(const char* ID) : OBDescriptor(ID){};
  virtual const char* Description(){return "Total Spin Multiplicity";};

  double Predict(OBBase* pOb)
  {
    OBMol* pmol = dynamic_cast<OBMol*> (pOb);
    if(!pmol)
      return 0;
    return pmol->GetTotalSpinMultiplicity();
  }
};

//Make a global instance
SpinFilter theTitleSpin("spinMult");


//**************************************************************

class InChIFilter : public OBDescriptor
{
public:
  InChIFilter(const char* ID) : OBDescriptor(ID){};
  virtual const char* Description(){return "IUPAC InChI identifier";};
  virtual bool Compare(OBBase* pOb, istream& optionText, bool noEval);
  virtual double GetStringValue(OBBase* pOb, std::string& svalue);
};

double InChIFilter::GetStringValue(OBBase* pOb, std::string& svalue)
{
  OBConversion conv;
  if(conv.SetOutFormat("inchi"))
    svalue = conv.WriteString(pOb);
  else
    obErrorLog.ThrowError(__FUNCTION__, "InChIFormat is not loaded" , obError);
  Trim(svalue);
  return std::numeric_limits<double>::quiet_NaN();
}

bool InChIFilter::Compare(OBBase* pOb, istream& optionText, bool noEval)
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

//Make a global instance
InChIFilter theInChIFilter("InChI");

}//namespace

