/**********************************************************************
filters.cpp - Some classes derived from OBDescriptor

Copyright (C) 2007 by Chris Morley

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
#include <openbabel/descriptor.h>
#include <openbabel/fingerprint.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <openbabel/parsmart.h>

using namespace std;
namespace OpenBabel {
//**************************************************************

class MWFilter : public OBDescriptor {
public:
  MWFilter(const char *ID) : OBDescriptor(ID){};
  virtual const char *Description() { return "Molecular Weight filter"; };
  virtual double Predict(OBBase *pOb, string *param = NULL) {
    OBMol *pmol = dynamic_cast<OBMol *>(pOb);
    if (!pmol)
      return 0;
    return pmol->GetMolWt();
  }
};
// Make a global instance
MWFilter theMWFilter("MW");

class RotatableBondsFilter : public OBDescriptor {
public:
  RotatableBondsFilter(const char *ID) : OBDescriptor(ID){};
  virtual const char *Description() { return "Rotatable bonds filter"; };
  virtual double Predict(OBBase *pOb, string *param = NULL) {
    OBMol *pmol = dynamic_cast<OBMol *>(pOb);
    if (!pmol)
      return 0;
    return pmol->NumRotors();
  }
};
// Make a global instance
RotatableBondsFilter theRBFilter("rotors");

//**************************************************************

/// Class to wrap SMARTS comparisons for use with --filter option
class SmartsFilter : public OBDescriptor {
public:
  SmartsFilter(const char *ID) : OBDescriptor(ID){};
  virtual const char *Description() { return "SMARTS filter"; };
  virtual bool Compare(OBBase *pOb, istream &optionText, bool noEval,
                       std::string *param = NULL);
};

/// For interpreting conditions like s!=c1ccccc1CN
/** The descriptor name can be s or smarts and is case independent
    The operator to return true for a match can be:
    one or more spaces, =, ==,  or nothing if the SMARTS string
    starts with a letter.
    To return true for a mismatch the operator is !=
    A space or tab should follow the SMARTS string.
 **/
bool SmartsFilter::Compare(OBBase *pOb, istream &optionText, bool noEval,
                           std::string *) {
  OBMol *pmol = dynamic_cast<OBMol *>(pOb);
  if (!pmol)
    return false;

  string smarts;
  bool matchornegate = ReadStringFromFilter(optionText, smarts);
  if (noEval)
    return false;
  OBSmartsPattern sp;
  if (!sp.Init(smarts))
    return false; // can't initialize the SMARTS, so fail gracefully

  bool ret = sp.Match(*pmol, true); // single match
  if (!matchornegate)
    ret = !ret;
  return ret;
}

// Make a global instances with alternative IDs
SmartsFilter firstSmartsFilter("smarts");
SmartsFilter secondSmartsFilter("s");

//**************************************************************

/// Class to filter on molecule title
class TitleFilter : public OBDescriptor {
public:
  TitleFilter(const char *ID) : OBDescriptor(ID){};
  virtual const char *Description() {
    return "For comparing a molecule's title";
  };
  virtual bool Compare(OBBase *pOb, istream &optionText, bool noEval,
                       std::string *param = NULL);
  virtual double GetStringValue(OBBase *pOb, std::string &svalue,
                                std::string *param = NULL);
  virtual bool LessThan(OBBase *pOb1, OBBase *pOb2);
};

bool TitleFilter::Compare(OBBase *pOb, istream &optionText, bool noEval,
                          std::string *) {
  OBMol *pmol = dynamic_cast<OBMol *>(pOb);
  if (!pmol)
    return false;

  string title(pmol->GetTitle());
  return CompareStringWithFilter(optionText, title, noEval);
}

double TitleFilter::GetStringValue(OBBase *pOb, std::string &svalue,
                                   std::string *) {
  OBMol *pmol = dynamic_cast<OBMol *>(pOb);
  if (pmol)
    svalue = pmol->GetTitle();
  return std::numeric_limits<double>::quiet_NaN();
}

bool TitleFilter::LessThan(OBBase *pOb1, OBBase *pOb2) {
  OBMol *pmol1 = dynamic_cast<OBMol *>(pOb1);
  OBMol *pmol2 = dynamic_cast<OBMol *>(pOb2);
  if (pmol1 == NULL || pmol2 == NULL)
    return false; // as a default to prevent dereferencing NULL pointers

  return strcmp(pmol1->GetTitle(), pmol2->GetTitle()) < 0;
}

// Make a global instance
TitleFilter theTitleFilter("title");

//**************************************************************
class FormulaDescriptor : public OBDescriptor {
public:
  FormulaDescriptor(const char *ID) : OBDescriptor(ID){};
  virtual const char *Description() { return "Chemical formula"; };

  virtual double GetStringValue(OBBase *pOb, std::string &svalue,
                                std::string *param = NULL) {
    OBMol *pmol = dynamic_cast<OBMol *>(pOb);
    if (pmol)
      svalue = pmol->GetSpacedFormula(1, ""); // actually unspaced
    return std::numeric_limits<double>::quiet_NaN();
  }

  virtual bool Compare(OBBase *pOb, istream &optionText, bool noEval,
                       std::string *) {
    string svalue;
    GetStringValue(pOb, svalue);
    return CompareStringWithFilter(optionText, svalue, noEval);
  }
};

FormulaDescriptor TheFormulaDescriptor("formula");

//******************************************************************
/* This descriptor uses a parameter, e.g. popcount(FP4)
class FPCount : public OBDescriptor
{
public:
  FPCount(const char* ID) : OBDescriptor(ID){};
  virtual const char* Description(){return "Count bits set in fingerprint whose
ID is in the parameter";}; virtual double Predict(OBBase* pOb, string*
param=NULL)
  {
    OBMol* pmol = dynamic_cast<OBMol*> (pOb);
    if(!pmol)
      return std::numeric_limits<double>::quiet_NaN();

    Trim(*param);
    OBFingerprint* pFP = OBFingerprint::FindFingerprint(param ? param->c_str() :
NULL); if(!pFP)
    {
      obErrorLog.ThrowError(__FUNCTION__, "Fingerprint type not available",
obError, onceOnly); return std::numeric_limits<double>::quiet_NaN();
    }
    vector<unsigned> FP;
    pFP->GetFingerprint(pmol,FP);
    int count=0;
    for (unsigned i=0;i<FP.size();++i)
    {
      int w = FP[i];
      for(;w;w=w<<1)
        if(w<0)
          ++count;
    }
    return count;
  }
};
//Make a global instance
FPCount theFPCount("popcount");
*/
//**************************************************************
} // namespace OpenBabel
