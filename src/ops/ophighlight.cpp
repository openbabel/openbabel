/**********************************************************************
ophighlight.cpp
Copyright (C) 2012 by Noel O'Boyle
Based on opisomorph.cpp Copyright (C) 2010 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include<openbabel/op.h>
#include <openbabel/generic.h>
#include <openbabel/bond.h>
#include <openbabel/mol.h>
#include <algorithm>

namespace OpenBabel
{

using namespace std;




//*****************************************************

class OpHighlight : public OBOp
{
public:
  OpHighlight(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return
    "<param> Highlight substructures in 2D depictions\n"
    "Usage: --highlight \"SMARTS1 color1 [SMARTS2 color2 ...]\"\n"
    "\n"
    "Valid colors are black, gray, white, red, green, blue, yellow,\n"
    "                 cyan, purple, teal and olive.\n"
    "Additional colors may be specified as hexadecimal RGB values\n"
    "preceded by #.\n\n"
    "The following will color the phenyl group green, and the\n"
    "carboxyl group orange:\n"
    "  obabel -:\"c1ccccc1CCC(=O)O\" -O mol.svg\n"
    "     --highlight \"c1ccccc1 green C(=O)O #FFA500\"";
  }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText=NULL, OpMap* pOptions=NULL, OBConversion* pConv=NULL);
private:
  bool AddDataToSubstruct(OBMol* pmol,
        const std::vector<int>& atomIdxs,
        const std::string& attribute,
        const std::string& value);
};

/////////////////////////////////////////////////////////////////
OpHighlight theOpHighlight("highlight"); //Global instance

//////////////////////////////////////////////////////////////////
bool OpHighlight::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;
  std::vector<std::string> vec;
  tokenize(vec, OptionText);
  
  // Jump in twos over the parameters! 
  // "SMARTS1 color1 SMARTS2 color2 ..."
  for(unsigned int vecIdx = 0; vecIdx < vec.size(); vecIdx += 2)
  {
    
    std::string smarts = vec[vecIdx];
    if(vecIdx + 1 == vec.size())
    {
      string msg = "No color specified for SMARTS string: " + smarts;
      obErrorLog.ThrowError(__FUNCTION__, msg, obError, onceOnly);
      delete pmol;
      pmol = NULL;
      pConv->SetOneObjectOnly(); //stop conversion
      return false;
    }
    std::string color = vec[vecIdx + 1];

    bool match = false;
    //These are a vector of each mapping, each containing atom indxs.
    vector<vector<int> > vecatomvec;
    vector<vector<int> >* pMappedAtoms = NULL;
    OBSmartsPattern sp;

    // Explicit H in SMARTS requires explicit H in the molecule.
    // Calling AddHydrogens() on a copy of the molecule  is done in parsmart.cpp
    // only when SMARTS contains [H]. Doing more has complications with atom typing,
    // so AddHydrogens here on the molecule (not a copy) when #1 detected.
    bool addHydrogens = (smarts.find("#1]")!=string::npos);

    if(!sp.Init(smarts))
    {
      string msg = smarts + " cannot be interpreted as a valid SMARTS ";
      obErrorLog.ThrowError(__FUNCTION__, msg, obError, onceOnly);
      delete pmol;
      pmol = NULL;
      pConv->SetOneObjectOnly(); //stop conversion
      return false;
    }

    if(addHydrogens)
      pmol->AddHydrogens(false,false);

    if( (match = sp.Match(*pmol)) ) // extra parens to indicate truth value
      pMappedAtoms = &sp.GetMapList();

    if(match)
    {
      vector<vector<int> >::iterator iter;
      for(iter=pMappedAtoms->begin();iter!=pMappedAtoms->end();++iter)//each match
         AddDataToSubstruct(pmol, *iter, "color", color);
    }
  }

  return true;
}

// Taken from opisomorph.cpp
bool OpHighlight::AddDataToSubstruct(OBMol* pmol,
        const std::vector<int>& atomIdxs,
        const std::string& attribute,
        const std::string& value)
{
  //Add data to atoms
  for(unsigned int j = 0; j < atomIdxs.size(); ++j)
  {
    OBAtom* pAtom = pmol->GetAtom(atomIdxs[j]);
    if(!pAtom)
      continue;
    OBPairData* dp = new OBPairData;
    dp->SetAttribute(attribute);
    dp->SetValue(value);
    pAtom->SetData(dp);
  }

  OBBond* pBond;
  vector<OBBond*>::iterator i;
  for(pBond = pmol->BeginBond(i); pBond; pBond = pmol->NextBond(i))
  {
    //Add data to bond if it joins two atoms in list
    if(count(atomIdxs.begin(), atomIdxs.end(), pBond->GetBeginAtomIdx())
        && count(atomIdxs.begin(), atomIdxs.end(), pBond->GetEndAtomIdx()))
    {
      OBPairData* dp = new OBPairData;
      dp->SetAttribute(attribute);
      dp->SetValue(value);
      pBond->SetData(dp);
    }
  }
  return true;
}

}//namespace

