/**********************************************************************
substructdata.cpp - Color substructures
Copyright (C) 2008 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>

namespace OpenBabel
{

using namespace std;
class OpSubstructData : public OBOp
{
public:
  OpSubstructData(const char* ID) : OBOp(ID, false){};
  const char* Description()
  { return 
  "<SMARTS> Match and color substructure\n"
    "Similar to -s option, but makes matched atoms and bonds green.\n"
    "OpSubstructureData::AddDataToSubstructure() can add any OBPairData\n"
    "to any set of atoms and the bonds which join them.\n";
  }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*);
  static  bool AddDataToSubstruct(OBMol* pMol, const std::vector<int>& atomIndxs, 
                           const std::string& attribute, const std::string& value);
};

/////////////////////////////////////////////////////////////////
OpSubstructData theOpSubstructData("ss"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpSubstructData::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;
  OBSmartsPattern sp;
  sp.Init(OptionText);
  if(!sp.Match(*pmol))
    return false; //do not convert a non-matching mol

  //color the substructure green
  return AddDataToSubstruct(pmol, sp.GetMapList()[0], "color", "green");
}

/**
Adds an OBPairData object to each atom and bond in a substructure.
The substructure's atoms are specified in an input parameter, a
vector of atom indx; the bonds are those in the molecule that join
these atoms. The attribute and value of the OBPairObject (the same
for all the added objects) are specified as parameters. Useful for
adding color to a substructure. Understood by cmlformat and svgformat.

To access this function:
OBOp* pOp = OBOp::FindType("ss");
if(pOp)
  dynamic_cast<OpSubstructData*>(pOp)->AddDataToSubstruct(mol, atomvec, "color", "green");
**/
bool OpSubstructData::AddDataToSubstruct(OBMol* pmol,
        const std::vector<int>& atomIdxs, 
        const std::string& attribute,
        const std::string& value)
{
  //Add data to atoms
  for(int j=0; j<atomIdxs.size(); ++j)
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