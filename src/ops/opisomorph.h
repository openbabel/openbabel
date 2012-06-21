/**********************************************************************
opisomorph.cpp - Enhanced -s option
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
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/query.h>
#include <vector>
#include <string>

namespace OpenBabel
{
/**
@since version 2.3
Adds an OBPairData object to each atom and bond in a substructure.
The substructure's atoms are specified in an input parameter, a
vector of atom indx; the bonds are those in the molecule that join
these atoms. The attribute and value of the OBPairObject (the same
for all the added objects) are specified as parameters.
**/
extern bool AddDataToSubstruct(OBMol* pmol, const std::vector<int>& atomIdxs,
        const std::string& attribute, const std::string& value);

/**
@since version 2.3
Deletes all atoms except those in @p atomIndxs
**/
extern bool ExtractSubstruct(OBMol* pmol, const std::vector<int>& atomIdxs);

extern bool MakeQueriesFromMolInFile(std::vector<OBQuery*>& queries
                         , const std::string& filename, int* pnAtoms, bool noH);

//*****************************************************
class OpNewS : public OBOp
{
public:
  OpNewS(const char* ID) : OBOp(ID, false){}
  const char* Description();
  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*);
  std::vector<int> GetMatchAtoms(){ return firstmatch; }
  
private:
  std::vector<std::string> vec;
  bool inv;
  int nPatternAtoms; //non-zero for exact matches
  OBQuery* query;
  std::vector<OBQuery*> queries;
  std::vector<int> firstmatch; //Idxes of first match by SMARTS or OBIsomorphismMapper
  bool showAll;
};

} //namespace

