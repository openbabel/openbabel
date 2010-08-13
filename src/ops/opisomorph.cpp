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
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>

namespace OpenBabel
{

using namespace std;

/**
@since version 2.3
Adds an OBPairData object to each atom and bond in a substructure.
The substructure's atoms are specified in an input parameter, a
vector of atom indx; the bonds are those in the molecule that join
these atoms. The attribute and value of the OBPairObject (the same
for all the added objects) are specified as parameters. 
**/
bool AddDataToSubstruct(OBMol* pmol,
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

//*****************************************************

OBQuery* MakeQueryFromMolInFile(const std::string& filename, int* pnAtoms)
{
    OBMol patternMol;
    patternMol.SetIsPatternStructure();
    OBConversion patternConv;
    OBFormat* pFormat;
    //Need to distinguish between filename and SMARTS. Not infallable...
    if( filename.empty() ||
        filename.find('.')==string::npos ||
        !(pFormat = patternConv.FormatFromExt(filename.c_str())) ||
        !patternConv.SetInFormat(pFormat) ||
        !patternConv.ReadFile(&patternMol, filename) ||
        patternMol.NumAtoms()==0)
      return NULL;

    *pnAtoms = patternMol.NumHvyAtoms();
    return CompileMoleculeQuery(&patternMol);
}

//*****************************************************
class OpNewS : public OBOp
{
public:
  OpNewS(const char* ID) : OBOp(ID, false){}
  const char* Description()
  { 
    return "Isomorphism filter (-s option replacement)(not displayed in GUI)\n"
      "This enhanced version can take a SMARTS parameter, for example:\n"
      "      babel in.smi -s \"c1ccccc1[#6] green\" out.cml \n"
      "Only molecules matching the SMARTS are converted. The optional second\n"
      "parameter causes the matched substructure to be colored. The coloring is\n"
      "recognized by SVGFormat and CMLFormat.\n"
      "If the second parameter is ``exact`` only exact matches are converted.\n"
      "If the SMARTS starts with a ``~`` character, only non-matching molecules \n"
      "are converted, e.g.\n"
      "    -s ~c1ccccc1[#6]\n"
      "This is an alternative to the -v option.\n\n"
      
      "The first parameter can also be a filename with an extension that\n"
      "can be interpreted as a file format:\n"
      "    -s \"pattern.mol exact\"\n"
      "The molecule in the file is used in an isomorphism test with the default\n"
      "matching: bonds by aromaticity or order, atoms only by atomic number.\n"
      "Explicit hydrogen atoms in this molecule are matched like any other atom.\n"
      "The test can be negated with a ~ before the file name.\n\n"

      "With the babel commandline interface, unless the option is at the end of\n"
      "a line, it is necessary to enclose all the parameters together in quotes,\n"
      "as in the first example above, because the -s option is expecting a\n"
      "single parameter. With obabel and the GUI this is not necessary.\n";
  }
  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*);
};

/////////////////////////////////////////////////////////////////
OpNewS theOpNewS("s"); //Global instance

bool OpNewS::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  // The SMARTS and any other parameters are extracted on the first molecule
  // and stored in the static variable vec. The parameter is cleared so that:
  // (a) the original -s option in transform.cpp is inactive, and
  // (b) the parsing does not have to be done again for multi-molecule files
  string txt(pmap->find("s")->second);
  static vector<string> vec;
  static bool inv;
  static int nPatternAtoms;
  static OBQuery* query;
  if(!txt.empty())
  {
    //Set up on first call
    tokenize(vec, txt);
    inv = false;
    if(vec[0][0]=='~')
    {
      inv = true;
      vec[0].erase(0,1);
    }

    //Interpret as a filename if possible
    query = MakeQueryFromMolInFile(vec[0], &nPatternAtoms);

    if(vec.size()>1 && vec[1]=="exact")
    {
      if(!query)
      {
        //Convert SMARTS to SMILES to count number of atoms
        OBConversion conv;
        OBMol patmol;
        if(!conv.SetInFormat("smi") || !conv.ReadString(&patmol, vec[0]))
        {
          obErrorLog.ThrowError(__FUNCTION__, "Cannot read the parameter of -s option, "
          "which has to be valid SMILES when the exact option is used.", obError, onceOnly);
          delete pmol;
          pConv->SetOneObjectOnly(); //stop conversion
          return false;
        }
        nPatternAtoms = patmol.NumHvyAtoms();   
      }
    }
    else
      nPatternAtoms = 0;

    pConv->AddOption("s", OBConversion::GENOPTIONS, "");
  }

  bool match;
  //These are a vector of each mapping, each containing atom indxs.
  vector<vector<int> > vecatomvec;
  vector<vector<int> >* pMappedAtoms = NULL;
  OBSmartsPattern sp;

  if(nPatternAtoms)
    if(pmol->NumHvyAtoms() != nPatternAtoms)
      return false;

  if(query) //filename supplied
  {
    OBIsomorphismMapper* mapper = OBIsomorphismMapper::GetInstance(query);
    OBIsomorphismMapper::Mappings mappings = mapper->MapUnique(pmol);
    if(match = !mappings.empty())
    {
      OBIsomorphismMapper::Mappings::iterator ita;
      OBIsomorphismMapper::Mapping::iterator itb;
      for(ita=mappings.begin(); ita!=mappings.end();++ita)//each mapping
      {
        vector<int> atomvec;
        for(itb=ita->begin(); itb!=ita->end();++itb)//each atom index
          atomvec.push_back(itb->second+1);
        vecatomvec.push_back(atomvec);
        atomvec.clear();
      }
      pMappedAtoms = &vecatomvec;
    }
  }
  else //SMARTS supplied
  {
    if(!sp.Init(vec[0]))
    {
      string msg = vec[0] + " cannot be interpreted as either valid SMARTS "
        "or the name of a file with an extension known to OpenBabel "
        "that contains a pattern molecule.";
      obErrorLog.ThrowError(__FUNCTION__, msg, obError, onceOnly);
      delete pmol;
      pmol = NULL;
      pConv->SetOneObjectOnly(); //stop conversion
      return false;
    }

    if(match = sp.Match(*pmol))
      pMappedAtoms = &sp.GetMapList();
  }

  if((!match && !inv) || (match && inv))
  {
    //delete a non-matching mol
    delete pmol;
    pmol = NULL;
    return false; 
  }

  if(!inv && vec.size()>=2 && !vec[1].empty()) // color the substructure
  {
    vector<vector<int> >::iterator iter;
    for(iter=pMappedAtoms->begin();iter!=pMappedAtoms->end();++iter)//each match
       AddDataToSubstruct(pmol, *iter, "color", vec[1]);
    return true;
  }

  if(pConv && pConv->IsLast())
    delete query;

  return true;
}

}//namespace
