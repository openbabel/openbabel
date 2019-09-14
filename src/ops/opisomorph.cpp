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
#include <openbabel/parsmart.h>
#include <openbabel/isomorphism.h>
#include "opisomorph.h"
#include <openbabel/generic.h>
#include <cstdlib>
#include <algorithm>

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
  for(unsigned int j=0; j<atomIdxs.size(); ++j)
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

/**
@since version 2.3
Deletes all atoms except those in @p atomIndxs
**/
bool ExtractSubstruct(OBMol* pmol, const std::vector<int>& atomIdxs)
{
  //Erase from the top to avoid invalidating the remaining ones
  for(int i = pmol->NumAtoms(); i; --i)
    if(find(atomIdxs.begin(),atomIdxs.end(), i)==atomIdxs.end())
      pmol->DeleteAtom(pmol->GetAtom(i));
  return true;
}

//*****************************************************

bool MakeQueriesFromMolInFile(vector<OBQuery*>& queries, const std::string& filename, int* pnAtoms, bool noH)
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
      return false;

    if(noH)
      patternMol.DeleteHydrogens();

    do
    {
      *pnAtoms = patternMol.NumHvyAtoms();
      queries.push_back(CompileMoleculeQuery(&patternMol));
    }while(patternConv.Read(&patternMol));
    return true;
}

const char* OpNewS::Description()
{
  return "Isomorphism filter(-s, -v options replacement)(not displayed in GUI)\n"
    "This enhanced version can take a SMARTS parameter, for example:\n"
    "      babel in.smi -s \"c1ccccc1[#6] green\" out.cml \n"
    "With -s, only molecules matching the SMARTS are converted.\n"
    "With -v, only molecules NOT matching the SMARTS are converted.\n"
    "The optional second parameter causes the matched substructure to be\n"
    "colored if it is a color name like ``green`` or a hex value like\n"
    "``#8dcb70``. The coloring is recognized by SVGFormat and CMLFormat.\n \n"

    "The first parameter can also be a filename with an extension that\n"
    "can be interpreted as a file format:\n"
    "    -s \"pattern.mol exact\"\n"
    "A molecule in the file is used in an isomorphism test with the default\n"
    "matching: bonds by aromaticity or order, atoms only by atomic\n"
    "number. Explicit hydrogen atoms in this molecule are matched like\n"
    "any other atom, unless there is a parameter ``noH``.\n"
    "If the pattern file contains more than one molecule, the test is\n"
    "an OR of them, i.e. with -s, a molecule is converted (and with -v\n"
    "is excluded) if it matches ANY of the pattern molecules.\n"
    "Multiple color parameters can be specified and the coloring in the\n"
    "converted molecule corresponds to the first pattern molecule matched,\n"
    "or the last color if there are fewer colors than pattern molecules.\n \n"

    "If the last parameter is ``showall``, all molecules are shown, even if\n"
    "they do not match. This allows the -s option to be used for highlighting.\n \n"

    "If the second parameter is ``exact`` only exact matches are converted.\n"
    "If the second parameter is ``extract`` all the atoms in the converted\n"
    "molecule are deleted except for those matched. Since these retain their\n"
    "coordinates, this can be used to prepare display templates.\n\n"

    "With SMARTS matching only, the number of unique occurrences in a molecule\n"
    "can be specified in the second parameter, e.g.\n"
    "    -s c1ccccc1 2   which matches if there are exactly two benzene rings\n"
    " or -s c1ccccc1 >2  which matches if there are more than two.\n"
    "(<2 also works.) The color of the substructure can be in the 3rd parameter.\n \n"

    "In the GUI (or on the commandline as an alternative to using -v) the test\n"
    "can be negated with a ~ before the SMARTS string or file name.\n \n"

    "With the ``babel`` commandline interface, unless the option is at the end\n"
    "of a line, it is necessary to enclose all the parameters together in quotes,\n"
    "as in the first example above, because the -s and -v options are\n"
    "expecting a single parameter. With obabel and the GUI this is not necessary.\n"
    "A command must not have more than a single -s or single -v option.\n"
    "The ``--filter`` option is more flexible.\n\n";
}

/////////////////////////////////////////////////////////////////
OpNewS theOpNewS("s"); //Global instances
OpNewS theOpNewV("v");

//////////////////////////////////////////////////////////////////
bool OpNewS::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;

  // The SMARTS and any other parameters are extracted on the first molecule
  // and stored in the member variables. The parameter is cleared so that
  // the original -s option in transform.cpp is inactive

  //string txt(pmap->find(GetID())->second); // ID can be "s" or "v"

  vector<OBQuery*>::iterator qiter;
  if(OptionText && *OptionText)//(!pConv || pConv->IsFirstInput())
  {
    //Set up on first call
    queries.clear();
    query=NULL;
    nPatternAtoms=0;
    inv=false;

    tokenize(vec, OptionText);
    inv = GetID()[0]=='v';
    if(vec[0][0]=='~')
    {
      inv = true;
      vec[0].erase(0,1);
    }

    //Do not filter out any molecules if there is a parameter "showall";
    //allows -s option to be used for highlighting substructures (--highlight also does this)
    vector<string>::iterator it = std::remove(vec.begin(), vec.end(),"showall");
    showAll = it != vec.end();
    if(showAll)
      vec.erase(it);

    //Store the number of matches required, if as a number in the second parameter, else 0.
    nmatches = 0;
    comparechar = '\0';
    if(vec.size()>1)
    {
      comparechar = vec[1][0];
      if(comparechar=='>' || comparechar=='<')
        vec[1].erase(0,1);
      else
        comparechar = '\0';
      nmatches = atoi(vec[1].c_str());
      if(nmatches) //remove this parameter to still allow coloring
        vec.erase(vec.begin()+1);
    }

    //Interpret as a filename if possible
    MakeQueriesFromMolInFile(queries, vec[0], &nPatternAtoms, strstr(OptionText,"noH"));
    vec.erase(remove(vec.begin(),vec.end(),"noH"),vec.end());//to prevent "noH2" being seen as a color
    
     
    if(queries.empty())
    {
      //SMARTS supplied
    
      // Explicit H in SMARTS requires explicit H in the molecule.
      // Calling AddHydrogens() on a copy of the molecule  is done in parsmart.cpp
      // only when SMARTS contains [H]. Doing more has complications with atom typing,
      // so AddHydrogens here on the molecule (not a copy) when #1 detected.
      addHydrogens = (vec[0].find("#1]")!=string::npos);

      // If extra target mols have been supplied, make a composite SMARTS
      // to test for any of the targets.
      if(ExtraMols.size()>0)
      {
        for(unsigned i=0;i<ExtraMols.size();++i)
        {
          OBConversion extraConv;
          extraConv.AddOption("h");
          if(!extraConv.SetOutFormat("smi"))
            return false;
          // Add option which avoids implicit H being added to the SMARTS.
          // The parameter must be present but can be anything.
          extraConv.AddOption("h",OBConversion::OUTOPTIONS, "X");
          xsmarts += ",$(" + extraConv.WriteString(ExtraMols[i], true) + ")";
        }
      }

      string ysmarts = xsmarts.empty() ? vec[0] : "[$(" + vec[0] + ")" + xsmarts +"]";
      xsmarts.clear();
      if(!sp.Init(ysmarts))
      {
        string msg = ysmarts + " cannot be interpreted as either valid SMARTS "
          "or the name of a file with an extension known to OpenBabel "
          "that contains one or more pattern molecules.";
        obErrorLog.ThrowError(__FUNCTION__, msg, obError, onceOnly);
        delete pmol;
        pmol = NULL;
        pConv->SetOneObjectOnly(); //stop conversion
        return false;
      }
    }
    else
    {
      // Target is in a file. Add extra targets if any supplied
      for(unsigned i=0;i<ExtraMols.size();++i)
        queries.push_back(CompileMoleculeQuery(static_cast<OBMol*>(ExtraMols[i])));
      ExtraMols.clear();
    }

    if(vec.size()>1 && vec[1]=="exact")
    {
      if(queries.empty())
      {
        //Convert SMARTS to SMILES to count number of atoms
        OBConversion conv;
        OBMol patmol;
        if(!conv.SetInFormat("smi") || !conv.ReadString(&patmol, vec[0]))
        {
          obErrorLog.ThrowError(__FUNCTION__, "Cannot read the parameter of -s option, "
          "which has to be valid SMILES when the exact option is used.", obError, onceOnly);
          delete pmol;
          if(pConv)
            pConv->SetOneObjectOnly(); //stop conversion
          return false;
        }
        nPatternAtoms = patmol.NumHvyAtoms();
      }
    }
    else
      nPatternAtoms = 0;

    //disable old versions
    if(pConv)
      pConv->AddOption(GetID(), OBConversion::GENOPTIONS, "");
  }

  bool match = false;
  //These are a vector of each mapping, each containing atom indxs.
  vector<vector<int> > vecatomvec;
  vector<vector<int> >* pMappedAtoms = NULL;

  if(nPatternAtoms)
    if(pmol->NumHvyAtoms() != nPatternAtoms)
      return false;

  unsigned int imol=0; //index of mol in pattern file
  if(!queries.empty()) //filename supplied
  {
    //match is set true if any of the structures match - OR behaviour
    for(qiter=queries.begin();qiter!=queries.end();++qiter, ++imol)
    {
      OBIsomorphismMapper* mapper = OBIsomorphismMapper::GetInstance(*qiter);
      OBIsomorphismMapper::Mappings mappings;
      mapper->MapUnique(pmol, mappings);
      if( (match = !mappings.empty()) ) // extra parens to indicate truth value
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
        break;
      }
    }
  }
  else //SMARTS supplied
  {

    if(addHydrogens)
      pmol->AddHydrogens(false,false);

    if( (match = sp.Match(*pmol)) ) // extra parens to indicate truth value
    {
      pMappedAtoms = &sp.GetMapList();
      if(nmatches!=0)
      {
        int n = sp.GetUMapList().size();
        if(comparechar=='>')      match = (n > nmatches);
        else if(comparechar=='<') match = (n < nmatches);
        else                      match = (n == nmatches);
      }
    }
  }

  if((!showAll && (!match && !inv)) || (match && inv))
  {
    //delete a non-matching mol
    delete pmol;
    pmol = NULL;
    return false;
  }

  if(match)
    //Copy the idxes of the first match to a member variable so that it can be retrieved from outside
    firstmatch.assign(pMappedAtoms->begin()->begin(), pMappedAtoms->begin()->end());
  else
    firstmatch.clear();

  if(match && !inv && vec.size()>=2 && !vec[1].empty() && !nPatternAtoms)
  {
    vector<vector<int> >::iterator iter;

    if (vec[1]=="extract" || (vec.size()>3 && vec[2]=="extract"))
    {
      //Delete all unmatched atoms. Use only the first match
      ExtractSubstruct(pmol, *pMappedAtoms->begin());
      return true;
    }

    // color the substructure if there is a second parameter which is not "exact" or "extract" or "noH"
    // with multiple color parameters use the one corresponding to the query molecule, or the last
    if(imol>vec.size()-2)
      imol = vec.size()-2;
    for(iter=pMappedAtoms->begin();iter!=pMappedAtoms->end();++iter)//each match
       AddDataToSubstruct(pmol, *iter, "color", vec[imol+1]);
    return true;
  }

  if(pConv && pConv->IsLast())
  {
    for(qiter=queries.begin();qiter!=queries.end();++qiter)
      delete *qiter;
    queries.clear();
  }
  return true;
}

bool OpNewS::ProcessVec(std::vector<OBBase*>& Extravec)
{
  //Adds extra target molecules (see FastSearchFormat)
  ExtraMols = Extravec;
  return true;
}

}//namespace
