/**********************************************************************
alias.cpp - implementation of OBGenericData class to hold alias information on atoms
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
#include <sstream>
#include <string>
#include <openbabel/alias.h>
#include <openbabel/obconversion.h>

using namespace std;
namespace OpenBabel
{

  bool AliasData::Expand(OBMol& mol, const unsigned int atomindex)
  {
    /*
    Interprets the alias text and adds atom(s) as appropriate to mol.
    Tries the following in turn until one is sucessful:
    1) If starts with number treat as isotope+element e.g. 2H
    2) Looks up alias in superatom.txt e.g. COOH Pr
    3) Parse as simple formula
    Returns false if none are successful.
    */

    //parse as isotopic atom
    if(isdigit(_alias[0]))
    {
      std::stringstream ss(_alias);
      int iso;
      std::string el;
      ss >> iso >>el;
      if(etab.GetAtomicNum(el.c_str())>0)
      {
        OBAtom* pAtom = mol.GetAtom(atomindex);
        if(!pAtom)
          return false;
        pAtom->SetIsotope(iso);
        pAtom->SetAtomicNum(etab.GetAtomicNum(el.c_str(),iso));
        return true;
      }
    }

    if(FromNameLookup(mol, atomindex))
      return true;

    //Crude implementation of formula parse
    //This should really not alter the molecule until the parsing has been
    //seen to be successful - another method may be better. But at present
    //there is no other method.
    //Only single character element symbols are handled
    //Atom which replaces atomindex is the first non-H 
    //Will parse ND2 DS CH-

    //(Adapt to use old code)
    char* txt = new char[_alias.size()+1];
    strcpy(txt, _alias.c_str());

    if(*txt=='?') //Assume that it is harmless to ignore this alias
    {
      delete[] txt;
      return true;
    }
    if(!isalpha(*txt)) //first char is the element that replaces atomindex
    {
      return false;
      delete[] txt;
    }
    //Swaps any leading H isotope with the first non-H atom
    if(*txt=='H' || *txt=='D' || *txt=='T')
    {
      char* p =txt+1;
      while(*p && *p=='H' && *p=='D' && *p=='T')p++;
      if(*p)
        std::swap(*p, *txt);
    }
    char symb[2];
    symb[0]=*(txt++);
    symb[1]='\0';
    OBAtom* pAtom = mol.GetAtom(atomindex);
    if(!pAtom)
      return false;
    int iso = 0;
    pAtom->SetAtomicNum(etab.GetAtomicNum(symb,iso));
    if(iso)
      pAtom->SetIsotope(iso);
    _expandedatoms.push_back(atomindex);

    while(*txt)
    {
      if(isspace(*txt)) {
        ++txt;
        continue;
      }
      int chg=0;
      if(*txt=='-')
        chg = -1;
      else if(*txt=='+')
        chg = 1;
      if(chg)
      {
        pAtom->SetFormalCharge(pAtom->GetFormalCharge()+chg);//put on central atom e.g. CH-
        ++txt;
        continue;
      }
      if(!isalpha(*txt))
        return false;
      symb[0]=*txt;
      int rep = atoi(++txt);
      if(rep)
        ++txt;
      do //for each rep
      {
        OBAtom* newAtom = mol.NewAtom();
        _expandedatoms.push_back(mol.NumAtoms());
        iso = 0;
        newAtom->SetAtomicNum(etab.GetAtomicNum(symb,iso));
        if(iso)
          newAtom->SetIsotope(iso);

        if (!mol.AddBond(atomindex,mol.NumAtoms(),1,0)) return false;
      }while(--rep>0);
    }
    return true;
  }

bool AliasData::FromNameLookup(OBMol& mol, const unsigned int atomindex)
{
  OBAtom* Xxatom = mol.GetAtom(atomindex);
  if(Xxatom->GetValence()>1)
  {
    obErrorLog.ThrowError(__FUNCTION__, _alias + " is multivalent, which is currently not supported.", obError);
    return false;
  }

  static std::map<std::string, std::string> table;
  if(table.empty())
    LoadFile(table);
  
  map<std::string, std::string>::iterator pos = table.find(_alias);
  if(pos==table.end())
  {
    obErrorLog.ThrowError(__FUNCTION__, "Alias " + _alias + " was not recognized.\n Output may not be correct.", obError, onceOnly);
    return false;
  }
  else if(mol.Has2D())//mol.GetDimension()>=2)
    obErrorLog.ThrowError(__FUNCTION__,
    "The current implementation does not assign coordinates to the "
    "atoms in the expanded alias, and so is usable only for 0D molecules", obWarning, onceOnly);


  //Convert SMILES of alias
  OBConversion conv;
  OBMol obAlias;
  obAlias.SetIsPatternStructure();
  if(conv.SetInFormat("smi"))
    conv.ReadString(&obAlias, pos->second);
  
  
  //Find index of atom to which XxAtom is attached, and the eventual index of first atom in fragment
  OBBondIterator bi;
  unsigned mainAttachIdx = (Xxatom->BeginNbrAtom(bi))->GetIdx();
  unsigned    aliasindex = mol.NumAtoms()+1;

 //Combine with main molecule 
  mol += obAlias;

  //remove trailing _ which += operator has added
  string title(mol.GetTitle());
  if(title[title.size()-1]=='_')
  {
    title.erase(title.size()-1);  
    mol.SetTitle(title);
  }

  //Remove bond to Xx atom
  mol.DeleteBond(*Xxatom->BeginBonds());

  //Connect to SMILES fragment
  mol.AddBond(mainAttachIdx, aliasindex, 1);

  //delete original Xx atom
  mol.DeleteAtom(Xxatom); 

  return true;
}

bool AliasData::LoadFile(std::map<std::string, std::string>& table)
{
  ifstream ifs;
  if (OpenDatafile(ifs, "superatom.txt").length() == 0)
  {
    obErrorLog.ThrowError(__FUNCTION__, "Cannot open superatom.txt", obError);
    return false;
  }
  string ln;
  while(getline(ifs, ln))
  {
    if (ln[0]=='#' || ln.empty())
      continue;
    std::vector<string> vec;
    if(tokenize(vec, ln) && vec.size()>=2)
      table[ vec[0] ] = vec[1];
  }
  return true;
}
}//namespace

//! \file alias.cpp
//! \brief OBGenericData class to for atom alias data (e.g., in 2D drawing programs for "COOH")
