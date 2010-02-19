/**********************************************************************
alias.cpp - implementation of an OBGenericData class to hold alias information on atoms
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
#include <openbabel/op.h>
#include <openbabel/builder.h>
#include <openbabel/parsmart.h>

using namespace std;
namespace OpenBabel
{
  std::string AliasData::GetAlias(bool rightAligned)const
  {
    if(rightAligned)
      return table().find(_alias)->second.right_form;//.first;
    else
      return _alias;
  }

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
  /*Converts an alias name (like COOH) to real chemistry:
    looks up in a table loaded from superatom.txt;
    converts the SMILES of the fragment and adds it to the molecule.
    If the molecule already has atom coordinates, uses builder to generate
    coordinates for the new atoms.
  */

  OBAtom* XxAtom = mol.GetAtom(atomindex);
  if(XxAtom->GetValence()>1)
  {
    obErrorLog.ThrowError(__FUNCTION__, _alias + " is multivalent, which is currently not supported.", obError);
    return false;
  }

  SuperAtomTable::iterator pos = table().find(_alias);
  if(pos==table().end())
  {
    obErrorLog.ThrowError(__FUNCTION__, "Alias " + _alias + " was not recognized.\n Output may not be correct.", obError, onceOnly);
    return false;
  }

  int dimension=0;
  if(mol.Has3D())
    dimension=3;
  else if(mol.Has2D())
    dimension=2;
  mol.SetDimension(dimension);

  //Convert SMILES of alias
  OBConversion conv;
  OBMol obFrag;
  obFrag.SetIsPatternStructure();
  if(conv.SetInFormat("smi"))
  {
    conv.ReadString(&obFrag, pos->second.smiles);//.second);
    _right_form = pos->second.right_form;
    _color      = pos->second.color;
  }

  //Find index of atom to which XxAtom is attached, and the eventual index of first atom in fragment
  OBBondIterator bi;
  unsigned mainAttachIdx = (XxAtom->BeginNbrAtom(bi))->GetIdx();
  unsigned    newFragIdx = mol.NumAtoms()+1;
  
  //obFrag.SetDimension(dimension);
  OBBuilder builder;
  //Give the fragment appropriate coordinates
  if(dimension!=0)
    builder.Build(obFrag);

 //Combine with main molecule 
  mol += obFrag;

  //delete original Xx atom
  mol.DeleteAtom(XxAtom, false);//delay deletion of the OBAtom object because this is attached to it
  //Correct indices for the deletion
  --newFragIdx;
  if(atomindex<mainAttachIdx)
    --mainAttachIdx;

  //Attach the fragment to the main molecule
  if(dimension==0)
    mol.AddBond(mainAttachIdx, newFragIdx, 1);
  else
    builder.Connect(mol, mainAttachIdx, newFragIdx);
  
  //Store the ids of the atoms which replace the alias (the last atoms in the combined molecule).
  //The ids do not change when other atoms are deleted.
  for(unsigned i=obFrag.NumAtoms();i;--i)
    _expandedatoms.push_back(mol.GetAtom(mol.NumAtoms()-i +1)->GetId());

  //Make a copy of this AliasData object (currently attached to XxAtom)
  //and attach it to the first atom of the fragment.
  mol.GetAtom(newFragIdx)->CloneData(this);

  delete(XxAtom);
  return true;
}

bool AliasData::LoadFile(SuperAtomTable& table)
{
  //In table: key=alias left-form; value=pair<alias right-form, SMILES>
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
    if(tokenize(vec, ln) && vec.size()>=3)
    {
      //table[ vec[0] ] = make_pair(vec[1], vec[2]);
      AliasItem item;
      item.right_form = vec[1];
      item.smiles     = vec[2];
      item.color      = vec.size()>=4 ? vec[3] : "";
      table[ vec[0] ] = item;
    }
  }
  return true;
}
#ifdef HAVE_SHARED_POINTER
bool AliasData::LoadFile(SmartsTable& smtable)
{
  //Re-parse the datafile. Seems simpler than trying to extract from the map.
  ifstream ifs;
  if (OpenDatafile(ifs, "superatom.txt").length() == 0)
  {
    obErrorLog.ThrowError(__FUNCTION__, "Cannot open superatom.txt", obError);
    return false;
  }
  string ln;
  while(getline(ifs, ln))
  {
    if ((ln[0]=='#' && ln[1]!='#') || ln.empty())
      continue;
    if(ln[0]=='#') //stop reading at line starting with ##
      break;
    std::vector<string> vec;
    if(tokenize(vec, ln) && vec.size()>=3)
    {
      //Convert SMILES with implicit H to SMARTS with explicit H.
      //Converting into and out of OBMol is a bit heavy, but saves 
      //worrying about edge cases in a string parse.
      stringstream ss('*'+vec[2]),// '*' added to SMILES because the superatom has to be attached
                   ssmarts; 
      OBConversion conv(&ss, &ssmarts);
      conv.AddOption("h",OBConversion::GENOPTIONS);//add explicit Hs...
      conv.AddOption("h");//...and output them to ensure the superatom itself is not substituted
      if(conv.SetInAndOutFormats("smi","smi"))
        conv.Convert();
      if(!ssmarts.str().empty())
      {
        //OBSmartsPattern objects are not copyable without complications,
        //so reference semantics used.

        shared_ptr<OBSmartsPattern> psp(new OBSmartsPattern);
        psp->Init(ssmarts.str());
        smtable.push_back(make_pair(vec[0], psp));
      }
    }
  }
  return true;
}
#endif

void AliasData::AddExpandedAtom(int id) { _expandedatoms.push_back(id); };

void AliasData::DeleteExpandedAtoms(OBMol& mol)
{
  //The atom that carries the AliasData object remains as an Xx atom;
  //the others are deleted.
  for(unsigned i=0;i<_expandedatoms.size();++i)
  {
    OBAtom* at = mol.GetAtomById(_expandedatoms[i]);
    if(!at)
      continue;
    if(at->HasData(AliasDataType))
      at->SetAtomicNum(0);
    else
      mol.DeleteAtom(at);
  }
  _expandedatoms.clear();
}

void AliasData::RevertToAliasForm(OBMol& mol)
{
  //Deleting atoms invalidates the iterator, so start again
  //and continue until all no unexpanded aliases are found in molecule. 
  bool acted;
  do
  {
    FOR_ATOMS_OF_MOL(a, mol)
    {
      acted=false;
      AliasData* ad = NULL; 
      if((ad = (static_cast<AliasData*>(a->GetData(AliasDataType)))) && ad->IsExpanded())
      {
        ad->DeleteExpandedAtoms(mol);
        acted = true;
        break;
      }
    }
  }while(acted);
}

#ifdef HAVE_SHARED_POINTER
bool AliasData::AddAliases(OBMol* pmol)
{
  static SmartsTable smtable;
  if(smtable.empty())
    LoadFile(smtable);
  set<int> AllExAtoms;
  SmartsTable::iterator iter;
  for(iter=smtable.begin();iter!=smtable.end();++iter)
  {
    if((*iter).second->Match(*pmol))
    {
      vector<std::vector<int> > mlist = (*iter).second->GetUMapList();      
      for(unsigned imatch=0;imatch<mlist.size();++imatch) //each match
      {
        AliasData* ad  = new AliasData;
        ad->SetAlias((*iter).first);
        //iatom==0 is the * that was added to the front of the SMILES, so start at 1
        for(unsigned iatom=1; iatom<mlist[imatch].size();++iatom)//each atom in match
        {
          int idx = mlist[imatch][iatom];
          if(AllExAtoms.count(idx))
          {
            //atom already appears in an alias so abandon this (smaller) alias
            delete ad;
            ad = NULL;
            break;
          }
          else
          {
            AllExAtoms.insert(idx);        
            int id  =(pmol->GetAtom(idx))->GetId();
            ad->AddExpandedAtom(id);
          }
        }
        if(ad)
          pmol->GetAtom(mlist[imatch][1])->SetData(ad);//attach alias to first expanded atom
      }
    }
  }
  return true;
}

//OpGenAlias is a wrapper for AddAliases. Use like:
//babel infile.xxx outfile.yyy --genalias
class OpGenAlias : public OBOp
{
public:
  OpGenAlias(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return "Generate aliases as an alternative representation."; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*);
};

/////////////////////////////////////////////////////////////////
OpGenAlias theOpGenAlias("genalias"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpGenAlias::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion*)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;
  return AliasData::AddAliases(pmol);
}
#endif // HAVE_SHARED_POINTER

}//namespace

//! \file alias.cpp
//! \brief OBGenericData class to for atom alias data (e.g., in 2D drawing programs for "COOH")
