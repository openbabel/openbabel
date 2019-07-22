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
#include <cstdlib>
#include <sstream>
#include <string>
#include <openbabel/alias.h>
#include <openbabel/obconversion.h>
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/builder.h>
#include <openbabel/obiter.h>
#include <openbabel/parsmart.h>
#include <openbabel/mcdlutil.h>
#include <openbabel/shared_ptr.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#define MARK_UNUSED(x) (void)(x)

using namespace std;
namespace OpenBabel
{
  std::string AliasData::GetAlias(bool rightAligned)const
  {
    if(rightAligned)
    {
      if(!_right_form.empty())
        return _right_form;
      if(table().find(_alias)!=table().end())
        return table().find(_alias)->second.right_form;
    }
    return _alias;
  }

  bool AliasData::Expand(OBMol& mol, const unsigned int atomindex)
  {
    /*
    Interprets the alias text and adds atom(s) as appropriate to mol.
    Tries the following in turn until one is sucessful:
    1) If starts with number treat as isotope+element e.g. 2H
    2) Looks up alias in superatom.txt e.g. COOH Pr
    3) If of the form Rn stored as a * atom with Atom Class data
    Returns false if none are successful.
    */

    //parse as isotopic atom
    if(isdigit(_alias[0]))
    {
      std::stringstream ss(_alias);
      int iso;
      std::string el;
      ss >> iso >>el;
      unsigned int elemno = OBElements::GetAtomicNum(el.c_str());
      if(elemno > 0)
      {
        OBAtom* pAtom = mol.GetAtom(atomindex);
        if(!pAtom)
          return false;
        pAtom->SetIsotope(iso);
        pAtom->SetAtomicNum(elemno);
        return true;
      }
    }

    if(FromNameLookup(mol, atomindex))
      return true;

    // Rn is stored as an atom with 0 atomic number and atomclass = n
    // R', R'' etc. are treated as R1, R2
    // Note that if the name contains anything after the number it is ignored.
    if(_alias[0]=='R' && (_alias[1]=='\'' || isdigit(_alias[1])))
    {
      unsigned int n = 1;
      if(_alias[1]=='\'')
        while(n<_alias.size()-1 && _alias[n]==_alias[n+1]) n++;
      else
        n = atoi(_alias.c_str()+1);

      OBPairInteger *atomclass = new OBPairInteger();
      atomclass->SetAttribute("Atom Class");
      atomclass->SetValue(n);
      mol.GetAtom(atomindex)->SetData(atomclass);

      if(atomindex <= mol.NumAtoms()) //needed for Rn aliases in mdlformat
        mol.GetAtom(atomindex)->SetAtomicNum(0);

      _right_form = _alias;
      return true;
    }

    obErrorLog.ThrowError(__FUNCTION__, "Alias " + _alias +
      " was not chemically interpreted\n", obWarning, onceOnly);
    return false;
  }

bool AliasData::FromNameLookup(OBMol& mol, const unsigned int atomindex)
{
  /*Converts an alias name (like COOH) to real chemistry:
    looks up in a table loaded from superatom.txt;
    converts the SMILES of the fragment and adds it to the molecule.
    If the molecule already has atom coordinates, generates coordinates
    for the new atoms, using builder for 3D and MCDL for 2D.
  */

  OBAtom* XxAtom = mol.GetAtom(atomindex);
/*  if(XxAtom->GetExplicitDegree()>1)
  {
    obErrorLog.ThrowError(__FUNCTION__, _alias + " is multivalent, which is currently not supported.", obWarning);
    return false;
  }
*/
  SuperAtomTable::iterator pos = table().find(_alias);
  if(pos==table().end())
    return false;

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
    conv.ReadString(&obFrag, '*' + pos->second.smiles);//Add dummy atom to SMILES
    _right_form = pos->second.right_form;
    _color      = pos->second.color;
  }
  obFrag.SetDimension(dimension);//will be same as parent

  //Find index of *first* atom to which XxAtom is attached (could be NULL)
  OBBondIterator bi;
  OBAtom* firstAttachAtom = XxAtom->BeginNbrAtom(bi);
  unsigned mainAttachIdx = firstAttachAtom ? firstAttachAtom->GetIdx() : 0;
  unsigned int firstAttachFlags = 0;
  unsigned int firstAttachOrder = 1;
  if (firstAttachAtom) {
    firstAttachFlags = mol.GetBond(XxAtom, firstAttachAtom)->GetFlags();
    firstAttachOrder = mol.GetBond(XxAtom, firstAttachAtom)->GetBondOrder();
  }
  //++Make list of other attachments* of XxAtom
  // (Added later so that the existing bonding of the XXAtom are retained)
  vector<pair<OBAtom*, unsigned> > otherAttachments;
  OBAtom* pAttach;
  while(firstAttachAtom && (pAttach = XxAtom->NextNbrAtom(bi)) ) // extra parentheses to minimize warnings
    otherAttachments.push_back(make_pair(pAttach, (*bi)->GetBondOrder()));

  //Copy coords of XxAtom to the first real atom in the fragment
  //so that the connecting bond is well defined for 2D case
  obFrag.GetAtom(2)->SetVector( XxAtom->GetVector());

  //delete original Xx atom
  mol.DeleteAtom(XxAtom, false);//delay deletion of the OBAtom object because this is attached to it
  //Correct indices for the deletion
  if(atomindex<mainAttachIdx)
    --mainAttachIdx;

  //Find the eventual index of first atom in fragment
  unsigned newFragIdx = mol.NumAtoms()+1;

  //Give the fragment appropriate coordinates
  if(dimension==3)
  {
    OBBuilder builder;

    builder.Build(obFrag);
    obFrag.DeleteAtom(obFrag.GetAtom(1));//remove dummy atom
    mol += obFrag; //Combine with main molecule
    if(mainAttachIdx) {
      builder.Connect(mol, mainAttachIdx, newFragIdx,XxAtom->GetVector(),firstAttachOrder);
    }
  }
  else // 0D, 2D
  {
    obFrag.DeleteAtom(obFrag.GetAtom(1));//remove dummy atom
    mol += obFrag; //Combine with main molecule and connect
    if(mainAttachIdx)
      mol.AddBond(mainAttachIdx, newFragIdx, 1, firstAttachFlags);
  }

  if(dimension==2)//Use MCDL
    groupRedraw(&mol, mol.NumBonds()-1, newFragIdx, true);

  //++Add bonds from list to newFragIdx
  while(!otherAttachments.empty())
  {
    mol.AddBond(otherAttachments.back().first->GetIdx(), newFragIdx, otherAttachments.back().second);
    otherAttachments.pop_back();
  }

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

        obsharedptr<OBSmartsPattern> psp(new OBSmartsPattern);
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
  //The atom that carries the AliasData object remains as an Xx atom with no charge;
  //the others are deleted. All the attached hydrogens are also deleted.
  for(unsigned i=0;i<_expandedatoms.size();++i)
  {
    OBAtom* at = mol.GetAtomById(_expandedatoms[i]);
    if(!at)
      continue;
    mol.DeleteHydrogens(at);
    if(at->HasData(AliasDataType))
    {
      at->SetAtomicNum(0);
      at->SetFormalCharge(0);
      at->SetSpinMultiplicity(0);
    }
    else
      mol.DeleteAtom(at);
  }
  _expandedatoms.clear();
}

void AliasData::RevertToAliasForm(OBMol& mol)
{
  //Deleting atoms invalidates the iterator, so start again
  //and continue until all no unexpanded aliases are found in molecule.
  bool acted = false;
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
  // Mark variables as unused to avoid warnings
  MARK_UNUSED(OptionText);
  MARK_UNUSED(pmap);

  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;
  return AliasData::AddAliases(pmol);
}
#endif // HAVE_SHARED_POINTER

}//namespace

//! \file alias.cpp
//! \brief OBGenericData class to for atom alias data (e.g., in 2D drawing programs for "COOH")
