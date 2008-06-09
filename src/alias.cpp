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
#include <openbabel/alias.h>

namespace OpenBabel
{

  bool AliasData::Expand(OBMol& mol, const unsigned int atomindex)
  {
    /*
    Interprets the alias text and adds atom as appropriate to mol.
    Tries the following in turn until one is sucessful:
    1) Looks up alias in fragments.txt (not implemented yet)
    2) Parse as simple formula
    Returns false if none are successful.
    */

    //Look up alias name here, add atoms as appropriate and return true if successful
    //else...

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
      return true;
    if(!isalpha(*txt)) //first char is the element that replaces atomindex
      return false;
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

}//namespace

//! \file alias.cpp
//! \brief OBGenericData class to for atom alias data (e.g., in 2D drawing programs for "COOH")
