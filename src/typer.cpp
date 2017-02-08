/**********************************************************************
typer.cpp - Open Babel atom typer.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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

#include <openbabel/mol.h>
#include <openbabel/typer.h>
#include <openbabel/obiter.h>

// private data header with default parameters
#include "atomtyp.h"

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

using namespace std;

namespace OpenBabel
{

  OBAromaticTyper  aromtyper;
  OBAtomTyper      atomtyper;

  /*! \class OBAtomTyper typer.h <openbabel/typer.h>
    \brief Assigns atom types, hybridization, implicit valence and formal charges

    The OBAtomTyper class is designed to read in a list of atom typing
    rules and apply them to molecules. The code that performs atom
    typing is not usually used directly as atom typing, hybridization
    assignment, implicit valence assignment and charge are all done
    automatically when their corresponding values are requested of
    atoms:
    \code
    atom->GetType();
    atom->GetFormalCharge();
    atom->GetHyb();
    \endcode
  */
  OBAtomTyper::OBAtomTyper()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "atomtyp.txt";
    _subdir = "data";
    _dataptr = AtomTypeData;
  }

  void OBAtomTyper::ParseLine(const char *buffer)
  {
    vector<string> vs;
    OBSmartsPattern *sp;

    if (EQn(buffer,"INTHYB",6))
      {
        tokenize(vs,buffer);
        if (vs.size() < 3)
          {
            obErrorLog.ThrowError(__FUNCTION__, " Could not parse INTHYB line in atom type table from atomtyp.txt", obInfo);
            return;
          }

        sp = new OBSmartsPattern;
        if (sp->Init(vs[1]))
          _vinthyb.push_back(pair<OBSmartsPattern*,int> (sp,atoi((char*)vs[2].c_str())));
        else
          {
            delete sp;
            sp = NULL;
            obErrorLog.ThrowError(__FUNCTION__, " Could not parse INTHYB line in atom type table from atomtyp.txt", obInfo);
            return;
          }
      }
    else if (EQn(buffer,"IMPVAL",6))
      {
        tokenize(vs,buffer);
        if (vs.size() < 3)
          {
            obErrorLog.ThrowError(__FUNCTION__, " Could not parse IMPVAL line in atom type table from atomtyp.txt", obInfo);
            return;
          }

        sp = new OBSmartsPattern;
        if (sp->Init(vs[1]))
          _vimpval.push_back(pair<OBSmartsPattern*,int> (sp,atoi((char*)vs[2].c_str())));
        else
          {
            obErrorLog.ThrowError(__FUNCTION__, " Could not parse IMPVAL line in atom type table from atomtyp.txt", obInfo);
            delete sp;
            sp = NULL;
            return;
          }
      }
    else if (EQn(buffer,"EXTTYP",6))
      {
        tokenize(vs,buffer);
        if (vs.size() < 3)
          {
            obErrorLog.ThrowError(__FUNCTION__, " Could not parse EXTTYP line in atom type table from atomtyp.txt", obInfo);
            return;
          }
        sp = new OBSmartsPattern;
        if (sp->Init(vs[1]))
          _vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
        else
          {
            delete sp;
            sp = NULL;
            obErrorLog.ThrowError(__FUNCTION__, " Could not parse EXTTYP line in atom type table from atomtyp.txt", obInfo);
            return;
          }
      }
  }

  OBAtomTyper::~OBAtomTyper()
  {
    vector<pair<OBSmartsPattern*,int> >::iterator i;
    for (i = _vinthyb.begin();i != _vinthyb.end();++i)
      {
        delete i->first;
        i->first = NULL;
      }
    for (i = _vimpval.begin();i != _vimpval.end();++i)
      {
        delete i->first;
        i->first = NULL;
      }

    vector<pair<OBSmartsPattern*,string> >::iterator j;
    for (j = _vexttyp.begin();j != _vexttyp.end();++j)
      {
        delete j->first;
        j->first = NULL;
      }

  }

  void OBAtomTyper::AssignTypes(OBMol &mol)
  {
    if (!_init)
      Init();

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::AssignTypes", obAuditMsg);

    mol.SetAtomTypesPerceived();

    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;

    for (i = _vexttyp.begin();i != _vexttyp.end();++i)
      if (i->first->Match(mol))
        {
          _mlist = i->first->GetMapList();
          for (j = _mlist.begin();j != _mlist.end();++j)
            mol.GetAtom((*j)[0])->SetType(i->second);
        }

    // Special cases
    vector<OBAtom*>::iterator a;
    OBAtom* atom;
    for (atom = mol.BeginAtom(a); atom; atom = mol.NextAtom(a)) {
      // guanidinium. Fixes PR#1800964
      if (strncasecmp(atom->GetType(),"C2", 2) == 0) {
        int guanidineN = 0;
        OBAtom *nbr;
        vector<OBBond*>::iterator k;
        for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k)) {
          if (strncasecmp(nbr->GetType(),"Npl", 3) == 0 ||
              strncasecmp(nbr->GetType(),"N2", 2) == 0 ||
              strncasecmp(nbr->GetType(),"Ng+", 3) == 0)
            ++guanidineN;
        }
        if (guanidineN == 3)
          atom->SetType("C+");

      } // end C2 carbon for guanidinium

    } // end special cases
  }

  void OBAtomTyper::AssignHyb(OBMol &mol)
  {
    if (!_init)
      Init();

    aromtyper.AssignAromaticFlags(mol);

    mol.SetHybridizationPerceived();
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::AssignHybridization", obAuditMsg);

    OBAtom *atom;
    vector<OBAtom*>::iterator k;
    for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k))
      atom->SetHyb(0);

    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,int> >::iterator i;

    for (i = _vinthyb.begin();i != _vinthyb.end();++i)
      if (i->first->Match(mol))
        {
          _mlist = i->first->GetMapList();
          for (j = _mlist.begin();j != _mlist.end();++j)
            mol.GetAtom((*j)[0])->SetHyb(i->second);
        }
  }

  void OBAtomTyper::AssignImplicitValence(OBMol &mol, bool CanBeLessThanActual)
  {
    // FF Make sure that valence has not been perceived
    if(mol.HasImplicitValencePerceived())
      return;

    if (!_init)
      Init();

    mol.SetImplicitValencePerceived();
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::AssignImplicitValence", obAuditMsg);

    // FF Ensure that the aromatic typer will not be called
    int oldflags = mol.GetFlags(); // save the current state flags
    mol.SetAromaticPerceived();    // and set the aromatic perceived flag on

    OBAtom *atom;
    vector<OBAtom*>::iterator k;
    for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k))
      atom->SetImplicitValence(atom->GetValence());

    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,int> >::iterator i;

    for (i = _vimpval.begin();i != _vimpval.end();++i)
      if (i->first->Match(mol))
        {
          _mlist = i->first->GetMapList();
          for (j = _mlist.begin();j != _mlist.end();++j)
            mol.GetAtom((*j)[0])->SetImplicitValence(i->second);
        }

    if (!mol.HasAromaticCorrected())
      CorrectAromaticNitrogens(mol);

    if(!CanBeLessThanActual)
      for (atom = mol.BeginAtom(k);atom;atom = mol.NextAtom(k))
        {
          if (atom->GetImplicitValence() < atom->GetValence())
            atom->SetImplicitValence(atom->GetValence());
        }

    // FF Come back to the initial flags
    mol.SetFlags(oldflags);

    return;
  }

  //! Currently sets OBMol::SetAromaticCorrected and returns.
  //! \deprecated Currently unused for anything significant.
  void OBAtomTyper::CorrectAromaticNitrogens(OBMol &mol)
  {
    if (!_init)
      Init();

    if (mol.HasAromaticCorrected())
      return;
    mol.SetAromaticCorrected();

    FOR_ATOMS_OF_MOL(atom, mol) {
      if (atom->IsNitrogen() && atom->IsAromatic()) {
        atom->SetHyb(2);
        atom->SetType("Nar");
        if (atom->HasDoubleBond()) {
          atom->SetImplicitValence(2 + atom->GetFormalCharge());
        } else {
          if (atom->GetImplicitValence() == 2)
            atom->SetImplicitValence(3 + atom->GetFormalCharge());
        }
      }
    }

    return;
  }

  /*! \class OBRingTyper typer.h <openbabel/typer.h>
    \brief Assigns ring types

    The OBRingTyper class is designed to read in a list of ring typing
    rules and apply them to molecules. The code that performs ring
    typing is not usually used directly as ring typing is done
    automatically when the ring type is requested of rings:
    \code
    vector<OBRing*>::iterator i;
    vector<OBRing*> rlist = mol.GetSSSR();

    for (i = rlist.begin();i != rlist.end();++i)
    cout << "ring type = " << (*i)->GetType() << endl;
    \endcode
  */
  OBRingTyper::OBRingTyper()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "ringtyp.txt";
    _subdir = "data";
    //_dataptr = RingTypeData;
  }

  void OBRingTyper::ParseLine(const char *buffer)
  {
    vector<string> vs;
    OBSmartsPattern *sp;

    if (EQn(buffer,"RINGTYP",7)) {
      tokenize(vs,buffer);
      if (vs.size() < 3) {
        obErrorLog.ThrowError(__FUNCTION__, " Could not parse RING line in ring type table from ringtyp.txt", obInfo);
        return;
      }
      sp = new OBSmartsPattern;
      if (sp->Init(vs[2]))
        _ringtyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[1]));
      else {
        delete sp;
        sp = NULL;
        obErrorLog.ThrowError(__FUNCTION__, " Could not parse RING line in ring type table from ringtyp.txt", obInfo);
        return;
      }
    }
  }

  OBRingTyper::~OBRingTyper()
  {
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    for (i = _ringtyp.begin();i != _ringtyp.end();++i) {
      delete i->first;
      i->first = NULL;
    }
  }

  void OBRingTyper::AssignTypes(OBMol &mol)
  {
    if (!_init)
      Init();

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OBRing::AssignTypes", obAuditMsg);

    mol.SetRingTypesPerceived();

    vector<vector<int> >::iterator j2;
    vector<pair<OBSmartsPattern*,string> >::iterator i2;

    vector<OBRing*>::iterator i;
    vector<int>::iterator j;
    vector<OBRing*> rlist = mol.GetSSSR();

    unsigned int member_count;
    for (i2 = _ringtyp.begin();i2 != _ringtyp.end();++i2) { // for each ring type
      if (i2->first->Match(mol)) {
        _mlist = i2->first->GetMapList();
        for (j2 = _mlist.begin();j2 != _mlist.end();++j2) { // for each found match

          for (i = rlist.begin();i != rlist.end();++i) { // for each ring
            member_count = 0;

            for(j = j2->begin(); j != j2->end(); ++j) { // for each atom in the match
              if ((*i)->IsMember(mol.GetAtom(*j)))
                member_count++;
            }

            if ((*i)->Size() == member_count)
              (*i)->SetType(i2->second);
          }
        }
      }
    }
  }


  /*! \class OBAromaticTyper typer.h <openbabel/typer.h>
    \brief Assigns aromatic typing to atoms and bonds

    The OBAromaticTyper class applies a set of
    aromatic perception rules to molecules. The code
    that performs typing is not usually used directly -- it is usually
    done automatically when their corresponding values are requested of atoms
    or bonds.
    \code
    atom->IsAromatic();
    bond->IsAromatic();
    bond->IsDouble(); // needs to check aromaticity and define Kekule structures
    \endcode
  */
  OBAromaticTyper::OBAromaticTyper()
  {
    _init = true;
    // The following values are no longer used by this class
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "aromatic.txt";
    _subdir = "data";
    _dataptr = (const char*)0; // AromaticData;
  }

  void OBAromaticTyper::ParseLine(const char *buffer)
  {
    // This function is no longer used by this class
  }

  OBAromaticTyper::~OBAromaticTyper()
  {
  }

  // Start of predicates for AssignOBAromaticityModel
  static bool HasRingDblBond(OBAtom *atm)
  {
    FOR_BONDS_OF_ATOM(bond, atm) {
      if (bond->GetBondOrder() == 2 && bond->IsInRing())
        return true;
    }
    return false;
  }
  static bool HasAcyclicDblBondToOxygen(OBAtom *atm)
  {
    FOR_BONDS_OF_ATOM(bond, atm) {
      if (bond->GetBondOrder() == 2 && !bond->IsInRing() && bond->GetNbrAtom(atm)->GetAtomicNum() == 8)
        return true;
    }
    return false;
  }
  static bool HasTwoRingOxygens(OBAtom *atm)
  {
    unsigned int totRingOxygens = 0;
    FOR_BONDS_OF_ATOM(bond, atm) {
      if (bond->GetBondOrder()==1) {
        OBAtom *nbr = bond->GetNbrAtom(atm);
        if (nbr->IsInRing() && nbr->GetAtomicNum() == 8) {
          if (totRingOxygens == 1)
            return true;
          totRingOxygens++;
        }
      }
    }
    return false;
  }
  static unsigned int ThreeNitrogenTestsInOne(OBAtom *atm)
  {
    // Return a bitwise combination of 
    //   1 -> HasAcyclicDblBondToOxygen
    //   2 -> HasTwoSingleRingBonds
    //   4 -> HasRingDblBond

    unsigned int ans = 0;
    unsigned int totSingleRingBonds = 0;

    FOR_BONDS_OF_ATOM(bond, atm) {
      OBAtom *nbr = bond->GetNbrAtom(atm);
      unsigned int nbr_elem = nbr->GetAtomicNum();
      if (nbr_elem == 1) continue;
      unsigned int bo = bond->GetBondOrder();
      if (bond->IsInRing()) {
        switch (bo) {
        case 2:
          ans |= 4; break; // HasRingDblBond
        case 1:
          totSingleRingBonds++; break;
        }
      } else { // acyclic
        if (bo == 2 && nbr_elem == 8)
          ans |= 1; // HasAcyclicDblBondToOxygen
      }
    }
    
    if (totSingleRingBonds >= 2)
      ans |= 2;

    return ans;
  }
  static unsigned int ThreeCarbonTestsInOne(OBAtom *atm)
  {
    // Return a bitwise combination of 
    //   1 -> HasAcyclicDblBond
    //   2 -> HasAcyclicBondToOxygen
    //   or just return 
    //   99 -> HasAcyclicDblBondToNonCarbon

    unsigned int ans = 0;

    FOR_BONDS_OF_ATOM(bond, atm) {
      if (bond->IsInRing()) continue;
      unsigned int nbr_elem = bond->GetNbrAtom(atm)->GetAtomicNum();
      unsigned int bo = bond->GetBondOrder();
      switch (bo) {
      case 2:
        if (nbr_elem != 6)
          return 99; // Can return straight-away as has priority
        ans |= 1;
        // no break - fall through
      case 1:
        if (nbr_elem == 8)
          ans |= 2;
      }
    }

    return ans;
  }
  static bool IsSP2(OBAtom * atm, unsigned int elem)
  {
    // phosphorus and sulfur may be initially typed as sp3
    switch (elem) {
    case 6:
    case 7:
    case 8:
      if (atm->GetHyb() != 2)
        return false;
    }
    return true;
  }
  // End of predicates for AssignOBAromaticityModel

  static bool AssignOBAromaticityModel(OBAtom *atm, int &min, int &max)
  {
    // The Open Babel aromaticity model
    //
    // Return minimum and maximum pi-electrons contributed to an aromatic system
    // The return value indicates a potentially aromatic atom (i.e. any patterns matched)
    //
    // The original code used SMARTS patterns organised in increasing order of
    // prioirity (i.e. later matches overrode earlier ones). These SMARTS patterns
    // are now implemented in the code below, but are included in the comments
    // for reference (Case 1->22).

    unsigned int elem = atm->GetAtomicNum();
    if (!atm->IsInRing() || atm->GetImplicitValence() > 3 || !IsSP2(atm, elem)) {
      min = 0; max = 0;
      return false;
    }
    
    unsigned int deg = atm->GetHvyValence();
    unsigned int allatomdeg = atm->GetImplicitValence();
    int chg = atm->GetFormalCharge();
    unsigned int threeNTests, threeCTests;
    switch (elem) {
    case 6: // carbon
      switch (chg) {
      case 0:
        if (HasRingDblBond(atm)) {
          min = 1; max = 1; return true; // Case 4 [#6r+0]=@*
        }
        switch (deg) {
        case 2:
          min = 1; max = 1; return true; // Case 1 [#6rD2+0]
        case 3:
          threeCTests = ThreeCarbonTestsInOne(atm);
          switch (threeCTests) {
          case 99: // HasAcyclicDblBondToNonCarbon (priority 1)
            // external double bonds to hetero atoms contribute no electrons to the
            // aromatic systems -- quinoid systems are non - aromatic, e.g. 1, 4 - benzoquinone
            min = 0; max = 0; return true; // Case 6 [#6rD3+0]=!@[!#6]
          case 1: // HasAcyclicDblBond (priority 2)
          case 3:
            min = 1; max = 1; return true; // Case 5 [#6rD3+0]=!@*
          case 2: // HasAcyclicBondToOxygen (priority 3)
            min = 0; max = 1; return true; // Case 2 [#6rD3+0]~!@[#8] (exo ketone or alcohol -- don't know which)
          }
        }
        break;
      case 1:
        if (allatomdeg == 2 || allatomdeg == 3) {
          min = 1; max = 1; return true; // Case 3 [#6rD2+,#6rD3+]
        }
        break;
      case -1:
        if (allatomdeg == 3) {
          min = 2; max = 2; return true; // Case 7 [#6rD3-]
        }
      }
      break;
    case 7: // nitrogen
      switch (chg) {
      case 0:
        threeNTests = ThreeNitrogenTestsInOne(atm);
        if (deg == 3 && (threeNTests & 1)) { // HasAcyclicDblBondToOxygen(atm))
          min = 1; max = 1; return true; // Case 13 [#7rD3+0]=O
        }
        if (threeNTests & 2) { // HasTwoSingleRingBonds(atm))
          min = 1; max = 2; return true; // Case 10 [#7r+0](-@*)-@*
        }
        switch (deg) {
        case 2:
          if (threeNTests & 4) { // HasRingDblBond(atm))
            min = 1; max = 1; return true; // Case 11 [#7rD2+0]=@*
          }
          else {
            min = 1; max = 2; return true; // Case 8 [#7rD2+0]
          }
          break;
        case 3:
          min = 1; max = 2; return true; // Case 9 [#7rD3+0]
        }
        break;
      case -1:
        if (deg == 2) {
          min = 2; max = 2; return true; // Case 14 [#7rD2-]
        }
        break;
      case 1:
        if ((deg==2 && (atm->ImplicitHydrogenCount() + atm->BOSum())==4) || deg == 3) { // tried deg 2 or 3 but fails on emolecules 6884346
          min = 1; max = 1; return true; // Case 12 [#7rD3+]
        }
      }
      break;
    case 8: // oxygen
      switch (chg) {
      case 0:
        min = 2; max = 2; return true; // Case 15 [#8r+0]
      case 1:
        min = 1; max = 1; return true; // Case 16 [#8r+]
      }
      break;
    case 15: // phosphorus
      // Accounts Chem Res 1978 11 p. 153
      // phosphole, phosphabenzene(not v.aromatic)
      if (deg == 3) {
        if (HasTwoRingOxygens(atm)) {
          min = 0; max = 0; return true; // Case 21 [#15rD3]([#8r])[#8r]
        }
        if (chg == 0) {
          min = 2; max = 2; return true; // Case 20 [#15rD3+0]
        }
      }
      break;
    case 16: // sulfur
      switch (chg) {
      case 0:
        switch (deg) {
        case 2:
          min = 2; max = 2; return true; // Case 17 [#16rD2+0]
        case 3:
          if (HasAcyclicDblBondToOxygen(atm)) {
            min = 2; max = 2; return true; // Case 19 [#16rD3+0]=!@O
          }
        }
        break;
      case 1:
        if (allatomdeg == 2) {
          min = 1; max = 1; return true; // Case 18 [#16rD2+]
        }
      }
      break;
    case 34: // selenium (selenophene)
      if (chg == 0 && deg == 2) {
        min = 2; max = 2; return true; // Case 22 [#34rD2+0]
      }
    }

    // Considered and excluded:
    // arsabenzene, etc. (*really* not v.aromatic)
    // [#33rD3+0]		2	2
    // tellurophene, etc. (*really* not v.aromatic)
    // [#52rD2+0]		2	2
    // stilbabenzene, etc. (very little aromatic character)
    // [#51rD3+0]		2	2

    min = 0; max = 0; // nothing matched
    return false;
  }

  void OBAromaticTyper::AssignAromaticFlags(OBMol &mol)
  {
    if (!_init)
      Init();

    if (mol.HasAromaticPerceived())
      return;
    mol.SetAromaticPerceived();
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::AssignAromaticFlags", obAuditMsg);

    _vpa.clear();
    _vpa.resize(mol.NumAtoms()+1);
    _velec.clear();
    _velec.resize(mol.NumAtoms()+1);
    _root.clear();
    _root.resize(mol.NumAtoms()+1);

    OBBond *bond;
    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;

    //unset all aromatic flags
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      atom->UnsetAromatic();
    for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
      bond->UnsetAromatic();

    // New code using lookups instead of SMARTS patterns
    FOR_ATOMS_OF_MOL(atom, mol) {
      unsigned int idx = atom->GetIdx();
      _vpa[idx] = AssignOBAromaticityModel(&(*atom), _velec[idx].first, _velec[idx].second);
      }

    //propagate potentially aromatic atoms
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      if (_vpa[atom->GetIdx()])
        PropagatePotentialAromatic(atom);

    //select root atoms
    SelectRootAtoms(mol);

    ExcludeSmallRing(mol); //remove 3 membered rings from consideration

    //loop over root atoms and look for aromatic rings
    _visit.clear();
    _visit.resize(mol.NumAtoms()+1);
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      if (_root[atom->GetIdx()])
        CheckAromaticity(atom,14);

    //for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    //	  if (atom->IsAromatic())
    //		  cerr << "aro = " <<atom->GetIdx()  << endl;

    //for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
    //if (bond->IsAromatic())
    //cerr << bond->GetIdx() << ' ' << bond->IsAromatic() << endl;
  }

  /** \brief Traverse a potentially aromatic cycle starting at @p root.
      \return  True if the cycle is likely aromatic
      \param root  The initial, "root" atom in traversing this ring
      \param atom  The current atom to visit and check
      \param prev  The bond traversed in moving to this @p atom
      \param er    The min and max number of pi electrons for this ring
      \param depth The maximum number of atoms to visit in a ring (e.g., 6)

      This method traverses a potentially aromatic ring, adding up the possible
      pi electrons for each atom. At the end (e.g., when @p atom == @p root)
      the Huekel 4n+2 rule is checked to see if there is a possible electronic
      configuration which corresponds to aromaticity.
  **/
  bool OBAromaticTyper::TraverseCycle(OBAtom *root, OBAtom *atom, OBBond *prev,
                                      std::pair<int,int> &er,int depth)
  {
    if (atom == root)
      {
        int i;
        for (i = er.first;i <= er.second;++i)
          if (i%4 == 2 && i > 2)
            return(true);

        return(false);
      }

    if (!depth || !_vpa[atom->GetIdx()] || _visit[atom->GetIdx()])
      return(false);

    bool result = false;

    depth--;
    er.first  += _velec[atom->GetIdx()].first;
    er.second += _velec[atom->GetIdx()].second;

    _visit[atom->GetIdx()] = true;
    OBAtom *nbr;
    vector<OBBond*>::iterator i;
    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
      if (*i != prev && (*i)->IsInRing() && _vpa[nbr->GetIdx()])
        {
          if (TraverseCycle(root,nbr,(OBBond*)(*i),er,depth))
            {
              result = true;
              ((OBBond*) *i)->SetAromatic();
            }
        }

    _visit[atom->GetIdx()] = false;
    if (result)
      atom->SetAromatic();

    er.first  -= _velec[atom->GetIdx()].first;
    er.second -= _velec[atom->GetIdx()].second;

    return(result);
  }

  void OBAromaticTyper::CheckAromaticity(OBAtom *atom,int depth)
  {
    OBAtom *nbr;
    vector<OBBond*>::iterator i;

    pair<int,int> erange;
    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
      if ((*i)->IsInRing()) // check all rings, regardless of assumed aromaticity
        {
          erange = _velec[atom->GetIdx()];

          if (TraverseCycle(atom,nbr,(OBBond *)*i,erange,depth-1))
            {
              atom->SetAromatic();
              ((OBBond*) *i)->SetAromatic();
            }
        }
  }

  void OBAromaticTyper::PropagatePotentialAromatic(OBAtom *atom)
  {
    int count = 0;
    OBAtom *nbr;
    vector<OBBond*>::iterator i;

    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
      if ((*i)->IsInRing() && _vpa[nbr->GetIdx()])
        count++;

    if (count < 2)
      {
        _vpa[atom->GetIdx()] = false;
        if (count == 1)
          for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
            if ((*i)->IsInRing() && _vpa[nbr->GetIdx()])
              PropagatePotentialAromatic(nbr);
      }
  }

  /**
   * \brief Select the root atoms for traversing atoms in rings.
   *
   * Picking only the begin atom of a closure bond can cause
   * difficulties when the selected atom is an inner atom
   * with three neighbour ring atoms. Why ? Because this atom
   * can get trapped by the other atoms when determining aromaticity,
   * because a simple visited flag is used in the
   * OBAromaticTyper::TraverseCycle() method.
   *
   * Ported from JOELib, copyright Joerg Wegner, 2003 under the GPL version 2
   * Improved by Fabian (fab5) in 2009 -- PR#2889708
   *
   * @param mol the molecule
   * @param avoidInnerRingAtoms inner closure ring atoms with more than 2 neighbours will be avoided
   *
   */
  void OBAromaticTyper::SelectRootAtoms(OBMol &mol, bool avoidInnerRingAtoms)
  {
    OBBond *bond;
    OBAtom *atom, *nbr, *nbr2;
    OBRing *ring;
    //    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j, l, nbr2Iter;
    vector<OBRing*> sssRings = mol.GetSSSR();
    vector<OBRing*>::iterator k;

    int rootAtom;
    int ringNbrs;
    int heavyNbrs;
    int newRoot = -1;
    vector<int> tmpRootAtoms;
    vector<int> tmp;

    vector<OBBond*> cbonds;
    vector< vector<OBRing*> > ringAtoms; // store ring pointers on an atom basis

    //generate list of closure bonds
    for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
      {
        if( bond->IsClosure() )
          {
            cbonds.push_back(bond);
            if(avoidInnerRingAtoms)
              {
                tmpRootAtoms.push_back(bond->GetBeginAtomIdx());
              }
          }
      }

    if(avoidInnerRingAtoms)
      {
        //for every atom fill vector with ring pointer it's associated with
        ringAtoms.resize(mol.NumAtoms()+1);
        for (k = sssRings.begin();k != sssRings.end();++k)
          {
            tmp = (*k)->_path;
            for (unsigned int j (0),j_end(tmp.size()); j < j_end; ++j)
              {
                ringAtoms[tmp[j]].push_back(*k);
              }
          }
      }


    //loop over closure bonds
    for(OBBondIterator bd(cbonds.begin()),bd_end(cbonds.end());bd!=bd_end;++bd)
      {
        bond = *bd;

        // BASIC APPROACH
        // pick beginning atom at closure bond
        // this is really ready, isn't it ! ;-)
        rootAtom = bond->GetBeginAtomIdx();
        _root[rootAtom] = true;

        // EXTENDED APPROACH
        if (avoidInnerRingAtoms)
          {
            // count the number of neighbor ring atoms
            atom = mol.GetAtom(rootAtom);
            ringNbrs = heavyNbrs = 0;

            for (nbr = atom->BeginNbrAtom(l);nbr;nbr = atom->NextNbrAtom(l))
              {
                // we can get this from atom->GetHvyValence()
                // but we need to find neighbors in rings too
                // so let's save some time
                if (!nbr->IsHydrogen())
                  {
                    heavyNbrs++;
                    if (nbr->IsInRing())
                      ringNbrs++;
                  }

                // if this atom has more than 2 neighbor ring atoms
                // we could get trapped later when traversing cycles
                // which can cause aromaticity false detection
                newRoot = -1;

                if (ringNbrs > 2)
                  {
                    // try to find another root atom
                    // only loop over rings which contain rootAtom
                    for(k = ringAtoms[rootAtom].begin() ; k != ringAtoms[rootAtom].end(); ++k)
                      {
                        ring = (*k);
                        tmp = ring->_path;

                        bool checkThisRing = false;
                        int rootAtomNumber=0;
                        int idx=0;
                        // avoiding two root atoms in one ring !
                        for (unsigned int j = 0; j < tmpRootAtoms.size(); ++j)
                          {
                            idx= tmpRootAtoms[j];
                            if(ring->IsInRing(idx))
                              {
                                rootAtomNumber++;
                                if(rootAtomNumber>=2)
                                  break;
                              }
                          }
                        if(rootAtomNumber<2)
                          {
                            for (unsigned int j = 0; j < tmp.size(); ++j)
                              {
                                // find critical ring
                                if (tmp[j] == rootAtom)
                                  {
                                    checkThisRing = true;
                                  }
                                else
                                  {
                                    // second root atom in this ring ?
                                    if (_root[tmp[j]] == true)
                                      {
                                        // when there is a second root
                                        // atom this ring can not be
                                        // used for getting an other
                                        // root atom
                                        checkThisRing = false;

                                        break;
                                      }
                                  }
                              }
                          }

                        // check ring for getting another
                        // root atom to avoid aromaticity typer problems
                        if (checkThisRing)
                          {
                            // check if we can find another root atom
                            for (unsigned int m = 0; m < tmp.size(); ++m)
                              {
                                ringNbrs = heavyNbrs = 0;
                                for (nbr2 = (mol.GetAtom(tmp[m]))->BeginNbrAtom(nbr2Iter);
                                     nbr2;nbr2 = (mol.GetAtom(tmp[m]))->NextNbrAtom(nbr2Iter))
                                  {
                                    if (!nbr2->IsHydrogen())
                                      {
                                        heavyNbrs++;

                                        if (nbr2->IsInRing())
                                          ringNbrs++;
                                      }
                                  }

                                // if the number of neighboured heavy atoms is also
                                // the number of neighboured ring atoms, the aromaticity
                                // typer could be stuck in a local traversing trap
                                if (ringNbrs <= 2 && ring->IsInRing((mol.GetAtom(tmp[m])->GetIdx())))
                                  {
                                    newRoot = tmp[m];
                                  }
                              }
                          }
                      }

                    if ((newRoot != -1) && (rootAtom != newRoot))
                      {
                        // unset root atom
                        _root[rootAtom] = false;

                        // pick new root atom
                        _root[newRoot] = true;
                      }
                  } // if (ringNbrs > 2)

              } // end for
          } // if (avoid)
      } // end for(closure bonds)
  }

  void OBAromaticTyper::ExcludeSmallRing(OBMol &mol)
  {
    OBAtom *atom,*nbr1,*nbr2;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j,k;

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      if (_root[atom->GetIdx()])
        for (nbr1 = atom->BeginNbrAtom(j);nbr1;nbr1 = atom->NextNbrAtom(j))
          if ((*j)->IsInRing() && _vpa[nbr1->GetIdx()])
            for (nbr2 = nbr1->BeginNbrAtom(k);nbr2;nbr2 = nbr1->NextNbrAtom(k))
              if (nbr2 != atom && (*k)->IsInRing() && _vpa[nbr2->GetIdx()])
                if (atom->IsConnected(nbr2))
                  _root[atom->GetIdx()] = false;
  }

} //namespace OpenBabel;

//! \file typer.cpp
//! \brief Open Babel atom and aromaticity typer.
