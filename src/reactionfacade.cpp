/**********************************************************************
reactionfacade.cpp - A facade class to help interrogate and manipulate
                     reactions

Copyright (C) 2018 by Noel M. O'Boyle

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

#include <openbabel/mol.h>
#include <openbabel/reactionfacade.h>

#include <set>

namespace OpenBabel
{
  void OBReactionFacade::AssignComponentIds(bool wipe)
  {
    unsigned int compid = 1;

    OBMolAtomDFSIter iter(_mol);
    while(iter) { // for each connected component
      do { // for each atom in connected component
        if (wipe || !(iter->HasData("rxncomp")))
          SetComponentId(&*iter, compid);
      } while ((iter++).next());
      compid++;
    }
  }

  int OBReactionFacade::GetId(const char* idtype, OBAtom *atom)
  {
    int idval = 0;
    OBGenericData *data = atom->GetData(idtype);
    if (data) {
      OBPairInteger *pi = (OBPairInteger*)data;
      idval = pi->GetGenericValue();
    }
    return idval;
  }

  void OBReactionFacade::SetId(const char* idtype, OBAtom* atom, int idval)
  {
    // Note: this replaces any existing data
    OBGenericData *data = atom->GetData(idtype);
    OBPairInteger *pi;
    if (data) {
      pi = (OBPairInteger*)data;
      pi->SetValue(idval);
    }
    else {
      pi = new OBPairInteger();
      pi->SetAttribute(idtype);
      pi->SetValue(idval);
      atom->SetData(pi);
    }
  }

  bool OBReactionFacade::IsValid()
  {
    if (!_mol->IsReaction()) {
      obErrorLog.ThrowError(__FUNCTION__, "The molecule is not marked as a reaction. Use SetIsReaction().", obWarning);
      return false;
    }
    FOR_ATOMS_OF_MOL(atom, _mol) {
      if (!atom->HasData("rxncomp")) {
        obErrorLog.ThrowError(__FUNCTION__, "The molecule contains an atom that is missing reaction component information. Use SetComponentId().", obWarning);
        return false;
      }
      if (!atom->HasData("rxnrole")) {
        obErrorLog.ThrowError(__FUNCTION__, "The molecule contains an atom that is missing reaction role information. Use SetRole().", obWarning);
        return false;
      }
    }
    return true;
  }

  void OBReactionFacade::FindComponents()
  {
    std::set<unsigned int> reactant_components;
    std::set<unsigned int> product_components;
    std::set<unsigned int> agent_components;
    std::set<unsigned int> unassigned_components;

    FOR_ATOMS_OF_MOL(atom, _mol) {
      OBAtom* atm = &*atom;
      unsigned int component = GetComponentId(&*atom);
      switch (GetRole(&*atom)) {
      case REACTANT:
        reactant_components.insert(component);
        break;
      case PRODUCT:
        product_components.insert(component);
        break;
      case AGENT:
        agent_components.insert(component);
        break;
      default:
        unassigned_components.insert(component);
      }
    }
    // Convert to vector so we have random access - note: these will be in sorted order
    for(std::set<unsigned int>::iterator sit = reactant_components.begin(); sit != reactant_components.end(); ++sit)
      _reactant_components.push_back(*sit);
    for (std::set<unsigned int>::iterator sit = product_components.begin(); sit != product_components.end(); ++sit)
      _product_components.push_back(*sit);
    for (std::set<unsigned int>::iterator sit = agent_components.begin(); sit != agent_components.end(); ++sit)
      _agent_components.push_back(*sit);
    for (std::set<unsigned int>::iterator sit = unassigned_components.begin(); sit != unassigned_components.end(); ++sit)
      _unassigned_components.push_back(*sit);

    _found_components = true;
  }

  std::vector<unsigned int>* OBReactionFacade::GetComponentIds(OBReactionRole rxnrole)
  {
    if (!_found_components)
      FindComponents();

    switch (rxnrole) {
    case REACTANT:
      return &_reactant_components;
    case PRODUCT:
      return &_product_components;
    case AGENT:
      return &_agent_components;
    case NO_REACTIONROLE:
      return &_unassigned_components;
    }
    return (std::vector<unsigned int>*)0;
  }

  unsigned int OBReactionFacade::NumComponents(OBReactionRole rxnrole)
  {
    std::vector<unsigned int> *component_ids = GetComponentIds(rxnrole);
    return (unsigned int) component_ids->size();
  }

  bool OBReactionFacade::GetComponent(OBMol* mol, OBReactionRole rxnrole, unsigned int num)
  {
    std::vector<unsigned int> *component_ids = GetComponentIds(rxnrole);
    if (num >= component_ids->size())
      return false;
    OBBitVec atoms;
    unsigned int componentId = (*component_ids)[num];
    FOR_ATOMS_OF_MOL(atom, _mol) {
      if (GetRole(&*atom) == rxnrole && GetComponentId(&*atom) == componentId)
        atoms.SetBitOn(atom->GetIdx());
    }
    bool ok = _mol->CopySubstructure(*mol, &atoms);
    return ok;
  }

  bool OBReactionFacade::ReassignComponent(OBReactionRole oldrole, unsigned int num, OBReactionRole newrole)
  {
    std::vector<unsigned int> *component_ids = GetComponentIds(oldrole);
    if (num >= component_ids->size())
      return false;
    unsigned int componentId = (*component_ids)[num];
    FOR_ATOMS_OF_MOL(atom, _mol) {
      if (GetRole(&*atom) == oldrole && GetComponentId(&*atom) == componentId) {
        SetRole(&*atom, newrole);
      }
    }
    // remove the entry from the original component_ids and update the new one
    component_ids->erase(component_ids->begin() + num);
    GetComponentIds(newrole)->push_back(componentId);
    return true;
  }

  void OBReactionFacade::AddComponent(OBMol* mol, OBReactionRole rxnrole)
  {
    if (!_found_components)
      FindComponents();

    // Get the current maximum component id
    unsigned int max_compid = 0;
    if (!_product_components.empty())
      max_compid = _product_components.back();
    if (!_agent_components.empty() && _agent_components.back()>max_compid)
      max_compid = _agent_components.back();
    if (!_reactant_components.empty() && _reactant_components.back()>max_compid)
      max_compid = _reactant_components.back();
    if (!_unassigned_components.empty() && _unassigned_components.back()>max_compid)
      max_compid = _unassigned_components.back();

    int new_compid = max_compid+1;
    if (new_compid == 0)
      new_compid = 1;

    FOR_ATOMS_OF_MOL(atom, mol) {
      SetRole(&*atom, rxnrole);
      SetComponentId(&*atom, new_compid);
    }
    *_mol += *mol;

    GetComponentIds(rxnrole)->push_back(new_compid);
  }

} // namespace OpenBabel


//! \file reactionfacade.cpp
//! \brief A facade class to help interrogate and manipulate reactions
