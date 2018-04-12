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
    unsigned int compid = 0;

    OBMolAtomDFSIter iter(_mol);
    while(iter) { // for each connected component
      do { // for each atom in connected component
        if (wipe || !(iter->HasData("rxncomp")))
          SetComponentId(&*iter, compid);
      } while ((iter++).next());
      compid++;
    }
  }

  unsigned int OBReactionFacade::GetId(const char* idtype, OBAtom *atom)
  {
    unsigned int idval = 0;
    OBGenericData *data = atom->GetData(idtype);
    if (data) {
      OBPairInteger *pi = (OBPairInteger*)data;
      idval = pi->GetGenericValue();
    }
    return idval;
  }

  void OBReactionFacade::SetId(const char* idtype, OBAtom* atom, unsigned int idval)
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

  void OBReactionFacade::FindComponents()
  {
    std::set<unsigned int> reactant_components;
    std::set<unsigned int> product_components;
    std::set<unsigned int> agent_components;

    FOR_ATOMS_OF_MOL(atom, _mol) {
      unsigned int component = GetComponentId(&*atom);
      switch (GetRole(&*atom)) {
      case 0: // reactant
        reactant_components.insert(component);
        break;
      case 1: // product
        product_components.insert(component);
        break;
      case 2: // agent
        agent_components.insert(component);
        break;
      }
    }
    // Convert to vector so we have random access - note: these will be in sorted order
    for(std::set<unsigned int>::iterator sit = reactant_components.begin(); sit != reactant_components.end(); ++sit)
      _reactant_components.push_back(*sit);
    for (std::set<unsigned int>::iterator sit = product_components.begin(); sit != product_components.end(); ++sit)
      _product_components.push_back(*sit);
    for (std::set<unsigned int>::iterator sit = agent_components.begin(); sit != agent_components.end(); ++sit)
      _agent_components.push_back(*sit);

    _found_components = true;
  }

  std::vector<unsigned int>* OBReactionFacade::GetComponentIds(unsigned int component_type)
  {
    if (!_found_components)
      FindComponents();

    switch (component_type) {
    case 0:
      return &_reactant_components;
    case 1:
      return &_product_components;
    case 2:
      return &_agent_components;
    }
    return (std::vector<unsigned int>*)0;
  }

  unsigned int OBReactionFacade::NumComponents(unsigned int component_type)
  {
    std::vector<unsigned int> *component_ids = GetComponentIds(component_type);
    return (unsigned int) component_ids->size();
  }

  bool OBReactionFacade::GetComponent(unsigned int component_type, OBMol *mol, unsigned int num)
  {
    std::vector<unsigned int> *component_ids = GetComponentIds(component_type);
    if (num >= component_ids->size())
      return false;
    OBBitVec atoms;
    unsigned int componentId = (*component_ids)[num];
    FOR_ATOMS_OF_MOL(atom, _mol) {
      if (GetRole(&*atom) == 0 && GetComponentId(&*atom) == componentId)
        atoms.SetBitOn(atom->GetIdx());
    }
    bool ok = true; // _mol->CopySubstructure(mol, atoms);
    return ok;
  }

} // namespace OpenBabel


//! \file reactionfacade.cpp
//! \brief A facade class to help interrogate and manipulate reactions
