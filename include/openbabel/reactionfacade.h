/**********************************************************************
reactionfacade.h - A facade class to help interrogate and manipulate
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

#ifndef OB_REACTIONFACADE_H
#define OB_REACTIONFACADE_H

#include <openbabel/babelconfig.h>

namespace OpenBabel
{
  enum OBReactionRole { NO_REACTIONROLE, REACTANT, AGENT, PRODUCT };

  class OBAPI OBReactionFacade
  {
  public:
    OBReactionFacade(OBMol *mol): _mol(mol), _found_components(false)
    {
    };
    OBReactionRole GetRole(OBAtom *atom)
    {
      int rxnrole = GetId("rxnrole", atom);
      switch(rxnrole) {
      default: case 0: return NO_REACTIONROLE;
      case 1: return REACTANT;
      case 2: return AGENT;
      case 3: return PRODUCT;
      }
    }
    void SetRole(OBAtom* atom, OBReactionRole rxnrole)
    {
      SetId("rxnrole", atom, rxnrole);
    }
    unsigned int GetComponentId(OBAtom *atom)
    {
      return GetId("rxncomp", atom);
    }
    void SetComponentId(OBAtom* atom, unsigned int compid)
    {
      SetId("rxncomp", atom, compid);
    }
    bool IsValid();
    bool GetComponent(OBMol *mol, OBReactionRole rxnrole, unsigned int num);
    unsigned int NumComponents(OBReactionRole rxnrole);
    bool ReassignComponent(OBReactionRole oldrole, unsigned int num, OBReactionRole newrole);
    void AddComponent(OBMol *mol, OBReactionRole rxnrole); 

    void AssignComponentIds(bool wipe=true);
    void ClearInternalState()
    {
      _found_components = false;
    }

  private:
    OBMol* _mol;
    bool _found_components;
    std::vector<unsigned int> _unassigned_components;
    std::vector<unsigned int> _reactant_components;
    std::vector<unsigned int> _product_components;
    std::vector<unsigned int> _agent_components;

    int GetId(const char* idtype, OBAtom *atom);
    void SetId(const char* idtype, OBAtom* atom, int idval);
    void FindComponents();
    std::vector<unsigned int>* GetComponentIds(OBReactionRole rxnrole);
  };
} // namespace OpenBabel
#endif // OB_REACTIONFACADE_H

//! \file reactionfacade.h
//! \brief A facade class to help interrogate and manipulate reactions
