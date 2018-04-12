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
  class OBAPI OBReactionFacade
  {
  public:
    OBReactionFacade(OBMol *mol): _mol(mol), _found_components(false)
    {
    };
    unsigned int GetRole(OBAtom *atom)
    {
      return GetId("rxnrole", atom);
    }
    void SetRole(OBAtom* atom, unsigned int rxnrole)
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
    unsigned int NumReactants()
    {
      return NumComponents(0);
    }
    unsigned int NumProducts()
    {
      return NumComponents(1);
    }
    unsigned int NumAgents()
    {
      return NumComponents(2);
    }
    bool GetReactant(OBMol *mol, unsigned int num)
    {
      return GetComponent(0, mol, num);
    }
    bool GetProduct(OBMol *mol, unsigned int num)
    {
      return GetComponent(1, mol, num);
    }
    bool GetAgent(OBMol *mol, unsigned int num)
    {
      return GetComponent(2, mol, num);
    }
    void AssignComponentIds(bool wipe=true);

  private:
    OBMol* _mol;
    bool _found_components;
    std::vector<unsigned int> _reactant_components;
    std::vector<unsigned int> _product_components;
    std::vector<unsigned int> _agent_components;

    unsigned int GetId(const char* idtype, OBAtom *atom);
    void SetId(const char* idtype, OBAtom* atom, unsigned int idval);
    void FindComponents();
    std::vector<unsigned int>* GetComponentIds(unsigned int component_type);
    unsigned int NumComponents(unsigned int component_type);
    bool GetComponent(unsigned int component_type, OBMol *mol, unsigned int num);
  };
} // namespace OpenBabel
#endif // OB_REACTIONFACADE_H

//! \file reactionfacade.h
//! \brief A facade class to help interrogate and manipulate reactions
