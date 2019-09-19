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

namespace OpenBabel
{
  class OBMol;
  class OBAtom;
  class OBReactionFacadePrivate;

  /**
  * The various roles a reaction component can have
  * @sa OBReactionFacade
  */
  enum OBReactionRole {
    NO_REACTIONROLE, //!< no reaction role - useful for temporarily hiding a component
    REACTANT,        //!< reactant
    AGENT,           //!< agent, a term that includes solvents and catalysts 
    PRODUCT          //!< product
  };
  /**
   * \class OBReactionFacade reactionfacade.h <openbabel/reactionfacade.h>
   * \brief Facade to simplify manipulation of reactions stored as OBMol objects
   *
   * All of the information defining a reaction is stored as part of an OBMol.
   * This information is:
   * - a flag indicating that the OBMol represents a reaction (see OBMol::SetIsReaction(),
   *   and OBMol::IsReaction())
   * - #OBPairInteger data attached to each atom indicating the reaction role ("rxnrole")
   * - #OBPairInteger data attached to each atom indicating the reaction component Id ("rxncomp")
   *
   * Everything that may need to be done with reactions can be implemented using this
   * information. The purpose of this class is simply to provide convenience functions
   * that cover a range of likely use cases. If your use case is not included, it should
   * still be straightforward for you to manipulate the underlying information to
   * achieve your goal.
   *
   * @since version 3.0
   */
  class OBAPI OBReactionFacade
  {
  public:
    /**
     * @name Constructor
     * @{
     */
    /**
     * The constructor requires an OBMol. Member functions retrieve or manipulate reaction
     * data stored in this molecule.
     */
    OBReactionFacade(OBMol *mol);
    ~OBReactionFacade();
    /**
     * @}
     * @name Low-level methods
     * These are convenience functions that set/get reaction-related data on individual atoms.
     * The actual data is stored in the underlying OBMol as OBPairInteger data with the
     * keys "rxnrole" or "rxncomp" (reaction role, reaction component).
     * The other methods of this class are implemented using these methods. Care should be
     * taken using the setter functions, as (unlike the other
     * methods) their use can result in an inconsistent state. IsValid() may be used to check
     * whether this is the case.
     * @{
     */
    /**
     * Assigns a component Id to every atom in the molecule based on connected components.
     * If @p wipe is \c false, then any atoms that already have a component Id are skipped.
     */
    void AssignComponentIds(bool wipe = true);
    /**
     * Return the reaction role of an atom.
     */
    OBReactionRole GetRole(OBAtom *atom);
    /**
     * Return the component Id of an atom.
     */
    unsigned int GetComponentId(OBAtom *atom);
    /**
     * Set the reaction role of an atom.
     */
    void SetRole(OBAtom* atom, OBReactionRole rxnrole);
    /**
     * Set the component Id of an atom.
     */
    void SetComponentId(OBAtom* atom, unsigned int compid);
    /**
     * @}
     * @name Check validity of reaction data
     * @{
     */
    /**
     * Check whether this OBMol is in a valid state to represent a reaction.
     *
     * Note, first of all, that this does not check whether a reaction is
     * chemically valid. It simply checks whether the OBMol has the following
     * properties:
     * - (a) it is marked as a reaction
     * - (b) every atom is marked as having a reaction role
     * - (c) every atom is marked as belonging to a particular reaction component
     * - (d) every atom in a connected component is marked as belong to the same
     *       reaction component and having the same reaction role
     *
     * If (a) is not true, while it does not affect any of the methods of the
     * OBReactionFacade, the reaction may not be written appropriately by an
     * output format.
     *
     * If (b), (c) or (d) is not true, then the behavior of the toolkit is undefined
     * when handling this molecule.
     *
     * While not checked by this function, it is also a good idea to avoid reusing 
     * the same component Id for different reaction roles, e.g. having a reactant
     * with component Id 1 as well as a product with component Id 1. If you reassign
     * roles, then the duplication of component Ids may cause two components to
     * be merged into one.
     *
     */
    bool IsValid();
    /**
     * @}
     * @name High-level methods
     * These are a set of methods that manipulate reactions at the component level.
     * If you use these high-level methods, and subsequently modify the underlying
     * OBMol or its reaction data, you may need to clear cached information on the
     * reaction components via ClearInternalState().
     * @{
     */
    /**
    * Add a new component to the reaction
    */
    void AddComponent(OBMol *mol, OBReactionRole rxnrole);
    /**
    * Clear the list of found components.
    *
    * On first
    * use any of the high-level methods,
    * the number of components for each reaction role is
    * determined and each is assigned a number from 0 to the number of components.
    * This information is cached for future accesses. If you modify the underlying
    * OBMol in a way that invalidates this information (e.g. using the low-level
    * methods), you should call ClearInternalState() or create a new OBReactionFacade.
    */
    void ClearInternalState();
    /**
    * Copy a component from the reaction into the provided OBMol.
    
    * Note that all
    * reaction data is copied, and so this may be useful for building up a new
    * reaction composed of just some of the components of the original reaction.
    * If copying components from multiple reactions into a single one, to avoid
    * duplicate component Ids it may be easier to copy into an empty OBMol, and
    * then use AddComponent() to add this OBMol into the new reaction.
    *
    * \return A boolean indicating whether the reaction has a component with that number.
    */
    bool GetComponent(OBMol *mol, OBReactionRole rxnrole, unsigned int num);
    /**
    * Return the number of reaction components with a particular reaction role.
    */
    unsigned int NumComponents(OBReactionRole rxnrole);
    /**
    * Reassign the reaction role of a given component. If assigned to #NO_REACTIONROLE,
    * this may be used to hide the reaction component when writing to a reaction
    * output format.
    * 
    * It should be noted that reassigning a component changes the numbers assigned.
    * For example,
    * in a reaction with 2 reactants and 2 products (each numbered from 0 to 1), if you
    * reassign reactant number
    * 0 to be a product, it will now be product number 2 while former reactant number 1
    * will now be reactant number 0.
    * 
    * \return A boolean indicating whether the reaction has a component with that number.
    */
    bool ReassignComponent(OBReactionRole oldrole, unsigned int num, OBReactionRole newrole);
    /**
     * @}
     */

  private:
    OBReactionFacadePrivate *d;
  };
} // namespace OpenBabel
#endif // OB_REACTIONFACADE_H

//! \file reactionfacade.h
//! \brief A facade class to help interrogate and manipulate reactions
