/**********************************************************************
alias.h - OBGenericData class to hold alias information on atoms
Copyright (C) Copyright (C) 2007 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_ALIAS_H
#define OB_ALIAS_H

#include <vector>
#include <openbabel/mol.h>

// This macro is used in DLL builds. If it has not
// been set in babelconfig.h, define it as nothing.
#ifndef OBAPI
  #define OBAPI
#endif

namespace OpenBabel
{
  /// @addtogroup data Data classes 
  //@{

  /** @class OBAliasData alias.h <openbabel/alias.h>
      @brief Indicate atoms as aliases for larger functional groups
      @since version 2.2

      An object of this class can be attached to an OBAtom if it is considered
      to be a placeholder for an alias, such as Ph, Ser, etc. 

      If the alias has not been interpreted chemically, the element type of the
      placeholder atom should be set to Xx so that the molecule is not interpreted 
      incorrectly by formats which do not consider this class.

      If the alias has been interpreted chemically, the alias may remain as extra
      information or as a hint for an alternative representation, for example to
      a chemical drawing program. The m_expandedatoms vector would then contains
      the indices of the atoms to which the alias is an alternative.
  */
  class OBAPI OBAliasData : public OBGenericData
  {
    private:
      std::string               m_alias;         //!< the alias
      std::vector<unsigned int> m_expandedatoms; //!< the atoms for this alias, if any
    public:
      virtual OBGenericData* Clone(OBBase* /*parent*/) const 
      { return new OBAliasData(*this); }
      //! @name OBAliasData methods
      //@{
      /**
       * Constructor.
       */
      OBAliasData(): OBGenericData("Alias", OBGenericDataType::AliasData)
      { 
      }
      /**
       * Set the alias string.
       * @param alias The alias string.
       */
      void SetAlias(const std::string& alias) { m_alias = alias; }
      /**
       * Set the alias string.
       * @param alias The alias string.
       */
      void SetAlias(const char* alias) { m_alias = alias;}
      /**
       * Get the alias string.
       * @return The alias string.
       */
      std::string GetAlias() const { return m_alias; }
      /**
       * Set the atom indices of the atoms that could or do replace the alias.
       * @param atoms The atom indices.
       */
      void SetExpandedAtoms(std::vector<unsigned int>& atoms) { m_expandedatoms = atoms; }
      /**
       * Return the indices of the atoms that the alias could or does replace
       */
      std::vector<unsigned int> GetExpandedAtoms() const { return m_expandedatoms; }
      /**
       * @return True if the atom indices are expanded.
       */
      bool IsExpanded() const { return !m_expandedatoms.empty(); }
      /**
       * This functions tries to expand the alias by interpreting it.
       * @param mol The molecule to expand.
       * @param atomindex The atom index for the atom with this alias. The molecule 
       * will be expanded starting at this atom.
       *
       * @warning Experimental code!
       */
      bool Expand(OBMol& mol, const unsigned int atomindex);
      //@}
  };

  //@} group

} // namespace

#endif // OB_ALIAS_H

//! @file alias.h
//! @brief OBGenericData class to for atom alias data (e.g., in 2D drawing programs for "COOH")
