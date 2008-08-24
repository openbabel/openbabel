/**********************************************************************
Copyright (C) 2007 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_ATOMCLASS_H
#define OB_ATOMCLASS_H

#include <openbabel/base.h>
#include <vector>

#ifdef HAVE_SSTREAM
#include <sstream>
#endif

// This macro is used in DLL builds. If it has not
// been set in babelconfig.h, define it as nothing.
#ifndef OBAPI
  #define OBAPI
#endif

namespace OpenBabel
{
  /// @addtogroup data Data classes 
  //@{
  
  /** @class OBAtomClassData atomclass.h <openbabel/atomclass.h>
      @brief Handle atom classes in reaction SMILES/SMIRKS
      @since version 2.2

      Class for attaching to OBMol to hold the info from SMILES like [C:2]
      Useful for reaction SMILES (SMIRKS). It influences the atom id attribute in CML.
      Not all atoms need have an atom class. 
      The atom class can be any positive or negative integer.
  */
  class OBAPI OBAtomClassData : public OBGenericData
  {
    protected:
      std::map<unsigned int, int> m_map; //!< index is atom index; value is class
    public:
      virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBAtomClassData(*this); }
      //! @name OBAtomClassData member functions
      //@{
      /**
       * Constructor.
       */
      OBAtomClassData(): OBGenericData("Atom Class", OBGenericDataType::AtomClassData) 
      {
      }
      /**
       * Erase contents.
       */
      void Clear(){ m_map.clear(); }
      /**
       * Add an individual value
       * @param index The atom index.
       * @param cl The atom class for index @p index.
       */
      void Add(unsigned int index, int cl) { m_map[index] = cl; }
      /**
       * @param index The atom index.
       * @return True if there is an atom class entry for @p index.
       */
      bool HasClass(unsigned int index) const { return m_map.find(index) != m_map.end(); }
      /**
       * @param index The atom index.
       * @return Value of class index (Test with HasClass first).
       */
      int GetClass(unsigned int index) const 
      {
        std::map<unsigned int, int>::const_iterator pos = m_map.find(index);
        if(pos != m_map.end())
          return pos->second;
        return -9999;
      }
      /**
       * If there is an entry for @p index, return ":n" where n is the atomclass value;
       * otherwise return an empty string.
       * @param index The atom index.
       */
      std::string GetClassString(unsigned int index)
      {
        std::stringstream ss;
        std::map<unsigned int, int>::const_iterator pos = m_map.find(index);
        if(pos != m_map.end())
          ss << ':' << pos->second;
        return ss.str();
      }
      /**
       * @return The number of atoms for which an atom class is set.
       */
      unsigned int Size() { return m_map.size(); }
      //@}
  };

  //@} group

} //namespace

#endif // OB_ATOMCLASS_H

//! @file atomclass.h
//! @brief Handle atom classes in reaction SMILES/SMIRKS
