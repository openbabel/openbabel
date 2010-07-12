/**********************************************************************
  query.h - OBQuery, OBQueryAtom & OBQueryBond classes

  Copyright (C) 2010 by Tim Vandermeersch
 
  This file is part of the Open Babel project.
  For more information, see <http://openbabel.sourceforge.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/
#ifndef OB_QUERY_H
#define OB_QUERY_H

#include <openbabel/mol.h>

namespace OpenBabel {

  class OBQueryBond;

  /**
   * The OBQueryAtom class defines an interface for query atoms. The class provides
   * some general methods and properties to access the topology information. The Matches
   * method can be reimplemented in subclasses to get custom matching behavior.
   *
   * The default Matches implementation only checks the atomic number.
   *
   * @sa OBQuery OBQueryBond OBIsomorphismMapper
   * @since version 2.3
   */
  class OBAPI OBQueryAtom
  {
    public:
      friend class OBQuery;
      friend class OBQueryBond;
      /**
       * Constructor.
       * @param atomicNum The atomic number for this query atom.
       */
      OBQueryAtom(int atomicNum = 6) : m_atomicNum(atomicNum) {}
      /**
       * Get the index for this query atom. Atoms are indexed starting from 0.
       * This method is used by OBIsomorphismMapper implementations.
       */
      unsigned int GetIndex() const 
      {
        return m_index;
      }
      /**
       * Get the query bonds for this atom.
       * This method is used by OBIsomorphismMapper implementations.
       */
      const std::vector<OBQueryBond*>& GetBonds() const
      {
        return m_bonds;
      }
      /**
       * Get the neighbor query atoms.
       * This method is used by OBIsomorphismMapper implementations.
       */
      const std::vector<OBQueryAtom*>& GetNbrs() const
      {
        return m_nbrs;
      }
      /**
       * This is the match method to verify if an OBQueryAtom and OBAtom class match.
       * The default implementation only checks if the atomic numbers match. Reimplement
       * this method in a subclass for more advances matching.
       * This method is used by OBIsomorphismMapper implementations.
       * @param atom The OBAtom object to compare this OBQueryAtom with.
       */
      virtual bool Matches(const OBAtom *atom) const
      {
        if (atom->GetAtomicNum() != m_atomicNum)
          return false;
        return true;
      }
    protected:
      unsigned int m_index;
      int m_atomicNum;
      std::vector<OBQueryBond*> m_bonds;
      std::vector<OBQueryAtom*> m_nbrs;
  };

  /**
   * The OBQueryBond class defines an interface for query bonds. The class provides
   * some general methods and properties to access the topology information. The Matches
   * method can be reimplemented in subclasses to get custom matching behavior.
   *
   * The default Matches implementation only checks if the bonds are both aromatic,
   * otherwise the bond orders are compared.
   *
   * @sa OBQuery OBQueryAtom OBIsomorphismMapper
   * @since version 2.3
   */
  class OBAPI OBQueryBond
  {
    public:
      friend class OBQuery;
      /**
       * Constructor.
       */
      OBQueryBond(OBQueryAtom *begin, OBQueryAtom *end, int order = 1, bool aromatic = false) : 
          m_begin(begin), m_end(end), m_order(order), m_aromatic(aromatic) 
      {
        m_begin->m_bonds.push_back(this);
        m_end->m_bonds.push_back(this);
        m_begin->m_nbrs.push_back(m_end);
        m_end->m_nbrs.push_back(m_begin);
      }
      /**
       * Get the index for this query bonds. Query bonds are indexed starting from 0.
       */
      unsigned int GetIndex() const 
      {
        return m_index;
      }
      /**
       * Get the begin atom.
       */
      OBQueryAtom* GetBeginAtom() const { return m_begin; }
      /**
       * Get the end atom.
       */
      OBQueryAtom* GetEndAtom() const { return m_end; }
      /**
       * This is the match method to verify if an OBQueryBond and OBBond class match.
       * The default implementation checks if both bonds are aromatic and compares the 
       * bond orders otherwise. Reimplement this method in a subclass for more 
       * advances matching.
       * This method is used by OBIsomorphismMapper implementations.
       * @param atom The OBBond object to compare this OBQueryBond with.
       */
      virtual bool Matches(const OBBond *bond) const
      {
        if (m_aromatic)
          return bond->IsAromatic();
        return bond->GetBondOrder() == m_order;      
      }
    protected:
      unsigned int m_index;
      OBQueryAtom *m_begin, *m_end;
      int m_order;
      bool m_aromatic;
  };

  /**
   * @since version 2.3
   */
  class OBAPI OBQuery
  {
    public:
      unsigned int NumAtoms() const 
      { 
        return m_atoms.size();
      }
      unsigned int NumBonds() const
      {
        return m_bonds.size();
      }
      const std::vector<OBQueryAtom*>& GetAtoms() const 
      {
        return m_atoms;
      }
      const std::vector<OBQueryBond*>& GetBonds() const
      {
        return m_bonds;
      }
      OBQueryBond* GetBond(OBQueryAtom *begin, OBQueryAtom *end) const
      {
        for (unsigned int i = 0; i < begin->GetBonds().size(); ++i)
          if (begin->GetNbrs()[i] == end)
            return begin->GetBonds()[i];        
        return 0;
      }
      void AddAtom(OBQueryAtom *atom)
      {
        atom->m_index = m_atoms.size();
        m_atoms.push_back(atom);
      }
      void AddBond(OBQueryBond *bond)
      {
        bond->m_index = m_bonds.size();
        m_bonds.push_back(bond);
      }
    protected:
      std::vector<OBQueryAtom*> m_atoms;
      std::vector<OBQueryBond*> m_bonds;
  };
  
  /**
   * @since version 2.3
   */
  OBAPI OBQuery* CompileMoleculeQuery(OBMol *mol);
  
  /**
   * @since version 2.3
   */
  OBAPI OBQuery* CompileSmilesQuery(const std::string &smiles);

}


#endif
