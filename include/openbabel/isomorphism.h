/**********************************************************************
  isomophism.h - OBIsomorphismMapper class for finding isomorphisms

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
#ifndef OB_ISOMORPHISM_H
#define OB_ISOMORPHISM_H

#include <openbabel/mol.h>

namespace OpenBabel {

  class OBQuery;

  ///@addtogroup substructure Substructure Searching
  ///@{

  /**
   * @since version 2.3
   * @class OBIsomorphismMapper
   * @brief Abstract class defining interface for isomorphism (i.e. substructure) searches.
   *
   * The OBIsomorphism class is an abstract class which defines an interface for
   * performing isomorphism (i.e. substructure) searches. It uses a OBQuery and
   * tries to map this onto a queried OBMol. A single mapping is represented by
   * a OBIsomorphismMapper::Mapping which is a std::map mapping query indexes to
   * queried indexes. Both query and queried indexes in the map start from 0.
   * Multiple mappings can be stored in a OBIsomorphismMapper::Mappings object
   * which is a std::vector of OBIsomorphismMapper objects.
   *
   * Since this is an abstract class with pure virtual methods, this class can't
   * be instantiated directly. To get a pointer to a subclass, the GetInstance()
   * method can be used which also sets the query. Once an instance is obtained,
   * the desired mapping function can be used to perform the mapping (i.e. MapFirst(),
   * MapUnique() or MapAll()).
   *
   * A typical example:
   * @code
   * OBMol *queried;
   * // ... initialize queried ...
   * OBQuery *query = CompileSmilesQuery("c1ccccc1");
   * OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
   * OBIsomorphismMapper::Mappings maps = mapper->MapUnique(mol);
   *
   * std::cout << "found " << maps.size() << " unique mappings" << std::endl;
   *
   * delete mapper;
   * delete query;
   * @endcode
   *
   * All mapping methods take an optional mask parameter. This can be used to
   * restrict the search to a part of the queried OBMol. The masked atoms in
   * the OBBitVec are  indexed from 1. A special case of isomorphism search is
   * an automorphism search where the query and queried molecule are the same.
   * Automorphism searches can be done using the MapAll method but an additional
   * FindAutomorphisms() function is provided for convenience.
   */
  class OBAPI OBIsomorphismMapper
  {
    public:
      /**
       * Type for an individual mapping.
       */
      typedef std::vector< std::pair<unsigned int,unsigned int> > Mapping;
      /**
       * Type for a collection (std::vector) of Mapping objects.
       */
      typedef std::vector<Mapping> Mappings;

      /**
       * Constructor. OBIsomorphismMapper is an abstract class, use GetInstance() to
       * get an instance of a derived class.
       * @param query The search query.
       */
      OBIsomorphismMapper(OBQuery *query);
      virtual ~OBIsomorphismMapper();

      /**
       * Get a pointer to an instance of the specified @p algorithm. This pointer
       * has to be delted when the instance is no longer needed.
       * @param query The search query to be mapped.
       * @param algorithm The algorithm for the mapper.
       * @return OBIsomorphismMapper instance or 0 if there is no subclass implementing
       * the specified @p algorithm.
       */
      static OBIsomorphismMapper* GetInstance(OBQuery *query, const std::string &algorithm = std::string("VF2"));

      /**
       * Find a single mapping in @p queried.
       * @param queried The molecule to search.
       * @param mask A mask to restrict the search to a part of the queried molecule.
       * The default empty mask will result in all atoms being considered. The mask
       * indexes start from 1 (i.e. OBAtom::GetIdx()).
       */
      virtual void MapFirst(const OBMol *queried, Mapping &map, const OBBitVec &mask = OBBitVec()) = 0;
      /**
       * Find all unique mappings in @p queried. A mapping is unique when there is no previous
       * mapping covering the same queried atoms. For two mappings, some overlap is allowed but
       * at least one atom should be different.
       * @param queried The molecule to search.
       * @param mask A mask to restrict the search to a part of the queried molecule.
       * The default empty mask will result in all atoms being considered. The mask
       * indexes start from 1 (i.e. OBAtom::GetIdx()).
       */
      virtual void MapUnique(const OBMol *queried, Mappings &maps, const OBBitVec &mask = OBBitVec()) = 0;
      /**
       * Find all mappings in @p queried. This function is used by FindAutomorphisms()
       * with a query that is a copy of the queried molecule (taking the mask into
       * account).
       * @param queried The molecule to search.
       * @param mask A mask to restrict the search to a part of the queried molecule.
       * The default empty mask will result in all atoms being considered. The mask
       * indexes start from 1 (i.e. OBAtom::GetIdx()).
       */
      virtual void MapAll(const OBMol *queried, Mappings &maps, const OBBitVec &mask = OBBitVec()) = 0;

      /**
       * Set the timeout in seconds.
       */
      void SetTimeout(unsigned int seconds) { m_timeout = seconds; }

      std::size_t GetMemory() const { return m_memory; }
      void SetMemory(std::size_t memory) { m_memory = memory; }
      std::size_t GetMaxMemory() const { return m_maxMemory; }
      void SetMaxMemory(std::size_t maxMemory) { m_maxMemory = maxMemory; }

    protected:
      OBQuery *m_query; //!< The search query.
      unsigned int m_timeout; //!< The timeout in seconds
      std::size_t m_memory, m_maxMemory;
  };


  typedef OBIsomorphismMapper::Mapping Automorphism;
  typedef OBIsomorphismMapper::Mappings Automorphisms;

  inline bool MapsTo(const OBIsomorphismMapper::Mapping &map, unsigned int queryIndex, unsigned int &queriedIndex)
  {
    OBIsomorphismMapper::Mapping::const_iterator i;
    for (i = map.begin(); i != map.end(); ++i)
      if (i->first == queryIndex) {
        queriedIndex = i->second;
        return true;
      }

    return false;
  }

  /**
   * Find the automorphisms of a molecule by using a OBIsomorphismMapper.
   * @since version 2.3
   */
  OBAPI bool FindAutomorphisms(OBMol *mol, Automorphisms &aut, const OBBitVec &mask = OBBitVec());
  OBAPI bool FindAutomorphisms(OBMol *mol, Automorphisms &aut, const std::vector<unsigned int> &symmetry_classes,
      const OBBitVec &mask = OBBitVec());



  /**
   * @page substructure Substructure Search
   *
   * Substructure searching is finding a mapping for a query to a target molecule.
   * Such a mapping is also known as a graph isomorphism. A graph isomorphism maps
   * the vertexes (i.e. atoms) from the query to vertexes in the molecule such that
   * two vertexes adjacent in the query are also adjacent in the target molecule.
   * In other words, no bonds are broken and no new bonds are formed.
   *
   * @section smarts SMARTS Substructure Search
   * Smarts is an extension of smiles to create powerful queries. Smarts substructure
   * search has been available in OpenBabel for many years. It is also used for many
   * of OpenBabel's algorithms. Although smarts is only a syntax for queries, the
   * implementation has it's own matching algorithm. For many purposes smarts are the
   * easiest way to do substructure searches. See the OBSmartsPattern documentation
   * for details on how to use smarts.
   *
   * @section query Queries
   * Since OpenBabel version 2.3, there are some classes for representing generic
   * queries. The OBQuery, OBQueryAtom and OBQueryBond class define interfaces that
   * can be reimplemented to get custom behavior. The classes also contain some
   * methods to access topological information which are used by the mapping
   * algorithms. The default implementations allow very simple exact substructure
   * searches to be performed but subclassing allows very advanced queries to be
   * used (e.g. smarts).
   *
   * While it is possible to construct these queries manually, "compilers" are
   * provided to convert a query representation to a OBQuery object. Currently,
   * only two exact substructure search compilers exist. The first is
   * CompileMoleculeQuery which converts an OBMol object to an OBQuery object.
   * The second is CompileSmilesQuery and converts a smiles string to an OBQuery
   * object.
   *
   * @code
   * OBMol *mol = new OBMol;
   *
   * // ... read molecule ...
   *
   * OBQuery *query;
   * query = CompileMoleculeQuery(mol);
   * query = CompileSmilesQuery("c1ccccc1CC(=O)O");
   * @endcode
   *
   * @section mapping Mapping Isomorphisms
   * The OBIsomorphismMapper class defined an interface for mapping queries to
   * target molecules. Multiple implementations can be added but they all do the
   * same. The MapFirst, MapUnique and MapAll methods are used for gettings the
   * map(s).
   *
   * @subsection MapFirst
   * This method returns the first map found. The main reason for getting only
   * one map is improved performance since it is considerably faster than
   * MapUnique and MapAll. However, depending on the use case a single map is
   * all that is needed. For example, to check if a molecule in a database
   * contains a phenyl ring, a single mapping is enough.
   *
   * @subsection MapUnique
   * MapUnique returns all unique maps. A map is considered unique if there is
   * no other map covering exactly the same atoms in the target molecule. For
   * example, when a phenyl query is performed on a molecule with 2 phenyl rings,
   * MapUnique will return 2 maps. These 2 maps are selected from the 24 found
   * non-duplicate maps (6 atoms to start from * 2 directions (CW/CCW) * 2 rings).
   *
   * @subsection MapAll
   * MapAll returns all non-duplicate maps. For example, when a phenyl query is
   * performed on a molecule with 2 phenyl rings, MapAll will return 24 maps
   * (6 atoms to start from * 2 directions (CW/CCW) * 2 rings).
   *
   * @section automorphisms Automorphisms
   * The automorphisms of a graph or molecule are a group of isomorphism mappings
   * of the molecule onto itself (i.e. the query and target are the same). The
   * automorphisms make it easy to take symmetry into account. See FindAutomorphisms
   * for detials.
   *
   *
   */

  ///@}

}


#endif
