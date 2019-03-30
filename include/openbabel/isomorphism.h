/**********************************************************************
  isomophism.h - OBIsomorphismMapper class for finding isomorphisms.

  Copyright (C) 2010 by Tim Vandermeersch

  This file is part of the Open Babel project.
  For more information, see <http://openbabel.org/>

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

#include <openbabel/bitvec.h>

namespace OpenBabel {

  class OBQuery;
  class OBMol;

  ///@addtogroup substructure Substructure Searching
  ///@{

  /**
   * @class OBIsomorphismMapper isomorphism.h <openbabel/isomorphism.h>
   * @since version 2.3
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
       * @typedef std::vector< std::pair<unsigned int,unsigned int> > Mapping
       * Type for an individual mapping.
       */
      typedef std::vector< std::pair<unsigned int,unsigned int> > Mapping;
      /**
       * @typedef std::vector<OBIsomorphismMapper::Mapping> Mappings
       * Type for a collection (std::vector) of Mapping objects.
       */
      typedef std::vector< Mapping > Mappings;

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
       * @param map Reference to the object to store the result in.
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
       * @param maps Reference to the object to store the results in.
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
       * @param maps Reference to the object to store the results in.
       * @param mask A mask to restrict the search to a part of the queried molecule.
       * @param maxMemory Memory limit for the @p maps object in bytes. Default is 300MB.
       * The default empty mask will result in all atoms being considered. The mask
       * indexes start from 1 (i.e. OBAtom::GetIdx()).
       */
      virtual void MapAll(const OBMol *queried, Mappings &maps, const OBBitVec &mask = OBBitVec(), std::size_t maxMemory = 3000000) = 0;

      /**
       * @class Functor isomorphism.h <openbabel/isomorphism.h>
       * @brief Functor base class to be used in combination with MapGeneric.
       * @see @ref MapGeneric
       * @since 2.3
       */
      class Functor
      {
        public:
          virtual ~Functor() {}
          /**
           * This function is called every time an isomorphism is discovered.
           * Returing true, will abort the mapping process. The map is passed
           * as non-const reference and may be modified (e.g. swap).
           *
           * @see @ref MapGeneric
           */
          virtual bool operator()(Mapping &map) = 0;
      };
      /**
       * Find all mappings in @p queried. The functor will be called when a mapping is found.
       * @param functor The functor to handle found mappings.
       * @param queried The molecule to search.
       * @param mask A mask to restrict the search to a part of the queried molecule.
       * The default empty mask will result in all atoms being considered. The mask
       * indexes start from 1 (i.e. OBAtom::GetIdx()).
       */
      virtual void MapGeneric(Functor &functor, const OBMol *queried, const OBBitVec &mask = OBBitVec()) = 0;


      /**
       * Set the timeout in seconds.
       */
      void SetTimeout(unsigned int seconds) { m_timeout = seconds; }

    protected:
      OBQuery *m_query; //!< The search query.
      unsigned int m_timeout; //!< The timeout in seconds
  };

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
   * @typedef OBIsomorphismMapper::Mapping Automorphism
   * @brief A single automorphic permutation.
   * @since 2.3
   */
  typedef OBIsomorphismMapper::Mapping Automorphism;
  /**
   * @typedef OBIsomorphismMapper::Mappings Automorphisms
   * @brief A group of automorphic permutations.
   * @since 2.3
   */
  typedef OBIsomorphismMapper::Mappings Automorphisms;

  /**
   * Find the automorphisms of a molecule by using an OBIsomorphismMapper. This
   * function is a wrapper for FindAutomorphisms with a functor to store all
   * automorphisms.
   *
   * @param mol The molecule for which to find the automorphisms.
   * @param aut The result will be stored here.
   * @param symmetry_classes The graph invariants to use. See OBGraphSym or use
   * the FindAutomorphisms function that computes these for you below.
   * @param mask A bit vector specifying which atoms to consider. An empty mask
   * will consider all atoms. The bits are indexed from 1 (i.e. OBAtom::GetIdx()).
   * @param maxMemory Maximum memory limit for @p aut. The number of automorphisms
   * for certain graphs can be large. The default is 300MB, consider using a functor
   * to process automorphisms when they are found.
   *
   * @since version 2.3
   */
  OBAPI bool FindAutomorphisms(OBMol *mol, std::vector<OBIsomorphismMapper::Mapping> &aut, const std::vector<unsigned int> &symmetry_classes,
      const OBBitVec &mask = OBBitVec(), std::size_t maxMemory = 3000000);
  /**
   * Find the automorphisms of a molecule by using an OBIsomorphismMapper. This
   * function will first find the graph invariants (i.e. symmetry_classes) using
   * the mask. This function is a wrapper for FindAutomorphisms with a functor to
   * store all automorphisms.
   *
   * @since version 2.3
   */
  OBAPI bool FindAutomorphisms(OBMol *mol, std::vector<OBIsomorphismMapper::Mapping> &aut, const OBBitVec &mask = OBBitVec(),
      std::size_t maxMemory = 3000000);

  /**
   * Find the automorphisms of a molecule by using an OBIsomorphismMapper. This
   * is the main implementation for finding automorphisms and uses an
   * OBIsomorphismMapper::Functor to process found isomorphisms. Wrapper functions
   * are provided to find and store all automorphisms but the number of automorphisms
   * can be large for certain graphs making it not feasible to store all automorphisms
   * in memory (RAM).
   *
   * @see  @ref MapGeneric
   * @since version 2.3
   */
  OBAPI void FindAutomorphisms(OBIsomorphismMapper::Functor &functor, OBMol *mol,
      const std::vector<unsigned int> &symmetry_classes, const OBBitVec &mask = OBBitVec());

  /**
   * @page substructure Substructure Search
   * @since version 2.3
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
   * #include <openbabel/query.h>
   * using namespace OpenBabel;
   *
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
   * @subsection MapGeneric
   * MapGeneric takes a functor object and calls the functor to handle found
   * isomorphisms. This allows for custom mapping results to be obtained by
   * filtering or avoid storing all results. To implement a custom functor,
   * a simple class that inherits OBIsomorphismMapper::Functor and implements
   * the required operator().
   *
   * @code
   * #include <openbabel/isomorphism.h>
   * using namespace OpenBabel;
   *
   * class MyCustomFunctor : public OBIsomorphismMapper::Functor
   * {
   *   private:
   *     // store all mappings in m_data
   *     std::vector<OBIsomorphismMapper::Mapping> &m_data;
   *   public:
   *     MyCustomFunctor(std::vector<OBIsomorphismMapper::Mapping> &data) : m_data(data) {}
   *     bool operator()(OBIsomorphismMapper::Mapping &map)
   *     {
   *       // code to handle found isomorphism...
   *       // examples: - store the mapping
   *       //           - filter mappings
   *       //           - use the found map in some way
   *
   *       m_data.push_back(map);
   *
   *       // continue mapping
   *       return false;
   *     }
   * }
   * @endcode
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

/// @file isomorphism.h
/// @brief OBIsomorphismMapper class for finding isomorphisms.
