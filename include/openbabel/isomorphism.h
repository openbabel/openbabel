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
   */
  class OBAPI OBIsomorphismMapper
  {
    public:
      typedef std::map<unsigned int,unsigned int> Mapping;
      typedef std::vector<Mapping> Mappings;

      OBIsomorphismMapper(OBQuery *query) : m_query(query) {}

      static OBIsomorphismMapper* GetInstance(OBQuery *query, const std::string &algorithm = std::string("bfs"));

      virtual Mapping MapFirst(const OBMol *queried) = 0;
      virtual Mappings MapUnique(const OBMol *queried) = 0;
      virtual Mappings MapAll(const OBMol *queried) = 0;
    protected:
      OBQuery *m_query;
  };


  /**
   * Find the automorphisms of a molecule by using a OBIsomorphismMapper.
   * @since version 2.3
   */
  OBAPI OBIsomorphismMapper::Mappings FindAutomorphisms(OBMol *mol);


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
