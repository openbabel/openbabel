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



}


#endif
