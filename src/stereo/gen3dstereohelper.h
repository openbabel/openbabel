/**********************************************************************
  gen3dstereohelper.h - Helper class for 3D coordinate generation

  Copyright (C) 2020 by Tim Vandermeersch

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
#ifndef OB_GEN3DSTEREOHELPER_H
#define OB_GEN3DSTEREOHELPER_H

#include <openbabel/babelconfig.h>
#include <string>
#include <vector>

namespace OpenBabel {

  class OBMol;

  /**
   * @brief Helper class for 3D coordinate generation.
   *
   * This class can be used to check stereochemistry when generating 3D
   * coordinates. This class also keeps track of unspecified stereochemistry
   * to ensure this information is not lost when calling StereoFrom3D().
   */
  class OBAPI OBGen3DStereoHelper
  {
    public:
      /**
       * @brief Store stereochemical information for later comparison.
       */
      void Setup(OBMol *mol);
      /**
       * @brief Check the stereochemistry.
       *
       * This function will perceive stereochemistry from 3D and compare this
       * with the stereochemistry that was stored when Setup() was called.
       *
       * @return True if the stereochemistry is correct.
       */
      bool Check(OBMol *mol);
    private:
      std::string m_inputSmiles;
      std::vector<unsigned long> m_unspecifiedTetrahedral;
      std::vector<unsigned long> m_unspecifiedCisTrans;
  };

} // namespace OpenBabel

#endif
