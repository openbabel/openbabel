/**********************************************************************
  tetraplanar.h - OBTetraPlanarStereo

  Copyright (C) 2009 by Tim Vandermeersch

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
#ifndef OB_TETRAPLANAR_H
#define OB_TETRAPLANAR_H

#include "stereo.h"
#include <algorithm> // std::rotate

namespace OpenBabel {

  /**
   * @class OBTetraPlanarStereo tetraplanar.h <openbabel/stereo/tetraplanar.h>
   * @brief Base class for handling and storing planar stereochemistry with 4 reference atoms.
   *
   * @image html tetraplanar.png
   *
   * @section Combinations
   * The four reference ids can be treated like a sequence of 4 numbers. Each element can
   * only occur once. This means there are 4! = 24 combinations.
   *
   * These are the 24 possible combinations or permutations.
   * @code
   * 1234   2134   3124   4123
   * 1243   2143   3142   4132
   * 1324   2314   3214   4213
   * 1342   2341   3241   4231
   * 1423   2413   3412   4312
   * 1432   2431   3421   4321
   * @endcode
   *
   * Based on which reference ids are on opposite sides (180°), it is possible
   * to divide these 24 combinations in three sets. For this it is also needed
   * to make use of a shape to map the reference id indexes on the points in
   * the plane. Without these shapes or a fixed meaning, these sequences have
   * no meaning. The use of these shapes (U, Z & 4) is illustrated in the figure
   * below:
   * @image html SPshapes.png
   *
   * In the figure, it can be seen the OBStereo::ShapeU implies the 1st
   * reference id in the sequence is opposite to the 3th and the 2nd is opposite
   * to the 4th. The same can be done for OBStereo::ShapeZ and OBStereo::Shape
   *
   * Grouped by opposite reference ids using OBStereo::Shape4
   * @code
   * 1-2, 3-4 : 1234, 1243, 2134, 2143, 3412, 3421, 4312, 4321
   * 1-3, 2-4 : 1324, 1242, 2413, 2431, 3124, 3142, 4213, 4231
   * 1-4, 2-3 : 1423, 1432, 2314, 2341, 3214, 3241, 4123, 4132
   * @endcode
   *
   * Internally the reference ids are stored in a U shape. The OBTetraPlanarStereo::ToInternal()
   * and OBTetraPlanarStereo::ToShape() methods convert between the internally used U shape
   * and the other shapes.
   *
   * Like all stereo classes, errors, warnings or info is reported using OBMessageHandler.
   *
   * @note U shaped ordering can also be considered circular and is the only shape
   * that can be rotated lexicographically.
   *
   * @since version 2.3
   */

  class OBAPI OBTetraPlanarStereo : public OBStereoBase
  {
    public:
      OBTetraPlanarStereo(OBMol *mol);
      virtual ~OBTetraPlanarStereo();

      template <typename ConfigType>
      static ConfigType ToConfig(const ConfigType &cfg, unsigned long start,
          OBStereo::Shape shape = OBStereo::ShapeU)
      {
        ConfigType result = cfg;
        result.shape = shape;

        // convert from U/Z/4 to U shape
        switch (cfg.shape) {
          case OBStereo::ShapeU:
            break;
          case OBStereo::ShapeZ:
            OBStereo::Permutate(result.refs, 2, 3); // convert to U shape
            break;
          case OBStereo::Shape4:
            OBStereo::Permutate(result.refs, 1, 2); // convert to U shape
            break;
        }


        // since refs are U shaped we can rotate the refs lexicographically
        for (int i = 0; i < 4; ++i) {
          std::rotate(result.refs.begin(), result.refs.begin() + 1, result.refs.end());
          // start at refs[0]?
          if (result.refs.at(0) == start)
            break;
        }

        // convert from U to desired U/Z/4
        // (don't change refs[0]!)
        switch (shape) {
          case OBStereo::ShapeU:
            break;
          case OBStereo::ShapeZ:
            OBStereo::Permutate(result.refs, 2, 3); // convert to Z shape
            break;
          case OBStereo::Shape4:
            OBStereo::Permutate(result.refs, 1, 2); // convert to 4 shape
            break;
        }

        return result;
      }
  
  };

}

#endif

/// @file tetraplanar.h
/// @brief Base class for handling and storing planar stereochemistry with 4 reference atoms.
