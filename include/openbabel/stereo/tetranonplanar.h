/**********************************************************************
  tetranonplanar.h - OBTetraNonPlanarStereo

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
#ifndef OB_TETRANONPLANAR_H
#define OB_TETRANONPLANAR_H

#include <openbabel/stereo/stereo.h>
#include <openbabel/oberror.h>

namespace OpenBabel {

  ///@addtogroup stereo Stereochemistry
  ///@{

  /**
   * @class OBTetraNonPlanarStereo tetranonplanar.h <openbabel/stereo/tetranonplanar.h>
   * @brief Base class for handling and storing non-planar stereochemistry with 4 reference atom ids.
   *
   * @image html tetranonplanar.png
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
   * However, when dealing with tetrahedral stereochemistry, it is often easier
   * to visualize by viewing from/towards one of the reference atoms to/from the
   * center atom. This reduces the 24 possible combinations to 3! = 6.
   *
   * @code
   * 123   321
   * 231   213
   * 312   132
   * @endcode
   *
   * These can be devided in 2 groups: clockwise or anti-clockwise
   *
   * @code
   * clockwise: 123, 231, 312
   * anti-clockwise: 321, 213, 132
   * @endcode
   *
   * Since subclass ConfigType structs accept refs viewing from/towards any atom,
   * it is needed to have some rules for converting.
   *
   * A single permutation of two consecutive elements in a sequence of 3
   * changes the winding. All permutations can be expressed as a combination
   * of consecutive permutations. The number of consecutive permutations can
   * be calculated from the difference in inversions (NumInversions()).
   *
   * If we exchange the from atom with another atom in the sequence, the oddness
   * of the difference in inversions between the 2 sequences is calculated. If this
   * is even, no extra permutation is needed. If this is odd, an extra permutation
   * is needed.
   *
   * Switching between viewing from and viewing towards reverses the winding.
   *
   * Like all stereo classes, errors, warnings or info is reported using OBMessageHandler.
   *
   * @since version 2.3
   */
  class OBAPI OBTetraNonPlanarStereo : public OBStereoBase
  {
    public:
      /**
       * Constructor
       * @param mol The molecule.
       */
      OBTetraNonPlanarStereo(OBMol *mol);
      /**
       * Destructor.
       */
      virtual ~OBTetraNonPlanarStereo();

      /**
       * Convert a @p ConfigType struct from any View/Winding to the
       * desired representation.
       *
       * This is a template method which works on ConfigType structs from
       * OBTetraNonPlanar subclasses. The subclasses can decide what data
       * member are needed to store the stereogenic unit (i.e. 1 atom for
       * tetrahedral, 3 for allene like, ...) and still use this generic
       * method to handle the real stereochemistry.
       *
       * A ConfigType struct should at least have the following data members:
       * @code
       * class SomeNonPlanarStereo : public TetraNonPlanarStereo
       * {
       *   public:
       *     struct Config
       *     {
       *       // constructor(s) are recommended!
       *
       *       // subclass specific stereogenic unit
       *       ...
       *
       *       union {
       *         unsigned long from;
       *         unsigned long towards;
       *       };
       *       OBStereo::Refs refs;
       *       OBStereo::Winding winding;
       *       OBStereo::View view;
       *     };
       * };
       * @endcode
       *
       * @param cfg A ConfigType struct from a OBTetraNonPlanar subclass.
       * @param from_or_towards The desired from/towards reference id (see @p view)
       * @param winding The desired winding.
       * @param view The desired viewing direction.
       *
       * @return The ConfigType struct with desired from/towards, winding and view.
       */
      template <typename ConfigType>
      static ConfigType ToConfig(const ConfigType &cfg, unsigned long from_or_towards,
          OBStereo::Winding winding = OBStereo::Clockwise,
          OBStereo::View view = OBStereo::ViewFrom)
      {
        if (cfg.from == OBStereo::NoRef) {
          obErrorLog.ThrowError(__FUNCTION__,
              "OBTetraNonPlanarStereo::ToConfig : Invalid from in ConfigType struct.", obError);
          return ConfigType();
        }
        if (cfg.refs.size() != 3) {
          obErrorLog.ThrowError(__FUNCTION__,
              "OBTetraNonPlanarStereo::ToConfig : Invalid refs size.", obError);
          return ConfigType();
        }

        // copy the internal refs
        ConfigType result = cfg;
        result.from = from_or_towards;
        result.winding = winding;
        result.view = view;

        // keep track of the permuations by using the oddness
        bool odd = false;

        // find id
        if (cfg.from != from_or_towards) {
          // move id to front and remove it = 1 permutation
          for (int i = 0; i < 3; ++i) {
            if (cfg.refs.at(i) == from_or_towards) {
              result.refs[i] = cfg.from;
              break;
            }
          }
          // 1 permutation perfromed --> odd = true
          odd = true;
        }

        // clockwise <-> anti-clockwise : odd = true
        if (winding == cfg.winding)
          odd = !odd;
        // ViewFrom <-> ViewTowards : odd = true
        if (view == cfg.view)
          odd = !odd;

        // make sure we actually found id
        if (result.refs.size() == 3) {
          if (odd)
            OBStereo::Permutate(result.refs, 1, 2);
          return result;
        }

        obErrorLog.ThrowError(__FUNCTION__,
            "OBTetraNonPlanarStereo::ToConfig : Parameter id not found in internal refs.", obError);
        return result;
      }
      /**
       * Change the winding of the ConfigType struct while maintaining the stereochemistry.
       */
      template <typename ConfigType>
      static void ChangeWinding(ConfigType &cfg)
      {
        cfg.winding = (cfg.winding == OBStereo::Clockwise) ? OBStereo::AntiClockwise : OBStereo::Clockwise;
        OBStereo::Permutate(cfg.refs, 1, 2);
      }
      /**
       * Change the view of the ConfigType struct while maintaining the stereochemistry.
       */
      template <typename ConfigType>
      static void ChangeView(ConfigType &cfg)
      {
        cfg.view = (cfg.view == OBStereo::ViewFrom) ? OBStereo::ViewTowards : OBStereo::ViewFrom;
        OBStereo::Permutate(cfg.refs, 1, 2);
      }
      /**
       * Invert the stereochemistry of the ConfigType struct.
       */
      template <typename ConfigType>
      static void Invert(ConfigType &cfg)
      {
        OBStereo::Permutate(cfg.refs, 1, 2);
      }
 };

  ///@} // addtogroup

}

#endif

/// @file tetranonplanar.h
/// @brief Base class for handling and storing non-planar stereochemistry with 4 reference atom ids.
