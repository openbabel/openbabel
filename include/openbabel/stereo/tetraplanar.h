#ifndef OB_TETRAPLANAR_H
#define OB_TETRAPLANAR_H

#include "stereo.h"

namespace OpenBabel {

  /**
   * @class OBTetraPlanarStereo
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
   * Bases on which reference ids are on opposite sides (180°), it is possible
   * to devide these 24 combinations in three sets. For this it is also needed
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
   */

  class OBTetraPlanarStereo : public OBStereoBase
  {
    public:
      OBTetraPlanarStereo(OBMol *mol);
      virtual ~OBTetraPlanarStereo();

      /**
       * Subclasses must implement the SetRefs method.
       */
      virtual void SetRefs(const std::vector<unsigned long> &refs, 
        OBStereo::Shape shape = OBStereo::ShapeU) = 0;
      /**
       * Subclasses must implement the GetRefs method.
       */
      virtual std::vector<unsigned long> GetRefs(OBStereo::Shape shape = OBStereo::ShapeU) const = 0;
      /**
       * Subclasses must implement the Compare method to check 
       * if @p refs match the stored stereochemistry.
       */
      virtual bool Compare(const std::vector<unsigned long> &refs, 
          OBStereo::Shape shape) const = 0;
      /**
       * Convert a sequence of reference ids from U, Z or 4 shape to 
       * internal U shape.
       * @note this method does nothing if a U shape is given as input.
       */
      static std::vector<unsigned long> ToInternal(const std::vector<unsigned long> &refs, 
          OBStereo::Shape shape);
      /**
       * Convert a sequence of reference ids from internal U shape 
       * to U, Z or 4 shape.
       * @note this method does nothing if a U shape is given as input.
       */
      static std::vector<unsigned long> ToShape(const std::vector<unsigned long> &refs, 
          OBStereo::Shape shape);
  };

}

#endif
