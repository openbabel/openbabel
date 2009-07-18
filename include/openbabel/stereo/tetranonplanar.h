#ifndef OB_TETRANONPLANAR_H
#define OB_TETRANONPLANAR_H

#include <openbabel/stereo/stereo.h>

namespace OpenBabel {

  /**
   * @class OBTetraNonPlanarStereo
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
   * Since SetRefs and GetRefs accept refs viewing from/towards any atom in 
   * the sequence, it is needed to have some rules for converting.
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
   *
   */
  class OBTetraNonPlanarStereo : public OBStereoBase
  {
    public:
      OBTetraNonPlanarStereo(OBMol *mol);
      virtual ~OBTetraNonPlanarStereo();

      /**
       * Subclasses must implement the SetRefs method.
       */
      virtual void SetRefs(const std::vector<unsigned long> &refs, unsigned long id,
          OBStereo::Winding winding = OBStereo::Clockwise, 
          OBStereo::View view = OBStereo::ViewFrom) = 0;
      /**
       * Subclasses must implement the GetRefs method.
       */
      virtual std::vector<unsigned long> GetRefs(unsigned long id,
          OBStereo::Winding winding = OBStereo::Clockwise, 
          OBStereo::View view = OBStereo::ViewFrom) const = 0;
      /**
       * Subclasses must implement the Compare function to check 
       * if @p refs match the stored stereochemistry.
       */
      virtual bool Compare(const std::vector<unsigned long> &refs, unsigned long id,
          OBStereo::Winding winding, OBStereo::View view) const = 0;
      /**
       * Convert a sequence of reference ids from any View/Winding to 
       * internal clockwise representation.
       */
      static std::vector<unsigned long> ToInternal(const std::vector<unsigned long> &refs, 
          unsigned long id, OBStereo::Winding winding, OBStereo::View view);
      /**
       * Convert a sequence of reference ids from internal clockwise to any View/Winding.
       */
      static std::vector<unsigned long> ToView(const std::vector<unsigned long> &refs, 
          unsigned long id, OBStereo::Winding winding, OBStereo::View view);
      /**
       * Compute the inversion vector for @p refs and return the sum of it's 
       * elements. The ith element in the inversion vector is the number of 
       * element to the right of element i with a lower value.
       *
       * The number of inversions is the same as the number of interchanges
       * of consecutive elements. 
       *
       * When working with 3 refs from a tetrahedral configuration:
       * @code
       * permutation   inversion vector    sum
       * -------------------------------------
       * 123           0 0 0               0 (even) -> clockwise
       * 132           0 1 0               1 (odd)  -> anti-clockwise
       * 213           1 0 0               1 (odd)  -> anti-clockwise
       * 231           1 1 0               2 (even) -> clockwise
       * 312           2 0 0               2 (even) -> clockwise
       * 321           2 1 0               3 (odd)  -> anti-clockwise
       * @endcode
       */
      static int NumInversions(const std::vector<unsigned long> &refs);
      /**
       * @param refs The sequence with N elements to permutate.
       * @param i Element i (0...N-1) will be mutated to j and vice versa.
       * @param j Element j (0...N-1) will be mutated to i and vice versa.
       *
       * @note This method does nothing if i equals j.
       */
      static std::vector<unsigned long> Permutate(const std::vector<unsigned long> &refs, 
          int i, int j);
  };

}

#endif
