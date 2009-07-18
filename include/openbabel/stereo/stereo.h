/**********************************************************************
  stereo.h - OBStereo & OBStereoBase

  Copyright (C) 2009 by Tim Vandermeersch
 
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
#ifndef OB_STEREO_H
#define OB_STEREO_H

#include <openbabel/base.h> // OBGenericData
#include <vector>
#include <climits> // UINT_MAX

namespace OpenBabel {

  struct OBAPI OBStereo 
  {
    /**
     * The various types of stereochemistry
     *
     */
    enum Type {
      None                = (1<<0), //!< no stereochemistry
      Unknown             = (1<<2), //!< unknown
      Unspecified         = (1<<3), //!< not specified
      CisTrans            = (1<<4), //!< cis/trans double bond
      ExtendedCisTrans    = (1<<5), //!< allene, biphenyl, ...
      SquarePlanar        = (1<<6), //!< Square-planar stereochemistry
      TetraPlanar         = CisTrans | ExtendedCisTrans | SquarePlanar, //!< planar configurations of four atoms
      Tetrahedral         = (1<<7), //!< tetrahedral
      ExtendedTetrahedral = (1<<8), //!< extended tetrahedral
      TetreNonPlanar      = Tetrahedral | ExtendedTetrahedral
      /*
      TrigonalBipyramidal = 4, //!< Trigonal-bipyramidal stereochemistry
      Octahedral          = 5, //!< Octahedral stereochemistry
      DoubleBond          = 6, //!< Double bond stereochemistry
      */
    };

    /**
     * Shapes used by OBTetraPlanarStereo subclasses for 
     * setting/getting reference ids.
     */
    enum Shape {
      ShapeU = 1,
      ShapeZ = 2,
      Shape4 = 3
    };

    /**
     * Views used by OBTetraNonPlanarStereo subclasses for
     * setting/getting reference ids.
     */
    enum View
    {
      ViewFrom = 1, //!< view from the atom (id parameter) towards the center atom
      ViewTowards = 2, //!< view from center atom towards the atom (id paramater)
    };

    /**
     * Windings used by OBTetraNonPlanar subclasses for 
     * setting/getting reference ids.
     */
    enum Winding {
      Clockwise = 1,     //!< Clockwise winding
      AntiClockwise = 2  //!< AntiClockiwe winding (or CounterClockwise
    };

    /**
     * Some useful predefined ids. 
     */
    enum {
      NoId = UINT_MAX,       //!< no id
      ImplicitId = UINT_MAX - 1  //!< implicit id (hydrogen, N lone pair, ..)
    };

    typedef std::vector<unsigned long> Refs;
    typedef std::vector<unsigned long>::iterator RefIter;
    typedef std::vector<unsigned long>::const_iterator ConstRefIter;

    /**
     * Create a std::vector<unsigned long> filled with @p id1, @p id2, @p id3 & @p id4.
     */
    static Refs MakeRefs(unsigned long id1, unsigned long id2,
        unsigned long id3, unsigned long id4 = NoId)
    {
      Refs refs(3);
      refs[0] = id1;
      refs[1] = id2;
      refs[2] = id3;
      if (id4 != NoId)
        refs.push_back(id4);
      return refs;
    }

    /**
     * Check if @p refs1 and @p refs2 contain the same ids regardless 
     * of their order.
     *
     * @code
     * OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), OBStereo::MakeRefs(2, 3, 1)) // true
     * OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), OBStereo::MakeRefs(3, 2, 1)) // true
     * OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), OBStereo::MakeRefs(3, 4, 1)) // false
     * @endcode
     */
    static bool ContainsSameRefs(const Refs &refs1, const Refs &refs2);

    /**
     * @return True if @p refs contains @p id.
     */
    static bool ContainsRef(const Refs &refs, unsigned long id);

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
    static int NumInversions(const OBStereo::Refs &refs);

    /**
     * Permutate element @p i with @p j in @p refs.
     *
     * @param refs The sequence with N elements to permutate.
     * @param i Element i (0...N-1) will be mutated to j and vice versa.
     * @param j Element j (0...N-1) will be mutated to i and vice versa.
     *
     * @note This method does nothing if i equals j.
     */
    static void Permutate(OBStereo::Refs &refs, int i, int j);
    /**
     * Get @p refs with element @p i and @p j permutated.
     *
     * @param refs The sequence with N elements to permutate.
     * @param i Element i (0...N-1) will be mutated to j and vice versa.
     * @param j Element j (0...N-1) will be mutated to i and vice versa.
     *
     * @return @p refs with @i and @j permutated.
     *
     * @note This method does nothing if i equals j.
     */
    static OBStereo::Refs Permutated(const OBStereo::Refs &refs, int i, int j);
 
  };

  class OBMol;
  class OBStereoBase : public OBGenericData
  {
    public:
      OBStereoBase(OBMol *mol) : 
        OBGenericData("StereoData", OBGenericDataType::StereoData, perceived),
        m_mol(mol) 
      {
      }
      virtual ~OBStereoBase() { m_mol = 0; }
      /**
       * Get the molecule. This can be used by subclasses when more
       * information is needed (e.g. OBCisTransStereo::GetCisRef, ...).
       */
      OBMol* GetMolecule() const { return m_mol; }
      /**
       * Reimplemented by subclasses to return type.
       */
      virtual OBStereo::Type GetType() const = 0;
    private:
      OBMol *m_mol; //!< the parent molecule
  };

  OBAPI void StereoFrom2D(OBMol *);
  OBAPI void StereoFrom3D(OBMol *);

}

#endif
