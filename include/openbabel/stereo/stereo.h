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
#include <map>
#include <climits> // UINT_MAX

namespace OpenBabel {

  struct OBAPI OBStereo 
  {
    /**
     * Bond directions used by StereoFrom0D to translate to
     * internal CisTransStereo representation
     */
    enum BondDirection { // Values taken from MDL format
      NotStereo =   0,
      UpBond =      1,
      DownBond =    6,
      UnknownDir =  4
    };
    /**
     * The various types of stereochemistry
     *
     */
    enum Type {
      CisTrans            = (1<<0), //!< cis/trans double bond
      ExtendedCisTrans    = (1<<1), //!< allene, biphenyl, ...
      SquarePlanar        = (1<<2), //!< Square-planar stereochemistry
      Tetrahedral         = (1<<3), //!< tetrahedral
      ExtendedTetrahedral = (1<<4), //!< extended tetrahedral
      TrigonalBipyramidal = (1<<5), //!< Trigonal-bipyramidal stereochemistry
      Octahedral          = (1<<6)  //!< Octahedral stereochemistry
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
  class OBAPI OBStereoBase : public OBGenericData
  {
    public:
      /**
       * Constructor. By default, the stereochemistry is specified. Use 
       * SetSpecified(false) for unspecified/unknown stereochemistry.
       *
       * @param mol The molecule.
       */
      OBStereoBase(OBMol *mol) : 
        OBGenericData("StereoData", OBGenericDataType::StereoData, perceived),
        m_mol(mol), m_specified(true)
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
      /**
       * Set wether the stereochemistry is specified. Comparing a specified 
       * OBStereoBase derived class (or it's Config struct) with an unspecified 
       * one, always returns true.
       */
      void SetSpecified(bool specified) { m_specified = specified; }
      /**
       * @return True if the stereochemistry is specified.
       */
      bool IsSpecified() const { return m_specified; }
    private:
      OBMol *m_mol; //!< the parent molecule
      bool m_specified; //!< true if the stereochemistry is specified, false if unknown/unspecified
  };

  class OBTetrahedralStereo;
  class OBCisTransStereo;
  class OBAPI OBStereoFacade
  {
    public:
      /**
       * Constructor with @p mol and @p perceive parameter.
       *
       * @param perceive If true, PerceiveStereo will be called if the 
       * OBMol::HasChiralityPerceived() flag is not set. (default is true)
       */
      OBStereoFacade(OBMol *mol, bool perceive = true) : 
          m_mol(mol), m_init(false), m_perceive(perceive)
      {
      }

      /**
       * Check if atom with @p id is a tetrahedral center. 
       * @return True if the atom with @p id has tetrahedral stereochemistry.
       */
      bool HasTetrahedralStereo(unsigned long atomId);
      /**
       * 
       */
      bool HasCisTransStereo(unsigned long bondId);

      OBTetrahedralStereo* GetTetrahedralStereo(unsigned long atomId);
      OBCisTransStereo* GetCisTransStereo(unsigned long bondId);
    

    private:
      inline void EnsureInit() { if (!m_init) InitMaps(); }
      void InitMaps();

      OBMol *m_mol;
      bool m_init;
      bool m_perceive;
      std::map<unsigned long, OBTetrahedralStereo*> m_tetrahedralMap;
      std::map<unsigned long, OBCisTransStereo*> m_cistransMap;
  };
  class OBBond;
  /**
   * Convert 2D/3D coordinates to OBStereo objects.
   * 
   * @code
   * if (mol->Has3D())
   *   StereoFrom3D(mol, force);
   * else
   *   StereoFrom2D(mol, force);
   * @endcode
   */
  OBAPI void PerceiveStereo(OBMol *mol, bool force = false); 
  /**
   * Convert the 2D depiction of molecule @p mol to OBStereo objects. 
   *
   * First, symmetry analysis taking stereochemistry into account is 
   * performed iterativily (see OBGraphSym). Next the 2D coordinates, 
   * OBBond::Wedge, OBBond::Hash, OBBond::WedgeOrHash and OBBond::CisOrTrans 
   * are used to construct OBStereoBase derived objects to store the 
   * stereochemistry.
   *
     @verbatim
     Reference:
     [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the 
     Unambiguous Identification of the Stereochemical Characteristics of 
     Compounds During Their Registration in Databases. Molecules 2000, 6,
     915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
     @endverbatim
   *
   * @param mol The molecule containing 2D coordinates.
   * @param force Force to run the perception even if the results are cached.
   */
  OBAPI void StereoFrom2D(OBMol *mol, bool force = false);
  /**
   * Convert the 3D coordinates of molecule @p mol to OBStereo objects.
   *
   * @param mol The molecule containing 3D coordinates.
   * @param force Force to run the perception even if the results are cached.
   */
  OBAPI void StereoFrom3D(OBMol *mol, bool force = false);
  /**
   * Add missing OBStereo objects. Unlike StereoFrom2D and StereoFrom2D, this
   * method only adds objects for previously unidentified objects since we 
   * don't want to loose any information. 
   *
   * @param mol The molecule.
   * @param updown A map of OBStereo::BondDirection for cis/trans bonds
   */
  OBAPI void StereoFrom0D(OBMol *mol,
      const std::map<OBBond*, enum OBStereo::BondDirection> *updown = NULL);

  /**
   * Get a vector with all OBTetrahedralStereo objects for the molecule. This 
   * function is used by StereoFrom3D with the @p addToMol parameter is set 
   * to true. 
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBTetrahedralStereo objects will be added 
   * to the molecule using OBBase::SetData().
   */
  OBAPI std::vector<OBTetrahedralStereo*> TetrahedralFrom3D(OBMol *mol, 
      const std::vector<unsigned int> symClasses, bool addToMol = true);
  /**
   * Get a vector with all OBTetrahedralStereo objects for the molecule. This 
   * function is used by StereoFrom2D with the @p addToMol parameter is set 
   * to true. 
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBTetrahedralStereo objects will be added 
   * to the molecule using OBBase::SetData().
   */
  OBAPI std::vector<OBTetrahedralStereo*> TetrahedralFrom2D(OBMol *mol, 
      const std::vector<unsigned int> symClasses, bool addToMol = true);
  /**
   * Get a vector with all OBTetrahedralStereo objects for the molecule. This 
   * function is used by StereoFrom0D with the @p addToMol parameter is set 
   * to true. 
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBTetrahedralStereo objects will be added 
   * to the molecule using OBBase::SetData().
   */
  OBAPI std::vector<OBTetrahedralStereo*> TetrahedralFrom0D(OBMol *mol, 
      const std::vector<unsigned int> symClasses, bool addToMol = true);
  /**
   * Get a vector with all OBCisTransStereo objects for the molecule. This 
   * function is used by StereoFrom3D with the @p addToMol parameter is set 
   * to true. 
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBCisTransStereo objects will be added 
   * to the molecule using OBBase::SetData().
   */
  OBAPI std::vector<OBCisTransStereo*> CisTransFrom3D(OBMol *mol, 
      const std::vector<unsigned int> &symClasses, bool addToMol = true);
  /**
   * Get a vector with all OBCisTransStereo objects for the molecule. This 
   * function is used by StereoFrom2D with the @p addToMol parameter is set 
   * to true. 
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBCisTransStereo objects will be added 
   * to the molecule using OBBase::SetData().
   */
  OBAPI std::vector<OBCisTransStereo*> CisTransFrom2D(OBMol *mol, 
      const std::vector<unsigned int> &symClasses, bool addToMol = true);
  /**
   * Get a vector with all OBCisTransStereo objects for the molecule. This 
   * function is used by StereoFrom0D with the @p addToMol parameter is set 
   * to true. 
   *
   * This function is also used for symmetry analysis to handle cases where 
   * there are two atoms in the same symmetry class that don't have the same 
   * stereochemistry. In this situation, the @p addToMol parameter is set to 
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param symClasses The symmetry classes to use.
   * @param addToMol If true, the OBCisTransStereo objects will be added 
   * to the molecule using OBBase::SetData().
   */
  OBAPI std::vector<OBCisTransStereo*> CisTransFrom0D(OBMol *mol, 
      const std::vector<unsigned int> &symClasses,
      std::map<OBBond*, OBStereo::BondDirection> *updown = NULL,
      bool addToMol = true);
 
  /**
   * Find all tetrahedral centers using the symmetry classes in the 
   * @p symClasses map.
   *
   * @return The atom ids for all tetrahedral centers.
   */
  OBAPI std::vector<unsigned long> FindTetrahedralAtoms(OBMol *mol, 
      const std::vector<unsigned int> &symClasses);
  /**
   * Find all cis/trans bonds using the symmetry classes in the 
   * @p symClasses map.
   *
   * @return The bond ids for all cis/trans bonds.
   */
  OBAPI std::vector<unsigned long> FindCisTransBonds(OBMol *mol, 
      const std::vector<unsigned int> &symClasses);
  
  void CisTransFromUpDown(OBMol &mol,
      const std::vector<unsigned int> &ctbonds,
      const std::map<OBBond*, OBStereo::BondDirection> &updown);
}

#endif
