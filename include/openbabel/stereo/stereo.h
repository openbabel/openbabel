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

  /**
   * @brief Placeholder for enums & Ref/Refs related functions.
   *
   * The OBStereo struct contains a number of enums with predefined values. 
   * These are OBStereo::BondDirection, OBStereo::Type, OBStereo::Shape, 
   * OBStereo::View, OBStereo::Winding. There are enums 
   * which only apply to certain types of stereochemistry but having them 
   * in 1 place makes it easier to remember. 
   *
   * The OBStereo struct also contains typedefs and functions which 
   * are crucial to fully understand how to use the OBStereoBase derived 
   * classes (i.e. OBTetrahedralStereo, OBCisTransStereo, ...). Ref variables
   * and Refs lists are a way to uniquely reference atoms in the molecule. In
   * most cases these Ref variables are the same as the unique atom ids 
   * (OBAtom::GetId). However, 2 special cases are provided:
   *
   * - OBStereo::NoRef: An initial value for Ref variables. The constructors of 
   *   the various Config structs set all refs to NoRef (Refs lists will remain 
   *   empty). This value is considered invalid when comparing stereochemistry.
   * - OBStereo::ImplicitRef: Can be used to replace implicit hydrogen ids. 
   *   Even with explicit hydrogens, it is still valid to replace the Ref of 
   *   the hydrogen with ImplicitRef. This flexibility also applies to 
   *   comparing stereochemistry. It is possible to compare stereochemistry
   *   (e.g. OBTetrahedral::Config::operator==) between a Config struct with 
   *   an ImplicitRef and a Config struct where the Ref is the atom id of the 
   *   explicit hydrogen. See actual documentation for details about operator==.
   *
   * There are utility functions which make it easier to handle Refs. The most
   * frequently used one is OBStereo::MakeRefs to create lists containing 3 or
   * 4 Ref values. Formats and library use normally doesn't need the other 
   * functions.
   *
   * @sa OBStereoBase OBStereoFacade
   */
  struct OBAPI OBStereo 
  {
    /**
     * The various types of stereochemistry
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
     * Bond directions used by StereoFrom0D to translate to
     * internal CisTransStereo representation.
     */
    enum BondDirection { // Values taken from MDL format
      NotStereo =   0,
      UpBond =      1,
      DownBond =    6,
      UnknownDir =  4
    };

    /**
     * Shapes used by OBTetraPlanarStereo subclasses for 
     * setting/getting reference ids.
     * @sa OBTetraPlanarStereo
     */
    enum Shape {
      ShapeU = 1,
      ShapeZ = 2,
      Shape4 = 3
    };

    /**
     * Views used by OBTetraNonPlanarStereo subclasses for
     * setting/getting reference ids.
     * @sa OBTetraNonPlanarStereo
     */
    enum View
    {
      ViewFrom = 1, //!< view from the atom (id parameter) towards the center atom
      ViewTowards = 2, //!< view from center atom towards the atom (id paramater)
    };

    /**
     * Windings used by OBTetraNonPlanar subclasses for 
     * setting/getting reference ids.
     * @sa OBTetraNonPlanar
     */
    enum Winding {
      Clockwise = 1,     //!< Clockwise winding
      AntiClockwise = 2  //!< AntiClockiwe winding (or CounterClockwise
    };

    ///@name Ref & Refs types 
    //@{
    /**
     * All stereo classes work with variables of the type Ref to uniquely
     * identify atoms. In most cases these Ref variables are the same as the
     * unique atom ids. 
     *
     * @sa OBAtom::GetId() OBStereo
     */
    typedef unsigned long Ref;
    /**
     * Special case Ref values.
     */
    enum {
      NoRef = UINT_MAX,       //!< No Ref set (invalid Ref)
      ImplicitRef = UINT_MAX - 1  //!< Implicit Ref (i.e. hydrogen, N lone pair, ...).
    };
    /**
     * A list (std::vector) of Ref variables.
     */
    typedef std::vector<Ref> Refs;
    /**
     * Iterator (std::iterator) for a Refs list.
     */
    typedef Refs::iterator RefIter;
    /**
     * Iterator (std::iterator) for a const Refs list.
     */ 
    typedef Refs::const_iterator ConstRefIter;
    //@}

    ///@name Refs utility functions
    //@{
    /**
     * Create a Refs list filled with @p ref1, @p ref2, @p ref3 & @p ref4.
     * @p ref4 is not added to the returned Refs if it is equal to NoRef.
     *
     * @return A Refs list containing the specified Ref values. 
     */
    static Refs MakeRefs(Ref ref1, Ref ref2, Ref ref3, Ref ref4 = NoRef)
    {
      Refs refs(3);
      refs[0] = ref1;
      refs[1] = ref2;
      refs[2] = ref3;
      if (ref4 != NoRef)
        refs.push_back(ref4);
      return refs;
    }
    /**
     * Check if @p refs1 and @p refs2 contain the same Ref values regardless 
     * of their order.
     *
     * @code
     * OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), OBStereo::MakeRefs(2, 3, 1)) // true
     * OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), OBStereo::MakeRefs(3, 2, 1)) // true
     * OBStereo::ContainsSameRefs(OBStereo::MakeRefs(1, 2, 3), OBStereo::MakeRefs(3, 4, 1)) // false
     * @endcode
     *
     * @return True if @p refs1 and @p refs2 contain the same Ref values.
     */
    static bool ContainsSameRefs(const Refs &refs1, const Refs &refs2);
    /**
     * @return True if @p refs contains @p ref.
     */
    static bool ContainsRef(const Refs &refs, unsigned long ref);
    //@}

    ///@name Low-level functions used by implementation.
    //@{
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
    static int NumInversions(const Refs &refs);
    /**
     * Permutate element @p i with @p j in @p refs.
     *
     * @param refs The sequence with N elements to permutate.
     * @param i Element i (0...N-1) will be mutated to j and vice versa.
     * @param j Element j (0...N-1) will be mutated to i and vice versa.
     *
     * @note This method does nothing if i equals j.
     */
    static void Permutate(Refs &refs, int i, int j);
    /**
     * Get @p refs with element @p i and @p j permutated.
     *
     * @param refs The sequence with N elements to permutate.
     * @param i Element @p i (0...N-1) will be mutated to @p j and vice versa.
     * @param j Element @p j (0...N-1) will be mutated to @p i and vice versa.
     *
     * @return @p refs with elements @p i and @p j permutated.
     *
     * @note This method does nothing if @p i equals @p j.
     */
    static Refs Permutated(const Refs &refs, int i, int j);
    //@}
 
  };

  // fwd decl
  class OBMol;
  /**
   * @brief Base class for all stereochemistry classes.
   *
   * All stereochemistry classes are derived from OBStereoBase. This class
   * inherits from OBGenericData which allows the objects to be stored in
   * the molecule. The attribute (OBGenericData::GetAttribute) is set to
   * "StereoData" and the data type is OBGenericDataType::StereoData. The 
   * pure virtual OBStereoBase::GetType function must be implemented by
   * derived classes to return a type defined in OBStereo::Type.
   *
   * Use the OBStereoFacade for easy access to the derived classes. 
   *
   * OBStereoBase keeps track of the OBMol object. This must always be
   * a valid (not 0 or deleted) pointer and can only be set using the 
   * constructor. Subclasses can use this to get more information on bonding 
   * for example. Finally, OBStereoBase also keeps track of the specified
   * flag. By default, this is always set to true.
   *
   * @sa OBStereo OBStereoFacade
   */
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
      /**
       * Destructor.
       */
      virtual ~OBStereoBase() { m_mol = 0; }
      /**
       * Get the molecule. This can be used by subclasses when more
       * information is needed (e.g. OBCisTransStereo::GetCisRef, ...).
       */
      OBMol* GetMolecule() const { return m_mol; }
      /**
       * Reimplemented by subclasses to return the type defined in OBStereo::Type.
       */
      virtual OBStereo::Type GetType() const = 0;
      /**
       * Set whether the stereochemistry is specified. Comparing a specified 
       * OBStereoBase derived class (or it's Config struct) with an unspecified 
       * one, always returns true.
       */
      void SetSpecified(bool specified) { m_specified = specified; }
      /**
       * @return True if the stereochemistry is specified.
       */
      bool IsSpecified() const { return m_specified; }
    private:
      OBMol *m_mol; //!< The parent molecule.
      bool m_specified; //!< True if the stereochemistry is specified, false if unknown/unspecified.
  };

  // fwd decl
  class OBTetrahedralStereo;
  class OBCisTransStereo;
  /**
   * @brief Facade to simplify retrieval of OBStereoBase derived objects.
   *
   * The OBStereoFacade helps with retrieving OBStereoBase derived objects 
   * (i.e. OBTetrahedralStereo, OBCisTransStereo, ...) from an OBMol. This 
   * is done by iterating over all OBGenericData objects with data type 
   * OBGenericDataType::StereoData and checking the OBStereo::Type using 
   * OBStereoBase::GetType.
   *
   * @sa OBStereo OBStereoBase
   */
  class OBAPI OBStereoFacade
  {
    public:
      /**
       * Constructor with @p mol and @p perceive parameter.
       *
       * @param mol The molecule.
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
       * Check if bond with @p id is a stereogenic cis/trans double bond. 
       * @return True if the bond with @p id has cis/trans stereochemistry.
       */
      bool HasCisTransStereo(unsigned long bondId);

      /**
       * Get the OBTetrahedralStereo object with @p atomId as center. This 
       * function returns 0 if there is no OBTetrahedralStereo object found 
       * with the specified center.
       */
      OBTetrahedralStereo* GetTetrahedralStereo(unsigned long atomId);
      /**
       * Get the OBTetrahedralStereo object with @p atomId as center. This 
       * function returns 0 if there is no OBTetrahedralStereo object found 
       * with the specified center.
       */
      OBCisTransStereo* GetCisTransStereo(unsigned long bondId);
    

    private:
      /**
       * Ensure the maps are initialized and initialize them only once.
       */
      inline void EnsureInit() { if (!m_init) InitMaps(); }
      /**
       * Initialize @p m_tetrahedralMap and m_cistransMap to contain the 
       * data objects. If @p m_perceive is true and chirality isn't perceived
       * yet, PerceiveStereo will be called.
       */
      void InitMaps();

      OBMol *m_mol;
      bool m_init;
      bool m_perceive;
      std::map<unsigned long, OBTetrahedralStereo*> m_tetrahedralMap;
      std::map<unsigned long, OBCisTransStereo*> m_cistransMap;
  };

  // fwd decl
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
  OBAPI void StereoFrom2D(OBMol *mol, bool tetfrom0D = false, bool force = false);
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
      std::map<OBBond*, enum OBStereo::BondDirection> *updown = NULL);

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
 
  /**
   * Create and fill OBCisTransStereo objects using the specified
   * @p ctbonds (bond ids) and map containing directions 
   * (OBStereo::BondDirection). This function is intended to be used
   * by 0D formats. The OBCisTransStereo objects will be stored in the 
   * OBMol.
   *
   * @param mol The molecule.
   * @param ctbonds std::vector containing bond ids for cis/trans bonds
   * @param updown std::map containing up to four bond directions per bond 
   *        id in @p ctbonds.
   */
  OBAPI void CisTransFromUpDown(OBMol *mol,
      const std::vector<unsigned long> &ctbonds,
      std::map<OBBond*, OBStereo::BondDirection> *updown);
}

#endif
