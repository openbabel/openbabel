/**********************************************************************
  stereo.h - OBStereo & OBStereoBase

  Copyright (C) 2009-2010 by Tim Vandermeersch

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
#ifndef OB_STEREO_H
#define OB_STEREO_H

#include <openbabel/base.h> // OBGenericData
#include <openbabel/isomorphism.h> // Automorphisms
#include <vector>
#include <map>
#include <set>
#include <climits> // UINT_MAX

namespace OpenBabel {

  ///@addtogroup stereo Stereochemistry
  ///@{

  /**
   * @class OBStereo stereo.h <openbabel/stereo/stereo.h>
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
   * @since version 2.3
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
     * Bond directions used by StereoFrom0D() to translate to
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
     *
     * @image html SPshapes.png
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
      ViewTowards = 2 //!< view from center atom towards the atom (id parameter)
    };

    /**
     * Windings used by OBTetraNonPlanar subclasses for
     * setting/getting reference ids.
     * @sa OBTetraNonPlanar
     */
    enum Winding {
      Clockwise = 1,     //!< Clockwise winding
      AntiClockwise = 2, //!< AntiClockwise winding (or CounterClockwise)
      UnknownWinding = 3 //!< The configuration is specified as unknown (squiggly line in depiction)
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
    static void Permutate(Refs &refs, unsigned int i, unsigned int j);
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
    static Refs Permutated(const Refs &refs, unsigned int i, unsigned int j);
    //@}

  };

  /**
   * @struct OBStereoUnit stereo.h <openbabel/stereo/stereo.h>
   * @brief Struct representing a single stereogenic unit.
   * @since 2.3
   */
  struct OBStereoUnit
  {
    /**
     * Default constructor creating an invalid OBStereoUnit. The unit can be
     * made valid by setting the type, id and para data members.
     */
    OBStereoUnit() : type(static_cast<OBStereo::Type>(0)), id(OBStereo::NoRef), para(false)
    {
    }

    /**
     * Constructor specifying all data to create a valid OBStereoUnit.
     */
    OBStereoUnit(OBStereo::Type _type, unsigned long _id, bool _para = false) :
        type(_type), id(_id), para(_para)
    {
    }

    OBStereo::Type type; //!< the type for this stereogenic unit
    unsigned long id; //! the atom/bond (depends on type) unique id
    bool para; //! para- (=ressemble) or true-stereocenter
  };
  /**
   * @brief A single set of OBStereoUnit objects.
   *
   * This type can be used to represent all stereogenic units in a molecule and
   * is used as return type of FinStereogenicUnits(). This set is also the input
   * for many functions requiring this information (e.g. StereoFrom2D, ...).
   */
  typedef std::vector<OBStereoUnit> OBStereoUnitSet;
  /**
   * @brief A set of sets of OBStereoUnit objects.
   *
   * This type is used for cases where there is some relationship between
   * individual OBStereoUnit objects.
   */
  typedef std::vector<OBStereoUnitSet> OBStereoUnitSetOfSets;


  // fwd decl
  class OBMol;
  /**
   * @class OBStereoBase stereo.h <openbabel/stereo/stereo.h>
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
   * @since version 2.3
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
      virtual ~OBStereoBase() { m_mol = nullptr; }

      ///@name Geniric (for all OBStereo::Type) stereochemistry
      //@{
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
      //@}
    private:
      OBMol *m_mol; //!< The parent molecule.
      bool m_specified; //!< True if the stereochemistry is specified, false if unknown/unspecified.
  };

  // fwd decl
  class OBTetrahedralStereo;
  class OBCisTransStereo;
  class OBSquarePlanarStereo;
  /**
   * @class OBStereoFacade stereo.h <openbabel/stereo/stereo.h>
   * @brief Facade to simplify retrieval of OBStereoBase derived objects.
   *
   * The OBStereoFacade helps with retrieving OBStereoBase derived objects
   * (i.e. OBTetrahedralStereo, OBCisTransStereo, ...) from an OBMol. This
   * is done by iterating over all OBGenericData objects with data type
   * OBGenericDataType::StereoData and checking the OBStereo::Type using
   * OBStereoBase::GetType.
   *
   * @sa OBStereo OBStereoBase
   * @since version 2.3
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

      ///@name Tetrahedral stereochemistry
      ///@{
      /**
       * Get the number of tetrahedral stereocenters.
       */
      unsigned int NumTetrahedralStereo();
      /**
       * Get all the OBTetrahedralStereo objects.
       */
      std::vector<OBTetrahedralStereo*> GetAllTetrahedralStereo();
      /**
       * Check if atom with @p id is a tetrahedral center.
       * @return True if the atom with @p id has tetrahedral stereochemistry.
       */
      bool HasTetrahedralStereo(unsigned long atomId);
      /**
       * Get the OBTetrahedralStereo object with @p atomId as center. This
       * function returns 0 if there is no OBTetrahedralStereo object found
       * with the specified center.
       */
      OBTetrahedralStereo* GetTetrahedralStereo(unsigned long atomId);
      ///@}

      ///@name Cis/Trans stereochemistry
      ///@{
      /**
       * Get the number of cis/trans stereocenters.
       */
      unsigned int NumCisTransStereo();
      /**
       * Get all the OBCisTransStereo objects.
       */
      std::vector<OBCisTransStereo*> GetAllCisTransStereo();
      /**
       * Check if bond with @p id is a stereogenic cis/trans double bond.
       * @return True if the bond with @p id has cis/trans stereochemistry.
       */
      bool HasCisTransStereo(unsigned long bondId);
      /**
       * Get the OBTetrahedralStereo object with @p bondId as double bond.
       * This function returns 0 if there is no OBCisTransStereo object found
       * with the specified bond.
       */
      OBCisTransStereo* GetCisTransStereo(unsigned long bondId);
      ///@}

      ///@name SquarePlanar stereochemistry
      ///@{
      /**
       * Get the number of square-planar stereocenters.
       */
      unsigned int NumSquarePlanarStereo();
      /**
       * Get all the OBSquarePlanarStereo objects.
       */
      std::vector<OBSquarePlanarStereo*> GetAllSquarePlanarStereo();
      /**
       * Check if atom with @p id is a stereogenic square-planar atom.
       * @return True if the atom with @p id has square-planar stereochemistry.
       */
      bool HasSquarePlanarStereo(unsigned long atomId);
      /**
       * Get the OBSquarePlanarStereo object with @p atomId as center. This
       * function returns 0 if there is no OBSquarePlanarStereo object found
       * with the specified center.
       */
      OBSquarePlanarStereo* GetSquarePlanarStereo(unsigned long atomId);
      ///@}

      template<int StereoType>
      bool HasStereo(unsigned long id);
      template<typename T>
      T* GetStereo(unsigned long id);


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
      std::map<unsigned long, OBSquarePlanarStereo*> m_squarePlanarMap;
  };

  // fwd decl
  class OBBond;
  ///@name High level functions
  ///@{
  /**
   * Convert 0D/2D/3D coordinates to OBStereo objects. The right function will
   * be selected based on the molecule's dimensionality
   * (i.e. OBMol::GetDimension()).
   *
   * @sa StereoFrom3D StereoFrom2D StereoFrom0D
   * @since version 2.3
   */
  OBAPI void PerceiveStereo(OBMol *mol, bool force = false);
  /**
   * Convert the 2D depiction of molecule @p mol to OBStereo objects.
   * This function makes use of the lower level functions
   * TetrahedralFrom2D(), CisTransFrom2D(), SquarePlanarFrom2D(), ...
   *
   * First, symmetry analysis taking stereochemistry into account is
   * performed iteratively (see OBGraphSym). Next the 2D coordinates,
   * OBBond::Wedge, OBBond::Hash, OBBond::WedgeOrHash and OBBond::CisOrTrans
   * are used to construct OBStereoBase derived objects to store the
   * stereochemistry. These objects will be added to @p mol.
   *
   * Unless perception is forced, this function does nothing if stereochemistry
   * has already been perceived (i.e. OBMol::HasChiralityPerceived()). Before
   * doing the actual perception, any data of the OBGenericDataType::StereoData
   * type will be deleted.
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
   * @param updown A map of OBStereo::BondDirection for cis/trans bonds
   * @param force Force to run the perception even if the results are cached.
   *
   * @sa StereoFrom3D StereoFrom0D PerceiveStereo
   * @since version 2.3
   */
  OBAPI void StereoFrom2D(OBMol *mol,
    std::map<OBBond*, enum OBStereo::BondDirection> *updown = nullptr, bool force = false);
  /**
   * Convert the 3D coordinates of molecule @p mol to OBStereo objects. This
   * function makes use of the lower level functions TetrahedralFrom3D(),
   * CisTransFrom3D(), SquarePlanarFrom3D(), ...
   *
   * Unless perception is forced, this function does nothing if stereochemistry
   * has already been perceived (i.e. OBMol::HasChiralityPerceived()). Before
   * doing the actual perception, any data of the OBGenericDataType::StereoData
   * type will be deleted.
   *
   * @param mol The molecule containing 3D coordinates.
   * @param force Force to run the perception even if the results are cached.
   *
   * @sa StereoFrom3D StereoFrom0D PerceiveStereo
   * @since version 2.3
   */
  OBAPI void StereoFrom3D(OBMol *mol, bool force = false);
  /**
   * Add missing OBStereo objects. Unlike StereoFrom3D() and StereoFrom2D(), this
   * method only adds objects for previously unidentified objects since we
   * don't want to loose any information. The Config::specified flag for the
   * newly added structs is always set to false.
   *
   * For example, a smiles is read which has two tetrahedral centers. Only one has
   * stereochemisrty specified using a '@' character. StereoFrom0D() will detect the
   * second tetrahedral atom and add an OBTetrahedralStereo object to the molecule.
   *
   * @param mol The molecule.
   *
   * @sa StereoFrom3D StereoFrom2D PerceiveStereo
   * @since version 2.3
   */
  OBAPI void StereoFrom0D(OBMol *mol);
  ///@}

  ///@name Low level functions
  ///@{
  /**
   * Get a vector with all OBTetrahedralStereo objects for the molecule. This
   * function is used by StereoFrom3D() with the @p addToMol parameter is set
   * to true.
   *
   * The algorithm to convert the 3D coordinates to OBTetrahedralStereo object
   * uses the sign of the volume described by the 4 center atom neighbors. Given
   * 4 points \f$a\f$, \f$b\f$, \f$c\f$ and \f$d\f$, the signed volume \f$S_v\f$
   * is defined as:
   *
     \f[ S_v = \left| \begin{array}{ccc}
     x_b - x_a & y_b - y_a & z_b - z_a \\
     x_c - x_a & y_c - y_a & z_c - z_a \\
     x_d - x_a & y_d - y_a & z_d - z_a
     \end{array} \right| \f]
   *
   * The sign of \f$S_v\f$ changes when any of the points cross the plane defined
   * by the other 3 points. To make this less abstract one could say that
   * a change of sign is equal to inverting the tetrahedral stereochemistry.
   *
   * In case there are only 3 neighbor atoms for the tetrahedral center, the
   * center atom itself is used as 4th point. This only changes the magnitude
   * and not the sign of \f$S_v\f$ because the center atom is still on the same
   * side of the plane.
   *
   * This function is also used for symmetry analysis to handle cases where
   * there are two atoms in the same symmetry class that don't have the same
   * stereochemistry. In this situation, the @p addToMol parameter is set to
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param stereoUnits The stereogenic units.
   * @param addToMol If true, the OBTetrahedralStereo objects will be added
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom3D FindStereogenicUnits
   * @since version 2.3
   */
  OBAPI std::vector<OBTetrahedralStereo*> TetrahedralFrom3D(OBMol *mol,
      const OBStereoUnitSet &stereoUnits, bool addToMol = true);
  /**
   * Get a vector with all OBTetrahedralStereo objects for the molecule. This
   * function is used by StereoFrom2D() with the @p addToMol parameter is set
   * to true.
   *
   * The algorithm to convert the 2D coordinates and bond properties
   * (i.e. OBBond::Wedge, OBBond::Hash, OBBond::WedgeOrHash and OBBond::CisOrTrans)
   * uses the sign of a triangle. Given 3 points \f$a\f$, \f$b\f$ and \f$c\f$, the
   * sign of the trianle \f$S_t\f$ is defined as:
   *
     \f[ S_t = (x_a - x_c) (y_b - y_c) - (y_a - y_c) (x_b - x_c) \f]
   *
   * This is equation 6 from on the referenced web page. The 3 points used
   * to calculate the triangle sign always remain in the same plane (i.e. z = 0).
   * The actual meaning of \f$S_t\f$ (i.e. assignment of OBStereo::Winding) depends
   * on the 4th atom. When the atom is in front of the plane, the sign should be
   * changed to have the same absolute meaning for an atom behind the plane and the
   * same triangle. It is important to note that none of the z coordinates is ever
   * changed, the molecule always stays 2D (unlike methods which set a pseudo-z
   * coordinate).
   *
   * @todo document bond property interpretation!
   *
   * This function is also used for symmetry analysis to handle cases where
   * there are two atoms in the same symmetry class that don't have the same
   * stereochemistry. In this situation, the @p addToMol parameter is set to
   * false and the returned objects will need to be deleted explicitly.
   *
     @verbatim
     Reference:
     [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the
     Unambiguous Identification of the Stereochemical Characteristics of
     Compounds During Their Registration in Databases. Molecules 2000, 6,
     915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
     @endverbatim
   *
   * @param mol The molecule.
   * @param stereoUnits The stereogenic units.
   * @param addToMol If true, the OBTetrahedralStereo objects will be added
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom2D FindStereogenicUnits
   * @since version 2.3
   */
  OBAPI std::vector<OBTetrahedralStereo*> TetrahedralFrom2D(OBMol *mol,
      const OBStereoUnitSet &stereoUnits, bool addToMol = true);
  /**
   * Get a vector with all OBTetrahedralStereo objects for the molecule. This
   * function is used by StereoFrom0D() with the @p addToMol parameter is set
   * to true. There is no algorithm used here, all specified flags will be
   * set to false.
   *
   * This function is also used for symmetry analysis to handle cases where
   * there are two atoms in the same symmetry class that don't have the same
   * stereochemistry. In this situation, the @p addToMol parameter is set to
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param stereoUnits The stereogenic units.
   * @param addToMol If true, the OBTetrahedralStereo objects will be added
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom0D FindStereogenicUnits
   * @since version 2.3
   */
  OBAPI std::vector<OBTetrahedralStereo*> TetrahedralFrom0D(OBMol *mol,
      const OBStereoUnitSet &stereoUnits, bool addToMol = true);

  /**
   * Get a vector with all OBCisTransStereo objects for the molecule. This
   * function is used by StereoFrom3D() with the @p addToMol parameter is set
   * to true.
   *
   * The algorithm to convert the 3D coordinates to OBCisTransStereo objects
   * considers the signed distance between the attached atoms and the plane
   * through the double bond at right angles to the plane of the attached
   * atoms. Bonds on the same side (cis) will share the same sign for the
   * signed distance.
   *
   * Missing atom coordinates (OBStereo::ImplicitRef) and their bond
   * vectors will be computed if needed.
   *
   @verbatim
         0      3     Get signed distance of 0 and 2 to the plane
          \    /      that goes through the double bond and is at
           C==C       right angles to the stereo bonds.
          /    \
         1      2     If the two signed distances have the same sign
                      then they are cis; if not, then trans.
   @endverbatim
   *
   * This function is also used for symmetry analysis to handle cases where
   * there are two atoms in the same symmetry class that don't have the same
   * stereochemistry. In this situation, the @p addToMol parameter is set to
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param stereoUnits The stereogenic units.
   * @param addToMol If true, the OBCisTransStereo objects will be added
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom3D FindStereogenicUnits
   * @since version 2.3
   */
  OBAPI std::vector<OBCisTransStereo*> CisTransFrom3D(OBMol *mol,
      const OBStereoUnitSet &stereoUnits, bool addToMol = true);
  /**
   * Get a vector with all OBCisTransStereo objects for the molecule. This
   * function is used by StereoFrom2D() with the @p addToMol parameter is set
   * to true.
   *
   * This function is also used for symmetry analysis to handle cases where
   * there are two atoms in the same symmetry class that don't have the same
   * stereochemistry. In this situation, the @p addToMol parameter is set to
   * false and the returned objects will need to be deleted explicitly.
   *
   * The algorithm for converting the 2D coordinates uses the same triangle
   * sign as TetrahedralFrom2D(). Depending on sign of 2 triangles, the right
   * OBStereo::Shape is selected.
   @verbatim
      0      3
       \    /        2 triangles: 0-1-b & 2-3-a
        a==b    -->  same sign: U
       /    \        opposite sign: Z
      1      2
   @endverbatim
   *
   * @param mol The molecule.
   * @param stereoUnits The stereogenic units.
   * @param updown A map of OBStereo::BondDirection for cis/trans bonds
   * @param addToMol If true, the OBCisTransStereo objects will be added
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom2D FindStereogenicUnits
   * @since version 2.3
   */
  OBAPI std::vector<OBCisTransStereo*> CisTransFrom2D(OBMol *mol,
      const OBStereoUnitSet &stereoUnits,
      const std::map<OBBond*, enum OBStereo::BondDirection> *updown = nullptr, bool addToMol = true);
  /**
   * Convert a molecule's OBTetrahedralStereo objects to a series of hash or
   * wedge bonds. Note that the molecule itself is not modified; the result
   * is returned in the maps @p updown and @p from, which indicate
   * the origin and direction of each hash or wedge bond.
   *
   * When converting, the following guidelines are followed when trying to
   * find the best candidate bond to set up/down for each OBTetrahedralStereo
   * object:
   * -# Should not already be set
   * -# Should not be connected to a 2nd tet center
   *    (this is acceptable in theory as the wedge is only at one end, but
   *     in practice it may cause confusion and thus we avoid it)
   * -# Preferably is not in a cycle
   * -# Preferably is a terminal H
   *
   * If no bond can be found that matches rules 1 and 2 (and in theory this is possible)
   * then an error message is logged and the function returns false. (If you find an
   * example where this occurs, please file a bug.)
   *
   * @param mol The molecule.
   * @param updown A map of OBStereo::BondDirection for each hash/wedge bond
   * @param from A map of OBStereo::Ref indicating the origin of each hash/wedge bond
   * @return True or False depending on whether the conversion was successful
   * @since version 2.3
   */
  OBAPI bool TetStereoToWedgeHash(OBMol &mol,
      std::map<OBBond*, enum OBStereo::BondDirection> &updown,
      std::map<OBBond*, OBStereo::Ref> &from);
  /**
   * Return a set of double bonds corresponding to the OBCisTransStereo objects
   * for which the stereochemistry is undefined.
   *
   * Note that this functions just iterates over the existing OBCisTransStereo
   * objects - it does not try to identify new ones.
   *
   * @param mol The molecule
   * @return A set of bonds with unspecified cis/trans stereochemistry
   * @since version 2.3
   */
  OBAPI std::set<OBBond*> GetUnspecifiedCisTrans(OBMol& mol);
  /**
   * Convert any reference to @p atomId in a stereo object to an OBStereo::ImplicitRef.
   * This function is called from OBMol::DeleteHydrogens()
   * (via OBMol::DeleteHydrogen()) to remove any explicit references to a
   * hydrogen atom that has been deleted. However, the code is not specific
   * to hydrogen atoms and could be used for other atoms.
   *
   * @param mol The molecule
   * @param atomId The Id of the atom to be converted to an OBStereo::ImplicitRef
   * @since version 2.3
   */
  OBAPI void StereoRefToImplicit(OBMol& mol, OBStereo::Ref atomId);
  /**
   * Convert any reference to an OBStereo::ImplicitRef attached to @p centerId
   * in a stereo object to an explicit reference to @p newId.
   * This function is called from OBMol::AddHydrogens() and
   * OBMol::AddHydrogen() to convert any implicit references to a
   * hydrogen atom that has just been added. However, the code is not specific
   * to hydrogen atoms and could be used for other atoms.
   *
   * @param mol The molecule
   * @param centerId The Id of the atom to which the new explicit atom is attached
   * @param newId The Id of the atom which was previously an OBStereo::ImplicitRef
   * @since version 2.4
   */
  OBAPI void ImplicitRefToStereo(OBMol& mol, OBStereo::Ref centerId, OBStereo::Ref newId);
  /**
   * Get a vector with all OBCisTransStereo objects for the molecule. This
   * function is used by StereoFrom0D() with the @p addToMol parameter is set
   * to true. There is no algorithm used here, all specified flags will be
   * set to false.
   *
   * This function is also used for symmetry analysis to handle cases where
   * there are two atoms in the same symmetry class that don't have the same
   * stereochemistry. In this situation, the @p addToMol parameter is set to
   * false and the returned objects will need to be deleted explicitly.
   *
   * @param mol The molecule.
   * @param stereoUnits The stereogenic units.
   * @param addToMol If true, the OBCisTransStereo objects will be added
   * to the molecule using OBBase::SetData().
   *
   * @sa StereoFrom0D FindStereogenicUnits
   * @since version 2.3
   */
  OBAPI std::vector<OBCisTransStereo*> CisTransFrom0D(OBMol *mol,
      const OBStereoUnitSet &stereoUnits,
      bool addToMol = true);
  ///@}


  ///@name Stereogenic unit identification
  ///@{
  /**
   * Find the stereogenic units in a molecule using a set of rules.<sup>1</sup>
   *
   * The potential stereocenters are identified first. A potential tetrahedral
   * stereogenic atom is any atom meeting the following criteria:
   *
   * - sp3 hybridization
   * - at least 3 "heavy" neighbors
   *
   * Nitrogen is treated as a special case since the barrier of inversion is
   * low in many cases making the atom non-stereogenic. Only bridge-head
   * nitrogen atoms (i.e. nitrogen has 3 neighbors in rings) will be
   * considered stereogenic.
   *
   * Potential stereogenic double bonds are identified using another set of
   * simple criteria:
   *
   * - must be a double bond
   * - must not be in a ring
   * - both begin and end atom should have at least one single bond
   *
   * True stereocenters (i.e. stereocenters with topologically different
   * ligands) are identified first. For tetrahedral stereocenters, true
   * stereocenters will have 4 different neighbor atom symmetry classes
   * and this can be expressed using T1234 to classify these stereocenters.
   * For stereogenic bonds, a similar classification C12 can be used but
   * both begin and end atom have their own classification and the bond
   * is only a true stereocenter if both atoms are C12.
   *
   * Para stereocenters are all stereocenters where there are at least two
   * equivalent neighbor atom symmetry classes. These are T1123, T1112, T1111
   * and T1122 for tetrahedral stereocenters and C11 for double bonds. To
   * determine which of the remaining potential stereocenters really are
   * stereocenters, a set of rules is used.<sup>1</sup>
   *
   * Rule 1 is applied recusively:
   *
   * All rings are merged "mergedRings". A merged ring is simply a fragment consisting
   * of all atoms of a ring system (bridged, spiro, adjacent, ...). If two rings in the
   * SSSR set share an atom, they are merged.
   *
   * Each merged must at least have two para-stereocenters (or 1 true + 1 para) in order
   * for the para-stereocenter to be valid. This is repeated until no new stereocenters
   * are identified.
   *
   * rule 1a for double bonds:
   * - bond atom in ring has two identical symmetry classes for it's neighbor atoms (-> para)
   * - other bond atom:
   *   - has two different symmetry classes for it's neighbours -> new stereocenter
   *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
   *
   * rule 1b for tetracoord atoms:
   * - at least two neighbour symmetry classes are the same (-> para)
   * - other pair:
   *   - has two different symmetry classes for it's neighbours -> new stereocenter
   *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
   *
   * Rules 2 and 3 are applied sequential (i.e. only once).
   *
   * Rule 2a for tetracoordinate carbon:
   * - 1 or 2 pair identical ligands
   * - each ligand contains at least 1 true-stereocenter or 2 para-stereocenters (from rule 1)
   *
   * Rule 2b for tetracoordinate carbon:
   * - 3 or 4 identical ligands with at least
   *   - 2 true-stereocenters
   *   - 2 separate assemblies of para-stereocenters (from rule 1)
   *
   * Rule 3 for double bonds:
   * - 1 or 2 pair identical ligands (on begin and end atom)
   * - each pair contains at least 1 true-stereocenter or 2 para-stereocenters (from rule 1)
   *
   *
   * @verbatim
     Reference:
     [1] M. Razinger, K. Balasubramanian, M. Perdih, M. E. Munk, Stereoisomer
     Generation in Computer-Enhanced Structure Elucidation, J. Chem. Inf.
     Comput. Sci. 1993, 33, 812-825
     @endverbatim
   */
  OBAPI OBStereoUnitSet FindStereogenicUnits(OBMol *mol,
      const std::vector<unsigned int> &symClasses);
  /**
   * @brief Find the stereogenic units in a molecule making use of the automorphisms.
   *
   * The potential stereocenters are identified first. A potential tetrahedral
   * stereogenic atom is any atom meeting the following criteria:
   *
   * - sp3 hybridization
   * - at least 3 "heavy" neighbors
   *
   * Nitrogen is treated as a special case since the barrier of inversion is
   * low in many cases making the atom non-stereogenic. Only bridge-head
   * nitrogen atoms (i.e. nitrogen has 3 neighbors in rings) will be
   * considered stereogenic.
   *
   * Potential stereogenic double bonds are identified using another set of
   * simple criteria:
   *
   * - must be a double bond
   * - must not be in a ring
   * - both begin and end atom should have at least one single bond
   *
   * Once the potential stereocenters are found, the automorphisms are the key
   * to identifying real stereogenic units. Automorphisms can be seen as
   * permutations that permutate a graph back to the same graph. Such a
   * permutation can only exchange atoms with the same symmetry class and it
   * follows that the use of automorphisms takes symmetry into account. The
   * definitions below use a concept where the automorphisms cause inversions
   * of configuration to potential stereocenters. Such an inversion occurs
   * whenever an automorphism exchanges two equivalent (i.e. with the same
   * symmetry class) neighbor atoms attached to the potential stereogenic unit.
   *
   * @par Definition for tetrahedral stereocenters:
   * A potential stereocenter really is a stereocenter if there exists no automorphic
   * permutation causing an inversion of the configuration of only the potential
   * stereogenic unit under consideration.
   * If there exists at least one automorphic permutation causing an inversion of
   * the configuration, then the potential stereogenic center can be a stereogenic
   * center if the number of topologically equivalent neighbors (ligands) of the
   * potential stereogenic center is less than or equal to the number of different
   * configurations of these ligands.<sup>1</sup>
   *
   * The actual number of configurations needed for the ligands depends on the
   * classification (i.e. T1234, T1123, ...) of the stereo center. These classes
   * reflect the symmetry classes of the neighbor atoms of the center.
   *
   * - T1123: 1 true stereocenter OR 2 para stereocenters
   * - T1122: 1 true stereocenter OR 2 para stereocenters (for both)
   * - T1112: 2 true stereocenters OR 2 para stereocenter assemblies
   * - T1111: 2 true stereocenters OR 2 para stereocenter assemblies
   *
   * @par Definition for double bond stereocenters:
   * A potential stereogenic double bond really is a stereogenic bond if there
   * exists no automorphic permutation causing an inversion of the configuration
   * of only the potential stereogenic unit under consideration. The bond can still
   * be a stereogenic bond if there exists such an automorphism when the number of
   * configurations of the pair of topologically equivalent geminal ligands, which
   * are exchanged by the automorphism, is greater than or equal to two (i.e. the
   * number of topologically equivalent geminal ligands.<sup>1</sup>
   *
   * For stereogenic bonds, there is only one case but both begin and end atom
   * have to be checked.
   *
   * - C11: 1 true stereocenter OR 1 para stereocenter
   *
   * These criteria are analogous to the rules from the Razinger paper on
   * stereoisomer generation. Since the existence of stereocenters can depend
   * on the existence of other stereocenters (in the ligands), the stereocenters
   * are found by iterating until no new stereocenters are found.
   *
   * @verbatim
     Reference:
     [1] M. Perdih, M. Razinger, Stereochemistry and Sequence Rules:
     A Proposal for Modification of Cahn-Ingold-Prelog System,
     Tetrahedron: Asymmetry, 1994, Vol. 5, No. 5, 835-861
     @endverbatim
   */
  OBAPI OBStereoUnitSet FindStereogenicUnits(OBMol *mol,
      const std::vector<unsigned int> &symClasses,
      const Automorphisms &automorphisms);
  ///@}

  /**
   * @page Stereochemistry
   * @section overview Overview of classes
   *
   * There are many molecules which contain stereogenic elements. However,
   * certain cases (i.e. tetrahedral, cis/trans) are more common than others
   * (i.e. allene, biphenyl, octrahedral, ...). For the common stereogenic
   * units, classes are provided. The inheritance of these classes resembles
   * the way they are split into groups.
   *
   * - OBStereoBase
   *   - OBTetraNonPlanarStereo
   *     - OBTetrahedralStereo
   *     - OBExtendedTetrahedralStereo
   *   - OBTetraPlanarStereo
   *     - OBCisTransStereo
   *     - OBExtendedCisTransStereo
   *     - OBSquarePlanarStereo
   *   - OBAxialStereo
   *     - OBTrigonalBipyrimidalStereo
   *     - OBOctahedralStereo
   *
   * @image html tetranonplanar.png
   * @image html tetraplanar.png
   *
   * All specific classes (i.e. OBTetrahedralStereo, ...) have embedded Config
   * structs which define the actual stereochemistry. All these Config structs
   * use OBStereo::Ref values to reference or uniquely identify atoms. Make sure
   * to read about OBStereo::Ref and the related functions (in OBStereo). OBStereo
   * is also a placeholder for various enums with predefined values for parameters
   * etc. These enums are used throughout the different stereo classes but having
   * these enums in a single location makes it easier to remember. When working
   * with stereo classes, you normally don't need to use any of the parent classes
   * directly. Only OBStereo and the specific class are needed.
   *
   * @section usage Basic usage
   *
   * The OBStereoFacade hides the complexity of working with stereochemistry. When
   * using openbabel as a library, this is by far the easiest way to access
   * stereochemistry information.
   * The header for the specific OBStereo::Type type is all you need to include.
   * These are:
   * - @em openbabel/stereo/tetrahedral.h
   * - @em openbabel/stereo/cistrans.h
   * - @em openbabel/stereo/squareplanar.h
   *
   * All these headers also include @em openbabel/stereo/stereo.h providing
   * declarations for OBStereo & OBStereoFacade.
   *
     @code
     #include <iostream>
     #include <openbabel/mol.h>
     #include <openbabel/obconversion.h>

     #include <openbabel/stereo/tetrahedral.h>

     using namespace OpenBabel;

     int main()
     {
       OBMol mol;
       OBConversion conv;
       conv.SetInFormat("smi");
       conv.ReadString(&mol, "C[C@H](Cl)Br");

       OBStereoFacade facade(&mol);

       FOR_ATOMS_OF_MOL(atom, mol) {
         if (facade.HasTetrahedralStereo(atom->GetId()))
           std::cout << facade.GetTetrahedralStereo(atom->GetId()) << std::endl;
       }
     }
     @endcode
   *
   * All specific stereo classes and their embedded Config struct have an
   * operator<< function which allows them to be used with std::ostream objects
   * (e.g. std::cout, std::err, ...). These functions are often useful when
   * debugging code.
   *
   * @section details Details on implementation
   *
   * The detection of stereogenic units start with symmetry analysis. However, a
   * complete symmetry analysis also needs to take stereochemistry into account.
   * In practice, this means stereochemistry will be found iteratively. At each
   * iteration, the current atom symmetry classes are used to identify stereogenic
   * units. The details about how the symmetry classes are used depends on the type
   * (OBStereo::Type) of stereogenic unit. For tetrahedral centers, having 3 heavy
   * atom neighbors with different symmetry classes or 4 neighbors with different
   * symmetry classes means the atom is chiral. See FindStereogenicUnits() for
   * details.
   *
   * After identifying the stereogenic units, Config structs with all the
   * information on the spacial arrangement of the groups still have to be
   * created. This involves interpreting various ways to represent
   * stereochemisrty:
   *
   * - 3D coordinates: StereoFrom3D()
   * - 2D coordinates: StereoFrom2D()
   * - 0D coordinates: StereoFrom0D()
   *
   * Both StereoFrom3D() and StereoFrom2D() delete all existing stereochemistry objects
   * before adding new ones. For molecules with 3D coordinates, it is evident that
   * all information is specified by the coordinates itself. However, if a file format
   * uses stereo parity flags, Config structs must be constructed using lower level
   * functions and StereoFrom3D() should not be called. In these cases information
   * could be lost by calling StereoFrom3D() after reading the file (the stereo flag might have
   * indicated the stereochemistry was unspecified or the flag might not match the
   * coordinates). In the case of 2D molecules, the coordinates together with bond
   * properties (OBBond::Hash, OBBond::Wedge, OBBond::WedgeOrHash and
   * OBBond::CisOrTrans) define the stereochemistry. Again, lower level functions
   * can be used when stereo flags need to be used.
   *
   * StereoFrom0D() works slightly different than 3D/2D. Here, deleting the
   * stereochemistry would always result in lost information. Instead StereoFrom0D()
   * only adds new objects for stereogenic units which were previously not found.
   * For example, a smiles is read which has two tetrahedral centers. Only one has
   * stereochemistry specified using a '@' character. StereoFrom0D() will detect the
   * second tetrahedral atom and add an OBTetrahedralStereo object to the molecule.
   * The Config::specified flag for the newly added structs is always set to false.
   *
   * Assuming the format code has correctly set the molecule dimensions (OBMol::GetDimesions),
   * PerceiveStereo() will automatically select the correct function to call.
   * When StereoFrom3D(), StereoFrom2D() or StereoFrom0D() are not used, make sure to always
   * set OBMol::HasChiralityPerceived() before returning from the format's ReadMolecule().
   *
   *
   * @section formats Guidelines for formats
   * @subsection input Reading files
   *
   * - Read the section above
   * - The MDL format (mdlformat.cpp) is a good example for 2D/3D formats with or
   *   without parity flags.
   * - The SMILES format (smilesformat.cpp) is a good example for 0D formats.
   *
   * @subsection output Writing files
   *
   * For many file formats no additional code is needed. For example, if a 3D format
   * doesn't require stereo parity flags, writing the coordinates is enough. For 2D
   * file formats it will often suffice to write the coordinates and bond properties.
   * If parity flags are needed, the OBStereoFacade class can be used to retrieve the
   * objects for all types of stereochemistry supported by the file format.
   *
   *
   *
   *
   *
   * @since version 2.3
   */

  ///@}  addtogroup
}

#endif

//! \file stereo.h
//! \brief Process molecular stereochemistry information.
