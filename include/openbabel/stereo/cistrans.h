/**********************************************************************
  cistrans.h - Class for handling and storing cis/trans stereochemistry.

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
#ifndef OB_CISTRANS_H
#define OB_CISTRANS_H

#include <openbabel/stereo/tetraplanar.h>
#include <vector>

namespace OpenBabel {

///@addtogroup stereo Stereochemistry
///@{
/**
 * @class OBCisTransStereo cistrans.h <openbabel/stereo/cistrans.h>
 * @brief Class for handling and storing cis/trans stereochemistry.
 *
 * @image html cistrans.png
 *
 * The OBCisTransStereo class is used to represent cis/trans stereochemistry.
 * Like all OBTetraPlanarStereo subclasses, it uses the OBStereo::Shape parameters
 * to set/get the OBStereo::Ref values in the Config struct. However, since the
 * orientation of the double bond matters, all methods in the "query methods"
 * section check the bonding by actually checking the atoms in the molecule
 * provided through OBStereoBase's constructor. If OBMol::GetAtomById(id) returns
 * 0 for a single id, it will be considered a deleted or implicit hydrogen if
 * the valences confirm this.
 *
 * The use of OBStereo::Shape is illustarted in the image below.
 * @image html SPshapes.png
 *
 * @code
 * //
 * // 3       6        3       6
 * //  \     /         |       |
 * //   0===1      =   | 0   1 |
 * //  /     \         |       |
 * // 4       5        4-------5
 * //
 * OBCisTransStereo ct;
 * ct.SetCenters(0, 1);
 * ct.SetRefs(OBStereo::MakeRefs(3, 4, 5, 6), OBStereo::ShapeU);
 *
 * ct.IsTrans(3, 5); // true
 * ct.IsTrans(3, 4); // false
 * ct.IsTrans(3, 6); // false
 *
 * ct.IsCis(3, 6); // true
 * ct.IsCis(3, 4); // false
 * ct.IsCis(3, 5); // false
 * @endcode
 *
 * @since version 2.3
 */
class OBAPI OBCisTransStereo : public OBTetraPlanarStereo
{
  public:
    /**
     * \struct Config cistrans.h <openbabel/stereo/cistrans.h>
     * \brief Stereochemical configuration for double-bond cis/trans stereochemistry
     *
     * The config struct represents the stereochemistry in a well defined way.
     * For cis/trans stereo bonds, the following data members define the spacial
     * arrengement of the atoms.
     *
     * - OBStereo::Ref @p begin: The begin atom for the double bond.
     * - OBStereo::Ref @p end: The end atom for the double bond.
     * - OBStereo::Refs @p refs: The 4 atoms connected to the double bond.
     * - OBStereo::Shape @p shape: The shape formed by the @p refs by connecting them
     *   in the same order as they occur in @p refs.
     *
     * @image html cistrans.png
     * @image html SPshapes.png
     *
     * Only @p begin and @p end are specific for OBCisTransStereo::Config. The other
     * data members occur in all OBTetraPlanarStereo derived classes.
     */
#ifndef SWIG
    struct OBAPI Config
    {
      /**
       * Default constructor. Initializes @p begin and @p end to OBStereo::NoRef
       * and @p shape to OBStereo::ShapeU.
       */
      Config() : begin(OBStereo::NoRef), end(OBStereo::NoRef), shape(OBStereo::ShapeU),
          specified(true)
      {  }
      /**
       * Constructor with all parameters.
       *
       * @param _begin The double bond begin atom id.
       * @param _end The double bond end atom id.
       * @param _refs The 4 reference ids.
       * @param _shape The shape for the 4 reference ids.
       */
      Config(unsigned long _begin, unsigned long _end, const OBStereo::Refs &_refs,
          OBStereo::Shape _shape = OBStereo::ShapeU) : begin(_begin), end(_end),
          refs(_refs), shape(_shape), specified(true)
      {  }
      /**
       * Equal to operator. Comparing OBCisTransStereo::Config structs
       * is done using the information stored in the struct's data members
       * (i.e. begin, end, refs and shape).
       *
       * There are a number of cases resuling in false being returned:
       * - @p begin and @p end don't match (is checked using the 2 combinations)
       * - One of the Refs lists does not contain 4 elements.
       * - 2 or more OBStereo::ImplicitRef values in a single Config struct
       * - (The two @p refs don't share a single common element)
       *
       * In the simplest case where both @p refs contain exactly the same elements
       * (OBStereo::ContainsSameRefs()), coould include OBStereo::ImplicitRef), both Config
       * struct are normalized to OBStereo::ShapeU starting with the same element.
       * After this normalization, there are two possible orientations to overlay the
       * shape on the double bond. From the illustration below, it can be seen only
       * @p refs[2] has to be checked in order to conclude both Config structs
       * have the same stereochemistry.
       *
         @verbatim
         1      4    1      4    1------4
          \    /     |      |           |
           C==C      |      |           |
          /    \     |      |           |
         2      3    2------3    2------3

                     1 2 3 4     1 2 3 4
                     |   |       |   |      <- in any case, refs[0] & refs[2] remain unchanged
                     1 2 3 4     1 4 3 2
        @endverbatim
       *
       * When comparing a Config struct with explicit hydrogen(s) to one with
       * implicit hydrogen(s), both @p refs are also normalized to OBStereo::ShapeU
       * starting with the same common element. This shared element cannot be
       * OBStereo::ImplicitRef. Depending on the position of the OBStereo::ImplicitRef
       * element(s) in the @p refs, 3 cases are possible:
       *
        @verbatim

         refs[2] != OBStereo::ImplicitId:

           (analog to the case above where they contained the same elements )

           1 2 3 4
           |   |      <- refs[0] & refs[2] remain unchanged
           1 H 3 H

         else:

           1 2 3 4
           |     |    <- refs[0] & refs[3] remain unchanged
           1 H H 4

           1 2 3 4
           | |        <- refs[0] & refs[1] remain unchanged
           1 2 H H
        @endverbatim
       *
       * In each case, the orientation of the U shape is also defined since
       * there can be only one OBStereo::ImplicitRef for each side of the
       * double bond.
       *
       * @return True if both Config structs represent the stereochemistry.
       */
      bool operator==(const Config &other) const;
      /**
       * Not equal to operator. This is the inverse of the Equal to operator==.
       *
       * @return True if the two Config structs represent a different stereochemistry.
       */
      bool operator!=(const Config &other) const
      {
        return !(*this == other);
      }

      /**
       * @name Data members defining stereochemistry.
       * @{
       */
      unsigned long begin, end; //<! The double bond begin and end ids.
      OBStereo::Refs refs; //!< The 4 reference ids.
      OBStereo::Shape shape; //!< The shape of the 4 reference ids.
      bool specified; //!< True if the stereochemistry is specified. When false, the described
                      //!< special orientation is only accidental (i.e. unspecified).
      //@}
    };
#endif
    /**
     * Constructor.
     */
    OBCisTransStereo(OBMol *mol);
    /**
     * Destructor.
     */
    virtual ~OBCisTransStereo();

    ///@name Cis/Trans stereochemistry
    ///@{
    /**
     * Get the OBStereo::Type for this object.
     * @return OBStereo::CisTrans
     */
    OBStereo::Type GetType() const { return OBStereo::CisTrans; }
    /**
     * @return True if this object is valid. This object is valid if all these
     * conditions are met:
     * - @p begin != OBStereo::NoRef
     * - @p end != OBStereo::NoRef
     * - @p refs contains 4 elements
     */
    bool IsValid() const;

    /**
     * Set the configuration using a Config struct.
     */
#ifndef SWIG
    void SetConfig(const Config &config);
    /**
     * Get the configuration as Config struct.
     */
    Config GetConfig(OBStereo::Shape shape = OBStereo::ShapeU) const;
    /**
     * Get the configuration as Config struct and ensure refs[0] is
     * equal to @p start.
     */
    Config GetConfig(unsigned long start,
        OBStereo::Shape shape = OBStereo::ShapeU) const;
#endif
    /**
     * Compare the stereochemistry stored in the Config struct with the
     * stereochemistry specified in the Config struct from @p other.
     *
     * @copydoc Config::operator==()
     */
    bool operator==(const OBCisTransStereo &other) const;
    /**
     * Not equal to operator. This is the inverse of the Equal to operator==.
     *
     * @return True if the two Config structs represent a different stereochemistry.
     */
    bool operator!=(const OBCisTransStereo &other) const
    {
      return !(*this == other);
    }
    ///@}

    /*
     * Implement OBGenericData::Clone().
     */
    OBGenericData* Clone(OBBase *mol) const;


    //! @name Query methods to compare stereochemistry.
    //@{
    /**
     * Check if the two atoms for @p id1 & @p id2 are bonded to the same
     * atom. If the atoms for one of the atoms doesn't exist (anymore),
     * the valence of the begin and end atom is checked. If the exising
     * atom is bonded to the begin atom and end->GetExplicitDegree() == 2, the
     * ids are considered to be on different atoms. The reasoning behind
     * this is that hydrogens may be deleted. However, you can also use
     * OBStereo::ImplicitRef explicitly in code like:
     *
     * @code
     * //
     * //  F       F      F      F     0      3
     * //   \     /       |      |     |      |
     * //    C===C        | C  C |     | 1  2 |
     * //   /     \       |      |     |      |
     * // (H)     (H)    (H)----(H)    H------H
     * //
     * reading smiles F/C=C\F cis-difluorethene
     * OBCisTransStereo ct(mol);
     * ct.SetCenters(1, 2);
     * ct.SetRefs(OBStereo::MakeRefs(0, OBStereo::ImplicitRef, OBStereo::ImplicitRef, 3));
     * ...
     * @endcode
     *
     * @return True if @p id1 and @p id2 are bonded to the same atom
     * taking implicit hydrogens into account.
     */
    bool IsOnSameAtom(unsigned long id1, unsigned long id2) const;
    /**
     * @return True if the two reference ids are placed trans configuration.
     */
    bool IsTrans(unsigned long id1, unsigned long id2) const;
    /**
     * @return True if the two reference ids are placed in a cis configuration.
     */
    bool IsCis(unsigned long id1, unsigned long id2) const;
    /**
     * @image html gettransref.png
     * Get the reference id trans from reference @p id.
     */
    unsigned long GetTransRef(unsigned long id) const;
    /**
     * @image html getcisref.png
     * Get the reference id cis from reference @p id.
     */
    unsigned long GetCisRef(unsigned long id) const;
    //@}

  private:
    Config m_cfg; //!< internal configuration
    // The following function sits behind GetCisRef and GetTransRef
    unsigned long GetCisOrTransRef(unsigned long id, bool getcisref) const;
};
///@}

} // namespace OpenBabel

#ifndef SWIG
namespace std {

///@addtogroup stereo Stereochemistry
///@{
/**
 * @code
 * OBCisTransStereo::Config cfg;
 * cfg.begin = 0;
 * cfg.end = 1;
 * cfg.refs = OBStereo::MakeRefs(2, 3, 4, 5);
 * cfg.shape = OBStereo::ShapeU;
 *
 * OBCisTransStereo ct(mol);
 * ct.SetConfig(cfg)
 *
 * cout << "ct = " << ct << endl;
 *
 * // output
 * OBCisTransStereo(begin = 0, end = 1, refs = 2 3 4 5, shape = U)
 * @endcode
 */
OBAPI ostream& operator<<(ostream &out, const OpenBabel::OBCisTransStereo &ct);
/**
 * @code
 * OBCisTransStereo::Config cfg;
 * cfg.begin = 0;
 * cfg.end = 1;
 * cfg.refs = OBStereo::MakeRefs(2, 3, 4, 5);
 * cfg.shape = OBStereo::ShapeU;
 *
 * cout << "cfg = " << cfg << endl;
 *
 * // output
 * OBCisTransStereo::Config(begin = 0, end = 1, refs = 2 3 4 5, shape = U)
 * @endcode
 */
OBAPI ostream& operator<<(ostream &out, const OpenBabel::OBCisTransStereo::Config &cfg);

///@}

} // namespace std
#endif // Not SWIG

#endif

//! \file cistrans.h
//! \brief Store and convert cis/trans double-bond stereochemistry
