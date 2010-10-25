/**********************************************************************
  tetrahedral.h - OBTetrahedralStereo

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
#ifndef OB_TETRAHEDRAL_H
#define OB_TETRAHEDRAL_H

#include <openbabel/stereo/tetranonplanar.h>

namespace OpenBabel {

///@addtogroup stereo Stereochemistry
///@{
/**
 * @class OBTetrahedralStereo tetrahedral.h <openbabel/stereo/tetrahedral.h>
 * @brief Class for handling and storing tetrahedral atom stereochemistry.
 *
 * The OBTetrahedralStereo class is used to represent tetrahedral atom
 * stereo chemistry. Before continuing reading, it is recommeded to
 * read the OBTetraNonPlanar documentation first. Since this
 * class inherits OBTetraNonPlanerStereo it should be clear that the class
 * stores the spacial arrangement of 4 non-planer atoms. However, for the
 * tetrahedral case, one additional OBStereo::Ref is needed: the atom id for
 * the center atom.
 *
 * The stereochemistry is set, retrieved and internally stored using the
 * OBTetrahedralStereo::Config struct. OBTetrahedralStereo functions as a
 * wrapper around the Config struct and provides functions to get a copy of
 * the Config struct converted to any view direction or atom and winding
 * (OBTetrahedralStereo::GetConfig).
 *
 * @image html tetrahedral.png
 *
 * Like all stereo classes, errors, warnings or info is reported using OBMessageHandler.
 *
 * @sa OBStereo OBStereoBase OBTetraNonPlanarStereo OBStereoFacade
 * @since version 2.3
 */
class OBAPI OBTetrahedralStereo : public OBTetraNonPlanarStereo
{
  public:
#ifndef SWIG
    /**
     * \struct Config tetrahedral.h <openbabel/stereo/tetrahedral.h>
     * \brief Stereochemical configuration for tetrahedral stereocenters
     *
     * The config struct represents the stereochemistry in a well defined way. For
     * tetrahedral stereo centers, the following data members define the special
     * orientation of the atoms:
     *
     * - OBStereo::Ref @p center: Atom id of the stereogenic center atom.
     * - OBStereo::Ref @p from/towards: Atom id (or OBStereo::ImplicitRef) for the
     *   atom to view from/towards.
     * - OBStereo::Refs @p refs: The three remaining atom ids (may also contain one
     *   OBStereo::NoRef element if from/towards is set to a real atom id).
     * - OBStereo::View @p view: Specify the viewing from or towards the atom with
     *   @p from/towards id.
     * - OBStereo::Winding @p winding: Clockwise or AntiClockwise (order in the Refs @p refs list)
     *
     * @image html tetrahedral.png
     *
     * Only @p center is specific for OBTetrahedralStereo::Config. The other
     * data members occur in all OBTetraNonPlanarStereo derived classes.
     */
    struct OBAPI Config
    {
      /**
       * Default constructor. Initializes the @p from/torards and @p center to
       * OBStereo::NoRef, the @p winding to OBStereo::Clockwise and @p view to
       * OBStereo::ViewFrom.
       */
      Config() : center(OBStereo::NoRef), from(OBStereo::NoRef),
          winding(OBStereo::Clockwise), view(OBStereo::ViewFrom),
          specified(true)
      {  }
      /**
       * Constructor with all parameters.
       *
       * @param _center The center (chiral) atom id.
       * @param from_or_towards The atom id from which to view or view towards (see @p view).
       * @param _refs The 3 reference ids.
       * @param _winding The winding for the 3 ids in @p _refs.
       * @param _view Specify viewing from or towards the atom with @p from_or_towards id.
       */
      Config(unsigned long _center, unsigned long from_or_towards,
          const OBStereo::Refs &_refs, OBStereo::Winding _winding = OBStereo::Clockwise,
          OBStereo::View _view = OBStereo::ViewFrom) : center(_center),
          from(from_or_towards), refs(_refs), winding(_winding), view(_view),
          specified(true)
      {  }
      /**
       * Equal to operator. Comparing OBTetrahedralStereo::Config structs
       * is done using the information stored in the struct's data members
       * (i.e. @p view, @p winding, @p from/towards and @p refs).
       *
       * There are a number of cases resuling in false being returned:
       * - The centers don't match.
       * - One of the Refs lists does not contain 3 elements.
       * - 2 or more OBStereo::ImplicitRef values in a single Config struct
       *
       * When either Config struct is unspecified (i.e. the stereochemistry
       * implied is accidental), true is returned.
       *
       * It doesn't matter if the two Config structs use the same view, same
       * from/towards Ref or the same winding. All needed conversions will be
       * carried out automatically (see OBTetraNonPlanerStereo::ToConfig). These
       * conversions ensure the spacial orientation of the 4 groups remains
       * unchanged.
       *
       * Another key feature is the ability to comapre Config structs regardless
       * of implicit (OBStereo::ImplicitRef) or explicit hydrogens. This is best
       * illustrated with some examples. In these examples the same ref has
       * already been selected as from/towards atom and both use the same winding
       * and view direction. We will focus on how the three remaining refs are
       * interpreted.
       *
       @verbatim
         234 == 234 // true
         2H4 == 234 // 3 is missing, must be the implicit --> 234 == 234 // true
         2H4 == 243 // same as above, but now 234 == 243 // false
         234 == H34 // 2 is missing, must be implicit --> 234 == 234 // true
       @endverbatim
       *
       * By comparing the second and third example above, it can be clearly seen
       * that the value of 1 Ref can actually be ignored. It's position in the
       * sequence (or the winding) is defined by the two explicit Ref values.
       *
       * @return True if both Config structs represent the same stereochemistry.
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
      unsigned long center; //<! The center (chiral) atom id.
      /**
       * This anonymous union helps to keep code clean. Both the @p from and
       * @p towards data members contain the same OBStereo::Ref (same memory
       * address) but can be used interchangeably to match the context of
       * the code. The real viewing direction is specified by the @p view
       * data member.
       */
      union {
        unsigned long from; //<! The viewing from atom id.
        unsigned long towards; //<! The viewing towards id.
      };
      OBStereo::Refs refs; //!< The 3 reference ids.
      OBStereo::Winding winding; //<! The winding for the 3 reference ids.
      OBStereo::View view; //!< Specify viewing from or towards the atom with @p from/towards id.
      bool specified; //!< True if the stereochemistry is specified. When false, the described
                      //!< special orientation is only accidental (i.e. unspecified).
      //@}
    };
#endif
    /**
     * Constructor.
     */
    OBTetrahedralStereo(OBMol *mol);
    /**
     * Destructor.
     */
    virtual ~OBTetrahedralStereo();

    ///@name Tetrahedral stereochemistry
    ///@{
    /**
     * Get the OBStereo::Type for this object.
     * @return OBStereo::Tetrahedral
     */
    OBStereo::Type GetType() const { return OBStereo::Tetrahedral; }
    /**
     * @return True if this object is valid. This object is valid if all (center, from
     * and ref) atom ids are set.
     */
    bool IsValid() const;
#ifndef SWIG
    /**
     * Set the configuration using a Config struct.
     */
    void SetConfig(const Config &config);
    /**
     * Get the configuration as Config struct.
     */
    Config GetConfig(OBStereo::Winding winding = OBStereo::Clockwise,
        OBStereo::View view = OBStereo::ViewFrom) const;
    /**
     * Get the configuration as Config struct viewing from/towards the specified id.
     */
    Config GetConfig(unsigned long from_or_towards,
        OBStereo::Winding winding = OBStereo::Clockwise,
        OBStereo::View view = OBStereo::ViewFrom) const;
#endif
    /**
     * Compare the stereochemistry stored in the Config struct with the
     * stereochemistry specified in the Config struct from @p other.
     *
     * @copydoc Config::operator==()
     */
    bool operator==(const OBTetrahedralStereo &other) const;
    /**
     * Not equal to operator. This is the inverse of the Equal to operator==.
     *
     * @return True if the two Config structs represent a different stereochemistry.
     */
    bool operator!=(const OBTetrahedralStereo &other) const
    {
      return !(*this == other);
    }
    ///@}

    /*
     * Implement OBGenericData::Clone().
     */
    OBGenericData* Clone(OBBase *mol) const;
  private:
    Config m_cfg; //!< internal configuration
};
///@}
// end addtogroup doxygen

} // namespace OpenBabel

#ifndef SWIG
namespace std {

///@addtogroup stereo Stereochemistry
///@{
/**
 * @code
 * OBTetrahedralStereo::Config cfg;
 * cfg.center = 0;
 * cfg.towards = 4;
 * cfg.refs = OBStereo::MakeRefs(1, 2, 3);
 * cfg.winding = OBStereo::AntiClockwise;
 * cfg.view = OBStereo::ViewTowards;
 *
 * OBTetrahedralStereo ts(mol);
 * ts.SetConfig(cfg)
 *
 * cout << "ts = " << ts << endl;
 *
 * // output
 * OBTetrahedralStereo(center = 0, viewTowards = 4, refs = 1 2 3, anti-clockwise)
 * @endcode
 */
OBAPI ostream& operator<<(ostream &out, const OpenBabel::OBTetrahedralStereo &ts);
/**
 * @code
 * OBTetrahedralStereo::Config cfg;
 * cfg.center = 0;
 * cfg.towards = 4;
 * cfg.refs = OBStereo::MakeRefs(1, 2, 3);
 * cfg.winding = OBStereo::AntiClockwise;
 * cfg.view = OBStereo::ViewTowards;
 *
 * cout << "cfg = " << cfg << endl;
 *
 * // output
 * OBTetrahedralStereo::Config(center = 0, viewTowards = 4, refs = 1 2 3, anti-clockwise)
 * @endcode
 */
OBAPI ostream& operator<<(ostream &out, const OpenBabel::OBTetrahedralStereo::Config &cfg);

///@}

} // namespace std
#endif // SWIG

#endif

//! \file tetrahedral.h
//! \brief Handle general non-planar tetrahedral stereochemistry
