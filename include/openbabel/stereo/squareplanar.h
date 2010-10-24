/**********************************************************************
  squareplanar.h - Class for handling and storing squareplanar stereochemistry.

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
#ifndef OB_SQUAREPLANAR_H
#define OB_SQUAREPLANAR_H

#include <openbabel/stereo/tetraplanar.h>
#include <vector>

namespace OpenBabel {

///@addtogroup stereo Stereochemistry
///@{
/**
 * @class OBSquarePlanarStereo squareplanar.h <openbabel/stereo/squareplanar.h>
 * @brief Class for handling and storing square planar stereochemistry.
 *
 * @image html squareplanar.png
 *
 * The OBSquarePlanarStereo class is used to represent square planar stereochemistry.
 * Like all OBTetraPlanarStereo subclasses, it uses the OBStereo::Shape parameters
 * to set/get the reference ids.
 *
 * This class works with reference ids only, it never uses the molecule
 * to get more information. Like all stereo classes, errors,
 * warnings or info is reported using OBMessageHandler.
 */
class OBAPI OBSquarePlanarStereo : public OBTetraPlanarStereo
{
  public:
    /**
     * \struct Config squareplanar.h <openbabel/stereo/squareplanar.h>
     * \brief Stereochemical configuration for square planar stereocenters
     *
     * The config struct represents the stereochemistry in a well defined way.
     * For squareplanar stereocenters, the following data members define the spacial
     * arrengement of the atoms.
     *
     * - OBStereo::Ref @p center: The central atom.
     * - OBStereo::Refs @p refs: The 4 atoms connected to the double bond.
     * - OBStereo::Shape @p shape: The shape formed by the @p refs by connecting them
     *   in the same order as they occur in @p refs.
     *
     * @image html squareplanar.png
     * @image html SPshapes.png
     *
     * Only @p center are specific for OBSquarePlanarStereo::Config. The other
     * data members occur in all OBTetraPlanarStereo derived classes.
     */
#ifndef SWIG
    struct OBAPI Config
    {
      /**
       * Default constructor. Initializes @p center to OBStereo::NoRef
       * and @p shape to OBStereo::ShapeU.
       */
      Config() : center(OBStereo::NoRef), shape(OBStereo::ShapeU),
          specified(true)
      {  }
      /**
       * Constructor with all parameters.
       *
       * @param _center The atom id for the central atom.
       * @param _refs The 4 reference ids.
       * @param _shape The shape for the 4 reference ids.
       */
      Config(unsigned long _center, const OBStereo::Refs &_refs,
          OBStereo::Shape _shape = OBStereo::ShapeU) : center(_center),
          refs(_refs), shape(_shape), specified(true)
      {  }
      /**
       * Equal to operator. Comparing OBSquarePlanarStereo::Config structs
       * is done using the information stored in the struct's data members
       * (i.e. center, refs and shape).
       *
       * There are a number of cases resuling in false being returned:
       * - @p center atom ids don't match
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
         1   4    1      4    1------4
          \ /     |      |           |
           C      |      |           |
          / \     |      |           |
         2   3    2------3    2------3

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
      unsigned long center; //<! The central atom id.
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
    OBSquarePlanarStereo(OBMol *mol);
    /**
     * Destructor.
     */
    virtual ~OBSquarePlanarStereo();

    ///@name SquarePlanar stereochemistry
    ///@{
    /**
     * Get the OBStereo::Type for this object.
     * @return OBStereo::SquarePlanar
     */
    OBStereo::Type GetType() const { return OBStereo::SquarePlanar; }
    /**
     * @return True if this object is valid. This object is valid if all (center and
     * and 4 reference) atom ids are set.
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
    bool operator==(const OBSquarePlanarStereo &other) const;
    /**
     * Not equal to operator. This is the inverse of the Equal to operator==.
     *
     * @return True if the two Config structs represent a different stereochemistry.
     */
    bool operator!=(const OBSquarePlanarStereo &other) const
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
     * Get the reference id cis from reference @p id.
     */
    std::vector<unsigned long> GetCisRefs(unsigned long id) const;
    //@}

  private:
    Config m_cfg; //!< internal configuration
    // The following function sits behind GetCisRef and GetTransRef
    unsigned long GetCisOrTransRef(unsigned long id, bool getcisref) const;
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
 * OBSquarePlanarStereo::Config cfg;
 * cfg.center = 0;
 * cfg.refs = OBStereo::MakeRefs(1, 2, 3, 4);
 * cfg.shape = OBStereo::ShapeU;
 *
 * OBSquarePlanarStereo ct(mol);
 * ct.SetConfig(cfg)
 *
 * cout << "ct = " << ct << endl;
 *
 * // output
 * OBSquarePlanarStereo(center = 0, refs = 1 2 3 4, shape = U)
 * @endcode
 */
OBAPI ostream& operator<<(ostream &out, const OpenBabel::OBSquarePlanarStereo &ct);
/**
 * @code
 * OBSquarePlanarStereo::Config cfg;
 * cfg.center = 0;
 * cfg.refs = OBStereo::MakeRefs(1, 2, 3, 4);
 * cfg.shape = OBStereo::ShapeU;
 *
 * cout << "cfg = " << cfg << endl;
 *
 * // output
 * OBSquarePlanarStereo::Config(center = 0, refs = 1 2 3 4, shape = U)
 * @endcode
 */
OBAPI ostream& operator<<(ostream &out, const OpenBabel::OBSquarePlanarStereo::Config &cfg);

///@}

} // namespace std
#endif // Not SWIG

#endif

//! \file squareplanar.h
//! \brief Store and convert square-planar stereochemistry
