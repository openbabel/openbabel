/**********************************************************************
  tetrahedral.h - OBTetrahedralStereo

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
#ifndef OB_TETRAHEDRAL_H
#define OB_TETRAHEDRAL_H

#include <openbabel/stereo/tetranonplanar.h>

namespace OpenBabel {

/**
 * @class OBTetrahedralStereo
 * @brief Handle and store tetrahedral atom stereo chemistry.
 *
 * @image html tetrahedral.png
 *
 * The OBTetrahedralStereo class is used to represent tetrahedral atom 
 * stereo chemistry. The stereochemistry is set, retrieved and internally
 * stored using the OBtetrahedralStereo::Config struct.
 *  
 * Like all stereo classes, errors, warnings or info is reported using OBMessageHandler.
 */
class OBAPI OBTetrahedralStereo : public OBTetraNonPlanarStereo
{
  public:
    /**
     * The config struct represents the stereochemistry in a well defined way.
     */
#ifndef SWIG
    struct OBAPI Config
    {
      /**
       * Default constructor.
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
       * (i.e. view, winding, from/towards and refs).
       *
       * If the centers don't match false is returned. If one of the Refs lists
       * does not contain 3 elements, false is returned. 2 or more 
       * OBStereo::ImplicitRef values in a single Config struct also result in
       * false being returned.
       *
       * It doesn't matter if the two Config structs use the same view, same 
       * from/towards Ref or the same winding. All needed conversions will be
       * carried out automatically.
       *
       * Another key feature is the ability to comapre Config structs regardless
       * of implicit (OBStereo::ImplicitRef) or explicit hydrogens. This is best 
       * illustrated with some examples. In these examples the same ref has 
       * already been selected as from/towards atom before. We will focus on how 
       * the three remaining refs are interpreted.
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
       * @code
       * OBTetrahedralStereo::Config cfg1, cfg2;
       * ...
       * if (cfg1 == cfg2) {
       *   // cfg1 and cfg2 represent the same stereochemistry
       *   ...
       * }
       * @endcode
       *
       * @return True if both Config structs represent the stereochemistry.
       */
      bool operator==(const Config &other) const;
      /**
       * Not equal to operator.
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
       * @p towards data members contain the same id (same memory address) 
       * but can be used interchangeably to match the context of the code. 
       * The real viewing direction is specified by the @p view data member.
       */
      union {
        unsigned long from; //<! The viewing from atom id.
        unsigned long towards; //<! The viewing towards id.
      };
      OBStereo::Refs refs; //!< The 3 reference ids.
      OBStereo::Winding winding; //<! The winding for the 3 reference ids.
      OBStereo::View view; //!< Specify viewing from or towards the atom with @p from/towards id.
      bool specified;
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
     * Compare the internally stored stereochemistry with the 
     * stereochemistry specified by @p other.
     *
     * @return True if both OBTetrahedralStereo objects represent the same
     * stereochemistry.
     */
    bool operator==(const OBTetrahedralStereo &other) const;
    bool operator!=(const OBTetrahedralStereo &other) const
    {
      return !(*this == other); 
    }
    
    /*
     * Implement OBGenericData::Clone().
     */
    OBGenericData* Clone(OBBase *mol) const;
  private:
    Config m_cfg; //!< internal configuration 
};

} // namespace OpenBabel
#ifndef SWIG
namespace std {

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

} // namespace std
#endif // SWIG

#endif
