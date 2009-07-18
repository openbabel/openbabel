#ifndef OB_CISTRANS_H
#define OB_CISTRANS_H

#include <openbabel/stereo/tetraplanar.h>
#include <vector>

namespace OpenBabel {

/**
 * @class OBCisTransStereo
 * @brief Class for handling and storing cis/trans stereochemistry.
 *
 * @image html cistrans.png 
 *
 * The OBCisTransStereo class is used to represent cis/trans stereochemistry.
 * Like all OBTetraPlanarStereo subclasses, it uses the OBStereo::Shape parameters
 * to set/get the reference ids. However, since the orientation of the double bond
 * matters, all methods in the "query methods" section check the bonding by actually
 * checking the atoms in the molecule provided through the constructor. If 
 * OBMol::GetAtomById(id) returns 0 for a single id, it will be considered a deleted 
 * hydrogen if the valences confirm this. This class uses the OBMessageHandler to 
 * report errors, warnings and info.
 *
 * An example:
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
 *
 *
 */
class OBAPI OBCisTransStereo : public OBTetraPlanarStereo
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
       * Equal to operator.
       *
       * @code
       * OBCisTransStereo::Config cfg1, cfg2;
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
      unsigned long begin, end; //<! The double bond begin and end ids.
      OBStereo::Refs refs; //!< The 4 reference ids.
      OBStereo::Shape shape; //!< The shape of the 4 reference ids.
      bool specified;
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

    /**
     * Get the OBStereo::Type for this object.
     * @return OBStereo::CisTrans
     */
    OBStereo::Type GetType() const { return OBStereo::CisTrans; }
    /**
     * @return True if this object is valid. This object is valid if all (center and
     * and reference) atom ids are set.
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
     * Compare the internally stored stereochemistry with the 
     * stereochemistry specified by @p other.
     *
     * @return True if both OBTetrahedralStereo objects represent the same
     * stereochemistry.
     */
    bool operator==(const OBCisTransStereo &other) const;
    bool operator!=(const OBCisTransStereo &other) const
    {
      return !(*this == other); 
    }
    
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
     * atom is bonded to the begin atom and end->GetValence() == 2, the 
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
};

} // namespace OpenBabel
#ifndef SWIG
namespace std {

/**
 * @code
 * OBTetrahedralStereo::Config cfg;
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

} // namespace std
#endif // Not SWIG

#endif
