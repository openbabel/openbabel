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
class OBCisTransStereo : public OBTetraPlanarStereo
{
  public:
    OBCisTransStereo(OBMol *mol);
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
     * Set the central, double bonded atoms
     */
    void SetCenters(unsigned long begin, unsigned long end);
    /**
     * Set the begin atom for the double bond
     */
    void SetBegin(unsigned long begin);
    /**
     * Set the end atom for the double bond
     */
    void SetEnd(unsigned long end);
    /**
     * @return The double bond begin atom.
     */
    unsigned long GetBegin() const;
    /**
     * @return The double bond begin atom.
     */
    unsigned long GetEnd() const;

    //! @name Methods to get and set the reference ids.
    //@{
    /**
     * Set the 4 reference ids. The @p shape parameter specifies how the 
     * reference ids are layed out.
     *
     * @param refs The 4 reference ids.
     * @param shape The reference id order in the returned list.
     * These are all examples of the same configuration:
     * @image html SPshapes.png
     */
    void SetRefs(const std::vector<unsigned long> &refs, 
        OBStereo::Shape shape = OBStereo::ShapeU);
    /**
     * Get a list of the 4 reference ids. The ids occur in the sequence 
     * following the specified shape.
     *
     * @param shape The reference id order in the returned list.
     * These are all examples of the same configuration:
     * @image html SPshapes.png
     */
    std::vector<unsigned long> GetRefs(OBStereo::Shape shape = OBStereo::ShapeU) const;
    //@}

    //! @name Query methods to compare stereochemistry.
    //@{
    /**
     * Check if the two atoms for @p id1 & @p id2 are bonded to the same 
     * atom. If the atoms for one of the atoms doesn't exist (anymore), 
     * the valence of the begin and end atom is checked. If the exising 
     * atom is bonded to the begin atom and end->GetValence() == 2, the 
     * ids are considered to be on different atoms. The reasoning behind
     * this is that hydrogens may be deleted. However, you can also use
     * OBStereo::HydrogenId explicitly in code like:
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
     * ct.SetRefs(OBStereo::MakeRefs(0, OBStereo::HydrogenId, OBStereo::HydrogenId, 3));
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
    /**
     * This function checks to see if the internal reference ids match with @p refs.
     * This is done by checking which reference ids are trans. See the OBTetraPlanarStereo
     * class documentation for more information about the 3 possible configurations.
     * \return True if the stored reference ids match the configuration 
     * of @p refs using @p shape. 
     */
    bool Compare(const std::vector<unsigned long> &refs, OBStereo::Shape shape) const;
    //@}
 

  private:

    unsigned long m_begin; //!< the center/chiral atom
    unsigned long m_end; //!< the center/chiral atom
    std::vector<unsigned long> m_refs; //!< the 4 clockwise atom refs
};

}

#endif
