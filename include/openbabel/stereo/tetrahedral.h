#ifndef OB_TETRAHEDRAL_H
#define OB_TETRAHEDRAL_H

#include <openbabel/stereo/tetranonplanar.h>
#include <vector>

namespace OpenBabel {

/**
 * @class OBTetrahedralStereo
 * @brief Handle and store tetrahedral atom stereo chemistry.
 *
 * @image html tetrahedral.png
 *
 * The OBTetrahedralStereo class is used to represent tetrahedral atom 
 * stereo chemisrty. It uses an internal representation similar to the
 * SMILES model. When viewing from the 'from' atom towards the center atom,
 * the three reference atoms are stored clockwise. However, it is important to 
 * note that this class doesn't take any ranking of the atoms (by index, 
 * CIP rules, ...) into account. All references to atoms are stored using 
 * unique ids which are guaranteed to remain unchanged during an OBMol
 * object's lifespan (even when indexes change after canonicalization for 
 * example). In other words, the stereo chemistry around the central atom 
 * is defined in a local way. This is different from CIP rankings or canonical 
 * label based rankings which may change when parts of the molecule get 
 * deleted or new parts are added (both are based on the molecule as a 
 * complete graph). Another benefit is that this class can correctly be
 * used in substructure search.
 * 
 * Methods are provided to set the atom ids in various ways. When using
 * these, the input will be converted to the internal representation.
 *
 *
 *
 *
 */
class OBAPI OBTetrahedralStereo : public OBTetraNonPlanarStereo
{
  public:
    OBTetrahedralStereo(OBMol *mol);
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
    /**
     * Set the central atom. This is the chiral atom.
     */
    void SetCenter(unsigned long id);
    /**
     * @return The central atom. This is the chiral atom.
     */
    unsigned long GetCenter() const;

    //! @name Methods using internal representation.
    //@{

    /**
     * Set the 4 reference atoms. It is possible to give the ids in clockwise or 
     * anti-clockwise winding as long as you specify the right one using the @p winding 
     * parameter. The ids will automatically be converted to the internal 
     * representation (clockwise). The @p view parameter specifies if the @p refs
     * are @p winding from or towards the atom with id @p id.
     *
     * @param refs std::vector containing the 3 reference ids. This vector should always be size 3.
     * @param id The ViewFrom or ViewTowards atom.
     * @param winding Specify the winding for the 3 reference ids.
     * @param view OBStereo::ViewFrom to OBStereo::ViewTowards for viewing from or 
     * towards the atom specified by @p id.
     */
    void SetRefs(const std::vector<unsigned long> &refs, unsigned long id, 
        OBStereo::Winding winding = OBStereo::Clockwise,
        OBStereo::View view = OBStereo::ViewFrom);
    /**
     * Get the 4 reference atoms. It is possible to retrieve the ids in clockwise or 
     * anticlockwise order using the @p winding parameter..
     *
     * @note We look from the from atom (GetFromAtom()) to the center atom (GetCenterAtom())
     *
     * @note We look from the ViewFrom atom towards the center atom.
     * @note We look towards the ViewTowards# atom from the center atom.
     *
     * @param id The ViewFrom or ViewTowards atom.
     * @param winding Specify the winding for the returned reference ids.
     * @param view OBStereo::ViewFrom to OBStereo::ViewTowards for viewing from or 
     * towards the atom specified by @p id.
     *
     * @return std::vector containing the 3 reference ids in the requested form. This vector 
     * is empty if this object is not valid (IsValid()).
     */
    std::vector<unsigned long> GetRefs(unsigned long id,
        OBStereo::Winding winding = OBStereo::Clockwise,
        OBStereo::View view = OBStereo::ViewFrom) const;
    //@}

    //! @name Query methods to compare stereochemistry.
    //@{
    /**
     * Compare the internally stored stereochemistry with the 
     * stereochemistry specified by the parameters.
     */
    bool Compare(const std::vector<unsigned long> &refs, unsigned long id, 
        OBStereo::Winding winding, OBStereo::View view) const;
    //@}

    OBGenericData* Clone(OBBase *mol) const;
  private:
    unsigned long m_center; //!< the center/chiral atom
    std::vector<unsigned long> m_refs; //!< the 4 clockwise atom refs
};

}

#endif
