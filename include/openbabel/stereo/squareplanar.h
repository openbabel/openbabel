#ifndef OB_SQUAREPLANAR_H
#define OB_SQUAREPLANAR_H

#include <openbabel/stereo/tetraplanar.h>
#include <vector>

namespace OpenBabel {

///@addtogroup stereo Stereochemistry
///@{
/**
 * @class OBSquarePlanarStereo
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
///@}
{
  public:
    OBSquarePlanarStereo(OBMol *mol);
    virtual ~OBSquarePlanarStereo();

    /**
     * Get the OBStereo::Type for this object.
     * @return OBStereo::SquarePlanar
     */
    OBStereo::Type GetType() const { return OBStereo::SquarePlanar; }
    /**
     * @return True if this object is valid. This object is valid if all (center and
     * and reference) atom ids are set.
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
    /**
     * Get a list of the 4 reference ids. The first id will be the @p start reference id.
     *
     * @param start The first id in the sequence. The other ids will be adjusted.
     * @param shape The reference id order in the returned list.
     * These are all examples of the same configuration:
     * @image html SPshapes.png
     */
    std::vector<unsigned long> GetRefs(unsigned long start, 
        OBStereo::Shape shape = OBStereo::ShapeU) const;
    //@}

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
     * Get the reference id trans from reference @p id.
     */
    unsigned long GetTransRef(unsigned long id) const;
    /**
     * Get the reference id cis from reference @p id.
     */
    std::vector<unsigned long> GetCisRefs(unsigned long id) const;
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
    unsigned long m_center; //!< the center/chiral atom
    std::vector<unsigned long> m_refs; //!< the 4 clockwise atom refs
};

}

#endif
