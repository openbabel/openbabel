#ifndef OB_STEREO_H
#define OB_STEREO_H

#include <openbabel/base.h> // OBGenericData
#include <vector>
#include <algorithm>

namespace OpenBabel {

  struct OBStereo 
  {
    /**
     * The various types of stereochemistry
     *
     */
    enum Type {
      None                = (1<<0), //!< no stereochemistry
      Unknown             = (1<<2), //!< unknown
      Unspecified         = (1<<3), //!< not specified
      CisTrans            = (1<<4), //!< cis/trans double bond
      ExtendedCisTrans    = (1<<5), //!< allene, biphenyl, ...
      SquarePlanar        = (1<<6), //!< Square-planar stereochemistry
      TetraPlanar         = CisTrans | ExtendedCisTrans | SquarePlanar, //!< planar configurations od four atoms
      Tetrahedral         = (1<<7), //!< tetrahedral
      ExtendedTetrahedral = (1<<8), //!< extended tetrahedral
      TetreNonPlanar      = Tetrahedral | ExtendedTetrahedral
      /*
      TrigonalBipyramidal = 4, //!< Trigonal-bipyramidal stereochemistry
      Octahedral          = 5, //!< Octahedral stereochemistry
      DoubleBond          = 6, //!< Double bond stereochemistry
      */
    };

    /**
     * Shapes used by OBTetraPlanarStereo subclasses for 
     * setting/getting reference ids.
     */
    enum Shape {
      ShapeU = 1,
      ShapeZ = 2,
      Shape4 = 3
    };

    /**
     * Views used by OBTetraNonPlanarStereo subclasses for
     * setting/getting reference ids.
     */
    enum View
    {
      ViewFrom = 1, //!< view from the atom (id parameter) towards the center atom
      ViewTowards = 2, //!< view from center atom towards the atom (id paramter)
    };

    /**
     * Windings used by OBTetraNonPlanar subclasses for 
     * setting/getting reference ids.
     */
    enum Winding {
      Clockwise = 1,     //!< Clockwise winding
      AntiClockwise = 2  //!< AntiClockiwe winding (or CounterClockwise
    };

    /**
     * Some useful predefined ids. 
     */
    enum {
      NoId = -1,       //!< no id
      HydrogenId = -2  //!< hydrogen, unknown id
    };

    typedef std::vector<unsigned long> Refs;
    typedef std::vector<unsigned long>::iterator RefIter;

    /**
     * Create a std::vector<unsigned long> filled with @p id1, @p id2, @p id3 & @p id4.
     */
    static std::vector<unsigned long> MakeRefs(unsigned long id1, unsigned long id2,
        unsigned long id3, unsigned long id4 = NoId)
    {
      std::vector<unsigned long> refs(3);
      refs[0] = id1;
      refs[1] = id2;
      refs[2] = id3;
      if (id4 != NoId)
        refs.push_back(id4);
      return refs;
    }

  };

  class OBMol;
  class OBStereoBase : public OBGenericData
  {
    public:
      OBStereoBase(OBMol *mol) : 
        OBGenericData("StereoData", OBGenericDataType::StereoData, perceived),
        m_mol(mol) 
      {
      }
      virtual ~OBStereoBase() { m_mol = 0; }
      /**
       * Get the molecule. This can be used by subclasses when more
       * information is needed (e.g. OBCisTransStereo::GetCisRef, ...).
       */
      OBMol* GetMolecule() const { return m_mol; }
      /**
       * Reimplemented by subclasses to return type.
       */
      virtual OBStereo::Type GetType() const = 0;
    private:
      OBMol *m_mol; //!< the parent molecule
  };


}

#endif
