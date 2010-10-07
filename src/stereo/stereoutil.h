#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>

namespace OpenBabel {

  /**
   * These are the different types of stereocenter classification used throughout
   * this file.
   */
  enum NeighborSymmetryClasses
  {
    // Tetrahedral
    T1234 = 1234, // 4 different symmetry classes
    T1123 = 1123, // 3 different symmetry classes, 1 class duplicated (2 times)
    T1122 = 1122, // 2 different symmetry classes, 1 class duplicated (3 times)
    T1112 = 1112, // 2 different symmetry classes, each class duplicated (2 times)
    T1111 = 1111, // 1 symmetry class, duplictaed 4 times
    // CisTrans
    C12 = 12, // 2 different symmetry classes
    C11 = 11 // the same symmetry class
  };


// defined in src/stereo/perception.cpp
  int classifyTetrahedralNbrSymClasses(const std::vector<unsigned int> &symClasses, OBAtom *atom);
  int classifyCisTransNbrSymClasses(const std::vector<unsigned int> &symClasses, OBBond *doubleBond, OBAtom *atom);
  OBStereoUnitSet orderSetBySymmetryClasses(OBMol *mol, const OBStereoUnitSet &set,
      const std::vector<unsigned int> &symmetry_classes);
  unsigned int findDuplicatedSymmetryClass(OBAtom *atom, const std::vector<unsigned int> &symClasses);
  void findDuplicatedSymmetryClasses(OBAtom *atom, const std::vector<unsigned int> &symClasses,
      unsigned int &duplicated1, unsigned int &duplicated2);
  OBBitVec getFragment(OBAtom *atom, OBAtom *skip);
  bool isTetrahedral(OBAtom *atom, const OBStereoUnitSet &units);
  bool isCisTrans(OBBond *bond, const OBStereoUnitSet &units);
  unsigned int getSymClass(const OBStereoUnit &unit, OBMol *mol, const std::vector<unsigned int> &symmetry_classes);

// graphsym.cpp
  void IdsToSymClasses(OBMol *mol, OBCisTransStereo::Config &config,
      const std::vector<unsigned int> &symClasses);
  void IdsToSymClasses(OBMol *mol, OBTetrahedralStereo::Config &config,
      const std::vector<unsigned int> &symClasses);
  bool setsContainSameSymmetryClasses(OBMol *mol, const OBStereoUnitSet &set1,
      const OBStereoUnitSet &set2, const std::vector<unsigned int> &symmetry_classes);
  void orderSetByFragmentRecursive(OBStereoUnitSet &ordered, OBAtom *atom, OBAtom *skip,
      const OBStereoUnitSet &set, OBBitVec &fragment);
  OBStereoUnitSet orderSetByFragment(OBAtom *atom, OBAtom *skip, const OBStereoUnitSet &set);
  OBStereoUnitSet orderSetByRing(OBMol *mol, const OBBitVec &fragment, const OBStereoUnitSet &set);
  int findDescriptor(const OBTetrahedralStereo::Config &config);
  int findDescriptor(const OBCisTransStereo::Config &config);
  int findDescriptor(OBAtom *center, const std::vector<unsigned int> &symmetry_classes);
  bool isUnitResolved(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symmetry_classes);
  int findDescriptor(OBBond *bond, const std::vector<unsigned int> &symmetry_classes);
  std::vector<int> findDescriptorVector(OBMol *mol, const OBBitVec &fragment,
      const OBStereoUnitSet &orderedUnits, const std::vector<unsigned int> &symmetry_classes);
  int findDescriptorVectorValue(OBMol *mol, const OBBitVec &fragment,
      const OBStereoUnitSet &orderedUnits, const std::vector<unsigned int> &symmetry_classes);
  OBStereoUnitSet findSetContainingUnit(const OBStereoUnitSetOfSets &sets, const OBStereoUnit &unit);

  bool isUnitInFragment(OBMol *mol, const OBStereoUnit &unit, const OBBitVec &fragment);


  template<int StereoType> struct OBStereoTypeTraits;
  template<> struct OBStereoTypeTraits<OBStereo::Tetrahedral>
  {
    enum {
      Type = OBStereo::Tetrahedral // the stereo type
    };
    typedef OBAtom Center; // the stereogenic center type
    typedef OBTetrahedralStereo Stereo; // The stereo data object
    typedef OBTetrahedralStereo::Config Config; // the config struct
    static Center* GetCenter(OBMol *mol, unsigned long id) { return mol->GetAtomById(id); }
  };
  template<> struct OBStereoTypeTraits<OBStereo::CisTrans>
  {
    enum {
      Type = OBStereo::CisTrans,
    };
    typedef OBBond Center;
    typedef OBCisTransStereo Stereo;
    typedef OBCisTransStereo::Config Config;
    static Center* GetCenter(OBMol *mol, unsigned long id) { return mol->GetBondById(id); }
  };

  template<typename OBAtomOrOBBond> struct OBStereoCenterTraits;
  template<> struct OBStereoCenterTraits<OBAtom> : OBStereoTypeTraits<OBStereo::Tetrahedral> {};
  template<> struct OBStereoCenterTraits<OBBond> : OBStereoTypeTraits<OBStereo::CisTrans> {};

  template<typename OBAtomOrOBBond>
  bool isInSet(OBAtomOrOBBond *atomOrBond, const OBStereoUnitSet &units)
  {
    for (std::size_t i = 0; i < units.size(); ++i) {
      const OBStereoUnit &unit = units[i];
      if (unit.type != static_cast<OBStereo::Type>(OBStereoCenterTraits<OBAtomOrOBBond>::Type))
        continue;
      if (unit.id == atomOrBond->GetId())
        return true;
    }
    return false;
  }

  int findDescriptor(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symmetry_classes, OBStereoFacade &stereoFacade);
  int classifyNbrSymClasses(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symmetry_classes);
  bool isResolved(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symClasses);
  bool isSpecified(OBMol *mol, const OBStereoUnit &unit, OBStereoFacade &stereoFacade);
  OBStereoUnitSet removeUnspecifiedUnits(OBMol *mol, const OBStereoUnitSet &units, OBStereoFacade &stereoFacade);

}
