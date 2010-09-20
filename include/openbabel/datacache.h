#include <openbabel/bitvec.h>
#include <openbabel/isomorphism.h>
#include <openbabel/stereo/stereo.h>

namespace OpenBabel {

  class OBMol;

  class OBAPI DataCache
  {
    public:
      DataCache(OBMol *mol);

      void GetSymmetryClasses(std::vector<unsigned int> &symmetry_classes, const OBBitVec &mask = OBBitVec(), bool optimizeMask = true);

      void GetAutomorphisms(Automorphisms &automorphisms, const OBBitVec &mask = OBBitVec());

      void GetStereogenicUnits(OBStereoUnitSet &units, const OBBitVec &mask = OBBitVec(), bool optimizeMask = true);

    private:
      void RemoveHydrogensFromMask(OBBitVec &output, const OBBitVec &input);
      void OptimizeStereoMask(OBBitVec &output, const OBBitVec &input);
      void Validate();

      bool m_validated;
      OBMol *m_mol;
  };

}
