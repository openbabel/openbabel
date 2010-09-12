#include "stereoutil.h"




namespace OpenBabel {

  template<int StereoType>
  int findDescriptorTemplate(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symmetry_classes)
  {
    typedef OBStereoTypeTraits<StereoType> Traits;

    OBStereoFacade stereoFacade(mol, false);

    if (stereoFacade.HasStereo<StereoType>(unit.id)) {
      typename Traits::Config symConfig = stereoFacade.GetStereo<typename Traits::Stereo>(unit.id)->GetConfig();
      IdsToSymClasses(mol, symConfig, symmetry_classes);
      return findDescriptor(symConfig);
    }

    return 0;
  }

  int findDescriptor(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symmetry_classes)
  {
    switch (unit.type) {
      case OBStereo::Tetrahedral:
        return findDescriptorTemplate<OBStereo::Tetrahedral>(mol, unit, symmetry_classes);
      case OBStereo::CisTrans:
        return findDescriptorTemplate<OBStereo::CisTrans>(mol, unit, symmetry_classes);
    }
    return 0;
  }

  template<int StereoType>
  bool isSpecifiedTemplate(OBMol *mol, const OBStereoUnit &unit)
  {
    typedef OBStereoTypeTraits<StereoType> Traits;

    OBStereoFacade stereoFacade(mol, false);
    if (stereoFacade.HasStereo<StereoType>(unit.id))
      return stereoFacade.GetStereo<typename Traits::Stereo>(unit.id)->GetConfig().specified;

    return false;
  }
 
  bool isSpecified(OBMol *mol, const OBStereoUnit &unit)
  {
    switch (unit.type) {
      case OBStereo::Tetrahedral:
        return isSpecifiedTemplate<OBStereo::Tetrahedral>(mol, unit);
      case OBStereo::CisTrans:
        return isSpecifiedTemplate<OBStereo::CisTrans>(mol, unit);
    }
    return false;
  }

  OBStereoUnitSet removeUnspecifiedUnits(OBMol *mol, const OBStereoUnitSet &units)
  {
    OBStereoUnitSet result;
    for (std::size_t i = 0; i < units.size(); ++i)
      if (isSpecified(mol, units[i]))
        result.push_back(units[i]);
    return result;
  }

     
  int classifyNbrSymClasses(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symClasses)
  {
    switch (unit.type) {
      case OBStereo::Tetrahedral:
        {
          OBAtom *center = mol->GetAtomById(unit.id);
          if (!center)
            return 0;
          return classifyTetrahedralNbrSymClasses(symClasses, center);
        }
      case OBStereo::CisTrans:
        {
          OBBond *bond = mol->GetBondById(unit.id);
          if (!bond)
            return 0;
          if (classifyCisTransNbrSymClasses(symClasses, bond, bond->GetBeginAtom()) == C11)
            return C11;
          return classifyCisTransNbrSymClasses(symClasses, bond, bond->GetEndAtom());
        }
    }
    return 0;
  }

  bool isResolved(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symClasses)
  {
    switch (unit.type) {
      case OBStereo::Tetrahedral:
        return (classifyNbrSymClasses(mol, unit, symClasses) == T1234);
      case OBStereo::CisTrans:
        return (classifyNbrSymClasses(mol, unit, symClasses) == C12);
    }
    return false;  
  }

}
