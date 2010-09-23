#include "stereoutil.h"
#include <openbabel/graphsym.h>




namespace OpenBabel {

  /*
  template<int StereoType>
  int findDescriptorTemplate(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symmetry_classes, OBStereoFacade &stereoFacade)
  {
    typedef OBStereoTypeTraits<StereoType> Traits;

    if (stereoFacade.HasStereo<StereoType>(unit.id)) {
      typename Traits::Config symConfig = stereoFacade.GetStereo<typename Traits::Stereo>(unit.id)->GetConfig();
      IdsToSymClasses(mol, symConfig, symmetry_classes);
      return findDescriptor(symConfig);
    }

    return 0;
  }

  int findDescriptor(OBMol *mol, const OBStereoUnit &unit, const std::vector<unsigned int> &symmetry_classes, OBStereoFacade &stereoFacade)
  {
    switch (unit.type) {
      case OBStereo::Tetrahedral:
        return findDescriptorTemplate<OBStereo::Tetrahedral>(mol, unit, symmetry_classes, stereoFacade);
      case OBStereo::CisTrans:
        return findDescriptorTemplate<OBStereo::CisTrans>(mol, unit, symmetry_classes, stereoFacade);
    }
    return 0;
  }
  */

  /**
   * Translate atom ids to symmetry classes for a CisTrans config struct.
   */
  void IdsToSymClasses(OBMol *mol, OBCisTransStereo::Config &config,
      const std::vector<unsigned int> &symClasses)
  {
    OBAtom *atom;
    // begin
    atom = mol->GetAtomById(config.begin);
    if (atom) {
      if (atom->IsHydrogen())
        config.begin = OBGraphSym::NoSymmetryClass;
      else
        config.begin = symClasses.at(atom->GetIndex());
    }
    // end
    atom = mol->GetAtomById(config.end);
    if (atom) {
      if (atom->IsHydrogen())
        config.end = OBGraphSym::NoSymmetryClass;
      else
        config.end = symClasses.at(atom->GetIndex());
    }
    // refs
    for (std::size_t i = 0; i < config.refs.size(); ++i) {
      atom = mol->GetAtomById(config.refs.at(i));
      if (atom) {
        if (atom->IsHydrogen())
          config.refs[i] = OBGraphSym::NoSymmetryClass;
        else
          config.refs[i] = symClasses.at(atom->GetIndex());
      }
    }
  }

  /**
   * Translate atom ids to symmetry classes for a Tetrahedral config struct.
   */
  void IdsToSymClasses(OBMol *mol, OBTetrahedralStereo::Config &config,
      const std::vector<unsigned int> &symClasses)
  {
    OBAtom *atom;
    // center
    atom = mol->GetAtomById(config.center);
    if (atom) {
      if (atom->IsHydrogen())
        config.center = OBGraphSym::NoSymmetryClass;
      else
        config.center = symClasses.at(atom->GetIndex());
    }
    // from/towards
    atom = mol->GetAtomById(config.from);
    if (atom) {
      if (atom->IsHydrogen())
        config.from = OBGraphSym::NoSymmetryClass;
      else
        config.from = symClasses.at(atom->GetIndex());
    }
    // refs
    for (std::size_t i = 0; i < config.refs.size(); ++i) {
      atom = mol->GetAtomById(config.refs.at(i));
      if (atom) {
        if (atom->IsHydrogen())
          config.refs[i] = OBGraphSym::NoSymmetryClass;
        else
          config.refs[i] = symClasses.at(atom->GetIndex());
      }
    }
  }


  /**
   * Find a descriptor for the specified tetrahedral config struct. The
   * returned descriptor had a value of 0 or 1.
   *
   * @note These descriptors depend on the values in the config struct but are
   * not related to CIP descriptors for example. However, if the values in the
   * config structs are replaced with CIP priorities, this function could be
   * be used.
   */
  /*
  int findDescriptor(const OBTetrahedralStereo::Config &config)
  {
    std::vector<unsigned long> refs = config.refs;
    refs.insert(refs.begin(), config.from);
    if (OBStereo::NumInversions(refs) % 2)
      return 1;
    else
      return 0;
  }
  */

  /**
   * Find a descriptor for the specified cis/trans config struct. The
   * returned descriptor had a value of 0 or 1.
   *
   * @note These descriptors depend on the values in the config struct but are
   * not related to CIP descriptors for example. However, if the values in the
   * config structs are replaced with CIP priorities, this function could be
   * be used.
   */
  /*
  int findDescriptor(const OBCisTransStereo::Config &config)
  {
    std::vector<unsigned long> refs1(2), refs2(2);
    refs1[0] = config.refs[0];
    refs1[1] = config.refs[1];
    refs2[0] = config.refs[2];
    refs2[1] = config.refs[3];
    if ((OBStereo::NumInversions(refs1) % 2 + OBStereo::NumInversions(refs2) % 2) % 2)
      return 1;
    else
      return 0;
  }
  */

  /**
   * Find a descriptor for the specified atom.
   */
  /*
  int findDescriptor(OBAtom *center, const std::vector<unsigned int> &symmetry_classes, OBStereoFacade &stereoFacade)
  {
    OBMol *mol = center->GetParent();
    if (!stereoFacade.HasTetrahedralStereo(center->GetId()))
      return -1;
    OBTetrahedralStereo::Config symConfig = stereoFacade.GetTetrahedralStereo(center->GetId())->GetConfig();
    IdsToSymClasses(mol, symConfig, symmetry_classes);
    return findDescriptor(symConfig);
  }
  */

  /**
   * Find a descriptor for the specified bond.
   */
  /*
  int findDescriptor(OBBond *bond, const std::vector<unsigned int> &symmetry_classes, OBStereoFacade &stereoFacade)
  {
    OBMol *mol = bond->GetParent();
    if (!stereoFacade.HasCisTransStereo(bond->GetId()))
      return -1;
    OBCisTransStereo::Config symConfig = stereoFacade.GetCisTransStereo(bond->GetId())->GetConfig();
    IdsToSymClasses(mol, symConfig, symmetry_classes);
    return findDescriptor(symConfig);
  }
  */

  // GENERIC
  /**
   * Find the descriptor vector for all stereocenters in the specified fragment.
   * The individual descriptors will be ordered in the same way are specified in
   * the @p orderedUnits parameter. This method works by replacing the unique
   * atom ids in the config structs by symmetry classes to obtain canonical
   * results. This only works for resolved (T1234 & C12) stereogenic units and
   * unresolved units are excluded from the vector. The elements in the returned
   * vector are 0 or 1 and obtained using the findDescriptor functions above.
   */
  /*
  std::vector<int> findDescriptorVector(OBMol *mol, const OBBitVec &fragment,
      const OBStereoUnitSet &orderedUnits, const std::vector<unsigned int> &symmetry_classes, OBStereoFacade &stereoFacade)
  {
    OBStereoUnitSet units = removeUnspecifiedUnits(mol, orderedUnits, stereoFacade);
    std::vector<int> v;
    for (unsigned int i = 0; i < units.size(); ++i) {
      const OBStereoUnit &unit = units[i];

      bool isInFragment = false;
      if (unit.type == OBStereo::Tetrahedral) {
        if (fragment.BitIsOn(mol->GetAtomById(unit.id)->GetIdx()))
          isInFragment = true;
      } else if(unit.type == OBStereo::CisTrans) {
        OBBond *bond = mol->GetBondById(unit.id);
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();
        if (fragment.BitIsOn(begin->GetIdx()) || fragment.BitIsOn(end->GetIdx()))
          isInFragment = true;
      }
      if (!isInFragment)
        continue;

      v.push_back(findDescriptor(mol, unit, symmetry_classes, stereoFacade));
    }

    return v;
  }
  */

  /**
   * Find the value for a descriptor vector (see findDescriptorVector()). A
   * descriptor vectors contains 0 or 1 for each stereogenic unit in the fragment
   * in the same order specified in @p orderedUnits. This function reads these
   * vectors as binary numbers and returns the decimal value. A single value
   * for a descriptor vector, makes it easier to compare these vectors.
   */
  /*
  int findDescriptorVectorValue(OBMol *mol, const OBBitVec &fragment,
      const OBStereoUnitSet &orderedUnits, const std::vector<unsigned int> &symmetry_classes, OBStereoFacade &stereoFacade)
  {
    std::vector<int> v = findDescriptorVector(mol, fragment, orderedUnits, symmetry_classes, stereoFacade);

    //cout << "dv_ref = ";
    int value = 0;
    for (unsigned int i = 0; i < v.size(); ++i) {
      //cout << v[i];
      if (!v[i])
        continue;
      int power = v.size() - i - 1;
      // bit shift is equivalent to 2^power
      // i.e., 1 << 0 == 1, 1 << 1 == 2, etc.
      value += 1 << power;
    }
    //cout << endl;

    return value;
  }
  */

  /**
   * Get the symmetry class for the specified stereogenic unit. For tetrahedral
   * stereogenic units, the symmetry class is simply the symmetry class for the
   * center atom. For double bond stereogenic units, the lowest symmetry class
   * for the begin and end atom is returned.
   */
  unsigned int getSymClass(const OBStereoUnit &unit, OBMol *mol, const std::vector<unsigned int> &symmetry_classes)
  {
    if (unit.type == OBStereo::Tetrahedral) {
      OBAtom *center = mol->GetAtomById(unit.id);
      return symmetry_classes[center->GetIndex()];
    } else
    if (unit.type == OBStereo::CisTrans) {
      OBBond *bond = mol->GetBondById(unit.id);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();
      return std::min(symmetry_classes[begin->GetIndex()], symmetry_classes[end->GetIndex()]);
    }
    return OBGraphSym::NoSymmetryClass; // FIXME
  }


  template<int StereoType>
  bool isSpecifiedTemplate(OBMol *mol, const OBStereoUnit &unit, OBStereoFacade &stereoFacade)
  {
    typedef OBStereoTypeTraits<StereoType> Traits;

    if (stereoFacade.HasStereo<StereoType>(unit.id))
      return stereoFacade.GetStereo<typename Traits::Stereo>(unit.id)->GetConfig().specified;

    return false;
  }

  bool isSpecified(OBMol *mol, const OBStereoUnit &unit, OBStereoFacade &stereoFacade)
  {
    switch (unit.type) {
      case OBStereo::Tetrahedral:
        return isSpecifiedTemplate<OBStereo::Tetrahedral>(mol, unit, stereoFacade);
      case OBStereo::CisTrans:
        return isSpecifiedTemplate<OBStereo::CisTrans>(mol, unit, stereoFacade);
    }
    return false;
  }

  OBStereoUnitSet removeUnspecifiedUnits(OBMol *mol, const OBStereoUnitSet &units, OBStereoFacade &stereoFacade)
  {
    OBStereoUnitSet result;
    for (std::size_t i = 0; i < units.size(); ++i)
      if (isSpecified(mol, units[i], stereoFacade))
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
