#include <openbabel/datacache.h>
#include <openbabel/graphsym.h>
#include <openbabel/isomorphism.h>
#include <openbabel/stereo/stereo.h>

namespace OpenBabel {

  typedef OBPairTemplate< std::vector<unsigned int> > OBVectorUIntData;

  DataCache::DataCache(OBMol *mol) :
      m_validated(false), m_mol(mol)
  {
  }

  template<typename T>
  void SetData(OBMol *mol, const std::string &attribute, const T &value)
  {
    OBPairTemplate<T> *data = new OBPairTemplate<T>;
    data->SetAttribute(attribute);
    data->SetValue(value);
    mol->SetData(data);
  }

  template<typename T>
  OBPairTemplate<T>* GetData(OBMol *mol, const std::string &attribute)
  {
    return dynamic_cast<OBPairTemplate<T>*>(mol->GetData(attribute));
  }

  void DataCache::Validate()
  {
    bool isValid = false;
    OBVectorUIntData *data = GetData< std::vector<unsigned int> >(m_mol, "GraphInvariants");

    std::vector<unsigned int> invariants;
    m_mol->GetGIVector(invariants);

    if (data) {
      if (invariants == data->GetGenericValue())
        isValid = true;
    }

    if (!isValid) {
      m_mol->DeleteData("GraphInvariants");
      m_mol->DeleteData("SymmetryClasses");
      m_mol->DeleteData("Automorphisms");

      SetData(m_mol, "GraphInvariants", invariants);
    }

    m_validated = true;
  }

  void DataCache::GetSymmetryClasses(std::vector<unsigned int> &symmetry_classes, const OBBitVec &mask, bool optimizeMask)
  {
    if (!m_validated)
      Validate();

    OBBitVec maskCopy;
    if (optimizeMask)
      RemoveHydrogensFromMask(maskCopy, mask);
    else
      maskCopy = mask;

    OBVectorUIntData *data = GetData< std::vector<unsigned int> >(m_mol, "SymmetryClasses");

    if (data) {
      OBPairTemplate<OBBitVec> *maskData = GetData<OBBitVec>(m_mol, "SymmetryClassesMask");
      if (maskData && maskData->GetGenericValue() == maskCopy) {
        symmetry_classes = data->GetGenericValue();
        return;
      }
    }

    if (!maskCopy.IsEmpty()) {
      OBGraphSym gs(m_mol, &maskCopy);
      gs.GetSymmetry(symmetry_classes);
    } else {
      OBGraphSym gs(m_mol);
      gs.GetSymmetry(symmetry_classes);
    }

    m_mol->DeleteData("SymmetryClasses");
    m_mol->DeleteData("SymmetryClassesMask");

    SetData(m_mol, "SymmetryClasses", symmetry_classes);
    SetData(m_mol, "SymmetryClassesMask", maskCopy);
  }

  void DataCache::GetAutomorphisms(Automorphisms &automorphisms, const OBBitVec &mask)
  {
    if (!m_validated)
      Validate();

    OBPairTemplate<Automorphisms> *data = GetData<Automorphisms>(m_mol, "Automorphisms");

    if (data) {
      OBPairTemplate<OBBitVec> *maskData = GetData<OBBitVec>(m_mol, "AutomorphismsMask");
      if (maskData && maskData->GetGenericValue() == mask) {
        automorphisms = data->GetGenericValue();
        return;
      }
    }

    std::vector<unsigned int> symmetry_classes;
    GetSymmetryClasses(symmetry_classes, mask);

    FindAutomorphisms(m_mol, automorphisms, symmetry_classes, mask);

    m_mol->DeleteData("Automorphisms");
    m_mol->DeleteData("AutomorphismsMask");

    SetData(m_mol, "Automorphisms", automorphisms);
    SetData(m_mol, "AutomorphismsMask", mask);
  }

  void DataCache::GetStereogenicUnits(OBStereoUnitSet &units, const OBBitVec &mask, bool optimizeMask)
  {
    if (!m_validated)
      Validate();

    OBBitVec stereomask;
    OptimizeStereoMask(stereomask, mask);

    OBPairTemplate<OBStereoUnitSet> *data = GetData<OBStereoUnitSet>(m_mol, "StereogenicUnits");

    if (data) {
      OBPairTemplate<OBBitVec> *maskData = GetData<OBBitVec>(m_mol, "StereogenicUnitsMask");
      if (maskData && maskData->GetGenericValue() == stereomask) {
        units = data->GetGenericValue();
        return;
      }
    }

    std::vector<unsigned int> symmetry_classes;
    GetSymmetryClasses(symmetry_classes, mask);

    Automorphisms automorphisms;
    GetAutomorphisms(automorphisms, stereomask);

    units = FindStereogenicUnits(m_mol, symmetry_classes, automorphisms);

    m_mol->DeleteData("StereogenicUnits");
    m_mol->DeleteData("StereogenicUnitsMask");

    SetData(m_mol, "StereogenicUnits", units);
    SetData(m_mol, "StereogenicUnitsMask", stereomask);
  }

  void DataCache::RemoveHydrogensFromMask(OBBitVec &output, const OBBitVec &input)
  {
    FOR_ATOMS_OF_MOL (atom, m_mol) {
      if (!atom->IsHydrogen())
        output.SetBitOn(atom->GetIdx());
    }

    if (!input.IsEmpty())
      output &= input;
  }

  void DataCache::OptimizeStereoMask(OBBitVec &output, const OBBitVec &input)
  {
    FOR_ATOMS_OF_MOL (atom, m_mol) {
      if (!atom->IsHydrogen() && atom->GetValence() < 7)
        output.SetBitOn(atom->GetIdx());
    }

    FOR_ATOMS_OF_MOL (atom, m_mol) {
      std::vector<OBAtom*> fluorines, chlorines, bromines, iodines;

      int heavyValence = 0;
      FOR_NBORS_OF_ATOM (nbr, &*atom) {
        if (!atom->IsHydrogen())
          heavyValence++;

        switch (nbr->GetAtomicNum()) {
          case 9: // F
            fluorines.push_back(&*nbr);
            break;
          case 17: // Cl
            chlorines.push_back(&*nbr);
            break;
          case 35: // Br
            bromines.push_back(&*nbr);
            break;
          case 53: // I
            iodines.push_back(&*nbr);
            break;
          default:
            break;
        }

        if (heavyValence == 3 || heavyValence == 4) {
          if (fluorines.size() >= 2)
            for (std::size_t i = 0; i < fluorines.size(); ++i)
              output.SetBitOff(fluorines[i]->GetIdx());
          if (chlorines.size() >= 2)
            for (std::size_t i = 0; i < chlorines.size(); ++i)
              output.SetBitOff(chlorines[i]->GetIdx());
          if (bromines.size() >= 2)
            for (std::size_t i = 0; i < bromines.size(); ++i)
              output.SetBitOff(bromines[i]->GetIdx());
          if (iodines.size() >= 2)
            for (std::size_t i = 0; i < iodines.size(); ++i)
              output.SetBitOff(iodines[i]->GetIdx());
        } else if (heavyValence > 4) {
          if (fluorines.size() == heavyValence)
            for (std::size_t i = 0; i < fluorines.size(); ++i)
              output.SetBitOff(fluorines[i]->GetIdx());
          if (chlorines.size() == heavyValence)
            for (std::size_t i = 0; i < chlorines.size(); ++i)
              output.SetBitOff(chlorines[i]->GetIdx());
          if (bromines.size() == heavyValence)
            for (std::size_t i = 0; i < bromines.size(); ++i)
              output.SetBitOff(bromines[i]->GetIdx());
          if (iodines.size() == heavyValence)
            for (std::size_t i = 0; i < iodines.size(); ++i)
              output.SetBitOff(iodines[i]->GetIdx());
        }
      }
    }

    if (!input.IsEmpty())
      output &= input;
  }


}
