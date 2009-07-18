#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/mol.h>
#include <cassert>

using namespace std;

namespace OpenBabel {

  OBTetrahedralStereo::OBTetrahedralStereo(OBMol *mol) : OBTetraNonPlanarStereo(mol),
      m_center(-1)
  {
  }

  OBTetrahedralStereo::~OBTetrahedralStereo()
  {
  }

  bool OBTetrahedralStereo::IsValid() const
  {
    if (m_center == -1)
      return false;
    if (m_refs.size() != 4)
      return false;
    return true;
  }

  void OBTetrahedralStereo::SetCenter(unsigned long id)
  {
    m_center = id;
  }

  unsigned long OBTetrahedralStereo::GetCenter() const
  {
    return m_center;
  }

  void OBTetrahedralStereo::SetRefs(const std::vector<unsigned long> &refs,
     unsigned long id, OBStereo::Winding winding, OBStereo::View view)
  {
    assert(refs.size() == 3);
    m_refs = OBTetraNonPlanarStereo::ToInternal(refs, id, winding, view);
  }

  std::vector<unsigned long> OBTetrahedralStereo::GetRefs(unsigned long id, 
      OBStereo::Winding winding, OBStereo::View view) const
  {
    if (m_refs.empty())
      return m_refs;
    return OBTetraNonPlanarStereo::ToView(m_refs, id, winding, view);
  }

  bool OBTetrahedralStereo::Compare(const std::vector<unsigned long> &refs, 
      unsigned long id, OBStereo::Winding winding, OBStereo::View view) const
  {
    if (!IsValid() || (refs.size() != 3))
      return false;

    std::vector<unsigned long> clockwiseFrom1 = OBTetraNonPlanarStereo::ToInternal(refs, id, winding, view);
    int Ni1 = OBTetraNonPlanarStereo::NumInversions(refs);
    cout << "Ni1 = " << Ni1 << endl;
    int Ni2 = OBTetraNonPlanarStereo::NumInversions(clockwiseFrom1);
    cout << "Ni2 = " << Ni2 << endl;

    return ((Ni1 + Ni2) % 2 == 0);
    /*
    if (m_refs.at(0) == clockwiseFrom1.at(1)) {
      clockwiseFrom1[1] = clockwiseFrom1.at(2);
      clockwiseFrom1[2] = clockwiseFrom1.at(0);
      clockwiseFrom1[0] = m_refs.at(0);      // use m_refs[0] rather than tmp var since it's got what we need
    } else if (m_refs.at(0) == clockwiseFrom1.at(2)) {
      clockwiseFrom1[2] = clockwiseFrom1.at(1);
      clockwiseFrom1[1] = clockwiseFrom1.at(0);
      clockwiseFrom1[0] = m_refs.at(0);      // ditto re tmp usage
    }

    return (m_refs.at(2) == clockwiseFrom1.at(2));
    */
  }

  OBGenericData* OBTetrahedralStereo::Clone(OBBase *mol) const
  {
    OBTetrahedralStereo *data = new OBTetrahedralStereo(static_cast<OBMol*>(mol));
    data->SetCenter(m_center);
    data->SetRefs(OBStereo::MakeRefs(m_refs[1], m_refs[2], m_refs[3]), m_refs[0]);
    return data;      
  }


}

