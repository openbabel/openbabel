#include <openbabel/stereo/squareplanar.h>
#include <cassert>

namespace OpenBabel {

  OBSquarePlanarStereo::OBSquarePlanarStereo(OBMol *mol) : OBTetraPlanarStereo(mol), m_center(OBStereo::NoId)
  {
  }

  OBSquarePlanarStereo::~OBSquarePlanarStereo()
  {
  }

  // testIsValid
  bool OBSquarePlanarStereo::IsValid() const
  {
    if (m_center == OBStereo::NoId)
      return false;
    if (m_refs.size() != 4)
      return false;
    return true;
  }

  // testCenter
  void OBSquarePlanarStereo::SetCenter(unsigned long id)
  {
    m_center = id; 
  }
    
  // testCenter
  unsigned long OBSquarePlanarStereo::GetCenter() const
  {
    return m_center;
  }

  // testRefs1
  void OBSquarePlanarStereo::SetRefs(const std::vector<unsigned long> &refs, 
      OBStereo::Shape shape)
  {
    assert( refs.size() == 4);
    m_refs = OBTetraPlanarStereo::ToInternal(refs, shape);
  }
  
  // testRefs1
  std::vector<unsigned long> OBSquarePlanarStereo::GetRefs(OBStereo::Shape shape) const
  {
    if (m_refs.empty())
      return m_refs;
    return OBTetraPlanarStereo::ToShape(m_refs, shape);
  }
  
  std::vector<unsigned long> OBSquarePlanarStereo::GetRefs(unsigned long start, 
      OBStereo::Shape shape) const
  {
    if (m_refs.empty())
      return m_refs;

    std::vector<unsigned long> refs(m_refs);
    // since m_refs are U shaped we can rotate the refs lexicographically
    for (int i = 0; i < 4; ++i) {
      std::rotate(refs.begin(), refs.begin() + 2, refs.end());
      if (refs.at(0) == start)
        return OBTetraPlanarStereo::ToShape(refs, shape);
    }

    // start not found, return empty sequence
    refs.clear();
    return refs;
  }

  // testCisTrans
  bool OBSquarePlanarStereo::IsTrans(unsigned long id1, unsigned long id2) const
  {
    return (GetTransRef(id1) == id2);
  }

  // testCisTrans
  bool OBSquarePlanarStereo::IsCis(unsigned long id1, unsigned long id2) const
  {
    if (m_refs.size() != 4)
      return false;

    std::vector<unsigned long> cis = GetCisRefs(id1);
    if (cis.size() != 2)
      return false;

    if ((cis.at(0) == id2) || (cis.at(1) == id2))
      return true;

    return false;
  }
 
  // testCisTrans  
  unsigned long OBSquarePlanarStereo::GetTransRef(unsigned long id) const
  {
    if (m_refs.size() != 4)
      return false;

    // find id1
    for (int i = 0; i < 4; ++i) {
      if (m_refs.at(i) == id) {
        // use it's index to compare id2 with the opposite reference id
        int j = (i > 1) ? i - 2 : i + 2;
        return m_refs.at(j);
      }
    }

    // id not found
    return OBStereo::NoId;
  }

  // testCisTrans  
  std::vector<unsigned long> OBSquarePlanarStereo::GetCisRefs(unsigned long id) const
  {
    std::vector<unsigned long> refs;
    if (m_refs.size() != 4)
      return refs;

    // find id
    for (int i = 0; i < 4; ++i) {
      if (m_refs.at(i) == id) {
        // use it's index to get the left/right reference ids
        int j = (i > 0) ? i - 1 : 3;
        int k = (i < 3) ? i + 1 : 0;
        refs.push_back(m_refs.at(j));
        refs.push_back(m_refs.at(k));
        return refs;
      }
    }

    // id not found
    return refs;  
  }
    
  // testCompare    
  bool OBSquarePlanarStereo::Compare(const std::vector<unsigned long> &refs, OBStereo::Shape shape) const
  {
    // 1-2, 3-4
    // 1-3, 2-4
    // 1-4, 2-3
  
    if (!IsValid() || (refs.size() != 4))
      return false;

    std::vector<unsigned long> u = OBTetraPlanarStereo::ToInternal(refs, shape);
    unsigned long a1 = u.at(0);
    unsigned long b1 = u.at(2);

    unsigned long a2 = GetTransRef(b1);
    if (a1 == a2)
      return true;

    return false;
  }


}
