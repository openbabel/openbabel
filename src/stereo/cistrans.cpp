#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>
#include <openbabel/oberror.h>
#include <cassert>

using namespace std;

namespace OpenBabel {

  OBCisTransStereo::OBCisTransStereo(OBMol *mol) : OBTetraPlanarStereo(mol), 
      m_begin(OBStereo::NoId), m_end(OBStereo::NoId)
  {
  }

  OBCisTransStereo::~OBCisTransStereo()
  {
  
  }

  bool OBCisTransStereo::IsValid() const
  {
    if ((m_begin == OBStereo::NoId) || (m_end == OBStereo::NoId))
      return false;
    if (m_refs.size() != 4)
      return false;
    return true;
  }

  void OBCisTransStereo::SetCenters(unsigned long begin, unsigned long end)
  {
    m_begin = begin;
    m_end = end;
  }
    
  unsigned long OBCisTransStereo::GetBegin() const
  {
    return m_begin;
  }
 
  void OBCisTransStereo::SetBegin(unsigned long id)
  {
    m_begin = id;
  }
  
  unsigned long OBCisTransStereo::GetEnd() const
  {
    return m_end;
  }

  void OBCisTransStereo::SetEnd(unsigned long id)
  {
    m_end = id;
  }
 
  void OBCisTransStereo::SetRefs(const std::vector<unsigned long> &refs, 
      OBStereo::Shape shape)
  {
    assert( refs.size() == 4 );
    m_refs = OBTetraPlanarStereo::ToInternal(refs, shape);
  }
  
  std::vector<unsigned long> OBCisTransStereo::GetRefs(OBStereo::Shape shape) const
  {
    if (m_refs.empty())
      return m_refs;
    return OBTetraPlanarStereo::ToShape(m_refs, shape);
  }
  
  bool OBCisTransStereo::IsTrans(unsigned long id1, unsigned long id2) const
  {
    return (GetTransRef(id1) == id2);
  }
  
  bool OBCisTransStereo::IsCis(unsigned long id1, unsigned long id2) const
  {
    return (GetCisRef(id1) == id2);
  }
 
  bool OBCisTransStereo::Compare(const std::vector<unsigned long> &refs, OBStereo::Shape shape) const
  {
    if (!IsValid() || (refs.size() != 4))
      return false;

    std::vector<unsigned long> u = OBTetraPlanarStereo::ToInternal(refs, shape);
    unsigned long a1 = u.at(0);
    unsigned long b1 = u.at(2);

    if ((a1 == OBStereo::HydrogenId) && (b1 == OBStereo::HydrogenId)) {
      a1 = u.at(1);
      b1 = u.at(3);
    }

    if (b1 != OBStereo::HydrogenId)
      if (a1 == GetTransRef(b1))
        return true;
    if (a1 != OBStereo::HydrogenId)
      if (b1 == GetTransRef(a1))
        return true;

    return false;
  }

  unsigned long OBCisTransStereo::GetTransRef(unsigned long id) const
  {
    if (!IsValid())
      return OBStereo::NoId;

    if (id == OBStereo::HydrogenId)
      return OBStereo::NoId;

    // find id1
    for (int i = 0; i < 4; ++i) {
      if (m_refs.at(i) == id) {
        // use it's index to compare id2 with the opposite reference id
        int j = (i > 1) ? i - 2 : i + 2;
        // make sure they are not bonded to the same atom
        unsigned long transId = m_refs.at(j);
        if (transId == OBStereo::HydrogenId)
          return OBStereo::HydrogenId;
        if (IsOnSameAtom(id, transId)) {
          obErrorLog.ThrowError(__FUNCTION__, 
              "OBCisTransStereo::GetTransRef : References don't match bond orientation", obError);
          return OBStereo::NoId;
        }
        return transId;
      }
    }

    // id not found
    return OBStereo::NoId;
  }

  unsigned long OBCisTransStereo::GetCisRef(unsigned long id) const
  {
    if (!IsValid())
      return OBStereo::NoId;

    if (id == OBStereo::HydrogenId)
      return OBStereo::NoId;

    // find id
    for (int i = 0; i < 4; ++i) {
      if (m_refs.at(i) == id) {
        // use it's index to get the left/right reference ids
        int j = (i > 0) ? i - 1 : 3;
        int k = (i < 3) ? i + 1 : 0;
        // make sure they are not bonded to the same atom
        if (m_refs.at(j) != OBStereo::HydrogenId)
          if (!IsOnSameAtom(id, m_refs.at(j)))
            return m_refs.at(j);
        if (m_refs.at(k) != OBStereo::HydrogenId)
          if (!IsOnSameAtom(id, m_refs.at(k)))
            return m_refs.at(k);

        if ((m_refs.at(j) == OBStereo::HydrogenId) && (m_refs.at(k) == OBStereo::HydrogenId)) {
          return OBStereo::HydrogenId;       
        }

        obErrorLog.ThrowError(__FUNCTION__, 
            "OBCisTransStereo::GetTransRef : References don't match bond orientation", obError);
        return OBStereo::NoId;
      }
    }

    // id not found
    return OBStereo::NoId;
  }
    
  bool OBCisTransStereo::IsOnSameAtom(unsigned long id1, unsigned long id2) const
  {
    const OBMol *mol = GetMolecule();
    if (!mol) {
      obErrorLog.ThrowError(__FUNCTION__, "OBCisTransStereo::IsOnSameAtom : No valid molecule set", obError);
      return false;
    }

    OBAtom *begin = mol->GetAtomById(m_begin);
    if (!begin) {
      obErrorLog.ThrowError(__FUNCTION__, "OBCisTransStereo::IsOnSameAtom : Begin reference id is not valid.", obError);
      return false;
    }
    OBAtom *end = mol->GetAtomById(m_end);
    if (!end) {
      obErrorLog.ThrowError(__FUNCTION__, "OBCisTransStereo::IsOnSameAtom : End reference id is not valid.", obError);
      return false;
    }

    OBAtom *a = mol->GetAtomById(id1);
    OBAtom *b = mol->GetAtomById(id2);
    
    if (a && b) {
      // both on begin atom?
      if (a->IsConnected(begin) && b->IsConnected(begin))
          return true;
      // both on end atom?
      if (a->IsConnected(end) && b->IsConnected(end))
          return true;
      return false;
    } else {
      if (a) {
        // b atom not found, could be a deleted hydrogen...
        if (a->IsConnected(begin)) {
          // a is connected to begin. if this is the atom missing a hydrogen, return false
          if (begin->GetValence() == 2)
            return true;
          // check if the end atom really is missing an atom
          if (end->GetValence() != 2) {
            obErrorLog.ThrowError(__FUNCTION__, 
                "OBCisTransStereo::IsOnSameAtom : id2 is not valid and is not a missing hydrogen.", obError);
            return false;
          }
          // inform user we are treating id2 as deleted hydrogen 
          obErrorLog.ThrowError(__FUNCTION__, 
              "OBCisTransStereo::IsOnSameAtom : Atom with id2 doesn't exist anymore, must be a (deleted) hydrogen.", obInfo);
        } else if (a->IsConnected(end)) {
          // a is connected to end. again, if this is the atom missing a hydrogen, return false
          if (end->GetValence() == 2)
            return true;
          // check if the begin atom really is missing an atom
          if (begin->GetValence() != 2) {
            obErrorLog.ThrowError(__FUNCTION__, 
                "OBCisTransStereo::IsOnSameAtom : id2 is not valid and is not a missing hydrogen.", obError);
            return true;
          }
          // inform user we are treating id2 as deleted hydrogen 
          obErrorLog.ThrowError(__FUNCTION__, 
              "OBCisTransStereo::IsOnSameAtom : Atom with id2 doesn't exist, must be a (deleted) hydrogen.", obInfo);
 
        } else {
          obErrorLog.ThrowError(__FUNCTION__, 
              "OBCisTransStereo::IsOnSameAtom : Atom with id1 isn't connected to the begin or end atom.", obError);
          return true;      
        }
      } else if (b) {
        // a atom not found, could be a deleted hydrogen...
        if (b->IsConnected(begin)) {
          // b is connected to begin. if this is the atom missing a hydrogen, return false
          if (begin->GetValence() == 2)
            return true;
          // check if the end atom really is missing an atom
          if (end->GetValence() != 2) {
            obErrorLog.ThrowError(__FUNCTION__, 
                "OBCisTransStereo::IsOnSameAtom : id1 is not valid and is not a missing hydrogen.", obError);
            return true;
          }
          // inform user we are treating id1 as deleted hydrogen 
          obErrorLog.ThrowError(__FUNCTION__, 
              "OBCisTransStereo::IsOnSameAtom : Atom with id1 doesn't exist, must be a (deleted) hydrogen.", obInfo);
        } else if (b->IsConnected(end)) {
          // a is connected to end. again, if this is the atom missing a hydrogen, return false
          if (end->GetValence() == 2)
            return true;
          // check if the begin atom really is missing an atom
          if (begin->GetValence() != 2) {
            obErrorLog.ThrowError(__FUNCTION__, 
                "OBCisTransStereo::IsOnSameAtom : id1 is not valid and is not a missing hydrogen.", obError);
            return true;
          }
          // inform user we are treating id2 as deleted hydrogen 
          obErrorLog.ThrowError(__FUNCTION__, 
              "OBCisTransStereo::IsOnSameAtom : Atom with id1 doesn't exist, must be a (deleted) hydrogen.", obInfo);
        } else {
          obErrorLog.ThrowError(__FUNCTION__, 
              "OBCisTransStereo::IsOnSameAtom : Atom with id1 isn't connected to the begin or end atom.", obError);
          return true;      
        }
      } else {
        OBAtom *c = 0, *d = 0;
        // no a & b, check the remaining ids which will reveal same info
        for (int i = 0; i < 4; ++i) {
          if ((m_refs.at(i) == id1) || (m_refs.at(i) == id2))
            continue;
          if (!c) {
            c = mol->GetAtomById(m_refs.at(i));
          } else {
            d = mol->GetAtomById(m_refs.at(i));
          }
        }
        if (!c || !d) {
          obErrorLog.ThrowError(__FUNCTION__, "OBCisTransStereo::IsOnSameAtom : invalid stereochemistry!", obError);
          return true;
        }
        if ((begin->GetValence() != 2) || (end->GetValence() != 2)) {
          obErrorLog.ThrowError(__FUNCTION__, "OBCisTransStereo::IsOnSameAtom : invalid stereochemistry!", obError);
          return true;
        }
        obErrorLog.ThrowError(__FUNCTION__, 
            "OBCisTransStereo::IsOnSameAtom : Atoms with id1 & id2 don't exist, must be a (deleted) hydrogens.", obInfo);
        return IsOnSameAtom(c->GetId(), d->GetId());
      }
    }

    return false;  
  }

}
