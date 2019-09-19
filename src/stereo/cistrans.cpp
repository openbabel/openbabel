#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/oberror.h>

using namespace std;

namespace OpenBabel {

  //
  // OBCisTransStereo::Config struct
  //

  bool OBCisTransStereo::Config::operator==(const Config &other) const
  {
    if ((begin != other.begin) && (begin != other.end))
      return false;
    if ((end != other.begin) && (end != other.end))
      return false;
    if ((refs.size() != 4) || (other.refs.size() != 4))
      return false;

    Config u1, u2;
    if (!OBStereo::ContainsSameRefs(refs, other.refs)) {
      // find a ref that occurs in both
      for (OBStereo::ConstRefIter i = refs.begin(); i != refs.end(); ++i)
        if (OBStereo::ContainsRef(other.refs, *i)) {
          u1 = OBTetraPlanarStereo::ToConfig(*this, *i, OBStereo::ShapeU); // refs[0] = u1.refs[0]
          u2 = OBTetraPlanarStereo::ToConfig(other, *i, OBStereo::ShapeU); // refs[0] = u2.refs[0]
        }

      // check if they actualy share an id...
      if (u1.refs.empty())
        return false;
    } else {
      // normalize the other Config struct
      u1 = OBTetraPlanarStereo::ToConfig(*this, refs.at(0), OBStereo::ShapeU); // refs[0] = u1.refs[0]
      u2 = OBTetraPlanarStereo::ToConfig(other, refs.at(0), OBStereo::ShapeU); // refs[0] = u2.refs[0]
      // both now start with the same ref
      //
      // 2 possiblilities:
      //
      //   1 2 3 4      1 2 3 4
      //   |   |        |   |      <- in any case, refs[0] & refs[2] remain unchanged
      //   1 2 3 4      1 4 3 2
      //
      return (u1.refs[2] == u2.refs[2]);
    }

    // possibilities:
    //
    //   1 2 3 4
    //   |   |      <- refs[0] & refs[2] remain unchanged
    //   1 H 3 H
    //
    //   1 2 3 4
    //   |     |    <- refs[0] & refs[3] remain unchanged
    //   1 H H 4
    //
    //   1 2 3 4
    //   | |        <- refs[0] & refs[1] remain unchanged
    //   1 2 H H
    if ((u1.refs[2] == OBStereo::ImplicitRef) || (u2.refs[2] == OBStereo::ImplicitRef)) {
      // 1 2 H 4
      if ((u1.refs[3] == OBStereo::ImplicitRef) || (u2.refs[3] == OBStereo::ImplicitRef)) {
        return (u1.refs[1] == u2.refs[1]); // 1 2 H H
      } else {
        return (u1.refs[3] == u2.refs[3]); // 1 H H 4
      }
    } else
      return (u1.refs[2] == u2.refs[2]); // 1 2 3 4  &  1 H 3 4  &  1 2 3 H

    return false;
  }

  //
  // OBCisTransStereo class
  //

  OBCisTransStereo::OBCisTransStereo(OBMol *mol) : OBTetraPlanarStereo(mol)
  {
  }

  OBCisTransStereo::~OBCisTransStereo()
  {
  }

  bool OBCisTransStereo::IsValid() const
  {
    if ((m_cfg.begin == OBStereo::NoRef) || (m_cfg.end == OBStereo::NoRef))
      return false;
    if (m_cfg.refs.size() != 4)
      return false;
    return true;
  }

  void OBCisTransStereo::SetConfig(const Config &config)
  {
    if (config.begin == OBStereo::NoRef) {
      obErrorLog.ThrowError(__FUNCTION__,
          "OBCisTransStereo::SetConfig : double bond begin id is invalid.", obError);
      m_cfg = Config();
      return;
    }
    if (config.end == OBStereo::NoRef) {
      obErrorLog.ThrowError(__FUNCTION__,
          "OBCisTransStereo::SetConfig : double bond end id is invalid.", obError);
      m_cfg = Config();
      return;
    }
    if (config.refs.size() != 4) {
      std::stringstream ss;
      ss << "OBCisTransStereo::SetConfig : found " << config.refs.size();
      ss << " reference ids, should be 4.";
      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
      m_cfg = Config();
      return;
    }

    // store using U shape
    m_cfg = OBTetraPlanarStereo::ToConfig(config, config.refs.at(0), OBStereo::ShapeU);
  }

  OBCisTransStereo::Config OBCisTransStereo::GetConfig(OBStereo::Shape shape) const
  {
    if (!IsValid())
      return Config();

    return OBTetraPlanarStereo::ToConfig(m_cfg, m_cfg.refs.at(0), shape);
  }

  OBCisTransStereo::Config OBCisTransStereo::GetConfig(unsigned long start,
      OBStereo::Shape shape) const
  {
    if (!IsValid())
      return Config();

    return OBTetraPlanarStereo::ToConfig(m_cfg, start, shape);
  }

  bool OBCisTransStereo::IsTrans(unsigned long id1, unsigned long id2) const
  {
    return (GetTransRef(id1) == id2);
  }

  bool OBCisTransStereo::IsCis(unsigned long id1, unsigned long id2) const
  {
    return (GetCisRef(id1) == id2);
  }

  bool OBCisTransStereo::operator==(const OBCisTransStereo &other) const
  {
    if (!IsValid() || !other.IsValid())
      return false;

    Config u = OBTetraPlanarStereo::ToConfig(other.GetConfig(),
        m_cfg.refs.at(0), OBStereo::ShapeU);
    unsigned long a1 = u.refs.at(0);
    unsigned long b1 = u.refs.at(2);

    if ((a1 == OBStereo::ImplicitRef) && (b1 == OBStereo::ImplicitRef)) {
      a1 = u.refs.at(1);
      b1 = u.refs.at(3);
    }

    if (b1 != OBStereo::ImplicitRef)
      if (a1 == GetTransRef(b1))
        return true;
    if (a1 != OBStereo::ImplicitRef)
      if (b1 == GetTransRef(a1))
        return true;

    return false;
  }

  unsigned long OBCisTransStereo::GetCisOrTransRef(unsigned long id, bool getcisref) const
  {
    if (!IsValid())
      return OBStereo::NoRef;

    if (id == OBStereo::ImplicitRef)
      return OBStereo::NoRef;

    // find id
    for (int i = 0; i < 4; ++i) {
      if (m_cfg.refs.at(i) == id) {
        // Use its index to find the index of the cis (or trans) atom
        int j;
        if (getcisref) // GetCisRef
          j = 3 - i; // Convert 0 to 3, and 3 to 0
        else // GetTransRef
          j = (i > 1) ? i - 2 : i + 2;

        unsigned long refId = m_cfg.refs.at(j);
        return refId;
      }
    }

    // id not found
    return OBStereo::NoRef;
  }

  unsigned long OBCisTransStereo::GetTransRef(unsigned long id) const
  {
    return GetCisOrTransRef(id, false);
  }

  unsigned long OBCisTransStereo::GetCisRef(unsigned long id) const
  {
    return GetCisOrTransRef(id, true);
  }

  bool OBCisTransStereo::IsOnSameAtom(unsigned long id1, unsigned long id2) const
  {
    const OBMol *mol = GetMolecule();
    if (!mol) {
      obErrorLog.ThrowError(__FUNCTION__, "OBCisTransStereo::IsOnSameAtom : No valid molecule set", obError);
      return false;
    }

    OBAtom *begin = mol->GetAtomById(m_cfg.begin);
    if (!begin) {
      obErrorLog.ThrowError(__FUNCTION__, "OBCisTransStereo::IsOnSameAtom : Begin reference id is not valid.", obError);
      return false;
    }
    OBAtom *end = mol->GetAtomById(m_cfg.end);
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
          if (begin->GetExplicitDegree() == 2)
            return true;
          // check if the end atom really is missing an atom
          if (end->GetExplicitDegree() != 2) {
            obErrorLog.ThrowError(__FUNCTION__,
                "OBCisTransStereo::IsOnSameAtom : id2 is not valid and is not a missing hydrogen.", obError);
            return false;
          }
          // inform user we are treating id2 as deleted hydrogen
          obErrorLog.ThrowError(__FUNCTION__,
              "OBCisTransStereo::IsOnSameAtom : Atom with id2 doesn't exist anymore, must be a (deleted) hydrogen.", obInfo);
        } else if (a->IsConnected(end)) {
          // a is connected to end. again, if this is the atom missing a hydrogen, return false
          if (end->GetExplicitDegree() == 2)
            return true;
          // check if the begin atom really is missing an atom
          if (begin->GetExplicitDegree() != 2) {
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
          if (begin->GetExplicitDegree() == 2)
            return true;
          // check if the end atom really is missing an atom
          if (end->GetExplicitDegree() != 2) {
            obErrorLog.ThrowError(__FUNCTION__,
                "OBCisTransStereo::IsOnSameAtom : id1 is not valid and is not a missing hydrogen.", obError);
            return true;
          }
          // inform user we are treating id1 as deleted hydrogen
          obErrorLog.ThrowError(__FUNCTION__,
              "OBCisTransStereo::IsOnSameAtom : Atom with id1 doesn't exist, must be a (deleted) hydrogen.", obInfo);
        } else if (b->IsConnected(end)) {
          // a is connected to end. again, if this is the atom missing a hydrogen, return false
          if (end->GetExplicitDegree() == 2)
            return true;
          // check if the begin atom really is missing an atom
          if (begin->GetExplicitDegree() != 2) {
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
          if ((m_cfg.refs.at(i) == id1) || (m_cfg.refs.at(i) == id2))
            continue;
          if (!c) {
            c = mol->GetAtomById(m_cfg.refs.at(i));
          } else {
            d = mol->GetAtomById(m_cfg.refs.at(i));
          }
        }
        if (!c || !d) {
          obErrorLog.ThrowError(__FUNCTION__, "OBCisTransStereo::IsOnSameAtom : invalid stereochemistry!", obError);
          return true;
        }
        if ((begin->GetExplicitDegree() != 2) || (end->GetExplicitDegree() != 2)) {
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

  OBGenericData* OBCisTransStereo::Clone(OBBase *mol) const
  {
    OBCisTransStereo *data = new OBCisTransStereo(static_cast<OBMol*>(mol));
    data->SetConfig(m_cfg);
    return data;
  }

} // namespace OpenBabel

namespace std {

  ostream& operator<<(ostream &out, const OpenBabel::OBCisTransStereo &ct)
  {
    OpenBabel::OBCisTransStereo::Config cfg = ct.GetConfig();
    out << "OBCisTransStereo(begin = " << cfg.begin;
    out << ", end = " << cfg.end;

    out << ", refs = ";
    for (OpenBabel::OBStereo::Refs::iterator i = cfg.refs.begin(); i != cfg.refs.end(); ++i)
      if (*i != OpenBabel::OBStereo::ImplicitRef)
        out << *i << " ";
      else
        out << "H ";

    switch (cfg.shape) {
      case OpenBabel::OBStereo::ShapeU:
        out << ", shape = U)";
        break;
      case OpenBabel::OBStereo::ShapeZ:
        out << ", shape = Z)";
        break;
      case OpenBabel::OBStereo::Shape4:
        out << ", shape = 4)";
        break;
    }

    return out;
  }

  ostream& operator<<(ostream &out, const OpenBabel::OBCisTransStereo::Config &cfg)
  {
    out << "OBCisTransStereo::Config(begin = " << cfg.begin;
    out << ", end = " << cfg.end;

    out << ", refs = ";
    for (OpenBabel::OBStereo::Refs::const_iterator i = cfg.refs.begin(); i != cfg.refs.end(); ++i)
      if (*i != OpenBabel::OBStereo::ImplicitRef)
        out << *i << " ";
      else
        out << "H ";

    switch (cfg.shape) {
      case OpenBabel::OBStereo::ShapeU:
        out << ", shape = U)";
        break;
      case OpenBabel::OBStereo::ShapeZ:
        out << ", shape = Z)";
        break;
      case OpenBabel::OBStereo::Shape4:
        out << ", shape = 4)";
        break;
    }

    return out;
  }

} // namespace std

