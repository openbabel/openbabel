#include <openbabel/stereo/squareplanar.h>
#include <openbabel/mol.h>
#include <openbabel/oberror.h>
#include <algorithm> // std::rotate
#include <cassert>

namespace OpenBabel {

  //
  // OBSquarePlanarStereo::Config struct
  //

  bool OBSquarePlanarStereo::Config::operator==(const Config &other) const
  {
    if (center != other.center)
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
  // OBSquarePlanar class
  //

  OBSquarePlanarStereo::OBSquarePlanarStereo(OBMol *mol) : OBTetraPlanarStereo(mol)
  {
  }

  OBSquarePlanarStereo::~OBSquarePlanarStereo()
  {
  }

  bool OBSquarePlanarStereo::IsValid() const
  {
    if (m_cfg.center == OBStereo::NoRef)
      return false;
    if (m_cfg.refs.size() != 4)
      return false;
    return true;
  }

  void OBSquarePlanarStereo::SetConfig(const Config &config)
  {
    if (config.center == OBStereo::NoRef) {
      obErrorLog.ThrowError(__FUNCTION__,
          "OBSquarePlanarStereo::SetConfig : center id is invalid.", obError);
      m_cfg = Config();
      return;
    }
    if (config.refs.size() != 4) {
      std::stringstream ss;
      ss << "OBSquarePlanarStereo::SetConfig : found " << config.refs.size();
      ss << " reference ids, should be 4.";
      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
      m_cfg = Config();
      return;
    }

    // store using U shape
    m_cfg = OBTetraPlanarStereo::ToConfig(config, config.refs.at(0), OBStereo::ShapeU);
  }

  OBSquarePlanarStereo::Config OBSquarePlanarStereo::GetConfig(OBStereo::Shape shape) const
  {
    if (!IsValid())
      return Config();

    return OBTetraPlanarStereo::ToConfig(m_cfg, m_cfg.refs.at(0), shape);
  }

  OBSquarePlanarStereo::Config OBSquarePlanarStereo::GetConfig(unsigned long start,
      OBStereo::Shape shape) const
  {
    if (!IsValid())
      return Config();

    return OBTetraPlanarStereo::ToConfig(m_cfg, start, shape);
  }

  bool OBSquarePlanarStereo::operator==(const OBSquarePlanarStereo &other) const
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

  bool OBSquarePlanarStereo::IsTrans(unsigned long id1, unsigned long id2) const
  {
    return (GetTransRef(id1) == id2);
  }

  bool OBSquarePlanarStereo::IsCis(unsigned long id1, unsigned long id2) const
  {
    if (m_cfg.refs.size() != 4)
      return false;

    std::vector<unsigned long> cis = GetCisRefs(id1);
    if (cis.size() != 2)
      return false;

    if ((cis.at(0) == id2) || (cis.at(1) == id2))
      return true;

    return false;
  }

  unsigned long OBSquarePlanarStereo::GetTransRef(unsigned long id) const
  {
    if (m_cfg.refs.size() != 4)
      return false;

    // find id1
    for (int i = 0; i < 4; ++i) {
      if (m_cfg.refs.at(i) == id) {
        // use it's index to compare id2 with the opposite reference id
        int j = (i > 1) ? i - 2 : i + 2;
        return m_cfg.refs.at(j);
      }
    }

    // id not found
    return OBStereo::NoRef;
  }

  std::vector<unsigned long> OBSquarePlanarStereo::GetCisRefs(unsigned long id) const
  {
    std::vector<unsigned long> refs;
    if (m_cfg.refs.size() != 4)
      return refs;

    // find id
    for (int i = 0; i < 4; ++i) {
      if (m_cfg.refs.at(i) == id) {
        // use it's index to get the left/right reference ids
        int j = (i > 0) ? i - 1 : 3;
        int k = (i < 3) ? i + 1 : 0;
        refs.push_back(m_cfg.refs.at(j));
        refs.push_back(m_cfg.refs.at(k));
        return refs;
      }
    }

    // id not found
    return refs;
  }

  OBGenericData* OBSquarePlanarStereo::Clone(OBBase *mol) const
  {
    OBSquarePlanarStereo *data = new OBSquarePlanarStereo(static_cast<OBMol*>(mol));
    data->SetConfig(m_cfg);
    return data;
  }

} // namespace OpenBabel

namespace std {

  ostream& operator<<(ostream &out, const OpenBabel::OBSquarePlanarStereo &ct)
  {
    OpenBabel::OBSquarePlanarStereo::Config cfg = ct.GetConfig();
    out << "OBSquarePlanarStereo(center = " << cfg.center;

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

  ostream& operator<<(ostream &out, const OpenBabel::OBSquarePlanarStereo::Config &cfg)
  {
    out << "OBSquarePlanarStereo::Config(center = " << cfg.center;

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

