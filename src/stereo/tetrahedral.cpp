/**********************************************************************
  tetrahedral.h - OBTetrahedralStereo

  Copyright (C) 2009 by Tim Vandermeersch

  This file is part of the Open Babel project.
  For more information, see <http://openbabel.org/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/mol.h>

namespace OpenBabel {

  //
  // OBTetrahedralStereo::Config struct
  //

  bool OBTetrahedralStereo::Config::operator==(const OBTetrahedralStereo::Config &other) const
  {
    if (center != other.center)
      return false;
    if ((refs.size() != 3) || (other.refs.size() != 3))
      return false;
    // return true if either is unspecified (i.e. accidental)
    if (!specified || !other.specified)
      return true;

    // Convert both Config's refs to same from, winding and view while
    // avoiding having an ImplicitRef in the 'from' position of either
    Config thisConfig, otherConfig;
    if (from == OBStereo::ImplicitRef) {
      thisConfig = OBTetraNonPlanarStereo::ToConfig(*this, refs[0], winding, view);
      otherConfig = OBTetraNonPlanarStereo::ToConfig(other, thisConfig.from, winding, view);
    }
    else if (other.from == OBStereo::ImplicitRef) {
      otherConfig = OBTetraNonPlanarStereo::ToConfig(other, other.refs[0], winding, view);
      thisConfig = OBTetraNonPlanarStereo::ToConfig(*this, otherConfig.from, winding, view);
    }
    else {
      thisConfig = *this;
      otherConfig = OBTetraNonPlanarStereo::ToConfig(other, thisConfig.from, winding, view);
    }

    if (!OBStereo::ContainsSameRefs(thisConfig.refs, otherConfig.refs)) {
      if (OBStereo::ContainsRef(thisConfig.refs, OBStereo::ImplicitRef)) {
        // if both refs already contain ImplicitRef, return false
        if (OBStereo::ContainsRef(otherConfig.refs, OBStereo::ImplicitRef))
          return false;

        // example: *this       = 23H
        //          otherConfig = 234 --> 23H

        // for each ref in otherConfig
        for (unsigned int i = 0; i < otherConfig.refs.size(); ++i) {
          bool found = false;
          for (OBStereo::RefIter j = thisConfig.refs.begin(); j != thisConfig.refs.end(); ++j)
            if (otherConfig.refs.at(i) == *j)
              found = true;

          if (!found) {
            // the ref from otherConfig is not found in this config
            otherConfig.refs[i] = OBStereo::ImplicitRef;
            break;
          }
        }
      } else
      if (OBStereo::ContainsRef(otherConfig.refs, OBStereo::ImplicitRef)) {
        // if both refs already contain ImplicitRef, return false
        if (OBStereo::ContainsRef(thisConfig.refs, OBStereo::ImplicitRef))
          return false;

        // example: *this       = 234
        //          otherConfig = 23H --> 234

        // for each ref in *this
        for (unsigned int i = 0; i < thisConfig.refs.size(); ++i) {
          bool found = false;
          // for each refs in otherConfig
          for (OBStereo::RefIter j = otherConfig.refs.begin(); j != otherConfig.refs.end(); ++j)
            if (thisConfig.refs.at(i) == *j)
              found = true;

          if (!found) {
            for (OBStereo::RefIter j = otherConfig.refs.begin(); j != otherConfig.refs.end(); ++j)
              if (*j == OBStereo::ImplicitRef)
                *j = thisConfig.refs.at(i);
            break;
          }
        }
      }
    }

    int Ni1 = OBStereo::NumInversions(thisConfig.refs);
    int Ni2 = OBStereo::NumInversions(otherConfig.refs);
    return ((Ni1 + Ni2) % 2 == 0);
  }

  //
  // OBTetrahedralStereo class
  //

  OBTetrahedralStereo::OBTetrahedralStereo(OBMol *mol) :
      OBTetraNonPlanarStereo(mol), m_cfg(Config())
  {
  }

  OBTetrahedralStereo::~OBTetrahedralStereo()
  {
  }

  bool OBTetrahedralStereo::IsValid() const
  {
    if (m_cfg.center == OBStereo::NoRef)
      return false;
    if (m_cfg.from == OBStereo::NoRef)
      return false;
    if (m_cfg.refs.size() != 3)
      return false;
    return true;
  }

  void OBTetrahedralStereo::SetConfig(const Config &config)
  {
    if (config.center == OBStereo::NoRef) {
      obErrorLog.ThrowError(__FUNCTION__,
          "OBTetrahedralStereo::SetConfig : center atom id is invalid.", obError);
      m_cfg = Config();
      return;
    }
    if (config.from == OBStereo::NoRef) {
      obErrorLog.ThrowError(__FUNCTION__,
          "OBTetrahedralStereo::SetConfig : from/towards atom id is invalid.", obError);
      m_cfg = Config();
      return;
    }
    if (config.refs.size() != 3) {
      std::stringstream ss;
      ss << "OBTetrahedralStereo::SetConfig : found " << config.refs.size();
      ss << " reference ids, should be 3.";
      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
      m_cfg = Config();
      return;
    }

    m_cfg = config;
  }

  OBTetrahedralStereo::Config OBTetrahedralStereo::GetConfig(
      OBStereo::Winding winding, OBStereo::View view) const
  {
    if (!IsValid())
      return Config();

    if (m_cfg.winding != OBStereo::UnknownWinding)
      return OBTetraNonPlanarStereo::ToConfig(m_cfg, m_cfg.from, winding, view);
    else
      return OBTetraNonPlanarStereo::ToConfig(m_cfg, m_cfg.from, OBStereo::UnknownWinding, view);
  }

  OBTetrahedralStereo::Config OBTetrahedralStereo::GetConfig(unsigned long from_or_towards,
        OBStereo::Winding winding, OBStereo::View view) const
  {
    if (!IsValid())
      return Config();

    if (m_cfg.winding != OBStereo::UnknownWinding)
      return OBTetraNonPlanarStereo::ToConfig(m_cfg, from_or_towards, winding, view);
    else
      return OBTetraNonPlanarStereo::ToConfig(m_cfg, from_or_towards, OBStereo::UnknownWinding, view);
  }

  bool OBTetrahedralStereo::operator==(const OBTetrahedralStereo &other) const
  {
    if (!IsValid() || !other.IsValid())
      return false;

    if (m_cfg == other.GetConfig())
      return true;

    return false;
  }

  OBGenericData* OBTetrahedralStereo::Clone(OBBase *mol) const
  {
    OBTetrahedralStereo *data = new OBTetrahedralStereo(static_cast<OBMol*>(mol));
    data->SetConfig(m_cfg);
    return data;
  }

} // namespace OpenBabel

namespace std {

  using namespace OpenBabel;

  ostream& operator<<(ostream &out, const OBTetrahedralStereo &ts)
  {
    OBTetrahedralStereo::Config cfg = ts.GetConfig();
    out << "OBTetrahedralStereo(center = " << cfg.center;
    if (cfg.view == OBStereo::ViewFrom)
      out << ", viewFrom = ";
    else
      out << ", viewTowards = ";

    if (cfg.from == OBStereo::ImplicitRef)
      out << "H";
    else
      out << cfg.from;

    out << ", refs = ";
    for (OBStereo::Refs::iterator i = cfg.refs.begin(); i != cfg.refs.end(); ++i)
      if (*i != OBStereo::ImplicitRef)
        out << *i << " ";
      else
        out << "H ";

    if (!cfg.specified)
      out << ", unspecified)";
    else {
      if (cfg.winding == OBStereo::Clockwise)
        out << ", clockwise)";
      else
        out << ", anti-clockwise)";
    }

    return out;
  }

  ostream& operator<<(ostream &out, const OBTetrahedralStereo::Config &cfg)
  {
    out << "OBTetrahedralStereo::Config(center = " << cfg.center;
    if (cfg.view == OBStereo::ViewFrom)
      out << ", viewFrom = ";
    else
      out << ", viewTowards = ";

    if (cfg.from == OBStereo::ImplicitRef)
      out << "H";
    else
      out << cfg.from;

    out << ", refs = ";
    for (OBStereo::Refs::const_iterator i = cfg.refs.begin(); i != cfg.refs.end(); ++i)
      if (*i != OBStereo::ImplicitRef)
        out << *i << " ";
      else
        out << "H ";

    if (!cfg.specified)
      out << ", unspecified)";
    else {
      if (cfg.winding == OBStereo::Clockwise)
        out << ", clockwise)";
      else
        out << ", anti-clockwise)";
    }

    return out;
  }

} // namespace std

