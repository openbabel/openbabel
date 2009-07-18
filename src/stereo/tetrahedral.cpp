/**********************************************************************
  tetrahedral.h - OBTetrahedralStereo

  Copyright (C) 2009 by Tim Vandermeersch
 
  This file is part of the Open Babel project.
  For more information, see <http://openbabel.sourceforge.net/>

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

  bool OBTetrahedralStereo::Config::operator==(const Config &other) const
  {
    if (center != other.center)
      return false;
    if ((refs.size() != 3) || (other.refs.size() != 3))
      return false;
    
    Config otherConfig = OBTetraNonPlanarStereo::ToConfig(other, from, winding, view);


    if (!OBStereo::ContainsSameRefs(refs, otherConfig.refs))
      return false;

    // normalize the other Config struct
    int Ni1 = OBStereo::NumInversions(refs);
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
    if (m_cfg.center == OBStereo::NoId)
      return false;
    if (m_cfg.from == OBStereo::NoId)
      return false;
    if (m_cfg.refs.size() != 3)
      return false;
    return true;
  }

  void OBTetrahedralStereo::SetConfig(const Config &config)
  {
    if (config.center == OBStereo::NoId) {
      obErrorLog.ThrowError(__FUNCTION__, 
          "OBTetrahedralStereo::SetConfig : center atom id is invalid.", obError);
      m_cfg = Config();
      return;
    }
    if (config.from == OBStereo::NoId) {
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

    return OBTetraNonPlanarStereo::ToConfig(m_cfg, m_cfg.from, winding, view);
  }
  
  OBTetrahedralStereo::Config OBTetrahedralStereo::GetConfig(unsigned long from_or_towards, 
        OBStereo::Winding winding, OBStereo::View view) const
  {
    if (!IsValid())
      return Config();

    return OBTetraNonPlanarStereo::ToConfig(m_cfg, from_or_towards, winding, view);
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

  ostream& operator<<(ostream &out, const OpenBabel::OBTetrahedralStereo &ts)
  {
    OpenBabel::OBTetrahedralStereo::Config cfg = ts.GetConfig();
    out << "OBTetrahedralStereo(center = " << cfg.center;
    if (cfg.view == OpenBabel::OBStereo::ViewFrom)
      out << ", viewFrom = " << cfg.from;
    else
      out << ", viewTowards = " << cfg.towards;
 
    out << ", refs = ";
    for (OpenBabel::OBStereo::Refs::iterator i = cfg.refs.begin(); i != cfg.refs.end(); ++i)
      out << *i << " ";

    if (cfg.winding == OpenBabel::OBStereo::Clockwise)
      out << ", clockwise)";
    else
      out << ", anti-clockwise)";

    return out;
  }

  ostream& operator<<(ostream &out, const OpenBabel::OBTetrahedralStereo::Config &cfg)
  {
    out << "OBTetrahedralStereo::Config(center = " << cfg.center;
    if (cfg.view == OpenBabel::OBStereo::ViewFrom)
      out << ", viewFrom = " << cfg.from;
    else
      out << ", viewTowards = " << cfg.towards;

    out << ", refs = ";
    for (OpenBabel::OBStereo::Refs::const_iterator i = cfg.refs.begin(); i != cfg.refs.end(); ++i)
      out << *i << " ";

    if (cfg.winding == OpenBabel::OBStereo::Clockwise)
      out << ", clockwise)";
    else
      out << ", anti-clockwise)";

    return out;
  }

} // namespace std

