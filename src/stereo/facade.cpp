/**********************************************************************
  stereofacade.cpp - OBStereoFacade

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
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/squareplanar.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>

namespace OpenBabel {

  unsigned int OBStereoFacade::NumTetrahedralStereo()
  {
    EnsureInit();
    return static_cast<unsigned int> (m_tetrahedralMap.size());
  }

  unsigned int OBStereoFacade::NumCisTransStereo()
  {
    EnsureInit();
    return static_cast<unsigned int> (m_cistransMap.size());
  }

  unsigned int OBStereoFacade::NumSquarePlanarStereo()
  {
    EnsureInit();
    return static_cast<unsigned int> (m_squarePlanarMap.size());
  }

  bool OBStereoFacade::HasTetrahedralStereo(unsigned long atomId)
  {
    EnsureInit();
    if (m_tetrahedralMap.find(atomId) != m_tetrahedralMap.end())
      return true;
    return false;
  }

  bool OBStereoFacade::HasCisTransStereo(unsigned long bondId)
  {
    EnsureInit();
    if (m_cistransMap.find(bondId) != m_cistransMap.end())
      return true;
    return false;
  }

  bool OBStereoFacade::HasSquarePlanarStereo(unsigned long atomId)
  {
    EnsureInit();
    if (m_squarePlanarMap.find(atomId) != m_squarePlanarMap.end())
      return true;
    return false;
  }

  OBTetrahedralStereo* OBStereoFacade::GetTetrahedralStereo(unsigned long atomId)
  {
    if (!HasTetrahedralStereo(atomId))
      return 0;
    return m_tetrahedralMap[atomId];
  }

  OBCisTransStereo* OBStereoFacade::GetCisTransStereo(unsigned long bondId)
  {
    if (!HasCisTransStereo(bondId))
      return 0;
    return m_cistransMap[bondId];
  }

  OBSquarePlanarStereo* OBStereoFacade::GetSquarePlanarStereo(unsigned long atomId)
  {
    if (!HasSquarePlanarStereo(atomId))
      return 0;
    return m_squarePlanarMap[atomId];
  }

  void OBStereoFacade::InitMaps()
  {
    if (m_perceive && !m_mol->HasChiralityPerceived())
      PerceiveStereo(m_mol);

    std::vector<OBGenericData *> stereoData = m_mol->GetAllData(OBGenericDataType::StereoData);

    std::vector<OBGenericData*>::iterator data;
    for (data = stereoData.begin(); data != stereoData.end(); ++data) {
      OBStereo::Type type = ((OBStereoBase*)*data)->GetType();
      if (type == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo::Config config = ts->GetConfig();
        if (config.center == OBStereo::NoRef)
          continue;
        m_tetrahedralMap[config.center] = ts;
      } else
      if (type == OBStereo::SquarePlanar) {
        OBSquarePlanarStereo *sp = dynamic_cast<OBSquarePlanarStereo*>(*data);
        OBSquarePlanarStereo::Config config = sp->GetConfig();
        if (config.center == OBStereo::NoRef)
          continue;
        m_squarePlanarMap[config.center] = sp;
      } else
      if (type == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        OBCisTransStereo::Config config = ct->GetConfig();
        // find the bond id from begin & end atom ids
        unsigned long id = OBStereo::NoRef;
        OBAtom *a = m_mol->GetAtomById(config.begin);
        if (!a)
          continue;
        FOR_BONDS_OF_ATOM (bond, a) {
          unsigned long beginId = bond->GetBeginAtom()->GetId();
          unsigned long endId = bond->GetEndAtom()->GetId();
          if ((beginId == config.begin && endId == config.end) ||
              (beginId == config.end && endId == config.begin)) {
            id = bond->GetId();
            break;
          }
        }
        if (id == OBStereo::NoRef)
          continue;
        m_cistransMap[id] = ct;
      }
    }

    m_init = true;
  }

  template<>
  bool OBStereoFacade::HasStereo<OBStereo::Tetrahedral>(unsigned long id)
  {
    return HasTetrahedralStereo(id);
  }

  template<>
  bool OBStereoFacade::HasStereo<OBStereo::CisTrans>(unsigned long id)
  {
    return HasCisTransStereo(id);
  }

  template<>
  OBTetrahedralStereo* OBStereoFacade::GetStereo<OBTetrahedralStereo>(unsigned long id)
  {
    return GetTetrahedralStereo(id);
  }

  template<>
  OBCisTransStereo* OBStereoFacade::GetStereo<OBCisTransStereo>(unsigned long id)
  {
    return GetCisTransStereo(id);
  }




}

