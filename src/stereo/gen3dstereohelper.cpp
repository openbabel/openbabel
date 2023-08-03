/**********************************************************************
  gen3dstereohelper.cpp - Helper class for 3D coordinate generation

  Copyright (C) 2020 by Tim Vandermeersch

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

#include "gen3dstereohelper.h"
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>

namespace OpenBabel {

  void OBGen3DStereoHelper::Setup(OBMol *mol)
  {
    m_unspecifiedTetrahedral.clear();
    m_unspecifiedCisTrans.clear();

    // Store canonical SMILES of original molecule
    OBConversion conv;
    conv.SetOutFormat("can");
    m_inputSmiles = conv.WriteString(mol, true);

    // Keep track of unspecified stereochemistry
    OBStereoFacade facade(mol);

    std::vector<OBTetrahedralStereo*> tetrahedral = facade.GetAllTetrahedralStereo();
    for (std::size_t i = 0; i < tetrahedral.size(); ++i) {
      OBTetrahedralStereo::Config cfg = tetrahedral[i]->GetConfig();
      if (!cfg.specified)
        m_unspecifiedTetrahedral.push_back(cfg.center);
    }

    std::vector<OBCisTransStereo*> cistrans = facade.GetAllCisTransStereo();
    for (std::size_t i = 0; i < cistrans.size(); ++i) {
      OBCisTransStereo::Config cfg = cistrans[i]->GetConfig();
      OBAtom *begin = mol->GetAtomById(cfg.begin);
      OBAtom *end = mol->GetAtomById(cfg.end);
      if (!begin || !end)
        continue;
      OBBond *bond = mol->GetBond(begin, end);
      if (!bond)
        continue;
      if (!cfg.specified)
        m_unspecifiedCisTrans.push_back(bond->GetId());
    }
  }

  bool OBGen3DStereoHelper::Check(OBMol *mol)
  {
    // Perceive stereo from 3D coords
    StereoFrom3D(mol, true); // true  = force

    // Make sure to respect previously unspecifed stereochemistry
    OBStereoFacade facade(mol);

    for (std::size_t i = 0; i < m_unspecifiedTetrahedral.size(); ++i) {
      OBTetrahedralStereo *ts = facade.GetTetrahedralStereo(m_unspecifiedTetrahedral[i]);
      if (!ts)
        continue;
      OBTetrahedralStereo::Config cfg = ts->GetConfig();
      cfg.specified = false;
      ts->SetConfig(cfg);
    }

    for (std::size_t i = 0; i < m_unspecifiedCisTrans.size(); ++i) {
      OBCisTransStereo *ct = facade.GetCisTransStereo(m_unspecifiedCisTrans[i]);
      if (!ct)
        continue;
      OBCisTransStereo::Config cfg = ct->GetConfig();
      cfg.specified = false;
      ct->SetConfig(cfg);
    }

    // Generate canonical SMILES with stereochemistry perceived from 3D coords.
    OBConversion conv;
    conv.SetOutFormat("can");
    std::string predictedSmiles = conv.WriteString(mol, true);

    return m_inputSmiles == predictedSmiles;
  }

} // namespace OpenBabel
