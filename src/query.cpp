/**********************************************************************
  query.cpp - OBQuery, OBQueryAtom & OBQueryBond classes.

  Copyright (C) 2010 by Tim Vandermeersch

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

#include <openbabel/query.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>
#include <algorithm>

using namespace std;

namespace OpenBabel {

  OBQuery::~OBQuery()
  {
    std::for_each(m_atoms.begin(), m_atoms.end(), DeleteObject());
    std::for_each(m_bonds.begin(), m_bonds.end(), DeleteObject());
  }


  OBQuery* CompileMoleculeQuery(OBMol *mol, const OBBitVec &mask)
  {
    // set all atoms to 1 if the mask is empty
    OBBitVec mask2 = mask;
    if (!mask2.CountBits())
      for (unsigned int i = 0; i < mol->NumAtoms(); ++i)
        mask2.SetBitOn(i + 1);

    OBQuery *query = new OBQuery;
    unsigned int offset = 0;
    std::vector<unsigned int> indexes;
    FOR_ATOMS_OF_MOL (obatom, mol) {
      indexes.push_back(obatom->GetIndex() - offset);
      if (!mask2.BitIsSet(obatom->GetIndex() + 1)) {
        offset++;
        continue;
      }
      query->AddAtom(new OBQueryAtom(obatom->GetAtomicNum(), obatom->IsInRing(), obatom->IsAromatic()));
    }
    FOR_BONDS_OF_MOL (obbond, mol) {
      unsigned int beginIndex = obbond->GetBeginAtom()->GetIndex();
      unsigned int endIndex = obbond->GetEndAtom()->GetIndex();
      if (!mask2.BitIsSet(beginIndex + 1) || !mask2.BitIsSet(endIndex + 1))
        continue;

      query->AddBond(new OBQueryBond(query->GetAtoms()[indexes[beginIndex]], query->GetAtoms()[indexes[endIndex]],
            obbond->GetBondOrder(), obbond->IsAromatic()));
    }

    return query;
  }

  OBQuery* CompileSmilesQuery(const std::string &smiles, const OBBitVec &mask)
  {
    OBConversion conv;
    conv.SetInFormat("smi");
    OBMol mol;
    conv.ReadString(&mol, smiles);
    return CompileMoleculeQuery(&mol, mask);
  }

}

/// @file query.cpp
/// @brief OBQuery, OBQueryAtom & OBQueryBond classes.
