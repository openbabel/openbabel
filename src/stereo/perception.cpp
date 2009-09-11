/**********************************************************************
  perception.cpp - Stereochemistry perception

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
#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>
#include <openbabel/graphsym.h>
#include <openbabel/oberror.h>

namespace OpenBabel {
 
  ////////////////////////////////////////////////////////////////////////////
  //
  //  General
  //
  ////////////////////////////////////////////////////////////////////////////

  void PerceiveStereo(OBMol *mol, bool force)
  {
    switch (mol->GetDimension()) {
      case 3:
        StereoFrom3D(mol, force);
        break;
      case 2:
        StereoFrom2D(mol, force);
        break;
      default:
        StereoFrom0D(mol);
        break;
    }
    
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::PerceiveStereo", obAuditMsg);
  }

  std::vector<OBBitVec> mergeRings(OBMol *mol)
  {
    std::vector<OBRing*> rings = mol->GetSSSR();

    std::vector<std::vector<OBRing*> > bridgedRings;
    for (int n = 0; rings.size(); ++n) {
      cout << "rings.size = " << rings.size() << endl;
      OBRing *ring = rings[0];
      std::vector<OBRing*> bridge;
      // add the first ring
      bridge.push_back(ring);

      // add rings that share an atom with the current first ring
      for (unsigned int j = 1; j < rings.size(); ++j) {
        for (unsigned int k = 0; k < ring->_path.size(); ++k)
          if (std::find(rings[j]->_path.begin(),
                        rings[j]->_path.end(),
                        ring->_path.at(k)) != rings[j]->_path.end())
            bridge.push_back(rings[j]);
      }

      bridgedRings.push_back(bridge);
      cout << "bridge.size = " << bridge.size() << endl;

      // erase new rings from original list
      for (unsigned int i = 0; i < bridge.size(); ++i) {
        std::vector<OBRing*>::iterator newEnd = std::remove(rings.begin(), rings.end(), bridge[i]);
        rings.erase(newEnd, rings.end());
      }

    }

    // store the merged rings in OBBitVecs
    std::vector<OBBitVec> result;
    for (unsigned int i = 0; i < bridgedRings.size(); ++i) {
      OBBitVec bits;
      for (unsigned int j = 0; j < bridgedRings.at(i).size(); ++j)
        for (unsigned int k = 0; k < bridgedRings.at(i).at(j)->_path.size(); ++k)
          bits.SetBitOn( bridgedRings.at(i).at(j)->_path[k] );

      result.push_back(bits);
    }
    cout << "# free ring systems = " << result.size() << endl;

    return result;
  }

  bool isInSameMergedRing(const std::vector<OBBitVec> &mergedRings, unsigned int idx1, unsigned int idx2)
  {
    std::vector<OBBitVec>::const_iterator bits;
    for (bits = mergedRings.begin(); bits != mergedRings.end(); ++bits)
      if ((*bits).BitIsSet( idx1 ) && (*bits).BitIsSet( idx2 ))
        return true;
    return false;  
  }

  std::vector<unsigned long> FindTetrahedralAtoms(OBMol *mol, const std::vector<unsigned int> &symClasses)
  {
    std::vector<unsigned long> centers;

    // do quick test to see if there are any possible chiral centers
    bool mayHaveChiralCenter = false;
    OBAtom *atom, *nbr;
    std::vector<OBAtom*>::iterator i;
    for (atom = mol->BeginAtom(i); atom; atom = mol->NextAtom(i))
      if (atom->GetHyb() == 3 && atom->GetHvyValence() >= 3) {
        mayHaveChiralCenter = true;
        break;
      }

    if (!mayHaveChiralCenter)
      return centers;
    if (symClasses.size() != mol->NumAtoms())
      return centers;

    std::vector<unsigned int> plist; // para-stereocenters canditates

    bool ischiral;
    for (atom = mol->BeginAtom(i); atom; atom = mol->NextAtom(i)) {
      if (atom->IsNitrogen() || atom->IsPhosphorus() || atom->IsSulfur())
        continue;
      if (atom->GetHyb() == 3 && atom->GetHvyValence() >= 3) {
        // list containing neighbor symmetry classes
        std::vector<unsigned int> tlist; 
        ischiral = true;

        // check neighbors to see if this atom is stereogenic
        std::vector<OBBond*>::iterator j;
        for (nbr = atom->BeginNbrAtom(j); nbr; nbr = atom->NextNbrAtom(j)) {
          // check if we already have a neighbor with this symmetry class
          std::vector<unsigned int>::iterator k;
          for (k = tlist.begin(); k != tlist.end(); ++k)
            if (symClasses[nbr->GetIndex()] == *k) {
              ischiral = false;
              // if so, might still be a para-stereocenter
              plist.push_back(atom->GetIdx());
            }

          if (ischiral)
            // keep track of all neighbors, so we can detect duplicates
            tlist.push_back(symClasses[nbr->GetIndex()]);
          else
            break;
        }

        if (ischiral) {
          // true-stereocenter found
          centers.push_back(atom->GetId());
        }
      }
    }

    // find para-stereocenters
    if (plist.size()) {
      std::vector<OBBitVec> mergedRings = mergeRings(mol);

      for (unsigned int i = 0; i < plist.size(); ++i) {
        for (unsigned int j = 0; j < plist.size(); ++j) {
          if (i == j)
            continue;
          if (isInSameMergedRing(mergedRings, plist.at(i), plist.at(j))) {
            centers.push_back(mol->GetAtom(plist.at(i))->GetId());
            centers.push_back(mol->GetAtom(plist.at(j))->GetId());
          }
        }
      }
    }

    // make sure we don't return duplicates
    std::sort(centers.begin(), centers.end());
    std::vector<unsigned long>::iterator newEnd = std::unique(centers.begin(), centers.end());
    centers.erase(newEnd, centers.end());

    return centers;
  }
  
  std::vector<unsigned long> FindCisTransBonds(OBMol *mol, const std::vector<unsigned int> &symClasses)
  {
    std::vector<unsigned long> bonds;
          
    //do quick test to see if there are any possible chiral centers
    bool mayHaveCisTransBond = false;
    std::vector<OBBond*>::iterator i;
    for (OBBond *bond = mol->BeginBond(i); bond; bond = mol->NextBond(i))
      if (bond->GetBO() == 2 && !bond->IsInRing()) {
        mayHaveCisTransBond = true;
        break;
      }

    if (!mayHaveCisTransBond)
      return bonds;
    if (symClasses.size() != mol->NumAtoms())
      return bonds;

    bool isCisTrans;
    for (OBBond *bond = mol->BeginBond(i); bond; bond = mol->NextBond(i)) {
      if (bond->IsInRing())
        continue;

      if (bond->GetBO() == 2) {
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();
        if (!begin || !end) 
          continue;

        // Needs to have at least one explicit single bond at either end
        if (!begin->HasSingleBond() || !end->HasSingleBond())
          continue;
          
        isCisTrans = true;
        std::vector<OBBond*>::iterator j;
         
        if (begin->GetValence() == 2) {
          // Begin atom has two explicit neighbors. One is the end atom. The other should 
          // be a heavy atom - this is what we test here.
          // (There is a third, implicit, neighbor which is either a hydrogen
          // or a lone pair.)
          if (begin->ExplicitHydrogenCount() == 1)
            isCisTrans = false;
        } else if (begin->GetValence() == 3) {
          std::vector<unsigned int> tlist;
          
          for (OBAtom *nbr = begin->BeginNbrAtom(j); nbr; nbr = begin->NextNbrAtom(j)) {
            // skip end atom
            if (nbr->GetId() == end->GetId())
              continue;
            // do we already have an atom with this symmetry class?
            if (tlist.size()) {
              // compare second with first
              if (symClasses[nbr->GetIndex()] == tlist.at(0))
                isCisTrans = false;
              break;
            }
              
            // save first summetry class
            tlist.push_back(symClasses[nbr->GetIndex()]);
          }
        } else {
          // Valence is not 2 or 3, for example SR3=NR
          isCisTrans = false;
        }

        if (!isCisTrans) 
          continue;

        if (end->GetValence() == 2) {
          // see comment above for begin atom
          if (end->ExplicitHydrogenCount() == 1)
            isCisTrans = false;
        } else if (end->GetValence() == 3) { 
          std::vector<unsigned int> tlist;
          
          for (OBAtom *nbr = end->BeginNbrAtom(j); nbr; nbr = end->NextNbrAtom(j)) {
            // skip end atom
            if (nbr->GetId() == begin->GetId())
              continue;
            // do we already have an atom with this symmetry class?
            if (tlist.size()) {
              // compare second with first
              if (symClasses[nbr->GetIndex()] == tlist.at(0))
                isCisTrans = false;
              break;
            }
                
            // save first summetry class
            tlist.push_back(symClasses[nbr->GetIndex()]);
          }
        } else {
          // Valence is not 2 or 3, for example SR3=NR
          isCisTrans = false;
        }

        if (isCisTrans)
          bonds.push_back(bond->GetId());
      }
    }

    return bonds;
  }
  
  void FindParaStereocenters(OBMol *mol, const std::vector<unsigned int> &symClasses, 
      const std::vector<unsigned long> &tetrahedralAtomIds, const std::vector<unsigned long> &cistransBondItds,
      std::vector<unsigned long> &paraAtomIds, std::vector<unsigned long> &paraAtomIds)
  {
  
  }

  /**
   * Perform symmetry analysis.
   *
   * @return vector containing symmetry classes index by OBAtom::GetIndex().
   */
  std::vector<unsigned int> FindSymmetry(OBMol *mol)
  {
    OBGraphSym symmetry(mol);
    std::vector<unsigned int> symClasses;
    symmetry.GetSymmetry(symClasses);
    return symClasses;
  }

  ////////////////////////////////////////////////////////////////////////////
  //
  //  From0D
  //
  ////////////////////////////////////////////////////////////////////////////

  void StereoFrom0D(OBMol *mol, std::map <OBBond*, OBStereo::BondDirection> *updown)
  {
    if (mol->HasChiralityPerceived())
      return;
     
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom0D", obAuditMsg);

    std::vector<unsigned int> symClasses = FindSymmetry(mol);
    TetrahedralFrom0D(mol, symClasses);
    CisTransFrom0D(mol, symClasses, updown);
    mol->SetChiralityPerceived();
  }

  std::vector<OBTetrahedralStereo*> TetrahedralFrom0D(OBMol *mol, 
      const std::vector<unsigned int> symClasses, bool addToMol)
  {
    std::vector<OBTetrahedralStereo*> configs;
    // make sure the number of atoms matches the symClasses' size
    if (symClasses.size() != mol->NumAtoms())
      return configs;
    
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::TetrahedralFrom0D", obAuditMsg);

    // find all tetrahedral centers
    std::vector<unsigned long> centers = FindTetrahedralAtoms(mol, symClasses);

    // Delete any existing stereo objects that are not a member of 'centers'
    // and make a map of the remaining ones
    std::map<unsigned long, OBTetrahedralStereo*> existingMap;
    std::vector<OBGenericData*>::iterator data;
    std::vector<OBGenericData*> stereoData = mol->GetAllData(OBGenericDataType::StereoData);
    for (data = stereoData.begin(); data != stereoData.end(); ++data) {
      if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        unsigned long center = ts->GetConfig().center;
        if (std::find(centers.begin(), centers.end(), center) == centers.end()) {
          // According to OpenBabel, this is not a tetrahedral stereo
          obErrorLog.ThrowError(__FUNCTION__, "Removed spurious TetrahedralStereo object", obAuditMsg);
          mol->DeleteData(ts);
        }
        else {
          existingMap[center] = ts;
          configs.push_back(ts);
        }
      }
    }

    std::vector<unsigned long>::iterator i;
    for (i = centers.begin(); i != centers.end(); ++i) {
      // if there already exists a OBTetrahedralStereo object for this 
      // center, continue
      if (existingMap.find(*i) != existingMap.end())
        continue;

      OBAtom *center = mol->GetAtomById(*i);
 
      OBTetrahedralStereo::Config config;
      config.specified = false;
      config.center = *i;
      FOR_NBORS_OF_ATOM(nbr, center) {
        if (config.from == OBStereo::NoRef)
          config.from = nbr->GetId();
        else
          config.refs.push_back(nbr->GetId());
      }

      if ((config.refs.size() == 2))
        config.refs.push_back(OBStereo::ImplicitRef); // need to add largest number on end to work

      OBTetrahedralStereo *th = new OBTetrahedralStereo(mol);
      th->SetConfig(config);
      
      configs.push_back(th);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(th);
    }

    return configs;    
  }

  std::vector<OBCisTransStereo*> CisTransFrom0D(OBMol *mol, 
      const std::vector<unsigned int> &symClasses, 
      std::map<OBBond*, enum OBStereo::BondDirection> *updown, bool addToMol)
  {
    std::vector<OBCisTransStereo*> configs;
    if (symClasses.size() != mol->NumAtoms())
      return configs;

    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::CisTransFrom0D", obAuditMsg);
 
    // find all cis/trans bonds
    std::vector<unsigned long> bonds = FindCisTransBonds(mol, symClasses);

    // Add any CisTransStereo objects indicated by the 'updown' map
    if (updown && updown->size() > 0)
        CisTransFromUpDown(mol, bonds, updown);

    // Delete any existing stereo objects that are not a member of 'bonds'
    // and make a map of the remaining ones
    std::map<unsigned long, OBCisTransStereo*> existingMap;
    std::vector<OBGenericData*>::iterator data;
    std::vector<OBGenericData*> stereoData = mol->GetAllData(OBGenericDataType::StereoData);
    for (data = stereoData.begin(); data != stereoData.end(); ++data) {
      if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        OBCisTransStereo::Config config = ct->GetConfig();
        // find the bond id from begin & end atom ids
        unsigned long id = OBStereo::NoRef;
        OBAtom *a = mol->GetAtomById(config.begin);
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

        if (std::find(bonds.begin(), bonds.end(), id) == bonds.end()) {
          // According to OpenBabel, this is not a cis trans stereo
          obErrorLog.ThrowError(__FUNCTION__, "Removed spurious CisTransStereo object", obAuditMsg);
          mol->DeleteData(ct);
        }
        else {
          existingMap[id] = ct;
          configs.push_back(ct);
        }
      }
    }

    std::vector<unsigned long>::iterator i;
    for (i = bonds.begin(); i != bonds.end(); ++i) {
      // If there already exists a OBCisTransStereo object for this 
      // bond, leave it alone
      if (existingMap.find(*i) != existingMap.end())
        continue;

      OBBond *bond = mol->GetBondById(*i);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();

      OBCisTransStereo::Config config;
      config.specified = false;
      // begin
      config.begin = begin->GetId();
      FOR_NBORS_OF_ATOM (nbr, begin) {
        if (nbr->GetId() == end->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
      }
      if (config.refs.size() == 1) {
        config.refs.push_back(OBStereo::ImplicitRef);
      }
      // end
      config.end = end->GetId();
      FOR_NBORS_OF_ATOM (nbr, end) {
        if (nbr->GetId() == begin->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
      }
      if (config.refs.size() == 3) {
        config.refs.push_back(OBStereo::ImplicitRef);
      }

      OBCisTransStereo *ct = new OBCisTransStereo(mol);
      ct->SetConfig(config);
      
      configs.push_back(ct);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(ct);
    }

    return configs;
  }



  ////////////////////////////////////////////////////////////////////////////
  //
  //  From3D
  //
  ////////////////////////////////////////////////////////////////////////////

  void StereoFrom3D(OBMol *mol, bool force)
  {
    if (mol->HasChiralityPerceived() && !force)
      return;
     
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom3D", obAuditMsg);

    mol->DeleteData(OBGenericDataType::StereoData);
    std::vector<unsigned int> symClasses = FindSymmetry(mol);
    TetrahedralFrom3D(mol, symClasses);
    CisTransFrom3D(mol, symClasses);
    mol->SetChiralityPerceived();
  }

  //! Calculate the "sign of a volume" given by a set of 4 coordinates
  double VolumeSign(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d)
  {
    vector3 A, B, C;
    A = b - a;
    B = c - a;
    C = d - a;
    matrix3x3 m(A, B, C);
    return m.determinant();
  }

  std::vector<OBTetrahedralStereo*> TetrahedralFrom3D(OBMol *mol, 
      const std::vector<unsigned int> symClasses, bool addToMol)
  {
    std::vector<OBTetrahedralStereo*> configs;
    // make sure the number of atoms matches the symClasses' size
    if (symClasses.size() != mol->NumAtoms())
      return configs;
 
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::TetrahedralFrom3D", obAuditMsg);

    // find all tetrahedral centers
    std::vector<unsigned long> centers = FindTetrahedralAtoms(mol, symClasses);
      
    std::vector<unsigned long>::iterator i;
    for (i = centers.begin(); i != centers.end(); ++i) {
      OBAtom *center = mol->GetAtomById(*i);
 
      // make sure we have at least 3 heavy atom neighbors
      // timvdm 28 Jun 2009: This is already checked in FindTetrahedralAtoms
      if (center->GetHvyValence() < 3) {
        std::stringstream errorMsg;
        errorMsg << "Cannot calculate a signed volume for an atom with a heavy atom valence of " 
                 << center->GetHvyValence() << std::endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        continue;
      }
        
      OBTetrahedralStereo::Config config;
      config.center = *i;
      FOR_NBORS_OF_ATOM(nbr, center) {
        if (config.from == OBStereo::NoRef)
          config.from = nbr->GetId();
        else
          config.refs.push_back(nbr->GetId());
      }

      bool use_central_atom = false;
           
      // Create a vector with the coordinates of the neighbor atoms
      std::vector<vector3> nbrCoords;
      OBAtom *from = mol->GetAtomById(config.from);
      nbrCoords.push_back(from->GetVector());
      for (OBStereo::RefIter id = config.refs.begin(); id != config.refs.end(); ++id) {
        OBAtom *nbr = mol->GetAtomById(*id);
        nbrCoords.push_back(nbr->GetVector());
      }
    
        // Checks for a neighbour having 0 co-ords (added hydrogen etc)
        /* FIXME: needed? if the molecule has 3D coords, additional 
         * hydrogens will get coords using OBAtom::GetNewBondVector
        for (std::vector<vector3>::iterator coord = nbrCoords.begin(); coord != nbrCoords.end(); ++coord) { 
          // are the coordinates zero to 6 or more significant figures
          if (coord->IsApprox(VZero, 1.0e-6)) {
            if (!use_central_atom) {
              use_central_atom = true;
            } else {
              obErrorLog.ThrowError(__FUNCTION__, 
                  "More than 2 neighbours have 0 co-ords when attempting 3D chiral calculation", obInfo);
            }
          }
        }
        */

      // If we have three heavy atoms we can use the chiral center atom itself for the fourth
      // will always give same sign (for tetrahedron), magnitude will be smaller.
      if ((config.refs.size() == 2) || use_central_atom) {
        nbrCoords.push_back(center->GetVector());
        config.refs.push_back(OBStereo::ImplicitRef); // need to add largest number on end to work
      }

      double sign = VolumeSign(nbrCoords[0], nbrCoords[1], nbrCoords[2], nbrCoords[3]);
      if (sign < 0.0)
        config.winding = OBStereo::AntiClockwise;

      OBTetrahedralStereo *th = new OBTetrahedralStereo(mol);
      th->SetConfig(config);
      
      configs.push_back(th);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(th);
    }

    return configs;    
  }

  std::vector<OBCisTransStereo*> CisTransFrom3D(OBMol *mol, 
      const std::vector<unsigned int> &symClasses, bool addToMol)
  {
    std::vector<OBCisTransStereo*> configs;
    if (symClasses.size() != mol->NumAtoms())
      return configs;
    
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::CisTransFrom3D", obAuditMsg);
 
    // find all cis/trans bonds
    std::vector<unsigned long> bonds = FindCisTransBonds(mol, symClasses);
    
    std::vector<unsigned long>::iterator i;
    for (i = bonds.begin(); i != bonds.end(); ++i) {
      OBBond *bond = mol->GetBondById(*i);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();

      // Create a vector with the coordinates of the neighbor atoms
      std::vector<vector3> bondVecs;
      OBCisTransStereo::Config config;
      // begin
      config.begin = begin->GetId();
      FOR_NBORS_OF_ATOM (nbr, begin) {
        if (nbr->GetId() == end->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
        bondVecs.push_back(nbr->GetVector() - begin->GetVector());
      }
      if (config.refs.size() == 1) {
        config.refs.push_back(OBStereo::ImplicitRef);
        vector3 pos;
        mol->GetAtomById(config.refs.at(0))->GetNewBondVector(pos, 1.0);
        bondVecs.push_back(pos - begin->GetVector());
      }
      // end
      config.end = end->GetId();
      FOR_NBORS_OF_ATOM (nbr, end) {
        if (nbr->GetId() == begin->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
        bondVecs.push_back(nbr->GetVector() - end->GetVector());
      }
      if (config.refs.size() == 3) {
        config.refs.push_back(OBStereo::ImplicitRef);
        vector3 pos;
        mol->GetAtomById(config.refs.at(2))->GetNewBondVector(pos, 1.0);
        bondVecs.push_back(pos - end->GetVector());
      }

      // 0      3       0   3
      //  \    /         \ /      angle 0-*-3 & 1-*-2: 60 degrees (cis)
      //   C==C    -->    *       angle 0-*-1 & 2-*-3: 120 degrees (same bond atom)
      //  /    \         / \      angle 0-*-2 & 1-*-3: 180 degrees (trans)
      // 1      2       1   2
      double angle = vectorAngle(bondVecs[0], bondVecs[2]);
      if (IsNear(angle, 60.0, 10.0))
        config.shape = OBStereo::ShapeZ;
      if (IsNear(angle, 180.0, 10.0))
        config.shape = OBStereo::ShapeU;

      OBCisTransStereo *ct = new OBCisTransStereo(mol);
      ct->SetConfig(config);
      
      configs.push_back(ct);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(ct);
    }

    return configs;
  }

  ////////////////////////////////////////////////////////////////////////////
  //
  //  From2D
  //
  //  Reference:
  //  [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the 
  //  Unambiguous Identification of the Stereochemical Characteristics of 
  //  Compounds During Their Registration in Databases. Molecules 2000, 6,
  //  915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
  ////////////////////////////////////////////////////////////////////////////

  void StereoFrom2D(OBMol *mol, bool tetfrom0D, bool force)
  {
    if (mol->HasChiralityPerceived() && !force)
      return;
      
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom2D", obAuditMsg);

    std::vector<unsigned int> symClasses = FindSymmetry(mol);
    if (!tetfrom0D) {
      mol->DeleteData(OBGenericDataType::StereoData);
      TetrahedralFrom2D(mol, symClasses);
    }
    else { // In MDL format, the "s" option means that tetstereo is set from 0D info
      TetrahedralFrom0D(mol, symClasses);
      std::vector<OBGenericData*>::iterator data;
      std::vector<OBGenericData*> stereoData = mol->GetAllData(OBGenericDataType::StereoData);
      for (data = stereoData.begin(); data != stereoData.end(); ++data)
        if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::CisTrans)
          mol->DeleteData(*data);
    }
    CisTransFrom2D(mol, symClasses);
    mol->SetChiralityPerceived();
  }
 
  //! Calculate the "sign of a triangle" given by a set of 3 2D coordinates
  double TriangleSign(const vector3 &a, const vector3 &b, const vector3 &c)
  {
    // equation 6 from [1]
    return (a.x() - c.x()) * (b.y() - c.y()) - (a.y() - c.y()) * (b.x() - c.x());
  }

  std::vector<OBTetrahedralStereo*> TetrahedralFrom2D(OBMol *mol, 
      const std::vector<unsigned int> symClasses, bool addToMol)
  {
    std::vector<OBTetrahedralStereo*> configs;
    // make sure the number of atoms matches the symClasses' size
    if (symClasses.size() != mol->NumAtoms())
      return configs;
    
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::TetrahedralFrom2D", obAuditMsg);
 
    // find all tetrahedral centers
    std::vector<unsigned long> centers = FindTetrahedralAtoms(mol, symClasses);
      
    std::vector<unsigned long>::iterator i;
    for (i = centers.begin(); i != centers.end(); ++i) {
      OBAtom *center = mol->GetAtomById(*i);
 
      // make sure we have at least 3 heavy atom neighbors
      if (center->GetHvyValence() < 3) {
        std::stringstream errorMsg;
        errorMsg << "Cannot calculate a signed volume for an atom with a heavy atom valence of " 
                 << center->GetHvyValence() << std::endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        continue;
      }
        
      OBTetrahedralStereo::Config config;
      config.center = *i;
        
      // find the hash, wedge and 2 plane atoms
      std::vector<OBAtom*> planeAtoms;
      std::vector<OBAtom*> wedgeAtoms;
      std::vector<OBAtom*> hashAtoms;
      FOR_BONDS_OF_ATOM(bond, center) {
        OBAtom *nbr = bond->GetNbrAtom(center);
        // hash bonds
        if (bond->IsHash()) {
          if (bond->GetBeginAtom()->GetId() == center->GetId()) {
            // this is a 'real' hash bond going from center to nbr
            hashAtoms.push_back(nbr);
          } else {
            // this is an 'inverted' hash bond going from nbr to center
            wedgeAtoms.push_back(nbr);
          }
        } else if (bond->IsWedge()) {
          // wedge bonds
          if (bond->GetBeginAtom()->GetId() == center->GetId()) {
            // this is a 'real' wedge bond going from center to nbr
            wedgeAtoms.push_back(nbr);
          } else {
            // this is an 'inverted' hash bond going from nbr to center
            hashAtoms.push_back(nbr);
          }
        } else if (bond->IsWedgeOrHash()) {
          config.specified = false;
          break;
        } else { 
          // plane bonds
          planeAtoms.push_back(nbr);
        }
      }

      bool success = true;
      
      using namespace std;
      if (!config.specified) {
        // unspecified
        FOR_NBORS_OF_ATOM (nbr, center)
          if (config.from == OBStereo::NoRef)
            config.from = nbr->GetId();
          else
            config.refs.push_back(nbr->GetId());
        while (config.refs.size() < 3)
          config.refs.push_back(OBStereo::ImplicitRef);
      } else
      if (planeAtoms.size() == 2) {
        if (hashAtoms.size() == 1 && wedgeAtoms.size() == 1) {
          // plane1 + plane2, hash, wedge
          config.from = wedgeAtoms[0]->GetId();
          config.refs.resize(3);
          config.refs[0] = planeAtoms[0]->GetId();
          config.refs[1] = planeAtoms[1]->GetId();
          config.refs[2] = hashAtoms[0]->GetId();
          double sign = TriangleSign(planeAtoms[0]->GetVector(), 
              planeAtoms[1]->GetVector(), hashAtoms[0]->GetVector());
          if (sign > 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else if (hashAtoms.size() == 1 && wedgeAtoms.size() == 0) {
          // plane1 + plane2, hash
          config.from = OBStereo::ImplicitRef;
          config.refs.resize(3);
          config.refs[0] = planeAtoms[0]->GetId();
          config.refs[1] = planeAtoms[1]->GetId();
          config.refs[2] = hashAtoms[0]->GetId();
          double sign = TriangleSign(planeAtoms[0]->GetVector(), 
              planeAtoms[1]->GetVector(), hashAtoms[0]->GetVector());
          if (sign > 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else if (hashAtoms.size() == 0 && wedgeAtoms.size() == 1) {
          // plane1 + plane2, wedge
          config.towards = OBStereo::ImplicitRef;
          config.view = OBStereo::ViewTowards;
          config.refs.resize(3);
          config.refs[0] = planeAtoms[0]->GetId();
          config.refs[1] = planeAtoms[1]->GetId();
          config.refs[2] = wedgeAtoms[0]->GetId();
          double sign = TriangleSign(planeAtoms[0]->GetVector(), 
              planeAtoms[1]->GetVector(), wedgeAtoms[0]->GetVector());
          if (sign > 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else {
          success = false;
        }
      } else if (planeAtoms.size() == 3) {
        if (hashAtoms.size() == 1 && wedgeAtoms.size() == 0) {
          // plane1 + plane2 + plane3, hash
          config.towards = hashAtoms[0]->GetId();
          config.view = OBStereo::ViewTowards;
          config.refs.resize(3);
          config.refs[0] = planeAtoms[0]->GetId();
          config.refs[1] = planeAtoms[1]->GetId();
          config.refs[2] = planeAtoms[2]->GetId();
          double sign = TriangleSign(planeAtoms[0]->GetVector(), 
              planeAtoms[1]->GetVector(), planeAtoms[2]->GetVector());
          if (sign < 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else if (hashAtoms.size() == 0 && wedgeAtoms.size() == 1) {
          // plane1 + plane2, wedge
          config.from = wedgeAtoms[0]->GetId();
          config.view = OBStereo::ViewFrom;
          config.refs.resize(3);
          config.refs[0] = planeAtoms[0]->GetId();
          config.refs[1] = planeAtoms[1]->GetId();
          config.refs[2] = planeAtoms[2]->GetId();
          double sign = TriangleSign(planeAtoms[0]->GetVector(), 
              planeAtoms[1]->GetVector(), planeAtoms[2]->GetVector());
          if (sign < 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else {
          success = false;
        }
      
      } else {
        success = false;
      }

      if (!success) {
//         std::stringstream errorMsg;
//         errorMsg << "Symmetry analysis found atom with id " << center->GetId() 
//             << " to be a tetrahedral atom but the wedge/hash bonds can't be interpreted." << std::endl
//             << " # in-plane bonds = " << planeAtoms.size() << std::endl
//             << " # wedge bonds = " << wedgeAtoms.size() << std::endl
//             << " # hash bonds = " << hashAtoms.size() << std::endl
//             << std::endl;
//         obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        continue;
      }
 

      OBTetrahedralStereo *th = new OBTetrahedralStereo(mol);
      th->SetConfig(config);

      configs.push_back(th);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(th);
    }
   
    return configs;
  }

  std::vector<OBCisTransStereo*> CisTransFrom2D(OBMol *mol, 
      const std::vector<unsigned int> &symClasses, bool addToMol)
  {
    std::vector<OBCisTransStereo*> configs;
    // make sure the number of atoms matches the symClasses' size
    if (symClasses.size() != mol->NumAtoms())
      return configs;

    // find all cis/trans bonds
    std::vector<unsigned long> bonds = FindCisTransBonds(mol, symClasses);
    
    std::vector<unsigned long>::iterator i;
    for (i = bonds.begin(); i != bonds.end(); ++i) {
      OBBond *bond = mol->GetBondById(*i);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();

      // Create a vector with the coordinates of the neighbor atoms
      std::vector<vector3> bondVecs;
      OBCisTransStereo::Config config;
      // begin
      config.begin = begin->GetId();
      FOR_NBORS_OF_ATOM (nbr, begin) {
        if (nbr->GetId() == end->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
        bondVecs.push_back(nbr->GetVector());
      }
      if (config.refs.size() == 1) {
        config.refs.push_back(OBStereo::ImplicitRef);
        vector3 pos;
        begin->GetNewBondVector(pos, 1.0);
        bondVecs.push_back(pos);
      }
      // end
      config.end = end->GetId();
      FOR_NBORS_OF_ATOM (nbr, end) {
        if (nbr->GetId() == begin->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
        bondVecs.push_back(nbr->GetVector());
      }
      if (config.refs.size() == 3) {
        config.refs.push_back(OBStereo::ImplicitRef);
        vector3 pos;
        end->GetNewBondVector(pos, 1.0);
        bondVecs.push_back(pos);
      }

      // 0      3       
      //  \    /        2 triangles: 0-1-b & 2-3-a
      //   a==b    -->  same sign: U
      //  /    \        opposite sign: Z
      // 1      2       
      /*
      double sign1 = TriangleSign(begin->GetVector(), end->GetVector(), bondVecs[0]);
      double sign2 = TriangleSign(begin->GetVector(), end->GetVector(), bondVecs[2]);
      */
      double sign1 = TriangleSign(bondVecs[0], bondVecs[1], end->GetVector());
      double sign2 = TriangleSign(bondVecs[2], bondVecs[3], begin->GetVector());
      double sign = sign1 * sign2;

      if (sign < 0.0) // opposite sign
        config.shape = OBStereo::ShapeZ;

      OBCisTransStereo *ct = new OBCisTransStereo(mol);
      ct->SetConfig(config);

      configs.push_back(ct);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(ct);
    }

    return configs;
  }
  void CisTransFromUpDown(OBMol *mol, const std::vector<unsigned long> &ctbonds,
    std::map<OBBond*, OBStereo::BondDirection> *updown)
  {
    // Create a vector of CisTransStereo objects for the molecule

    // Loop across the known cistrans bonds
    vector<unsigned long>::const_iterator bondId_it;
    for (bondId_it = ctbonds.begin(); bondId_it != ctbonds.end(); bondId_it++) {
      OBBond* dbl_bond = mol->GetBondById(*bondId_it);
      
      OBAtom *a1 = dbl_bond->GetBeginAtom();
      OBAtom *a2 = dbl_bond->GetEndAtom();

      // Get the bonds of neighbors of atom1 and atom2
      OBBond *a1_b1 = NULL, *a1_b2 = NULL, *a2_b1 = NULL, *a2_b2 = NULL;
      OBStereo::BondDirection a1_stereo, a2_stereo;

      FOR_BONDS_OF_ATOM(bi, a1) {
        OBBond *b = &(*bi);
        if (b == dbl_bond) continue;  // skip the double bond we're working on
        if (a1_b1 == NULL && updown->find(b) != updown->end())
        {
          a1_b1 = b;    // remember a stereo bond of Atom1
          a1_stereo = (*updown)[b];
        }
        else
          a1_b2 = b;    // remember a 2nd bond of Atom1
      }

      FOR_BONDS_OF_ATOM(bi, a2) {
        OBBond *b = &(*bi);
        if (b == dbl_bond) continue;
        if (a2_b1 == NULL && updown->find(b) != updown->end())
        {
          a2_b1 = b;    // remember a stereo bond of Atom2
          a2_stereo = (*updown)[b];
        }
        else
          a2_b2 = b;    // remember a 2nd bond of Atom2
      }
      
      if (a1_b1 == NULL || a2_b1 == NULL) continue; // No cis/trans
      
      // a1_b2 and/or a2_b2 will be NULL if there are bonds to implicit hydrogens
      unsigned int second = (a1_b2 == NULL) ? OBStereo::ImplicitRef : a1_b2->GetNbrAtom(a1)->GetId();
      unsigned int fourth = (a2_b2 == NULL) ? OBStereo::ImplicitRef : a2_b2->GetNbrAtom(a2)->GetId();

      // If a1_stereo==a2_stereo, this means cis for a1_b1 and a2_b1.
      OBCisTransStereo *ct = new OBCisTransStereo(mol);
      OBCisTransStereo::Config cfg;
      cfg.begin = a1->GetId();
      cfg.end = a2->GetId();

      

      if (a1_stereo == a2_stereo)
        cfg.refs = OBStereo::MakeRefs(a1_b1->GetNbrAtom(a1)->GetId(), second,
                                      fourth, a2_b1->GetNbrAtom(a2)->GetId());
      else
        cfg.refs = OBStereo::MakeRefs(a1_b1->GetNbrAtom(a1)->GetId(), second,
                                      a2_b1->GetNbrAtom(a2)->GetId(), fourth);
      if (a1_stereo == OBStereo::UnknownDir || a2_stereo == OBStereo::UnknownDir)
        cfg.specified = false;

      ct->SetConfig(cfg);
      // add the data to the atom
      mol->SetData(ct);
    }
  } 
}

