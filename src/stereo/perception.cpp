/**********************************************************************
  perception.cpp - Stereochemistry perception

  Copyright (C) 2009-2010 by Tim Vandermeersch

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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <openbabel/obutil.h>
#include <openbabel/obiter.h>
#include <openbabel/graphsym.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/canon.h>
#include <openbabel/oberror.h>
#include <openbabel/elements.h>
#include <cassert>

#include "stereoutil.h"

#include <cmath>
#include <limits>
#include <set>
#include <iterator>
#include <functional>

#define DEBUG 0
#define DEBUG_INVERSIONS 0
#define IMPLICIT_CIS_RING_SIZE 8
#define DELTA_ANGLE_FOR_OVERLAPPING_BONDS 4.0 // In degrees

using namespace std;

// debug function
template<typename T>
void print_vector(const std::string &label, const std::vector<T> &v)
{
  std::cout << label << ": ";
  for (std::size_t i = 0; i < v.size(); ++i)
    std::cout << v[i] << " ";
  std::cout << endl;
}

namespace OpenBabel {

  OBAtom* findAtomWithSymmetryClass(OBAtom *atom, unsigned int symClass, const std::vector<unsigned int> &symClasses);
  bool containsAtLeast_1true_2para(OBAtom *ligandAtom, OBAtom *atom, const OBStereoUnitSet &units);
  bool containsAtLeast_2true_2paraAssemblies(OBAtom *ligandAtom, OBAtom *atom, const OBStereoUnitSet &units, const std::vector<OBBitVec> &mergedRings);

  //////////////////////////////////////////////////////////////////////////////
  //
  //  General
  //
  //////////////////////////////////////////////////////////////////////////////

  void PerceiveStereo(OBMol *mol, bool force)
  {
    switch (mol->GetDimension()) {
      case 3:
        StereoFrom3D(mol, force);
        break;
      case 2:
        StereoFrom2D(mol, 0, force);
        break;
      default:
        StereoFrom0D(mol);
        break;
    }
    if (obErrorLog.GetOutputLevel() >= obAuditMsg)
      obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::PerceiveStereo", obAuditMsg);
  }

  /**
   * Perform a quick check for tetrahedral stereo centers. Used by
   * FindStereogenicUnits to return quickly if there are no stereogenic units.
   */
  bool mayHaveTetrahedralCenter(OBMol *mol)
  {
    std::vector<OBAtom*>::iterator ia;
    for (OBAtom *atom = mol->BeginAtom(ia); atom; atom = mol->NextAtom(ia))
      if (atom->GetHyb() == 3 && atom->GetHvyDegree() >= 3) {
        return true;
      }
    return false;
  }

  /**
   * Perform a quick check for stereogenic bonds. Used by FindStereogenicUnits
   * to return quickly if there are no stereogenic units.
   */
  bool mayHaveCisTransBond(OBMol *mol)
  {
    std::vector<OBBond*>::iterator ib;
    for (OBBond *bond = mol->BeginBond(ib); bond; bond = mol->NextBond(ib))
      if (bond->GetBondOrder() == 2) {
        return true;
      }
    return false;
  }

  /**
   * Check if the specified atom is a potential stereogenic atom.
   *
   * Criteria:
   * - sp3 hybridization (or P and sp3d hybridization)
   * - not connected to more than 4 atoms
   * - at least 3 "heavy" neighbors
   *
   * Nitrogen (neutral) is treated as a special case since the barrier of inversion is
   * low in many cases making the atom non-stereogenic. Only bridge-head
   * nitrogen atoms (i.e. nitrogen has 3 neighbors in rings) will be
   * considered stereogenic.
   */
  bool isPotentialTetrahedral(OBAtom *atom)
  {
    // consider only potential steroecenters
    if ((atom->GetHyb() != 3 && !(atom->GetHyb() == 5 && atom->GetAtomicNum() == OBElements::Phosphorus))
        || atom->GetTotalDegree() > 4 || atom->GetHvyDegree() < 3 || atom->GetHvyDegree() > 4)
      return false;
    // skip non-chiral N
    if (atom->GetAtomicNum() == OBElements::Nitrogen && atom->GetFormalCharge()==0) {
      int nbrRingAtomCount = 0;
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (nbr->IsInRing())
          nbrRingAtomCount++;
      }
      if (nbrRingAtomCount < 3)
        return false;
    }
    if (atom->GetAtomicNum() == OBElements::Carbon) {
      if (atom->GetFormalCharge())
        return false;
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (nbr->GetAtomicNum() == 26 && nbr->GetExplicitDegree() > 7)
          return false;
      }
    }

    return true;
  }

  /**
   * Check if the specified bond is a potential stereogenic bond.
   *
   * Criteria:
   * - must be a double bond
   * - must not be in a ring
   * - both begin and end atom should have at least one single bond
   */
  bool isPotentialCisTrans(OBBond *bond)
  {
    if (bond->GetBondOrder() != 2)
      return false;
    if (bond->IsInRing())
      return false;
    if (!bond->GetBeginAtom()->HasSingleBond() || !bond->GetEndAtom()->HasSingleBond())
      return false;
    if (bond->GetBeginAtom()->GetHvyDegree() == 1 || bond->GetEndAtom()->GetHvyDegree() == 1)
      return false;
    if (bond->GetBeginAtom()->GetHvyDegree() > 3 || bond->GetEndAtom()->GetHvyDegree() > 3)
      return false;
    return true;
  }










  ////////////////////////////////////////////////////////////////////////////////


  /**
   * Check if the specified stereogenic unit is in a fragment.
   */
  bool isUnitInFragment(OBMol *mol, const OBStereoUnit &unit, const OBBitVec &fragment)
  {
    if (unit.type == OBStereo::Tetrahedral) {
      if (fragment.BitIsSet(unit.id))
        return true;
    } else if(unit.type == OBStereo::CisTrans) {
      OBBond *bond = mol->GetBondById(unit.id);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();
      if (fragment.BitIsSet(begin->GetId()) || fragment.BitIsSet(end->GetId()))
        return true;
    }
    return false;
  }









  //////////////////////////////////////////////////////////////////////////////////





  /**
   * Check if the specified atom is a tetrahedral center (i.e. there is a Tetrahedral
   * OBStereoUnit in units with the same id)
   */
  bool isTetrahedral(OBAtom *atom, const OBStereoUnitSet &units)
  {
    for (std::size_t i = 0; i < units.size(); ++i) {
      const OBStereoUnit &unit = units[i];
      if (unit.type != OBStereo::Tetrahedral)
        continue;
      if (unit.id == atom->GetId())
        return true;
    }
    return false;
  }

  /**
   * Check if the specified bond is a double bond stereocenter (i.e. there is a CisTrans
   * OBStereoUnit in units with the same id)
   */
  bool isCisTrans(OBBond *bond, const OBStereoUnitSet &units)
  {
    for (std::size_t i = 0; i < units.size(); ++i) {
      const OBStereoUnit &unit = units[i];
      if (unit.type != OBStereo::CisTrans)
        continue;
      if (unit.id == bond->GetId())
        return true;
    }
    return false;
  }


  /**
   * Classify the tetrahedral atom using the NeighborSymmetryClasses types.
   */
  int classifyTetrahedralNbrSymClasses(const std::vector<unsigned int> &symClasses, OBAtom *atom)
  {
    std::vector<unsigned int> nbrClasses, nbrClassesCopy, uniqueClasses;
    FOR_NBORS_OF_ATOM (nbr, atom)
      nbrClasses.push_back(symClasses.at(nbr->GetIndex()));
    // add an implicit ref if there are only 3 explicit
    if (nbrClasses.size() == 3)
      nbrClasses.push_back(OBStereo::ImplicitRef);

    // use some STL to work out the number of unique classes
    nbrClassesCopy = nbrClasses; // keep copy for count below
    std::sort(nbrClasses.begin(), nbrClasses.end());
    std::vector<unsigned int>::iterator endLoc = std::unique(nbrClasses.begin(), nbrClasses.end());
    std::copy(nbrClasses.begin(), endLoc, std::back_inserter(uniqueClasses));

    switch (uniqueClasses.size()) {
      case 4:
        return T1234; // e.g. 1 2 3 4
      case 3:
        return T1123; // e.g. 1 1 2 3
      case 2:
        // differentiate between T1122 and T1112
        if (std::count(nbrClassesCopy.begin(), nbrClassesCopy.end(), uniqueClasses.at(0)) == 2)
          return T1122; // e.g. 1 1 2 2
        else
          return T1112; // e.g. 1 1 1 2
      case 1:
	  default:
        return T1111; // e.g. 1 1 1 1
    }
  }

  /**
   * Classify the cis/trans bond using the NeighborSymmetryClasses types.
   */
  int classifyCisTransNbrSymClasses(const std::vector<unsigned int> &symClasses, OBBond *doubleBond, OBAtom *atom)
  {
    std::vector<unsigned int> nbrClasses, uniqueClasses;
    FOR_NBORS_OF_ATOM (nbr, atom) {
      if (nbr->GetIdx() != doubleBond->GetNbrAtom(atom)->GetIdx())
        nbrClasses.push_back(symClasses.at(nbr->GetIndex()));
    }

    if (nbrClasses.size() == 1)
      nbrClasses.push_back(OBStereo::ImplicitRef);

    if (nbrClasses.at(0) == nbrClasses.at(1))
      return C11; // e.g. 1 1
    else
      return C12; // e.g. 1 2
  }

  /**
   * Merge the rings in a molecule and return the result as OBBitVec objects.
   * Rings are merged if they share at least one atom (e.g. bridged, spiro,
   * adjacent, ...).
   */
  std::vector<OBBitVec> mergeRings(OBMol *mol, const std::vector<unsigned int> &symClasses)
  {
    std::vector<OBRing*> rings = mol->GetSSSR();

    std::vector<OBBitVec> result;
    for (std::size_t i = 0; i < rings.size(); ++i) {
      // check if ring shares atom with previously found ring
      bool found = false;
      for (std::size_t j = 0; j < result.size(); ++j) {
        std::vector<unsigned int> shared;
        // foreach ring atom
        for (std::size_t k = 0; k < rings[i]->_path.size(); ++k) {
          // check if the ring atom is in the current result bitvec
          if (result[j].BitIsSet(rings[i]->_path[k])) {
            shared.push_back(rings[i]->_path[k]);
          }
        }

        if (shared.size() > 1) {
          found = true;
        } else if (shared.size() == 1) {
          int classification = classifyTetrahedralNbrSymClasses(symClasses, mol->GetAtom(shared[0]));
          if (classification == T1122 || classification == T1111)
            found = true;
        }

        if (found) {
          // add bits for the atoms in the ring
          for (std::size_t l = 0; l < rings[i]->_path.size(); ++l)
            result[j].SetBitOn(rings[i]->_path[l]);
          break;
        }
      }

      // add the ring as a new bitvec if it shares no atom with a previous ring
      if (!found) {
        OBBitVec r;
        for (std::size_t l = 0; l < rings[i]->_path.size(); ++l)
          r.SetBitOn(rings[i]->_path[l]);
        result.push_back(r);
      }
    }

    return result;
  }

  /*
  bool isInSameMergedRing(const std::vector<OBBitVec> &mergedRings, unsigned int idx1, unsigned int idx2)
  {
    std::vector<OBBitVec>::const_iterator bits;
    for (bits = mergedRings.begin(); bits != mergedRings.end(); ++bits)
      if ((*bits).BitIsSet( idx1 ) && (*bits).BitIsSet( idx2 ))
        return true;
    return false;
  }
  */

  /**
   * Helper function for getFragment below.
   */
  void addNbrs(OBBitVec &fragment, OBAtom *atom, OBAtom *skip)
  {
    FOR_NBORS_OF_ATOM (nbr, atom) {
      // don't pass through skip
      if (nbr->GetId() == skip->GetId())
        continue;
      // skip visited atoms
      if (fragment.BitIsSet(nbr->GetId()))
        continue;
      // add the neighbor atom to the fragment
      fragment.SetBitOn(nbr->GetId());
      // recurse...
      addNbrs(fragment, &*nbr, skip);
    }
  }

  /**
   * Create an OBBitVec objects with bets set for the fragment consisting of all
   * atoms for which there is a path to atom without going through skip. These
   * fragment bitvecs are indexed by unique id (i.e. OBAtom::GetId()).
   */
  OBBitVec getFragment(OBAtom *atom, OBAtom *skip)
  {
    OBBitVec fragment;
    fragment.SetBitOn(atom->GetId());
    // start the recursion
    addNbrs(fragment, atom, skip);
    return fragment;
  }


  struct StereoRing
  {
    struct ParaAtom
    {
      typedef OBAtom CenterType;

      ParaAtom(unsigned long _id, unsigned int idx) : id(_id), inIdx(idx) {}
      OBAtom* GetCenter(OBMol *mol) const { return mol->GetAtomById(id); }
      bool isInRing(const StereoRing &ring) const
      {
        for (std::size_t i = 0; i < ring.paraAtoms.size(); ++i)
          if (ring.paraAtoms[i].inIdx == inIdx)
            return true;
        return false;
      }

      unsigned long id;
      union {
        unsigned int inIdx, outIdx;
      };
      std::vector<OBAtom*> insideNbrs, outsideNbrs;
    };
    struct ParaBond
    {
      typedef OBBond CenterType;
      ParaBond(unsigned long _id, unsigned int _inIdx, unsigned int _outIdx) : id(_id), inIdx(_inIdx), outIdx(_outIdx) {}
      OBBond* GetCenter(OBMol *mol) const { return mol->GetBondById(id); }
      bool isInRing(const StereoRing &ring) const
      {
        for (std::size_t i = 0; i < ring.paraBonds.size(); ++i)
          if (ring.paraBonds[i].inIdx == inIdx)
            return true;
        return false;
      }

      unsigned long id;
      unsigned int inIdx, outIdx;
      std::vector<OBAtom*> insideNbrs, outsideNbrs;
    };

    StereoRing() : trueCount(0) {}

    std::vector<ParaAtom> paraAtoms;
    std::vector<ParaBond> paraBonds;
    unsigned int trueCount;
  };

  template<typename Type>
  bool checkLigands(const Type &currentPara, const OBStereoUnitSet &units)
  {
    if (currentPara.outsideNbrs.size() == 1) {
      //cout << "OK: " << __LINE__ << endl;
      return true;
    }
    OBMol *mol = currentPara.insideNbrs[0]->GetParent();
    assert(mol->GetAtom(currentPara.outIdx));
    OBBitVec ligand = getFragment(currentPara.outsideNbrs[0], mol->GetAtom(currentPara.outIdx));
    for (OBStereoUnitSet::const_iterator u2 = units.begin(); u2 != units.end(); ++u2) {
      if (isUnitInFragment(mol, *u2, ligand)) {
        //cout << "OK: " << __LINE__ << endl;
        return true;
      }
    }
    //cout << "NOT OK: " << __LINE__ << endl;
    return false;
  }


  template<typename Type>
  bool ApplyRule1(const Type &currentPara, const std::vector<unsigned int> &symmetry_classes,
      const std::vector<StereoRing> &rings, std::vector<bool> &visitedRings, const OBStereoUnitSet &units,
      std::vector<unsigned int> stereoAtoms)
  {
    bool foundRing = false;
    unsigned int idx = currentPara.inIdx;

    /*
    for (std::size_t i = 0; i < visitedRings.size(); ++i)
      if (visitedRings[i])
        cout << "  ";
    cout << "ApplyRule1(" << currentPara.inIdx << ", " << currentPara.outIdx << ", outside = " << currentPara.outsideNbrs.size() << ")" << endl;
    */

    for (std::size_t i = 0; i < rings.size(); ++i) {
      // skip visited rings
      if (visitedRings[i])
        continue;

      // Check if currentPara is in this ring
      if (!currentPara.isInRing(rings[i]))
        continue;

      //
      // A new ring containing currentPara is found
      //
      foundRing = true;

      // if there are one or more true stereo centers, currentPara is a stereo center
      if (rings[i].trueCount) {
        //cout << "OK: " << __LINE__ << endl;
        return true;
      }

      // check if there is at least one other potential atom
      for (std::size_t j = 0; j < rings[i].paraAtoms.size(); ++j) {
        const StereoRing::ParaAtom &paraAtom = rings[i].paraAtoms[j];
        // skip idx
        if (paraAtom.inIdx == idx)
          continue;
        // there is another atom already identified as stereo atom
        if (std::find(stereoAtoms.begin(), stereoAtoms.end(), paraAtom.inIdx) != stereoAtoms.end()) {
          //cout << "OK: " << __LINE__ << endl;
          return true;
        }

        if (paraAtom.outsideNbrs.size() == 1) {
          // only 1 ring substituent, the other is implicit H -> topologically different
          //cout << "OK: " << __LINE__ << endl;
          return true;
        } else {
          if (paraAtom.outsideNbrs.size() != 2)
            return false;
          // two ring substituents, need to check for topological difference
          if (symmetry_classes[paraAtom.outsideNbrs[0]->GetIndex()] != symmetry_classes[paraAtom.outsideNbrs[1]->GetIndex()]) {
            // they are different
            //cout << "OK: " << __LINE__ << endl;
            return true;
          } else {
            // they are the same and they might also be in a ring -> apply rule 1 recursive
            visitedRings[i] = true;
            if (ApplyRule1(paraAtom, symmetry_classes, rings, visitedRings, units, stereoAtoms)) {
              //cout << "OK: " << __LINE__ << endl;
              return true;
            }
          }
        }
      }
      // check if there is at least one other potential bond
      for (std::size_t j = 0; j < rings[i].paraBonds.size(); ++j) {
        const StereoRing::ParaBond &paraBond = rings[i].paraBonds[j];
        // skip idx
        if (paraBond.inIdx == idx)
          continue;
        // there is another atom already identified as stereo atom
        if (std::find(stereoAtoms.begin(), stereoAtoms.end(), paraBond.inIdx) != stereoAtoms.end()) {
          //cout << "OK: " << __LINE__ << endl;
          return true;
        }

        if (paraBond.outsideNbrs.size() == 1) {
          // only 1 ring substituent, the other is implicit H -> topologically different
          //cout << "OK: " << __LINE__ << endl;
          return true;
        } else {
          if (paraBond.outsideNbrs.size() != 2)
            continue;
          // two ring substituents, need to check for topological difference
          if (symmetry_classes[paraBond.outsideNbrs[0]->GetIndex()] != symmetry_classes[paraBond.outsideNbrs[1]->GetIndex()]) {
            // they are different
            //cout << "OK: " << __LINE__ << endl;
            return true;
          } else {
            // they are the same and they might also be in a ring -> apply rule 1 recursive
            visitedRings[i] = true;
            if (ApplyRule1(paraBond, symmetry_classes, rings, visitedRings, units, stereoAtoms)) {
              //cout << "OK: " << __LINE__ << endl;
              return true;
            }
          }
        }
      }

    }

    // if a non-visited ring was found and true was not returned -> it does not
    // contain any stereocenters other than idx
    if (foundRing) {
      //cout << "NOT OK: " << __LINE__ << endl;
      return false;
    }

    //cout << "NOT OK: " << __LINE__ << endl;
    return false;
  }

  void StartRule1(const std::vector<unsigned int> &symmetry_classes, const std::vector<StereoRing> &rings,
      OBStereoUnitSet &units, std::vector<unsigned int> &stereoAtoms)
  {
    for (std::size_t i = 0; i < rings.size(); ++i) {
      //cout << "Checking ring: " << i << endl;

      // tetrahedral atoms
      for (std::size_t j = 0; j < rings[i].paraAtoms.size(); ++j) {
        const StereoRing::ParaAtom &paraAtom = rings[i].paraAtoms[j];
        // skip the atom if it is already in stereoAtoms
        if (std::find(stereoAtoms.begin(), stereoAtoms.end(), paraAtom.inIdx) != stereoAtoms.end())
          continue;

        std::vector<bool> visitedRings(rings.size(), false);
        //visitedRings[i] = true;
        if (ApplyRule1(paraAtom, symmetry_classes, rings, visitedRings, units, stereoAtoms)) {
          bool isStereoUnit = false;
          if (paraAtom.outsideNbrs.size() == 1)
            isStereoUnit = true;
          if (paraAtom.outsideNbrs.size() == 2) {
            if (symmetry_classes[paraAtom.outsideNbrs[0]->GetIndex()] == symmetry_classes[paraAtom.outsideNbrs[1]->GetIndex()]) {
              // check for spiro atom
              bool isSpiro = false;
              for (std::size_t k = 0; k < rings[i].paraAtoms.size(); ++k) {
                const StereoRing::ParaAtom &paraAtom2 = rings[i].paraAtoms[k];
                if (paraAtom.inIdx == paraAtom2.outIdx && paraAtom.insideNbrs == paraAtom2.outsideNbrs) {
                  isSpiro = true;
                  if (ApplyRule1(paraAtom2, symmetry_classes, rings, visitedRings, units, stereoAtoms))
                    isStereoUnit = true;
                }
              }
              if (!isSpiro)
                isStereoUnit = checkLigands(paraAtom, units);
              //cout << "isStereoUnit = " << isStereoUnit << endl;
            } else {
              isStereoUnit = true;
            }
          }

          if (isStereoUnit) {
            stereoAtoms.push_back(paraAtom.inIdx);
            OBAtom *atom = paraAtom.insideNbrs[0]->GetParent()->GetAtomById(paraAtom.id);
            units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId(), true));
          }
        }

      }

      // cistrans bonds
      for (std::size_t j = 0; j < rings[i].paraBonds.size(); ++j) {
        const StereoRing::ParaBond &paraBond = rings[i].paraBonds[j];
        // skip the atom if it is already in stereoAtoms
        if (std::find(stereoAtoms.begin(), stereoAtoms.end(), paraBond.inIdx) != stereoAtoms.end())
          continue;

        std::vector<bool> visitedRings(rings.size(), false);
        //visitedRings[i] = true;
        if (ApplyRule1(paraBond, symmetry_classes, rings, visitedRings, units, stereoAtoms)) {
          bool isStereoUnit = false;
          if (paraBond.outsideNbrs.size() == 1)
            isStereoUnit = true;
          if (paraBond.outsideNbrs.size() == 2) {
            if (symmetry_classes[paraBond.outsideNbrs[0]->GetIndex()] == symmetry_classes[paraBond.outsideNbrs[1]->GetIndex()]) {
              // check for spiro bond
              bool isSpiro = false;
              for (std::size_t k = 0; k < rings[i].paraBonds.size(); ++k) {
                const StereoRing::ParaBond &paraBond2 = rings[i].paraBonds[k];
                if (paraBond.inIdx == paraBond2.outIdx && paraBond.insideNbrs == paraBond2.outsideNbrs) {
                  isSpiro = true;
                  if (ApplyRule1(paraBond2, symmetry_classes, rings, visitedRings, units, stereoAtoms))
                    isStereoUnit = true;
                }
              }
              if (!isSpiro)
                isStereoUnit = checkLigands(paraBond, units);
              //cout << "isStereoUnit = " << isStereoUnit << endl;
            } else {
              isStereoUnit = true;
            }
          }

          if (isStereoUnit) {
            stereoAtoms.push_back(paraBond.inIdx);
            stereoAtoms.push_back(paraBond.outIdx);
            OBBond *bond = paraBond.insideNbrs[0]->GetParent()->GetBondById(paraBond.id);
            units.push_back(OBStereoUnit(OBStereo::CisTrans, bond->GetId(), true));
          }
        }

      }


    }

  }

  /**
   * Find the stereogenic units in a molecule using a set of rules.
   *
   * This is a public function: see header for details.
   */
  OBStereoUnitSet FindStereogenicUnits(OBMol *mol, const std::vector<unsigned int> &symClasses)
  {
    OBStereoUnitSet units;

    // do quick test to see if there are any possible stereogenic units
    if (!mayHaveTetrahedralCenter(mol) && !mayHaveCisTransBond(mol))
      return units;

    // make sure we have symmetry classes for all atoms
    if (symClasses.size() != mol->NumAtoms())
      return units;

    // para-stereocenters candidates
    std::vector<unsigned int> stereoAtoms; // Tetrahedral = idx, CisTrans = begin & end idx
    std::vector<unsigned int> paraAtoms;
    std::vector<unsigned int> paraBonds;

    /**
     * true Tetrahedral stereocenters:
     * - have four different symmetry classes for the ligands to the central atom
     */
    bool ischiral;
    std::vector<OBAtom*>::iterator ia;
    for (OBAtom *atom = mol->BeginAtom(ia); atom; atom = mol->NextAtom(ia)) {
      if (!isPotentialTetrahedral(atom))
        continue;

      // list containing neighbor symmetry classes
      std::vector<unsigned int> tlist;
      ischiral = true;

      // check neighbors to see if this atom is stereogenic
      std::vector<OBBond*>::iterator j;
      for (OBAtom *nbr = atom->BeginNbrAtom(j); nbr; nbr = atom->NextNbrAtom(j)) {
        // check if we already have a neighbor with this symmetry class
        std::vector<unsigned int>::iterator k;
        for (k = tlist.begin(); k != tlist.end(); ++k)
          if (symClasses[nbr->GetIndex()] == *k) {
            ischiral = false;
            // if so, might still be a para-stereocenter
            paraAtoms.push_back(atom->GetIdx());
          }

        if (ischiral)
          // keep track of all neighbors, so we can detect duplicates
          tlist.push_back(symClasses[nbr->GetIndex()]);
        else
          break;
      }

      if (ischiral) {
        // true-stereocenter found
        stereoAtoms.push_back(atom->GetIdx());
        units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId()));
      }
    }

    /**
     * true CisTrans stereocenters:
     * - each terminal has two different symmetry classes for it's ligands
     */
    bool isCisTransBond;
    std::vector<OBBond*>::iterator ib;
    for (OBBond *bond = mol->BeginBond(ib); bond; bond = mol->NextBond(ib)) {
      if (bond->IsInRing() && bond->IsAromatic())
        continue; // Exclude C=C in phenyl rings for example

      if (bond->GetBondOrder() == 2) {
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();
        if (!begin || !end)
          continue;

        if (begin->GetTotalDegree() > 3 || end->GetTotalDegree() > 3)
          continue; // e.g. C=Ru where the Ru has four substituents

        // Needs to have at least one explicit single bond at either end
        // FIXME: timvdm: what about C=C=C=C
        if (!begin->HasSingleBond() || !end->HasSingleBond())
          continue;

        isCisTransBond = true;
        std::vector<OBBond*>::iterator j;

        if (begin->GetExplicitDegree() == 2) {
          // Begin atom has two explicit neighbors. One is the end atom. The other should
          // be a heavy atom - this is what we test here.
          // (There is a third, implicit, neighbor which is either a hydrogen
          // or a lone pair.)
          if (begin->ExplicitHydrogenCount() == 1)
            isCisTransBond = false;
        } else if (begin->GetExplicitDegree() == 3) {
          std::vector<unsigned int> tlist;

          for (OBAtom *nbr = begin->BeginNbrAtom(j); nbr; nbr = begin->NextNbrAtom(j)) {
            // skip end atom
            if (nbr->GetId() == end->GetId())
              continue;
            // do we already have an atom with this symmetry class?
            if (tlist.size()) {
              // compare second with first
              if (symClasses[nbr->GetIndex()] == tlist.at(0)) {
                isCisTransBond = false;
                // if same, might still be a para-stereocenter
                paraBonds.push_back(bond->GetIdx());
              }
              break;
            }

            // save first symmetry class
            tlist.push_back(symClasses[nbr->GetIndex()]);
          }
        } else {
          // Valence is not 2 or 3, for example SR3=NR
          isCisTransBond = false;
        }

        if (!isCisTransBond)
          continue;

        if (end->GetExplicitDegree() == 2) {
          // see comment above for begin atom
          if (end->ExplicitHydrogenCount() == 1)
            isCisTransBond = false;
        } else if (end->GetExplicitDegree() == 3) {
          std::vector<unsigned int> tlist;

          for (OBAtom *nbr = end->BeginNbrAtom(j); nbr; nbr = end->NextNbrAtom(j)) {
            // skip end atom
            if (nbr->GetId() == begin->GetId())
              continue;
            // do we already have an atom with this symmetry class?
            if (tlist.size()) {
              // compare second with first
              if (symClasses[nbr->GetIndex()] == tlist.at(0)) {
                // if same, might still be a para-stereocenter
                paraBonds.push_back(bond->GetIdx());
                isCisTransBond = false;
              }
              break;
            }

            // save first symmetry class
            tlist.push_back(symClasses[nbr->GetIndex()]);
          }
        } else {
          // Valence is not 2 or 3, for example SR3=NR
          isCisTransBond = false;
        }

        if (isCisTransBond)
          // true-stereocenter found
          units.push_back(OBStereoUnit(OBStereo::CisTrans, bond->GetId()));
      }
    }

    /**
     * Apply rule 1 from the Razinger paper recusively:
     *
     * All rings are merged "mergedRings". A merged ring is simply a fragment consisting
     * of all atoms of a ring system (bridged, spiro, adjacent, ...). If two rings in the
     * SSSR set share an atom, they are merged.
     *
     * Each merged must at least have two para-stereocenters (or 1 true + 1 para) in order
     * for the para-stereocenter to be valid. This is repeated until no new stereocenters
     * are identified.
     *
     * rule 1a for double bonds:
     * - bond atom in ring has two identical symmetry classes for it's neighbor atoms (-> para)
     * - other bond atom:
     *   - has two different symmetry classes for it's neighbours -> new stereocenter
     *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
     *
     * rule 1b for tetracoord atoms:
     * - at least two neighbour symmetry classes are the same (-> para)
     * - other pair:
     *   - has two different symmetry classes for it's neighbours -> new stereocenter
     *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
     *
     * NOTE: there must always be at least 2 new stereocenters (or one existing + 1 newly found) in order for them to be valid
     */
    std::vector<OBRing*> lssr = mol->GetLSSR();
    std::vector<StereoRing> rings;

    //cout << "=====================================================" << endl;
    for (std::size_t i = 0; i < lssr.size(); ++i) {
      rings.push_back(StereoRing());
      StereoRing &ring = rings.back();


      for (std::size_t j = 0; j < stereoAtoms.size(); ++j)
        if (lssr[i]->_pathset.BitIsSet(stereoAtoms[j]))
          ring.trueCount++;

      //cout << "StereoRing: trueCount = " << ring.trueCount << endl;
      for (std::size_t j = 0; j < paraAtoms.size(); ++j) {
        if (lssr[i]->_pathset.BitIsSet(paraAtoms[j])) {
          OBAtom *atom = mol->GetAtom(paraAtoms[j]);
          ring.paraAtoms.push_back(StereoRing::ParaAtom(atom->GetId(), paraAtoms[j]));

          FOR_NBORS_OF_ATOM (nbr, mol->GetAtom(paraAtoms[j])) {
            if (lssr[i]->_pathset.BitIsSet(nbr->GetIdx()))
              ring.paraAtoms.back().insideNbrs.push_back(&*nbr);
            else
              ring.paraAtoms.back().outsideNbrs.push_back(&*nbr);
          }

          //cout << "  ParaAtom(idx = " << ring.paraAtoms.back().inIdx << ", outside = " << ring.paraAtoms.back().outsideNbrs.size() << ")" << endl;
          if (ring.paraAtoms.back().insideNbrs.size() != 2)
            ring.paraAtoms.pop_back();
        }
      }

      for (std::size_t j = 0; j < paraBonds.size(); ++j) {
        OBBond *bond = mol->GetBond(paraBonds[j]);
        unsigned int beginIdx = bond->GetBeginAtomIdx();
        unsigned int endIdx = bond->GetEndAtomIdx();

        if (lssr[i]->_pathset.BitIsSet(beginIdx)) {
          ring.paraBonds.push_back(StereoRing::ParaBond(bond->GetId(), beginIdx, endIdx));

          FOR_NBORS_OF_ATOM (nbr, bond->GetBeginAtom()) {
            if (nbr->GetIdx() == endIdx)
              continue;
            ring.paraBonds.back().insideNbrs.push_back(&*nbr);
          }
          FOR_NBORS_OF_ATOM (nbr, bond->GetEndAtom()) {
            if (nbr->GetIdx() == beginIdx)
              continue;
            ring.paraBonds.back().outsideNbrs.push_back(&*nbr);
          }

          //cout << "  ParaBond(inIdx = " << beginIdx << ", outIdx = " << endIdx << ", outside = " << ring.paraBonds.back().outsideNbrs.size() << ")" << endl;
          if (ring.paraBonds.back().insideNbrs.size() != 2)
            ring.paraBonds.pop_back();
        }

        if (lssr[i]->_pathset.BitIsSet(endIdx)) {
          ring.paraBonds.push_back(StereoRing::ParaBond(bond->GetId(), endIdx, beginIdx));

          FOR_NBORS_OF_ATOM (nbr, bond->GetEndAtom()) {
            if (nbr->GetIdx() == beginIdx)
              continue;
            ring.paraBonds.back().insideNbrs.push_back(&*nbr);
          }
          FOR_NBORS_OF_ATOM (nbr, bond->GetBeginAtom()) {
            if (nbr->GetIdx() == endIdx)
              continue;
            ring.paraBonds.back().outsideNbrs.push_back(&*nbr);
          }

          //cout << "  ParaBond(inIdx = " << endIdx << ", outIdx = " << beginIdx << ", outside = " << ring.paraBonds.back().outsideNbrs.size() << ")" << endl;
          if (ring.paraBonds.back().insideNbrs.size() != 2)
            ring.paraBonds.pop_back();
        }

      }

      if (ring.paraAtoms.size() + ring.paraBonds.size() == 1) {
        ring.paraAtoms.clear();
        ring.paraBonds.clear();
      }

    }
    //cout << "=====================================================" << endl;

    unsigned int numStereoUnits;
    do {
      numStereoUnits = units.size();
      StartRule1(symClasses, rings, units, stereoAtoms);
    } while (units.size() > numStereoUnits);


    std::vector<OBBitVec> mergedRings = mergeRings(mol, symClasses);
    /**
     * Apply rule 2a for tetracoordinate carbon:
     * - 1 or 2 pair identical ligands
     * - each pair contains at least 1 true-stereocenter or 2 para-stereocenters
     *
     * Apply rule 2b for tetracoordinate carbon:
     * - 3 or 4 identical ligands with at least
     *   - 2 true-stereocenters
     *   - 2 separate assemblies of para-stereocenters
     */
    for (std::vector<unsigned int>::iterator idx = paraAtoms.begin(); idx != paraAtoms.end(); ++idx) {
      OBAtom *atom = mol->GetAtom(*idx);
      // make sure we didn't add this atom already from rule 1
      bool alreadyAdded = false;
      for (OBStereoUnitSet::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
        if ((*u2).type == OBStereo::Tetrahedral)
          if (atom->GetId() == (*u2).id) {
            alreadyAdded = true;
          }
      }
      if (alreadyAdded)
        continue;


      int classification = classifyTetrahedralNbrSymClasses(symClasses, atom);
      switch (classification) {
        case T1123:
          // rule 2a with 1 pair
          {
            unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses);
            OBAtom *ligandAtom = findAtomWithSymmetryClass(atom, duplicatedSymClass, symClasses);
            if (containsAtLeast_1true_2para(ligandAtom, atom, units)) {
              units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId(), true));
            }
          }
          break;
        case T1122:
          // rule 2a with 2 pairs
          {
            unsigned int duplicatedSymClass1, duplicatedSymClass2;
            findDuplicatedSymmetryClasses(atom, symClasses, duplicatedSymClass1, duplicatedSymClass2);
            OBAtom *ligandAtom1 = findAtomWithSymmetryClass(atom, duplicatedSymClass1, symClasses);
            OBAtom *ligandAtom2 = findAtomWithSymmetryClass(atom, duplicatedSymClass2, symClasses);
            if (containsAtLeast_1true_2para(ligandAtom1, atom, units) &&
                containsAtLeast_1true_2para(ligandAtom2, atom, units))
              units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId(), true));
          }
          break;
        case T1112:
          // rule 2b with 3 identical
          {
            unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses);
            OBAtom *ligandAtom = findAtomWithSymmetryClass(atom, duplicatedSymClass, symClasses);
            if (containsAtLeast_2true_2paraAssemblies(ligandAtom, atom, units, mergedRings))
              units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId(), true));
          }
          break;
        case T1111:
          // rule 2b with 4 identical
          {
            unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses);
            OBAtom *ligandAtom = findAtomWithSymmetryClass(atom, duplicatedSymClass, symClasses);
            if (containsAtLeast_2true_2paraAssemblies(ligandAtom, atom, units, mergedRings)) {
              units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId(), true));
            }
          }
          break;

      }

    }

    /**
     * Apply rule 3 for double bonds.
     * - 1 or 2 pair identical ligands (on begin and end atom)
     * - each pair contains at least 1 true-stereocenter or 2 para-stereocenters (from rule1)
     */
    for (std::vector<unsigned int>::iterator idx = paraBonds.begin(); idx != paraBonds.end(); ++idx) {
      OBBond *bond = mol->GetBond(*idx);

      // make sure we didn't add this atom already from rule 1
      bool alreadyAdded = false;
      for (OBStereoUnitSet::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
        if ((*u2).type == OBStereo::CisTrans)
          if (bond->GetId() == (*u2).id) {
            alreadyAdded = true;
          }
      }
      if (alreadyAdded)
        continue;

      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();

      int beginClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond->GetBeginAtom());
      bool beginValid = false;
      switch (beginClassification) {
        case C12:
          beginValid = true;
          break;
        case C11:
          {
            // find the ligand
            OBAtom *ligandAtom = 0;
            FOR_NBORS_OF_ATOM (nbr, begin) {
              if ((nbr->GetIdx() != bond->GetBeginAtomIdx()) && (nbr->GetIdx() != bond->GetEndAtomIdx())) {
                ligandAtom = &*nbr;
                break;
              }
            }

            OBBitVec ligand = getFragment(ligandAtom, begin);
            for (OBStereoUnitSet::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
              if ((*u2).type == OBStereo::Tetrahedral) {
                if (ligand.BitIsSet((*u2).id))
                  beginValid = true;
              } else if((*u2).type == OBStereo::CisTrans) {
                OBBond *bond = mol->GetBondById((*u2).id);
                OBAtom *begin = bond->GetBeginAtom();
                OBAtom *end = bond->GetEndAtom();
                if (ligand.BitIsSet(begin->GetId()) || ligand.BitIsSet(end->GetId()))
                  beginValid = true;
              }
            }
          }
          break;
      }

      if (!beginValid)
        continue;

      int endClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond->GetEndAtom());
      bool endValid = false;
      switch (endClassification) {
        case C12:
          endValid = true;
          break;
        case C11:
          {
            // find the ligand
            OBAtom *ligandAtom = 0;
            FOR_NBORS_OF_ATOM (nbr, end) {
              if ((nbr->GetIdx() != bond->GetBeginAtomIdx()) && (nbr->GetIdx() != bond->GetEndAtomIdx())) {
                ligandAtom = &*nbr;
                break;
              }
            }

            OBBitVec ligand = getFragment(ligandAtom, end);
            for (OBStereoUnitSet::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
              if ((*u2).type == OBStereo::Tetrahedral) {
                if (ligand.BitIsSet((*u2).id))
                  endValid = true;
              } else if((*u2).type == OBStereo::CisTrans) {
                OBBond *bond = mol->GetBondById((*u2).id);
                OBAtom *begin = bond->GetBeginAtom();
                OBAtom *end = bond->GetEndAtom();
                if (ligand.BitIsSet(begin->GetId()) || ligand.BitIsSet(end->GetId()))
                  endValid = true;
              }
            }
          }
          break;
      }

      if (endValid)
        units.push_back(OBStereoUnit(OBStereo::CisTrans, bond->GetId(), true));
    }

    if (DEBUG) {
      for (OBStereoUnitSet::iterator unit = units.begin(); unit != units.end(); ++unit) {
        if (unit->type == OBStereo::Tetrahedral)
          cout << "Tetrahedral(center = " << unit->id << ", para = " << unit->para << ")" << endl;
        if (unit->type == OBStereo::CisTrans)
          cout << "CisTrans(bond = " << unit->id << ", para = " << unit->para << ")" << endl;
        if (unit->type == OBStereo::SquarePlanar)
          cout << "SquarePlanar(bond = " << unit->id << ", para = " << unit->para << ")" << endl;
      }
    }


    return units;
  }





















// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX






  /**
   * Helper function for FindStereogenicUnits using automorphisms.
   *
   * Find the duplicated symmetry class for neighbors of atom. This method only works if there is
   * only one duplicated symmetry class (i.e. T1123, T1112, T1111).
   */
  unsigned int findDuplicatedSymmetryClass(OBAtom *atom, const std::vector<unsigned int> &symClasses)
  {
    // find the duplicated symmetry class
    unsigned int duplicatedSymClass = OBGraphSym::NoSymmetryClass; // FIXME
    std::vector<unsigned int> nbrSymClasses;
    FOR_NBORS_OF_ATOM (nbr, atom) {
      nbrSymClasses.push_back(symClasses.at(nbr->GetIndex()));
    }
    for (std::size_t i = 0; i < nbrSymClasses.size(); ++i) {
      if (std::count(nbrSymClasses.begin(), nbrSymClasses.end(), nbrSymClasses.at(i)) >= 2) {
        duplicatedSymClass = nbrSymClasses.at(i);
        break;
      }
    }
    return duplicatedSymClass;
  }

  /**
   * Helper function for FindStereogenicUnits using automorphisms.
   *
   * Find the duplicated symmetry classes for neighbors of atom. This method only works for the
   * T1122 case.
   */
  void findDuplicatedSymmetryClasses(OBAtom *atom, const std::vector<unsigned int> &symClasses,
      unsigned int &duplicated1, unsigned int &duplicated2)
  {
    std::vector<unsigned int> nbrSymClasses;
    FOR_NBORS_OF_ATOM (nbr, atom)
      nbrSymClasses.push_back(symClasses.at(nbr->GetIndex()));
    std::sort(nbrSymClasses.begin(), nbrSymClasses.end());
    duplicated1 = nbrSymClasses[0];
    duplicated2 = nbrSymClasses[2];
  }

  /**
   * Helper function for FindStereogenicUnits using automorphisms.
   *
   * Find the duplicated symmetry classes for neighbors of atoms. This method works for all
   * cases (i.e. T1234, T1123, T1112, T1111 and T1122).
   */
  std::vector<unsigned int> findDuplicatedSymmetryClasses(OBAtom *atom, const std::vector<unsigned int> &symClasses)
  {
    std::vector<unsigned int> nbrSymClasses, result;
    FOR_NBORS_OF_ATOM (nbr, atom)
      nbrSymClasses.push_back(symClasses.at(nbr->GetIndex()));

    std::sort(nbrSymClasses.begin(), nbrSymClasses.end());
    for (std::size_t i = 0; i < nbrSymClasses.size(); ++i)
      if (std::count(nbrSymClasses.begin(), nbrSymClasses.end(), nbrSymClasses[i]) > 1)
        if (std::find(result.begin(), result.end(), nbrSymClasses[i]) == result.end())
          result.push_back(nbrSymClasses[i]);
    return result;
  }

  inline bool ComparePairSecond(const std::pair<unsigned int, unsigned int> &a,
      const std::pair<unsigned int, unsigned int> &b)
  {
    return (a.second < b.second);
  }



  /**
   * Helper functions for FindStereogenicUnits (using automorphisms).
   *
   * These functions determine if an automorphism permutation invert the
   * configuration of stereocenters by exchanging equivalent neighbor atoms
   * (i.e. neighbor atoms with the same topological symmetry class).
   *
   * @note: The molecule should be ordered by topological canonical labels.
   */
  struct StereoInverted {
    struct Entry {
      Automorphism p;
      std::vector<OBAtom*> invertedAtoms;
      std::vector<OBBond*> invertedBonds;
    };

    /**
     * Check if the specified automorphism causes an inversion of configuration
     * for the specified tetrahedral stereogenic center.
     */
    static bool permutationInvertsTetrahedralCenter(const Automorphism &p,
        OBAtom *center, const std::vector<unsigned int> &symmetry_classes,
        const std::vector<unsigned int> &canon_labels)
    {
      // Find the duplicated ligand symmetry class(es)
      std::vector<unsigned int> duplicatedSymClasses = findDuplicatedSymmetryClasses(center, symmetry_classes);

      if (DEBUG_INVERSIONS) {
        cout << "permutationInvertsTetrahedralCenter(" << center->GetIndex() << ")" << endl;
        print_vector("duplicatedSymClasses", duplicatedSymClasses);
      }

      std::vector< std::vector<OBAtom*> > duplicatedAtoms;

      int permutated = 0;
      for (std::size_t i = 0; i < duplicatedSymClasses.size(); ++i) {
        unsigned int duplicatedSymClass = duplicatedSymClasses[i];

        duplicatedAtoms.resize(duplicatedAtoms.size()+1);

        // Store the ligand indexes for the atoms with the duplicated symmetry class
        std::vector< std::pair<unsigned int, unsigned int> > tlist1;
        FOR_NBORS_OF_ATOM (nbr, center) {
          if (symmetry_classes[nbr->GetIndex()] == duplicatedSymClass) {
            tlist1.push_back(std::make_pair(nbr->GetIndex(), canon_labels[nbr->GetIndex()]));
            duplicatedAtoms.back().push_back(&*nbr);
          }
        }
        // Sort the indexes
        std::sort(tlist1.begin(), tlist1.end(), ComparePairSecond);

        //if (DEBUG_INVERSIONS) print_vector("tlist 1", tlist1);

        // Translate the sorted indexes using the automorphism
        std::vector<unsigned long> tlist2;
        for (std::size_t j = 0; j < tlist1.size(); ++j) {
          unsigned int t;
          if (MapsTo(p, tlist1[j].first, t))
            tlist2.push_back(canon_labels[t]);
        }

        if (DEBUG_INVERSIONS) print_vector("tlist 2", tlist2);

        // Permute the flag
        if (OBStereo::NumInversions(tlist2) % 2)
          //permutated = !permutated;
          permutated++;
      }

      if (permutated == 2) {
        std::vector<OBRing*> lssr = center->GetParent()->GetLSSR();
        assert( duplicatedAtoms.size() == 2 );
        assert( duplicatedAtoms[0].size() == 2 );
        assert( duplicatedAtoms[1].size() == 2 );
        for (std::size_t i = 0; i < lssr.size(); ++i) {
          if (lssr[i]->_pathset.BitIsSet(duplicatedAtoms[0][0]->GetIdx()) &&
              lssr[i]->_pathset.BitIsSet(duplicatedAtoms[0][1]->GetIdx()))
            return false;
          if (lssr[i]->_pathset.BitIsSet(duplicatedAtoms[1][0]->GetIdx()) &&
              lssr[i]->_pathset.BitIsSet(duplicatedAtoms[1][1]->GetIdx()))
            return false;
        }
        return true;
      }

      return permutated;
    }

    static bool permutationInvertsCisTransBeginOrEndAtom(const Automorphism &p, OBBond *bond, OBAtom *beginOrEnd,
        const std::vector<unsigned int> &canon_labels)
    {
      OBAtom *otherAtom = bond->GetNbrAtom(beginOrEnd);

      std::vector< std::pair<unsigned int, unsigned int> > tlist1;
      // Store the neighbor indexes in tlist1
      FOR_NBORS_OF_ATOM (nbr, beginOrEnd) {
        // skip the other double bond atom
        if (nbr->GetId() == otherAtom->GetId())
          continue;
        tlist1.push_back(std::make_pair(nbr->GetIndex(), canon_labels[nbr->GetIndex()]));
      }
      // Sort the indexes
      std::sort(tlist1.begin(), tlist1.end(), ComparePairSecond);

      // Translate the sorted indexes using the automorphism
      std::vector<unsigned long> tlist2;
      for (std::size_t j = 0; j < tlist1.size(); ++j) {
        unsigned int t;
        if (MapsTo(p, tlist1[j].first, t))
          tlist2.push_back(canon_labels[t]);
      }

      return (OBStereo::NumInversions(tlist2) % 2);
    }

    /**
     * Check if the specified automorphism causes an inversion of configuration
     * for the specfied stereogenic double bond.
     */
    static bool permutationInvertsCisTransCenter(const Automorphism &p, OBBond *bond,
        const std::vector<unsigned int> &canon_labels)
    {
      // begin atom
      bool beginInverted = permutationInvertsCisTransBeginOrEndAtom(p, bond, bond->GetBeginAtom(), canon_labels);
      // end atom
      bool endInverted = permutationInvertsCisTransBeginOrEndAtom(p, bond, bond->GetEndAtom(), canon_labels);

      // combine result using xor operation
      if (beginInverted ^ endInverted)
        return true;
      return false;
    }

    /**
     * Perform the computation.
     */
    static std::vector<Entry> compute(OBMol *mol, const std::vector<unsigned int> &symClasses,
        const Automorphisms &automorphisms)
    {
      if (DEBUG_INVERSIONS) cout << "ENTER StereoInverted::compute()" << endl;

      // We need topological canonical labels for this
      std::vector<unsigned int> canon_labels;
      CanonicalLabels(mol, symClasses, canon_labels, OBBitVec(), 5, true);

      // the result
      std::vector<Entry> result;

      // make a list of stereogenic centers inverted by the automorphism permutations
      for (std::size_t i = 0; i < automorphisms.size(); ++i) {
        Entry entry;
        entry.p = automorphisms[i];

        if (DEBUG_INVERSIONS) cout << "----> Checking automorphism " << i+1 << endl;

        // Check the atoms
        std::vector<OBAtom*>::iterator ia;
        for (OBAtom *atom = mol->BeginAtom(ia); atom; atom = mol->NextAtom(ia)) {
          // consider only potential stereo centers
          if (!isPotentialTetrahedral(atom))
            continue;
          // add the atom to the inverted list if the automorphism inverses it's configuration
          if (permutationInvertsTetrahedralCenter(automorphisms[i], atom, symClasses, canon_labels))
            entry.invertedAtoms.push_back(atom);
        }

        // Check the bonds
        std::vector<OBBond*>::iterator ib;
        for (OBBond *bond = mol->BeginBond(ib); bond; bond = mol->NextBond(ib)) {
          // consider only potential stereo centers
          if (!isPotentialCisTrans(bond))
            continue;
          // add the bond to the inverted list if the automorphism inverses it's configuration
          if (permutationInvertsCisTransCenter(entry.p, bond, canon_labels))
            entry.invertedBonds.push_back(bond);
        }

        if (DEBUG_INVERSIONS) {
          cout << "automorphism " << i+1 << "     ";
          for (std::size_t j = 0; j < mol->NumAtoms(); ++j) {
            unsigned int t;
            if (MapsTo(entry.p, j, t)) {
              if (t < 10) {
                cout << " " << t << " ";
              } else {
                cout << t << " ";
              }
            }
          }
          cout << endl;
          cout << "  invertedAtoms: ";
          for (std::size_t l = 0; l < entry.invertedAtoms.size(); ++l)
            cout << entry.invertedAtoms[l]->GetId() << " ";
          cout << endl;
          cout << "  invertedBonds: ";
          for (std::size_t l = 0; l < entry.invertedBonds.size(); ++l)
            cout << entry.invertedBonds[l]->GetId() << " ";
          cout << endl;
        }

        result.push_back(entry);
      }

      if (DEBUG_INVERSIONS) cout << "EXIT StereoInverted::compute()" << endl;

      return result;
    }


  };

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //
  //  FindStereogenicUnits using automorphisms
  //
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  /**
   * Find an atom with the specified symmetry class. The first found atom is returned
   * or 0 when there is no such atom. This function is intended to be used in cases
   * where any atom with the specified symmetry class can be used. For example, when
   * checking a fragments for stereocenters, the result will be the same for any atom
   * with a specified (duplicated) symmetry class.
   */
  OBAtom* findAtomWithSymmetryClass(OBAtom *atom, unsigned int symClass, const std::vector<unsigned int> &symClasses)
  {
    OBAtom *ligandAtom = 0;
    FOR_NBORS_OF_ATOM (nbr, atom)
      if (symClasses.at(nbr->GetIndex()) == symClass)
        ligandAtom = &*nbr;
    return ligandAtom;
  }

  /**
   * Helper function to determine if a stereogenic center with duplicated symmetry classes
   * really is a stereogenic center.
   *
   * Check if the ligandAtom's fragment (see getFragment()) contains at least one
   * true- or 1 para-stereocenter. This is rule 1 (a & b) in the Razinger paper on
   * stereoisomer generation.
   */
  bool containsAtLeast_1true_1para(OBAtom *ligandAtom, OBAtom *skip, const OBStereoUnitSet &units)
  {
    OBMol *mol = skip->GetParent();
    // create the fragment bitvec
    OBBitVec ligand = getFragment(ligandAtom, skip);
    for (OBStereoUnitSet::const_iterator u2 = units.begin(); u2 != units.end(); ++u2) {
      if (isUnitInFragment(mol, *u2, ligand))
        return true;
    }
    return false;
  }

  /**
   * Helper function to determine if a stereogenic center with duplicated symmetry classes
   * really is a stereogenic center.
   *
   * Check if the ligandAtom's fragment (see getFragment()) contains at least one
   * true- or 2 para-stereocenter. This is rule 2a and rule 3 in the Razinger
   * paper on stereoisomer generation.
   */
  bool containsAtLeast_1true_2para(OBAtom *ligandAtom, OBAtom *atom, const OBStereoUnitSet &units)
  {
    OBMol *mol = atom->GetParent();
    // check if ligand contains at least:
    // - 1 true-stereocenter
    // - 2 para-stereocenters
    OBBitVec ligand = getFragment(ligandAtom, atom);
    bool foundTrueStereoCenter = false;
    int paraStereoCenterCount = 0;
    for (OBStereoUnitSet::const_iterator u2 = units.begin(); u2 != units.end(); ++u2) {
      if (isUnitInFragment(mol, *u2, ligand)) {
        if ((*u2).para) {
          paraStereoCenterCount++;
        } else {
          foundTrueStereoCenter = true;
        }
      }
    }

    if (foundTrueStereoCenter || paraStereoCenterCount >= 2)
      return true;
    if (ligandAtom->IsInRing() && atom->IsInRing() && paraStereoCenterCount)
      return true;
    return false;
  }

  /**
   * Helper function to determine if a stereogenic center with duplicated symmetry classes
   * really is a stereogenic center.
   *
   * Check if the ligandAtom's fragment (see getFragment()) contains at least one
   * true- or 2 separate assemblies of at least 2 para-stereocenter. This is rule
   * 2b in the Razinger paper on stereoisomer generation.
   */
  bool containsAtLeast_2true_2paraAssemblies(OBAtom *ligandAtom, OBAtom *atom, const OBStereoUnitSet &units, const std::vector<OBBitVec> &mergedRings)
  {
    OBMol *mol = atom->GetParent();
    // check if ligand contains at least:
    // - 2 true-stereocenter
    // - 2 separate para-stereocenters assemblies
    OBBitVec ligand = getFragment(ligandAtom, atom);
    int trueStereoCenterCount = 0;
    std::vector<unsigned int> ringIndices;
    for (OBStereoUnitSet::const_iterator u2 = units.begin(); u2 != units.end(); ++u2) {
      if ((*u2).type == OBStereo::Tetrahedral) {
        if (ligand.BitIsSet((*u2).id)) {
          if ((*u2).para) {
            OBAtom *paraAtom = mol->GetAtomById((*u2).id);
            for (std::size_t ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
              if (mergedRings.at(ringIdx).BitIsSet(paraAtom->GetIdx()))
                if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end())
                  ringIndices.push_back(ringIdx);
            }
          } else
            trueStereoCenterCount++;
        }
      } else if((*u2).type == OBStereo::CisTrans) {
        OBBond *bond = mol->GetBondById((*u2).id);
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();
        if (ligand.BitIsSet(begin->GetId()) || ligand.BitIsSet(end->GetId())) {
          if ((*u2).para) {
            for (std::size_t ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
              if (mergedRings.at(ringIdx).BitIsSet(begin->GetIdx()) || mergedRings.at(ringIdx).BitIsSet(end->GetIdx())) {
                if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end()) {
                  ringIndices.push_back(ringIdx);
                }
              }
            }
          } else {
            trueStereoCenterCount++;
          }
        }
      }
    }

    if (trueStereoCenterCount >= 2 || ringIndices.size() >= 2)
      return true;
    return false;
  }

  /**
   * Find the stereogenic units in a molecule using automorphisms.
   *
   * This is a public function: see header for details.
   */
  OBStereoUnitSet FindStereogenicUnits(OBMol *mol,
      const std::vector<unsigned int> &symClasses, const Automorphisms &automorphisms)
  {
    OBStereoUnitSet units;

    // do quick test to see if there are any possible stereogenic units
    if (!mayHaveTetrahedralCenter(mol) && !mayHaveCisTransBond(mol))
      return units;

    // make sure we have symmetry classes for all atoms
    if (symClasses.size() != mol->NumAtoms())
      return units;

    // Compute which automorphisms cause inversion of configuration
    // for the stereogenic units
    std::vector<StereoInverted::Entry> inverted = StereoInverted::compute(mol, symClasses, automorphisms);

    std::vector<OBBitVec> mergedRings = mergeRings(mol, symClasses);

    std::vector<unsigned long> doneAtoms, doneBonds;
    unsigned int lastSize = units.size();
    while (true) {
      std::vector<OBAtom*>::iterator ia;
      for (OBAtom *atom = mol->BeginAtom(ia); atom; atom = mol->NextAtom(ia)) {
        if (std::find(doneAtoms.begin(), doneAtoms.end(), atom->GetId()) != doneAtoms.end())
          continue;
        // consider only potential steroecenters
        if (!isPotentialTetrahedral(atom))
          continue;

        // A potential stereocenter is really a stereocenter if there exists no automorphic
        // permutation causing an inversion of the configuration of only the potential
        // stereogenic unit under consideration.
        bool foundPermutation = false; // invert __only__ configuration of atom
        for (std::size_t i = 0; i < inverted.size(); ++i) {
          const std::vector<OBAtom*> &atoms = inverted[i].invertedAtoms;
          if (atoms.size() != 1)
            continue;
          const std::vector<OBBond*> &bonds = inverted[i].invertedBonds;
          if (bonds.size())
            continue;
          if (atoms[0] == atom) {
            foundPermutation = true;
            break;
          }
        }

        int classification = classifyTetrahedralNbrSymClasses(symClasses, atom);

        if (DEBUG_INVERSIONS)
          cout << "foundPermutation for id = " << atom->GetId() << ": " << foundPermutation << endl;

        if (!foundPermutation) {
          // true-stereocenter found
          bool isParaCenter = (classification == T1234) ? false : true;
          //cout << "found(2) " << atom->GetId() << endl;
          units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId(), isParaCenter));
          doneAtoms.push_back(atom->GetId());
        } else {
          // count ligand configurations:
          // If there exists at least one automorphic permutation causing the inversion of the
          // configuration of only the stereogenic unit under consideration, then the potential
          // stereocenter can be a stereocenter if the number of topologically equivalent neighbors
          // (ligands) of potential stereogenic is less than or equal to the number of configurations
          // of these ligands.
          //
          // In practise:
          //    T1123 -> 1 true stereocenter OR 2 para stereocenters
          //    T1122 -> 1 true stereocenter OR 2 para stereocenters (for both)
          //    T1112 -> 2 true stereocenters OR 2 para stereocenter assemblies
          //    T1111 -> 2 true stereocenters OR 2 para stereocenter assemblies
          switch (classification) {
            case T1123:
              {
                unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses);
                OBAtom *ligandAtom = findAtomWithSymmetryClass(atom, duplicatedSymClass, symClasses);
                if (containsAtLeast_1true_2para(ligandAtom, atom, units)) {
                  units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId(), true));
                  doneAtoms.push_back(atom->GetId());
                }
              }
              break;
            case T1122:
              {
                unsigned int duplicatedSymClass1, duplicatedSymClass2;
                findDuplicatedSymmetryClasses(atom, symClasses, duplicatedSymClass1, duplicatedSymClass2);
                OBAtom *ligandAtom1 = findAtomWithSymmetryClass(atom, duplicatedSymClass1, symClasses);
                OBAtom *ligandAtom2 = findAtomWithSymmetryClass(atom, duplicatedSymClass2, symClasses);
                if (containsAtLeast_1true_2para(ligandAtom1, atom, units) &&
                    containsAtLeast_1true_2para(ligandAtom2, atom, units)) {
                  units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId(), true));
                  doneAtoms.push_back(atom->GetId());
                }
              }
              break;
            case T1112:
            case T1111:
              {
                unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses);
                OBAtom *ligandAtom = findAtomWithSymmetryClass(atom, duplicatedSymClass, symClasses);
                if (containsAtLeast_2true_2paraAssemblies(ligandAtom, atom, units, mergedRings)) {
                  units.push_back(OBStereoUnit(OBStereo::Tetrahedral, atom->GetId(), true));
                  doneAtoms.push_back(atom->GetId());
                }
              }
              break;
          }
        }
      }

      std::vector<OBBond*>::iterator ib;
      for (OBBond *bond = mol->BeginBond(ib); bond; bond = mol->NextBond(ib)) {
        if (std::find(doneBonds.begin(), doneBonds.end(), bond->GetId()) != doneBonds.end())
          continue;
        if (!isPotentialCisTrans(bond))
          continue;

        // A double bond is a stereogenic bond if there exists no automorphic
        // permutation causing an inversion of the configuration of only the potential
        // stereogenic unit under consideration.
        bool foundPermutation = false; // invert __only__ configuration of atom
        for (std::size_t i = 0; i < inverted.size(); ++i) {
          const std::vector<OBAtom*> &atoms = inverted[i].invertedAtoms;
          // if any atoms are inverted, the bond can't be the only inverted stereocenter
          if (atoms.size())
            continue;
          const std::vector<OBBond*> &bonds = inverted[i].invertedBonds;
          // the bond should be the only inverted stereocenter
          if (bonds.size() != 1)
            continue;
          // check if it is this bond
          if (bonds[0] == bond) {
            foundPermutation = true;
            break;
          }
        }

        int beginClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond->GetBeginAtom());
        int endClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond->GetEndAtom());

        if (!foundPermutation) {
          // true-stereocenter found
          bool isParaCenter = (beginClassification == C12) && (endClassification == C12) ? false : true;
          units.push_back(OBStereoUnit(OBStereo::CisTrans, bond->GetId(), isParaCenter));
          doneBonds.push_back(bond->GetId());
        } else {
          // count ligand configurations:
          bool beginValid = false;
          switch (beginClassification) {
            case C12:
              beginValid = true;
              break;
            case C11:
              {
                // find the ligand
                OBAtom *ligandAtom = 0;
                FOR_NBORS_OF_ATOM (nbr, bond->GetBeginAtom()) {
                  if ((nbr->GetIdx() != bond->GetBeginAtomIdx()) && (nbr->GetIdx() != bond->GetEndAtomIdx())) {
                    ligandAtom = &*nbr;
                    break;
                  }
                }
                if (ligandAtom)
                  beginValid = containsAtLeast_1true_1para(ligandAtom, bond->GetBeginAtom(), units);
              }
              break;
          }

          if (!beginValid)
            continue;

          bool endValid = false;
          switch (endClassification) {
            case C12:
              endValid = true;
              break;
            case C11:
              {
                // find the ligand
                OBAtom *ligandAtom = 0;
                FOR_NBORS_OF_ATOM (nbr, bond->GetEndAtom()) {
                  if ((nbr->GetIdx() != bond->GetBeginAtomIdx()) && (nbr->GetIdx() != bond->GetEndAtomIdx())) {
                    ligandAtom = &*nbr;
                    break;
                  }
                }
                if (ligandAtom)
                  endValid = containsAtLeast_1true_1para(ligandAtom, bond->GetEndAtom(), units);
              }
              break;
          }

          if (endValid) {
            units.push_back(OBStereoUnit(OBStereo::CisTrans, bond->GetId(), true));
            doneBonds.push_back(bond->GetId());
          }
        }
      }


      if (units.size() == lastSize)
        break;
      lastSize = units.size();
    }

    if (DEBUG) {
      for (OBStereoUnitSet::iterator unit = units.begin(); unit != units.end(); ++unit) {
        if (unit->type == OBStereo::Tetrahedral)
          cout << "Tetrahedral(center = " << unit->id << ", para = " << unit->para << ")" << endl;
        if (unit->type == OBStereo::CisTrans)
          cout << "CisTrans(bond = " << unit->id << ", para = " << unit->para << ")" << endl;
        if (unit->type == OBStereo::SquarePlanar)
          cout << "SquarePlanar(bond = " << unit->id << ", para = " << unit->para << ")" << endl;
      }
    }

    return units;
  } // FindStereogenicUnits using automorphisms


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
  ////////////////////////////////////////////////////////////////////////////
  //
  //
  //  From0D
  //
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void StereoFrom0D(OBMol *mol)
  {
    if (mol->HasChiralityPerceived())
      return;

    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom0D", obAuditMsg);

    std::vector<unsigned int> symmetry_classes = FindSymmetry(mol);
    OBStereoUnitSet stereogenicUnits = FindStereogenicUnits(mol, symmetry_classes);

    TetrahedralFrom0D(mol, stereogenicUnits);
    CisTransFrom0D(mol, stereogenicUnits);
    mol->SetChiralityPerceived();
  }

  std::vector<OBTetrahedralStereo*> TetrahedralFrom0D(OBMol *mol,
      const OBStereoUnitSet &stereoUnits, bool addToMol)
  {
    std::vector<OBTetrahedralStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::TetrahedralFrom0D", obAuditMsg);

    // Delete any existing stereo objects that are not a member of 'centers'
    // and make a map of the remaining ones
    std::map<unsigned long, OBTetrahedralStereo*> existingMap;
    std::vector<OBGenericData*>::iterator data;
    std::vector<OBGenericData*> stereoData = mol->GetAllData(OBGenericDataType::StereoData);
    for (data = stereoData.begin(); data != stereoData.end(); ++data) {
      if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        unsigned long center = ts->GetConfig().center;
        // check if the center is really stereogenic
        bool isStereogenic = false;
        OBStereoUnitSet::const_iterator u;
        for (u = stereoUnits.begin(); u != stereoUnits.end(); ++u) {
          if ((*u).type == OBStereo::Tetrahedral)
            if ((*u).id == center)
              isStereogenic = true;
        }

        if (isStereogenic) {
          existingMap[center] = ts;
          configs.push_back(ts);
        } else {
          // According to OpenBabel, this is not a tetrahedral stereo
          obErrorLog.ThrowError(__FUNCTION__, "Removed spurious TetrahedralStereo object", obAuditMsg);
          mol->DeleteData(ts);
        }
      }
    }

    OBStereoUnitSet::const_iterator u;
    for (u = stereoUnits.begin(); u != stereoUnits.end(); ++u) {
      // skip non-tetrahedral units
      if ((*u).type != OBStereo::Tetrahedral)
        continue;
      // if there already exists a OBTetrahedralStereo object for this
      // center, continue
      if (existingMap.find((*u).id) != existingMap.end())
        continue;

      OBAtom *center = mol->GetAtomById((*u).id);

      OBTetrahedralStereo::Config config;
      config.specified = false;
      config.center = (*u).id;
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
      const OBStereoUnitSet &stereoUnits,
      bool addToMol)
  {
    std::vector<OBCisTransStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::CisTransFrom0D", obAuditMsg);

    std::vector<unsigned long> bonds;
    for (OBStereoUnitSet::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::CisTrans)
        bonds.push_back((*u).id);

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
      // bond, leave it alone unless it's in a ring of small size

      bool alreadyExists = (existingMap.find(*i) != existingMap.end());
      OBBond *bond = mol->GetBondById(*i);

      OBCisTransStereo *ct;
      OBCisTransStereo::Config config;
      if (alreadyExists)
      {
        ct = existingMap[*i];
        config = ct->GetConfig();
      }
      else
      {
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();

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

        ct = new OBCisTransStereo(mol);
        ct->SetConfig(config);
      }

      // For a double bond in a ring of size IMPLICIT_CIS_RING_SIZE or less
      // the stereochemistry is implicitly cis (in terms
      // of the ring atoms)
      OBRing* ring = bond->FindSmallestRing();
      if (ring && ring->Size() <= IMPLICIT_CIS_RING_SIZE) {

        // Find the ring atoms in the config.refs
        vector<unsigned int> ringrefs(2);
        for (int i = 0; i<2; ++i) {
          if (config.refs[i*2] != OBStereo::ImplicitRef && ring->IsMember(mol->GetAtomById(config.refs[i*2])))
            ringrefs[i] = config.refs[i*2];
          else
            ringrefs[i] = config.refs[i*2 + 1];
        }
        if (!ct->IsCis(ringrefs[0], ringrefs[1])) // Need to invert the stereo
          config.shape = OBStereo::ShapeZ;

        config.specified = true;
        ct->SetConfig(config);
      }

      configs.push_back(ct);
      // add the data to the molecule if needed
      if (addToMol && !alreadyExists)
        mol->SetData(ct);

    }

    return configs;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //
  //
  //  From3D
  //
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void StereoFrom3D(OBMol *mol, bool force)
  {
    if (mol->HasChiralityPerceived() && !force)
      return;

    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom3D", obAuditMsg);

    std::vector<unsigned int> symmetry_classes = FindSymmetry(mol);
    OBStereoUnitSet stereogenicUnits = FindStereogenicUnits(mol, symmetry_classes);

    mol->DeleteData(OBGenericDataType::StereoData);
    TetrahedralFrom3D(mol, stereogenicUnits);
    CisTransFrom3D(mol, stereogenicUnits);
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
      const OBStereoUnitSet &stereoUnits, bool addToMol)
  {
    std::vector<OBTetrahedralStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::TetrahedralFrom3D", obAuditMsg);

    // find all tetrahedral centers
    std::vector<unsigned long> centers;
    for (OBStereoUnitSet::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::Tetrahedral)
        centers.push_back((*u).id);

    std::vector<unsigned long>::iterator i;
    for (i = centers.begin(); i != centers.end(); ++i) {
      OBAtom *center = mol->GetAtomById(*i);

      // make sure we have at least 3 heavy atom neighbors
      // timvdm 28 Jun 2009: This is already checked in FindStereogenicUnits
      if (center->GetHvyDegree() < 3) {
        std::stringstream errorMsg;
        errorMsg << "Cannot calculate a signed volume for an atom with a heavy atom valence of "
                 << center->GetHvyDegree() << std::endl;
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
      // and check for a bond that indicates unspecified stereochemistry
      std::vector<vector3> nbrCoords;
      OBAtom *from = mol->GetAtomById(config.from);
      OBBond *bond = mol->GetBond(from, center);
      if (bond->IsWedgeOrHash() && bond->GetBeginAtom()==center)
        config.specified = false;

      nbrCoords.push_back(from->GetVector());
      for (OBStereo::RefIter id = config.refs.begin(); id != config.refs.end(); ++id) {
        OBAtom *nbr = mol->GetAtomById(*id);
        nbrCoords.push_back(nbr->GetVector());
        OBBond *bond = mol->GetBond(nbr, center);
        if (bond->IsWedgeOrHash() && bond->GetBeginAtom()==center)
          config.specified = false;
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
      const OBStereoUnitSet &stereoUnits, bool addToMol)
  {
    std::vector<OBCisTransStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::CisTransFrom3D", obAuditMsg);

    // find all cis/trans bonds
    std::vector<unsigned long> bonds;
    for (OBStereoUnitSet::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::CisTrans)
        bonds.push_back((*u).id);

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

      // 0      3     Get signed distance of 0 and 2 to the plane
      //  \    /      that goes through the double bond and is at
      //   C==C       right angles to the stereo bonds.
      //  /    \      If the two signed distances have the same sign
      // 1      2     then they are cis; if not, then trans.

      vector3 dbl_bond = end->GetVector() - begin->GetVector();
      vector3 above_plane = cross(dbl_bond, bondVecs[0]);
      double d0 = Point2PlaneSigned( mol->GetAtomById(config.refs[0])->GetVector(),
                               begin->GetVector(), end->GetVector(), above_plane);
      double d2 = Point2PlaneSigned( mol->GetAtomById(config.refs[2])->GetVector(),
                               begin->GetVector(), end->GetVector(), above_plane);

      if ((d0 > 0 && d2 > 0) || (d0 < 0 && d2 < 0))
        config.shape = OBStereo::ShapeZ;
      else
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
  ////////////////////////////////////////////////////////////////////////////
  //
  //  From2D
  //
  //  Reference:
  //  [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the
  //  Unambiguous Identification of the Stereochemical Characteristics of
  //  Compounds During Their Registration in Databases. Molecules 2000, 6,
  //  915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
  //
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void StereoFrom2D(OBMol *mol, std::map<OBBond*, enum OBStereo::BondDirection> *updown, bool force)
  {
    if (mol->HasChiralityPerceived() && !force)
      return;

    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom2D", obAuditMsg);

    std::vector<unsigned int> symmetry_classes = FindSymmetry(mol);
    OBStereoUnitSet stereogenicUnits = FindStereogenicUnits(mol, symmetry_classes);

    mol->DeleteData(OBGenericDataType::StereoData);
    TetrahedralFrom2D(mol, stereogenicUnits);
    CisTransFrom2D(mol, stereogenicUnits, updown);
    mol->SetChiralityPerceived();
  }
  //! Calculate the "sign of a triangle" given by a set of 3 2D coordinates
  double TriangleSign(const vector3 &a, const vector3 &b, const vector3 &c)
  {
    // equation 6 from [1]
    return (a.x() - c.x()) * (b.y() - c.y()) - (a.y() - c.y()) * (b.x() - c.x());
  }
  //! Calculate whether three vectors are arranged in order of increasing
  //! angle anticlockwise (true) or clockwise (false) relative to a central point.
  bool AngleOrder(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &center)
  {
    vector3 t, u, v;
    t = a - center;
    t.normalize();
    u = b - center;
    u.normalize();
    v = c - center;
    v.normalize();
    return TriangleSign(t, u, v) > 0;
  }
  //! Get the angle between three atoms (from -180 to +180)
  //! Note: OBAtom.GetAngle just returns 0->180
  double GetAngle(OBAtom *a, OBAtom *b, OBAtom *c)
  {
   vector3 v1,v2;

    v1 = a->GetVector() - b->GetVector();
    v2 = c->GetVector() - b->GetVector();
    if (IsNearZero(v1.length(), 1.0e-3)
      || IsNearZero(v2.length(), 1.0e-3)) {
        return(0.0);
    }

    double angle = (atan2(v2.y(),v2.x()) - atan2(v1.y(),v1.x())) * RAD_TO_DEG;
    while (angle < -180.0) angle += 360.0;
    while (angle > 180.0) angle -= 360.0;
    return angle;
  }
  std::vector<OBTetrahedralStereo*> TetrahedralFrom2D(OBMol *mol,
      const OBStereoUnitSet &stereoUnits, bool addToMol)
  {
    std::vector<OBTetrahedralStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::TetrahedralFrom2D", obAuditMsg);

    // find all tetrahedral centers
    std::vector<unsigned long> centers;
    for (OBStereoUnitSet::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::Tetrahedral)
        centers.push_back((*u).id);


    std::vector<unsigned long>::iterator i;
    for (i = centers.begin(); i != centers.end(); ++i) {
      OBAtom *center = mol->GetAtomById(*i);

      // make sure we have at least 3 heavy atom neighbors
      if (center->GetHvyDegree() < 3) {
        std::stringstream errorMsg;
        errorMsg << "Cannot calculate a signed volume for an atom with a heavy atom valence of "
                 << center->GetHvyDegree() << std::endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        continue;
      }

      OBTetrahedralStereo::Config config;
      config.center = *i;

      // We assume the 'tip-only' convention. That is, wedge or hash bonds only
      // determine the stereochemistry at their thin end (the BeginAtom)
      bool tiponly = true;

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
            if (tiponly)
              planeAtoms.push_back(nbr);
            else
              wedgeAtoms.push_back(nbr);
          }
        } else if (bond->IsWedge()) {
          // wedge bonds
          if (bond->GetBeginAtom()->GetId() == center->GetId()) {
            // this is a 'real' wedge bond going from center to nbr
            wedgeAtoms.push_back(nbr);
          } else {
            // this is an 'inverted' wedge bond going from nbr to center
            if (tiponly)
              planeAtoms.push_back(nbr);
            else
              hashAtoms.push_back(nbr);
          }
        } else if (bond->IsWedgeOrHash()) {
          if (!tiponly || (tiponly && bond->GetBeginAtom()->GetId() == center->GetId())) {
            config.specified = true;
            config.winding = OBStereo::UnknownWinding;
            break;
          }
          else
            planeAtoms.push_back(nbr);
        } else {
          // plane bonds
          planeAtoms.push_back(nbr);
        }
      }

      // Handle the case of a tet center with four plane atoms or
      //        3 plane atoms with the fourth bond implicit
      if (planeAtoms.size() == 4 || (planeAtoms.size() == 3 && center->GetExplicitDegree()==3))
        config.specified = false;

      bool success = true;

      using namespace std;
      if (!config.specified || (config.specified && config.winding==OBStereo::UnknownWinding)) {
        // unspecified or specified as unknown
        FOR_NBORS_OF_ATOM (nbr, center)
          if (config.from == OBStereo::NoRef)
            config.from = nbr->GetId();
          else
            config.refs.push_back(nbr->GetId());
        while (config.refs.size() < 3)
          config.refs.push_back(OBStereo::ImplicitRef);
      } else {

        // config.specified

        if (hashAtoms.size() == 4 || wedgeAtoms.size() == 4)
        {
          success = false;
        } else if (planeAtoms.size() + hashAtoms.size() + wedgeAtoms.size() == 4)
        {
          // Handle all explicit tetra with at least one stereobond
          vector<OBAtom*> order;

          // First of all, handle the case of three wedge (or three hash) and one other bond
          //          by converting it into a single hash (or single wedge) and three planes
          if (wedgeAtoms.size() == 3 || hashAtoms.size() == 3) {
            vector<OBAtom*> *pwedge, *phash;
            if (wedgeAtoms.size() == 3) {
              pwedge = &wedgeAtoms; phash = &hashAtoms;
            }
            else {
              phash = &wedgeAtoms; pwedge = &hashAtoms;
            }
            if (planeAtoms.size() == 0) { // Already has the hash bond
              planeAtoms.insert(planeAtoms.end(), pwedge->begin(), pwedge->end());
              pwedge->clear();
            }
            else { // Does not already have the hash bond
              phash->push_back(planeAtoms[0]);
              planeAtoms.clear();
              pwedge->clear();
            }
          }

          // Pick a stereobond on which to base the stereochemistry:
          bool wedge = wedgeAtoms.size() > 0;
          order.push_back(wedge?wedgeAtoms[0]:hashAtoms[0]);
          vector<OBAtom*> nbrs;
          FOR_NBORS_OF_ATOM(nbr, center) {
            if (&*nbr != order[0])
              nbrs.push_back(&*nbr);
          }
          // Add "nbrs" to "order" in order of anticlockwise stereo
          order.push_back(nbrs[0]);
          if (AngleOrder(order[0]->GetVector(), order[1]->GetVector(), nbrs[1]->GetVector(), center->GetVector()))
            order.push_back(nbrs[1]);
          else
            order.insert(order.begin()+1, nbrs[1]);
          if (AngleOrder(order[0]->GetVector(), order[2]->GetVector(), nbrs[2]->GetVector(), center->GetVector()))
            order.push_back(nbrs[2]);
          else {
            if (AngleOrder(order[0]->GetVector(), order[1]->GetVector(), nbrs[2]->GetVector(), center->GetVector()))
              order.insert(order.begin()+2, nbrs[2]);
            else
              order.insert(order.begin()+1, nbrs[2]);
          }

          // Handle the case of two planes with a wedge and hash bond opposite each other.
          // This is handled as in the InChI TechMan (Figure 9) by marking it ambiguous if
          // the (small) angle between the plane bonds is > 133, and basing the stereo on
          // the 'inner' bond otherwise. This is commonly used for stereo in rings.
          // See also Get2DTetrahedralAmbiguity() in ichister.c (part of InChI)
          if (planeAtoms.size() == 2 && wedgeAtoms.size() == 1) { // Two planes, 1 wedge, 1 hash
            if (order[2] == hashAtoms[0]) { // The wedge and hash are opposite
              double angle = GetAngle(order[1], center, order[3]); // The anticlockwise angle between the plane atoms
              if (angle > -133 && angle < 133) { // This value is from the InChI TechMan Figure 9
                if (angle > 0) { // Change to three planes and the hash bond
                  std::rotate(order.begin(), order.begin() + 2, order.end()); // Change the order so that it begins with the hash bond
                  wedge = false;
                  planeAtoms.push_back(wedgeAtoms[0]);
                  wedgeAtoms.clear();
                }
                else { // Change to three planes and the wedge bond (note: order is already correct)
                  planeAtoms.push_back(hashAtoms[0]);
                  hashAtoms.clear();
                }
              } // No need for "else" statement, as this will be picked up as ambiguous stereo below
            }
          }

          config.from = order[0]->GetId();
          config.refs.resize(3);
          for(int i=0; i<3; ++i)
            config.refs[i] = order[i+1]->GetId();
          if (wedge)
            config.winding = OBStereo::AntiClockwise;

          // Check for ambiguous stereo based on the members of "order".
          // If the first is a wedge bond, then the next should be a plane/hash, then plane/wedge, then plane/hash
          // If not, then the stereo is considered ambiguous.
          vector<OBAtom*> *pwedge, *phash;
          if (wedge) {
            pwedge = &wedgeAtoms; phash = &hashAtoms;
          }
          else {
            phash = &wedgeAtoms; pwedge = &hashAtoms;
          }
          if (std::find(pwedge->begin(), pwedge->end(), order[1]) != pwedge->end() ||
              std::find(phash->begin(), phash->end(), order[2]) != phash->end()  ||
              std::find(pwedge->begin(), pwedge->end(), order[3]) != pwedge->end()) { // Ambiguous stereo
                success = false;
          }
        } else // 3 explicit bonds from here on
          if(hashAtoms.size() == 0 || wedgeAtoms.size() == 0) {
            // Composed of just wedge bonds and plane bonds, or just hash bonds and plane bonds

            // Pick a stereobond on which to base the stereochemistry:
            vector<OBAtom*> order;
            bool wedge = wedgeAtoms.size() > 0;
            order.push_back(wedge?wedgeAtoms[0]:hashAtoms[0]);
            vector<OBAtom*> nbrs;
            FOR_NBORS_OF_ATOM(nbr, center) {
              if (&*nbr != order[0])
                nbrs.push_back(&*nbr);
            }
            // Add "nbrs" to "order" in order of anticlockwise stereo
            order.push_back(nbrs[0]);
            if (AngleOrder(order[0]->GetVector(), order[1]->GetVector(), nbrs[1]->GetVector(), center->GetVector()))
              order.push_back(nbrs[1]);
            else
              order.insert(order.begin()+1, nbrs[1]);

            // Handle the case of two planes with a wedge/hash in the small angle between them.
            // This is handled similar to the InChI TechMan (Figure 10) by treating the stereo bond
            // as being in the large angle. This is consistent with Symyx Draw.
            if (planeAtoms.size() == 2) { // Two planes, 1 stereo
              double angle = GetAngle(order[1], center, order[2]); // The anticlockwise angle between the plane atoms
              if (angle < 0) // Invert the stereo of the stereobond
                wedge = !wedge;
            }

            config.from = OBStereo::ImplicitRef;
            config.refs.resize(3);
            for(int i=0; i<3; ++i)
              config.refs[i] = order[i]->GetId();
            if (!wedge)
              config.winding = OBStereo::AntiClockwise;
        }  else { // 3 explicit bonds with at least one hash and at least one wedge
          success = false;
        }
      } // end of config.specified

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
      const OBStereoUnitSet &stereoUnits,
      const std::map<OBBond*, enum OBStereo::BondDirection> *updown, bool addToMol)
  {
    std::vector<OBCisTransStereo*> configs;
    std::map<OBBond*, enum OBStereo::BondDirection>::const_iterator ud_cit;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::CisTransFrom2D", obAuditMsg);

    // find all cis/trans bonds
    std::vector<unsigned long> bonds;
    for (OBStereoUnitSet::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::CisTrans)
        bonds.push_back((*u).id);

    std::vector<unsigned long>::iterator i;
    for (i = bonds.begin(); i != bonds.end(); ++i) {
      OBBond *bond = mol->GetBondById(*i);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();

      // Create a vector with the coordinates of the neighbor atoms
      std::vector<vector3> bondVecs;
      OBCisTransStereo::Config config;
      config.specified = true;

      // begin
      config.begin = begin->GetId();
      FOR_NBORS_OF_ATOM (nbr, begin) {
        if (nbr->GetId() == end->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
        bondVecs.push_back(nbr->GetVector());

        // Check whether a single bond with unknown dir starts at the dbl bond (tip-only convention)
        OBBond *b = mol->GetBond(begin, &*nbr);
        if (updown) {
          ud_cit = updown->find(b);
          if (ud_cit!=updown->end() && ud_cit->second==OBStereo::UnknownDir && b->GetBeginAtom()==begin)
            config.specified = false;
        }
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

        // Check whether a single bond with unknown dir starts at the dbl bond (tip-only convention)
        OBBond *b = mol->GetBond(end, &*nbr);
        if (updown) {
          ud_cit = updown->find(b);
          if (ud_cit!=updown->end() && ud_cit->second==OBStereo::UnknownDir && b->GetBeginAtom()==end)
            config.specified = false;
        }
      }
      if (config.refs.size() == 3) {
        config.refs.push_back(OBStereo::ImplicitRef);
        vector3 pos;
        end->GetNewBondVector(pos, 1.0);
        bondVecs.push_back(pos);
      }

      // Handle the case where the dbl bond is marked as unknown stereo
      if (updown) {
        ud_cit = updown->find(bond);
        if (ud_cit!=updown->end() && ud_cit->second==OBStereo::UnknownDir)
            config.specified = false;
      }

      if (config.specified==true) { // Work out the stereochemistry
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

  bool TetStereoToWedgeHash(OBMol &mol,
      std::map<OBBond*, enum OBStereo::BondDirection> &updown,
      std::map<OBBond*, OBStereo::Ref> &from)
  {
    // Store the tetcenters for the second loop (below)
    std::set <unsigned long> tetcenters;
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo::Config cfg = ts->GetConfig();
        tetcenters.insert(cfg.center);
      }

    // This loop sets one bond of each tet stereo to up or to down (2D only)
    std::set <OBBond *> alreadyset;
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo::Config cfg = ts->GetConfig();

        if (cfg.specified) {
          OBBond* chosen = (OBBond*) NULL;
          OBAtom* center = mol.GetAtomById(cfg.center);

          // Find the two bonds closest in angle and remember them if
          // they are closer than DELTA_ANGLE_FOR_OVERLAPPING_BONDS
          vector<OBAtom*> nbrs;
          FOR_NBORS_OF_ATOM(a, center)
            nbrs.push_back(&*a);
          double min_angle = 359.0;
          OBBond *close_bond_a = (OBBond*) NULL;
          OBBond *close_bond_b = (OBBond*) NULL;
          for (unsigned int i=0; i<nbrs.size() - 1; ++i)
            for (unsigned int j=i+1; j<nbrs.size(); ++j) {
              double angle = abs(nbrs[i]->GetAngle(center, nbrs[j]));
              if (angle < min_angle) {
                min_angle = angle;
                close_bond_a = mol.GetBond(center, nbrs[i]);
                close_bond_b = mol.GetBond(center, nbrs[j]);
              }
            }

          if (min_angle > DELTA_ANGLE_FOR_OVERLAPPING_BONDS) {
            close_bond_a = (OBBond*) NULL;
            close_bond_b = (OBBond*) NULL;
          }

          // Find the best candidate bond to set to up/down
          // 1. **Should not already be set**
          // 2. Should not be connected to a 2nd tet center
          //    (this is acceptable, as the wedge is only at one end, but will only confuse things)
          // 3. Preferably is not in a cycle
	  // 4. Prefer neighbor with fewer bonds over neighbor with more bonds
          // 5. Preferably is a terminal H, C, or heteroatom (in that order)
          // 6. If two bonds are overlapping, choose one of these
          //    (otherwise the InChI code will mark it as ambiguous)

          int max_bond_score = 0;   // The test below (score > max_bond_score)
          // gave incorrect results when score < 0 and max_bond_score was an unsigned int
          // see https://stackoverflow.com/questions/5416414/signed-unsigned-comparisons#5416498
          FOR_BONDS_OF_ATOM(b, center) {
            if (alreadyset.find(&*b) != alreadyset.end()) continue;

            OBAtom* nbr = b->GetNbrAtom(center);
	    int nbr_nbonds = nbr->GetExplicitDegree();
            int score = 0;
            if (!b->IsInRing()) {
	      if (!nbr->IsInRing())
		score += 8;		// non-ring bond to non-ring atom is good
	      else
		score += 2;		// non-ring bond to ring atom is bad
	    }
            if (tetcenters.find(nbr->GetId()) == tetcenters.end()) // Not a tetcenter
              score += 4;
	    if (nbr_nbonds == 1)	// terminal atom...
		score += 8;		// strongly prefer terminal atoms
	    else
	      score -= nbr_nbonds - 2;	// bond to atom with many bonds is penalized
	    if (nbr->GetAtomicNum() == OBElements::Hydrogen)
	      score += 2;		// prefer H
	    else if (nbr->GetAtomicNum() == OBElements::Carbon)
	      score += 1;		// then C
            if (&*b==close_bond_a || &*b==close_bond_b)
              score += 16;

            if (score > max_bond_score) {
              max_bond_score = score;
              chosen = &*b;
            }
          }

          if (chosen==NULL) { // There is a remote possibility of this but let's worry about 99.9% of cases first
            obErrorLog.ThrowError(__FUNCTION__,
              "Failed to set stereochemistry as unable to find an available bond", obError);
            return false;
          }
          alreadyset.insert(chosen);

          // Determine whether this bond should be set hash or wedge (or indeed unknown)
          // (Code inspired by perception.cpp, TetrahedralFrom2D: plane1 + plane2 + plane3, wedge)
          OBStereo::BondDirection bonddir = OBStereo::UnknownDir;
          if (cfg.winding != OBStereo::UnknownWinding) {
            OBTetrahedralStereo::Config test_cfg = cfg;

            // If there is an implicit ref; let's make that the 'from' atom
            // otherwise use the atom on the chosen bond
            bool implicit = false;
            if (test_cfg.from != OBStereo::ImplicitRef) {
              OBStereo::RefIter ri = std::find(test_cfg.refs.begin(), test_cfg.refs.end(), (unsigned long) OBStereo::ImplicitRef);
              if (ri!=test_cfg.refs.end()) {
                test_cfg = OBTetrahedralStereo::ToConfig(test_cfg, OBStereo::ImplicitRef);
                implicit = true;
              }
            }
            else
              implicit = true;

            bool anticlockwise_order;
            bool useup;
            if (implicit) {
              // Put the ref for the stereo bond second
              while (test_cfg.refs[1] != chosen->GetNbrAtom(center)->GetId())
                std::rotate(test_cfg.refs.begin(), test_cfg.refs.begin() + 2, test_cfg.refs.end());
              anticlockwise_order = AngleOrder(mol.GetAtomById(test_cfg.refs[0])->GetVector(),
                mol.GetAtomById(test_cfg.refs[1])->GetVector(), mol.GetAtomById(test_cfg.refs[2])->GetVector(),
                center->GetVector());
              // Get the angle between the plane bonds
              double angle = GetAngle(mol.GetAtomById(test_cfg.refs[0]), center, mol.GetAtomById(test_cfg.refs[2]));
              if ((angle<0 && anticlockwise_order) || (angle>0 && !anticlockwise_order)) // Is the stereobond in the bigger angle?
                // If the bonds are in anticlockwise order, a clockwise angle (<180) between plane bonds
                // implies that the stereo bond is in the bigger angle. Otherwise it has the opposite meaning.
                useup = anticlockwise_order;
              else
                useup = !anticlockwise_order;
              }
            else {
              test_cfg = OBTetrahedralStereo::ToConfig(test_cfg, chosen->GetNbrAtom(center)->GetId());
              anticlockwise_order = AngleOrder(mol.GetAtomById(test_cfg.refs[0])->GetVector(),
                mol.GetAtomById(test_cfg.refs[1])->GetVector(), mol.GetAtomById(test_cfg.refs[2])->GetVector(),
                center->GetVector());
              if (anticlockwise_order)
                useup = false;
              else
                useup = true;
            }


            // Set to UpBond (filled wedge from cfg.center to chosen_nbr) or DownBond
            bonddir = useup ? OBStereo::UpBond : OBStereo::DownBond;
          }
          updown[chosen] = bonddir;
          from[chosen] = cfg.center;
        }
      }
      return true;
  }

  set<OBBond*> GetUnspecifiedCisTrans(OBMol& mol)
  {
    // Get double bonds with unspecified CisTransStereo
    set<OBBond*> unspec_ctstereo;
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        OBCisTransStereo::Config cfg = ct->GetConfig();
        if (!cfg.specified) {
          OBBond* dbl_bond = mol.GetBond(mol.GetAtomById(cfg.begin), mol.GetAtomById(cfg.end));
          unspec_ctstereo.insert(dbl_bond);
        }
      }
    return unspec_ctstereo;
  }

  void StereoRefToImplicit(OBMol& mol, OBStereo::Ref atomId) {
    // The following is for use in replace_if(...) below
    const std::binder1st<std::equal_to<OBStereo::Ref> > equal_to_atomId = std::bind1st (equal_to<OBStereo::Ref>(), atomId);

    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data) {
      OBStereo::Type datatype = ((OBStereoBase*)*data)->GetType();

      if (datatype != OBStereo::CisTrans && datatype != OBStereo::Tetrahedral) {
        // Maybe I should just unset the stereochemistry if this happens?
        obErrorLog.ThrowError(__FUNCTION__,
            "This function should be updated to handle additional stereo types.\nSome stereochemistry objects may contain explicit refs to hydrogens which have been removed.", obWarning);
        continue;
      }

      // Replace any references to atomId with ImplicitRef
      if (datatype == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        OBCisTransStereo::Config ct_cfg = ct->GetConfig();
        replace_if(ct_cfg.refs.begin(), ct_cfg.refs.end(), equal_to_atomId, (OBStereo::Ref) OBStereo::ImplicitRef);
        ct->SetConfig(ct_cfg);
      }
      else if (datatype == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo::Config ts_cfg = ts->GetConfig();
        if (ts_cfg.from == atomId) ts_cfg.from = OBStereo::ImplicitRef;
        replace_if(ts_cfg.refs.begin(), ts_cfg.refs.end(), equal_to_atomId, (OBStereo::Ref) OBStereo::ImplicitRef);
        ts->SetConfig(ts_cfg);
      }
    }
  }

  void ImplicitRefToStereo(OBMol& mol, OBStereo::Ref centerId, OBStereo::Ref newId) {
    // The following is for use in replace_if(...) below
    const std::binder1st<std::equal_to<OBStereo::Ref> > equal_to_implicitRef = std::bind1st (equal_to<OBStereo::Ref>(), (OBStereo::Ref) OBStereo::ImplicitRef);

    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data) {
      OBStereo::Type datatype = ((OBStereoBase*)*data)->GetType();

      if (datatype != OBStereo::CisTrans && datatype != OBStereo::Tetrahedral) {
        // Maybe I should just unset the stereochemistry if this happens?
        obErrorLog.ThrowError(__FUNCTION__,
            "This function should be updated to handle additional stereo types.\nSome stereochemistry objects may contain implicit refs to hydrogens which need to be converted to explicit.", obWarning);
        continue;
      }

      // Replace any references to ImplicitRef (attached to centerId) with newId
      if (datatype == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        OBCisTransStereo::Config ct_cfg = ct->GetConfig();
        if (ct_cfg.begin == centerId || ct_cfg.end == centerId) {
          // Assumption: the first two refs are on the begin atom, the last two on the end atom
          if (ct_cfg.begin == centerId)
            replace_if(ct_cfg.refs.begin(), ct_cfg.refs.begin()+2, equal_to_implicitRef, (OBStereo::Ref) newId);
          if (ct_cfg.end == centerId)
            replace_if(ct_cfg.refs.begin()+2, ct_cfg.refs.end(), equal_to_implicitRef, (OBStereo::Ref) newId);
          ct->SetConfig(ct_cfg);
        }
      }
      else if (datatype == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo::Config ts_cfg = ts->GetConfig();
        if (ts_cfg.center == centerId) {
          if (ts_cfg.from == OBStereo::ImplicitRef) ts_cfg.from = newId;
          replace_if(ts_cfg.refs.begin(), ts_cfg.refs.end(), equal_to_implicitRef, (OBStereo::Ref) newId);
          ts->SetConfig(ts_cfg);
        }
      }
    }
  }

}
