/**********************************************************************
kekulize.cpp - Alternate algorithm to kekulize a molecule.

Copyright (C) 2004-2006 by Fabien Fontaine
Some portions Copyright (C) 2005-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2009 by Craig A. James

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#define DEBUG 0
// Note about DEBUG: Please don't remove or comment out the "if
// (DEBUG)..."  statements.  When DEBUG is 0, all modern optimizing
// compilers will remove these from the object code completely.

#ifndef MAX_TIME
#define MAX_TIME 15
#endif
#ifndef MAX_DEPTH
#define MAX_DEPTH 30
#endif

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/ring.h>
#include <openbabel/obiter.h>
#include <openbabel/graphsym.h>

#include <sstream>

// Some names to make the code more readable
#define UNASSIGNED 0
#define SINGLE 1
#define DOUBLE 2

#define NOT_IN_RINGS      -1
#define DOUBLE_ASSIGNED    0
#define DOUBLE_ALLOWED     1
#define DOUBLE_PROHIBITED  2

using namespace std;

namespace OpenBabel
{

  //! @cond
  // Hide this from doxygen
  namespace Kekulize {
    // Allow ourselves to restrict recursive searching (we will fall back to LSSR assignment which sometimes makes mistakes)
    struct Timeout
    {
      Timeout(time_t _maxTime) : maxTime(_maxTime)
      {
        startTime = time(NULL);
      }
      time_t startTime, maxTime;
    };
  }
  //! @endcond

  using namespace Kekulize;

  // Modified internal-only version
  // Keep track of which rings contain *all* atoms in the cycle
  // This method essentially does a modified depth-first search to find
  //  large aromatic cycles
  int expand_cycle (OBMol *mol, OBAtom *atom, OBBitVec &avisit, OBBitVec &cvisit,
                    const OBBitVec &potAromBonds, int rootIdx, Timeout &timeout,
                    int prevAtomIdx = -1, int depth = MAX_DEPTH);

  // Modified internal-only version
  // Install a timeout to prevent slow recursion
  bool expandKekulize(OBMol *mol, int bond_idx,
                      std::vector<int> &atomState,
                      std::vector<int> &bondState, Timeout &timeout);

  // LSSR-based version for fallback after timeout
  bool expand_kekulize_lssr(OBMol *mol,
                            std::vector<int> &atomState,
                            std::vector<int> &bondState,
                            std::vector<OBRing*> &lssr,
                            std::vector<bool> &lssrAssigned,
                            std::vector<OBBond *> &bondsThisRing);

  // Count the number of bonds in one ring that are already assigned
  // (single or double) during Kekule assignment for a ring system.  For
  // example, in naphthalene, if we've already assigned single/double to
  // the first ring, then calling this function for the second ring will
  // return 1: one bond is already assigned.
  int count_assigned_bonds(std::vector<OBBond *> &bondsThisRing,
                           std::vector<int> &bondState)
  {
    int n = 0;
    for (unsigned int i = 0; i < bondsThisRing.size(); i++) {
      if (bondState[bondsThisRing[i]->GetIdx()] != UNASSIGNED) {
        n++;
      }
    }
    return n;
  }

  // Determine if this atom is possibly aromatic (to cut down on the search space
  // Only considers C, N, O, Si, P, S, As, Se, Te
  bool isPotentialAromaticAtom(OBAtom *atom)
  {
    switch (atom->GetAtomicNum()) {
    case 6:
      return (atom->GetHvyValence() < 4);
    case 7:
    case 8:
    case 14: // Si
    case 15: // P
    case 16: // S
    case 33: // As
    case 34: // Se
    case 52: // Te
      return true;
    default:
      return false;
    }
  }

  // Check the molecule for potential aromatic bonds
  // (e.g., a bond to a non-aromatic atom doesn't count)
  // Minimizes the search space for Kekule assignments
  void potentialAromaticBonds(OBMol *mol, OBBitVec &bonds)
  {
    std::vector<OBRing*> rings = mol->GetLSSR();

    for (std::size_t i = 0; i < rings.size(); ++i) {
      bool skip = false;
      bool allCarbon = true;
      // check if all atoms are potential aromatic atoms
      for (std::size_t j = 0; j < rings[i]->_path.size(); ++j) {
        OBAtom *atom = mol->GetAtom(rings[i]->_path[j]);
        if (!atom || !isPotentialAromaticAtom(atom)) { //!atom is temp crash prevention
          skip = true;
          break;
        }
        if (!atom->IsCarbon())
          allCarbon = false;
      }

      if (skip)
        continue;
      // a 4 membered ring containing only carbon can't be aromatic
      if (rings[i]->_path.size() == 4 && allCarbon)
        continue;

      // set the bond bits
      for (std::size_t j = 1; j < rings[i]->_path.size(); ++j)
        bonds.SetBitOn(mol->GetBond(rings[i]->_path[j-1], rings[i]->_path[j])->GetIdx());
      bonds.SetBitOn(mol->GetBond(rings[i]->_path[rings[i]->_path.size()-1], rings[i]->_path[0])->GetIdx());
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  //! Returns a vector of the bonds in a ring, in circular order around the ring,
  //
  // The circular order (the fact that the bonds are end-to-end in the
  // vector) is important when assigning single/double bonds to an aromatic
  // system, because it helps you "walk" around a ring assigning
  // alternating single/double bonds, and you can often get it right on the
  // first try.
  void get_bonds_of_ring(OBMol *mol, OBRing *ring, std::vector<OBBond *> &ring_bonds) {

    vector<bool> used_bonds;
    used_bonds.resize(mol->NumBonds());
    for (unsigned int i = 0; i < used_bonds.size(); i++)
      used_bonds[i] = false;

    ring_bonds.clear();
    int start_idx = ring->_path[0];
    OBAtom *start_atom = mol->GetAtom(start_idx);
    OBAtom *atom2 = start_atom;
    OBAtom *atom;
    do {
      atom = atom2;				// move to the next atom in the ring
      FOR_BONDS_OF_ATOM(b, &(*atom)) {		// loop over its bonds to find the next atom...
        int bond_idx = b->GetIdx();
        if (used_bonds[bond_idx])		// is this the previous neighbor atom?
          continue;				//    skip it ... already done.
        atom2 = b->GetNbrAtom(atom);
        if (!ring->IsMember(atom2))		// is this neighbor in the ring?
          continue;				//    not in aromatic system, ignore it
        ring_bonds.push_back(&(*b));		// found the bond we want
        used_bonds[bond_idx] = true;
        break;
      }
    } while (atom2 && atom2 != start_atom);
  }

  ///////////////////////////////////////////////////////////////////////////////
  //! \brief Kekulize aromatic rings without using implicit valence
  //!
  //! This new perceive kekule bonds function has been designed to
  //! handle molecule files without explicit hydrogens such as pdb or xyz.
  //! (It can, of course, easily handle explicit hydrogens too.)
  //! The function does not rely on GetImplicitValence function
  //! The function looks for groups of aromatic cycle
  //! For each group it tries to guess the number of electrons given by each atom
  //! in order to satisfy the huckel (4n+2) rule
  //! If the huckel rule cannot be satisfied the algorithm try with its best alternative guess
  //! Then it recursively walk on the atoms of the cycle and assign single and double bonds
  //! \deprecated Will no longer be a public, visible method
  void OBMol::NewPerceiveKekuleBonds()
  {

    if (HasKekulePerceived())  return;
    SetKekulePerceived();

    OBAtom *atom;
    int n, de, minde;
    std::vector<OBAtom*> cycle;
    OBBitVec avisit,cvisit;
    avisit.Resize(NumAtoms()+1);
    cvisit.Resize(NumAtoms()+1);
    OBBond *bond;
    std::vector<OBBond*>::iterator bi;
    std::vector<int> electron;
    int sume, orden;
    int bestatom = 1;
    int bestorden = 99;
    // Init the kekulized bonds
    unsigned i;
    FOR_BONDS_OF_MOL(b, *this)
      {
        switch (b->GetBO())
          {
          case 2: b->SetKDouble(); break;
          case 3: b->SetKTriple(); break;
          case 1: b->SetKSingle(); break;
          }
      }

    // Check if we should allow "fused ring" analysis -- we fall back to expanding as needed
    // We'll only do this if an atom is in 3 or more rings.
    bool fusedRings = false;
    OBSmartsPattern fused; fused.Init("[aR3]");
    if (fused.Match(*this))
      fusedRings = true;

    OBBitVec potAromBonds;
    potentialAromaticBonds(this, potAromBonds);

    // Find all the groups of aromatic cycle
    for(i=1; i<= NumAtoms(); i++ ) {
      atom = GetAtom(i);
      if (!isPotentialAromaticAtom(atom))
        continue;
      //      cout << "Checking for cycle at " << i << endl;
      if (atom->HasAromaticBond() && !cvisit[i]) { // is new aromatic atom of an aromatic cycle ?

        avisit.Clear();
        electron.clear();
        cycle.clear();

        avisit.SetBitOn(i);

        if (fusedRings)
          expandcycle(atom, avisit, potAromBonds);
        else {
          Timeout timeout(MAX_TIME);
          int depth = expand_cycle(this, atom, avisit, cvisit, potAromBonds, atom->GetIdx(), timeout);
          if (depth <= 0)
            continue; // no valid cycle from this atom
        }

        // Check to see that at least 2 bonds are included
        // e.g.           __
        //           ____/  \   bad because, we should start at a ring atom
        //               \__/
        std::vector<OBBond*>::iterator b;
        OBAtom *nbr;

        unsigned int bondCount = 0;
        for (nbr = atom->BeginNbrAtom(b);nbr;nbr = atom->NextNbrAtom(b)) {
          if (avisit[nbr->GetIdx()])
            bondCount++;
        }
        if (bondCount < 2) {
          //          cout << "rejected cycle" << endl;
          continue; // this atom isn't a root for the cycle
        }

        //        cout << " depth: " << depth << endl;
        //store the atoms of the cycle(s)
        unsigned int j;
        if (DEBUG)  cout << " cycle: " << NumAtoms() << " ";
        for(j=1; j<= NumAtoms(); ++j) {
          if ( avisit[j] ) {
            if (DEBUG)  cout << "\t" << j;
            atom = GetAtom(j);
            cycle.push_back(atom);
          }
        }
        if (DEBUG) cout << endl;

        // This isn't a real aromatic cycle -- give up
        // Fixes PR#1965566
        // Fixes PR#1784204
        if (cycle.size() < 3)
          continue;

        // At the beginning each atom give one electron to the cycle
        for(j=0; j< cycle.size(); ++j) {
          atom = cycle[j];
          if (atom->IsOxygen()
              && atom->GetFormalCharge() == 0
              && atom->GetValence() == 2) // fixes PR#3075065
            electron.push_back(2);
          else
            electron.push_back(1);
        }

        // remove one electron if the atom makes a double bond out of the cycle
        sume = 0;
        bool doubleBondInRing = false;
        for(j=0; j< cycle.size(); ++j) {
          atom = cycle[j];
          for(bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
            if ( bond->IsDouble()) {
              // monitor this, because these count for electrons, even if this atom doesn't take a double bond
              if (bond->IsInRing())
                doubleBondInRing = true;

              OBAtom *atom2 = bond->GetNbrAtom(atom);
              int fcharge = atom->GetFormalCharge();
              int fcharge2 = atom2->GetFormalCharge();
              if(atom->IsNitrogen() && atom2->IsOxygen()
                 && fcharge == 0 && fcharge2 == 0) { //n=O to [n+][O-]
                atom->SetFormalCharge(1);
                atom2->SetFormalCharge(-1);
                bond->SetKSingle();
                bond->SetBO(1);
              }
              else {
                electron[j] = 0;
              }
            }
          }
          // count the number of electrons
          if (!doubleBondInRing) {
            sume += electron[j];
          } else {
            sume++; // one electron contribution for explicit double bonds into another ring
          }

          if (atom->IsNitrogen() && atom->GetFormalCharge() == 0 && atom->GetValence() >= 3) {
            electron[j] = 0; // no double bonds, first try to kekulize with 1 electron, add two later if it doesn't work
          }
        }

        // Save the electron state in case huckel rule is not satisfied
        vector<int> previousElectron = electron;

        // find the ideal number of electrons according to the huckel 4n+2 rule
        minde=99;
        for (int x=1; 1; ++x) { // do NOT use "i": main for() loop counter
          n = 4 *x +2;
          de = n - sume;
          if (  de < minde )
            minde=de;
          else if ( minde < 0 )
            minde=de;
          else
            break;
        }

        stringstream errorMsg;

        if (DEBUG)  cout << "sume " << sume << " minde before:" << minde << endl;
        // if huckel rule not satisfied some atoms must give more electrons
        int currentState;
        while ( minde != 0 ) {
          bestorden=99;
          for(j=0; j< cycle.size(); ++j) {
            currentState = electron[j];
            if (cycle[j]->IsNitrogen() && cycle[j]->GetFormalCharge() == 0 && cycle[j]->GetValence() >= 3)
              currentState = 1; // might need to add another electron

            if (currentState == 1) {
              orden = getorden(cycle[j]);
              if (orden < bestorden) {
                bestorden = orden;
                bestatom = j;
                if (DEBUG) {cout << "atom " << cycle[j]->GetIdx() << ", bestorden =" << orden << endl;}
              }
            }
          }
          if (bestorden==99) {  // no electron giving atom found
            errorMsg << "Kekulize: Huckel rule not satisfied for molecule " << GetTitle() << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
            break;             // Huckel rule cannot be satisfied
          }                    // try to kekulize anyway
          else {
            electron[bestatom] += 1;
            minde--;
            if (DEBUG) {cout << "electron[atom:" << cycle[bestatom]->GetIdx() << "] = " << electron[bestatom] << endl;}
          }
        }

        if (bestorden == 99) { // Huckel rule not satisfied, just try to get an even number of electron before kekulizing

          electron = previousElectron; // restore electon's state

          int odd = sume % 2;
          //cout << "odd:" << odd << endl;
          if(odd) { // odd number of electrons try to add an electron to the best possible atom
            for(j=0; j< cycle.size(); ++j) {

              currentState = electron[j];
              if (cycle[j]->IsNitrogen() && cycle[j]->GetFormalCharge() == 0 && cycle[j]->GetValence() >= 3)
                currentState = 1; // might need to add another electron

              if (currentState == 1) {
                orden = getorden(cycle[j]);
                if (orden < bestorden) {
                  bestorden = orden;
                  bestatom = j;
                }
              }
            }
            if (bestorden==99) {  // no electron giving atom found
              errorMsg << "Kekulize: Cannot get an even number of electron for molecule " << GetTitle() << "\n";
              obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
              break;             // impossible to choose an atom to obtain an even number of electron
            }                    // try to kekulize anyway
            else {
              electron[bestatom] += 1;
            }
          }
        }

        if (DEBUG) {
          cout << "sume: " << sume << " minde after:" << minde <<endl;
          cout << "atom: ";
          for(j=0; j < cycle.size(); ++j) {
            OBAtom *cycleAtom = cycle[j];
            cout << "\t" << cycleAtom->GetIdx();
          }
          cout << endl;

          cout << "atom: ";
          for(j=0; j < electron.size(); ++j) {
            cout << "\t" << electron[j];
        	}
        	cout << endl;
        }

        // kekulize the cycle(s)
        start_kekulize(cycle,electron);

        // Set the kekulized cycle(s) as visited
        for(j=1; j<= NumAtoms(); ++j) {
          if (avisit[j])
            cvisit.SetBitOn(j);
        }

      }
    }
    // Double bond have been assigned, set the remaining aromatic bonds to single
    //std::cout << "Set not assigned single bonds\n";
    FOR_BONDS_OF_MOL(b, *this)
      {
        //std::cout << "bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << " ";
        if (b->GetBO()==5 ) {
          b->SetKSingle();
          b->SetBO(1);
          //std::cout << "single\n";
        }
        //else
        //  std::cout << "double\n";
      }

    return;
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //! \brief Find a consistent assignment of single/double bonds to a
  //! Kekule' ring or a set of fused aromatic rings.
  //!
  //! The initial electronic state indicates whether each atom can make a
  //! double bond or not.  The function works recursively to "walk" around
  //! the ring or rings and try all possible arrangements of single and
  //! double bonds.
  //! \deprecated Will no longer be a public, visible method
  void OBMol::start_kekulize( std::vector <OBAtom*> &cycle, std::vector<int> &electron) {

    // initialize atom and bond states in constructor
    std::vector<int> atomState(NumAtoms() + 1, NOT_IN_RINGS); // atom index issue
    std::vector<int> bondState(NumBonds(), SINGLE);
    int idx;
    OBAtom *atom;
    OBBond *bond;

    // Figure out which atoms are in this ring system and whether or not each
    // atom can donate an electron.
    for (unsigned int i = 0; i < cycle.size(); ++i) {
      atom = cycle[i];
      idx =  atom->GetIdx();
      if (electron[i] == 1) {
        atomState[idx] = DOUBLE_ALLOWED;	// It has an electron it can donate
      } else {
        atomState[idx] = DOUBLE_PROHIBITED;	// No electrons to contribute to aromatic system
      }

      if (atomState[idx] == DOUBLE_ALLOWED && atom->IsNitrogen() && atom->GetFormalCharge() == 0 && atom->GetValence() == 3) {
        // Correct N with three explicit bonds, if we haven't already
        atomState[idx] = DOUBLE_PROHIBITED;
        if (DEBUG) { cout << "atom " << idx << " rejected NR3 double bonds " << endl; }
      }
      if (DEBUG) {cout << "atom " << idx << ": initial state = " << atomState[idx] << endl;}
    }

    // Do a recursive walk around the aromatic system, and see if we can find
    // an assignment of double/single bonds that works.
    atom = cycle[0];
    Timeout timeout(MAX_TIME); // install a timeout to fallback to LSSR analysis
    if (!expandKekulize(this, 0, atomState, bondState, timeout) ) {
      // OK, now try using LSSR analysis

      // Get the LSSR, and find all rings that are completely in cycle[], that
      // is, every atom in the ring is in this aromatic ring system.
      std::vector<OBRing*> lssr_tmp = GetLSSR();
      std::vector<OBRing*> lssr;
      std::vector<bool>atom_in_cycle(NumAtoms()+1);
      for (unsigned int i = 0; i <= NumAtoms(); i++)
        atom_in_cycle[i] = false;
      for (std::vector<OBAtom*>::iterator ai = cycle.begin(); ai != cycle.end(); ai++)
        atom_in_cycle[(*ai)->GetIdx()] = true;
      for (std::vector<OBRing*>::iterator ri = lssr_tmp.begin(); ri != lssr_tmp.end(); ri++) {
        OBRing *ring = *ri;
        bool ok = true;
        for (std::vector<int>::iterator pi = ring->_path.begin(); pi != ring->_path.end(); pi++) {
          if (!atom_in_cycle[*pi]) {
            ok = false;
            break;
          }
        }
        if (ok)
          lssr.push_back(ring);
        if (DEBUG) {cout << "ring " << (ri - lssr_tmp.begin()) << (ok ? " is" : " is not") << " in cycle" << endl;}
      }

      // Initialize the "ring is assigned" vector
      std::vector<bool> lssrAssigned;
      lssrAssigned.resize(lssr.size());
      for (unsigned int i = 0; i < lssr.size(); i++)
        lssrAssigned[i] = false;

      // Initialize an empty vector that expand_kekulize2() uses.
      std::vector<OBBond*> bondsThisRing;
      bondsThisRing.clear();

      // Do a recursive walk across the aromatic system, and see if we can find
      // an assignment of double/single bonds that works.
      if (!expand_kekulize_lssr(this, atomState, bondState, lssr, lssrAssigned, bondsThisRing)) {
        stringstream errorMsg;
        errorMsg << "Kekulize Error for molecule " << GetTitle() << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        return;
    }
      // If we did successfully use the LSSR assignment, we don't need to worry -- no error!
    }

    // We found a successful assignment, now actually change the bonds to double

    if (DEBUG) {std::cout << "Set double bonds\n";}
    for (unsigned int i = 0; i < NumBonds(); ++i) {
      bond = GetBond(i);
      if (DEBUG) { std::cout << "  bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << " ";}
      if (bond->GetBO()==5 && bondState[i] == DOUBLE) {
        if (   bond->GetBeginAtom()->IsSulfur()
               && bond->GetEndAtom()->IsSulfur()) {
          // no double bonds between aromatic sulfur atoms -- PR#1504089
          continue;
        }

        bond->SetKDouble();
        bond->SetBO(2);
        if (DEBUG) {std::cout << "double\n";}
      }
      else {
        if (DEBUG) {std::cout << "single\n";}
      }
    }

    return;
  }


  /////////////////////////////////////////////////////////////////////////////////////////
  //! \brief Recursively assign single and double bonds according to the electronical state
  //! of the atoms.
  //
  // This can be thought of as a large boolean equation, with each bond of
  // the molecule a "boolean variable".  The bond can be double or single,
  // sort of like true/false, and the "boolean equation" is whether the
  // particular combination of single/double values results in a correct
  // electronic configuration.  This algorithm is pretty much a brute-force
  // test of every combination of single/double that is electronically
  // reasonable.
  //
  // The recursion is called with bond N.  It tries to make bond N double,
  // if that works, it recurses to bond N+1, and if that works, we're done.
  // Otherwise, it reverts to a single bond at N, and recurses to bond N+1.
  // If that works, we return success, if it fails, return failure.
  //
  // This recursive algorithm guarantees that if it succeeds (finds bond
  // assignments and returns true), then atomState[] will have the double
  // bonds marked, but if it fails, (returns false), then the atomState[]
  // and bondState[] will be unaltered.
  //
  // Note that the time to run this test grows as O(2^N), where N is the
  // number of double bonds that are needed to make the Kekule bonds.  It
  // will be very slow for large systems such as fullerenes.  The way to
  // make it more efficient would be to include a symmetry analysis and
  // divide the bonds into symmetrically-equivalent classes, and only try
  // combinations of single/double bonds that differ when symmetry is taken
  // into account.
  //! \deprecated Will no longer be a public, visible method
  bool OBMol::expand_kekulize(int bond_idx,
                              std::vector<int> &atomState,
                              std::vector<int> &bondState)
  {
    Timeout timeout(MAX_TIME);
    return expandKekulize(this, bond_idx, atomState, bondState, timeout);
  }

  // Used in old expandKekulize code (i.e., fused-ring analysis)
  bool has_leftover_electrons(OBMol *mol, std::vector<int> &atomState)
  {
    FOR_ATOMS_OF_MOL(a, mol) {
      int idx = a->GetIdx();
      if (atomState[idx] == DOUBLE_ALLOWED) { // we haven't assigned this atom yet
        if (DEBUG) {cout << "  failure, extra electron on atom " << idx << endl;}
        return true;
      }
    }
    return false;	// no extra electrons found
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Check for leftover electrons.  This is used during expand_kekulize_lssr()
  // to make sure all of the 4n+2 electrons that are available for bonding in the
  // aromatic ring system are actually used.
  //
  // Scans the aromatic atoms, and checks each one to see if all of its bonds
  // have been assigned.  If they have, then all of the atom's electrons have to
  // be used up, otherwise it's an invalid bond assignment.
  //
  // This is called repeatedly as the bonds are assigned, to help "fail fast."
  // As the bonds are assigned to a complex ring system, there's no point waiting
  // until the end to discover that you made an invalid assignment early on.
  bool has_leftover_electrons(OBMol *mol, std::vector<int> &atomState,
                              std::vector<int> &bondState)
  {
    FOR_ATOMS_OF_MOL(a, mol) {
      int idx = a->GetIdx();
      bool atom_ok = false;
      if (atomState[idx] == DOUBLE_ALLOWED) {
        FOR_BONDS_OF_ATOM(b, &(*a)) {
          OBAtom *neighbor = b->GetNbrAtom(&(*a));
          if (atomState[neighbor->GetIdx()] == NOT_IN_RINGS) {
            continue;	// don't care about non-ring bonds
          }
          if (bondState[b->GetIdx()] == UNASSIGNED) {
            atom_ok = true;	// if one in-ring bond is unassigned, then
            break;		// it's OK if the atom still has an extra electron
          }
        }
        // If no unassigned in-ring bonds were found, and the atom
        // has an extra electron, it can't be valid.
        if (!atom_ok) {
          if (DEBUG) {cout << "  failure, extra electron on atom " << idx << endl;}
          return true;
        }
      }
    }
    return false;	// no extra electrons found
  }

  //! Check for leftover electrons.  This is used during expand_kekulize()
  //! to make sure all of the 4n+2 electrons that were available for bonding in the
  //! aromatic ring system were actually used during the assignment of single
  //! and double bonds.
  //! \deprecated Will no longer be a public, visible method
  bool OBMol::has_no_leftover_electrons(std::vector<int> &atomState)
  {
    return !has_leftover_electrons(this, atomState);
  }

  bool expandKekulize(OBMol *mol, int bond_idx,
                      std::vector<int> &atomState,
                      std::vector<int> &bondState, Timeout &timeout)
  {
    // Check for timeout first
    if (time(NULL) - timeout.startTime > timeout.maxTime) {
      // Only issue the error once, not for each recurison
      if (bond_idx == 0) obErrorLog.ThrowError(__FUNCTION__, "maximum time exceeded...", obError);
      return !has_leftover_electrons(mol, atomState);
    }

    // If all bonds are assigned, check that this is a sensible combination.
    // This is the "end of the line" for the recursion: Did this bond assignment
    // work or not?
    if (bond_idx >= bondState.size())
      return !has_leftover_electrons(mol, atomState);

    // Get the bond and its atoms
    OBBond *bond = mol->GetBond(bond_idx);
    OBAtom *atom1 = bond->GetBeginAtom();
    OBAtom *atom2 = bond->GetEndAtom();
    int idx1 = atom1->GetIdx();
    int idx2 = atom2->GetIdx();

    // If this bond isn't part of the aromatic system, move on to the next bond.
    if (atomState[idx1] == NOT_IN_RINGS || atomState[idx2] == NOT_IN_RINGS)
      return expandKekulize(mol, bond_idx + 1, atomState, bondState, timeout);

    // Remember the current state so that we can backtrack if the attempt fails
    // Keep previous states on heap to avoid stack overflows for large molecules
    vector<int>* ppreviousState     = new vector<int>(atomState);  // Backup the atom states
    vector<int>* ppreviousBondState = new vector<int>(bondState);  // ... and the bond states

    // Is a double bond allowed here?  Both atoms have to have an extra electron,
    // and not have been assigned a double bond from a previous step in recursion.
    if (atomState[idx1] == DOUBLE_ALLOWED && atomState[idx2] == DOUBLE_ALLOWED) {

      // Assign a double bond.
      atomState[idx1]  = DOUBLE_ASSIGNED;
      atomState[idx2]  = DOUBLE_ASSIGNED;
      bondState[bond_idx] = DOUBLE;        // set bond to double
      if (DEBUG) {std::cout << "bond " << bond_idx << " (atoms " << idx1 << " to " << idx2 << ") double\n";}

      // Recursively try the next bond
      if (expandKekulize(mol, bond_idx + 1, atomState, bondState, timeout)) {
        delete ppreviousState;
        delete ppreviousBondState;
        return true;
      }

      // If the double bond didn't work, roll back the changes and try a single bond.
      atomState = *ppreviousState;
      bondState = *ppreviousBondState;
      if (DEBUG) {cout << "  double on bond " << bond_idx << " failed." << endl;}
    }

    // Double bond not allowed here, or double bond failed, just recurse with a single bond.
    if (DEBUG) {cout << "bond " << bond_idx << " (atoms " << idx1 << " to " << idx2 << ") single" << endl;}
    if (expandKekulize(mol, bond_idx + 1, atomState, bondState, timeout)) {
      delete ppreviousState;
      delete ppreviousBondState;
      return true;
    }

    // If it didn't work, roll back the changes we made and return failure.
    if (DEBUG) {cout << "bond " << bond_idx << " single failed, rolling back changes" << endl;}
    atomState = *ppreviousState;
    bondState = *ppreviousBondState;
    delete ppreviousState;
    delete ppreviousBondState;
    return false;
  }

  //! Give the priority to give two electrons instead of 1
  //! \deprecated Will no longer be a public, visible method
  int OBMol::getorden( OBAtom *atom)
  {
    if ( atom->IsSulfur() && atom->GetFormalCharge() == 0) return 0;
    if ( atom->IsOxygen() ) return 2;
    if ( atom->GetAtomicNum() == 34 || atom->GetAtomicNum() == 52 ) return 3;
    if ( atom->IsNitrogen() && atom->GetFormalCharge() == 0 && atom->GetValence() == 3) return 4;
    if ( atom->IsAmideNitrogen() ) return 5;
    if ( atom->IsNitrogen() && atom->GetFormalCharge() == -1) return 6;
    if ( atom->IsNitrogen() && atom->GetFormalCharge() == 0 && atom->IsInRingSize(5) ) return 7;
    if ( atom->IsNitrogen() && atom->GetFormalCharge() == 0 ) return 8;
    if ( atom->IsCarbon() && atom->GetFormalCharge() == -1) return 9;
    //if ( atom->IsCarbon() ) return 9;

    return (100); //no atom found
  }

  //! Recursively find the aromatic atoms with an aromatic bond to the current atom
  //! \deprecated Will no longer be a public, visible method
  bool OBMol::expandcycle (OBAtom *atom, OBBitVec &avisit, const OBBitVec &potAromBonds)
  {
    OBAtom *nbr;
    //  OBBond *bond;
    std::vector<OBBond*>::iterator i;
    int natom;
    //for each neighbour atom test if it is in the aromatic ring
    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
      {
        if (!potAromBonds.BitIsSet((*i)->GetIdx()))
          continue;
        natom = nbr->GetIdx();
        if (!avisit[natom] && ((OBBond*) *i)->GetBO()==5
            && ((OBBond*) *i)->IsInRing()) {
          avisit.SetBitOn(natom);
          expandcycle(nbr, avisit, potAromBonds);
        }
      }

    return true;
  }

  //! Recursively find the aromatic atoms with an aromatic bond to the current atom
  int expand_cycle (OBMol *mol, OBAtom *atom, OBBitVec &avisit, OBBitVec &cvisit,
                    const OBBitVec &potAromBonds, int rootIdx, Timeout &timeout,
                    int prevAtomIdx, int depth)
  {
    // early termination -- too deep recursion, or timeout
    if (depth < 0)
      return depth;

    if (time(NULL) - timeout.startTime > timeout.maxTime) {
      obErrorLog.ThrowError(__FUNCTION__, "maximum time exceeded...", obError);
      return depth;
    }

    if (DEBUG) cout << " expand_cycle: " << atom->GetIdx() << " depth " << depth << endl;

    OBAtom *nbr;
    std::vector<OBBond*>::iterator i;
    int natom;

    // OK, here's the plan:
    // - If a neighboring atom is non-aromatic, we ignore it
    // - If the atom is in cvisit, it's already be assigned. Ignore it.
    // - If the atom is the previous step in the path, ignore it.
    // Otherwise recurse: look for a large cycle back to the root
    int trialScore, bestScore = 1000;
    OBBitVec trialMatch, bestMatch; // the best path we've found so far
    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
      {
        if (!potAromBonds.BitIsSet((*i)->GetIdx()))
          continue;
        natom = nbr->GetIdx();
        if ((*i)->GetBO() != 5)
          continue; // this is a non-aromatic bond, skip it
        if (natom == prevAtomIdx) {
          // the previous step in our path
          continue;
        }

        if (avisit[natom] && natom != rootIdx) {
          continue; // skip this path, we should try to get to the root again
        }

        if (natom == rootIdx) {
          bestMatch = avisit;
          bestMatch.SetBitOn(natom);
          bestScore = depth;
          continue; // don't recurse further
        }

        trialMatch = avisit;
        trialMatch.SetBitOn(natom);
        trialScore = expand_cycle(mol, nbr, trialMatch, cvisit, potAromBonds, rootIdx, timeout, atom->GetIdx(), depth - 1);
        if (trialScore > 0 && trialScore < bestScore) { // we found a larger, valid cycle
          if (DEBUG) cout << " score: " << trialScore << endl;
          bestMatch = trialMatch;
          bestScore = trialScore;
        }
      } // check all neighbors

    if (bestScore <= 0 || bestScore == 1000)
      return -1; // e.g., we have no valid choices

    avisit = bestMatch; // remember our path
    return bestScore;
  }

  /// NEW ANALYSIS methods: rely on LSSR methods to only target one ring. Works much faster on fullerens, graphene, etc.

  /////////////////////////////////////////////////////////////////////////////////////////
  //! \brief Recursive function to find a sensible kekule assignment of
  //! single and double bonds for an aromatic ring system.
  //
  // Very large aromatic ring systems can be very computationally difficult.  In theory,
  // there are 2^N ways to try assigning single/double to assign bonds (where N is the
  // number of bonds in the aromatic system).  For large molecules such as fullerenes,
  // this number (e.g. 2^60) exceeds all computing power in the world.
  //
  // To avoid this, three strategies are used.
  //
  // 1. The algorithm proceeds by assigning bonds to each complete ring
  //    from the LSSR, rather than treating the whole aromatic system.
  //    That reduces the complexity to roughly O(2^R) where R is the number
  //    of rings in the LSSR.
  //
  // 2. Each time a ring is completed (a trial assignment of
  //    single/double), the whole tentative assignment is checked for
  //    sensibility.  If any atom has all of its bonds assigned, but
  //    still has a leftover electron, then the assignment fails before a
  //    lot more work is invested.  (Atoms connected to bonds that haven't
  //    been examined yet are ignored during this test.)
  //
  // 3. When one ring is finished, the next ring is selected by finding the
  // 	one "most connected" to previously-considered rings.  Say, for
  // 	example, we're working on a fullerene and have just finished one
  // 	ring.  The second ring will be one that's adjacent to the first
  // 	ring (one of its bonds is already assigned).  Now for the third
  // 	ring, we can pick one that's in a "corner" of the first two, that
  // 	is, a ring that has TWO bonds already assigned.
  //
  // Strategy #3 means that the completed area of the ring system (the
  // bonds we've already assigned) will tend to have a "minimal perimeter",
  // and any atom that has one bond assigned will quickly have all of its
  // bonds assigned.  That means that strategy #2 can be used as early as
  // possible to detect bond assignments that can never work, long before
  // the entire aromatic system has been examined.  Strategy #1 means that
  // we can typically assign alternating single/double bonds correctly to
  // any particular ring on the first try.
  //
  // Together, these three strategies together mean that even very large
  // aromatic systems are typically solved in milliseconds.
  bool expand_kekulize_lssr(OBMol *mol,
                            std::vector<int> &atomState,
                            std::vector<int> &bondState,
                            std::vector<OBRing*> &lssr,
                            std::vector<bool> &lssrAssigned,
                            std::vector<OBBond *> &bondsThisRing)
  {
    OBBond *bond = NULL;
    OBAtom *atom1, *atom2;
    int idx1, idx2, bond_idx;

    if (DEBUG) {cout << "---------- expand_kekulize_lssr:" << endl;}

    // Find the next unassigned bond on this ring
    for (std::vector<OBBond*>::iterator b = bondsThisRing.begin(); b != bondsThisRing.end(); b++) {
      if (bondState[(*b)->GetIdx()] == UNASSIGNED) {
        bond = *b;
        break;
      }
    }

    if (bond) {

      bond_idx = bond->GetIdx();
      atom1 = bond->GetBeginAtom();
      atom2 = bond->GetEndAtom();
      idx1 = atom1->GetIdx();
      idx2 = atom2->GetIdx();

      // Remember the current state so that we can backtrack if the attempt fails
      // Keep previous states on heap to avoid stack overflows for large molecules
      vector<int>* ppreviousState     = new vector<int>(atomState);  // Backup the atom states
      vector<int>* ppreviousBondState = new vector<int>(bondState);  // ... and the bond states

      // Does a double bond work here?
      if (   atomState[idx1] == DOUBLE_ALLOWED
             && atomState[idx2] == DOUBLE_ALLOWED) {

        // Assign a double bond.
        atomState[idx1]  = DOUBLE_ASSIGNED;
        atomState[idx2]  = DOUBLE_ASSIGNED;
        bondState[bond_idx] = DOUBLE;        // set bond to double
        if (DEBUG) {std::cout << "bond " << bond_idx << " (atoms " << idx1 << " to " << idx2 << ") double" << endl;}

        // Recursively try the next bond
        if (expand_kekulize_lssr(mol, atomState, bondState, lssr, lssrAssigned, bondsThisRing)) {
          delete ppreviousState;
          delete ppreviousBondState;
          return true;
        }

        // If the double bond didn't work, roll back the changes and try a single bond.
        atomState = *ppreviousState;
        bondState = *ppreviousBondState;
        if (DEBUG) {cout << "  double on bond " << bond_idx << " failed." << endl;}
      }

      // Either a double bond not allowed here, or double bond failed, try a single bond.
      bondState[bond_idx] = SINGLE;        // set bond to single
      if (DEBUG) {cout << "bond " << bond_idx << " (atoms " << idx1 << " to " << idx2 << ") single" << endl;}

      // Recursively try the next bond
      if (expand_kekulize_lssr(mol, atomState, bondState, lssr, lssrAssigned, bondsThisRing)) {
        delete ppreviousState;
        delete ppreviousBondState;
        return true;
      }

      // If it didn't work, roll back the changes we made and return failure.
      if (DEBUG) {cout << "bond " << bond_idx << " single failed, rolling back changes" << endl;}
      atomState = *ppreviousState;
      bondState = *ppreviousBondState;
      delete ppreviousState;
      delete ppreviousBondState;
      return false;
    }

    // If we get here, it means we successfully assigned all of the bonds on
    // the current ring.  Double check that we haven't made an impossible
    // bond assignment.

    if (has_leftover_electrons(mol, atomState, bondState))
      return false;

    // Now its time to pick the next ring to work on.
    //
    // Find the ring that has the most bonds already assigned.  This will
    // tend to "localize" the search, so that the algorithm doesn't wander.
    // For example, imagine a large sheet of graphite ("chickenwire").
    // This algorithm will make it fill in a local area, rather than
    // stretching across the sheet.  You can think of it like minimizing
    // the perimeter of the area of assigned bonds.
    //
    // This is critical to the performance of this algorithm, because it
    // "fails early, fails fast."

    int best_rnum = -1;
    int best_nbonds_assigned = -1;
    std::vector<OBBond *> best_bondsThisRing;

    if (DEBUG) {cout << "Select next ring: " << lssr.size() << " rings to try" << endl;}

    for (unsigned int rnum = 0; rnum < lssr.size(); rnum++) {
      if (lssrAssigned[rnum])
        continue;
      OBRing *ring = lssr[rnum];
      get_bonds_of_ring(mol, ring, bondsThisRing);
      int nbonds_assigned = count_assigned_bonds(bondsThisRing, bondState);
      if (DEBUG) {cout << "Trying ring " << rnum << ": " << nbonds_assigned << " bonds assigned" << endl;}
      if (nbonds_assigned > best_nbonds_assigned) {
        best_rnum = rnum;
        best_nbonds_assigned = nbonds_assigned;
        best_bondsThisRing = bondsThisRing;
        if (DEBUG) {cout << "    new best ring: " << rnum << "\n";}
      }
    }

    // If we didn't find any unassigned ring at all, then we're done,
    // the Kekule assignment was a success.
    if (best_rnum == -1) {
      if (DEBUG) {cout << "    No more rings, single/double bond assignment succeeded." << endl;}
      return true;
    }

    // We found a new ring to assign.  Recursive call to assign its bonds.
    lssrAssigned[best_rnum] = true;
    if (!expand_kekulize_lssr(mol, atomState, bondState, lssr, lssrAssigned, best_bondsThisRing)) {
      lssrAssigned[best_rnum] = false;
      return false;
    } else {
      return true;
    }
  }

} // end namespace OpenBabel

//! \file kekulize.cpp
//! \brief Alternate algorithm to kekulize a molecule (OBMol::NewPerceiveKekuleBonds()).
