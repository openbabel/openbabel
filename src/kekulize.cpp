/**********************************************************************
kekulize.cpp - Alternate algorithm to kekulize a molecule.

Copyright (C) 2004-2006 by Fabien Fontaine
Some portions Copyright (C) 2005-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2009 by Craig A. James
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

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

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/ring.h>
#include <openbabel/obiter.h>
#include <openbabel/graphsym.h>

#include <sstream>

// Some names to make the code more readable
#define UNASSIGNED 0
#define SINGLE     1
#define DOUBLE     2

#define NOT_IN_RINGS      -1
#define DOUBLE_ASSIGNED    0
#define DOUBLE_ALLOWED     1
#define DOUBLE_PROHIBITED  2

using namespace std;

namespace OpenBabel 
{

  // Modified internal-only version
  // Keep track of which rings contain *all* atoms in the cycle
  // This method essentially does a modified depth-first search to find
  //  large aromatic cycles
  int expand_cycle (OBMol *mol, OBAtom *atom, OBBitVec &avisit, OBBitVec &cvisit, 
                    int rootIdx, int prevAtomIdx = -1, int depth = 19);

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

    // Find all the groups of aromatic cycle
    for(i=1; i<= NumAtoms(); i++ ) {
      atom = GetAtom(i);
      if (DEBUG) {cout << "Checking for cycle at " << i << endl;}
      if (atom->HasAromaticBond() && !cvisit[i]) { // is new aromatic atom of an aromatic cycle ?

        avisit.Clear();
        electron.clear();
        cycle.clear();
      
        avisit.SetBitOn(i);

        if (fusedRings)
          expandcycle(atom, avisit);
        else {
          int depth = expand_cycle(this, atom, avisit, cvisit, atom->GetIdx());
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
          if (DEBUG) {cout << "rejected cycle" << endl;}
          continue; // this atom isn't a root for the cycle
        }

        //store the atoms of the cycle(s)
        unsigned int j;
        if (DEBUG) {cout << " cycle: " << NumAtoms() << " ";}
        for(j=1; j<= NumAtoms(); ++j) {
          if ( avisit[j] ) {
            if (DEBUG) {cout << "\t" << j;}
            atom = GetAtom(j);
            cycle.push_back(atom);
          }
        }
        if (DEBUG) {cout << endl;}

        // This isn't a real aromatic cycle -- give up
        // Fixes PR#1965566
        // Fixes PR#1784204
        if (cycle.size() < 3)
          continue;

        // At the beginning each atom give one electron to the cycle
        for(j=0; j< cycle.size(); ++j) {
          electron.push_back(1);
        }
      
        // Remove one electron if the atom makes a double bond out of the cycle
        sume =0;
        for(j=0; j< cycle.size(); ++j) {
          atom = cycle[j];
          for(bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
            if ( bond->IsDouble() ) {
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
          sume += electron[j];
        }

        // Save the electron state in case huckel rule is not satisfied
        vector<int> previousElectron = electron;

        // find the ideal number of electrons according to the huckel 4n+2 rule
        minde=99;
        for (i=1; 1; ++i) {
          n = 4 *i +2;
          de = n - sume;
          if (  de < minde )
            minde=de;
          else if ( minde < 0 )
            minde=de;
          else
            break;
        }
      
        stringstream errorMsg;

        if (DEBUG) {cout << "minde before:" << minde << endl;}
        // if huckel rule not satisfied some atoms must give more electrons
        while ( minde != 0 ) {
          bestorden=99;
          for(j=0; j< cycle.size(); ++j) {
            if (electron[j] == 1) {
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

        if (bestorden == 99) {	// Huckel rule not satisfied, just try to get an
				// even number of electron before kekulizing
	
          electron = previousElectron; // restore electon's state
	
          int odd = sume % 2; 
          if (DEBUG) {cout << "odd:" << odd << endl;}
          if(odd) { // odd number of electrons try to add an electron to the best possible atom
            for(j=0; j< cycle.size(); ++j) {
              if (electron[j] == 1) {
                orden = getorden(cycle[j]);
                if (orden < bestorden) {
                  bestorden = orden;
                  bestatom = j;
                }
              }
            }
            if (bestorden==99) {  // no electron giving atom found
              errorMsg << "Kekulize: Cannot get an even number of electron for molecule " << GetTitle() << endl;
              obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
              break;             // impossible to choose an atom to obtain an even number of electron
            }                    // try to kekulize anyway
            else {
              electron[bestatom] += 1;
            }
          }
        }

        /*
          cout << "minde after:" << minde <<endl;
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
        */

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

  void OBMol::start_kekulize( std::vector <OBAtom*> &cycle, std::vector<int> &electron) {
  
    std::vector<int> atomState;
    std::vector<int> bondState;
    int idx;
    OBAtom *atom;
    OBBond *bond;

    if (DEBUG) {cout << "-------------------- start_kekulize:" << endl;}

    // Initially mark all bonds "unassigned" (not single or double)
    bondState.resize(NumBonds());
    for (int i = 0; i < NumBonds(); ++i) {
      bondState[i] = UNASSIGNED;
    }
  
    // Figure out which atoms are in this ring system and whether or not each
    // atom can donate an electron.

    atomState.resize(NumAtoms()+1);
    for (int i = 0; i < NumAtoms()+1; ++i) {
      atomState[i] = NOT_IN_RINGS;
    }
    for (int i = 0; i < cycle.size(); ++i) {
      atom = cycle[i];
      idx =  atom->GetIdx();
      if (electron[i] == 1) {
        atomState[idx] = DOUBLE_ALLOWED;	// It has an electron it can donate
      } else {
        atomState[idx] = DOUBLE_PROHIBITED;	// No electrons to contribute to aromatic system
      }
      if (DEBUG) {cout << "atom " << idx << ": initial state = " << atomState[idx] << endl;}
    }

    // Get the LSSR, and find all rings that are completely in cycle[], that
    // is, every atom in the ring is in this aromatic ring system.

    std::vector<OBRing*> lssr_tmp = GetLSSR();
    std::vector<OBRing*> lssr;
    std::vector<bool>atom_in_cycle(NumAtoms()+1);
    for (int i = 0; i <= NumAtoms(); i++)
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
    for (int i = 0; i < lssr.size(); i++)
      lssrAssigned[i] = false;

    // Initialize an empty vector that expand_kekulize2() uses.
    std::vector<OBBond*> bondsThisRing;
    bondsThisRing.clear();

    // Do a recursive walk across the aromatic system, and see if we can find
    // an assignment of double/single bonds that works.
    if (!expand_kekulize2(atomState, bondState, lssr, lssrAssigned, bondsThisRing)) {
      stringstream errorMsg;
      errorMsg << "Kekulize Error for molecule " << GetTitle() << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return;
    }

    // We found a successful assignment, now actually change the bonds to double

    if (DEBUG) {std::cout << "Set double bonds" << endl;}
    for (int i = 0; i < NumBonds(); ++i) {
      bond = GetBond(i);    
      if (bond->GetBO() == 5) {
	if (bondState[i] == DOUBLE) {
	  if (   bond->GetBeginAtom()->IsSulfur()
	      && bond->GetEndAtom()->IsSulfur()) {
	    continue;	    // no double bonds between aromatic sulfur atoms
	  }
	  bond->SetKDouble();
	  bond->SetBO(2);
	  if (DEBUG) { std::cout << "  bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << " double" << endl;}
	}
	else {
	  if (DEBUG) {	// a lot of complication just to avoid printing confusing debug info...
	    if (   atomState[bond->GetBeginAtom()->GetIdx()] != NOT_IN_RINGS
		&& atomState[bond->GetEndAtom()->GetIdx()] != NOT_IN_RINGS) {
	      if (DEBUG) { std::cout << "  bond " << bond->GetBeginAtomIdx() << " " << bond->GetEndAtomIdx() << " single" << endl;}
	    }
	  }
	}
      }
    }
    return;
  }


  // deprecated!!
  bool OBMol::expand_kekulize(int bond_idx,
                              std::vector<int> &atomState,
                              std::vector<int> &bondState) {}
  // deprecated!!
  bool OBMol::has_no_leftover_electrons(std::vector<int> &atomState) {}

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

  bool OBMol::expand_kekulize2(std::vector<int> &atomState,
			      std::vector<int> &bondState,
			      std::vector<OBRing*> &lssr,
			      std::vector<bool> &lssrAssigned,
			      std::vector<OBBond *> &bondsThisRing)
  {
    OBBond *bond = NULL;
    OBAtom *atom1, *atom2;
    int idx1, idx2, bond_idx;

    if (DEBUG) {cout << "---------- expand_kekulize2:" << endl;}

    // Remember the current state so that we can backtrack if the attempt fails
    vector<int> previousState         = atomState;	// Backup the atom states
    vector<int> previousBondState     = bondState;	// ... and the bond states

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

      // Does a double bond work here?
      if (   atomState[idx1] == DOUBLE_ALLOWED
	  && atomState[idx2] == DOUBLE_ALLOWED) {

	// Assign a double bond.
	atomState[idx1]  = DOUBLE_ASSIGNED;
	atomState[idx2]  = DOUBLE_ASSIGNED;
	bondState[bond_idx] = DOUBLE;        // set bond to double
	if (DEBUG) {std::cout << "bond " << bond_idx << " (atoms " << idx1 << " to " << idx2 << ") double" << endl;}

	// Recursively try the next bond
	if (expand_kekulize2(atomState, bondState, lssr, lssrAssigned, bondsThisRing))
	  return true;

	// If the double bond didn't work, roll back the changes and try a single bond.
	atomState = previousState;
	bondState = previousBondState;
	if (DEBUG) {cout << "  double on bond " << bond_idx << " failed." << endl;}
      }

      // Either a double bond not allowed here, or double bond failed, try a single bond.
      bondState[bond_idx] = SINGLE;        // set bond to single
      if (DEBUG) {cout << "bond " << bond_idx << " (atoms " << idx1 << " to " << idx2 << ") single" << endl;}

      // Recursively try the next bond
      if (expand_kekulize2(atomState, bondState, lssr, lssrAssigned, bondsThisRing))
	return true;

      // If it didn't work, roll back the changes we made and return failure.
      if (DEBUG) {cout << "bond " << bond_idx << " single failed, rolling back changes" << endl;}
      atomState = previousState;
      bondState = previousBondState;
      return false;
    }

    // If we get here, it means we successfully assigned all of the bonds on
    // the current ring.  Double check that we haven't made an impossible
    // bond assignment.

    if (has_leftover_electrons(atomState, bondState)) {
      atomState = previousState;
      bondState = previousBondState;
      return false;
    }

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

    for (int rnum = 0; rnum < lssr.size(); rnum++) {
      if (lssrAssigned[rnum])
	continue;
      OBRing *ring = lssr[rnum];
      get_bonds_of_ring(ring, bondsThisRing);
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
    if (!expand_kekulize2(atomState, bondState, lssr, lssrAssigned, best_bondsThisRing)) {
      lssrAssigned[best_rnum] = false;
      return false;
    } else {
      return true;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Check for leftover electrons.  This is used during expand_kekulize2() above
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

  bool OBMol::has_leftover_electrons(std::vector<int> &atomState,
				     std::vector<int> &bondState)
  {
    FOR_ATOMS_OF_MOL(a, this) {
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


  ////////////////////////////////////////////////////////////////////////////////
  //! Returns a vector of the bonds in a ring, in circular order around the ring,
  //
  // The circular order (the fact that the bonds are end-to-end in the
  // vector) is important when assigning single/double bonds to an aromatic
  // system, because it helps you "walk" around a ring assigning
  // alternating single/double bonds, and you can often get it right on the
  // first try.

  bool OBMol::get_bonds_of_ring(OBRing *ring, std::vector<OBBond *> &ring_bonds) {

    vector<bool> used_bonds;
    used_bonds.resize(NumBonds());
    for (int i = 0; i < used_bonds.size(); i++)
      used_bonds[i] = false;

    ring_bonds.clear();
    int start_idx = ring->_path[0];
    OBAtom *start_atom = GetAtom(start_idx);
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

  ////////////////////////////////////////////////////////////////////////////////
  // Count the number of bonds in one ring that are already assigned
  // (single or double) during Kekule assignment for a ring system.  For
  // example, in naphthalene, if we've already assigned single/double to
  // the first ring, then calling this function for the second ring will
  // return 1: one bond is already assigned.

  int OBMol::count_assigned_bonds(std::vector<OBBond *> &bondsThisRing,
				   std::vector<int> &bondState)
  {
    int n = 0;
    for (int i = 0; i < bondsThisRing.size(); i++) {
      if (bondState[bondsThisRing[i]->GetIdx()] != UNASSIGNED) {
	n++;
      }
    }
    return n;
  }


  //! Give the priority to give two electrons instead of 1

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
  bool OBMol::expandcycle (OBAtom *atom, OBBitVec &avisit, OBAtom *, int)
  {
    OBAtom *nbr;
    //  OBBond *bond;
    std::vector<OBBond*>::iterator i;
    int natom;
    //for each neighbour atom test if it is in the aromatic ring
    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
      {
        natom = nbr->GetIdx();
        // if (!avisit[natom] && nbr->IsAromatic() && ((OBBond*) *i)->IsAromatic()) {
        if (!avisit[natom] && ((OBBond*) *i)->GetBO()==5 
            && ((OBBond*) *i)->IsInRing()) {
          avisit.SetBitOn(natom);
          expandcycle(nbr, avisit);
        }
      }

    return true;
  }

  //! Recursively find the aromatic atoms with an aromatic bond to the current atom
  int expand_cycle (OBMol *mol, OBAtom *atom, OBBitVec &avisit, OBBitVec &cvisit, 
                    int rootIdx, int prevAtomIdx, int depth)
  {
    // early termination
    if (depth < 0)
      return depth;

    if (DEBUG) {cout << " expand_cycle: " << atom->GetIdx() << " depth " << depth << endl;}

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
        natom = nbr->GetIdx();
        if (DEBUG) {cout << " checking: " << natom << " bo: " << (*i)->GetBO() << endl;}
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
        trialScore = expand_cycle(mol, nbr, trialMatch, cvisit, rootIdx, atom->GetIdx(), depth - 1);
        if (trialScore > 0 && trialScore < bestScore) { // we found a larger, valid cycle
          if (DEBUG) {cout << " score: " << trialScore << endl;}
          bestMatch = trialMatch;
          bestScore = trialScore;
        }
      } // check all neighbors

    if (bestScore <= 0 || bestScore == 1000)
      return -1; // e.g., we have no valid choices

    avisit = bestMatch; // remember our path
    return bestScore;
  }

} // end namespace OpenBabel

//! \file kekulize.cpp
//! \brief Alternate algorithm to kekulize a molecule (OBMol::NewPerceiveKekuleBonds()).
