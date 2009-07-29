/**********************************************************************
kekulize.cpp - Alternate algorithm to kekulize a molecule.

Copyright (C) 2004-2006 by Fabien Fontaine
Some portions Copyright (C) 2005-2009 by Geoffrey R. Hutchison
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

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/ring.h>
#include <openbabel/obiter.h>
// #include <openbabel/graphsym.h>

#include <sstream>

// Some names to make the code more readable
#define SINGLE 1
#define DOUBLE 2

#define NOT_IN_RINGS      -1
#define DOUBLE_ASSIGNED    0
#define DOUBLE_ALLOWED     1
#define DOUBLE_PROHIBITED  2

using namespace std;

namespace OpenBabel 
{

  // Modified internal-only versions
  int expand_cycle (OBMol *mol, OBAtom *atom, OBBitVec &avisit, OBBitVec &cvisit, 
                    int rootIdx, int prevAtomIdx = -1, int depth = 19);
  bool expandkekulize(OBMol *mol, int bond_idx,
                      std::vector<int> &atomState,
                      std::vector<int> &bondState);
  bool has_no_leftover_electrons(OBMol *mol, std::vector<int> &atomState);

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
      //      cout << "Checking for cycle at " << i << endl;
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
          //          cout << "rejected cycle" << endl;
          continue; // this atom isn't a root for the cycle
        }

        //        cout << " depth: " << depth << endl;
        //store the atoms of the cycle(s)
        unsigned int j;
        //        cout << " cycle: " << NumAtoms() << " ";
        for(j=1; j<= NumAtoms(); ++j) {
          if ( avisit[j] ) {
            //            cout << "\t" << j;
            atom = GetAtom(j);
            cycle.push_back(atom);
          }
        }
        //        cout << endl;

        // This isn't a real aromatic cycle -- give up
        // Fixes PR#1965566
        // Fixes PR#1784204
        if (cycle.size() < 3)
          continue;

        // At the beginning each atom give one electron to the cycle
        for(j=0; j< cycle.size(); ++j) {
          electron.push_back(1);
        }

        // remove one electron if the atom make a double bond out of the cycle
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

        //        cout << "minde before:" << minde << endl;
        // if huckel rule not satisfied some atoms must give more electrons
        while ( minde != 0 ) {
          bestorden=99;
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
            errorMsg << "Kekulize: Huckel rule not satisfied for molecule " << GetTitle() << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
            break;             // Huckel rule cannot be satisfied
          }                    // try to kekulize anyway
          else {
            electron[bestatom] += 1;
            minde--;
          }
        }

        if (bestorden == 99) { // Huckel rule not satisfied, just try to get an even number of electron before kekulizing

          electron = previousElectron; // restore electon's state

          int odd = sume % 2; 
          //cout << "odd:" << odd << endl;
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
              errorMsg << "Kekulize: Cannot get an even number of electron for molecule " << GetTitle() << "\n";
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

    std::vector<int> initAtomState;
    std::vector<int> atomState;
    std::vector<int> initBondState;
    std::vector<int> bondState;
    int Idx;
    OBAtom *atom, *atom2;
    OBBond *bond;

    // Initialize the atom and arrays
    initAtomState.resize(NumAtoms()+1);
    atomState.resize(NumAtoms()+1);
    for (int i = 0; i < NumAtoms()+1; ++i) {
      initAtomState[i] = NOT_IN_RINGS;
      atomState[i]     = NOT_IN_RINGS;
    }

    // Initialize the bond arrays with single bonds
    initBondState.resize(NumBonds());
    bondState.resize(NumBonds());
    for (int i = 0; i < NumBonds(); ++i) {
      initBondState[i] = SINGLE;
      bondState[i] = SINGLE;
    }

    // Figure out which atoms are in this ring system and whether or not each
    // atom can donate an electron.

    for (int i = 0; i < cycle.size(); ++i) {
      atom = cycle[i];
      Idx =  atom->GetIdx();
      if (electron[i] == 1) {
        initAtomState[Idx] = DOUBLE_ALLOWED;	// It has an electron it can donate
      } else {
        initAtomState[Idx] = DOUBLE_PROHIBITED;	// No electrons to contribute to aromatic system
      }
      atomState[Idx] = initAtomState[Idx];	// initialize atoms' current state too
    }

    // Do a recursive walk around the aromatic system, and see if we can find
    // an assignment of double/single bonds that works.
    atom = cycle[0];
    if (!expandkekulize(this, 0, atomState, bondState) ) {
      stringstream errorMsg;
      errorMsg << "Kekulize Error for molecule " << GetTitle() << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return;
    }

    // We found a successful assignment, now actually change the bonds to double

    if (DEBUG) {std::cout << "Set double bonds\n";}
    for (int i = 0; i < NumBonds(); ++i) {
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

  bool expandkekulize(OBMol *mol, int bond_idx,
                       std::vector<int> &atomState,
                       std::vector<int> &bondState)
  {
    // If all bonds are assigned, check that this is a sensible combination.
    // This is the "end of the line" for the recursion: Did this bond assignment
    // work or not?
    if (bond_idx >= bondState.size())
      return has_no_leftover_electrons(mol, atomState);

    // Get the bond and its atoms
    OBBond *bond = mol->GetBond(bond_idx);
    OBAtom *atom1 = bond->GetBeginAtom();
    OBAtom *atom2 = bond->GetEndAtom();
    int idx1 = atom1->GetIdx();
    int idx2 = atom2->GetIdx();

    // If this bond isn't part of the aromatic system, move on to the next bond.
    if (atomState[idx1] == NOT_IN_RINGS || atomState[idx2] == NOT_IN_RINGS)
      return expandkekulize(mol, bond_idx + 1, atomState, bondState);

    // Remember the current state so that we can backtrack if the attempt fails
    vector<int> previousState         = atomState;	// Backup the atom states
    vector<int> previousBondState     = bondState;	// ... and the bond states

    // Is a double bond allowed here?  Both atoms have to have an extra electron,
    // and not have been assigned a double bond from a previous step in recursion.

    if (atomState[idx1] == DOUBLE_ALLOWED && atomState[idx2] == DOUBLE_ALLOWED) {

      // Assign a double bond.
      atomState[idx1]  = DOUBLE_ASSIGNED;
      atomState[idx2]  = DOUBLE_ASSIGNED;
      bondState[bond_idx] = DOUBLE;        // set bond to double
      if (DEBUG) {std::cout << "bond " << bond_idx << " double\n";}

      // Recursively try the next bond
      if (expandkekulize(mol, bond_idx + 1, atomState, bondState))
        return true;

      // If the double bond didn't work, roll back the changes and try a single bond.
      atomState = previousState;
      bondState = previousBondState;
      if (DEBUG) {cout << "  double on bond " << bond_idx << " failed." << endl;}
    }

    // Double bond not allowed here, or double bond failed, just recurse with a single bond.
    if (DEBUG) {cout << "bond " << bond_idx << " single" << endl;}
    if (expandkekulize(mol, bond_idx + 1, atomState, bondState))
      return true;

    // If it didn't work, roll back the changes we made and return failure.
    if (DEBUG) {cout << "bond " << bond_idx << " single failed, rolling back changes" << endl;}
    atomState = previousState;
    bondState = previousBondState;
    return false;
  }

  // Check for leftover electrons.  This is used during expandkekulize() above
  // to make sure all of the 4n+2 electrons that were available for bonding in the
  // aromatic ring system were actually used during the assignment of single
  // and double bonds.
  bool has_no_leftover_electrons(OBMol *mol, std::vector<int> &atomState)
  {
    FOR_ATOMS_OF_MOL(a, mol) {
      int idx = a->GetIdx();
      if (atomState[idx] == DOUBLE_ALLOWED) {
        if (DEBUG) {cout << "  failure, extra electron on atom " << a->GetIdx() << endl;}
        return false;
      }
    }
    return true;	// no extra electrons found
  }

  int OBMol::expand_kekulize(OBAtom *atom1, OBAtom *atom2, std::vector<int> &currentState,
                             std::vector<int> &initState, std::vector<int> &bcurrentState,
                             std::vector<int> &binitState, std::vector<bool> &mark)
  {
    return 0; // deprecated over expandkekulize above
  }

  //! \brief Give the priority to give two electrons instead of 1
  //!
  //! Higher numbers will have the precedence to give two electrons (and get a double bond)
  //! instead of one
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

  // Keep track of which rings contain *all* atoms in the cycle
  // This method essentially does a modified depth-first search to find
  //  large aromatic cycles
  //! Recursively find the aromatic atoms with an aromatic bond to the current atom
  int expand_cycle (OBMol *mol, OBAtom *atom, OBBitVec &avisit, OBBitVec &cvisit, 
                    int rootIdx, int prevAtomIdx, int depth)
  {
    // early termination
    if (depth < 0)
      return depth;

    //    cout << " expand_cycle: " << atom->GetIdx() << " depth " << depth << endl;

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
        //        cout << " checking: " << natom << " bo: " << (*i)->GetBO() << endl;
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
          //          cout << " score: " << trialScore << endl;
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
