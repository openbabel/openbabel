/**********************************************************************
kekulize.cpp - Kekulize a molecule

Copyright (C) 2016 Noel M. O'Boyle

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

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/kekulize.h>

namespace OpenBabel
{
  class Kekulizer
  {
  public:
    Kekulizer(OBMol* mol) : m_mol(mol) {}
    bool GreedyMatch();
  private:
    void AssignDoubleBonds(OBBitVec &doubleBonds);
    OBMol* m_mol;
  };

  static unsigned int TotalNumberOfBonds(OBAtom* atom)
  {
    return atom->GetImplicitHydrogen() + atom->GetValence();
  }

  static bool IsSpecialCase(OBAtom* atom, OBAtom* nbr)
  {
    switch (nbr->GetAtomicNum()) {
    case 8:
      // e.g. pyridine N-oxide as the double bond form
      if (atom->GetAtomicNum() == 7 && TotalNumberOfBonds(atom) == 3 &&
          atom->GetFormalCharge() == 0)
        return true;
      break;
    case 16:
      // ?? TODO ??
      if (TotalNumberOfBonds(atom) == 4)
        return true;
      break;
    }
    return false;
  }

  static bool NeedsDoubleBond(OBAtom* atom)
  {
    if (!atom->IsAromatic())
      return false;

    // Note to self:
    //   Bonds not in a ring should not be marked as aromatic - this should be done by the SMILES reader.
    //   If we have written the SMILES ourselves, and we ensure that we always write out the bond
    //   symbol between aromatic atoms, then bonds with a bond symbol cannot be aromatic - again
    //   this should be handled by the SMILES reader.

    // Does it already have an explicit double bond?
    FOR_BONDS_OF_ATOM(bond, atom) {
      if (bond->IsAromatic()) continue;
      OBAtom *nbr = bond->GetNbrAtom(atom);
      switch (bond->GetBondOrder()) {
      case -1: case 0: case 1:
        continue;
      case 2:
        if (IsSpecialCase(atom, nbr))
          return true;
        return false;
      default: // bond order > 2
        return false;
      }
    }

    // Is it one of the cases where we know that it only has single bonds?
    int chg = atom->GetFormalCharge();
    int deg = TotalNumberOfBonds(atom);
    switch (atom->GetAtomicNum()) {
    case 6:
      if (deg == 3 && (chg == 1 || chg == -1))
        return false;
      break;
    case 7: case 15: case 33: case 51:
      switch (chg) {
      case 0: // e.g. a pyrrole-type nitrogen
        if (deg == 3 || deg > 4)
          return false;
        break;
      case -1:
        if (deg == 2)
          return false;
        break;
      case 1:
        if (deg > 3)
          return false;
      }
      break;
    case 8: case 16: case 34: case 52:
      switch (chg) {
      case 0:
        if (deg == 2 || deg == 4 || deg > 5)
          return false;
        break;
      case -1: case 1:
        if (deg == 3 || deg == 5 || deg > 6)
          return false;
      }
    }

    return true; // It needs a double bond
  }

  static unsigned int GetMaxAtomIdx(OBMol* mol)
  {
    OBAtom* lastatom = (OBAtom*)0;
    FOR_ATOMS_OF_MOL(atom, mol) {
      lastatom = &(*atom);
    }
    return lastatom ? lastatom->GetIdx() : 0;
  }

  static unsigned int GetMaxBondIdx(OBMol* mol)
  {
    OBBond* lastbond = (OBBond*)0;
    FOR_BONDS_OF_MOL(bond, mol) {
      lastbond = &(*bond);
    }
    return lastbond ? lastbond->GetIdx() : 0;
  }

  class NodeIterator
  {
  public:
    NodeIterator(unsigned int *&degrees, unsigned int atomArraySize) :
      m_degrees(degrees), m_atomArraySize(atomArraySize),
      m_counter(0), finishedDegTwo(false)
    { }
    unsigned int next()
    {
      m_counter++;

      if (!finishedDegTwo) { // return deg 2 nodes first
        for (; m_counter < m_atomArraySize; ++m_counter) {
          if (m_degrees[m_counter] == 2) {
            return m_counter;
          }
        }
        finishedDegTwo = true;
        m_counter = 1; // first atom has idx 1
      }

      // return nodes with degree > 2
      for (; m_counter < m_atomArraySize; ++m_counter) {
        if (m_degrees[m_counter] > 2) {
          return m_counter;
        }
      }

      // Finished - return 0 signalling the end of iteration
      return 0;
    }
  private:
    unsigned int *&m_degrees; 
    unsigned int m_atomArraySize;
    unsigned int m_counter;
    bool finishedDegTwo;
  };

  void Kekulizer::AssignDoubleBonds(OBBitVec &doubleBonds)
  {
    int bit;
    for (bit = doubleBonds.FirstBit(); bit != doubleBonds.EndBit(); bit = doubleBonds.NextBit(bit)) {
      m_mol->GetBond(bit)->SetBondOrder(2);
    }
  }

  bool Kekulizer::GreedyMatch()
  {
    unsigned int atomArraySize = GetMaxAtomIdx(m_mol) + 1;
    unsigned int bondArraySize = GetMaxBondIdx(m_mol) + 1;

    // What atoms need a double bond? The job of kekulization is
    // to give all of these atoms a single double bond.
    OBBitVec needs_dbl_bond(atomArraySize); // defaults to all False
    FOR_ATOMS_OF_MOL(atom, m_mol) {
      if (NeedsDoubleBond(&*atom))
        needs_dbl_bond.SetBitOn(atom->GetIdx());
    }
    std::vector<bool> tmp;
    for (int i = 0; i < atomArraySize; ++i) {
      tmp.push_back(needs_dbl_bond.BitIsOn(i));
    }

    // Create lookup of degrees
    unsigned int *degrees = (unsigned int*)malloc(sizeof(unsigned int)*atomArraySize);
    std::vector<OBAtom*> degreeOneAtoms;
    FOR_ATOMS_OF_MOL(atom, m_mol) {
      if (!needs_dbl_bond.BitIsOn(atom->GetIdx())) continue;
      unsigned int mdeg = 0;
      FOR_BONDS_OF_ATOM(bond, &*atom) {
        if (!bond->IsAromatic()) continue;
        OBAtom *nbr = bond->GetNbrAtom(&*atom);
        if (needs_dbl_bond.BitIsOn(nbr->GetIdx()))
          mdeg++;
      }
      degrees[atom->GetIdx()] = mdeg;
      if (mdeg == 1)
        degreeOneAtoms.push_back(&*atom);
    }
    
    // Location of assigned double bonds
    OBBitVec doubleBonds(bondArraySize); // defaults to all False

    bool finished = false;
    while (true) { // Main loop
      
      // Complete all of the degree one nodes
      while (!degreeOneAtoms.empty()) {
        OBAtom* atom = degreeOneAtoms.back();
        degreeOneAtoms.pop_back();
        // some nodes may already have been handled
        if (!needs_dbl_bond.BitIsOn(atom->GetIdx())) continue;
        FOR_BONDS_OF_ATOM(bond, atom) {
          if (!bond->IsAromatic()) continue;
          OBAtom *nbr = bond->GetNbrAtom(&*atom);
          if (!needs_dbl_bond.BitIsOn(nbr->GetIdx())) continue;
          // create a double bond from atom -> nbr
          doubleBonds.SetBitOn(bond->GetIdx());
          needs_dbl_bond.SetBitOff(atom->GetIdx());
          needs_dbl_bond.SetBitOff(nbr->GetIdx());
          // now update degree information for nbr's neighbors
          FOR_BONDS_OF_ATOM(nbrbond, nbr) {
            if (&(*nbrbond) == &(*bond) || !nbrbond->IsAromatic()) continue;
            OBAtom* nbrnbr = nbrbond->GetNbrAtom(nbr);
            unsigned int nbrnbrIdx = nbrnbr->GetIdx();
            if (!needs_dbl_bond.BitIsOn(nbrnbrIdx)) continue;
            degrees[nbrnbrIdx]--;
            if (degrees[nbrnbrIdx] == 1)
              degreeOneAtoms.push_back(nbrnbr);
          }
          // only a single double bond can be made to atom so we can break here
          break;
        }
      }
      
      if (needs_dbl_bond.IsEmpty()) {
        finished = true;
        break;
      }

      // Now handle any remaining degree 2 or 3 nodes
      // We handle deg 2 nodes first and then 3, and the iteration over these nodes
      // is abstracted away. Once a double-bond is added that generates more
      // degree one nodes, then the iterator is exited
      NodeIterator iterator(degrees, atomArraySize);
      bool change = false;
      while (unsigned int atomIdx = iterator.next()) {
        if (!needs_dbl_bond.BitIsOn(atomIdx)) continue;
        // The following is almost identical to the code above for deg 1 atoms
        // except for handling the variable 'change'
        OBAtom *atom = m_mol->GetAtom(atomIdx);
        FOR_BONDS_OF_ATOM(bond, atom) {
          if (!bond->IsAromatic()) continue;
          OBAtom *nbr = bond->GetNbrAtom(&*atom);
          if (!needs_dbl_bond.BitIsOn(nbr->GetIdx())) continue;
          // create a double bond from atom -> nbr
          doubleBonds.SetBitOn(bond->GetIdx());
          needs_dbl_bond.SetBitOff(atomIdx);
          needs_dbl_bond.SetBitOff(nbr->GetIdx());
          // now update degree information for both atom's and nbr's neighbors
          for(int N=0; N<2; N++) {
            OBAtom *ref = N == 0 ? atom : nbr;
            FOR_BONDS_OF_ATOM(nbrbond, ref) {
              if (&(*nbrbond) == &(*bond) || !nbrbond->IsAromatic()) continue;
              OBAtom* nbrnbr = nbrbond->GetNbrAtom(ref);
              unsigned int nbrnbrIdx = nbrnbr->GetIdx();
              if (!needs_dbl_bond.BitIsOn(nbrnbrIdx)) continue;
              degrees[nbrnbrIdx]--;
              if (degrees[nbrnbrIdx] == 1) {
                degreeOneAtoms.push_back(nbrnbr);
                change = true;
              }
            }
          }
          // only a single double bond can be made to atom so we can break here
          break;
        }
        if (change)
          break; // exit the iterator once we have actually set a double bond
      }

      // We exit if we are finished or if no degree 2/3 nodes can be set
      if (!change)
        break;
    }

    // Assign double bonds
    AssignDoubleBonds(doubleBonds);

    // Tidy up
    free(degrees);

    return finished;
  }

  bool OBKekulize(OBMol* mol)
  {
    Kekulizer kekulizer(mol);
    bool success = kekulizer.GreedyMatch();
    return success;
  }


} // end namespace OpenBabel

//! \file kekulize.cpp
//! \brief Algorithm to kekulize a molecule
