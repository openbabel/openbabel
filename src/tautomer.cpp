/**********************************************************************
tautomer.cpp - Tautomer support

  Copyright (C) 2011,2020 by Tim Vandermeersch

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

#include <openbabel/tautomer.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/graphsym.h>
#include <openbabel/canon.h>
#include <openbabel/obconversion.h>
#include <cassert>
#include <algorithm>

namespace OpenBabel {



//#define DEBUG 1

  /**
   * http://www.daylight.com/meetings/emug99/Delany/taut_html/index.htm
   */
  struct TautomerImpl
  {

#ifdef DEBUG
    static std::string indentation(int n)
    {
      return std::string(n * 4, ' ');
    }

    // Use RAII to keep track of enetering and leaving a function
    class EnterExit
    {
      public:
        explicit EnterExit(const std::string &function, int indent = 0) : m_function(function), m_indent(indent)
        {
          std::cout << ::OpenBabel::TautomerImpl::indentation(m_indent) << "Enter " << m_function << "()..." << std::endl;
        }

        ~EnterExit()
        {
          std::cout << ::OpenBabel::TautomerImpl::indentation(m_indent) << "Exit " << m_function << "()..." << std::endl;
        }

        std::string indentation() const
        {
          return ::OpenBabel::TautomerImpl::indentation(m_indent + 1);
        }

      private:
        std::string m_function;
        int m_indent;
    };
#endif

    /*   H
     *   |        nitrogen donor
     * --N--
     *
     * ==N--      nitrogen acceptor
     *
     * --O--H     oxygen(*) donor
     *
     * ==O        oxygen(*) acceptor
     *
     * (*) sulfur, selenium and tellurium behave the same as oxygen
     */
    enum AtomNumber {
      Carbon = 6,
      Nitrogen = 7,
      Oxygen = 8,
      Sulfur = 16,
      Selenium = 34,
      Tellurium = 52
    };

    enum Type {
      Donor,
      Acceptor,
      Hybridized,
      Other,
      Assigned,
      Unassigned,
      Single,
      Double
    };

    /**
     * Use RAII in the backtracking algoirthm to reduce bugs.
     */
    class AssignDonorRAII
    {
      public:
        AssignDonorRAII(std::vector<Type> &atomTypes, int &hydrogenCounter, OBAtom *atom) : m_atomTypes(atomTypes), m_atom(atom), m_hydrogenCounter(hydrogenCounter)
        {
          redo(atom);
        }

        ~AssignDonorRAII()
        {
          undo();
        }

        void redo(OBAtom *atom)
        {
#ifdef DEBUG
          std::cout << "AssignDonorRAII::redo(" << atom->GetIndex() << ")" << std::endl;
#endif
          assert(m_hydrogenCounter >= 1);
          m_atom = atom;
          m_atomTypes[m_atom->GetIndex()] = Donor;
          m_atom->SetImplicitHCount(m_atom->GetImplicitHCount() + 1);
          --m_hydrogenCounter;
        }

        void undo()
        {
          if (!m_atom)
            return;
#ifdef DEBUG
          std::cout << "AssignDonorRAII::undo(" << m_atom->GetIndex() << ")" << std::endl;
#endif
          m_atomTypes[m_atom->GetIndex()] = Unassigned;
          m_atom->SetImplicitHCount(m_atom->GetImplicitHCount() - 1);
          ++m_hydrogenCounter;
          m_atom = nullptr;
        }

        /**
         * For the canonical tautomer we do not use the functor and exit search
         * with canonical tautomer. Therefore, the action from redo (increment
         * implicit hydrogen count) should not be undone. In other words, we
         * release or "leak" the resource.
         */
        void release()
        {
          m_atom = nullptr;
        }

      private:
        std::vector<Type> &m_atomTypes;
        OBAtom *m_atom;
        int &m_hydrogenCounter;
    };

    /**
     * Use RAII in the backtracking algoirthm to reduce bugs.
     */
    class AssignAcceptorRAII
    {
      public:
        AssignAcceptorRAII(std::vector<Type> &atomTypes, OBAtom *atom) : m_atomTypes(atomTypes), m_atom(atom)
        {
          redo(atom);
        }

        ~AssignAcceptorRAII()
        {
          undo();
        }

        void redo(OBAtom *atom)
        {
#ifdef DEBUG
          std::cout << "AssignAcceptorRAII::redo(" << atom->GetIndex() << ")" << std::endl;
#endif
          m_atom = atom;
          m_atomTypes[m_atom->GetIndex()] = Acceptor;
        }

        void undo()
        {
          if (!m_atom)
            return;
#ifdef DEBUG
          std::cout << "AssignAcceptorRAII::undo(" << m_atom->GetIndex() << ")" << std::endl;
#endif
          m_atomTypes[m_atom->GetIndex()] = Unassigned;
          m_atom = nullptr;
        }

      private:
        std::vector<Type> &m_atomTypes;
        OBAtom *m_atom;
    };

    /**
     * Use RAII in the backtracking algoirthm to reduce bugs.
     */
    class PropagationRAII
    {
      public:
        PropagationRAII(std::vector<Type> &atomTypes, std::vector<Type> &bondTypes, int &hydrogenCounter)
          : m_atomTypes(atomTypes), m_bondTypes(bondTypes), m_hydrogenCounter(hydrogenCounter)
        {
        }

        ~PropagationRAII()
        {
          // donors
          for (std::size_t i = 0; i < m_donors.size(); ++i) {
            m_atomTypes[m_donors[i]->GetIndex()] = Unassigned;
            m_donors[i]->SetImplicitHCount(m_donors[i]->GetImplicitHCount() - 1);
          }
          m_hydrogenCounter += m_donors.size();
          // acceptors
          for (std::size_t i = 0; i < m_acceptors.size(); ++i)
            m_atomTypes[m_acceptors[i]->GetIndex()] = Unassigned;
          // bonds
          for (std::size_t i = 0; i < m_bonds.size(); ++i)
            m_bondTypes[m_bonds[i]->GetIdx()] = Unassigned;
        }

        void assignDonor(OBAtom *atom)
        {
          assert(m_hydrogenCounter >= 1);
          m_donors.push_back(atom);
          m_atomTypes[atom->GetIndex()] = Donor;
          atom->SetImplicitHCount(atom->GetImplicitHCount() + 1);
          --m_hydrogenCounter;
        }

        void assignAcceptor(OBAtom *atom)
        {
          m_acceptors.push_back(atom);
          m_atomTypes[atom->GetIndex()] = Acceptor;
        }

        void assignBond(OBBond *bond, Type type)
        {
          m_bonds.push_back(bond);
          m_bondTypes[bond->GetIdx()] = type;
        }

        // see AssignDonorRAII::release()
        void release()
        {
          std::copy(m_donors.begin(), m_donors.end(), std::back_inserter(m_acceptors));
          m_donors.clear();
        }

      private:
        std::vector<Type> &m_atomTypes;
        std::vector<Type> &m_bondTypes;
        std::vector<OBAtom*> m_acceptors; //!< propagated atoms
        std::vector<OBAtom*> m_donors; //!< propagated atoms
        std::vector<OBBond*> m_bonds; //!< propagated bonds
        int &m_hydrogenCounter; //!< global hydrogen counter
    };


    bool m_canonical, m_foundLeafNode;
    std::vector<OBAtom*> m_canonAtoms;

#ifdef DEBUG
    void PrintAtomTypes(const std::vector<Type> &atomTypes, int indent = 0)
    {
      std::cout << indentation(indent);
      for(std::size_t i = 0; i < atomTypes.size(); ++i) {
        std::cout << i;
        switch (atomTypes[i]) {
          case Donor:
            std::cout << "D ";
            break;
          case Acceptor:
            std::cout << "A ";
            break;
          case Hybridized:
            std::cout << "H ";
            break;
          case Other:
            std::cout << "O ";
            break;
          case Unassigned:
            std::cout << "? ";
            break;
          default:
            assert(0);
        }
      }
      std::cout << std::endl;
    }

    void PrintBondTypes(OBMol *mol, const std::vector<Type> &bondTypes, int indent = 0)
    {
      std::cout << indentation(indent);
      for(std::size_t i = 0; i < bondTypes.size(); ++i) {
        OBBond *bond = mol->GetBond(i);
        std::cout << bond->GetBeginAtomIdx() - 1;

        switch (bondTypes[i]) {
          case Assigned:
            switch (bond->GetBondOrder()) {
              case 1:
                std::cout << "-";
                break;
              case 2:
                std::cout << "=";
                break;
              case 3:
                std::cout << "#";
                break;
              default:
                assert(0);
            }
            break;
          case Single:
            std::cout << "-";
            break;
          case Double:
            std::cout << "=";
            break;
          case Unassigned:
            std::cout << "?";
            break;
          default:
            assert(0);
        }

        std::cout << bond->GetEndAtomIdx() - 1 << " ";
      }
      std::cout << std::endl;
    }
#endif

    std::vector<Type> InitializeAtomTypes(OBMol *mol)
    {
      std::vector<Type> types;

      FOR_ATOMS_OF_MOL (atom, mol) {
        switch (atom->GetAtomicNum()) {
          case Carbon:
            if (atom->HasDoubleBond())
              types.push_back(Hybridized);
            else
              types.push_back(Other);
            break;

          case Nitrogen:
            if (atom->HasDoubleBond() && atom->GetTotalValence() == 3)
              types.push_back(Acceptor);
            else if (atom->GetImplicitHCount() && atom->GetTotalValence() == 3)
              types.push_back(Donor);
            else
              types.push_back(Other);
            break;

          case Oxygen:
          case Sulfur:
          case Selenium:
          case Tellurium:
            if (atom->HasDoubleBond() && atom->GetTotalValence() == 2)
              types.push_back(Acceptor);
            else if (atom->GetImplicitHCount() && atom->GetTotalValence() == 2)
              types.push_back(Donor);
            else
              types.push_back(Other);
            break;

          default:
            types.push_back(Other);
        }
      }

      assert( types.size() == mol->NumAtoms() );
      return types;
    }

    std::vector<Type> InitializeBondTypes(OBMol *mol, const std::vector<Type> &atomTypes)
    {
      std::vector<Type> bondTypes;

      FOR_BONDS_OF_MOL (bond, mol) {
        assert( bond->GetBeginAtomIdx() <= mol->NumAtoms() );
        assert( bond->GetEndAtomIdx() <= mol->NumAtoms() );
        if (atomTypes[bond->GetBeginAtomIdx()-1] != Other && atomTypes[bond->GetEndAtomIdx()-1] != Other)
          bondTypes.push_back(Unassigned);
        else
          bondTypes.push_back(Assigned);
      }

      return bondTypes;
    }

    bool IsPropagationValid(OBMol *mol, const std::vector<Type> &atomTypes, const std::vector<Type> &bondTypes, int depth)
    {
      FOR_ATOMS_OF_MOL (atom, mol) {
        // Count atom's bond types
        int numSingleBond = 0, numDoubleBond = 0, numUnassigned = 0;
        FOR_BONDS_OF_ATOM (bond, &*atom) {
          if (bondTypes[bond->GetIdx()] == Single)
            numSingleBond++;
          if (bondTypes[bond->GetIdx()] == Double)
            numDoubleBond++;
          if (bondTypes[bond->GetIdx()] == Unassigned)
            numUnassigned++;
        }

        // A donor can not have any double bonds
        if (atomTypes[atom->GetIndex()] == Donor)
          if (numDoubleBond) {
#ifdef DEBUG
            std::cout << indentation(depth) << "invalid Donor" << std::endl;
#endif
            return false;
          }

        // An acceptor or hybridized atom should:
        if (atomTypes[atom->GetIndex()] == Acceptor || atomTypes[atom->GetIndex()] == Hybridized) {
          // have at least one unassigned or double bond
          if (!numUnassigned && !numDoubleBond) {
#ifdef DEBUG
            std::cout << indentation(depth) << "invalid Acceptor/Hybridized [no unassigned/double bond] for atom " << atom->GetIndex() << std::endl;
#endif
            return false;
          }
          // have no more than one double bond
          if (numDoubleBond > 1) {
#ifdef DEBUG
            std::cout << indentation(depth) << "invalid Acceptor/Hybridized [multiple double bonds] for atom " << atom->GetIndex() << std::endl;
#endif
            return false;
          }
        }
      }

      return true;
    }

    bool IsLeafNode(const std::vector<Type> &atomTypes) const
    {
      // Check to see if we are at a leaf node (i.e. no unassigned atoms left)
      for (std::size_t i = 0; i < atomTypes.size(); ++i)
        if (atomTypes[i] == Unassigned)
          return false;
      return true;
    }

    void AssignmentPropagation(OBMol *mol, std::vector<Type> &atomTypes, std::vector<Type> &bondTypes,
        int &numHydrogens, TautomerFunctor &functor, int depth)
    {
#ifdef DEBUG
      EnterExit ee("AssignmentPropagation", depth);
#endif
      PropagationRAII propagation(atomTypes, bondTypes, numHydrogens);

      bool changed = true;
      while (changed) {
        changed = false;

        // Assigning an atom Donor causes all of its bonds to be assigned single
        FOR_ATOMS_OF_MOL (atom, mol) {
          if (atomTypes[atom->GetIndex()] == Donor) {
            FOR_BONDS_OF_ATOM (bond, &*atom) {
              if (bondTypes[bond->GetIdx()] == Unassigned) {
                propagation.assignBond(&*bond, Single);
                changed = true;
#ifdef DEBUG
                std::cout << ee.indentation() << "-> Rule 1: Assign " << bond->GetBeginAtomIdx()-1 << "-" << bond->GetEndAtomIdx()-1 << " Single" << std::endl;
#endif
              }
            }
          }
        }

        // Any unassigned atom that has all of its bonds assigned single, must be assigned a donor
        if (numHydrogens) {
          FOR_ATOMS_OF_MOL (atom, mol) {
            if (atomTypes[atom->GetIndex()] != Unassigned)
              continue;
            bool allBondsSingle = true;
            FOR_BONDS_OF_ATOM (bond, &*atom) {
              if (bondTypes[bond->GetIdx()] != Single) {
                allBondsSingle = false;
                break;
              }
            }

            if (allBondsSingle) {
              propagation.assignDonor(&*atom);
              changed = true;
#ifdef DEBUG
              std::cout << ee.indentation() << "-> Rule 2: Assign " << atom->GetIndex() << " Donor" << std::endl;
#endif
            }
          }
        }

        // Any unassigned atom with an assigned double bond must be assigned an acceptor
        FOR_ATOMS_OF_MOL (atom, mol) {
          if (atomTypes[atom->GetIndex()] != Unassigned)
            continue;
          bool hasDoubleBond = false;
          FOR_BONDS_OF_ATOM (bond, &*atom) {
            if (bondTypes[bond->GetIdx()] == Double) {
              hasDoubleBond = true;
              break;
            }
          }

          if (hasDoubleBond) {
            propagation.assignAcceptor(&*atom);
            changed = true;
#ifdef DEBUG
            std::cout << ee.indentation() << "-> Rule 3: Assign " << atom->GetIndex() << " Acceptor" << std::endl;
#endif
          }
        }

        // A hybridized or acceptor atom with assigned double bond assigns all the remaining unassigned bonds single
        FOR_ATOMS_OF_MOL (atom, mol) {
          if (atomTypes[atom->GetIndex()] != Acceptor && atomTypes[atom->GetIndex()] != Hybridized)
            continue;
          bool hasDoubleBond = false;
          FOR_BONDS_OF_ATOM (bond, &*atom) {
            if (bondTypes[bond->GetIdx()] == Double) {
              hasDoubleBond = true;
              break;
            }
          }

          if (hasDoubleBond) {
            FOR_BONDS_OF_ATOM (bond, &*atom) {
              if (bondTypes[bond->GetIdx()] == Unassigned) {
                propagation.assignBond(&*bond, Single);
                changed = true;
#ifdef DEBUG
                std::cout << ee.indentation() << "-> Rule 4: Assign " << bond->GetBeginAtomIdx()-1 << "-" << bond->GetEndAtomIdx()-1 << " Single" << std::endl;
#endif
              }
            }
          }
        }

        // A hybridized or acceptor atom with a single unassigned bond assigns it double if it doesn't already have one
        FOR_ATOMS_OF_MOL (atom, mol) {
          if (atomTypes[atom->GetIndex()] != Acceptor && atomTypes[atom->GetIndex()] != Hybridized)
            continue;
          bool hasDoubleBond = false;
          int numUnassigned = 0;
          FOR_BONDS_OF_ATOM (bond, &*atom) {
            if (bondTypes[bond->GetIdx()] == Double)
              hasDoubleBond = true;
            if (bondTypes[bond->GetIdx()] == Unassigned)
              numUnassigned++;
          }


          if (!hasDoubleBond && numUnassigned == 1) {
            FOR_BONDS_OF_ATOM (bond, &*atom) {
              if (bondTypes[bond->GetIdx()] == Unassigned) {
                propagation.assignBond(&*bond, Double);
                changed = true;
#ifdef DEBUG
                std::cout << ee.indentation() << "-> Rule 5: Assign " << bond->GetBeginAtomIdx()-1 << "-" << bond->GetEndAtomIdx()-1 << " Double" << std::endl;
#endif
              }
            }
          }
        }

      } // while (changed)

      // An invalid atom type causes the search to backtrack
      if (!IsPropagationValid(mol, atomTypes, bondTypes, depth + 1))
        return;

      // Check to see if we are at a leaf node (i.e. no unassigned atoms left)
      if (IsLeafNode(atomTypes)) {

        // Check to make sure there are no unassigned bonds remaining
        // This happens for phenyl rings...
        bool unassignedBonds = false;
        for (std::size_t i = 0; i < bondTypes.size(); ++i)
          if (bondTypes[i] == Unassigned) {
            unassignedBonds = true;
            // Randomly assign a single bond, the other bonds will be propagated
            OBBond *bond = mol->GetBond(i);
            propagation.assignBond(bond, Single);
#ifdef DEBUG
            std::cout << ee.indentation() << "-> Unasigned bonds remaining, assigning " << bond->GetBeginAtomIdx()-1 << "-" << bond->GetEndAtomIdx()-1 << " Single" << std::endl;
#endif
            break;
          }

        if (unassignedBonds) {
          AssignmentPropagation(mol, atomTypes, bondTypes, numHydrogens, functor, depth + 1);
          if (m_canonical && m_foundLeafNode)
            propagation.release();
          return;
        }

        //
        if (!numHydrogens) {
          m_foundLeafNode = true;

#ifdef DEBUG
          std::cout << ee.indentation() << "  --> LeafNode reached..." << std::endl;
#endif
          // Copy information to mol and reset aromaticity
          FOR_ATOMS_OF_MOL (atom, mol)
            atom->SetAromatic(false);

          FOR_BONDS_OF_MOL (bond, mol) {
            if (bondTypes[bond->GetIdx()] == Single)
              bond->SetBondOrder(1);
            if (bondTypes[bond->GetIdx()] == Double)
              bond->SetBondOrder(2);
            bond->SetAromatic(false);
          }
          mol->SetAromaticPerceived(false);

          // Call functor with tautomer
          functor(mol);
        }

      } else { // IsLeafNode
        // Go to the next unassigned atom
        EnumerateRecursive(mol, atomTypes, bondTypes,numHydrogens, functor, depth + 1);
      }

      if (m_canonical && m_foundLeafNode)
        propagation.release();
    }

    OBAtom *SelectNextAtom(const std::vector<Type> &atomTypes) const
    {
      // Select next lowest canonical unassigned atom
      for (std::size_t i = 0; i < m_canonAtoms.size(); ++i)
        if (atomTypes[m_canonAtoms[i]->GetIndex()] == Unassigned)
          return m_canonAtoms[i];
      return nullptr;
    }

    bool HasTooManyHydrogensLeft(const std::vector<Type> &atomTypes, int numHydrogens) const
    {
      // Check to ensure there are at least enough unassigned atoms left to assign hydrogens to
      return std::count(atomTypes.begin(), atomTypes.end(), Unassigned) < numHydrogens;
    }

    void EnumerateRecursive(OBMol *mol, std::vector<Type> &atomTypes, std::vector<Type> &bondTypes,
        int numHydrogens, TautomerFunctor &functor, int depth)
    {
#ifdef DEBUG
      EnterExit ee("EnumerateRecursive", depth);

      std::cout << ee.indentation() << numHydrogens << " hydrogens left [1]..." << std::endl;
      PrintAtomTypes(atomTypes, depth+1);
      PrintBondTypes(mol, bondTypes, depth+1);
#endif

      if (m_canonical && m_foundLeafNode) {
#ifdef DEBUG
        std::cout << ee.indentation() << "--> Canonical and leaf node found, backtracking..." << std::endl;
#endif
        return;
      }

      // Select next lowest canonical unassigned atom
      OBAtom *nextAtom = SelectNextAtom(atomTypes);
      if (!nextAtom) {
#ifdef DEBUG
        std::cout << ee.indentation() << "--> No unassigned atoms left, backtracking..." << std::endl;
#endif
        return;
      }

      // Assign it a hydrogen if there are hydrogens left
      if (numHydrogens) {
        AssignDonorRAII donor(atomTypes, numHydrogens, nextAtom);
#ifdef DEBUG
        std::cout << ee.indentation() << "-> Assign atom " << nextAtom->GetIndex() << " Donor" << std::endl;
#endif

        // Propagate...
        AssignmentPropagation(mol, atomTypes, bondTypes, numHydrogens, functor, depth + 1);

#ifdef DEBUG
        std::cout << ee.indentation() << numHydrogens << " hydrogens left [1b]..." << std::endl;
        PrintAtomTypes(atomTypes, depth+1);
        PrintBondTypes(mol, bondTypes, depth+1);
#endif

        if (m_canonical && m_foundLeafNode) {
          donor.release();
          return; // early termination
        }
      }

#ifdef DEBUG
      std::cout << ee.indentation() << numHydrogens << " hydrogens left [2]..." << std::endl;
      PrintAtomTypes(atomTypes, depth+1);
      PrintBondTypes(mol, bondTypes, depth+1);
#endif

      // Assign it acceptor after backtracking
      AssignAcceptorRAII acceptor(atomTypes, nextAtom);
#ifdef DEBUG
      std::cout << ee.indentation() << "-> Assign atom " << nextAtom->GetIndex() << " Acceptor" << std::endl;
#endif

      // Check to ensure there are at least enough unassigned atoms left to assign hydrogens to
      if (!HasTooManyHydrogensLeft(atomTypes, numHydrogens)) {
        // Propagate...
        AssignmentPropagation(mol, atomTypes, bondTypes, numHydrogens, functor, depth + 1);
      } else {
#ifdef DEBUG
        std::cout << ee.indentation() << "--> Too many hydrogens left, backtracking..." << std::endl;
#endif
      }
    }

#ifdef DEBUG
    static void SanityCheckHydrogens(OBAtom *atom, int indent)
    {
      OBMol* mol = atom->GetParent();
      unsigned int totH = 0;
      FOR_ATOMS_OF_MOL(atm, mol) {
        totH += atm->GetImplicitHCount();
      }
      std::cout << indentation(indent) << "Total no. of hydrogens in molecule: " << totH << std::endl;
    }
#endif

    static void DecrementImplicitHCount(OBAtom* atom, int indent)
    {
      atom->SetImplicitHCount(atom->GetImplicitHCount() - 1);
#ifdef DEBUG
      std::cout << indentation(indent) << "Decremented " << atom->GetIndex() << " (" << atom->GetAtomicNum() << ") to " << atom->GetImplicitHCount() << std::endl;
      SanityCheckHydrogens(atom, indent);
#endif
    }

    void Enumerate(OBMol *mol, TautomerFunctor &functor, bool canonical)
    {
      m_canonical = canonical;
      m_foundLeafNode = false;
      // Make hydrogens implicit
      mol->DeleteHydrogens();

      // Assign atom types
      std::vector<Type> atomTypes = InitializeAtomTypes(mol);
      // Assign bond types
      std::vector<Type> bondTypes = InitializeBondTypes(mol, atomTypes);

      // If all bonds are assigned around an atom, mark it as Other (e.g. sulfon oxygen)
      // This avoids premature backtracking (see IsPropagationValid)
      FOR_ATOMS_OF_MOL (atom, mol) {
        int numUnassignedBonds = 0;
        FOR_BONDS_OF_ATOM (bond, &*atom)
          if (bondTypes[bond->GetIdx()] == Unassigned)
            ++numUnassignedBonds;
        if (!numUnassignedBonds)
          atomTypes[atom->GetIndex()] = Other;
      }

#ifdef DEBUG
      PrintAtomTypes(atomTypes);
      PrintBondTypes(mol, bondTypes);
#endif

      // Store original h counts
      std::vector<unsigned int> hcounts;
      FOR_ATOMS_OF_MOL(atom, mol)
        hcounts.push_back(atom->GetImplicitHCount());
#ifdef DEBUG
      SanityCheckHydrogens(mol->GetAtom(1), 0);
#endif

      // Count the number of hydrogens
      // Mark donor/acceptor atoms as unassigned
      int numHydrogens = 0;
      for (std::size_t i = 0; i < atomTypes.size(); ++i) {
        if (atomTypes[i] == Donor) {
          numHydrogens++;
          OBAtom* atm = mol->GetAtom(i + 1); // off-by-one
          DecrementImplicitHCount(atm, 0);
        }
        if (atomTypes[i] == Donor || atomTypes[i] == Acceptor)
          atomTypes[i] = Unassigned;
      }
#ifdef DEBUG
      PrintAtomTypes(atomTypes);
#endif

      // Store original bond orders
      std::vector<int> bondOrders;
      FOR_BONDS_OF_MOL (bond, mol)
        bondOrders.push_back(bond->GetBondOrder());

      if (canonical) {
        // Set all unassigned bonds to single
        FOR_BONDS_OF_MOL (bond, mol)
          if (bondTypes[bond->GetIdx()] == Unassigned)
            bond->SetBondOrder(1);
        mol->SetAromaticPerceived(false);

        // Compute canonical labels
        std::vector<unsigned int> symmetryClasses;
        OBGraphSym gs(mol);
        gs.GetSymmetry(symmetryClasses);

        std::vector<unsigned int> canonLabels;
        CanonicalLabels(mol, symmetryClasses, canonLabels);

        m_canonAtoms.resize(mol->NumAtoms());
        for (std::size_t i = 0; i < mol->NumAtoms(); ++i)
          m_canonAtoms[canonLabels[i]-1] = mol->GetAtom(i+1);
      } else {
        FOR_ATOMS_OF_MOL (atom, mol)
          m_canonAtoms.push_back(&*atom);
      }

      EnumerateRecursive(mol, atomTypes, bondTypes, numHydrogens, functor, 0);

      if (!canonical || !m_foundLeafNode) {
        // Restore original bond orders & aromaticity
        FOR_BONDS_OF_MOL (bond, mol) {
          bond->SetBondOrder(bondOrders[bond->GetIdx()]);
          bond->SetAromatic(false);
        }
        // Restore original implicit H counts & aromaticity
        FOR_ATOMS_OF_MOL(atom, mol) {
          atom->SetImplicitHCount(hcounts[atom->GetIndex()]);
          atom->SetAromatic(false);
        }
      }

      mol->SetAromaticPerceived(false);

      // If no leaf was found, call functor with input molecule
      if (!canonical && !m_foundLeafNode)
        functor(mol);
    }


  }; // TuatomerImpl

  void UniqueTautomerFunctor::operator()(OBMol *mol)
  {
    OpenBabel::OBConversion conv;
    conv.SetOutFormat("can");
    std::string smiles = conv.WriteString(mol, true);

    // Make sure the found tautomer is unique
    if (std::find(m_smiles.begin(), m_smiles.end(), smiles) != m_smiles.end())
      return;

    m_smiles.push_back(smiles);
    this->operator()(mol, smiles);
  }

  void EnumerateTautomers(OBMol *mol, TautomerFunctor &functor)
  {
    TautomerImpl impl;
    impl.Enumerate(mol, functor, false);
  }

  void CanonicalTautomer(OBMol *mol)
  {
    class Functor : public TautomerFunctor {
      void operator()(OBMol *) {}
    };
    Functor functor;
    TautomerImpl impl;
    impl.Enumerate(mol, functor, true);
  }

}

//! \file tautomer.cpp
//! \brief Tautomer support.
