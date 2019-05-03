/**********************************************************************
tautomer.cpp - Tautomer support

  Copyright (C) 2011 by Tim Vandermeersch

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

namespace OpenBabel {

//#define DEBUG 1

  /**
   * http://www.daylight.com/meetings/emug99/Delany/taut_html/index.htm
   */
  struct TautomerImpl
  {
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

    bool m_canonical, m_foundLeafNode;

    // debug
    void print_atom_types(const std::vector<Type> &types)
    {
      std::cout << "Atom Types:" << std::endl;
      for(std::size_t i = 0; i < types.size(); ++i) {
        std::cout << "  " << i << ": ";

        switch (types[i]) {
          case Assigned:
            std::cout << "Assigned" << std::endl;
            break;
          case Single:
            std::cout << "Single" << std::endl;
            break;
          case Double:
            std::cout << "Double" << std::endl;
            break;
          case Donor:
            std::cout << "Donor" << std::endl;
            break;
          case Acceptor:
            std::cout << "Acceptor" << std::endl;
            break;
          case Hybridized:
            std::cout << "Hybridized" << std::endl;
            break;
          case Other:
            std::cout << "Other" << std::endl;
            break;
          case Unassigned:
            std::cout << "Unassigned" << std::endl;
            break;
        }
      }
    }

    // debug
    void print_bond_types(const std::vector<Type> &types)
    {
      std::cout << "Bond Types:" << std::endl;
      for(std::size_t i = 0; i < types.size(); ++i) {
        std::cout << "  " << i << ": ";

        switch (types[i]) {
          case Assigned:
            std::cout << "Assigned" << std::endl;
            break;
          case Single:
            std::cout << "Single" << std::endl;
            break;
          case Double:
            std::cout << "Double" << std::endl;
            break;
          case Donor:
            std::cout << "Donor" << std::endl;
            break;
          case Acceptor:
            std::cout << "Acceptor" << std::endl;
            break;
          case Hybridized:
            std::cout << "Hybridized" << std::endl;
            break;
          case Other:
            std::cout << "Other" << std::endl;
            break;
          case Unassigned:
            std::cout << "Unassigned" << std::endl;
            break;
        }
      }
    }


    void AssignAtomTypes(OBMol *mol, std::vector<Type> &types)
    {
      FOR_ATOMS_OF_MOL (atom, mol) {
        switch (atom->GetAtomicNum()) {
          case Carbon:
            if (atom->HasDoubleBond())
              types.push_back(Hybridized);
            else
              types.push_back(Other);
            break;

          case Nitrogen:
            if (atom->HasDoubleBond())
              types.push_back(Acceptor);
            else if (atom->GetImplicitHCount())
              types.push_back(Donor);
            else
              types.push_back(Other);
            break;

          case Oxygen:
          case Sulfur:
          case Selenium:
          case Tellurium:
            if (atom->HasDoubleBond())
              types.push_back(Acceptor);
            else if (atom->GetImplicitHCount())
              types.push_back(Donor);
            else
              types.push_back(Other);
            break;

          default:
            types.push_back(Other);        
        }
      }

      assert( types.size() == mol->NumAtoms() );
    }

    void AssignBondTypes(OBMol *mol, const std::vector<Type> &atomTypes, std::vector<Type> &bondTypes)
    {
      FOR_BONDS_OF_MOL (bond, mol) {
        assert( bond->GetBeginAtomIdx() <= mol->NumAtoms() );
        assert( bond->GetEndAtomIdx() <= mol->NumAtoms() );
        if (atomTypes[bond->GetBeginAtomIdx()-1] != Other && atomTypes[bond->GetEndAtomIdx()-1] != Other)
          bondTypes.push_back(Unassigned);
        else 
          bondTypes.push_back(Assigned);
      }
    }

    struct Level 
    {
      OBAtom *assigned;
      std::vector<OBAtom*> propagatedAtoms;
      std::vector<OBBond*> propagatedBonds;
    };

    static void SanityCheckHydrogens(OBAtom *atom)
    {
      OBMol* mol = atom->GetParent();
      unsigned int totH = 0;
      FOR_ATOMS_OF_MOL(atm, mol) {
        totH += atm->GetImplicitHCount();
      }
      printf("Total no. of hydrogens in molecule: %d\n", totH);
    }

    static void DecrementImplicitHCount(OBAtom* atom)
    {
      atom->SetImplicitHCount(atom->GetImplicitHCount() - 1);
#ifdef DEBUG
      printf("Decremented %d (%d) to %d\n", atom->GetIndex(), atom->GetAtomicNum(), atom->GetImplicitHCount());
      SanityCheckHydrogens(atom);
#endif
    }

    static void IncrementImplicitHCount(OBAtom* atom)
    {
      atom->SetImplicitHCount(atom->GetImplicitHCount() + 1);
#ifdef DEBUG
      printf("Incremented %d (%d) to %d\n", atom->GetIndex(), atom->GetAtomicNum(), atom->GetImplicitHCount());
      SanityCheckHydrogens(atom);
#endif
    }

    void print_assigned(const std::vector<Level> &levels, const std::vector<Type> &atomTypes)
    {
      for (std::size_t i = 0; i < levels.size(); ++i)
        std::cout << levels[i].assigned->GetIndex() << " ";
      std::cout << std::endl;
      for (std::size_t i = 0; i < levels.size(); ++i)
        if (atomTypes[levels[i].assigned->GetIndex()] == Donor)
          std::cout << "D ";
        else
          std::cout << "A ";
      std::cout << std::endl;
    }

    void Backtrack(std::vector<Type> &atomTypes, std::vector<Type> &bondTypes, std::vector<Level> &levels, int &numHydrogens)
    {
      OBAtom *lastAtom = levels.back().assigned;
#ifdef DEBUG
      std::cout << "  Backtrack... " << lastAtom->GetIndex() << std::endl;
#endif
      if (atomTypes[lastAtom->GetIndex()] == Donor)
        DecrementImplicitHCount(lastAtom); // Noel: Why no need for accompanying numHydrogens++?
      atomTypes[lastAtom->GetIndex()] = Unassigned;

      std::vector<OBAtom*> &propagatedAtoms = levels.back().propagatedAtoms;
      for (std::size_t i = 0; i < propagatedAtoms.size(); ++i) {
        if (atomTypes[propagatedAtoms[i]->GetIndex()] == Donor) {
          DecrementImplicitHCount(propagatedAtoms[i]);
          numHydrogens++;
        }
        atomTypes[propagatedAtoms[i]->GetIndex()] = Unassigned;
      }

      std::vector<OBBond*> &propagatedBonds = levels.back().propagatedBonds;
      for (std::size_t i = 0; i < propagatedBonds.size(); ++i)
        bondTypes[propagatedBonds[i]->GetIdx()] = Unassigned;

      levels.pop_back();
    }

    void AssignmentPropagation(OBMol *mol, std::vector<Type> &atomTypes, std::vector<Type> &bondTypes, const std::vector<OBAtom*> &canonAtoms,
        int &numHydrogens, std::vector<Level> &levels, TautomerFunctor &functor)
    {
      bool changed = true;
      while (changed) {
        changed = false;

        // Assigning an atom Donor causes all of its bonds to be assigned single
        FOR_ATOMS_OF_MOL (atom, mol) {
          if (atomTypes[atom->GetIndex()] == Donor) {
            FOR_BONDS_OF_ATOM (bond, &*atom) {
              if (bondTypes[bond->GetIdx()] == Unassigned) {
                bondTypes[bond->GetIdx()] = Single;
                levels.back().propagatedBonds.push_back(&*bond);
                changed = true;
#ifdef DEBUG
                std::cout << "    -> Rule 1: Assign " << bond->GetBeginAtomIdx()-1 << "-" << bond->GetEndAtomIdx()-1 << " Single" << std::endl;
#endif
              }
            }
          }
        }

        // Any unassigned atom that has all of its bonds assigned single, must be assigned a donor
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
            atomTypes[atom->GetIndex()] = Donor; 
            IncrementImplicitHCount(&*atom);
            numHydrogens--;
            levels.back().propagatedAtoms.push_back(&*atom);
            changed = true;
#ifdef DEBUG
            std::cout << "    -> Rule 2: Assign " << atom->GetIndex() << " Donor" << std::endl;
#endif
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
            atomTypes[atom->GetIndex()] = Acceptor;
            levels.back().propagatedAtoms.push_back(&*atom);
            changed = true;
#ifdef DEBUG
            std::cout << "    -> Rule 3: Assign " << atom->GetIndex() << " Acceptor" << std::endl;
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
                bondTypes[bond->GetIdx()] = Single;
                levels.back().propagatedBonds.push_back(&*bond);
                changed = true;
#ifdef DEBUG
                std::cout << "    -> Rule 4: Assign " << bond->GetBeginAtomIdx()-1 << "-" << bond->GetEndAtomIdx()-1 << " Single" << std::endl;
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
                bondTypes[bond->GetIdx()] = Double;
                levels.back().propagatedBonds.push_back(&*bond);
                changed = true;
#ifdef DEBUG
                std::cout << "    -> Rule 5: Assign " << bond->GetBeginAtomIdx()-1 << "-" << bond->GetEndAtomIdx()-1 << " Double" << std::endl;
#endif
              }
            }
          }
        }

      } // while (changed)

      // An invalid atom type causes the search to backtrack
      bool invalid = false;
      FOR_ATOMS_OF_MOL (atom, mol) {
          int numSingleBond = 0, numDoubleBond = 0, numUnassigned = 0;
          FOR_BONDS_OF_ATOM (bond, &*atom) {
            if (bondTypes[bond->GetIdx()] == Single)
              numSingleBond++;
            if (bondTypes[bond->GetIdx()] == Double)
              numDoubleBond++;
            if (bondTypes[bond->GetIdx()] == Unassigned)
              numUnassigned++;
          }


        if (atomTypes[atom->GetIndex()] == Donor)
          if (numDoubleBond) {
            invalid = true;
#ifdef DEBUG
            std::cout << "invalid Donor" << std::endl;
#endif
          }

        if (atomTypes[atom->GetIndex()] == Acceptor || atomTypes[atom->GetIndex()] == Hybridized) {
          if (!numUnassigned && !numDoubleBond) {
            invalid = true;
#ifdef DEBUG
            std::cout << "invalid Acceptor/Hybridized 1" << std::endl;
#endif
          }
          if (numDoubleBond > 1) {
            invalid = true;
#ifdef DEBUG
            std::cout << "invalid Acceptor/Hybridized 2" << std::endl;
#endif
          }
        }
      }

      // Check to see if we are at a leaf node (i.e. no unassigned atoms left)
      bool leafNode = true;
      for (std::size_t i = 0; i < atomTypes.size(); ++i)
        if (atomTypes[i] == Unassigned) {
          leafNode = false;
          break;
        }

      if (leafNode) {
        // Check to make sure there are no unassigned bonds remaining
        // This happens for phenyl rings...
        bool unassignedBonds = false;
        for (std::size_t i = 0; i < bondTypes.size(); ++i)
        if (bondTypes[i] == Unassigned) {
          unassignedBonds = true;
          // Randomly assign a single bond, the other bonds will be propagated
          bondTypes[i] = Single;
          break;
        }

        if (unassignedBonds) {
          AssignmentPropagation(mol, atomTypes, bondTypes, canonAtoms, numHydrogens, levels, functor);
          return;
        }

      }

      if (!invalid) {
        if (leafNode) {
          FOR_BONDS_OF_MOL (bond, mol) {
            if (bondTypes[bond->GetIdx()] == Single)
              bond->SetBondOrder(1);
            if (bondTypes[bond->GetIdx()] == Double)
              bond->SetBondOrder(2);
          }
          if (!numHydrogens) {
#ifdef DEBUG
            std::cout << "  --> LeafNode reached..." << std::endl;
#endif
            m_foundLeafNode = true;
            
            mol->BeginModify();
            mol->EndModify();
            functor(mol);
          }
        } else
          // Go to the next unassigned atom
          EnumerateRecursive(mol, atomTypes, bondTypes, canonAtoms, numHydrogens, levels, functor);
      }

    }

    void EnumerateRecursive(OBMol *mol, std::vector<Type> &atomTypes, std::vector<Type> &bondTypes, const std::vector<OBAtom*> &canonAtoms,
        int numHydrogens, std::vector<Level> &levels, TautomerFunctor &functor)
    {
      if (m_canonical && m_foundLeafNode)
        return;
#ifdef DEBUG
      std::cout << "EnumerateRecursive" << std::endl;
#endif
      // Select next lowest canonical unassigned atom
      for (std::size_t i = 0; i < canonAtoms.size(); ++i)
        if (atomTypes[canonAtoms[i]->GetIndex()] == Unassigned) {
          if (numHydrogens) {
            // Assign it a hydrogen if there are hydrogens left
            atomTypes[canonAtoms[i]->GetIndex()] = Donor;
            IncrementImplicitHCount(canonAtoms[i]);
            numHydrogens--;
#ifdef DEBUG
            std::cout << "  Assigned " << canonAtoms[i]->GetIndex() << " Donor" << std::endl;
#endif
          } else {
            // Assign it acceptor otherwise
            atomTypes[canonAtoms[i]->GetIndex()] = Acceptor;
#ifdef DEBUG
            std::cout << "  Assigned " << canonAtoms[i]->GetIndex() << " Acceptor" << std::endl;
#endif
          }

          // store the atom for backtracking later
          levels.push_back(Level());
          levels.back().assigned = canonAtoms[i];
          break;
        }


      // Check to ensure there are at least enough unassigned atoms left to assign hydrogens to
      int numUnassigned = 0;
      for (std::size_t i = 0; i < canonAtoms.size(); ++i)
        if (atomTypes[canonAtoms[i]->GetIndex()] == Unassigned)
          numUnassigned++;
      if (numUnassigned < numHydrogens) {
        // Backtrack
#ifdef DEBUG
        std::cout << "  Backtrack..." << std::endl;
#endif
        Backtrack(atomTypes, bondTypes, levels, numHydrogens);
        return;
      }
 

      if (levels.size()) {
        AssignmentPropagation(mol, atomTypes, bondTypes, canonAtoms, numHydrogens, levels, functor);
        if (m_canonical && m_foundLeafNode)
          return; // early termination

#ifdef DEBUG
        std::cout << "Change?" << std::endl;
        print_assigned(levels, atomTypes);
#endif
        
        OBAtom *lastAtom = levels.back().assigned;
        if (atomTypes[lastAtom->GetIndex()] == Donor) {
          // Change atom to acceptor and continue the search
          atomTypes[lastAtom->GetIndex()] = Acceptor;
          DecrementImplicitHCount(lastAtom);
          numHydrogens++;
#ifdef DEBUG
          std::cout << "  Change " << lastAtom->GetIndex() << " to Acceptor" << std::endl;
#endif

          std::vector<OBAtom*> &propagatedAtoms = levels.back().propagatedAtoms;
          for (std::size_t i = 0; i < propagatedAtoms.size(); ++i) {
            if (atomTypes[propagatedAtoms[i]->GetIndex()] == Donor) {
              DecrementImplicitHCount(propagatedAtoms[i]);
              numHydrogens++;
            }
            atomTypes[propagatedAtoms[i]->GetIndex()] = Unassigned;
          }

          std::vector<OBBond*> &propagatedBonds = levels.back().propagatedBonds;
          for (std::size_t i = 0; i < propagatedBonds.size(); ++i)
            bondTypes[propagatedBonds[i]->GetIdx()] = Unassigned;


          AssignmentPropagation(mol, atomTypes, bondTypes, canonAtoms, numHydrogens, levels, functor);
          if (m_canonical && m_foundLeafNode)
            return; // early termination
        } 

      // Backtrack
#ifdef DEBUG
        print_assigned(levels, atomTypes);
#endif
        Backtrack(atomTypes, bondTypes, levels, numHydrogens);
      }
    }


    void Enumerate(OBMol *mol, TautomerFunctor &functor, bool canonical)
    {
      m_canonical = canonical;
      m_foundLeafNode = false;
      // Make hydrogens implicit
      mol->DeleteHydrogens();

      // Assign atom types
      std::vector<Type> atomTypes;
      AssignAtomTypes(mol, atomTypes);
#ifdef DEBUG
      print_atom_types(atomTypes);
#endif

      // Assign bond types
      std::vector<Type> bondTypes;
      AssignBondTypes(mol, atomTypes, bondTypes);
#ifdef DEBUG
      print_bond_types(bondTypes);
#endif
      // Store original h counts
      std::vector<unsigned int> hcounts;
      FOR_ATOMS_OF_MOL(atom, mol)
        hcounts.push_back(atom->GetImplicitHCount());
#ifdef DEBUG
      SanityCheckHydrogens(mol->GetAtom(1));
#endif

      // Count the number of hydrogens
      // Mark donor/acceptor atoms as unassigned
      int numHydrogens = 0;
      for (std::size_t i = 0; i < atomTypes.size(); ++i) {
        if (atomTypes[i] == Donor) {
          numHydrogens++;
          OBAtom* atm = mol->GetAtom(i + 1); // off-by-one
          DecrementImplicitHCount(atm);
        }
        if (atomTypes[i] == Donor || atomTypes[i] == Acceptor)
          atomTypes[i] = Unassigned;
      }
#ifdef DEBUG
      print_atom_types(atomTypes);
#endif

      // Store original bond orders
      std::vector<int> bondOrders;
      FOR_BONDS_OF_MOL (bond, mol)
        bondOrders.push_back(bond->GetBondOrder());

      // FIXME
      std::vector<unsigned int> canonLabels;
      if (canonical) {
        // Set all unassigned bonds to single
        FOR_BONDS_OF_MOL (bond, mol)
          if (bondTypes[bond->GetIdx()] == Unassigned)
            bond->SetBondOrder(1);
        mol->SetAromaticPerceived(false);
        
        std::vector<unsigned int> symmetryClasses;
        OBGraphSym gs(mol);
        gs.GetSymmetry(symmetryClasses);
        CanonicalLabels(mol, symmetryClasses, canonLabels);
      } else {
        for (std::size_t i = 0; i < mol->NumAtoms(); ++i)
          canonLabels.push_back(i+1);
      }

      std::vector<OBAtom*> canonAtoms(mol->NumAtoms());
      for (std::size_t i = 0; i < mol->NumAtoms(); ++i)
        canonAtoms[canonLabels[i]-1] = mol->GetAtom(i+1);
      
      std::vector<Level> levels;
      EnumerateRecursive(mol, atomTypes, bondTypes, canonAtoms, numHydrogens, levels, functor);

      if (!canonical) {
        // Restore original bond orders
        FOR_BONDS_OF_MOL (bond, mol)
          bond->SetBondOrder(bondOrders[bond->GetIdx()]);
        // Restore original implicit H counts
        FOR_ATOMS_OF_MOL(atom, mol)
          atom->SetImplicitHCount(hcounts[atom->GetIndex()]);
      }

    }

  
  }; // TuatomerImpl

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
