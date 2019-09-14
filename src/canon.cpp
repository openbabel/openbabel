/**********************************************************************
canon.h - Canonical labeling.

  Copyright (C) 2009-2010 by Tim Vandermeersch
  Copyright (C) 2005-2006, eMolecules, Inc. (www.emolecules.com)
                           Craig A. James

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

#include <openbabel/canon.h>
#include <openbabel/graphsym.h>
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/obutil.h>
#include <openbabel/elements.h>

#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>

#include <iterator> // std::istream_iterator
#include <cassert>
#include <algorithm>
#include <cmath>
#include <limits>

#include "stereo/stereoutil.h"

#define DEBUG 0

#define MAX_IDENTITY_NODES 50

using namespace std;

// debug function
template<typename T>
void print_vector(const std::string &label, const std::vector<T> &v)
{
  std::cout << label << ": ";
  for (std::size_t i = 0; i < v.size(); ++i)
    if (v[i] < 10)
      std::cout << " " << v[i] << " ";
    else
      std::cout << v[i] << " ";

  std::cout << endl;
}


namespace OpenBabel {

  static unsigned int TotalHydrogenCount(OBAtom* atom)
  {
    return atom->ExplicitHydrogenCount() + atom->GetImplicitHCount();
  }

  inline bool CompareBondPairSecond(const std::pair<OBBond*,unsigned int> &a,const std::pair<OBBond*,unsigned int> &b)
  {
    return(a.second < b.second);
  }

  /**
   * Helper function for getFragment below.
   */
  void addNbrs(OBBitVec &fragment, OBAtom *atom, const OBBitVec &mask, const std::vector<OBBond*> &metalloceneBonds)
  {
    FOR_NBORS_OF_ATOM (nbr, atom) {
      if (!mask.BitIsSet(nbr->GetIdx()))
        continue;
      // skip visited atoms
      if (fragment.BitIsSet(nbr->GetIdx()))
        continue;
      // skip mettalocene bonds
      if (std::find(metalloceneBonds.begin(), metalloceneBonds.end(),
            atom->GetParent()->GetBond(atom, &*nbr)) != metalloceneBonds.end())
        continue;
      // add the neighbor atom to the fragment
      fragment.SetBitOn(nbr->GetIdx());
      // recurse...
      addNbrs(fragment, &*nbr, mask, metalloceneBonds);
    }
  }

  /**
   * Create an OBBitVec objects with bets set for the fragment consisting of all
   * atoms for which there is a path to atom without going through skip. These
   * fragment bitvecs are indexed by atom idx (i.e. OBAtom::GetIdx()).
   */
  OBBitVec getFragment(OBAtom *atom, const OBBitVec &mask, const std::vector<OBBond*> &metalloceneBonds = std::vector<OBBond*>())
  {
    OBBitVec fragment;
    fragment.SetBitOn(atom->GetIdx());
    // start the recursion
    addNbrs(fragment, atom, mask, metalloceneBonds);
    return fragment;
  }

  OBBitVec getFragment(OBAtom *atom, OBAtom *skip, const OBBitVec &mask)
  {
    struct getFragmentImpl
    {
      static void addNbrs(OBBitVec &fragment, OBAtom *atom, OBAtom *skip, const OBBitVec &mask)
      {
        FOR_NBORS_OF_ATOM (nbr, atom) {
          // don't pass through skip
          if (nbr->GetIdx() == skip->GetIdx())
            continue;
          // skip visited atoms
          if (fragment.BitIsSet(nbr->GetIdx()))
            continue;
          if (!mask.BitIsSet(nbr->GetIdx()))
            continue;
          // add the neighbor atom to the fragment
          fragment.SetBitOn(nbr->GetIdx());
          // recurse...
          addNbrs(fragment, &*nbr, skip, mask);
        }
      }
    };

    OBBitVec fragment;
    fragment.SetBitOn(atom->GetIdx());
    // start the recursion
    getFragmentImpl::addNbrs(fragment, atom, skip, mask);
    return fragment;
  }



  bool isFerroceneBond(OBBond *bond)
  {
    if (bond->GetBondOrder() != 1)
      return false;

    OBAtom *Fe = 0, *C = 0;

    OBAtom *begin = bond->GetBeginAtom();
    if (begin->GetAtomicNum() == 26)
      Fe = begin;
    if (begin->GetAtomicNum() == 6)
      C = begin;

    OBAtom *end = bond->GetEndAtom();
    if (end->GetAtomicNum() == 26)
      Fe = end;
    if (end->GetAtomicNum() == 6)
      C = end;

    if (!Fe || !C)
      return false;

    if (Fe->GetExplicitDegree() < 10)
      return false;

    return C->HasDoubleBond() && C->IsInRing();
  }

  void findMetalloceneBonds(std::vector<OBBond*> &bonds, OBMol *mol, const std::vector<unsigned int> &symmetry_classes)
  {
    FOR_ATOMS_OF_MOL (atom, mol) {
      if (!atom->IsInRingSize(3))
        continue;
      std::vector<unsigned int> nbrSymClasses;
      FOR_NBORS_OF_ATOM (nbr, &*atom) {
        if (nbr->IsInRingSize(3))
          nbrSymClasses.push_back(symmetry_classes[nbr->GetIndex()]);
      }

      if (nbrSymClasses.size() < 8)
        continue;

      std::sort(nbrSymClasses.begin(), nbrSymClasses.end());
      unsigned int numUnique = std::unique(nbrSymClasses.begin(), nbrSymClasses.end()) - nbrSymClasses.begin();
      if (numUnique > 1)
        continue;

      FOR_NBORS_OF_ATOM (nbr, &*atom)
        bonds.push_back(mol->GetBond(&*atom, &*nbr));
    }
  }


  /**
   * @page canonical_code_algorithm Canonical Coding Algorithm
   *
   * @section canonical_introduction Introduction
   * The aim of a canonical coding algorithm is to assign unique labels to
   * atoms regardless of their input order. In practise, these labels or orders
   * are atom indexes resulting from the order in which the atoms are defined
   * in an input file. Although most chemical file formats could be used to
   * store canonical ordered molecules, canonical single line notations are used
   * more often since they allow two canonical molecules to be compared using a
   * string comparison. Two well known examples are canonical smiles and inchi.
   * While the canonical smiles for the same molecule should always be the same
   * when using a specific implementation (i.e. toolkits), these canonical
   * smiles are usually not transferable between implementations. There is only
   * one inchi implementation and the canonical code, although not formally
   * specified is always the same. A typical use case for canonical codes is to
   * determine if a molecule is already in a database.
   *
   * The rest of this page documents the OpenBabel canonical coding algorithm.
   *
   * @section Topological Symmetry
   * The starting point of the canonical coding algorithm is the topological
   * symmetry or graph symmetry. This symmetry can be expressed by assigning
   * symmetry classes to atoms. Symmetric atoms have the same symmetry class.
   * The term color is often used in graph theory as a synonym for symmetry
   * classes.
   *
   * Symmetry classes are used as graph invariants to order atoms. The symmetry
   * classes can be used in two ways:
   *
   * - lowest symmetry class = lowest label (i.e. 1)
   * - highest symmetry class = lowest label
   *
   * Both methods are valid. OpenBabel uses the second method. In the following
   * sections there are more of these decisions that have to be made. The first,
   * choice is used as example to illustrate that it doesn't make any real
   * difference as long as the same choice is used every time. When trying to
   * reproduce OpenBabel canonical smiles the same choices or conventions have
   * to be used.
   *
   * The used initial graph invariant is specified below and the iteration
   * algorithm is similar to the original Morgan algorithm. The OBGraphSym class
   * should be consulted for a detailed description.
     @verbatim
          0                   1                   2                   3
          1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2
     GI = D D D D D D D D D D D V V V V A R N N N N N N N B B B B C C C C

     GI: Graph invariant (32 bits)
     D: Graph theoretical distance (10 bits)
     V: Heavy valence (4 bits)
     A: Aromaticity (1 bit)
     R: Ring atom (1 bit)
     N: Atomic number (7 bits)
     B: Heavy bond sum (4 bits)
     C: Formal charge + 7 (4 bits)
     @endverbatim
   *
   * (FIXME: ugly smiles)
   *
   * The symmetry classes encode various attributes for each atom. Consult the
   * OBGraphSym documentation for details. (note: when trying to reproduce
   * OpenBabel canonical smiles, make sure these graph invariants are the same)
   *
   * @section canonical_labeling Labeling
   * The labeling starts by selecting the initial atom using a set of rules to
   * ensure it starts at terminal atoms. (FIXME) Label 1 is assigned to this atom
   * and the labeling continues using this atom as the current atom:
   *
   * -# Find the neighbor atoms of the current atom to label next. Already labeled
   *   and atoms not in the fragment are ignored. If there are no such atoms, go
   *   to 3.
   * -# Sort the neighbors by their symmetry classes. Each group of atoms with the
   *    same symmetry class is now handled in sequence. If there is only one atom
   *    in the group, this atom is labeled and the process continues with the next
   *    group. If there are multiple (n) atoms in a group, all n! permutations of
   *    the group are tried. (note: The Optimization section below mentions some
   *    exceptions which affect the final result. These optimizations are part of
   *    the algorithm)
   * -# When all neighbor atoms of the current atom are labeled, the next atom
   *    (i.e. the atom with current_label+1 label) is considered the current atom
   *    and step 1 through 3 are repeated. If there is no next atom, the labeling
   *    is done.
   *
   * The algorithm just described results in a set of labellings that are candidates
   * for the canonical labeling. For each labeling, a code is generated to select
   * the canonical code and the associated labels.
   *
   * @section canonical_code Canonical Code
   * A canonical code can be anything encoding the molecular graph with all the
   * properties that should be canonicalized. The code should be invariant to the
   * original atom order. At this point, it may seem that this is similar to the
   * graph invariants from above (i.e. symmetry classes). However, generating these
   * codes allows for additional attributes to be included. A simple example are
   * isotopes. More importantly, stereochemistry is included in these codes which
   * is an attribute that does not depend an isolated atom.
   *
   * The actual code is a list of numbers (unsigned int) constructed by joining
   * smaller codes. These smaller codes usually encode specific attributes.
   *
   * @subsection canonical_from FROM code
   * When labeling the molecule as described earlier, all but the initial atom in
   * a connected fragment have a "from" atom (i.e. current atom in the algorithm).
   * The FROM code contains the label of the "from" atom, for each atom (ignoring
   * the initial atom) ordered by their assigned label. The FROM code is a
   * spanning tree containing all bonds that are not ring closures.
   *
   * @subsection canonical_closure CLOSURE code
   * The ring closure bonds are easily identified while creating the FROM code.
   * The CLOSURE code contains two labels for each closure bond. The labels are
   * sorted in the following way:
   *
   * - [1 6] < [6 1] (i.e. 1 6 is used, individual bond sorting)
   * - [1 3][1 4] < [1 4][1 3] (i.e. 1 3 1 4 is used, sorting bonds)
   * - [1 3][5 7] < [5 7][1 3]
   *
   * The square brackets are only to indicate individual bonds, the code is just
   * a single list.
   *
   * @subsection canonical_atomtype ATOM-TYPES & ISOTOPES code
   * The ATOM-TYPES code is simply the atomic number for each atom ordered by
   * labels. If the molecule contains isotopes, an additional ISOTOPES code is
   * used to correctly canonicalize these molecules.
   *
   * @subsection canonical_bondtype BOND-TYPES code
   * The BOND-TYPES code contains, for each bond, the bond order or 5 if the
   * bond is aromatic. The bonds are ordered using the order obtained during
   * the FROM code construction. The closure bonds are last, ordered as
   * described in the CLOSURE code section.
   *
   * @subsection canonical_stereo STEREO code
   * The STEREO code contains a descriptor for each stereo center. The stereo
   * centers are ordered using the central atom(s) (i.e. center for tetrahedral
   * stereochemistry, lowest label for cistrans).
   *
   * @subsection canonical_canonical Canonical Code
   * The canonical code is the greatest code (elements are compared, the first
   * difference determines order). The complete canonical code can also be seen
   * as layers. The first layer is the FROM code together with the CLOSURE code.
   * Together these codes describe the graph without node (atom) or edge (bond)
   * properties. The ATOM-TYPES code is the layer adding node properties. The
   * optional ISOTOPES code also encodes node properties. Edge properties are
   * specified in the BOND-TYPES code. The previous codes together are enough
   * to specify the constitution of a molecule. Additional layers such as the
   * STEREO code could be ignored for certain use cases. For example, if two
   * molecules have the same canonical code up to the STEREO part, they are
   * stereoisomers.
   *
   * @section canonical_fragments Disconnected Fragments
   * Disconnected fragments are handled individually and the resulting canonical
   * codes are sorted. The final canonical code for the whole molecule, are the
   * sorted canonical codes joined together.
   *
   * @section canonical_optimization Optimizations
   * This section documents optimizations to reduce the number of states which
   * affect the final result. The aim of these optimizations is to avoid the n!
   * permutations in step 2 of the labeling algorithm. Some of the optimizations
   * affect the final result and should be implemented when trying to reproduce
   * the same canonical code.
   *
   * @subsection canonical_opt1 Optimization 1
   * If the current atom is not in a ring, it follows that the bonds to all it's
   * neighbors are not ring bonds. It also follows that there is no path between
   * two neighbor atoms with the same symmetry class without going through the
   * current atom. Therefore, these ligands (i.e. the fragments that would be
   * created by deleting the current atom) can be labeled individually. For
   * example, for an atom with 4 neighbors with the same symmetry class, the
   * number of permutations to label is reduced from 4! * 4 = 96 to 4. The times
   * 4 is because after generating the 24 permutations, the ligands are not
   * labeled yet.
   *
   * @note This optimization <b>affects the final result</b>.
   *
   * @subsection canonical_opt2 Optimization 2 (not implemented)
   * If the bond(s) to the neighbor atom(s) with the same symmetry classes is not
   * a ring bond, optimization 1 could be used. This is the more general version
   * of optimization 1.
   *
   * @note This optimization <b>affects the final result</b> but is <b>not
   * implemented</b>.
   *
   * @subsection canonical_opt3 Optimization 3 (not implemented)
   * If the current atom is a stereo center, only permutations with the highest
   * descriptor produce a greatest canonical code. This does not change the final
   * result and can be implemented without changing canonical smiles etc.
   *
   * @subsection canonical_opt4 Other optimizations
   * In essence, the canonical coding algorithm is very similar to the algorithm
   * from the well-known nauty package [1, 2]. However, the optimizations mentioned
   * above are far easier to implement, document and result in acceptable
   * performance for most structures. This enables anyone to reproduce OpenBabel
   * canonical smiles without the absolute requirement to implement the complex
   * optimizations from nauty. Since the optimizations from the nauty paper do
   * not have any affect on the final result, these can be implemented as required.
   *
   * The OpenBabel algorithm uses the molecule graph as the main object while
   * searching for a canonical code. Nauty starts from a partition of the set
   * of nodes which is refined. This is not the essential part though. The
   * similarity between the algorithms is that both algorithms share the
   * concept of a search tree. A completed canonical candidate code is a leaf
   * node in this search tree. If a second labeling is found with the same
   * canonical code, there exists a permutation mapping the first labels to
   * the second. This permutation is an automorphism that takes all attributes
   * of the canonical code into account. The found automorphism is not only
   * an automorphism of the molecule but also an isomorphism of the search tree.
   * These accidentally discovered automorphisms can be used to reduce the
   * number of states to consider in various ways.
   *
   * @subsubsection canonical_opt5 Optimization 4
   *
   *
   *
   *
   *
     @verbatim
     [1] Brendan D. McKay, Backtrack programming and isomorphism rejection on
         ordered subsets, ARS Combinatoria, Vol. 5, 1979, pp. 65-99
         http://cs.anu.edu.au/~bdm/papers/backtrack1978.pdf
     [2] Brendan D. McKay, Practical graph isomorphism, Congressus Numerantium,
         Vol. 30, 1981, pp. 45-87
         http://cs.anu.edu.au/~bdm/papers/pgi.pdf
     @endverbatim
   */

  struct CanonicalLabelsImpl
  {
    typedef std::vector<OBAtom*> Orbit;
    typedef std::vector<Orbit> Orbits;

    static void print_orbits(const Orbits &orbits)
    {
      for (std::size_t j = 0; j < orbits.size(); ++j) {
        cout << "( ";
        for (std::size_t k = 0; k < orbits[j].size(); ++k)
          cout << orbits[j][k]->GetIndex() << " ";
        cout << ") ";
      }
    }
    static void print_orbits(const std::string &label, const Orbits &orbits)
    {
      cout << label << ": ";
      print_orbits(orbits);
      cout << endl;
    }


    /**
     * Structure to represent a completed labeling and it's associated canonical
     * candidate code.
     */
    struct FullCode
    {
      /**
       * The atom labels starting from 1. All excluded atoms have label 0.
       */
      std::vector<unsigned int> labels;
      /**
       * The canonical candidate code resulting from the @p labels.
       */
      std::vector<unsigned short> code;

      /**
       * Default constructor, used to create initial empty bestCode object.
       */
      FullCode()
      {
      }

      /**
       * Constructor specifying both the @p labels and @p from (the FROM part
       * of the canonical code). The other parts of the canonical code are
       * added using the CompleteCode function.
       */
      FullCode(const std::vector<unsigned int> &_labels, const std::vector<unsigned short> &from)
          : labels(_labels), code(from)
      {
      }

      /**
       * Compare this object's canonical candidate code to the other object's
       * code and return true if this code is greater. This method is used to
       * select the final canonical code. In other words, the canonical labels
       * for a molecule are the labels with the greatest canonical code.
       */
      inline bool operator>(const FullCode &other) const
      {
        return (code > other.code);
      }
    };

    /**
     * Structure to represent a partial labeling and the associated FROM part
     * of the canonical code.
     */
    struct PartialCode
    {
      /**
       * The atoms in canonical order. An atom is added when a label is assigned
       * to it. When all atoms are labeled, this vector will contain all atoms
       * in the fragment.
       */
      std::vector<OBAtom*> atoms;
      /**
       * The bonds in canonical order. Since the atoms are labeled from another
       * atom ("from"), these are the bonds connecting the "from" atom with the
       * next labeled atom. The FROM list forms a spanning tree that does not
       * include the ring closure bonds. These ring closure bonds are added at
       * the end of this vector, in the correct order, in the CompleteCode
       * function.
       */
      std::vector<OBBond*> bonds;
      /**
       * The FROM part of the canonical code. This is a vector containing the
       * label of the "from" atom for the atoms in canonical order. The initial
       * atom of a bonded fragment is not included in this list since there is
       * no "from" atom.
       */
      std::vector<unsigned short> from;
      /**
       * The atom labels starting from 1. All atoms start with label 0.
       */
      std::vector<unsigned int> labels;

      /**
       * Constructor specifying the number of atoms in the molecule. The @p
       * numAtoms are used to initialize the labels vector to contain label
       * 0 for all atoms.
       */
      PartialCode(std::size_t numAtoms)
      {
        labels.resize(numAtoms, 0);
      }

      /**
       * Add an atom to atoms. Used to add the initial atom.
       */
      inline void add(OBAtom *atom) { atoms.push_back(atom); }
      /**
       * Add a bond to bonds. Used to add the ring closure bonds.
       */
      inline void add(OBBond *bond) { bonds.push_back(bond); }
      /**
       * Add an atom that is labeled "from" another atom. This function adds
       * the @p atom to atoms, adds bond between @p fromAtom and @p atom to
       * bonds and updates the @p from vector.
       */
      inline void add(OBAtom *fromAtom, OBAtom *atom)
      {
        from.push_back(labels[fromAtom->GetIndex()]);
        atoms.push_back(atom);
        bonds.push_back(atom->GetParent()->GetBond(fromAtom, atom));
      }

      /**
       * Compare the @p from vector with a FullCode's code and return true if
       * this code is lower. This function is used to avoid completing a
       * labeling that will never result in a greatest canonical code.
       */
      inline bool operator<(const FullCode &other) const
      {
        std::size_t numFrom = std::min(from.size(), other.code.size());
        for (std::size_t i = 0; i < numFrom; ++i) {
          if (from[i] > other.code[i])
            return false;
          if (from[i] < other.code[i])
            return true;
        }

        return false;
      }
    };

    /**
     * Structure to store and compute stereo center parities. Currently for
     * both tetrahedral and cistrans stereochemistry. Precomputing these
     * StereoCenter objects is done for performance reasons since parities
     * need to be computed for each canonical candidate labeling.
     */
    struct StereoCenter
    {
      /**
       * Use the stored nbrIndexes to access the @p labels and compute a parity
       * or descriptor. This function returns 0 or 1 for both tetrahedral and
       * cistrans stereo centers.
       */
      int getDescriptor(const std::vector<unsigned int> &symmetry_classes, const std::vector<unsigned int> &labels) const
      {
        // Unspecified stereo centers have their own descriptor.
        if (nbrIndexes1.empty())
          return 2;
        std::vector<unsigned long> refs1, refs2;
        for (std::size_t i = 0; i < nbrIndexes1.size(); ++i) {
          if (nbrIndexes1[i] < labels.size())
            refs1.push_back(labels[nbrIndexes1[i]]);
          else
            refs1.push_back(nbrIndexes1[i]);
        }
        for (std::size_t i = 0; i < nbrIndexes2.size(); ++i)
          if (nbrIndexes2[i] < labels.size())
            refs2.push_back(labels[nbrIndexes2[i]]);
          else
            refs2.push_back(nbrIndexes2[i]);
        if (indexes.size() == 2) {
          bool symOrder = symmetry_classes[indexes[0]] < symmetry_classes[indexes[1]];
          bool canOrder = labels[indexes[0]] < labels[indexes[1]];
          if (symOrder != canOrder && symmetry_classes[indexes[0]] != symmetry_classes[indexes[1]])
            std::swap(refs1[0], refs1[1]);
        }
        return ((OBStereo::NumInversions(refs1) % 2 + OBStereo::NumInversions(refs2) % 2) % 2);
      }
      /**
       * Index(es) for the stereo center. This is used to sort the stereo
       * centers by label. Tetrahedral stereo centers have one index in this
       * vector and cistrans stereocenters have both double bond atom indexes
       * in this vector.
       */
      std::vector<unsigned int> indexes;
      /**
       * Indexes for the stereo center neighbor atoms. Tetrahedral centers have
       * all neighbor atoms in nbrIndexes1. CisTrans stereo centers store the
       * neighbor atoms for each double bond atom seperatly.
       */
      std::vector<unsigned int> nbrIndexes1, nbrIndexes2;
    };

    /**
     * Sort StereoCenter objects by their label. Used in combination with
     * std::sort to create STEREO code.
     */
    struct SortStereoCenters
    {
      const std::vector<unsigned int> &labels;
      SortStereoCenters(const std::vector<unsigned int> &_labels) : labels(_labels) {}
      inline unsigned int getLabel(const StereoCenter &c) const
      {
        switch (c.indexes.size()) {
          case 2:
            return std::min(labels[c.indexes[0]], labels[c.indexes[1]]);
          default:
            return labels[c.indexes[0]];
        }
      }
      inline bool operator()(const StereoCenter &c1, const StereoCenter &c2) const
      {
        return (getLabel(c1) < getLabel(c2));
      }
    };

    /**
     * Sort FullCode objects by their code. Used in combination with
     * std::sort order disconnected fragments.
     */
    static bool SortCode(const FullCode &code1, const FullCode &code2)
    {
      return (code1.code < code2.code);
    }

    /**
     * Sort FullCode objects with associated int in std::pair. Used to sort
     * ligands (i.e. equivalent fragments connected to a symmetric atom).
     */
    static bool SortCode2(const std::pair<int, FullCode> &code1, const std::pair<int, FullCode> &code2)
    {
      return (code1.second.code > code2.second.code);
    }

    /**
     * Sort atoms by ascending ranks (e.g. labels).
     */
    struct SortAtomsAscending
    {
      SortAtomsAscending(const std::vector<unsigned int> &_ranks) : ranks(_ranks) {}
      const std::vector<unsigned int> &ranks;
      inline bool operator()(const OBAtom *a1, const OBAtom *a2) const
      {
        return ranks[a1->GetIndex()] < ranks[a2->GetIndex()];
      }
    };

    /**
     * Sort atoms by descending ranks (e.g. symmetry classes)
     */
    struct SortAtomsDescending
    {
      SortAtomsDescending(const std::vector<unsigned int> &_ranks) : ranks(_ranks) {}
      const std::vector<unsigned int> &ranks;
      inline bool operator()(const OBAtom *a1, const OBAtom *a2) const
      {
        return ranks[a1->GetIndex()] > ranks[a2->GetIndex()];
      }
    };

    /**
     * Structure used while labeling a single connected fragment. All input is
     * specified using the constructor and @p code is generated as a result of
     * assigning labels.
     */
    struct State
    {
      State(const std::vector<unsigned int> &_symmetry_classes,
            const OBBitVec &_fragment, std::vector<StereoCenter> &_stereoCenters,
            std::vector<FullCode> &_identityCodes, Orbits &_orbits, OBBitVec &_mcr,
            bool _onlyOne) : symmetry_classes(_symmetry_classes), fragment(_fragment),
          onlyOne(_onlyOne), stereoCenters(_stereoCenters),
          code(_symmetry_classes.size()), identityCodes(_identityCodes),
          backtrackDepth(0), orbits(_orbits), mcr(_mcr)
      {
        mcr.Clear();
        if (mcr.IsEmpty())
          for (std::size_t i = 0; i < symmetry_classes.size(); ++i)
            mcr.SetBitOn(i+1);
      }

      const std::vector<unsigned int> &symmetry_classes;
      /**
       * The connected fragment. This is a subset of the mask.
       */
      const OBBitVec &fragment;
      const bool onlyOne;

      /**
       * The pre-computed stereo centers. Non-const since it needs to be
       * sorted.
       */
      std::vector<StereoCenter> &stereoCenters;
      /**
       * The partial code with labels and FROM code.
       */
      PartialCode code;

      /**
       * Identity nodes of the search tree.
       */
      std::vector<FullCode> identityCodes;
      unsigned int backtrackDepth;
      Orbits orbits;
      OBBitVec &mcr;
    };

    /**
     * Simple struct to manage timeout.
     */
    struct Timeout
    {
      Timeout(time_t _maxTime) : maxTime(_maxTime)
      {
        startTime = time(NULL);
      }
      time_t startTime, maxTime;
    };


    /**
     * Given a complete labeling (@p state.code.labels) and FROM code (@p state.code.from),
     * construct a FullCode with complete code.
     */
    static void CompleteCode(OBMol *mol, FullCode &fullcode, State &state)
    {
      PartialCode &code = state.code;

      // initialize the FullCode object with the found labels and the from code
      fullcode = FullCode(code.labels, code.from);
      unsigned int numClosures = 0;

      //
      // the RING-CLOSURE list
      //
      unsigned int current_label = 1; // the second atom is always labeled from 1
      // iterate over all atoms in the order they are labeled
      // (this ensures [1 3] < [3 1])
      for (std::size_t j = 0; j < code.atoms.size(); ++j) {
        OBAtom *atom = code.atoms[j];
        // still need to sort [1 3] and [1 4]
        std::vector<std::pair<OBBond*, unsigned int> > closures;
        FOR_BONDS_OF_ATOM (bond, atom) {
          // skip atoms not in the fragment
          if (!state.fragment.BitIsSet(bond->GetNbrAtom(atom)->GetIdx()))
            continue;
          // a closure bond is a bond not found while generating the FROM spanning tree.
          if (std::find(code.bonds.begin(), code.bonds.end(), &*bond) == code.bonds.end()) {
            closures.push_back(std::make_pair(&*bond, code.labels[bond->GetNbrAtom(atom)->GetIndex()]));
          }
        }

        // do the sorting: [1 3] < [1 4]
        std::sort(closures.begin(), closures.end(), CompareBondPairSecond);

        for (std::size_t k = 0; k < closures.size(); ++k) {
          // add the closure bond to the code
          fullcode.code.push_back(current_label);
          fullcode.code.push_back(closures[k].second);
          // add the bond to the list (needed for BOND-TYPES below)
          code.add(closures[k].first);
          numClosures++;
        }

        current_label++;
      }

      // Isotopes are only considered if there are isotopes.
      bool hasIsotope = false;
      // Charges are only considered if there is a formal charge.
      bool hasCharge = false;

      //
      // the ATOM-TYPES list
      //
      OBStereoFacade facade(mol);
      for (std::size_t j = 0; j < code.atoms.size(); ++j) {
        OBAtom *atom = code.atoms[j];
        if (atom->GetIsotope())
          hasIsotope = true;
        if (atom->GetFormalCharge())
          hasCharge = true;

	      // Include all hydrogens
        int hydrogens_to_include = TotalHydrogenCount(atom);

        unsigned int c = 10000 * atom->GetSpinMultiplicity() +
                          1000 * hydrogens_to_include +
                                 atom->GetAtomicNum();

        // add the atomic number to the code
        fullcode.code.push_back(c);
      }

      //
      // the (optional) ISOTOPES list
      //
      if (hasIsotope)
        for (std::size_t j = 0; j < code.atoms.size(); ++j)
          fullcode.code.push_back(code.atoms[j]->GetIsotope());

      //
      // the (optional) CHARGES list
      //
      if (hasCharge)
        for (std::size_t j = 0; j < code.atoms.size(); ++j)
          fullcode.code.push_back(7 + code.atoms[j]->GetFormalCharge());

      //
      // the BOND-TYPES list
      //
      for (std::size_t j = 0; j < code.bonds.size(); ++j) {
        OBBond *bond = code.bonds[j];
        if (bond->IsAromatic())
          fullcode.code.push_back(5);
        else
          fullcode.code.push_back(bond->GetBondOrder());
      }

      // the STEREO flag
      if (state.stereoCenters.size()) {
        // sort the stereo centers
        std::sort(state.stereoCenters.begin(), state.stereoCenters.end(), SortStereoCenters(code.labels));

        for (std::size_t i = 0; i < state.stereoCenters.size(); ++i) {
          bool isInFragment = false;
          for (std::size_t j = 0; j < state.stereoCenters[i].indexes.size(); ++j)
            if (state.fragment.BitIsSet(state.stereoCenters[i].indexes[j]+1))
              isInFragment = true;
          // ignore stereo centers not in this fragment
          if (isInFragment)
            fullcode.code.push_back(state.stereoCenters[i].getDescriptor(state.symmetry_classes, code.labels));
        }
      }

      // backtrack
      for (unsigned int i = 0; i < numClosures; ++i)
        code.bonds.pop_back();
    }

    /**
     * This function implements optimization 1 as described above.
     */
    static void LabelFragments(OBAtom *current, std::vector<OBAtom*> &nbrs, unsigned int label, Timeout &timeout, FullCode &bestCode, State &state)
    {
      OBMol *mol = current->GetParent();
      PartialCode &code = state.code;

      // Sort the neighbor atoms according to symmetry class.
      std::sort(nbrs.begin(), nbrs.end(), SortAtomsDescending(state.symmetry_classes));

      // Count the number of unique neighbor symmetry classes.
      unsigned int numUnique = 1;
      unsigned int lastSymClass = state.symmetry_classes[nbrs[0]->GetIndex()];
      for (std::size_t i = 0; i < nbrs.size(); ++i) {
        unsigned int symClass = state.symmetry_classes[nbrs[i]->GetIndex()];
        if (symClass != lastSymClass)
          numUnique++;
        lastSymClass = state.symmetry_classes[nbrs[i]->GetIndex()];
      }

      if (numUnique < nbrs.size()) {
        // The canonical codes for the equivalent ligands
        std::vector< std::pair<int, CanonicalLabelsImpl::FullCode> > lcodes;

        std::vector<unsigned int> ligandSizes;
        for (std::size_t i = 0; i < nbrs.size(); ++i) {
          OBBitVec ligand = getFragment(nbrs[i], current, state.fragment);
          ligandSizes.push_back(ligand.CountBits());

          CanonicalLabelsImpl::FullCode lbestCode;
          if (ligandSizes.back() == 1) {
            // Avoid additional state creation if the neighbor is a single terminal atom.
            lbestCode.code.push_back(nbrs[i]->GetAtomicNum());
            lbestCode.labels.resize(state.symmetry_classes.size(), 0);
            lbestCode.labels[nbrs[i]->GetIndex()] = 1;
          } else if (ligandSizes.back() == 2) {
            // Avoid additional state creation if the neighbor is a fragment of 2 atoms.
            lbestCode.labels.resize(state.symmetry_classes.size(), 0);
            lbestCode.code.push_back(1); // FROM
            lbestCode.code.push_back(nbrs[i]->GetAtomicNum()); // ATOM-TYPES 1
            lbestCode.labels[nbrs[i]->GetIndex()] = 1;
            FOR_NBORS_OF_ATOM (nbr, nbrs[i]) {
              if (!state.fragment.BitIsSet(nbr->GetIdx()))
                continue;
              if (code.labels[nbr->GetIndex()])
                continue;
              lbestCode.code.push_back(nbr->GetAtomicNum()); // ATOM-TYPES 2
              OBBond *bond = mol->GetBond(nbrs[i], &*nbr);
              if (bond->IsAromatic())
                lbestCode.code.push_back(5);
              else
                lbestCode.code.push_back(bond->GetBondOrder()); // BOND-TYPES 1
              lbestCode.labels[nbr->GetIndex()] = 2;
            }
          } else {
            // Start labeling from the ligand atom.
            std::vector<CanonicalLabelsImpl::FullCode> identityCodes;
            Orbits orbits;
            OBBitVec mcr;
            State lstate(state.symmetry_classes, ligand, state.stereoCenters, identityCodes, orbits, mcr, state.onlyOne);
            lstate.code.add(nbrs[i]);
            lstate.code.labels[nbrs[i]->GetIndex()] = 1;
            CanonicalLabelsRecursive(nbrs[i], 1, timeout, lbestCode, lstate);
          }

          // Store the canonical code (and labels) for the ligand.
          lcodes.push_back(std::make_pair(i, lbestCode));
        }

        // Sort the codes for the fragments. Each neighbor symmetry class is sorted separately.
        unsigned int firstIndex = 0;
        lastSymClass = state.symmetry_classes[nbrs[0]->GetIndex()];
        for (std::size_t i = 1; i < nbrs.size(); ++i) {
          unsigned int symClass = state.symmetry_classes[nbrs[i]->GetIndex()];
          if (symClass != lastSymClass) {
            std::sort(lcodes.begin() + firstIndex, lcodes.begin() + i, SortCode2);
            firstIndex = i;
          }
          lastSymClass = state.symmetry_classes[nbrs[i]->GetIndex()];
        }
        // Make sure to sort the last set too.
        std::sort(lcodes.begin() + firstIndex, lcodes.end(), SortCode2);

        // Label the neighbor atoms by updating code.
        std::vector<OBAtom*> atoms;
        unsigned int nextLbl = label + 1;
        for (std::size_t l = 0; l < lcodes.size(); ++l) {
          //print_vector("LIG CODE", lcodes[l].second.code);
          OBAtom *atom = nbrs[lcodes[l].first];
          code.add(current, atom);
          code.labels[atom->GetIndex()] = nextLbl;
          atoms.push_back(atom);
          nextLbl++;
        }

        // Convert the labels from the ligands to labels in the whole fragment.
        unsigned int ligandSize = *std::max_element(ligandSizes.begin(), ligandSizes.end());
        for (unsigned int lbl = 1; lbl < ligandSize; ++lbl) {
          for (std::size_t l = 0; l < lcodes.size(); ++l) {
            if (lbl >= ligandSizes[lcodes[l].first])
              continue;

            OBAtom *atom = 0;
            for (std::size_t i = 0; i < mol->NumAtoms(); ++i)
              if (lcodes[l].second.labels[i] == lbl) {
                atom = mol->GetAtom(i+1);
                break;
              }

            if (!atom)
              continue;

            std::vector<OBAtom*> atomNbrs;
            FOR_NBORS_OF_ATOM (nbr, atom) {
              if (lcodes[l].second.labels[nbr->GetIndex()] > lbl) {
                if (std::find(atoms.begin(), atoms.end(), &*nbr) != atoms.end())
                  continue;
                atomNbrs.push_back(&*nbr);
              }
            }

            std::sort(atomNbrs.begin(), atomNbrs.end(), SortAtomsAscending(lcodes[l].second.labels));

            for (std::size_t i = 0; i < atomNbrs.size(); ++i) {
              code.add(atom, atomNbrs[i]);
              code.labels[atomNbrs[i]->GetIndex()] = nextLbl;
              atoms.push_back(atomNbrs[i]);
              nextLbl++;
            }
          }
        }

        // Recurse...
        CanonicalLabelsRecursive(current, nextLbl - 1, timeout, bestCode, state);

        // Backtrack.
        for (std::size_t j = 0; j < atoms.size(); ++j) {
          code.atoms.pop_back();
          code.bonds.pop_back();
          code.from.pop_back();
          code.labels[atoms[j]->GetIndex()] = 0;
        }
      } else {
        // There are no duplicated symmetry classes for the neighbor atoms.
        // Label each neighbor atom in sorted sequence.
        unsigned int lbl = label;
        for (std::size_t i = 0; i < nbrs.size(); ++i) {
          lbl++;
          code.add(current, nbrs[i]);
          code.labels[nbrs[i]->GetIndex()] = lbl;
        }

        // Recurse...
        CanonicalLabelsRecursive(current, lbl, timeout, bestCode, state);

        // Backtrack.
        for (std::size_t i = 0; i < nbrs.size(); ++i) {
          code.atoms.pop_back();
          code.bonds.pop_back();
          code.from.pop_back();
          code.labels[nbrs[i]->GetIndex()] = 0;
        }
      }
    }

    static void FindOrbits(Orbits &orbits, OBMol *mol, const std::vector<unsigned int> &labels1,
        const std::vector<unsigned int> &labels2)
    {
      Orbits newOrbits;
      std::vector<bool> visited(labels1.size(), false);

      //print_vector("labels1", labels1);
      //print_vector("labels2", labels2);

      // Construct the orbits for automorphic permutation labels1 -> labels2
      for (std::size_t i = 0; i < labels1.size(); ++i) {
        if (visited[i])
          continue;
        unsigned int vi = labels1[i];
        if (vi == labels2[i])
          continue;
        //cout << "( ";
        std::size_t j = i;
        unsigned int vj = labels2[j];
        //cout << j << " ";
        newOrbits.resize(newOrbits.size()+1);
        newOrbits.back().push_back(mol->GetAtom(j+1));
        visited[i] = true;
        while (vi != vj) {
          for (std::size_t k = i; k < labels1.size(); ++k) {
            if (vj == labels1[k]) {
              j = k;
              vj =  labels2[j];
              //cout << j << " ";
              newOrbits.back().push_back(mol->GetAtom(j+1));
              visited[j] = true;
              break;
            }
          }
        }
        //cout << ") ";
      }
      //cout << endl;

      //print_orbits("newOrbits", newOrbits);

      // Merge the orbits with previously found orbits.
      for (std::size_t j = 0; j < newOrbits.size(); ++j) {
        bool merge = false;
        for (std::size_t k = 0; k < orbits.size(); ++k) {

          for (std::size_t l = 0; l < newOrbits[j].size(); ++l) {
            if (std::find(orbits[k].begin(), orbits[k].end(), newOrbits[j][l]) != orbits[k].end())
              merge = true;
          }

          if (merge) {
            for (std::size_t l = 0; l < newOrbits[j].size(); ++l)
              if (std::find(orbits[k].begin(), orbits[k].end(), newOrbits[j][l]) == orbits[k].end())
                orbits[k].push_back(newOrbits[j][l]);
            break;
          }
        }

        if (!merge)
          orbits.push_back(newOrbits[j]);
      }

//      print_orbits("orbits", orbits);

      newOrbits.clear();
      std::vector<bool> vivisted(orbits.size(), false);
      for (std::size_t i = 0; i < orbits.size(); ++i) {
        if (visited[i])
          continue;
        visited[i] = true;

        Orbit newOrbit = orbits[i];

        for (std::size_t j = i; j < orbits.size(); ++j) {
          if (visited[j])
            continue;

          std::sort(newOrbit.begin(), newOrbit.end());
          std::sort(orbits[j].begin(), orbits[j].end());
          Orbit result;

          std::set_intersection(newOrbit.begin(), newOrbit.end(),
                                orbits[j].begin(), orbits[j].end(),
                                std::back_inserter(result));
          if (result.empty())
            continue;

          visited[j] = true;
          result.clear();
          std::set_union(newOrbit.begin(), newOrbit.end(),
                         orbits[j].begin(), orbits[j].end(),
                         std::back_inserter(result));

          newOrbit = result;
        }

        newOrbits.push_back(newOrbit);
      }

      orbits = newOrbits;
    }

    /**
     * Update the minimum cell representations (mcr).
     */
    static void UpdateMcr(OBBitVec &mcr, Orbits &orbits, const std::vector<unsigned int> &bestLabels)
    {
      //print_orbits("UpdateMcr", orbits);

      for (std::size_t i = 0; i < bestLabels.size(); ++i)
        mcr.SetBitOn(i+1);

      for (std::size_t j = 0; j < orbits.size(); ++j) {
        std::sort(orbits[j].begin(), orbits[j].end(), SortAtomsAscending(bestLabels));

        for (std::size_t k = 0; k < orbits[j].size(); ++k)
          if (k)
            mcr.SetBitOff(orbits[j][k]->GetIdx());
      }
    }


    /**
     * This is the recursive function implementing the labeling algorithm
     * outlined above (steps 1-3). This function works for single connected
     * fragments and uses @p state to hold most information. The @p bestCode
     * is the current greatest code (if any). given the @p current atom and
     * the last @p label, this function labels the neighbor atoms with
     * label+1, label+2, ...
     */
    static void CanonicalLabelsRecursive(OBAtom *current, unsigned int label, Timeout &timeout, FullCode &bestCode, State &state)
    {
      OBMol *mol = current->GetParent();
      PartialCode &code = state.code;

      if (state.backtrackDepth) {
        //std::cout << "backtrackDepth = " << state.backtrackDepth << std::endl;

        if (code.atoms.size() > state.backtrackDepth) {
          //std::cout << "BACKTRACKING" << std::endl;
          return;
        }
        if (code.atoms.size() == state.backtrackDepth) {
          //std::cout << "BACKTRACK DONE 1" << std::endl;
          state.backtrackDepth = 0;
          return;
        } else
        if (code.atoms.size() < state.backtrackDepth) {
          //std::cout << "BACKTRACK DONE 2" << std::endl;
          state.backtrackDepth = 0;
        }
      }


      // Check if there is a full mapping.
      if (label == state.fragment.CountBits()) {
        // Complete the canonical code.
        FullCode fullcode;
        CompleteCode(mol, fullcode, state);

        //print_vector("TERMINAL", fullcode.code);

        // Check previously found codes to find redundant subtrees.
        for (std::size_t i = state.identityCodes.size(); i > 0; --i)
          if (fullcode.code == state.identityCodes[i-1].code) {
            //
            // An explicit automorphism has been found.
            //
            std::vector<unsigned int> v1(fullcode.labels.size(), 0);
            for (std::size_t j = 0; j < fullcode.labels.size(); ++j)
              if (fullcode.labels[j])
                v1[fullcode.labels[j]-1] = j + 1;

            std::vector<unsigned int> v2(state.identityCodes[i-1].labels.size(), 0);
            for (std::size_t j = 0; j < state.identityCodes[i-1].labels.size(); ++j)
              if (state.identityCodes[i-1].labels[j])
                v2[state.identityCodes[i-1].labels[j]-1] = j + 1;

            state.backtrackDepth = 0;
            for (std::size_t j = 0; j < v1.size(); ++j) {
              if (v1[j] != v2[j]) {
                // Yield backtrackDepth.
                return;
              }

              if (v1[j])
                state.backtrackDepth++;
            }
          }

        if (fullcode.code == bestCode.code) {
          UpdateMcr(state.mcr, state.orbits, bestCode.labels);
          FindOrbits(state.orbits, mol, fullcode.labels, bestCode.labels);
        } else if (fullcode > bestCode) {
          // if fullcode is greater than bestCode, we have found a new greatest code
          bestCode = fullcode;
        }


        if (state.identityCodes.size() < MAX_IDENTITY_NODES) {
          state.identityCodes.push_back(FullCode());
          state.identityCodes.back().labels.swap(fullcode.labels);
          state.identityCodes.back().code.swap(fullcode.code);
        } else {
          state.identityCodes[MAX_IDENTITY_NODES-1].labels.swap(fullcode.labels);
          state.identityCodes[MAX_IDENTITY_NODES-1].code.swap(fullcode.code);
        }


        return;
      }

      // Avoid endless loops.
      if (time(NULL) - timeout.startTime > timeout.maxTime) {
        return;
      }

      // If there is a bestCode and only one labeling is required, return.
      if (state.onlyOne && !bestCode.code.empty())
        return;

      // Abort early if this will not lead to a greatest canonical code. The
      // code.from vector is compared with the elements in the bestCode.
      if (code < bestCode)
        return;

      // Find the neighbors of the current atom to assign the next label(s).
      std::vector<OBAtom*> nbrs;
      std::vector<unsigned int> nbrSymClasses;

      FOR_NBORS_OF_ATOM (nbr, current) {
        // Skip atoms not in the fragment.
        if (!state.fragment.BitIsSet(nbr->GetIdx()))
          continue;
        // Skip already labeled atoms.
        if (code.labels[nbr->GetIndex()])
          continue;

        OBBond *bond = mol->GetBond(current, &*nbr);
        // Ugly, but it helps...                                        // <--- not documented! need better approach to handle this (not only ferrocene)...
        if (!isFerroceneBond(bond)) {
          nbrSymClasses.push_back(state.symmetry_classes[nbr->GetIndex()]);
          nbrs.push_back(&*nbr);
        }
      }

      if (nbrs.empty()) {
        // If there are no neighbor atoms to label, recurse with the next
        // current atom.
        unsigned int nextLabel = code.labels[current->GetIndex()] + 1;
        for (std::size_t i = 0; i < code.labels.size(); ++i) {
          if (code.labels[i] == nextLabel) {
            CanonicalLabelsRecursive(mol->GetAtom(i+1), label, timeout, bestCode, state);
            return;
          }
        }
        return;
      }

      if (!current->IsInRing()) {
        // Optimization to avoid generating the n! permutations. If the current
        // atom is not in a ring, the neighbor atoms with the same symmetry
        // are not connected and it is possible to label without considering
        // the rest of the molecule. The found ligand labels can be sorted and
        // the labels can be calculated using offsets.
        //
        // This modifies the final canonical code since it produces a different
        // canonical order but although different, the code is still canonical.
        LabelFragments(current, nbrs, label, timeout, bestCode, state);
      } else {
        // Create all possible labelings for the neighbors.
        //
        // Example: The current atom has 4 neighbors to label.
        //
        // neighbor atom   : 1 2 3 4
        // symmetry classes: 3 2 2 1
        //
        // The first neighbor with symmetry class 3 is added to allOrderedNbrs.
        //
        // allOrderedNbrs = [ [1] ]
        //
        // The two neighbor atoms with symmetry class 2 are added next. All
        // permutations are added (2 in this case).
        //
        // allOrderedNbrs = [ [1,2,3], [1,3,2] ]
        //
        // The last atom is similar to the first one.
        //
        // allOrderedNbrs = [ [1,2,3,4], [1,3,2,4] ]
        //
        // Note: If there are atoms with a large number of neighbor atoms
        // with the same symmetry class (n), this can result in a large number
        // of permutations (n!).
        std::vector<std::vector<OBAtom*> > allOrderedNbrs(1);
        while (!nbrs.empty()) {

          // Select the next nbr atoms with highest symmetry classes.
          unsigned int maxSymClass = *std::min_element(nbrSymClasses.begin(), nbrSymClasses.end());
          std::vector<OBAtom*> finalNbrs;
          for (std::size_t i = 0; i < nbrs.size(); ++i) {
            if (nbrSymClasses[i] == maxSymClass)
              finalNbrs.push_back(nbrs[i]);
          }

          // Remove the selected atoms from nbrs and nbrSymClasses (this could be made more efficient)
          for (std::size_t i = 0; i < finalNbrs.size(); ++i) {
            nbrs.erase(std::find(nbrs.begin(), nbrs.end(), finalNbrs[i]));
            nbrSymClasses.erase(std::find(nbrSymClasses.begin(), nbrSymClasses.end(), maxSymClass));
          }


          if (finalNbrs.size() == 1) {
            // If there is only one atom with the same symmetry class, label it
            // and select the next group of neighbor atoms with the same symmetry
            // class.
            for (std::size_t i = 0; i < allOrderedNbrs.size(); ++i)
              allOrderedNbrs[i].push_back(finalNbrs[0]);
          } else {
            // Sort the atoms lexicographically.
            std::sort(finalNbrs.begin(), finalNbrs.end());

            // Copy the current labelings for the neighbor atoms.
            std::vector<std::vector<OBAtom*> > allOrderedNbrsCopy(allOrderedNbrs);

            // Add the first permutation for the neighbor atoms.
//            if (state.mcr.BitIsSet(finalNbrs[0]->GetIdx()))
              for (std::size_t j = 0; j < allOrderedNbrs.size(); ++j) {
                for (std::size_t i = 0; i < finalNbrs.size(); ++i) {
                  allOrderedNbrs[j].push_back(finalNbrs[i]);
                }
              }

            // Add the other permutations.
            while (std::next_permutation(finalNbrs.begin(), finalNbrs.end())) {
              if (state.mcr.BitIsSet(finalNbrs[0]->GetIdx()))
                for (std::size_t j = 0; j < allOrderedNbrsCopy.size(); ++j) {
                  allOrderedNbrs.push_back(allOrderedNbrsCopy[j]);
                  for (std::size_t i = 0; i < finalNbrs.size(); ++i)
                    allOrderedNbrs.back().push_back(finalNbrs[i]);
                }
            }

          } // finalNbrs.size() != 1
        } // while (!nbrs.empty())

        if (DEBUG) {
          cout << "allOrderedNbrs:" << endl;
          for (std::size_t i = 0; i < allOrderedNbrs.size(); ++i) {
            for (std::size_t j = 0; j < allOrderedNbrs[i].size(); ++j) {
              cout << allOrderedNbrs[i][j]->GetIndex() << " ";
            }
            cout << endl;
          }
        }

        for (std::size_t i = 0; i < allOrderedNbrs.size(); ++i) {
          // Convert the order stored in allOrderedNbrs to labels.
          unsigned int lbl = label;
          for (std::size_t j = 0; j < allOrderedNbrs[i].size(); ++j) {
            lbl++;
            code.add(current, allOrderedNbrs[i][j]);
            code.labels[allOrderedNbrs[i][j]->GetIndex()] = lbl;
          }

          // Recurse...
          CanonicalLabelsRecursive(current, lbl, timeout, bestCode, state);

          // Backtrack...
          for (std::size_t j = 0; j < allOrderedNbrs[i].size(); ++j) {
            code.atoms.pop_back();
            code.bonds.pop_back();
            code.from.pop_back();
            code.labels[allOrderedNbrs[i][j]->GetIndex()] = 0;
          }

          // Optimization
          if (state.backtrackDepth) {
            if (code.atoms.size() <= state.backtrackDepth) {
              //std::cout << "BACKTRACK DONE 3" << std::endl;
              state.backtrackDepth = 0;
            }
          }

        }

      } // current->IsInRing()

    }

    /**
     * Select an initial atom from a fragment to assign the first label.
     */
    static std::vector<OBAtom*> findStartAtoms(OBMol *mol, const OBBitVec &fragment, const std::vector<unsigned int> &symmetry_classes)
    {
      // find the a symmetry class in the fragment using criteria
      std::vector<unsigned int> ranks;
      for (std::size_t i = 0; i < mol->NumAtoms(); ++i) {
        if (!fragment.BitIsSet(i+1))
          continue;

        OBAtom *atom = mol->GetAtom(i+1);
        unsigned int rank = 10000 * symmetry_classes[i]  +
                             1000 * atom->GetSpinMultiplicity() +
                               10 * (atom->GetFormalCharge() + 7) +
                                    TotalHydrogenCount(atom);

        ranks.push_back(rank);
      }

      unsigned int lowestRank = *std::min_element(ranks.begin(), ranks.end());

      std::vector<OBAtom*> result;
      for (std::size_t i = 0; i < mol->NumAtoms(); ++i) {
        if (!fragment.BitIsSet(i+1))
          continue;

        OBAtom *atom = mol->GetAtom(i+1);
        unsigned int rank = 10000 * symmetry_classes[i]  +
                             1000 * atom->GetSpinMultiplicity() +
                               10 * (atom->GetFormalCharge() + 7) +
                                    TotalHydrogenCount(atom);

        if (rank == lowestRank)
          result.push_back(atom);
      }

      return result;
    }

    /**
     * Compute the canonical labels for @p mol. This function handles a whole molecule with
     * disconnected fragments. It finds a canonical code for each fragment, sorts these codes
     * and computes labels for the molecule as a whole.
     *
     * This is the CanonicalLabelsImpl entry point.
     */
    static void CalcCanonicalLabels(OBMol *mol, const std::vector<unsigned int> &symmetry_classes,
        std::vector<unsigned int> &canonical_labels, const OBStereoUnitSet &stereoUnits, const OBBitVec &mask,
        OBStereoFacade *stereoFacade, int maxSeconds, bool onlyOne = false)
    {
      // Handle some special cases.
      if (!mol->NumAtoms())
        return;
      if (mol->NumAtoms() == 1) {
        canonical_labels.resize(1, 1);
        return;
      }

      // Set all labels to 0.
      canonical_labels.clear();
      canonical_labels.resize(mol->NumAtoms(), 0);

      std::vector<OBBond*> metalloceneBonds;
      findMetalloceneBonds(metalloceneBonds, mol, symmetry_classes);

      // Find the (dis)connected fragments.
      OBBitVec visited;
      std::vector<OBBitVec> fragments;
      for (std::size_t i = 0; i < mol->NumAtoms(); ++i) {
        if (!mask.BitIsSet(i+1) || visited.BitIsSet(i+1))
          continue;
        fragments.push_back(getFragment(mol->GetAtom(i+1), mask, metalloceneBonds));
        visited |= fragments.back();
      }

      // Pre-compute the stereo center information. (See StereoCenter)
      std::vector<StereoCenter> stereoCenters;
      if (stereoFacade) {
        for (std::size_t i = 0; i < stereoUnits.size(); ++i) {
          const OBStereoUnit &unit = stereoUnits[i];

          if (unit.type == OBStereo::Tetrahedral) {
            OBAtom *atom = mol->GetAtomById(unit.id);
            if (!atom)
              continue;
            // Add the StereoCenter indexes.
            stereoCenters.resize(stereoCenters.size()+1);
            stereoCenters.back().indexes.push_back(atom->GetIndex());

            if (!stereoFacade->HasTetrahedralStereo(unit.id))
              continue;
            OBTetrahedralStereo::Config config = stereoFacade->GetTetrahedralStereo(unit.id)->GetConfig();
            if (!config.specified)
              continue;

            // Add the neighbor atom indexes.
            OBAtom *from = mol->GetAtomById(config.from);
            if (from && from->GetAtomicNum() != OBElements::Hydrogen)
              stereoCenters.back().nbrIndexes1.push_back(from->GetIndex());
            else
              stereoCenters.back().nbrIndexes1.push_back(std::numeric_limits<unsigned int>::max());
            for (std::size_t j = 0; j < config.refs.size(); ++j) {
              OBAtom *ref = mol->GetAtomById(config.refs[j]);
              if (ref && ref->GetAtomicNum() != OBElements::Hydrogen)
                stereoCenters.back().nbrIndexes1.push_back(ref->GetIndex());
              else
                stereoCenters.back().nbrIndexes1.push_back(std::numeric_limits<unsigned int>::max());
            }
          } else if (unit.type == OBStereo::CisTrans) {
            OBBond *bond = mol->GetBondById(unit.id);
            if (!bond || bond->IsAromatic())
              continue;
            OBAtom *begin = bond->GetBeginAtom();
            OBAtom *end = bond->GetEndAtom();
            if (!begin || !end)
              continue;
            // Add the StereoCenter indexes.
            stereoCenters.resize(stereoCenters.size()+1);
            stereoCenters.back().indexes.push_back(begin->GetIndex());
            stereoCenters.back().indexes.push_back(end->GetIndex());

            if (!stereoFacade->HasCisTransStereo(unit.id))
              continue;
            OBCisTransStereo::Config config = stereoFacade->GetCisTransStereo(unit.id)->GetConfig();
            if (!config.specified)
              continue;

            // Add the neighbor atom indexes.
            for (std::size_t j = 0; j < config.refs.size(); ++j) {
              OBAtom *ref = mol->GetAtomById(config.refs[j]);
              unsigned int r = (ref && ref->GetAtomicNum() != OBElements::Hydrogen) ? ref->GetIndex() : std::numeric_limits<unsigned int>::max();
              if (stereoCenters.back().nbrIndexes1.size() < 2)
                stereoCenters.back().nbrIndexes1.push_back(r);
              else
                stereoCenters.back().nbrIndexes2.push_back(r);
            }
          }
        }
      }

      // Find the canonical code for each fragment.
      std::vector<CanonicalLabelsImpl::FullCode> fcodes;
      for (std::size_t f = 0; f < fragments.size(); ++f) {
        const OBBitVec &fragment = fragments[f];

        // Select the first atom.
        std::vector<OBAtom*> startAtoms = findStartAtoms(mol, fragment, symmetry_classes);

        CanonicalLabelsImpl::Timeout timeout(maxSeconds);
        CanonicalLabelsImpl::FullCode bestCode;
        std::vector<CanonicalLabelsImpl::FullCode> identityCodes;
        Orbits orbits;
        OBBitVec mcr;

        for (std::size_t i = 0; i < startAtoms.size(); ++i) {
          OBAtom *atom = startAtoms[i];

          // Start labeling of the fragment.
          State state(symmetry_classes, fragment, stereoCenters, identityCodes, orbits, mcr, onlyOne);
          //if (!state.mcr.BitIsSet(atom->GetIdx()) && atom->IsInRing())
          //  continue;

          state.code.add(atom);
          state.code.labels[atom->GetIndex()] = 1;
          CanonicalLabelsRecursive(atom, 1, timeout, bestCode, state);
        }

        // Throw an error if the timeout is exceeded.
        if (time(NULL) - timeout.startTime > timeout.maxTime) {
          obErrorLog.ThrowError(__FUNCTION__, "maximum time exceeded...", obError);
        }

        // Store the canonical code for the fragment.
        fcodes.push_back(bestCode);
      }

      // Sort the codes for the fragments.
      std::sort(fcodes.begin(), fcodes.end(), SortCode);

      // Construct the full labeling from the sorted fragment labels.
      unsigned int offset = 0;
      for (std::size_t f = 0; f < fcodes.size(); ++f) {
        //print_vector("CODE", fcodes[f].code);
        //print_vector("code_labels", fcodes[f].labels);
        if (fcodes[f].labels.size() == 0)
          continue; // defensive programming

        unsigned int max_label = 0;
        for (std::size_t i = 0; i < mol->NumAtoms(); ++i) {
          if (fcodes[f].labels[i]) {
            canonical_labels[i] = fcodes[f].labels[i] + offset;
            max_label = std::max(max_label, canonical_labels[i]);
          }
        }
        offset = max_label;
      }

    }

  }; // CanonicalLabelsImpl

  /*
   * See header for detailed description of parameters.
   *
   * The main purpose of this function is calling CanonicalLabelsImpl::CalcCanonicalLabels
   * with the correct parameters regarding stereochemistry.
   */
  void CanonicalLabels(OBMol *mol, const std::vector<unsigned int> &symmetry_classes,
      std::vector<unsigned int> &canonical_labels, const OBBitVec &mask,
      int maxSeconds, bool onlyOne)
  {
    // make sure the mask is valid: no mask = all atoms
    OBBitVec maskCopy(mask);
    if (!maskCopy.CountBits())
      FOR_ATOMS_OF_MOL (atom, mol)
        maskCopy.SetBitOn(atom->GetIdx());

    if (onlyOne) {
      // Only one labeling requested. This results in canonical labels that do not
      // consider stereochemistry. Used for finding stereo centers with automorphisms.
      CanonicalLabelsImpl::CalcCanonicalLabels(mol, symmetry_classes, canonical_labels, OBStereoUnitSet(), maskCopy, 0, maxSeconds, true);
    } else {
      std::vector<OBBond*> metalloceneBonds;
      findMetalloceneBonds(metalloceneBonds, mol, symmetry_classes);

      // Check if there are specified stereo centers.
      bool hasAtLeastOneDefined = false;
      OBStereoFacade sf(mol, false);
      FOR_ATOMS_OF_MOL (atom, mol) {
        FOR_BONDS_OF_ATOM (bond, &*atom)
        if (std::find(metalloceneBonds.begin(), metalloceneBonds.end(), &*bond) != metalloceneBonds.end())
          continue;
        if (sf.HasTetrahedralStereo(atom->GetId())) {
          if (sf.GetTetrahedralStereo(atom->GetId())->GetConfig().specified) {
            hasAtLeastOneDefined = true;
            break;
          }
        }
      }
      FOR_BONDS_OF_MOL (bond, mol) {
        if (sf.HasCisTransStereo(bond->GetId())) {
          if (sf.GetCisTransStereo(bond->GetId())->GetConfig().specified) {
            hasAtLeastOneDefined = true;
            break;
          }
        }
      }
      if (!hasAtLeastOneDefined) {
        // If there are no specified stereo centers, we don't need to find stereogenic units.
        CanonicalLabelsImpl::CalcCanonicalLabels(mol, symmetry_classes, canonical_labels, OBStereoUnitSet(), maskCopy, 0, maxSeconds);
        return;
      }

      // Find the stereogenic units
      OBStereoUnitSet stereoUnits = FindStereogenicUnits(mol, symmetry_classes);

      // Mark all invalid stereo data as unspecified
      std::vector<OBGenericData*> stereoData = mol->GetAllData(OBGenericDataType::StereoData);
      std::vector<OBGenericData*>::iterator data;
      for (data = stereoData.begin(); data != stereoData.end(); ++data) {
        OBStereo::Type type = ((OBStereoBase*)*data)->GetType();
        if (type == OBStereo::Tetrahedral) {
          OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
          OBTetrahedralStereo::Config config = ts->GetConfig();
          bool valid = true;
          if (!ts->IsValid())
            valid = false;
          OBAtom *center = mol->GetAtomById(config.center);
          if (!center)
            valid = false;
          else if (!isTetrahedral(center, stereoUnits))
            valid = false;

          if (!valid) {
            config.specified = false;
            ts->SetConfig(config);
          }
        }
        if (type == OBStereo::CisTrans) {
          OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
          OBCisTransStereo::Config config = ct->GetConfig();
          bool valid = true;
          if (!ct->IsValid())
            valid = false;
          OBAtom *beginAtom = mol->GetAtomById(config.begin);
          OBAtom *endAtom = mol->GetAtomById(config.end);
          if (!beginAtom || !endAtom)
            valid = false;
          else {
            OBBond *bond = mol->GetBond(beginAtom, endAtom);
            if (!bond)
              valid = false;
            else if (!isCisTrans(bond, stereoUnits))
              valid = false;
          }

          if (!valid) {
            config.specified = false;
            ct->SetConfig(config);
          }
        }
      }

      // Determine stereochemistry from coordinates if needed
      if (!mol->HasChiralityPerceived()) {
        switch (mol->GetDimension()) {
          case 2:
            mol->DeleteData(OBGenericDataType::StereoData);
            TetrahedralFrom2D(mol, stereoUnits);
            CisTransFrom2D(mol, stereoUnits);
            break;
          case 3:
            mol->DeleteData(OBGenericDataType::StereoData);
            TetrahedralFrom3D(mol, stereoUnits);
            CisTransFrom3D(mol, stereoUnits);
            break;
          default:
            TetrahedralFrom0D(mol, stereoUnits);
            CisTransFrom0D(mol, stereoUnits);
            break;
        }
      }

      // The mol->DeleteData() above deleted some of the data that the OBStereoFacade sf()
      // is using.  We have to create a new one.  (CJ - fixes a bug.)
      OBStereoFacade newsf(mol, false);

      // Start the labeling process
      CanonicalLabelsImpl::CalcCanonicalLabels(mol, symmetry_classes, canonical_labels, stereoUnits, maskCopy, &newsf, maxSeconds);
    }

    // if the labeling failed, just return the identity labels to avoid craches
    if (canonical_labels.empty())
      for (std::size_t i = 0; i < symmetry_classes.size(); ++i)
        canonical_labels.push_back(i+1);
  }



/* -*-C++-*-
+======================================================================
|
|                       !!!DEPRECATED!!!
|
|
| AUTHOR: Craig A. James, eMolecules, Inc.
|
| DESCRIPTION: CANONICALIZATION OF SMILES
|
|       This is a specialized SMILES canonicalization algorithm.  Although
|       it can be applied in the standard fashion to a whole molecule,
|       its real job is to generate canonical SMILES for fragments, or
|       "subsets", of the atoms of a molecule.
|
|       For example, consider the first three atoms of Oc1ccccc1.  With
|       a "normal" SMILES canonicalizer, you couldn't generate a SMILES
|       for Occ, because it's not a valid molecule.  However, this system
|       can do exactly that, by taking both the whole molecule (which
|       retains the aromaticity), and a "subset" bitmap that specifies
|       which atoms are to be included in the SMILES.
|
|       Canonicalization is carried out per Weininger et al (J. Chem.
|       Inf. Comput. Sci., Vol. 29, No. 2, 1989, pp 97-101), with some
|       modifications to handle bond symmetries not foreseen by Weininger
|       in that paper.
|
|       WARNING - KNOWN BUG: These functions make use of a bitmap vector
|       to represent a "fragment" -- a subset of the atoms in a molecule.
|       But this means the bonds of the fragment are implicit, not explicit,
|       which is incorrect.  For example, if you want to break one bond of
|       cyclehexane (C1CCCCC1), all six atoms will still be there, so the
|       "fragment" will be cyclic.  This is relevant when generating fragment
|       SMILES for ring systems where breaking a bond can reduce the number
|       of ring without removing any atoms.  We need to add a pair of bit
|       vectors, the atoms AND the bonds, to represent a fragment.  (Note
|       that this is also an ambiguity in OpenBabel itself, which represents
|       a ring as a set of atoms. This is only valid if the ring is a member
|       of a SSSR.)
+======================================================================
*/

} // namespace OpenBabel

//! \file canon.cpp
//! \brief Canonical numbering of SMILES, molecules and fragments
