/**********************************************************************
  graphsym.cpp - Symmetry detection algorithm

  Copyright (C) 2009-2010 by Tim Vandermeersch
  Copyright (C) 2005-2006, eMolecules, Inc. (www.emolecules.com)

  This file is part of the Open Babel project.
  For more information, see <http://openbabel.org/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/

#include <openbabel/graphsym.h>
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <openbabel/obiter.h>
#include <openbabel/generic.h>

#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/elements.h>

#include <iterator> // std::istream_iterator
#include <cassert>
#include <algorithm>
#include <cmath>
#include <limits>

#include "stereo/stereoutil.h"

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

// debug function
void print_sym_classes(const std::string &label, const std::vector<std::pair<OpenBabel::OBAtom*, unsigned int> > &atom_sym_classes)
{
  cout << label << ": ";
  for (unsigned int i = 0; i < atom_sym_classes.size(); i++)
    cout << atom_sym_classes[i].second << " ";
  cout << endl;
}

namespace OpenBabel {

  class OBGraphSymPrivate
  {
    public:
      OBBitVec _frag_atoms;
      OBMol* _pmol;
      std::vector<unsigned int> _canonLabels;
      OBStereoUnitSet _stereoUnits;

      unsigned int GetHvyDegree(OBAtom *atom);
      unsigned int GetHvyBondSum(OBAtom *atom);
      void FindRingAtoms(OBBitVec &ring_atoms);
      void CreateNewClassVector(std::vector<std::pair<OBAtom*,unsigned int> > &vp1,
                                std::vector<std::pair<OBAtom*,unsigned int> > &vp2);
      static void CreateNewClassVector(OBMol *mol, std::vector<std::pair<OBAtom*,unsigned int> > &vp1,
                                       std::vector<std::pair<OBAtom*,unsigned int> > &vp2);
      void GetGIVector(std::vector<unsigned int> &vid);
      bool GetGTDVector(std::vector<int> &gtd);
      static void CountAndRenumberClasses(std::vector<std::pair<OBAtom*,unsigned int> > &vp, unsigned int &count);
      int ExtendInvariants(std::vector<std::pair<OBAtom*, unsigned int> > &symmetry_classes);
      int CalculateSymmetry(std::vector<unsigned int> &symmetry_classes);
      int Iterate(std::vector<unsigned int> &symmetry_classes);
      void CanonicalLabels(const std::vector<unsigned int> &symmetry_classes, std::vector<unsigned int> &canon_labels, int maxSeconds);
  };

  OBGraphSym::OBGraphSym(OBMol* pmol, const OBBitVec* frag_atoms) : d(new OBGraphSymPrivate)
  {
    d->_pmol = pmol;
    if (frag_atoms) {
      d->_frag_atoms = *frag_atoms;
    } else {
      d->_frag_atoms.Resize(d->_pmol->NumAtoms());
      FOR_ATOMS_OF_MOL(a, d->_pmol)
        d->_frag_atoms.SetBitOn(a->GetIdx());
    }
  }

  OBGraphSym::~OBGraphSym()
  {
    delete d;
  }

  const unsigned int OBGraphSym::NoSymmetryClass = 0x7FFFFFFF;

  /**
   * Functions for use by the sort() method of a vector.
   */
  inline bool CompareUnsigned(const unsigned int &a,const unsigned int &b)
  {
    return(a<b);
  }

  inline bool ComparePairFirst(const std::pair<OBAtom*,unsigned int> &a,const std::pair<OBAtom*,unsigned int> &b)
  {
    return(a.first->GetIdx() < b.first->GetIdx());
  }

  inline bool ComparePairSecond(const std::pair<OBAtom*,unsigned int> &a,const std::pair<OBAtom*,unsigned int> &b)
  {
    return(a.second < b.second);
  }


  /**
   * Like OBAtom::GetHvyDegree(): Counts the number non-hydrogen
   * neighbors, but doesn't count atoms not in the fragment.
   */
  unsigned int OBGraphSymPrivate::GetHvyDegree(OBAtom *atom)
  {
    unsigned int count = 0;
    OBBond *bond;
    OBAtom *nbr;

    vector<OBBond*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms.BitIsSet(nbr->GetIdx()) && nbr->GetAtomicNum() != OBElements::Hydrogen)
        count++;
    }

    return(count);
  }

  /**
   * Sums the bond order over the bonds from this atom to other atoms
   * in the fragment.  Single = 1, double = 2, triple = 3, aromatic = 1.6,
   * but sum is rounded to nearest integer.
   *
   * This is used for fragment symmetry perception instead of the "implicit
   * valence" used by the standard OpenBabel symmetry perception.  It
   * has the same effect, but we don't have to worry about hydrogen counts,
   * EXCEPT for aromatic N, where the difference between n and [nH] is
   * critical.
   */
  unsigned int OBGraphSymPrivate::GetHvyBondSum(OBAtom *atom)
  {
    float count = 0.0f;
    OBBond *bond;
    OBAtom *nbr;

    vector<OBBond*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms.BitIsSet(nbr->GetIdx()) && nbr->GetAtomicNum() != OBElements::Hydrogen) {
        if (bond->IsAromatic())
          count += 1.6f;
        else
          count += (float)bond->GetBondOrder();
      }
    }
    if (atom->GetAtomicNum() == 7 && atom->IsAromatic() && atom->GetTotalDegree() == 3) {
      count += 1.0f;         // [nH] - add another bond
    }
    return(int(count + 0.5));     // round to nearest int
  }

  /**
   * Calculates the graph theoretical distance of each atom.
   * Vector is indexed from zero.
   *
   * NOTE: "Indexed from zero" means it's one off from the atom->GetIdx()
   * that's used to index atoms inside the molecule!
   *
   * NOTE: This function is hard to decipher, and seems to be misnamed.
   * A "distance" should be be between two atoms, but there's more here
   * than that.  It seems to be doing a breadth-first search to find the
   * most-distant atom from each atom, and reporting the number of steps
   * (which happens to be the graph-theoretical distance) to that atom.
   * The name "Graph Theoretical Distance" is thus misleading.
   */
  bool OBGraphSymPrivate::GetGTDVector(vector<int> &gtd)
  {
    gtd.clear();
    gtd.resize(_pmol->NumAtoms());

    int gtdcount, natom;
    OBBitVec used, curr, next;
    OBAtom *atom, *atom1;
    OBBond *bond;
    vector<OBNodeBase*>::iterator ai;
    vector<OBBond*>::iterator j;

    next.Clear();

    for (atom = _pmol->BeginAtom(ai); atom; atom = _pmol->NextAtom(ai)) {

      int idx = atom->GetIdx();
      if (!_frag_atoms.BitIsSet(idx)) {     // Not in this fragment?
        gtd[idx-1] = OBGraphSym::NoSymmetryClass;
        continue;
      }

      gtdcount = 0;
      used.Clear();curr.Clear();
      used.SetBitOn(idx);
      curr.SetBitOn(idx);

      while (!curr.IsEmpty()) {
        next.Clear();
        for (natom = curr.NextBit(-1);natom != curr.EndBit();natom = curr.NextBit(natom)) {
          atom1 = _pmol->GetAtom(natom);
          if (!_frag_atoms.BitIsSet(atom1->GetIdx()))
            continue;
          for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j)) {
            int nbr_idx = bond->GetNbrAtomIdx(atom1);
            if (   _frag_atoms.BitIsSet(nbr_idx)
                && !used.BitIsSet(nbr_idx)
                && !curr.BitIsSet(nbr_idx)
                && bond->GetNbrAtom(atom1)->GetAtomicNum() != OBElements::Hydrogen)
              next.SetBitOn(nbr_idx);
          }
        }
        used |= next;
        curr = next;
        gtdcount++;
      }
      gtd[idx-1] = gtdcount;
    }

    return(true);
  }

  /**
   * Finds all atoms that are part of a ring in the current fragment.
   * We start with the whole molecule's rings, and eliminate any that
   * have atoms not in the subset.  For the rings that are left, mark
   * each atom of the ring as a ring atom.
   *
   * @return A bit vector where TRUE means it's a ring atom.
   */
  void OBGraphSymPrivate::FindRingAtoms(OBBitVec &ring_atoms)
  {
    vector<OBRing*> sssRings;
    vector<OBRing*>::iterator ri;

    ring_atoms.Resize(_pmol->NumAtoms());
    ring_atoms.Clear();

    sssRings = _pmol->GetSSSR();
    for (ri = sssRings.begin(); ri != sssRings.end(); ++ri) {
      OBRing *ring = *ri;
      OBBitVec bvtmp = _frag_atoms & ring->_pathset;      // intersection: fragment and ring
      if (bvtmp == ring->_pathset)                        // all ring atoms in fragment?
        ring_atoms |= ring->_pathset;                     //   yes - add this ring's atoms
    }
  }

  /**
   * Calculates a set of graph invariant indexes using the graph theoretical
   * distance, number of connected heavy atoms, aromatic boolean, ring
   * boolean, atomic number, and summation of bond orders connected to the
   * atom.
   *
   * We have to recalculate which atoms are in rings by taking the fragment's
   * atoms into account when we generate the graph invarients.
   *
   * Vector is indexed from zero (not one, like atom->GetIdx()).
   *
   * NOTE: This may need to be extended to include the bond-invariant properties,
   * particularly the size of all rings the bond is in (from a SSSR).
   */
  void OBGraphSymPrivate::GetGIVector(vector<unsigned int> &vid)
  {
    // Prepare the vector...
    vid.clear();
    vid.resize(_pmol->NumAtoms());

    // The "graph theoretical distance" for each atom (see comments in the function)
    vector<int> v;
    GetGTDVector(v);

    // Compute the ring atoms for this particular fragment (set of atoms)
    OBBitVec ring_atoms;
    FindRingAtoms(ring_atoms);

    int i;
    OBAtom *atom;
    vector<OBNodeBase*>::iterator ai;
    for (i=0, atom = _pmol->BeginAtom(ai); atom; atom = _pmol->NextAtom(ai)) {
      //    vid[i] = 0;
      vid[i] = OBGraphSym::NoSymmetryClass;
      if (_frag_atoms.BitIsSet(atom->GetIdx())) {
        vid[i] =
          v[i]                                                    // 10 bits: graph-theoretical distance
          | (GetHvyDegree(atom)                <<10)  //  4 bits: heavy valence
          | (((atom->IsAromatic()) ? 1 : 0)                <<14)  //  1 bit:  aromaticity
          | (((ring_atoms.BitIsSet(atom->GetIdx())) ? 1 : 0)<<15)  //  1 bit:  ring atom
          | (atom->GetAtomicNum()                          <<16)  //  7 bits: atomic number
          | (GetHvyBondSum(atom)               <<23)  //  4 bits: heavy bond sum
          | ((7 + atom->GetFormalCharge())                 <<27); //  4 bits: formal charge
      }
      i++;
    }
  }

  /**
   * Creates a new vector of symmetry classes based on an existing
   * vector.  (Helper routine to GetGIDVector.)  On return, vp2 will
   * have newly-extended connectivity sums, but the numbers (the class
   * IDs) are very large.
   *
   * (Comments by CJ) This appears to compute the "extended connectivity
   * sums" similar to those described by Weininger, Morgan, etc. It uses
   * vp1 as its starting point (the current connectivity sums), and puts
   * the new sums in vp2.  Note that vp1 is modified along the way.
   *
   * Note that, per Weininger's warning, this assumes the initial class
   * ID's are less than 100, which is a BAD assumption, e.g. OCC...CCN
   * would have more than 100 symmetry classes if the chain is more than
   * 98 carbons long.  Should change this to use Weininger's product of
   * corresponding primes.
   */
  void OBGraphSymPrivate::CreateNewClassVector(std::vector<std::pair<OBAtom*,unsigned int> > &vp1,
                                        std::vector<std::pair<OBAtom*,unsigned int> > &vp2)
  {
    int m,id;
    OBAtom *atom, *nbr;
    vector<OBBond*>::iterator nbr_iter;
    vector<unsigned int>::iterator k;
    vector<pair<OBAtom*,unsigned int> >::iterator vp_iter;

#if DEBUG2
    cout << "CreateNewClassVector: START\n";
    //print_vector_pairs("    ", vp1);
#endif

    // There may be fewer atoms than in the whole molecule, so we can't
    // index the vp1 array by atom->GetIdx().  Instead, create a quick
    // mapping vector of idx-to-index for vp1.
    vector<int> idx2index(_pmol->NumAtoms() + 1, -1);  // natoms + 1
    int index = 0;
    for (vp_iter = vp1.begin(); vp_iter != vp1.end(); ++vp_iter) {
      int idx = vp_iter->first->GetIdx();
      idx2index[idx] = index++;
    }

    // vp2 will hold the newly-extended symmetry classes
    vp2.resize(vp1.size());
    vp2.clear();

    // Loop over original atoms.
    // Create a new extended varient for each atom.  Get its neighbors' class ID's,
    // sort them into ascending order, and create a sum of (c0 + c1*10^2 + c2*10^4 + ...)
    // which becomes the new class ID (where c0 is the current classID).

    for (vp_iter = vp1.begin(); vp_iter != vp1.end(); ++vp_iter) {
      atom = vp_iter->first;
      id   = vp_iter->second;
      vector<unsigned int> vtmp;
      for (nbr = atom->BeginNbrAtom(nbr_iter); nbr; nbr = atom->NextNbrAtom(nbr_iter)) {
        int idx = nbr->GetIdx();
        if (_frag_atoms.BitIsSet(idx))
          vtmp.push_back(vp1[idx2index[idx]].second);
      }

      sort(vtmp.begin(),vtmp.end(),CompareUnsigned);
      for (m = 100, k = vtmp.begin(); k != vtmp.end(); ++k, m*=100)
        id += *k * m;
      vp2.push_back(pair<OBAtom*,unsigned int> (atom, id));
    }
#if DEBUG2
    cout << "CreateNewClassVector: FINISH\n";
    //print_vector_pairs("    ", vp2);
#endif
  }

  void OBGraphSymPrivate::CreateNewClassVector(OBMol *mol, std::vector<std::pair<OBAtom*,unsigned int> > &vp1,
      std::vector<std::pair<OBAtom*,unsigned int> > &vp2)
  {
    int m,id;
    OBAtom *atom, *nbr;
    vector<OBBond*>::iterator nbr_iter;
    vector<unsigned int>::iterator k;
    vector<pair<OBAtom*,unsigned int> >::iterator vp_iter;

#if DEBUG2
    cout << "CreateNewClassVector: START\n";
    //print_vector_pairs("    ", vp1);
#endif

    // There may be fewer atoms than in the whole molecule, so we can't
    // index the vp1 array by atom->GetIdx().  Instead, create a quick
    // mapping vector of idx-to-index for vp1.
    vector<int> idx2index(mol->NumAtoms() + 1, -1);  // natoms + 1
    int index = 0;
    for (vp_iter = vp1.begin(); vp_iter != vp1.end(); ++vp_iter) {
      int idx = vp_iter->first->GetIdx();
      idx2index[idx] = index++;
    }

    // vp2 will hold the newly-extended symmetry classes
    vp2.resize(vp1.size());
    vp2.clear();

    // Loop over original atoms.
    // Create a new extended varient for each atom.  Get its neighbors' class ID's,
    // sort them into ascending order, and create a sum of (c0 + c1*10^2 + c2*10^4 + ...)
    // which becomes the new class ID (where c0 is the current classID).

    for (vp_iter = vp1.begin(); vp_iter != vp1.end(); ++vp_iter) {
      atom = vp_iter->first;
      id   = vp_iter->second;
      vector<unsigned int> vtmp;
      for (nbr = atom->BeginNbrAtom(nbr_iter); nbr; nbr = atom->NextNbrAtom(nbr_iter)) {
        int idx = nbr->GetIdx();
        vtmp.push_back(vp1[idx2index[idx]].second);
      }

      sort(vtmp.begin(),vtmp.end(),CompareUnsigned);
      for (m = 100, k = vtmp.begin(); k != vtmp.end(); ++k, m*=100)
        id += *k * m;
      vp2.push_back(pair<OBAtom*,unsigned int> (atom, id));
    }
#if DEBUG2
    cout << "CreateNewClassVector: FINISH\n";
    //print_vector_pairs("    ", vp2);
#endif

  }

  /**
   * Counts the number of unique symmetry classes in a list.
   *
   * (NOTE: CJ -- It also appears to MODIFY the list.  It sorts it in order
   * of class ID, then renumbers the ID's zero through N-1.  See the comments
   * in CreateNewClassVector() about how it returns very large numbers for the
   * class IDs it creates.  These are replaced by lower, sequential numbers here.)
   */
  void OBGraphSymPrivate::CountAndRenumberClasses(std::vector<std::pair<OBAtom*,unsigned int> > &vp,
                                           unsigned int &count)
  {
    count = 1;
    vector<pair<OBAtom*,unsigned int> >::iterator k;

    sort(vp.begin(), vp.end(), ComparePairSecond);
    k = vp.begin();
    if (k != vp.end()) {
      unsigned int id = k->second;
      if (id) {
      k->second = 1;
      ++k;
      for (;k != vp.end(); ++k) {
        if (k->second != id) {
          id = k->second;
          k->second = ++count;
        } else {
          k->second = count;
        }
      }
      }
    }
  }

  /**
   * This is the core of symmetry analysis.  Starting with a set of
   * classes on each atom, it "spreads" them using a sum-of-invariants
   * of each atom's class and its neighbors' classes.  This iterates
   * until a stable solution is found (further spreading doesn't
   * change the answer).
   *
   * @return The number of distinct symmetry classes found.
   */
  int OBGraphSymPrivate::ExtendInvariants(std::vector<std::pair<OBAtom*, unsigned int> > &symmetry_classes)
  {
    unsigned int nclasses1, nclasses2;
    vector<pair<OBAtom*,unsigned int> > tmp_classes;

    // How many classes are we starting with?  (The "renumber" part isn't relevant.)
    CountAndRenumberClasses(symmetry_classes, nclasses1);

    unsigned int nfragatoms = _frag_atoms.CountBits();

    // LOOP: Do extended sum-of-invarients until no further changes are
    // noted.  (Note: This is inefficient, as it re-computes extended sums
    // and re-sorts the entire list each time.  You can save a lot of time by
    // only recomputing and resorting within regions where there is a tie
    // initially.  But it's a lot more code.)

    if (nclasses1 < nfragatoms) {
      for (int i = 0; i < 100;i++) {  //sanity check - shouldn't ever hit this number
        CreateNewClassVector(symmetry_classes, tmp_classes);
        CountAndRenumberClasses(tmp_classes, nclasses2);
        symmetry_classes = tmp_classes;
        if (nclasses1 == nclasses2) break;
        nclasses1 = nclasses2;
      }
    }

    CreateNewClassVector(symmetry_classes, tmp_classes);
    CountAndRenumberClasses(tmp_classes, nclasses2);


    if (nclasses1 != nclasses2) {
      symmetry_classes = tmp_classes;
      return ExtendInvariants(symmetry_classes);
    }


    return nclasses1;
  }

  /**
   * Calculates a set of canonical symmetry identifiers for a molecule.
   * Atoms with the same symmetry ID are symmetrically equivalent.  By
   * "canonical", we mean it generates a repeatable labelling of the
   * atoms, i.e. the same fragment will get the same symmetry labels in
   * any molecule in which it occurs.
   *
   * Vector is indexed from zero, corresponding to (atom->GetIdx() - 1).
   *
   * The bit vector "_frag_atoms" specifies a fragment of the molecule,
   * where each bit represents the presence or absence of the atom in
   * the fragment.  Symmetry is computed as though the fragment is the
   * only part that exists.
   */
  int OBGraphSymPrivate::CalculateSymmetry(std::vector<unsigned int> &atom_sym_classes)
  {
    vector<unsigned int> vgi;
    vector<OBNodeBase*>::iterator j;
    OBAtom *atom;

    // Get vector of graph invariants.  These are the starting "symmetry classes".
    GetGIVector(vgi);


    // Create a vector-of-pairs, associating each atom with its Class ID.
    std::vector<std::pair<OBAtom*, unsigned int> > symmetry_classes;
    for (atom = _pmol->BeginAtom(j); atom; atom = _pmol->NextAtom(j)) {
      int idx = atom->GetIdx();
      if (_frag_atoms.BitIsSet(idx))
        symmetry_classes.push_back(pair<OBAtom*, unsigned int> (atom, vgi[idx-1]));
      //else
      //  symmetry_classes.push_back(pair<OBAtom*, unsigned int> (atom, OBGraphSym::NoSymmetryClass));
    }


    // The heart of the matter: Do extended sum-of-invariants until no further
    // changes are noted.
    int nclasses = ExtendInvariants(symmetry_classes);

    // Convert to a vector indexed by Index
    // Atoms not in the fragment will have a value of OBGraphSym::NoSymmetryClass
    atom_sym_classes.clear();
    atom_sym_classes.resize(_pmol->NumAtoms(), OBGraphSym::NoSymmetryClass);
    for (unsigned int i = 0; i < symmetry_classes.size(); ++i) {
      atom_sym_classes[symmetry_classes.at(i).first->GetIndex()] = symmetry_classes.at(i).second;
    }

    // Store the symmetry classes in an OBPairData
    stringstream temp;
    vector<unsigned int>::iterator sym_iter = atom_sym_classes.begin();
    if (sym_iter != atom_sym_classes.end())
      temp << (*sym_iter++);
    for (; sym_iter != atom_sym_classes.end(); ++sym_iter)
      temp << " " << (*sym_iter);

    OBPairData *symData = new OBPairData;
    symData->SetAttribute("OpenBabel Symmetry Classes");
    symData->SetOrigin(local); //will not show as sdf or cml property
    symData->SetValue(temp.str());
    _pmol->SetData(symData);

    return nclasses;
  }

  int OBGraphSym::GetSymmetry(std::vector<unsigned int> &symmetry_classes)
  {
    ClearSymmetry(); // For the moment just recalculate the symmetry classes

    // Check to see whether we have already calculated the symmetry classes
    OBPairData *pd = dynamic_cast<OBPairData*>(d->_pmol->GetData("OpenBabel Symmetry Classes"));

    int nclasses = 0;
    if (!pd) {
      nclasses = d->CalculateSymmetry(symmetry_classes);
    } else {
      istringstream iss(pd->GetValue());
      symmetry_classes.clear();
      copy(istream_iterator<unsigned int>(iss),
           istream_iterator<unsigned int>(),
           back_inserter<vector<unsigned int> >(symmetry_classes));
      // Now find the number of unique elements
      vector<unsigned int> copy_sym = symmetry_classes;
      sort(copy_sym.begin(), copy_sym.end());
      vector<unsigned int>::iterator end_pos = unique(copy_sym.begin(), copy_sym.end()); // Requires sorted elements
      nclasses = end_pos - copy_sym.begin();
    }

    return nclasses;
  }

  int OBGraphSymPrivate::Iterate(vector<unsigned int> &symClasses)
  {
    // Create a vector-of-pairs, associating each atom with its Class ID.
    vector<OBAtom*>::iterator j;
    std::vector<std::pair<OBAtom*, unsigned int> > symmetry_classes;
    for (OBAtom *atom = _pmol->BeginAtom(j); atom; atom = _pmol->NextAtom(j)) {
      int idx = atom->GetIdx();
      if (_frag_atoms.BitIsSet(idx))
        symmetry_classes.push_back(pair<OBAtom*, unsigned int> (atom, symClasses[idx-1]));
    }

    // The heart of the matter: Do extended sum-of-invariants until no further
    // changes are noted.
    int nclasses = ExtendInvariants(symmetry_classes);

    // Convert to a vector indexed by Index
    // Atoms not in the fragment will have a value of OBGraphSym::NoSymmetryClass
    symClasses.clear();
    symClasses.resize(_pmol->NumAtoms(), OBGraphSym::NoSymmetryClass);
    for (unsigned int i = 0; i < symmetry_classes.size(); ++i) {
      symClasses[symmetry_classes.at(i).first->GetIndex()] = symmetry_classes.at(i).second;
    }

    return nclasses;
  }

  // Clears perceived symmetry
  void OBGraphSym::ClearSymmetry()
  {
    d->_pmol->DeleteData("OpenBabel Symmetry Classes");
  }
















} // namespace OpenBabel

//! \file graphsym.cpp
//! \brief Handle and perceive graph symmtery for canonical numbering.
