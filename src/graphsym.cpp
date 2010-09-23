/**********************************************************************
  graphsym.cpp - Symmetry detection algorithm

  Copyright (C) 2009-2010 by Tim Vandermeersch
  Copyright (C) 2005-2006, eMolecules, Inc. (www.emolecules.com)

  This file is part of the Open Babel project.
  For more information, see <http://openbabel.sourceforge.net/>

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

#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>

#include <iterator> // std::istream_iterator
#include <cassert>
#include <algorithm>
#include <cmath>
#include <limits>

#include "stereo/stereoutil.h"

#define DEBUG 0

  static const char *red    = "\033[1;31m";
  static const char *green  = "\033[1;32m";
  static const char *yellow = "\033[1;33m";
  static const char *blue   = "\033[1;34m";
  static const char *normal = "\033[0m";

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
  for (int i = 0; i < atom_sym_classes.size(); i++)
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

      unsigned int GetHvyValence(OBAtom *atom);
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

  inline bool CompareBondPairSecond(const std::pair<OBBond*,unsigned int> &a,const std::pair<OBBond*,unsigned int> &b)
  {
    return(a.second < b.second);
  }

  /**
   * Like OBAtom::GetHvyValence(): Counts the number non-hydrogen
   * neighbors, but doesn't count atoms not in the fragment.
   */
  unsigned int OBGraphSymPrivate::GetHvyValence(OBAtom *atom)
  {
    unsigned int count = 0;
    OBBond *bond;
    OBAtom *nbr;

    vector<OBEdgeBase*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms.BitIsSet(nbr->GetIdx()) && !(nbr->IsHydrogen()))
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

    vector<OBEdgeBase*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms.BitIsSet(nbr->GetIdx()) && !(nbr->IsHydrogen())) {
        if (bond->IsSingle())        count += 1.0f;
        else if (bond->IsDouble())   count += 2.0f;
        else if (bond->IsTriple())   count += 3.0f;
        else if (bond->IsAromatic()) count += 1.6f;
      }
    }
    if (atom->GetAtomicNum() == 7 && atom->IsAromatic() && atom->GetImplicitValence() == 3) {
      count += 1;         // [nH] - add another bond
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
    vector<OBEdgeBase*>::iterator j;

    next.Clear();

    for (atom = _pmol->BeginAtom(ai); atom; atom = _pmol->NextAtom(ai)) {

      int idx = atom->GetIdx();
      if (!_frag_atoms.BitIsOn(idx)) {     // Not in this fragment?
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
          if (!_frag_atoms.BitIsOn(atom1->GetIdx()))
            continue;
          for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j)) {
            int nbr_idx = bond->GetNbrAtomIdx(atom1);
            if (   _frag_atoms.BitIsOn(nbr_idx)
                && !used.BitIsOn(nbr_idx)
                && !curr.BitIsOn(nbr_idx)
                && !(bond->GetNbrAtom(atom1))->IsHydrogen())
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
    for (ri = sssRings.begin(); ri != sssRings.end(); ri++) {
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
      if (_frag_atoms.BitIsOn(atom->GetIdx())) {
        vid[i] =
          v[i]                                                    // 10 bits: graph-theoretical distance
          | (GetHvyValence(atom)                <<10)  //  4 bits: heavy valence
          | (((atom->IsAromatic()) ? 1 : 0)                <<14)  //  1 bit:  aromaticity
          | (((ring_atoms.BitIsOn(atom->GetIdx())) ? 1 : 0)<<15)  //  1 bit:  ring atom
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
    vector<OBEdgeBase*>::iterator nbr_iter;
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
    for (vp_iter = vp1.begin(); vp_iter != vp1.end(); vp_iter++) {
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

    for (vp_iter = vp1.begin(); vp_iter != vp1.end(); vp_iter++) {
      atom = vp_iter->first;
      id   = vp_iter->second;
      vector<unsigned int> vtmp;
      for (nbr = atom->BeginNbrAtom(nbr_iter); nbr; nbr = atom->NextNbrAtom(nbr_iter)) {
        int idx = nbr->GetIdx();
        if (_frag_atoms.BitIsOn(idx))
          vtmp.push_back(vp1[idx2index[idx]].second);
      }

      sort(vtmp.begin(),vtmp.end(),CompareUnsigned);
      for (m = 100, k = vtmp.begin(); k != vtmp.end(); k++, m*=100)
        id += *k * m;
      vp2.push_back(pair<OBAtom*,unsigned int> (atom, id));
    }
#if DEBUG2
    cout << "CreateNewClassVector: FINISH\n";
    //print_vector_pairs("    ", vp2);
#endif
  }

  void OBGraphSymPrivate::CreateNewClassVector(OBMol *mol, vector<pair<OBAtom*,unsigned int> > &vp1,
      vector<pair<OBAtom*,unsigned int> > &vp2)
  {
    int m,id;
    OBAtom *atom, *nbr;
    vector<OBEdgeBase*>::iterator nbr_iter;
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
    for (vp_iter = vp1.begin(); vp_iter != vp1.end(); vp_iter++) {
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

    for (vp_iter = vp1.begin(); vp_iter != vp1.end(); vp_iter++) {
      atom = vp_iter->first;
      id   = vp_iter->second;
      vector<unsigned int> vtmp;
      for (nbr = atom->BeginNbrAtom(nbr_iter); nbr; nbr = atom->NextNbrAtom(nbr_iter)) {
        int idx = nbr->GetIdx();
        vtmp.push_back(vp1[idx2index[idx]].second);
      }

      sort(vtmp.begin(),vtmp.end(),CompareUnsigned);
      for (m = 100, k = vtmp.begin(); k != vtmp.end(); k++, m*=100)
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

    int natoms = _pmol->NumAtoms();
    int nfragatoms = _frag_atoms.CountBits();

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
      if (_frag_atoms.BitIsOn(idx))
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
    if (DEBUG)
      cout << "ENTER GetSymmetry()" << endl;
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
      if (_frag_atoms.BitIsOn(idx))
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












  /**
   * Helper function for getFragment below.
   */
  void addNbrs(OBBitVec &fragment, OBAtom *atom, const OBBitVec &mask)
  {
    FOR_NBORS_OF_ATOM (nbr, atom) {
      if (!mask.BitIsSet(nbr->GetIdx()))
        continue;
      // skip visited atoms
      if (fragment.BitIsSet(nbr->GetIdx()))
        continue;
      // add the neighbor atom to the fragment
      fragment.SetBitOn(nbr->GetIdx());
      // recurse...
      addNbrs(fragment, &*nbr, mask);
    }
  }

  /**
   * Create an OBBitVec objects with bets set for the fragment consisting of all
   * atoms for which there is a path to atom without going through skip. These
   * fragment bitvecs are indexed by atom idx (i.e. OBAtom::GetIdx()).
   */
  OBBitVec getFragment(OBAtom *atom, const OBBitVec &mask)
  {
    OBBitVec fragment;
    fragment.SetBitOn(atom->GetIdx());
    // start the recursion
    addNbrs(fragment, atom, mask);
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

    if (Fe->GetValence() < 10)
      return false;

    return C->HasDoubleBond() && C->IsInRing();
  }



  struct CanonicalLabelsImpl
  {
    struct FullCode
    {
      std::vector<unsigned int> labels, code;

      FullCode()
      {
      }

      FullCode(const std::vector<unsigned int> &_labels, const std::vector<unsigned int> &from)
          : labels(_labels), code(from)
      {
      }

      inline bool operator>(const FullCode &other) const
      {
        return (code > other.code);
      }
    };

    struct PartialCode
    {
      std::vector<OBAtom*> atoms;
      std::vector<OBBond*> bonds;
      std::vector<unsigned int> labels, from;

      PartialCode(std::size_t numAtoms)
      {
        labels.resize(numAtoms, 0);
      }

      inline void add(OBAtom *atom) { atoms.push_back(atom); }
      inline void add(OBBond *bond) { bonds.push_back(bond); }
      inline void add(OBAtom *fromAtom, OBAtom *atom)
      {
        from.push_back(labels[fromAtom->GetIndex()]);
        atoms.push_back(atom);
        bonds.push_back(atom->GetParent()->GetBond(fromAtom, atom));
      }

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

    struct StereoCenter
    {
      int getDescriptor(const std::vector<unsigned int> &labels) const
      {
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
        return ((OBStereo::NumInversions(refs1) % 2 + OBStereo::NumInversions(refs2) % 2) % 2);
      }
      std::vector<unsigned int> indexes;
      std::vector<unsigned int> nbrIndexes1, nbrIndexes2;
    };

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

    static bool SortCode(const FullCode &code1, const FullCode &code2)
    {
      return (code1.code < code2.code);
    }

    static bool SortCode2(const std::pair<int, FullCode> &code1, const std::pair<int, FullCode> &code2)
    {
      return (code1.second.code > code2.second.code);
    }

    struct SortAtomsAscending
    {
      SortAtomsAscending(const std::vector<unsigned int> &_labels) : labels(_labels) {}
      const std::vector<unsigned int> &labels;
      inline bool operator()(const OBAtom *a1, const OBAtom *a2) const
      {
        return labels[a1->GetIndex()] < labels[a2->GetIndex()];
      }
    };

    struct SortAtomsDescending
    {
      SortAtomsDescending(const std::vector<unsigned int> &_labels) : labels(_labels) {}
      const std::vector<unsigned int> &labels;
      inline bool operator()(const OBAtom *a1, const OBAtom *a2) const
      {
        return labels[a1->GetIndex()] > labels[a2->GetIndex()];
      }
    };

    struct State
    {
      State(const std::vector<unsigned int> &_symmetry_classes,
            const OBBitVec &_fragment, std::vector<StereoCenter> &_stereoCenters,
            bool _onlyOne) : symmetry_classes(_symmetry_classes),
          fragment(_fragment), stereoCenters(_stereoCenters), onlyOne(_onlyOne),
          code(_symmetry_classes.size())
      {
      }

      const std::vector<unsigned int> &symmetry_classes;
      const OBBitVec &fragment;
      const bool onlyOne;

      std::vector<StereoCenter> &stereoCenters;
      PartialCode code;
    };

    struct Timeout
    {
      Timeout(unsigned int _maxTime) : maxTime(_maxTime)
      {
        startTime = time(NULL);
      }
      unsigned int startTime, maxTime;
    };


    static void CompleteCode(OBMol *mol, FullCode &bestCode, State &state)
    {
      PartialCode &code = state.code;

      FullCode fullcode(code.labels, code.from);
      unsigned int numClosures = 0;

      // the RING-CLOSURE list
      unsigned int current_label = 1;
      for (std::size_t j = 0; j < code.atoms.size(); ++j) {
        OBAtom *atom = code.atoms[j];
        std::vector<std::pair<OBBond*, unsigned int> > closures;
        FOR_BONDS_OF_ATOM (bond, atom) {
          if (!state.fragment.BitIsSet(bond->GetNbrAtom(atom)->GetIdx()))
            continue;
          if (std::find(code.bonds.begin(), code.bonds.end(), &*bond) == code.bonds.end()) {
            closures.push_back(std::make_pair(&*bond, code.labels[bond->GetNbrAtom(atom)->GetIndex()]));
          }
        }

        std::sort(closures.begin(), closures.end(), CompareBondPairSecond);
        for (std::size_t k = 0; k < closures.size(); ++k) {
          fullcode.code.push_back(current_label);
          fullcode.code.push_back(closures[k].second);
          code.add(closures[k].first);
          numClosures++;
        }

        current_label++;
      }

      bool hasIsotope = false;

      // the ATOM-TYPES list
      for (std::size_t j = 0; j < code.atoms.size(); ++j) {
        OBAtom *atom = code.atoms[j];
        if (atom->GetIsotope())
          hasIsotope = true;
        fullcode.code.push_back(atom->GetAtomicNum());
      }

      // the (optional) ISOTOPES list
      if (hasIsotope)
        for (std::size_t j = 0; j < code.atoms.size(); ++j)
          fullcode.code.push_back(code.atoms[j]->GetIsotope());

      // the BOND-TYPES list
      for (std::size_t j = 0; j < code.bonds.size(); ++j) {
        OBBond *bond = code.bonds[j];
        if (bond->IsAromatic())
          fullcode.code.push_back(5);
        else
          fullcode.code.push_back(bond->GetBondOrder());
      }

      // the STEREO flag
      if (state.stereoCenters.size()) {
        std::vector<int> dv;

        std::sort(state.stereoCenters.begin(), state.stereoCenters.end(), SortStereoCenters(code.labels));
        for (std::size_t i = 0; i < state.stereoCenters.size(); ++i) {
          bool isInFragment = false;
          for (std::size_t j = 0; j < state.stereoCenters[i].indexes.size(); ++j)
            if (state.fragment.BitIsSet(state.stereoCenters[i].indexes[j]+1))
              isInFragment = true;
          if (isInFragment)
            dv.push_back(state.stereoCenters[i].getDescriptor(code.labels));
        }

        int value = 0;
        for (unsigned int i = 0; i < dv.size(); ++i) {
          if (!dv[i])
            continue;
          int power = dv.size() - i - 1;
          // bit shift is equivalent to 2^power
          // i.e., 1 << 0 == 1, 1 << 1 == 2, etc.
          value += 1 << power;
        }

        fullcode.code.push_back(value);
      }

      if (fullcode > bestCode) {
        bestCode.labels.swap(fullcode.labels);
        bestCode.code.swap(fullcode.code);
      }

      for (unsigned int i = 0; i < numClosures; ++i)
        code.bonds.pop_back();
    }

    static void LabelFragments(OBAtom *current, std::vector<OBAtom*> &nbrs, unsigned int label, Timeout &timeout, FullCode &bestCode, State &state)
    {
      OBMol *mol = current->GetParent();
      PartialCode &code = state.code;

      std::sort(nbrs.begin(), nbrs.end(), SortAtomsDescending(state.symmetry_classes));

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
            // avoid additional state creation
            lbestCode.code.push_back(nbrs[i]->GetAtomicNum());
            lbestCode.labels.resize(state.symmetry_classes.size(), 0);
            lbestCode.labels[nbrs[i]->GetIndex()] = 1;
          } else if (ligandSizes.back() == 2) {
            // avoid additional state creation
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
            // start labeling from the ligand atom
            State lstate(state.symmetry_classes, ligand, state.stereoCenters, state.onlyOne);
            lstate.code.add(nbrs[i]);
            lstate.code.labels[nbrs[i]->GetIndex()] = 1;
            CanonicalLabelsRecursive(nbrs[i], 1, timeout, lbestCode, lstate);
          }

          lcodes.push_back(std::make_pair(i, lbestCode));
        }

        // sort the codes for the fragments
        unsigned int firstIndex = 0;
        lastSymClass = state.symmetry_classes[nbrs[0]->GetIndex()];
        for (std::size_t i = 1; i < nbrs.size(); ++i) {
          unsigned int symClass = state.symmetry_classes[nbrs[i]->GetIndex()];
          if (symClass != lastSymClass) {
            std::sort(lcodes.begin() + firstIndex, lcodes.begin() + i - 1, SortCode2);
            firstIndex = i;
          }
          lastSymClass = state.symmetry_classes[nbrs[i]->GetIndex()];
        }

        std::vector<OBAtom*> atoms;
        unsigned int nextLbl = label + 1;
        for (std::size_t l = 0; l < lcodes.size(); ++l) {
          OBAtom *atom = nbrs[lcodes[l].first];
          code.add(current, atom);
          code.labels[atom->GetIndex()] = nextLbl;
          atoms.push_back(atom);
          nextLbl++;
        }

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

        CanonicalLabelsRecursive(current, nextLbl - 1, timeout, bestCode, state);

        for (std::size_t j = 0; j < atoms.size(); ++j) {
          code.atoms.pop_back();
          code.bonds.pop_back();
          code.from.pop_back();
          code.labels[atoms[j]->GetIndex()] = 0;
        }
      } else {
        unsigned int lbl = label;
        for (std::size_t i = 0; i < nbrs.size(); ++i) {
          lbl++;
          code.add(current, nbrs[i]);
          code.labels[nbrs[i]->GetIndex()] = lbl;
        }

        CanonicalLabelsRecursive(current, lbl, timeout, bestCode, state);

        for (std::size_t i = 0; i < nbrs.size(); ++i) {
          code.atoms.pop_back();
          code.bonds.pop_back();
          code.from.pop_back();
          code.labels[nbrs[i]->GetIndex()] = 0;
        }
      }
    }


    static void CanonicalLabelsRecursive(OBAtom *current, unsigned int label, Timeout &timeout, FullCode &bestCode, State &state)
    {
      OBMol *mol = current->GetParent();
      PartialCode &code = state.code;

      // Check if there is a full mapping
      if (label == state.fragment.CountBits()) {
        CompleteCode(mol, bestCode, state);
        return;
      }

      // avoid endless loops
      if (time(NULL) - timeout.startTime > timeout.maxTime) {
        return;
      }

      if (state.onlyOne && !bestCode.code.empty())
        return;

      // abort early if this will not lead to a maximal canonical code
      if (code < bestCode)
        return;


      std::vector<OBAtom*> nbrs;
      std::vector<unsigned int> nbrSymClasses;
      FOR_NBORS_OF_ATOM (nbr, current) {
        if (!state.fragment.BitIsSet(nbr->GetIdx()))
          continue;
        if (code.labels[nbr->GetIndex()])
          continue;

        OBBond *bond = mol->GetBond(current, &*nbr);
        if (!isFerroceneBond(bond)) {
          nbrSymClasses.push_back(state.symmetry_classes[nbr->GetIndex()]);
          nbrs.push_back(&*nbr);
        }
      }

      if (nbrs.empty()) {
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
        LabelFragments(current, nbrs, label, timeout, bestCode, state);
      } else { // current->IsInRing


        std::vector<std::vector<OBAtom*> > allOrderedNbrs(1);
        while (!nbrs.empty()) {

          // select the next nbr atoms with highest symmetry classes
          unsigned int maxSymClass = *std::max_element(nbrSymClasses.begin(), nbrSymClasses.end());
          std::vector<OBAtom*> finalNbrs;
          for (std::size_t i = 0; i < nbrs.size(); ++i) {
            if (nbrSymClasses[i] == maxSymClass)
              finalNbrs.push_back(nbrs[i]);
          }

          // remove the selected atoms from nbrs and nbrSymClasses
          for (std::size_t i = 0; i < finalNbrs.size(); ++i) {
            nbrs.erase(std::find(nbrs.begin(), nbrs.end(), finalNbrs[i]));
            nbrSymClasses.erase(std::find(nbrSymClasses.begin(), nbrSymClasses.end(), maxSymClass));
          }


          if (finalNbrs.size() == 1) {
            for (std::size_t i = 0; i < allOrderedNbrs.size(); ++i)
              allOrderedNbrs[i].push_back(finalNbrs[0]);
          } else {
            std::sort(finalNbrs.begin(), finalNbrs.end());

            std::vector<std::vector<OBAtom*> > allOrderedNbrsCopy(allOrderedNbrs);

            for (std::size_t j = 0; j < allOrderedNbrsCopy.size(); ++j) {
              if (!allOrderedNbrsCopy[j].empty())
                allOrderedNbrs.push_back(allOrderedNbrsCopy[j]);
              for (std::size_t i = 0; i < finalNbrs.size(); ++i) {
                allOrderedNbrs.back().push_back(finalNbrs[i]);
              }
            }

            while (std::next_permutation(finalNbrs.begin(), finalNbrs.end())) {
              for (std::size_t j = 0; j < allOrderedNbrsCopy.size(); ++j) {
                allOrderedNbrs.push_back(allOrderedNbrsCopy[j]);
                for (std::size_t i = 0; i < finalNbrs.size(); ++i) {
                  allOrderedNbrs.back().push_back(finalNbrs[i]);
                }
              }
            }

          } // finalNbrs.size() != 1
        } // while (!nbrs.empty())

        if (DEBUG) {
          /*
             cout << "allOrderedNbrs:" << endl;
             for (std::size_t i = 0; i < allOrderedNbrs.size(); ++i) {
             for (std::size_t j = 0; j < allOrderedNbrs[i].size(); ++j) {
             cout << allOrderedNbrs[i][j]->GetIndex() << " ";
             }
             cout << endl;
             }
           */
        }

        for (std::size_t i = 0; i < allOrderedNbrs.size(); ++i) {
          unsigned int lbl = label;
          for (std::size_t j = 0; j < allOrderedNbrs[i].size(); ++j) {
            lbl++;
            code.add(current, allOrderedNbrs[i][j]);
            code.labels[allOrderedNbrs[i][j]->GetIndex()] = lbl;
          }

          CanonicalLabelsRecursive(current, lbl, timeout, bestCode, state);

          for (std::size_t j = 0; j < allOrderedNbrs[i].size(); ++j) {
            code.atoms.pop_back();
            code.bonds.pop_back();
            code.from.pop_back();
            code.labels[allOrderedNbrs[i][j]->GetIndex()] = 0;
          }

        }

      }

    }

    static void CalcCanonicalLabels(OBMol *mol, const std::vector<unsigned int> &symmetry_classes,
        std::vector<unsigned int> &canonical_labels, const OBStereoUnitSet &stereoUnits, const OBBitVec &mask,
        OBStereoFacade *stereoFacade, int maxSeconds, bool onlyOne = false)
    {
      if (!mol->NumAtoms())
        return;
      if (mol->NumAtoms() == 1) {
        canonical_labels.resize(1, 1);
        return;
      }

      canonical_labels.resize(mol->NumAtoms(), 0);

      // find the connected fragments
      OBBitVec visited;
      std::vector<OBBitVec> fragments;
      for (std::size_t i = 0; i < mol->NumAtoms(); ++i) {
        if (!mask.BitIsSet(i+1) || visited.BitIsSet(i+1))
          continue;
        fragments.push_back(getFragment(mol->GetAtom(i+1), mask));
        visited |= fragments.back();
      }

      // precompute the stereo center information
      std::vector<StereoCenter> stereoCenters;
      if (stereoFacade) {
        for (std::size_t i = 0; i < stereoUnits.size(); ++i) {
          const OBStereoUnit &unit = stereoUnits[i];
          if (unit.type == OBStereo::Tetrahedral) {
            if (!stereoFacade->HasTetrahedralStereo(unit.id))
              continue;
            OBTetrahedralStereo::Config config = stereoFacade->GetTetrahedralStereo(unit.id)->GetConfig();
            if (!config.specified)
              continue;
            OBAtom *atom = mol->GetAtomById(config.center);
            if (!atom)
              continue;
            stereoCenters.resize(stereoCenters.size()+1);
            stereoCenters.back().indexes.push_back(atom->GetIndex());
            OBAtom *from = mol->GetAtomById(config.from);
            if (from && !from->IsHydrogen())
              stereoCenters.back().nbrIndexes1.push_back(from->GetIndex());
            else
              stereoCenters.back().nbrIndexes1.push_back(std::numeric_limits<unsigned int>::max());
            for (std::size_t j = 0; j < config.refs.size(); ++j) {
              OBAtom *ref = mol->GetAtomById(config.refs[j]);
              if (ref && !ref->IsHydrogen())
                stereoCenters.back().nbrIndexes1.push_back(ref->GetIndex());
              else
                stereoCenters.back().nbrIndexes1.push_back(std::numeric_limits<unsigned int>::max());
            }
          } else if (unit.type == OBStereo::CisTrans) {
            if (!stereoFacade->HasCisTransStereo(unit.id))
              continue;
            OBCisTransStereo::Config config = stereoFacade->GetCisTransStereo(unit.id)->GetConfig();
            if (!config.specified)
              continue;
            OBAtom *begin = mol->GetAtomById(config.begin);
            OBAtom *end = mol->GetAtomById(config.end);
            if (!begin || !end)
              continue;
            stereoCenters.resize(stereoCenters.size()+1);
            stereoCenters.back().indexes.push_back(begin->GetIndex());
            stereoCenters.back().indexes.push_back(end->GetIndex());
            for (std::size_t j = 0; j < config.refs.size(); ++j) {
              OBAtom *ref = mol->GetAtomById(config.refs[j]);
              unsigned int r = (ref && !ref->IsHydrogen()) ? ref->GetIndex() : std::numeric_limits<unsigned int>::max();
              if (stereoCenters.back().nbrIndexes1.size() < 2)
                stereoCenters.back().nbrIndexes1.push_back(r);
              else
                stereoCenters.back().nbrIndexes2.push_back(r);
            }
          }
        }
      }

      // find the canonical code for each fragment
      std::vector<CanonicalLabelsImpl::FullCode> fcodes;
      for (std::size_t f = 0; f < fragments.size(); ++f) {
        const OBBitVec &fragment = fragments[f];

        // find the maximal symmetry class in the fragment
        std::vector<unsigned int> nextSymClasses;
        for (std::size_t i = 0; i < mol->NumAtoms(); ++i) {
          if (fragment.BitIsSet(i+1))
            nextSymClasses.push_back(symmetry_classes[i]);
        }
        unsigned int maxSymClass = *std::max_element(nextSymClasses.begin(), nextSymClasses.end());

        CanonicalLabelsImpl::Timeout timeout(maxSeconds);
        CanonicalLabelsImpl::FullCode bestCode;
        for (std::size_t i = 0; i < mol->NumAtoms(); ++i) {
          if (!fragment.BitIsSet(i+1))
            continue;
          OBAtom *atom = mol->GetAtom(i+1);

          if (symmetry_classes[atom->GetIndex()] == maxSymClass) {
            // start labeling from the atom with the highest symmetry class
            State state(symmetry_classes, fragment, stereoCenters, onlyOne);
            state.code.add(atom);
            state.code.labels[atom->GetIndex()] = 1;
            CanonicalLabelsRecursive(atom, 1, timeout, bestCode, state);
          }

        }

        if (time(NULL) - timeout.startTime > timeout.maxTime) {
          obErrorLog.ThrowError(__FUNCTION__, "maximum time exceeded...", obError);
        }

        fcodes.push_back(bestCode);
      }

      // sort the codes for the fragments
      std::sort(fcodes.begin(), fcodes.end(), SortCode);

      // construct the full labeling from the sorted fragment labelings
      unsigned int offset = 0;
      for (std::size_t f = 0; f < fcodes.size(); ++f) {
        unsigned int max_label = 0;
        for (std::size_t i = 0; i < mol->NumAtoms(); ++i) {
          if (fcodes[f].labels[i]) {
            canonical_labels[i] = fcodes[f].labels[i] + offset;
            max_label = canonical_labels[i];
          }
        }
        offset += max_label;
      }


    }

  }; // CanonicalLabelsImpl

  void OBGraphSym::CanonicalLabels(std::vector<unsigned int> &canon_labels, int maxSeconds)
  {
    std::vector<unsigned int> symmetry_classes;
    GetSymmetry(symmetry_classes);

    d->CanonicalLabels(symmetry_classes, canon_labels, maxSeconds);
    if (canon_labels.empty())
      for (std::size_t i = 0; i < symmetry_classes.size(); ++i)
        canon_labels.push_back(i+1);
  }

  void OBGraphSym::CanonicalLabels(OBMol *mol, const std::vector<unsigned int> &symmetry_classes,
      std::vector<unsigned int> &canonical_labels, const OBBitVec &mask)
  {
    OBBitVec maskCopy(mask);
    if (!maskCopy.CountBits())
      FOR_ATOMS_OF_MOL (atom, mol)
        maskCopy.SetBitOn(atom->GetIdx());

    CanonicalLabelsImpl::CalcCanonicalLabels(mol, symmetry_classes, canonical_labels, OBStereoUnitSet(), maskCopy, 0, 10, true);
    if (canonical_labels.empty())
      for (std::size_t i = 0; i < symmetry_classes.size(); ++i)
        canonical_labels.push_back(i+1);
  }


  void OBGraphSymPrivate::CanonicalLabels(const std::vector<unsigned int> &symmetry_classes, std::vector<unsigned int> &canon_labels, int maxSeconds)
  {
    bool hasAtLeastOneDefined = false;
    OBStereoFacade sf(_pmol, false);
    FOR_ATOMS_OF_MOL (atom, _pmol) {
      if (sf.HasTetrahedralStereo(atom->GetId())) {
        if (sf.GetTetrahedralStereo(atom->GetId())->GetConfig().specified) {
          hasAtLeastOneDefined = true;
          break;
        }
      }
    }
    FOR_BONDS_OF_MOL (bond, _pmol) {
      if (sf.HasCisTransStereo(bond->GetId())) {
        if (sf.GetCisTransStereo(bond->GetId())->GetConfig().specified) {
          hasAtLeastOneDefined = true;
          break;
        }
      }
    }
    if (!hasAtLeastOneDefined) {
      CanonicalLabelsImpl::CalcCanonicalLabels(_pmol, symmetry_classes, canon_labels, OBStereoUnitSet(), _frag_atoms, 0, maxSeconds);
      return;
    }


    if (!_stereoUnits.size()) {

      _stereoUnits = FindStereogenicUnits(_pmol, symmetry_classes);

      // Mark all invalid stereo data as unspecified
      std::vector<OBGenericData*> stereoData = _pmol->GetAllData(OBGenericDataType::StereoData);
      std::vector<OBGenericData*>::iterator data;
      for (data = stereoData.begin(); data != stereoData.end(); ++data) {
        OBStereo::Type type = ((OBStereoBase*)*data)->GetType();
        if (type == OBStereo::Tetrahedral) {
          OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
          OBTetrahedralStereo::Config config = ts->GetConfig();
          bool valid = true;
          if (!ts->IsValid())
            valid = false;
          OBAtom *center = _pmol->GetAtomById(config.center);
          if (!center)
            valid = false;
          else if (!isTetrahedral(center, _stereoUnits))
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
          OBAtom *beginAtom = _pmol->GetAtomById(config.begin);
          OBAtom *endAtom = _pmol->GetAtomById(config.end);
          if (!beginAtom || !endAtom)
            valid = false;
          else {
            OBBond *bond = _pmol->GetBond(beginAtom, endAtom);
            if (!bond)
              valid = false;
            else if (!isCisTrans(bond, _stereoUnits))
              valid = false;
          }

          if (!valid) {
            config.specified = false;
            ct->SetConfig(config);
          }
        }
      }

    }

    // Determine stereochemistry from coordinates if needed
    if (!_pmol->HasChiralityPerceived()) {
      if (DEBUG)
        cout << "Determining stereochemistry from coords..." << endl;
      switch (_pmol->GetDimension()) {
        case 2:
          _pmol->DeleteData(OBGenericDataType::StereoData);
          TetrahedralFrom2D(_pmol, _stereoUnits);
          CisTransFrom2D(_pmol, _stereoUnits);
          break;
        case 3:
          _pmol->DeleteData(OBGenericDataType::StereoData);
          TetrahedralFrom3D(_pmol, _stereoUnits);
          CisTransFrom3D(_pmol, _stereoUnits);
          break;
        default:
          TetrahedralFrom0D(_pmol, _stereoUnits);
          CisTransFrom0D(_pmol, _stereoUnits);
          break;
      }
    }

    CanonicalLabelsImpl::CalcCanonicalLabels(_pmol, symmetry_classes, canon_labels, _stereoUnits, _frag_atoms, &sf, maxSeconds);
  }


  /**
   * @page canonical_code Canonical Coding Algorithm
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
   * @section Terminology
   * - Topology:
   *
   * @section Topological Symmetry
   * The starting point of the canonical coding algorithm is the topological
   * symmetry or graph symmetry. This symmetry can be expressed by assigning
   * symmetry classes to atoms. Symmetric atoms have the same symmetry class.
   * The term color is often used in graph theory as a synonym for symmetry
   * classes.
   *
   *
   *
   *
   * @subsection Canonical Code for Asymmetric Molecules
   *
   *
   * @subsection Canonical Code for Symmetric Molecules
   *
   * @section Stereochemistry


   For a molecule containing N atoms, there are N! possible
   * ways to order the atoms. The actual complexity of the algorithm needed
   * depends on the molecule.
   *
   *

   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   */

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

//! \file graphsym.cpp
//! \brief XXXX
