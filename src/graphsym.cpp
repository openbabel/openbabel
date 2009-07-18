/* -*-C++-*-

**********************************************************************
Copyright (C) 2005-2006, eMolecules, Inc. (www.emolecules.com)
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************

+======================================================================
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

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/canon.h>
#include <openbabel/graphsym.h>

using namespace std;

namespace OpenBabel {

/***************************************************************************
* FUNCTION: CompareXXX
*
* DESCRIPTION:
*       Three functions for use by the sort() method of a vector.
***************************************************************************/

// NOTE: Copied from OpenBabel/mol.cpp

  bool OBGraphSym::CompareUnsigned(const unsigned int &a,const unsigned int &b)
  {
    return(a<b);
  }

  bool OBGraphSym::ComparePairFirst(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
  {
    return(a.first->GetIdx() < b.first->GetIdx());
  }

  bool OBGraphSym::ComparePairSecond(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
  {
    return(a.second < b.second);
  }

  bool OBGraphSym::CompareBondPairSecond(const pair<OBBond*,unsigned int> &a,const pair<OBBond*,unsigned int> &b)
  {
    return(a.second < b.second);
  }


/***************************************************************************
* FUNCTION: GetValence
*
* DESCRIPTION:
*       Like OBAtom::GetValence(): Counts the number of neighbors, but
*       doesn't count atoms not in the fragment.
***************************************************************************/

  unsigned int OBGraphSym::GetValence(OBAtom *atom)
  {
    unsigned int count = 0;
    OBBond *bond;
    OBAtom *nbr;

    vector<OBEdgeBase*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms->BitIsSet(nbr->GetIdx()))
        count++;
    }
    return(count);
  }

/***************************************************************************
* FUNCTION: GetHvyValence
*
* DESCRIPTION:
*       Like OBAtom::GetHvyValence(): Counts the number non-hydrogen
*       neighbors, but doesn't count atoms not in the fragment.
***************************************************************************/

  unsigned int OBGraphSym::GetHvyValence(OBAtom *atom)
  {
    unsigned int count = 0;
    OBBond *bond;
    OBAtom *nbr;

    vector<OBEdgeBase*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms->BitIsSet(nbr->GetIdx()) && !(nbr->IsHydrogen()))
        count++;
    }
    
    return(count);
  }

/***************************************************************************
* FUNCTION: GetHvyBondSum
*
* DESCRIPTION:
*       Sums the bond order over the bonds from this atom to other atoms
*       in the fragment.  Single = 1, double = 2, triple = 3, aromatic = 1.6,
*       but sum is rounded to nearest integer.
*
*       This is used for fragment symmetry perception instead of the "implicit
*       valence" used by the standard OpenBabel symmetry perception.  It
*       has the same effect, but we don't have to worry about hydrogen counts,
*       EXCEPT for aromatic N, where the difference between n and [nH] is
*       critical.
***************************************************************************/

  unsigned int OBGraphSym::GetHvyBondSum(OBAtom *atom)
  {
    float count = 0.0;
    OBBond *bond;
    OBAtom *nbr;

    vector<OBEdgeBase*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms->BitIsSet(nbr->GetIdx()) && !(nbr->IsHydrogen())) {
        if (bond->IsSingle())        count += 1.0;
        else if (bond->IsDouble())   count += 2.0;
        else if (bond->IsTriple())   count += 3.0;
        else if (bond->IsAromatic()) count += 1.6;
      }
    }
    if (atom->GetAtomicNum() == 7 && atom->IsAromatic() && atom->GetImplicitValence() == 3) {
      count += 1;         // [nH] - add another bond
    }
    return(int(count + 0.5));     // round to nearest int
  }


/***************************************************************************
* FUNCTION: GetGTDVector
*
* DESCRIPTION:
*       
*       Calculates the graph theoretical distance of each atom.
*       Vector is indexed from zero.
*
*       NOTE: Copied from OpenBabel/mol.cpp
*
*       NOTE: "Indexed from zero" means it's one off from the atom->GetIdx()
*       that's used to index atoms inside the molecule!
*
*       NOTE: This function is hard to decipher, and seems to be misnamed.
*       A "distance" should be be between two atoms, but there's more here
*       than that.  It seems to be doing a breadth-first search to find the
*       most-distant atom from each atom, and reporting the number of steps
*       (which happens to be the graph-theoretical distance) to that atom.
*       The name "Graph Theoretical Distance" is thus misleading.
***************************************************************************/

  bool OBGraphSym::CalcGTDVector()
{
  _gtd.clear();
  _gtd.resize(_pmol->NumAtoms());
  
  int gtdcount, natom;
  OBBitVec used, curr, next;
  OBAtom *atom, *atom1;
  OBBond *bond;
  vector<OBNodeBase*>::iterator ai;
  vector<OBEdgeBase*>::iterator j;

  next.Clear();

  for (atom = _pmol->BeginAtom(ai); atom; atom = _pmol->NextAtom(ai)) {

    int idx = atom->GetIdx();
    if (!_frag_atoms->BitIsOn(idx)) {     // Not in this fragment?
      _gtd[idx-1] = 0;
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
        if (!_frag_atoms->BitIsOn(atom1->GetIdx()))
          continue;
        for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j)) {
          int nbr_idx = bond->GetNbrAtomIdx(atom1);
          if (   _frag_atoms->BitIsOn(nbr_idx)
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
    _gtd[idx-1] = gtdcount;
  }

  return(true);
}

/***************************************************************************
* FUNCTION: FindRingAtoms
*
* DESCRIPTION:
*       Finds all atoms that are part of a ring in the current fragment.
*       We start with the whole molecule's rings, and eliminate any that
*       have atoms not in the subset.  For the rings that are left, mark
*       each atom of the ring as a ring atom.
*
*       Returns a bit vector where TRUE means it's a ring atom.
***************************************************************************/

  void OBGraphSym::FindRingAtoms(OBBitVec &ring_atoms)
{
  vector<OBRing*> sssRings;
  vector<OBRing*>::iterator ri;

  ring_atoms.Resize(_pmol->NumAtoms());
  ring_atoms.Clear();

  sssRings = _pmol->GetSSSR();
  for (ri = sssRings.begin(); ri != sssRings.end(); ri++) {
    OBRing *ring = *ri;
    OBBitVec bvtmp = *_frag_atoms & ring->_pathset;       // intersection: fragment and ring
    if (bvtmp == ring->_pathset)                        // all ring atoms in fragment?
      ring_atoms |= ring->_pathset;                     //   yes - add this ring's atoms 
  }
}


/***************************************************************************
* FUNCTION: GetGIVector
*
* DESCRIPTION:
*       NOTE: Copied from OpenBabel/mol.cpp
*
*       Calculates a set of graph invariant indexes using the graph theoretical
*       distance, number of connected heavy atoms, aromatic boolean, ring
*       boolean, atomic number, and summation of bond orders connected to the
*       atom.
*
*       We have to recalculate which atoms are in rings by taking the fragment's
*       atoms into account when we generate the graph invarients.
*
*       Vector is indexed from zero (not one, like atom->GetIdx()).
*
*       NOTE: This may need to be extended to include the bond-invariant properties,
*       particularly the size of all rings the bond is in (from a SSSR).
***************************************************************************/

  void OBGraphSym::GetGIVector(vector<unsigned int> &vid)
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
    vid[i] = 0;
    if (_frag_atoms->BitIsOn(atom->GetIdx())) {
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

#if DEBUG
  cout << "    GetGIVector:  ";
  vector<unsigned int>::iterator ii;
  for (ii = vid.begin(); ii != vid.end(); ii++)
    cout << *ii << " ";
  cout << "\n";
#endif

}


/***************************************************************************
* FUNCTION: CreateNewClassVector
*
* DESCRIPTION:
*       NOTE: Derived from OpenBabel/mol.cpp
*
*       Creates a new vector of symmetry classes based on an existing
*       vector.  (Helper routine to GetGIDVector.)  On return, vp2 will
*       have newly-extended connectivity sums, but the numbers (the class
*       IDs) are very large.
*
*       (Comments by CJ) This appears to compute the "extended connectivity
*       sums" similar to those described by Weininger, Morgan, etc. It uses
*       vp1 as its starting point (the current connectivity sums), and puts
*       the new sums in vp2.  Note that vp1 is modified along the way.
*
*       Note that, per Weininger's warning, this assumes the initial class
*       ID's are less than 100, which is a BAD assumption, e.g. OCC...CCN
*       would have more than 100 symmetry classes if the chain is more than
*       98 carbons long.  Should change this to use Weininger's product of
*       corresponding primes.
***************************************************************************/

  void OBGraphSym::CreateNewClassVector(vector<pair<OBAtom*,unsigned int> > &vp1,
                                        vector<pair<OBAtom*,unsigned int> > &vp2,
                                        int natoms)
{
  int m,id;
  OBAtom *atom, *nbr;
  vector<OBEdgeBase*>::iterator nbr_iter;
  vector<unsigned int>::iterator k;
  vector<pair<OBAtom*,unsigned int> >::iterator vp_iter;

#if DEBUG
  cout << "CreateNewClassVector: START\n";
  print_vector_pairs("    ", vp1);
#endif

  // There may be fewer atoms than in the whole molecule, so we can't
  // index the vp1 array by atom->GetIdx().  Instead, create a quick
  // mapping vector of idx-to-index for vp1.
  vector<int> idx2index(natoms+1, -1);  // natoms + 1
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
      if (_frag_atoms->BitIsOn(idx))
        vtmp.push_back(vp1[idx2index[idx]].second);
    }

    sort(vtmp.begin(),vtmp.end(),CompareUnsigned);
    for (m = 100, k = vtmp.begin(); k != vtmp.end(); k++, m*=100) 
      id += *k * m;
    vp2.push_back(pair<OBAtom*,unsigned int> (atom, id));
  }
#if DEBUG
  cout << "CreateNewClassVector: FINISH\n";
  print_vector_pairs("    ", vp2);
#endif

}

/***************************************************************************
* FUNCTION: CountAndRenumberClasses
*
* DESCRIPTION:
*       NOTE: Copied from OpenBabel/mol.cpp
*
*       Counts the number of unique symmetry classes in a list.
*
*       (NOTE: CJ -- It also appears to MODIFY the list.  It sorts it in order
*       of class ID, then renumbers the ID's zero through N-1.  See the comments
*       in CreateNewClassVector() about how it returns very large numbers for the
*       class IDs it creates.  These are replaced by lower, sequential numbers here.)
***************************************************************************/

  void OBGraphSym::CountAndRenumberClasses(vector<pair<OBAtom*,unsigned int> > &vp,
                                           unsigned int &count)
{
  count = 1;
  vector<pair<OBAtom*,unsigned int> >::iterator k;

  sort(vp.begin(),vp.end(),ComparePairSecond);
  k = vp.begin();
  if (k != vp.end()) {
    unsigned int id = k->second;
    k->second = 1;
    ++k;
    for (;k != vp.end(); ++k) {
      if (k->second != id) {
        id = k->second;
        k->second = ++count;
      }
      else
        k->second = count;
    }
  }
}

/***************************************************************************
* FUNCTION: ExtendInvariants
*
* DESCRIPTION:
*       This is the core of symmetry analysis.  Starting with a set of
*       classes on each atom, it "spreads" them using a sum-of-invariants
*       of each atom's class and its neighbors' classes.  This iterates
*       until a stable solution is found (further spreading doesn't
*       change the answer).
*       
* RETURNS: The number of distinct symmetry classes found.
***************************************************************************/

  int OBGraphSym::ExtendInvariants(vector<pair<OBAtom*, unsigned int> > &symmetry_classes,
                            int nfragatoms, int natoms)
{
  unsigned int nclasses1, nclasses2;
  vector<pair<OBAtom*,unsigned int> > tmp_classes;

  // How many classes are we starting with?  (The "renumber" part isn't relevant.)
  CountAndRenumberClasses(symmetry_classes, nclasses1);
        
  // LOOP: Do extended sum-of-invarients until no further changes are
  // noted.  (Note: This is inefficient, as it re-computes extended sums
  // and re-sorts the entire list each time.  You can save a lot of time by
  // only recomputing and resorting within regions where there is a tie
  // initially.  But it's a lot more code.)

  if (nclasses1 < nfragatoms) {
    for (int i = 0; i < 100;i++) {  //sanity check - shouldn't ever hit this number
      CreateNewClassVector(symmetry_classes, tmp_classes, natoms);
      CountAndRenumberClasses(tmp_classes, nclasses2);
      symmetry_classes = tmp_classes;
      if (nclasses1 == nclasses2) break;
      nclasses1 = nclasses2;
    }
  }
  return nclasses1;
}

/***************************************************************************
* FUNCTION: CalculateSymmetry
*
* DESCRIPTION:
*       Calculates a set of canonical symmetry identifiers for a molecule.
*       Atoms with the same symmetry ID are symmetrically equivalent.  By
*       "canonical", we mean it generates a repeatable labelling of the
*       atoms, i.e. the same fragment will get the same symmetry labels in
*       any molecule in which it occurs.
*
*       Vector is indexed from zero, corresponding to (atom->GetIdx() - 1).
*
*       The bit vector "_frag_atoms" specifies a fragment of the molecule,
*       where each bit represents the presence or absence of the atom in
*       the fragment.  Symmetry is computed as though the fragment is the
*       only part that exists.
***************************************************************************/

  int OBGraphSym::CalculateSymmetry(vector<pair<OBAtom*, unsigned int> > &symmetry_classes)
{
  vector<unsigned int> vgi;
  vector<OBNodeBase*>::iterator j;
  OBAtom *atom;

  // How many atoms, and how many do we care about?
  int natoms = _pmol->NumAtoms();
  int nfragatoms = _frag_atoms->CountBits();

#if DEBUG
  for (atom = _pmol->BeginAtom(j); atom; atom = _pmol->NextAtom(j)) {
    if (_frag_atoms[atom->GetIdx()])
      cout << etab.GetSymbol(atom->GetAtomicNum()) << " ";    
  }
  cout << "\n";
#endif

  // Get vector of graph invariants.  These are the starting "symmetry classes".
  GetGIVector(vgi);

  // Create a vector-of-pairs, associating each atom with its Class ID.
  for (atom = _pmol->BeginAtom(j); atom; atom = _pmol->NextAtom(j)) {
    int idx = atom->GetIdx();
    if (_frag_atoms->BitIsOn(idx))
      symmetry_classes.push_back(pair<OBAtom*, unsigned int> (atom, vgi[idx-1]));
  }

  // The heart of the matter: Do extended sum-of-invariants until no further
  // changes are noted. 
  int nclasses = ExtendInvariants(symmetry_classes, nfragatoms, natoms);

  return nclasses;
}

} // namespace OpenBabel

//! \file graphsym.cpp
//! \brief XXXX
