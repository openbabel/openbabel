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

using namespace std;

namespace OpenBabel {

#define DEBUG 0

/***************************************************************************
* FUNCTION: print_vector_pairs
*
* DESCRIPTION:
*       For debugging only: Prints a vector of pair<OBAtom, int>.
***************************************************************************/

#if DEBUG
static void print_vector_pairs(char *prefix, vector<pair<OBAtom*,unsigned int> > &vp) {
  vector<pair<OBAtom*,unsigned int> >::iterator vp_iter;
  for (vp_iter = vp.begin(); vp_iter != vp.end(); vp_iter++) {
    OBAtom *atom = vp_iter->first;
    int symclass = vp_iter->second;
    cout << prefix << "atom: " << etab.GetSymbol(atom->GetAtomicNum()) << ", idx: " << atom->GetIdx() << ", symclass: " << symclass << "\n";
  }
}
#endif
  

/***************************************************************************
* FUNCTION: CompareXXX
*
* DESCRIPTION:
*       Three functions for use by the sort() method of a vector.
***************************************************************************/

// NOTE: Copied from OpenBabel/mol.cpp

static bool CompareUnsigned(const unsigned int &a,const unsigned int &b)
{
  return(a<b);
}

static bool ComparePairFirst(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
{
  return(a.first->GetIdx() < b.first->GetIdx());
}

static bool ComparePairSecond(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
{
  return(a.second < b.second);
}

static bool CompareBondPairSecond(const pair<OBBond*,unsigned int> &a,const pair<OBBond*,unsigned int> &b)
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

static unsigned int GetValence(OBAtom *atom, OBBitVec &frag_atoms)
{
  unsigned int count = 0;
  OBBond *bond;
  OBAtom *nbr;

  vector<OBEdgeBase*>::iterator bi;
  for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
    nbr = bond->GetNbrAtom(atom);
    if (frag_atoms.BitIsSet(nbr->GetIdx()))
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

static unsigned int GetHvyValence(OBAtom *atom, OBBitVec &frag_atoms)
{
  unsigned int count = 0;
  OBBond *bond;
  OBAtom *nbr;

  vector<OBEdgeBase*>::iterator bi;
  for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
    nbr = bond->GetNbrAtom(atom);
    if (frag_atoms.BitIsSet(nbr->GetIdx()) && !(nbr->IsHydrogen()))
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

static unsigned int GetHvyBondSum(OBAtom *atom, OBBitVec &frag_atoms)
{
  float count = 0.0;
  OBBond *bond;
  OBAtom *nbr;

  vector<OBEdgeBase*>::iterator bi;
  for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
    nbr = bond->GetNbrAtom(atom);
    if (frag_atoms.BitIsSet(nbr->GetIdx()) && !(nbr->IsHydrogen())) {
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

static bool GetGTDVector(OBMol *pmol,
                         OBBitVec &frag_atoms,
                         vector<int> &gtd)
{
  gtd.clear();
  gtd.resize(pmol->NumAtoms());
  
  int gtdcount, natom;
  OBBitVec used, curr, next;
  OBAtom *atom, *atom1;
  OBBond *bond;
  vector<OBNodeBase*>::iterator ai;
  vector<OBEdgeBase*>::iterator j;

  next.Clear();

  for (atom = pmol->BeginAtom(ai); atom; atom = pmol->NextAtom(ai)) {

    int idx = atom->GetIdx();
    if (!frag_atoms.BitIsOn(idx)) {     // Not in this fragment?
      gtd[idx-1] = 0;
      continue;
    }

    gtdcount = 0;
    used.Clear();curr.Clear();
    used.SetBitOn(idx);
    curr.SetBitOn(idx);

    while (!curr.IsEmpty()) {
      next.Clear();
      for (natom = curr.NextBit(-1);natom != curr.EndBit();natom = curr.NextBit(natom)) {
        atom1 = pmol->GetAtom(natom);
        if (!frag_atoms.BitIsOn(atom1->GetIdx()))
          continue;
        for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j)) {
          int nbr_idx = bond->GetNbrAtomIdx(atom1);
          if (   frag_atoms.BitIsOn(nbr_idx)
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

#if DEBUG
  cout << "    GetGTDVector: ";
  vector<int>::iterator ii;
  for (ii = gtd.begin(); ii != gtd.end(); ii++)
    cout << *ii << " ";
  cout << "\n";
#endif

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

static void FindRingAtoms(OBMol *pmol, OBBitVec &frag_atoms, OBBitVec &ring_atoms)
{
  vector<OBRing*> sssRings;
  vector<OBRing*>::iterator ri;

  ring_atoms.Resize(pmol->NumAtoms());
  ring_atoms.Clear();

  sssRings = pmol->GetSSSR();
  for (ri = sssRings.begin(); ri != sssRings.end(); ri++) {
    OBRing *ring = *ri;
    OBBitVec bvtmp = frag_atoms & ring->_pathset;       // intersection: fragment and ring
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

static void GetGIVector(OBMol *pmol,
                        OBBitVec &frag_atoms,
                        vector<unsigned int> &vid)
{
  // Prepare the vector...
  vid.clear();
  vid.resize(pmol->NumAtoms());

  // The "graph theoretical distance" for each atom (see comments in the function)
  vector<int> v;
  GetGTDVector(pmol, frag_atoms, v);

  // Compute the ring atoms for this particular fragment (set of atoms)
  OBBitVec ring_atoms;
  FindRingAtoms(pmol, frag_atoms, ring_atoms);

  int i;
  OBAtom *atom;
  vector<OBNodeBase*>::iterator ai;
  for (i=0, atom = pmol->BeginAtom(ai); atom; atom = pmol->NextAtom(ai)) {
    vid[i] = 0;
    if (frag_atoms.BitIsOn(atom->GetIdx())) {
      vid[i] = 
        v[i]                                                    // 10 bits: graph-theoretical distance
        | (GetHvyValence(atom,frag_atoms)                <<10)  //  4 bits: heavy valence
        | (((atom->IsAromatic()) ? 1 : 0)                <<14)  //  1 bit:  aromaticity
        | (((ring_atoms.BitIsOn(atom->GetIdx())) ? 1 : 0)<<15)  //  1 bit:  ring atom
        | (atom->GetAtomicNum()                          <<16)  //  7 bits: atomic number
        | (GetHvyBondSum(atom, frag_atoms)               <<23)  //  4 bits: heavy bond sum
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

static void CreateNewClassVector(vector<pair<OBAtom*,unsigned int> > &vp1,
                                 vector<pair<OBAtom*,unsigned int> > &vp2,
                                 OBBitVec &frag_atoms, int natoms)
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
      if (frag_atoms.BitIsOn(idx))
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

static void CountAndRenumberClasses(vector<pair<OBAtom*,unsigned int> > &vp,
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
* FUNCTION: ExtendInvarients
*
* DESCRIPTION:
*       This is the core of symmetry analysis.  Starting with a set of
*       classes on each atom, it "spreads" them using a sum-of-invarients
*       of each atom's class and its neighbors' classes.  This iterates
*       until a stable solution is found (further spreading doesn't
*       change the answer).
*       
* RETURNS: The number of distinct symmetry classes found.
***************************************************************************/

static int ExtendInvarients(vector<pair<OBAtom*, unsigned int> > &symmetry_classes,
                            OBBitVec &frag_atoms, int nfragatoms, int natoms)
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
      CreateNewClassVector(symmetry_classes, tmp_classes, frag_atoms, natoms);
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
*       The bit vector "frag_atoms" specifies a fragment of the molecule,
*       where each bit represents the presence or absence of the atom in
*       the fragment.  Symmetry is computed as though the fragment is the
*       only part that exists.
***************************************************************************/

static int CalculateSymmetry(OBMol *pmol,
                             OBBitVec &frag_atoms,
                             vector<pair<OBAtom*, unsigned int> > &symmetry_classes)
{
  vector<unsigned int> vgi;
  vector<OBNodeBase*>::iterator j;
  OBAtom *atom;

  // How many atoms, and how many do we care about?
  int natoms = pmol->NumAtoms();
  int nfragatoms = frag_atoms.CountBits();

#if DEBUG
  for (atom = pmol->BeginAtom(j); atom; atom = pmol->NextAtom(j)) {
    if (frag_atoms[atom->GetIdx()])
      cout << etab.GetSymbol(atom->GetAtomicNum()) << " ";    
  }
  cout << "\n";
#endif

  // Get vector of graph invariants.  These are the starting "symmetry classes".
  GetGIVector(pmol, frag_atoms, vgi);

  // Create a vector-of-pairs, associating each atom with its Class ID.
  for (atom = pmol->BeginAtom(j); atom; atom = pmol->NextAtom(j)) {
    int idx = atom->GetIdx();
    if (frag_atoms.BitIsOn(idx))
      symmetry_classes.push_back(pair<OBAtom*, unsigned int> (atom, vgi[idx-1]));
  }

  // The heart of the matter: Do extended sum-of-invarients until no further
  // changes are noted. 
  int nclasses = ExtendInvarients(symmetry_classes, frag_atoms, nfragatoms, natoms);

  return nclasses;
}

/***************************************************************************
* FUNCTION: BreakChiralTies
*
* DESCRIPTION:
*       After the achiral symmetry analysis ChiralSymmetry() is done, but
*       before the "tie breaker" step (see CanonicalLabels(), below), there
*       may be two (or more) chiral centers in the same symmetry class that
*       actually have different chirality.  This function finds such chiral
*       centers, and compares their chirality.  If it finds that two atoms
*       in the same symmetry class have the same chirality, it leaves them
*       alone, but if they have opposite chirality, it breaks the tie
*       between the two atoms.
*
*       Actually, it's more subtle than that.  Suppose there are a bunch of
*       chiral atoms in the same symmetry class.  The the class is divided
*       into two classes: All atoms with one chirality go into one class,
*       and all atoms with the opposite chirality go in the other class.
*
* INPUTS:
*       pmol                the molecule
*       frag_atoms          atoms of the molecules in this fragment
*       atom_sym_classes    vector of atom/symclass pairs
*
***************************************************************************/

static void BreakChiralTies(OBMol *pmol,
                            OBBitVec &frag_atoms, int nfragatoms,
                            vector<pair<OBAtom*, unsigned int> > &atom_sym_classes)
{
  vector<pair<OBAtom*,unsigned int> > vp1, vp2;

  // for keeping track of atoms we've already considered
  OBBitVec used_atoms;
  used_atoms.Clear();
  used_atoms.Resize(pmol->NumAtoms()+1);

  // Convert the atom/class pairs to an array indexed by atom idx.
  // This is just for convenience in the next step.  Note that there
  // will be "holes" in this vector since it's a molecule fragment.
  vector<unsigned int> idx2sym_class;
  idx2sym_class.resize(pmol->NumAtoms()+1);
  vector<pair<OBAtom*,unsigned int> >::iterator api;
  for (api = atom_sym_classes.begin(); api < atom_sym_classes.end(); api++)
    idx2sym_class[api->first->GetIdx()] = api->second;


  // Loop over all atoms...
  vector<OBNodeBase*>::iterator ai, aj;
  for (OBAtom *atom = pmol->BeginAtom(ai); atom; atom = pmol->NextAtom(ai)) {
  
    int idx = atom->GetIdx();
    
    // We only want: unused, chiral, and part of this fragment.
    if (!frag_atoms[idx])
      continue;
    if (used_atoms[idx])
      continue;
    used_atoms.SetBitOn(idx);
    if (GetValence(atom, frag_atoms) < 4)       // Valence relative to this fragment
      continue;
    if (!atom->IsChiral())
      continue;

    // Get its symmetry class
    int symclass = idx2sym_class[idx];
    
    // Start the vector of "same class" atoms by adding this atom
    vector<OBAtom *> same_class;
    same_class.push_back(atom);

    // Inner Loop: Find all other atoms with the same symmetry class
    OBAtom *atom2;
    aj = ai;
    for (atom2 = pmol->NextAtom(aj) ; atom2; atom2 = pmol->NextAtom(aj)) {
      int idx2 = atom2->GetIdx();
      if (used_atoms[idx2])
        continue;
      if (idx2sym_class[idx2] == symclass) {
        same_class.push_back(atom2);
        used_atoms.SetBitOn(idx2);
      }
    }
    
    // Unless at least two atoms in the class, there are no ties to break
    if (same_class.size() < 2)
      continue;

#if DEBUG
    cout << "BreakChiralTies: same_class = ";
    vector<OBAtom*>::iterator ia;
    for (ia = same_class.begin(); ia != same_class.end(); ia++)
      cout << (*ia)->GetIdx() << " ";
    cout << "\n";
#endif

    // Create a vector of each atom's neighbors, and sort them
    // (the neighbors) by symmetry class.  When this step is through,
    // we'll be able to compare the 3D coordinates.

    vector<pair<OBAtom*, vector<OBAtom*> > > sorted_neighbors;
    vector<OBAtom*>::iterator iatom;
    for (iatom = same_class.begin(); iatom != same_class.end(); iatom++) {

      // Gather the neighbors, and sort them by symmetry class
      vector<pair<OBAtom*,unsigned int> > neighbors;    // pair: neighbor/symclass
      FOR_NBORS_OF_ATOM(i_nbr, (*iatom)) {
        OBAtom *nbr = &(*i_nbr);
        int idx = nbr->GetIdx();
        if (frag_atoms[idx])
          neighbors.push_back(pair<OBAtom*,unsigned int>(nbr, idx2sym_class[idx]));
      }
      sort(neighbors.begin(), neighbors.end(), ComparePairSecond);

      // Copy the sorted neighbors to a simple vector (we don't care
      // about the symclass any more now that they're sorted).
      vector<OBAtom*> sorted;
      vector<pair<OBAtom *, unsigned int> >::iterator ni;
      for (ni = neighbors.begin(); ni != neighbors.end(); ni++)
        sorted.push_back((*ni).first);

      // Now put the chiral atom, and its sorted neighbors, into a vector.
      sorted_neighbors.push_back(pair<OBAtom*, vector<OBAtom*> >(*iatom, sorted));
    }
    
    // See the comments in cansmilesformat.cpp::GetChiralStereo regarding
    // "torsion".  If two atoms' neighbors, sorted the same way, have torsion
    // angles with the same sign (both positive or both negative torsions),
    // then their chiralities match, otherwise they're opposite chiralities.
    // Using this fact, divide the atoms in this symmetry class into two groups.

    OBAtom *ref_atom              = sorted_neighbors[0].first;
    vector<OBAtom*> ref_neighbors = sorted_neighbors[0].second;

    double t1 = CalcTorsionAngle(ref_neighbors[0]->GetVector(),
                                 ref_neighbors[1]->GetVector(),
                                 ref_neighbors[2]->GetVector(),
                                 ref_neighbors[3]->GetVector());

    vector<OBAtom*> symclass1, symclass2;
    symclass1.push_back(ref_atom);

    for (int i = 1; i < sorted_neighbors.size(); i++) {
      OBAtom *atom             = sorted_neighbors[i].first;
      vector<OBAtom*>neighbors = sorted_neighbors[i].second;
      double t2 = CalcTorsionAngle(neighbors[0]->GetVector(),
                                   neighbors[1]->GetVector(),
                                   neighbors[2]->GetVector(),
                                   neighbors[3]->GetVector());

      if (t1*t2 >= 0.0)
        symclass1.push_back(atom);      // t1 & t2 have same signs ==> same chirality
      else
        symclass2.push_back(atom);      // t1 & t2 have opposite signs ==> opposite chirality
    }

    // If there's nothing in symclass2, then we don't have to split 
    // the symmetry class.
    if (symclass2.empty())
      continue;

#if DEBUG
    cout << "Time to break chiral ties...\n";
#endif

    // Time to break the class in two. Double all symmetry classes
    // Then, either add one or subtract one to all the atoms in symclass2,
    // the atoms that didn't match the chirality of the "reference" class above.
    // The add-or-subtract decision is based on t1, the "torsion angle" calculated
    // above; by doing this, we force a consistent choice of '@' or '@@' for the,
    // two chiral centers; otherwise it's arbitrary and you'll get two different SMILES.
    for (int i = 0; i < atom_sym_classes.size(); i++) {
      atom_sym_classes[i].second *= 2;
      for (int j = 0; j < symclass2.size(); j++) {
        if (symclass2[j] == atom_sym_classes[i].first) {
          if (t1 > 0)
            atom_sym_classes[i].second += 1;
          else
            atom_sym_classes[i].second -= 1;
        }
      }
    }
    
    // Now propagate the change across the whole molecule with the
    // extended sum-of-invariants.
    ExtendInvarients(atom_sym_classes, frag_atoms, nfragatoms, pmol->NumAtoms());
    
  }

}

/***************************************************************************
* FUNCTION: FindConjugatedEZBonds
*
* DESCRIPTION:
*       When generating a canonical SMILES for an E/Z double bond, the
*       system needs to choose between two equivalent representations,
*       e.g. C/C=C/C vs. C\C=C\C, and C/C=C\C vs. C\C=C/C.  The function
*       FixCisTransBonds() below does this by selecting the lowest atom
*       in the canonical ordering, and if it's going to be "\", it "flips"
*       all of the Up/Down designations around the double bond so that
*       the first symbol printed is always "/".
*
*       However, it's trickier than that.  Consider:
*
*               C\C=C/C=C/C
*
*       If you just flip the two up/down bonds around the first double
*       bond, you'll get the second double bond wrong.  So to flip the
*       first one, you need to identify all of the conjugated double bonds,
*       and flip them as a set.
*
*       This function "walks" outward from a double bond and identifies
*       all neighboring double bonds, continuing until either a single
*       bond is encountered with no up/down designation, or the neighboring
*       bond isn't a double bond.
*       
*       On return, "flip_bonds" will be a set of single bonds that must
*       all be "flipped" (Up/Down reversed) as a set.
***************************************************************************/

static void FindConjugatedEZBonds(OBAtom *atom,
                                  OBBitVec &flip_bonds,
                                  OBBitVec &visited_atoms)
{
  visited_atoms.SetBitOn(atom->GetIdx());

  FOR_BONDS_OF_ATOM(bi, atom) {

    OBBond *bond = &(*bi);

    // Filter out bonds that can't be part of 
    // conjugated double-bond chain
    if (!bond->IsSingle())                      // only care about single-bonded neighbors
      continue;
    if (!bond->IsUp() && !bond->IsDown())       // bond must be marked up or down
      continue;
    OBAtom *nbr = bond->GetNbrAtom(atom);
    if (visited_atoms[nbr->GetIdx()])           // been there already?
      continue;

    // Ok, it's marked and it's single: this bond is part of the "flip set"
    flip_bonds.SetBitOn(bond->GetIdx());

    // Now see if the neighbor atom also has a double bond, and if so 
    // recursively carry the markings forward.
    if (nbr->HasDoubleBond())                   // neighbor must also have double bond
      FindConjugatedEZBonds(nbr, flip_bonds, visited_atoms);
  }

  // Now recursively call this function on the other atom across
  // the double bond from this one.
  FOR_BONDS_OF_ATOM(bond, atom) {
    if (bond->IsDouble()) {
      OBAtom *nbr = bond->GetNbrAtom(atom);
      if (visited_atoms[nbr->GetIdx()])
        break;
      FindConjugatedEZBonds(nbr, flip_bonds, visited_atoms);
      break;
    }
  }

#if DEBUG
  cout << "FindConjugatedEZBonds: atom " << atom->GetIdx() << ": " << flip_bonds << "\n";
#endif
}


/***************************************************************************
* FUNCTION: FixCisTransBonds
*
* DESCRIPTION:
*
*       NOTE: THIS ALGORITHM IS NOT SUFFICIENT.  It canonicalizes certain
*       aspects of cis/trans chirality, but it misses some important 
*       cases.  In particular, all symmetry analysis is done without taking
*       cis/trans information into account; only after the symmetry analysis
*       is done that the cis/trans double bonds are "canonicalized" by either
*       erasing the cis/trans bonds (when symmetry makes it not cis/trans), or
*       by rearranging the up/down bonds into a canonical configuration.
*
*       Here's an example of a molecule that won't canonicalize properly:
*
*
*              X   Y           The canonicalizer will think this is
*               \ /            symmetrical because it doesn't take the
*         A      D             double bonds into account when determining
*          \    / \            the chirality of D.
*           B==C   E==F
*                      \
*                       G
*
*       To fix this, we need to do a "tie breaker" between A and G using their
*       cis/trans information, in the same fashion that we do tie-breaking of
*       chiral centers, above.
*
*       Other than that, here are the normal comments...
*
*       This function fixes two problems in a molecule:
*
*       1. Many molecules have cis/trans markings that are bogus - the
*          molecule is actually symmetrical at one or both ends of the
*          double bond.
*
*       2. Many molecules have only two or three of the cis/trans bonds
*          marked as "up" and "down" bonds, when all of them should be
*          marked.
*
*       Using the symmetry classes, this function fixes these problems.
***************************************************************************/

static void FixCisTransBonds(OBMol *pmol,
                             OBBitVec &frag_atoms,
                             vector<unsigned int> &symmetry_classes,
                             vector<unsigned int> &canonical_labels)
{
  // Bit vector records which bonds are discovered to actually be up/down.
  // Some bonds may be marked up/down but after symmetry analysis we learn
  // that they are not.
  OBBitVec cis_trans_bonds(pmol->NumBonds() + 1);
  cis_trans_bonds.Clear();

  FOR_BONDS_OF_MOL(dbi, pmol) {

    OBBond *dbl_bond = &(*dbi);

    // Not a double bond?
    if (!dbl_bond->IsDouble() || dbl_bond->IsAromatic())
      continue;

    bool is_cis_trans = true;

    // Find the single bonds around the atoms connected by the double bond.
    // While we're at it, note whether the pair of atoms on either end are
    // identical, in which case it's not cis/trans.

    OBAtom *a1 = dbl_bond->GetBeginAtom();
    OBAtom *a2 = dbl_bond->GetEndAtom();

#if DEBUG
    cout << "FixCisTransBonds: dbl_bond = " << dbl_bond->GetIdx() << ", atoms = " << a1->GetIdx() << ", " << a2->GetIdx() << ")\n";
    cout << "          canonical labels = " << canonical_labels[a1->GetIdx()-1] << ", " << canonical_labels[a2->GetIdx()-1] << "\n";
#endif


    // neighbors of atom1 and atom2
    OBBond *a1_b1 = NULL, *a1_b2 = NULL, *a2_b1 = NULL, *a2_b2 = NULL;

    int symclass = -1;
    FOR_BONDS_OF_ATOM(bi, a1) {
      OBBond *b = &(*bi);
      if ((b) == (dbl_bond)) continue;  // skip the double bond we're working on
      if (b->IsSingle() || b->IsAromatic()) {
        if (symclass >= 0 && symclass == symmetry_classes[b->GetNbrAtom(a1)->GetIdx()-1]) {
          is_cis_trans = false;                 // two atoms, and both atoms the same
        }
        symclass = symmetry_classes[b->GetNbrAtom(a1)->GetIdx()-1];
        if (NULL == a1_b1)
          a1_b1 = b;    // remember 1st bond of Atom1
        else
          a1_b2 = b;    // remember 2nd bond of Atom1
      } else {
        is_cis_trans = false;           // it's not a single bond, can't be cis/trans
      }
    }

    symclass = -1;
    FOR_BONDS_OF_ATOM(bi, a2) {
      OBBond *b = &(*bi);
      if (b == dbl_bond) continue;
      if (b->IsSingle() || b->IsAromatic()) {
        if (symclass >= 0 && symclass == symmetry_classes[b->GetNbrAtom(a2)->GetIdx()-1]) {
          is_cis_trans = false;         // second connection, but both atoms the same
        }
        symclass = symmetry_classes[b->GetNbrAtom(a2)->GetIdx()-1];
        if (NULL == a2_b1)
          a2_b1 = b;    // remember 1st bond of Atom2
        else
          a2_b2 = b;    // remember 2nd bond of Atom2
      } else {
        is_cis_trans = false;           // it's not a single bond, can't be cis/trans
      }
    }

    // Check that both atoms on the double bond have at least one
    // other neighbor;
    int v1 = GetValence(a1, frag_atoms);
    int v2 = GetValence(a2, frag_atoms);
    if (v1 < 2 || v1 > 3 || v2 < 2 || v2 > 3) {
      is_cis_trans = false;
    }

    // Now check that at least two are marked cis/trans.
    if (is_cis_trans) {
      int count = 0;
      if (a1_b1 && (a1_b1->IsUp() || a1_b1->IsDown())) count++;
      if (a1_b2 && (a1_b2->IsUp() || a1_b2->IsDown())) count++;
      if (a2_b1 && (a2_b1->IsUp() || a2_b1->IsDown())) count++;
      if (a2_b2 && (a2_b2->IsUp() || a2_b2->IsDown())) count++;
      if (count < 2) {
        is_cis_trans = false;
      }
    }


    if (is_cis_trans) {

      // It's cis/trans, so record the fact for each bond.  Later, we'll clear
      // out any Up/Down bonds that aren't in this cis_trans_bonds bitmap.
      if (a1_b1) cis_trans_bonds.SetBitOn(a1_b1->GetIdx());
      if (a1_b2) cis_trans_bonds.SetBitOn(a1_b2->GetIdx());
      if (a2_b1) cis_trans_bonds.SetBitOn(a2_b1->GetIdx());
      if (a2_b2) cis_trans_bonds.SetBitOn(a2_b2->GetIdx());

      // Now we know we have a double bond with two, three, or four
      // single-bond neighbors, of which at least two are marked as cis/trans.
      // In this case, make sure all bonds have appropriate up/down markings,
      // and change the up/down marks to a canonical form.

#define CT_MATCH(b1, b2) \
  b1->IsUp() ? b2->SetDown() : (b1->IsDown() ? b2->SetUp() : (b2->IsUp() ? b1->SetDown() : (b2->IsDown() ? b1->SetUp() : b1->SetDown())))

#if DEBUG
      cout << "  Is cis/trans - matching: \n";
#endif

      // Make sure all bonds are actually marked
      if (a1_b1 && a1_b2) {CT_MATCH(a1_b1, a1_b2);}
      if (a2_b1 && a2_b2) {CT_MATCH(a2_b1, a2_b2);}

      // Find the neighbor atom with the lowest canonical label
      vector<pair<OBBond*, int> > bsym;
      if (a1_b1) bsym.push_back(pair<OBBond*, int>(a1_b1, canonical_labels[a1_b1->GetNbrAtom(a1)->GetIdx() - 1]));
      if (a1_b2) bsym.push_back(pair<OBBond*, int>(a1_b2, canonical_labels[a1_b2->GetNbrAtom(a1)->GetIdx() - 1]));
      if (a2_b1) bsym.push_back(pair<OBBond*, int>(a2_b1, canonical_labels[a2_b1->GetNbrAtom(a2)->GetIdx() - 1]));
      if (a2_b2) bsym.push_back(pair<OBBond*, int>(a2_b2, canonical_labels[a2_b2->GetNbrAtom(a2)->GetIdx() - 1]));
      sort(bsym.begin(), bsym.end(), CompareBondPairSecond);


      // If the lowest-symmetry-class bond is Down, flip everything.
      if (bsym[0].first->IsUp()) {

        // Find the set of conjugated double/single bonds that all
        // need to be flipped as a unit.
        OBBitVec flip_bonds(pmol->NumBonds() + 1);
        OBBitVec visited_atoms(pmol->NumAtoms() + 1);
        FindConjugatedEZBonds(a1, flip_bonds, visited_atoms);

        // Flip them all.
        FOR_BONDS_OF_MOL(bi, pmol) {
          OBBond *b = &(*bi);
          if (flip_bonds[b->GetIdx()]) {
            if (b->IsUp())
              b->SetDown();
            else if (b->IsDown())
              b->SetUp();
          }
        }
      }
    }
  }

  // OK, we've determined which bonds really are cis/trans, and we've
  // got the up/down right on all of them.  The final step is to 
  // remove any Up/Down marks on bonds that we determined weren't 
  // really up or down.

  FOR_BONDS_OF_MOL(bi, pmol) {
    OBBond *b = &(*bi);
    int idx = b->GetIdx();
    if (!cis_trans_bonds[idx]) {
      b->UnsetUp();
      b->UnsetDown();
    }
  }
  
  
}

/***************************************************************************
* FUNCTION: CanonicalLabels
*
* DESCRIPTION:
*       Generates a canonical labeling of the atoms of a molecule, and as
*       a side benefit, returns the symmetry classes of the molecule.
*
*       To create a canonical labeling, we need every node to have a unique
*       label.  The canonical symmetry classes (see CalculateSymmetry(),
*       above) are a good start, but the atoms in each symmetry class are
*       still indistinguishable.  For writing a canonical string, we need
*       to create an arbitrary, but canonical (repeatable) distinction
*       between the atoms in each symmetry class -- "break the ties" in the
*       symmetry values.
*
*       To break ties, we sort into symetry-class order, double all class
*       IDs, then arbitrarily subtract one from the first repeated symmetry
*       class, thus breaking the tie (see Weininger et al).  With this new
*       set of symmetry classes, we repeat the extended-connectivity sums
*       to "spread" the broken symmetry class, and check again.  This is
*       repeated until all symmetry is gone and every atom has a unique
*       label.
*
* RETURNS:
*       symmetry_classes - A vector, indexed by [ OBAtom::GetIdx() - 1].
*       canonical_labels - A vector, indexed by [ OBAtom::GetIdx() - 1].
***************************************************************************/

void CanonicalLabels(OBMol *pmol,
                     OBBitVec &frag_atoms,
                     vector<unsigned int> &symmetry_classes,
                     vector<unsigned int> &canonical_labels)    // on input: symclasses
{
  vector<pair<OBAtom*,unsigned int> > atom_sym_classes, vp1, vp2;
  vector<OBNodeBase*>::iterator j;
  unsigned int nclass1, nclass2; //number of classes
  int i;

  int nfragatoms = frag_atoms.CountBits();
  int natoms = pmol->NumAtoms();

  // Start by calculating symmetry classes, without accounting for chiral centers
  nclass1 = CalculateSymmetry(pmol, frag_atoms, atom_sym_classes);

#if DEBUG
  cout << "BEFORE TieBreaker: nclass1 = " << nclass1 << ", nfragatoms = " << nfragatoms << "\n";
  print_vector_pairs("    ", atom_sym_classes);
#endif

  // Deal with chiral symmetry
  BreakChiralTies(pmol, frag_atoms, nfragatoms, atom_sym_classes);

  // The symmetry classes are the starting point for the canonical labels
  vp1 = atom_sym_classes;

  if (nclass1 < nfragatoms) {
    int tie_broken = 1;
    while (tie_broken) {
      tie_broken = 0;
      int last_rank = -1;
      for (i = 0; i < vp1.size(); i++) {
        vp1[i].second *= 2;             // Double symmetry classes
        if (vp1[i].second == last_rank && !tie_broken) {
          vp1[i-1].second -= 1;         // Break a tie
          tie_broken = 1;
        }
        last_rank = vp1[i].second;
      }
      if (tie_broken) {
        for (i = 0; i < 100;i++) {  //sanity check - shouldn't ever hit this number
          CreateNewClassVector(vp1, vp2, frag_atoms, natoms);
          CountAndRenumberClasses(vp2, nclass2);
          vp1 = vp2;
          if (nclass1 == nclass2) break;
          nclass1 = nclass2;
        }
      } else {
        CountAndRenumberClasses(vp1, nclass1);  // no more ties - undo the doublings
      }
    }
  }

#if DEBUG
  cout << "AFTER TieBreaker: nclass1 = " << nclass1 << ", nfragatoms = " << nfragatoms << "\n";
  print_vector_pairs("    ", vp1);
#endif

  // For return values, convert vectors of atom/int pairs into one-dimensional
  // vectors of int, indexed by atom->GetIdx().
  //
  // Since we're working with a molecular fragment, the symmetry-class and
  // canonical-label vectors are both shorter than the number of atoms in
  // the molecule.  To fix it, we append all of the non-fragment atoms to
  // the end, with a very large value, then, re-sort the vector by the
  // atom's GetIdx() number, and finally copy the result so that the vector
  // corresponds to the atoms' order in the molecule.

#define NO_SYMCLASS 0x7FFFFFFF

  // Add non-fragment atoms so vector is same length as natoms in molecule
  for (OBAtom *atom = pmol->BeginAtom(j); atom; atom = pmol->NextAtom(j))
    if (!frag_atoms.BitIsOn(atom->GetIdx())) {
      atom_sym_classes.push_back(pair<OBAtom*,unsigned int> (atom, NO_SYMCLASS));
      vp1.push_back(pair<OBAtom*,unsigned int> (atom, NO_SYMCLASS));
    }

  vector<pair<OBAtom*,unsigned int> >::iterator k;

  // Sort and copy symmetry classes
  symmetry_classes.clear();
  sort(atom_sym_classes.begin(),atom_sym_classes.end(),ComparePairFirst);
  for (k = atom_sym_classes.begin();k != atom_sym_classes.end();k++)
    symmetry_classes.push_back(k->second);

  // Sort and copy canonical labels
  canonical_labels.clear();
  sort(vp1.begin(),vp1.end(),ComparePairFirst);
  for (k = vp1.begin();k != vp1.end();k++)
    canonical_labels.push_back(k->second);

  // Deal with cis/trans double bonds.
  FixCisTransBonds(pmol, frag_atoms, symmetry_classes, canonical_labels);

}

} // namespace OpenBabel

//! \file canon.cpp
//! \brief Canonical numbering of SMILES, molecules and fragments
