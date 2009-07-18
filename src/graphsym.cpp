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
#include <openbabel/graphsym.h>

#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>

using namespace std;

namespace OpenBabel {

/***************************************************************************
* FUNCTION: CompareXXX
*
* DESCRIPTION:
*       Three functions for use by the sort() method of a vector.
***************************************************************************/

// NOTE: Copied from OpenBabel/mol.cpp

  // Destructor
  OBGraphSym::~OBGraphSym()
  {
    // Nothing to free as we only hold pointers
  }

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

  bool OBGraphSym::GetGTDVector(vector<int> &gtd)
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
    if (!_frag_atoms->BitIsOn(idx)) {     // Not in this fragment?
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
    gtd[idx-1] = gtdcount;
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

void IdsToSymClasses(OBMol *mol, OBCisTransStereo::Config &config, 
    const std::vector<unsigned int> &symClasses)
{
  OBAtom *atom;
  // begin
  atom = mol->GetAtomById(config.begin);
  if (atom)
    config.begin = symClasses.at(atom->GetIndex());
  // end
  atom = mol->GetAtomById(config.end);
  if (atom)
    config.end = symClasses.at(atom->GetIndex());
  // refs
  for (unsigned int i = 0; i < config.refs.size(); ++i) {
    atom = mol->GetAtomById(config.refs.at(i));
    if (atom)
      config.refs[i] = symClasses.at(atom->GetIndex());
  }
}

void IdsToSymClasses(OBMol *mol, OBTetrahedralStereo::Config &config, 
    const std::vector<unsigned int> &symClasses)
{
  OBAtom *atom;
  // center
  atom = mol->GetAtomById(config.center);
  if (atom)
    config.center = symClasses.at(atom->GetIndex());
  // from/towards
  atom = mol->GetAtomById(config.from);
  if (atom)
    config.from = symClasses.at(atom->GetIndex());
  // refs
  for (unsigned int i = 0; i < config.refs.size(); ++i) {
    atom = mol->GetAtomById(config.refs.at(i));
    if (atom)
      config.refs[i] = symClasses.at(atom->GetIndex());
  }
}
void OBGraphSym::BreakChiralTies(OBMol *pmol,
                            OBBitVec &frag_atoms, int nfragatoms,
                            vector<pair<OBAtom*, unsigned int> > &atom_sym_classes)
{
  vector<pair<OBAtom*,unsigned int> > vp1, vp2;

  // for keeping track of atoms we've already considered
  OBBitVec used_atoms;
  used_atoms.Clear();
  used_atoms.Resize(pmol->NumAtoms());

  // Convert the atom/class pairs to an array indexed by atom idx.
  // This is just for convenience in the next step.  Note that there
  // will be "holes" in this vector since it's a molecule fragment.
  vector<unsigned int> index2sym_class(pmol->NumAtoms());
  vector<pair<OBAtom*,unsigned int> >::iterator api;
  for (api = atom_sym_classes.begin(); api < atom_sym_classes.end(); api++)
    index2sym_class[api->first->GetIndex()] = api->second;

  // 
  // Tetrahedral atoms
  //
  vector<unsigned long> tetrahedral = FindTetrahedralAtoms(_pmol, index2sym_class);
  for (vector<unsigned long>::iterator id1 = tetrahedral.begin(); id1 != tetrahedral.end(); ++id1) {
    OBAtom *atom1 = _pmol->GetAtomById(*id1);
    int index1 = atom1->GetIndex();
    // We only want: unused and part of this fragment.
    if (!frag_atoms[index1])
      continue;
    if (used_atoms[index1])
      continue;
    used_atoms.SetBitOn(index1);
    //if (GetValence(atom) < 4) // Valence relative to this fragment
    //  continue;

    // Get its symmetry class
    int symclass = index2sym_class[index1];
 
    // Start the vector of "same class" atoms by adding this atom
    vector<OBAtom *> same_class;
    same_class.push_back(atom1);

    // Inner Loop: Find all other atoms with the same symmetry class
    for (vector<unsigned long>::iterator id2 = id1+1; id2 != tetrahedral.end(); ++id2) {
      OBAtom *atom2 = _pmol->GetAtomById(*id2);
      int index2 = atom2->GetIndex();
      if (used_atoms[index2])
        continue;
      if (index2sym_class[index2] == symclass) {
        same_class.push_back(atom2);
        used_atoms.SetBitOn(index2);
      }
    }
    
    // Unless at least two atoms in the class, there are no ties to break
    if (same_class.size() < 2)
      continue;

    /*
    cout << "BreakChiralTies: same_class = ";
    vector<OBAtom*>::iterator ia;
    for (ia = same_class.begin(); ia != same_class.end(); ia++)
      cout << (*ia)->GetIndex() << " ";
    cout << "\n";
    */

    // find tetrahedral atoms using current symmetry classes
    _pmol->DeleteData(OBGenericDataType::StereoData);
    if (_pmol->Has3D())
      TetrahedralFrom3D(_pmol, index2sym_class);
    else
      TetrahedralFrom2D(_pmol, index2sym_class);

    // false = don't perceive stereochemistry, we already did it explicitly above
    OBStereoFacade stereoFacade(_pmol, false); 
    // get the reference config
    OBTetrahedralStereo::Config refConfig = stereoFacade.GetTetrahedralStereo(*id1)->GetConfig();
    // this 'odd' variable will always be the same (see below...)
    bool odd = OBStereo::NumInversions(refConfig.refs) % 2; // FIXME
    // convert id --> symmetry classes to compare
    IdsToSymClasses(_pmol, refConfig, index2sym_class);

    // Devide the cis/trans in two groups based on C/T by comparing the Configs
    vector<OBAtom*> symclass1, symclass2;

    vector<OBAtom*>::iterator iatom;
    for (iatom = same_class.begin(); iatom != same_class.end(); iatom++) {
      OBTetrahedralStereo::Config otherConfig = stereoFacade.GetTetrahedralStereo((*iatom)->GetId())->GetConfig();
      IdsToSymClasses(_pmol, otherConfig, index2sym_class);

      // compare
      if (refConfig == otherConfig)
        symclass1.push_back(*iatom); // same, add it to symclass1
      else
        symclass2.push_back(*iatom); // different, add to symclass2
    }

    // If there's nothing in symclass2, then we don't have to split 
    // the symmetry class.
    if (symclass2.empty())
      continue;
 
    // Time to break the class in two. Double all symmetry classes
    // Then, either add one or subtract one to all the atoms in symclass2,
    // the atoms that didn't match the chirality of the "reference" class above.
    // The add-or-subtract decision is based on odd, the "odd" calculated
    // above; by doing this, we force a consistent choice of '@' or '@@' for the,
    // two chiral centers; otherwise it's arbitrary and you'll get two different SMILES.
    for (int i = 0; i < atom_sym_classes.size(); i++) {
      atom_sym_classes[i].second *= 2;
      for (int j = 0; j < symclass2.size(); j++) {
        if (symclass2[j] == atom_sym_classes[i].first) {
          if (odd)
            atom_sym_classes[i].second += 1;
          else
            atom_sym_classes[i].second -= 1;
        }
      }
    }
    // Now propagate the change across the whole molecule with the
    // extended sum-of-invariants.
    ExtendInvariants(atom_sym_classes, nfragatoms, pmol->NumAtoms());
  }

  // 
  // Cis/Trans bonds
  //
  vector<unsigned long> cistrans = FindCisTransBonds(_pmol, index2sym_class);
  for (vector<unsigned long>::iterator id1 = cistrans.begin(); id1 != cistrans.end(); ++id1) {
    OBBond *bond1 = _pmol->GetBondById(*id1);
    unsigned int begin1 = bond1->GetBeginAtom()->GetIndex();
    unsigned int end1 = bond1->GetEndAtom()->GetIndex();
    
    // We only want: unused and part of this fragment.
    if (!frag_atoms[begin1] || !frag_atoms[end1])
      continue;
    if (used_atoms[begin1] || used_atoms[end1])
      continue;
    used_atoms.SetBitOn(begin1);
    used_atoms.SetBitOn(end1);

    // Get its symmetry class
    int beginSymClass = index2sym_class[begin1];
    int endSymClass = index2sym_class[end1];
 
    // Start the vector of "same class" bonds by adding this bond
    vector<OBBond*> same_class;
    same_class.push_back(bond1);

    // Inner Loop: Find all other bonds with the same symmetry classes
    for (vector<unsigned long>::iterator id2 = id1+1; id2 != cistrans.end(); ++id2) {
      OBBond *bond2 = _pmol->GetBondById(*id2);
      unsigned int begin2 = bond2->GetBeginAtom()->GetIndex();
      unsigned int end2 = bond2->GetEndAtom()->GetIndex();
      if (used_atoms[begin2] || used_atoms[end2])
        continue;
      if (((index2sym_class[begin2] == beginSymClass) && (index2sym_class[end2] == endSymClass)) ||
          ((index2sym_class[begin2] == endSymClass) && (index2sym_class[end2] == beginSymClass))) {
        same_class.push_back(bond2);
        used_atoms.SetBitOn(begin2);
        used_atoms.SetBitOn(end2);
      }
    }
    
    // Unless at least two atoms in the class, there are no ties to break
    if (same_class.size() < 2)
      continue;

    /*
    cout << "BreakChiralTies: same_class = ";
    vector<OBBond*>::iterator ib;
    for (ib = same_class.begin(); ib != same_class.end(); ib++)
      cout << (*ib)->GetIdx() << " ";
    cout << "\n";
    */

    // find cis/trans bonds using current symmetry classes
    _pmol->DeleteData(OBGenericDataType::StereoData);
    if (_pmol->Has3D())
      CisTransFrom3D(_pmol, index2sym_class);
    else
      CisTransFrom2D(_pmol, index2sym_class);

    // false = don't perceive stereochemistry, we already did it explicitly above
    OBStereoFacade stereoFacade(_pmol, false); 
    // get the reference config
    OBCisTransStereo::Config refConfig = stereoFacade.GetCisTransStereo(*id1)->GetConfig();
    // this 'odd' variable will always be the same (see below...)
    bool odd = OBStereo::NumInversions(refConfig.refs) % 2; // FIXME
    // convert id --> symmetry classes to compare
    IdsToSymClasses(_pmol, refConfig, index2sym_class);

    // Devide the cis/trans in two groups based on C/T by comparing the Configs
    vector<OBBond*> symclass1, symclass2;

    vector<OBBond*>::iterator ibond;
    for (ibond = same_class.begin(); ibond != same_class.end(); ibond++) {
      OBCisTransStereo::Config otherConfig = stereoFacade.GetCisTransStereo((*ibond)->GetId())->GetConfig();
      IdsToSymClasses(_pmol, otherConfig, index2sym_class);

      // compare
      if (refConfig == otherConfig)
        symclass1.push_back(*ibond); // same, add it to symclass1
      else
        symclass2.push_back(*ibond); // different, add to symclass2
    }

    // If there's nothing in symclass2, then we don't have to split 
    // the symmetry class.
    if (symclass2.empty())
      continue;
    
    // Time to break the class in two. Double all symmetry classes
    // Then, either add one or subtract one to all the atoms in symclass2,
    // the atoms that didn't match the chirality of the "reference" class above.
    // The add-or-subtract decision is based on odd, the number of inversions of 
    // the U shaped reference config above; by doing this, we force a consistent 
    // choice of '\' or '/' for the, two double bonds; otherwise it's arbitrary 
    // and you'll get two different SMILES.
    for (int i = 0; i < atom_sym_classes.size(); i++) {
      atom_sym_classes[i].second *= 2;
      for (int j = 0; j < symclass2.size(); j++) {
        // changing the begin atom should be enough
        if (symclass2[j]->GetBeginAtom() == atom_sym_classes[i].first) {
          if (odd)
            atom_sym_classes[i].second += 1;
          else
            atom_sym_classes[i].second -= 1;
        }
      }
    }

    // Now propagate the change across the whole molecule with the
    // extended sum-of-invariants.
    ExtendInvariants(atom_sym_classes, nfragatoms, _pmol->NumAtoms());
  }
 
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
      BreakChiralTies(_pmol, *_frag_atoms, nfragatoms, symmetry_classes);
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
