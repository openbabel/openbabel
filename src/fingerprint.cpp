/**********************************************************************
Fingerprint.cpp: Handle fingerprints.

Copyright (C) 2004 by Fabien Fontaine.

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

#include "mol.h"
#include "fingerprint.h"

using namespace std;

namespace OpenBabel {

bool WriteFingerprint(ostream &ofs, OBMol &mol)
{
  fingerprint fpt(mol.GetTitle());
  fpt.HashMol(mol);
  fpt.printFingerprint(ofs);

  return(true);
}



void AtomRing::printRingClosure(void)
{
  vector<int>::iterator i;

  for(i=_ringClosure.begin(); i != _ringClosure.end(); i++) {
    cout << ":" << *i ;
  }
}

fragment::fragment(OBAtom *atom) 
{
  AtomRing ar(atom);
  _atoms.push_back(ar);
  _hash=0;
  int i;
  for(i=0; i<MAX_FRAGMENT_SIZE; i++) 
    _index[i]=0;
}

OBAtom *fragment::GetLastAtom(void) 
{
  if (_atoms.empty()) 
    return NULL;
  else
    return (_atoms.back()).GetAtom();
}

OBAtom* fragment::GetBeforeLastAtom(void)
{
  if (_atoms.size() < 2)
    return NULL;
  
  list<AtomRing>::iterator i = _atoms.end();
  i--;i--;
  return (*i).GetAtom();
}


void fragment::SetRingClosure(OBAtom *atom)
{
  list<AtomRing>::iterator i;
  int n, Idx = atom->GetIdx(), pos=0;
  for(i = _atoms.begin(), n=1; i != _atoms.end(); i++, n++) {
    if ( (*i).GetAtomIdx() == Idx)
      pos = n;      
  }
  (_atoms.back()).AddRingClosure(pos);
  return;
}

void fragment::printAtoms(void)
{
  list<AtomRing>::iterator i;
  for(i= _atoms.begin(); i != _atoms.end(); i++) {
    cout << (*i).GetAtomIdx();
    (*i).printRingClosure();
    cout << " ";
  }
}
int fragment::getIndex(int i)
{
  if (i<MAX_FRAGMENT_SIZE) 
    return _index[i];
  else {
    cerr << "Error: Index higher than fragment size\n";
    return 0;
  }
}

void fragment::setIndexes(void)
{
  list<AtomRing>::iterator i;
  int Idx, c;
  for(i= _atoms.begin(), c=0; i != _atoms.end(); i++, c++) {
    Idx = (*i).GetAtomIdx();
    if (c == MAX_FRAGMENT_SIZE) {
      cerr << "Error: Fragment bigger than " << MAX_FRAGMENT_SIZE << endl;
      return;
    }
    _index[c]=Idx;
  }

  return;
}

void fragment::sortIndexes(void)
{
  int i;
  for( i=0; i < MAX_FRAGMENT_SIZE; i++) {
    if (_index[i] == 0)
      break;
  }
  
  qsort (_index, i, sizeof (int), *compareInt);

  return;
}

int compareInt (const void *a, const void *b)
{
 int *v1 = (int *) a, *v2 = (int *) b;
 return (*v1 > *v2) - (*v1 < *v2);
}

void fragment::printIndexes(void)
{
  int i;
  
  for( i=0; i < MAX_FRAGMENT_SIZE; i++) {
    if (_index[i] == 0)
      break;
    cout << _index[i] << " ";
  }

  return;
}

void fragment::setHashNumber(void)
{
  long int anumber=0, bnumber=0, cnumber=0;
  list<AtomRing>::iterator i;
  list<AtomRing>::iterator j;
  OBAtom *atomi, *atomj;
  OBBond *bond;
  //set atom number
  for(i=_atoms.begin(); i != _atoms.end(); i++) {
    atomi = (*i).GetAtom();
    anumber *= 8;
    anumber += atomi->GetAtomicNum();
    anumber *= 2;
    if ( atomi->IsAromatic())
      anumber += 1;
  }
  //set bond number
  OBMol *mol;
  i=j=_atoms.begin();
  j++;
  if (!_atoms.empty()) {
    atomi = (*i).GetAtom();
    mol = (OBMol*) atomi->GetParent();
    for( ; j != _atoms.end(); i++, j++) {
      atomi = (*i).GetAtom();
      atomj = (*j).GetAtom();
      bond = mol->GetBond(atomi,atomj);
      if (!bond) { 
	cout << "Error: No bond between two atoms\n";
	return;
      }
       bnumber *= 4;
       if (bond->IsDouble())
	bnumber +=1;
      if (bond->IsTriple())
	bnumber +=2;
      if (bond->IsAromatic())
	bnumber +=3;     
    }
  }
  //set ring closure number
  vector<int>::iterator v;
  vector<int> rc;
  int ni, nj;
  for(i=_atoms.begin(), ni=1; i != _atoms.end(); i++, ni++) {
    rc = (*i).GetRingClosure();
    atomi = (*i).GetAtom();
    if ( rc.empty() )
      continue;
    for (v = rc.begin(); v != rc.end(); v++) {

      // get the ring closure atom
      for(j=_atoms.begin(), nj=1; j != _atoms.end(); j++, nj++) {
	if (nj== *v) {
	  atomj = (*j).GetAtom();
	  break;
	}
      }
      // get the ring closure bond
      mol = (OBMol*) atomi->GetParent();
      bond = mol->GetBond(atomi,atomj);
      if (!bond) { 
	cout << "Error: No bond between two atoms\n";
	return;
      }
      cnumber *= 8;
      cnumber += ni;
      cnumber *= 8;
      cnumber += nj;
      cnumber *= 4;
      if (bond->IsDouble())
	cnumber +=1;
      if (bond->IsTriple())
	cnumber +=2;
      if (bond->IsAromatic())
	cnumber +=3;
    }
  }
  _hash = anumber + bnumber + cnumber;
  return;
}

fingerprint::fingerprint(const char *name)
{ 
  _name=name; 
  _pf = NULL; 
  _pfsize = 0; 
  _fpt.Resize(1024);
}
void fingerprint::GetFragments(OBMol &mol)
{
  
  OBAtom *atom; 
  OBBitVec avisit;
  avisit.Resize(mol.NumAtoms()+1);

  vector<OBNodeBase*>::iterator i;
  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i)) {
    fragment f(atom);
    expandFragment(f, avisit,1);
  }
  return;
}

void fingerprint::expandFragment(fragment &f, OBBitVec &avisit, int n)
{

  OBAtom *atom = f.GetLastAtom();
  OBAtom *blatom = f.GetBeforeLastAtom();
  OBAtom *nbr;
  OBBond *bond;
  std::vector<OBEdgeBase*>::iterator i;
  int natom,Idx;
 
  if ( !atom ) {
    cerr << "Error: empty fragment\n";
    return;
  }
  Idx = atom->GetIdx();
  avisit.SetBitOn(Idx);

  // check for ring closure
  if (blatom)
    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i)) {
      natom = nbr->GetIdx();
      if (avisit[natom] ) {
	if (blatom->GetIdx() != natom) // ring closure
	  f.SetRingClosure(nbr);
      }
    }
  
  // save the fragment
  _flist.push_back(f);
  if (n < MAX_FRAGMENT_SIZE) {
    //next fragment
    n++;
    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
      {
	natom = nbr->GetIdx();
	if (!avisit[natom] ) {
	  f.addAtom(nbr);
	  expandFragment(f, avisit, n);
	}
	
      }
  }
  f.RemoveLastAtom();
  avisit.SetBitOff(Idx);

  return;
}

void fingerprint::printFragments(void)
{

  int i;
  for(i=0; i< _flist.size(); i++) {
    cout << "Fragment " << i << ": ";
    (_flist[i]).printAtoms();
    cout << ", ";
    (_flist[i]).printIndexes();
    cout << " hash:";
    (_flist[i]).printHash();
    cout << endl;
  }
  return;
}

void fingerprint::sortFragmentIdx(void)
{
  int i;
  for(i=0; i< _flist.size(); i++)
    (_flist[i]).sortIndexes();
  
  return;
}

void fingerprint::setFragmentIdx(void)
{
  int i;
  for(i=0; i< _flist.size(); i++)
    (_flist[i]).setIndexes();
  
  return;
}

int fingerprint::setPointers(void)
{
  if ( _pf ) {
    free (_pf);
    _pf = NULL;
  }
    //delete [] _pf;
    
  int size = _flist.size();
  vector<fragment>::iterator i;
  fragment *frag, **p;

  //_pf = new fragment*[size];
  _pf = (fragment **) calloc(size,sizeof(fragment *));
  if (!_pf) {
    cerr << "Error unable to alloc memory for pointers\n";
    return 0;
  }
    
  _pfsize = size;
  for(i= _flist.begin(),p=_pf; i != _flist.end(); i++, p++) {
    frag = &(*i);
    *p = frag;
  }

  return 1;
}


void fingerprint::printPointers(void)
{
  int i;
  fragment *p;
  for(i=0; i< _pfsize; i++) {
    p = _pf[i];
    p->printAtoms();
    cout << ", ";
    p->printIndexes();
    cout << " hash:";
    p->printHash();
    cout << endl;
  }
  return;
}
void fingerprint::sortPointers(void)
{
  qsort (_pf, _pfsize, sizeof (fragment *), *compareFragment);
  return;
}


int compareFragment (const void *a, const void *b)
{
 fragment **v1 = (fragment **) a, **v2 = (fragment **) b;
 fragment *f1 = *v1, *f2 = *v2;
 int i;
 for(i=0; i<MAX_FRAGMENT_SIZE; i++)
   {
     if ( f1->getIndex(i) < f2->getIndex(i))
       return -1;
     if ( f1->getIndex(i) > f2->getIndex(i))
       return 1;
   }
 
 return 0;
}

void fingerprint::setUniqueFragments(void)
{
  fragment **puf, **pf;
  int isdifferent, i;

  // Store the fragment atom numbers
  setFragmentIdx();

  // Sort these atom numbers for comparison
  sortFragmentIdx();

  // Store the pointers to each fragment
  setPointers();
  
  // Sort these pointer
  sortPointers();
  
  // compute the hash number
  for(i=0; i< _flist.size(); i++) {
    (_flist[i]).setHashNumber();
  }
  // printPointers();
  // Get the fragments with the lower hash number
  if (_pfsize == 0)
    return;
  if ( _pfsize == 1) {
    _uflist.push_back(**_pf);
    return;
  }
  
  puf = _pf;
  for (pf = _pf+1, i=1; i < _pfsize; pf++, i++) {
  
    isdifferent = compareFragment(puf, pf);
    if (isdifferent) {
      _uflist.push_back(**puf);
      puf=pf;
    }
    else {
      if ( (*pf)->getHash() < (*puf)->getHash()  )
	puf = pf;
    }
  }
  _uflist.push_back(**puf);

  //printUniqueFragments();
  return;
}

void fingerprint::printUniqueFragments(void)
{

  int i;
  for(i=0; i< _uflist.size(); i++) {
    cout << "Fragment " << i << ": ";
    (_uflist[i]).printAtoms();
    cout << ", ";
    (_uflist[i]).printIndexes();
    cout << " hash:";
    (_uflist[i]).printHash();
    cout << endl;
  }
  return;
}

void fingerprint::HashFragments(void)
{
  vector<fragment>::iterator i;
  long int h;
  int mod, test;
  for(i= _uflist.begin(); i != _uflist.end(); i++) {
    h = i->getHash();
    mod = h % 1021;
    _fpt.SetBitOn(mod);
  }
}

void fingerprint::printFingerprint(ostream &ofs)
{
  int i;
  ofs << ">" << _name << endl;
  for(i=0; i< _fpt.GetSize()*SETWORD; i++ ) {
    if (_fpt[i])
      ofs << "1";
    else 
      ofs << "0";
  }
  ofs << endl;
}

void fingerprint::HashMol(OBMol &mol)
{
  // Obtain all the linear and ring fragments
  GetFragments(mol);
  // Remove duplicate fragments and get a hash number for each fragment
  setUniqueFragments();
  // Obtain the fingerprint of the fragments
  HashFragments();

  return;
}



}
