/**********************************************************************
Fingerprint.h - Handle fingerprints
 
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
#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include "mol.h"
#include <list>
#include <set>
#include <vector>

#define MAX_FRAGMENT_SIZE 7

namespace OpenBabel
{


///////////////////////////////////////////////////////////////////////////////
//! \brief Atom + Ring closure in the atom list
class AtomRing
{
private:
    OBAtom *_atom;                 //!< atom
    std::vector<int> _ringClosure; //!< position of the closure atom in the fragment
public:
    AtomRing(OBAtom *atom)
    {
        _atom=atom;
    }
    int GetAtomIdx(void)
    {
        return _atom->GetIdx();
    }
    OBAtom* GetAtom(void)
    {
        return _atom;
    }
    void AddRingClosure(int n)
    {
        _ringClosure.push_back(n);
    }
    void printRingClosure(void);
    std::vector<int> GetRingClosure(void)
    {
        return _ringClosure;
    }
};

///////////////////////////////////////////////////////////////////////////////
//! \brief Linear or cyclic fragment
class fragment
{
private:
    unsigned int _hash;         //!< hash number
    std::list<AtomRing> _atoms;     //!< list of the atoms of the fragment
    int _index[MAX_FRAGMENT_SIZE];  //!< indexes of the atoms of the fragment
public:
    fragment(OBAtom *atom);
    void addAtom(AtomRing atom)
    {
        _atoms.push_back(atom);
    }
    void RemoveLastAtom()
    {
        if ( ! _atoms.empty() )
            _atoms.pop_back();
    }
    void SetRingClosure(OBAtom *atom);
    OBAtom * GetLastAtom(void);
    OBAtom * GetBeforeLastAtom(void);
    void printAtoms(void);
    int getIndex(int i);
    void setIndexes(void);
    void sortIndexes(void);
    void printIndexes(void);
    void setHashNumber(void);
    unsigned int getHash (void)
    {
        return _hash;
    }
    void printHash (void)
    {
        std::cout << (unsigned int)_hash; //This only prints the last 32 bits!!!
    }

};

///////////////////////////////////////////////////////////////////////////////
//! \brief Set of linear or cyclic fragment hashed into a bitstring
class fingerprint
{
private:
    std::string _name;            //!< name of the fingerprint
    std::vector<fragment> _flist; //!< list of all the fragments
    std::vector<fragment> _uflist;//!< list of unique fragments
    fragment **_pf;               //!< pointers to each fragment used for sorting
    int _pfsize;                  //!< size of array _pf
    OBBitVec _fpt;                //!< fingerprint bitstring
    void setFragmentIdx(void);
    void sortFragmentIdx(void);
    int setPointers(void);
    void sortPointers(void);
public:
    fingerprint(void)
    {
        _pf = NULL;
        _pfsize = 0;
        _fpt.Resize(1024);
    }
    fingerprint(const char *name);
    ~fingerprint()
    {
        if (_pf)
        {
            free(_pf);
            _pf=NULL;
        }
    }//if (_pf) delete [] _pf;
    OBBitVec &GetFpt(void)
    {
        return _fpt;
    }
    std::string &GetName(void)
    {
        return _name;
    }
    void GetFragments(OBMol &mol);
    void expandFragment(fragment &f, OBBitVec &avisit, int n);
    void setUniqueFragments(void);
    void HashFragments(void);
    void printFragments(void);
    void printPointers(void);
    void printUniqueFragments(void);
    void printFingerprint(std::ostream &ofs);
    void HashMol(OBMol &mol);
};

int compareInt (const void *a, const void *b);
int compareFragment (const void *a, const void *b);

//**********************************************************
class fingerprint2 //: public BaseFingerprint
{
public:
	//Calculates the fingerprint. OBBitVec needs to be 1024 bits 
	virtual void GetFingerprint(OBMol& mol, OBBitVec& fpt);
private:
	typedef std::set<std::vector<int> > Fset;

	void getFragments(Fset& fragset, std::vector<int>& levels, std::vector<int> curfrag, 
			int level, OBAtom* patom, OBBond* pbond);
	void RemoveDuplicates(Fset& fragset, std::vector<int>& frag);
	int CalcHash(std::vector<int>& frag);
	void PrintFpt(std::vector<int>& f, int hash=0);
};


//*********************************************************
///Global function for getting a fingerprint of specified type in vecwords.
///Puts fingerprint of specified type folded down to nwords of 32bits into vecwords 
bool GetFingerprint(OBMol& mol, std::vector<unsigned int>& vecwords, int nwords=0, int type=0);

/// Calculates Tanimoto coefficient (bits set after AND / bits set after OR) 
double Tanimoto(const std::vector<unsigned int> vec1, const std::vector<unsigned int> vec2); 

} // end namespace OpenBabel

#endif
