/**********************************************************************
obiter.h - STL-style iterators for Open Babel.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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

#ifndef OB_OBITER_H
#define OB_OBITER_H

#include <openbabel/babelconfig.h>
#include <openbabel/bitvec.h>

#include <vector>
#include <stack>
#include <queue>


namespace OpenBabel
{
  class OBMol;
  class OBAtom;
  class OBBond;
  class OBResidue;

  // more detailed descriptions and documentation in obiter.cpp

  //! \brief Iterate over all atoms in an OBMol
  class OBAPI OBMolAtomIter {
    std::vector<OBAtom*>::iterator _i;
    OBMol *_parent;
    OBAtom *_ptr;
  public:

    OBMolAtomIter() :_parent(nullptr), _ptr(nullptr) { }
    OBMolAtomIter(OBMol *mol);
    OBMolAtomIter(OBMol &mol);
    OBMolAtomIter(const OBMolAtomIter &ai);
    ~OBMolAtomIter() { }

    OBMolAtomIter& operator=(const OBMolAtomIter &ai);
    //! \return Whether the iterator can still advance (i.e., visit more atoms)
    operator bool() const        { return _ptr != nullptr; }
    //! Preincrement iterator -- advance to next atom and return
    OBMolAtomIter& operator++();
    //! Postincrement iterator -- return the current atom and advance
    OBMolAtomIter  operator++(int);
    //! \return a pointer to the current atom
    OBAtom* operator->() const   { return _ptr;      }
    //! \return a reference to the current atom
    OBAtom& operator*() const    { return *_ptr;     }
  };

  //! \brief Iterate over all atoms in an OBMol in a depth-first search (DFS)
  class OBAPI OBMolAtomDFSIter {
    OBMol               *_parent;
    OBAtom              *_ptr;
    OBBitVec             _notVisited;
    std::stack<OBAtom *> _stack;
  public:

    OBMolAtomDFSIter() : _parent(nullptr), _ptr(nullptr) { }
    OBMolAtomDFSIter(OBMol *mol, int StartIndex=1);
    OBMolAtomDFSIter(OBMol &mol, int StartIndex=1);
    OBMolAtomDFSIter(const OBMolAtomDFSIter &ai);
    ~OBMolAtomDFSIter() { }

    OBMolAtomDFSIter& operator=(const OBMolAtomDFSIter &ai);
    //! \return Whether the iterator can still advance (i.e., visit more atoms)
    operator bool() const        { return _ptr != nullptr; }
    //! Preincrement -- advance to the next atom in the DFS list and return
    OBMolAtomDFSIter& operator++();
    //! Postincrement -- return the current atom and advance to the next atom
    OBMolAtomDFSIter  operator++(int);
    //! \return a pointer to the current atom
    OBAtom* operator->() const   { return _ptr;      }
    //! \return a reference to the current atom
    OBAtom& operator*() const    { return *_ptr;     }
    /// \return NULL if at the last atom in a fragment, else the next atom
    OBAtom* next()
    {
      if(_stack.empty())
        return nullptr; //end of a disconnected fragment
      else
        return _stack.top(); //the next atom
    }
  };

  //! \brief Iterate over all atoms in an OBMol in a breadth-first search (BFS)
  class OBAPI OBMolAtomBFSIter {
    OBMol               *_parent;
    OBAtom              *_ptr;
    OBBitVec             _notVisited;
    std::queue<OBAtom *> _queue;
    std::vector<int>     _depth;
  public:

    OBMolAtomBFSIter(): _parent(nullptr), _ptr(nullptr) { }
    OBMolAtomBFSIter(OBMol *mol, int StartIndex = 1);
    OBMolAtomBFSIter(OBMol &mol, int StartIndex = 1);
    OBMolAtomBFSIter(const OBMolAtomBFSIter &ai);
    ~OBMolAtomBFSIter() { }

    OBMolAtomBFSIter& operator=(const OBMolAtomBFSIter &ai);
    //! \return Whether the iterator can still advance (i.e., visit more atoms)
    operator bool() const        { return _ptr != nullptr; }
    //! Preincrement -- advance to the next atom in the BFS list and return
    OBMolAtomBFSIter& operator++();
    //! Postincrement -- return the current atom and advance to the next atom
    OBMolAtomBFSIter  operator++(int);
    //! \return a pointer to the current atom
    OBAtom* operator->() const   { return _ptr;      }
    //! \return a reference to the current atom
    OBAtom& operator*() const    { return *_ptr;     }
    //! \return the current depth of the iterator
    //! \since version 2.2
    int CurrentDepth() const;
  };

  //! \brief Iterate over all bonds in an OBMol in a breadth-first search (BFS)
  class OBAPI OBMolBondBFSIter {
    OBMol               *_parent;
    OBBond              *_ptr;
    OBBitVec             _notVisited;
    std::queue<OBBond *> _queue;
    std::vector<int>     _depth;
  public:

    OBMolBondBFSIter(): _parent(nullptr), _ptr(nullptr) { }
    OBMolBondBFSIter(OBMol *mol, int StartIndex = 0);
    OBMolBondBFSIter(OBMol &mol, int StartIndex = 0);
    OBMolBondBFSIter(const OBMolBondBFSIter &ai);
    ~OBMolBondBFSIter() { }

    OBMolBondBFSIter& operator=(const OBMolBondBFSIter &ai);
    //! \return Whether the iterator can still advance (i.e., visit more atoms)
    operator bool() const        { return _ptr != nullptr; }
    //! Preincrement -- advance to the next atom in the BFS list and return
    OBMolBondBFSIter& operator++();
    //! Postincrement -- return the current atom and advance to the next atom
    OBMolBondBFSIter  operator++(int);
    //! \return a pointer to the current atom
    OBBond* operator->() const   { return _ptr;      }
    //! \return a reference to the current atom
    OBBond& operator*() const    { return *_ptr;     }
    //! \return the current depth of the iterator
    //! \since version 2.2
    int CurrentDepth() const;
  };

  //! \brief Iterate over all bonds in an OBMol
  class OBAPI OBMolBondIter {
    std::vector<OBBond*>::iterator _i;
    OBMol *_parent;
    OBBond *_ptr;
  public:

    OBMolBondIter() : _parent(nullptr), _ptr(nullptr) {}
    OBMolBondIter(OBMol *mol);
    OBMolBondIter(OBMol &mol);
    OBMolBondIter(const OBMolBondIter &bi);
    ~OBMolBondIter() { }

    OBMolBondIter& operator=(const OBMolBondIter &bi);
    //! \return Whether the iterator can still advance (i.e., visit more bonds)
    operator bool() const        { return _ptr != nullptr; }
    //! Preincrement -- advance to the next bond and return
    OBMolBondIter& operator++();
    //! Postincrement -- return the current bond and advance to the next
    OBMolBondIter  operator++(int);
    //! \return a pointer to the current bond
    OBBond* operator->() const   { return _ptr;      }
    //! \return a reference to the current bond
    OBBond& operator*() const    { return *_ptr;     }
  };

  //! \brief Iterate over all neighboring atoms to an OBAtom
  class OBAPI OBAtomAtomIter {
    std::vector<OBBond*>::iterator _i;
    OBAtom *_parent;
    OBAtom *_ptr;
  public:

    OBAtomAtomIter() : _parent(nullptr), _ptr(nullptr) { }
    OBAtomAtomIter(OBAtom *atm);
    OBAtomAtomIter(OBAtom &atm);
    OBAtomAtomIter(const OBAtomAtomIter &ai);
    ~OBAtomAtomIter() { }

    OBAtomAtomIter& operator=(const OBAtomAtomIter &ai);
     //! \return Whether the iterator can still advance (i.e., visit more neighbors)
    operator bool() const        { return _ptr != nullptr; }
    //! Preincrement -- advance to the next neighbor and return
    OBAtomAtomIter& operator++();
    //! Postincrement -- return the current neighbor and advance to the next
    OBAtomAtomIter  operator++(int);
    //! \return a pointer to the current atom
    OBAtom* operator->() const   { return _ptr;      }
    //! \return a reference to the current atom
    OBAtom& operator*() const    { return *_ptr;     }
  };

  //! \brief Iterate over all bonds on an OBAtom
  class OBAPI OBAtomBondIter {
    std::vector<OBBond*>::iterator _i;
    OBAtom *_parent;
    OBBond *_ptr;
  public:

    OBAtomBondIter(): _parent(nullptr), _ptr(nullptr) { }
    OBAtomBondIter(OBAtom *atm);
    OBAtomBondIter(OBAtom &atm);
    OBAtomBondIter(const OBAtomBondIter &bi);
    ~OBAtomBondIter() { }

    OBAtomBondIter& operator=(const OBAtomBondIter &bi);
    //! \return Whether the iterator can still advance (i.e., visit more bonds)
    operator bool() const        { return _ptr != nullptr; }
    //! Preincrement -- advance to the next bond and return
    OBAtomBondIter& operator++();
    //! Postincrement -- return the current state and advance to the next bond
    OBAtomBondIter  operator++(int);
    //! \return a pointer to the current bond
    OBBond* operator->() const   { return _ptr; }
    //! \return a reference to the current bond
    OBBond& operator*() const    { return *_ptr;}
  };

  //! \brief Iterate over all residues in an OBMol
  class OBAPI OBResidueIter {
    std::vector<OBResidue*>::iterator _i;
    OBResidue *_ptr;
    OBMol *_parent;
  public:

    OBResidueIter() : _ptr(nullptr), _parent(nullptr) { }
    OBResidueIter(OBMol *mol);
    OBResidueIter(OBMol &mol);
    OBResidueIter(const OBResidueIter &ri);
    ~OBResidueIter() { }

    OBResidueIter& operator=(const OBResidueIter &ri);
    //! \return Whether the iterator can still advance (i.e., visit more residues)
    operator bool() const        { return _ptr != nullptr; }
    //! Preincrement -- advance to the next residue and return
    OBResidueIter& operator++();
    //! Postincrement -- return the current state and advance to the next residue
    OBResidueIter  operator++(int);
    //! \return a pointer to the current residue
    OBResidue* operator->() const{ return _ptr; }
    //! \return a reference to the current residue
    OBResidue& operator*() const { return *_ptr;}
  };

  //! \brief Iterate over all atoms in an OBResidue
  class OBAPI OBResidueAtomIter {
    std::vector<OBAtom*>::iterator _i;
    OBResidue *_parent;
    OBAtom    *_ptr;
  public:

    OBResidueAtomIter() : _parent(nullptr), _ptr(nullptr) { }
    OBResidueAtomIter(OBResidue *res);
    OBResidueAtomIter(OBResidue &res);
    OBResidueAtomIter(const OBResidueAtomIter &ri);
    ~OBResidueAtomIter() { }

    OBResidueAtomIter &operator = (const OBResidueAtomIter &ri);
    //! \return Whether the iterator can still advance (i.e., visit more atoms in this residue)
    operator bool() const        { return _ptr != nullptr; }
    //! Preincrement -- advance to the next atom (if any) and return
    OBResidueAtomIter& operator++ ();
    //! Postincrement -- return the current state and advance to the next atom (if any)
    OBResidueAtomIter  operator++ (int);
    //! \return a pointer to the current atom
    OBAtom *operator->() const   { return _ptr; }
    //! \return a reference to the current atom
    OBAtom &operator*() const    { return *_ptr;}
  };

  //! \brief Iterate over all angles in an OBMol
  class OBAPI OBMolAngleIter {
    OBMol     *_parent;
    std::vector<std::vector<unsigned int> > _vangle;
    std::vector<std::vector<unsigned int> >::iterator _i;
    std::vector<unsigned int> _angle;
  public:

    OBMolAngleIter() :_parent(nullptr) { }
    OBMolAngleIter(OBMol *mol);
    OBMolAngleIter(OBMol &mol);
    OBMolAngleIter(const OBMolAngleIter &ai);
    ~OBMolAngleIter() { }

    OBMolAngleIter& operator=(const OBMolAngleIter &ai);
    //! \return Whether the iterator can still advance (i.e., visit more angles)
    operator bool() const        { return (_i != _vangle.end()); }
    //! Preincrement -- advance to the next angle and return
    OBMolAngleIter& operator++();
    //! \return A vector of three atom indexes specifying the angle
    //! \see OBAtom::GetIdx()
    std::vector<unsigned int> operator*() const    { return _angle;     }
  };

  //! \brief Iterate over all torsions in an OBMol
  class OBAPI OBMolTorsionIter {
    OBMol *_parent;
    std::vector<std::vector<unsigned int> > _vtorsion;
    std::vector<std::vector<unsigned int> >::iterator _i;
    std::vector<unsigned int> _torsion;
  public:

    OBMolTorsionIter() :_parent(nullptr) { }
    OBMolTorsionIter(OBMol *mol);
    OBMolTorsionIter(OBMol &mol);
    OBMolTorsionIter(const OBMolTorsionIter &ai);
    ~OBMolTorsionIter() { }

    OBMolTorsionIter& operator=(const OBMolTorsionIter &ai);
    //! \return Whether the iterator can still advance (i.e., visit more torsions)
    operator bool() const        { return (_i != _vtorsion.end()); }
    //! Preincrement -- advance to the next torsion and return
    OBMolTorsionIter& operator++();
    //! \return A vector of four atom indexes specifying the torsion
    //! \see OBAtom::GetIdx()
    std::vector<unsigned int> operator*() const    { return _torsion;     }
  };

  //! \brief Iterate over all pairs of atoms (>1-4) in an OBMol
  class OBAPI OBMolPairIter {
    std::vector<OBAtom*>::iterator _i;
    std::vector<OBAtom*>::iterator _j;
    OBMol *_parent;
    //std::vector<std::vector<unsigned int> > _vpair;
    //std::vector<std::vector<unsigned int> >::iterator _i;
    std::vector<unsigned int> _pair;

  public:

    OBMolPairIter() :_parent(nullptr) { }
    OBMolPairIter(OBMol *mol);
    OBMolPairIter(OBMol &mol);
    OBMolPairIter(const OBMolPairIter &ai);
    ~OBMolPairIter() { }

    OBMolPairIter& operator=(const OBMolPairIter &ai);
    //! \return Whether the iterator can still advance (i.e., visit more 1-4 atom pairs)
    operator bool() const        { return _pair.size()>0; }
    //! Preincrement -- advance to the next 1-4 atom pair and return
    OBMolPairIter& operator++();
    //! \return A vector of two atom indexes specifying the 1-4 atom pair
    //! \see OBAtom::GetIdx()
    std::vector<unsigned int> operator*() const    { return _pair;     }
  };

  class OBRing;
  class OBRingData;

  //! \brief Iterate over all rings in an OBMol
  class OBAPI OBMolRingIter {
    std::vector<OBRing*>::iterator _i;
    OBRing *_ptr;
    OBMol *_parent;
    OBRingData *_rings;
  public:

    OBMolRingIter() : _ptr(nullptr), _parent(nullptr), _rings(nullptr) { }
    OBMolRingIter(OBMol *mol);
    OBMolRingIter(OBMol &mol);
    OBMolRingIter(const OBMolRingIter &ri);
    ~OBMolRingIter() { }

    OBMolRingIter& operator=(const OBMolRingIter &ri);
    //! \return Whether the iterator can advance (i.e., there are more rings)
    operator bool()      const { return _ptr != nullptr; }
    //! Preincrement -- advance to the next ring (if any) and return
    OBMolRingIter& operator++();
    //! Postincrement -- return the current state and advance to the next ring
    OBMolRingIter  operator++(int);
    //! \return A pointer to the current ring (if any)
    OBRing* operator->() const { return _ptr; }
    //! \return A reference to the current ring (if any)
    OBRing& operator*()  const { return *_ptr;}
  };

#define FOR_ATOMS_OF_MOL(a,m)     for( OpenBabel::OBMolAtomIter     a(m); a; ++a )
#define FOR_BONDS_OF_MOL(b,m)     for( OpenBabel::OBMolBondIter     b(m); b; ++b )
#define FOR_NBORS_OF_ATOM(a,p)    for( OpenBabel::OBAtomAtomIter    a(p); a; ++a )
#define FOR_BONDS_OF_ATOM(b,p)    for( OpenBabel::OBAtomBondIter    b(p); b; ++b )
#define FOR_RESIDUES_OF_MOL(r,m)  for( OpenBabel::OBResidueIter     r(m); r; ++r )
#define FOR_ATOMS_OF_RESIDUE(a,r) for( OpenBabel::OBResidueAtomIter a(r); a; ++a )
#define FOR_DFS_OF_MOL(a,m)       for( OpenBabel::OBMolAtomDFSIter  a(m); a; ++a )
#define FOR_BFS_OF_MOL(a,m)       for( OpenBabel::OBMolAtomBFSIter  a(m); a; ++a )
#define FOR_BONDBFS_OF_MOL(b,m)   for( OpenBabel::OBMolBondBFSIter  b(m); b; ++b )
#define FOR_RINGS_OF_MOL(r,m)     for( OpenBabel::OBMolRingIter     r(m); r; ++r )
#define FOR_ANGLES_OF_MOL(a,m)    for( OpenBabel::OBMolAngleIter    a(m); a; ++a )
#define FOR_TORSIONS_OF_MOL(t,m)  for( OpenBabel::OBMolTorsionIter  t(m); t; ++t )
#define FOR_PAIRS_OF_MOL(p,m)     for( OpenBabel::OBMolPairIter     p(m); p; ++p )

} // namespace OpenBabel
#endif // OB_OBITER_H

//! \file obiter.h
//! \brief STL-style iterators for Open Babel.
