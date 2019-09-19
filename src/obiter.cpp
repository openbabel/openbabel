/**********************************************************************
obiter.cpp - STL-style iterators for Open Babel

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
#include <openbabel/babelconfig.h>

#include <vector>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/generic.h>
#include <openbabel/obiter.h>

using namespace std;

namespace OpenBabel
{

  /** \class OBMolAtomIter obiter.h <openbabel/obiter.h>

      To facilitate iteration through all atoms in a molecule, without resorting
      to atom indexes (which <strong>will</strong> change in the future), a
      variety of iterator methods are provided.

      This has been made significantly easier by a series of macros in the
      obiter.h header file:

      \code
      \#define FOR_ATOMS_OF_MOL(a,m)     for( OBMolAtomIter     a(m); a; ++a )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      double exactMass = 0.0;
      FOR_ATOMS_OF_MOL(a, mol)
      {
         // The variable a behaves like OBAtom* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*a

         exactMass +=  a->GetExactMass();
      }
      \endcode
  **/

  OBMolAtomIter::OBMolAtomIter(OBMol *mol)
  {
    _parent = mol;
    _ptr = _parent->BeginAtom(_i);
  }

  OBMolAtomIter::OBMolAtomIter(OBMol &mol)
  {
    _parent = &mol;
    _ptr = _parent->BeginAtom(_i);
  }

  OBMolAtomIter::OBMolAtomIter(const OBMolAtomIter &ai)
  {
    _parent = ai._parent;
    _ptr = ai._ptr;
    _i = ai._i;
  }

  OBMolAtomIter& OBMolAtomIter::operator=(const OBMolAtomIter &ai)
  {
    if (this != &ai)
      {
        _parent = ai._parent;
        _ptr = ai._ptr;
        _i = ai._i;
      }
    return *this;
  }

  OBMolAtomIter& OBMolAtomIter::operator++()
  {
    _ptr = _parent->NextAtom(_i);
    return *this;
  }

  OBMolAtomIter OBMolAtomIter::operator++(int)
  {
    OBMolAtomIter tmp(*this);
    operator++();
    return tmp;
  }

  /** \class OBMolAtomDFSIter obiter.h <openbabel/obiter.h>

      \since version 2.1

      To facilitate iteration through all atoms in a molecule, without resorting
      to atom indexes (which <strong>will</strong> change in the future), a
      variety of iterator methods are provided.

      This class provides a depth-first search ordering of atoms. When one
      connected component is exhausted, the iterator will start at another until
      all atoms are visited. No guarantee is made as to the ordering of
      iteration through connected components.

      The iterator maintains an internal stack and list of visited
      atoms. As such it may not be appropriate for memory-constrained
      situations when iterating through large molecules.

      Use of this iterator has been made significantly easier by a series
      of macros in the obiter.h header file:

      \code
      \#define FOR_DFS_OF_MOL(a,m)     for( OBMolAtomDFSIter     a(m); a; ++a )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      FOR_DFS_OF_MOL(a, mol)
      {
         // The variable a behaves like OBAtom* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*a

      }
      \endcode
  **/

  OBMolAtomDFSIter::OBMolAtomDFSIter(OBMol *mol, int StartIndex):
    _parent(mol), _ptr(_parent->GetAtom(StartIndex))
  {
    if (!_ptr) return;

    _notVisited.Resize(_parent->NumAtoms());
    _notVisited.SetRangeOn(0, _parent->NumAtoms() - 1);
    
    _notVisited.SetBitOff(_ptr->GetIdx() - 1);

    vector<OBBond*>::iterator i;
    OBAtom *a;

    for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
      {
        _stack.push(a);
        _notVisited.SetBitOff(a->GetIdx() - 1);
      }
  }

  OBMolAtomDFSIter::OBMolAtomDFSIter(OBMol &mol, int StartIndex):
    _parent(&mol), _ptr(_parent->GetAtom(StartIndex))
  {
    if (!_ptr) return;

    _notVisited.Resize(_parent->NumAtoms());
    _notVisited.SetRangeOn(0, _parent->NumAtoms() - 1);
    
    _notVisited.SetBitOff(_ptr->GetIdx() - 1);

    vector<OBBond*>::iterator i;
    OBAtom *a;

    for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
      {
        _stack.push(a);
        _notVisited.SetBitOff(a->GetIdx() - 1);
      }
  }

  OBMolAtomDFSIter::OBMolAtomDFSIter(const OBMolAtomDFSIter &ai)
  {
    _parent = ai._parent;
    _ptr = ai._ptr;
    _notVisited = ai._notVisited;
    _stack = ai._stack;
  }

  OBMolAtomDFSIter& OBMolAtomDFSIter::operator=(const OBMolAtomDFSIter &ai)
  {
    if (this != &ai)
      {
        _parent = ai._parent;
        _ptr = ai._ptr;
        _notVisited = ai._notVisited;
        _stack = ai._stack;
      }
    return *this;
  }

  OBMolAtomDFSIter& OBMolAtomDFSIter::operator++()
  {
    if (!_stack.empty())
      {
        _ptr = _stack.top();
        _stack.pop();
      }
    else // are there any disconnected subgraphs?
      {
        int next = _notVisited.FirstBit();
        if (next != _notVisited.EndBit())
          {
            _ptr = _parent->GetAtom(next + 1);
            _notVisited.SetBitOff(next);
          }
        else
          _ptr = NULL;
      }

    if (_ptr)
      {
        vector<OBBond*>::iterator i;
        OBAtom *a;

        for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
          if (_notVisited[a->GetIdx() - 1])
            {
              _stack.push(a);
              _notVisited.SetBitOff(a->GetIdx() - 1);
            }
      }

    return *this;
  }

  OBMolAtomDFSIter OBMolAtomDFSIter::operator++(int)
  {
    OBMolAtomDFSIter tmp(*this);
    operator++();
    return tmp;
  }

  /** \class OBMolAtomBFSIter obiter.h <openbabel/obiter.h>

      \since version 2.1

      To facilitate iteration through all atoms in a molecule, without resorting
      to atom indexes (which <strong>will</strong> change in the future), a
      variety of iterator methods are provided.

      This class provides a breadth-first search ordering of atoms. When one
      connected component is exhausted, the iterator will start at another until
      all atoms are visited. No guarantee is made as to the ordering of
      iteration through connected components.

      The iterator maintains an internal queue and list of visited
      atoms. As such it may not be appropriate for memory-constrained
      situations when iterating through large molecules.

      Use of this iterator has been made significantly easier by a series
      of macros in the obiter.h header file:

      \code
      \#define FOR_BFS_OF_MOL(a,m)     for( OBMolAtomBFSIter     a(m); a; ++a )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      FOR_BFS_OF_MOL(a, mol)
      {
         // The variable a behaves like OBAtom* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*a

      }
      \endcode
  **/

  OBMolAtomBFSIter::OBMolAtomBFSIter(OBMol *mol, int StartIndex):
    _parent(mol), _ptr(_parent->GetAtom(StartIndex))
  {
    if (!_ptr) return;

    _notVisited.Resize(_parent->NumAtoms());
    _notVisited.SetRangeOn(0, _parent->NumAtoms() - 1);
    
    _notVisited.SetBitOff(_ptr->GetIdx() - 1);

    // Set up storage for the depths
    _depth.resize(_parent->NumAtoms() + 1, 0);
    _depth[_ptr->GetIdx()] = 1;

    vector<OBBond*>::iterator i;
    OBAtom *a;

    for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
      {
        _queue.push(a);
        _depth[a->GetIdx()] = 2;
        _notVisited.SetBitOff(a->GetIdx() - 1);
      }
  }

  OBMolAtomBFSIter::OBMolAtomBFSIter(OBMol &mol, int StartIndex):
    _parent(&mol), _ptr(_parent->GetAtom(StartIndex))
  {
    if (!_ptr) return;

    _notVisited.Resize(_parent->NumAtoms());
    _notVisited.SetRangeOn(0, _parent->NumAtoms() - 1);
    
    _notVisited.SetBitOff(_ptr->GetIdx() - 1);

    // Set up storage for the depths
    _depth.resize(_parent->NumAtoms() + 1, 0);
    _depth[_ptr->GetIdx()] = 1;

    vector<OBBond*>::iterator i;
    OBAtom *a;

    for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
      {
        _queue.push(a);
        _depth[a->GetIdx()] = 2;
        _notVisited.SetBitOff(a->GetIdx() - 1);
      }
  }

  OBMolAtomBFSIter::OBMolAtomBFSIter(const OBMolAtomBFSIter &ai)
  {
    _parent = ai._parent;
    _ptr = ai._ptr;
    _notVisited = ai._notVisited;
    _queue = ai._queue;
    _depth = ai._depth;
  }

  OBMolAtomBFSIter& OBMolAtomBFSIter::operator=(const OBMolAtomBFSIter &ai)
  {
    if (this != &ai)
      {
        _parent = ai._parent;
        _ptr = ai._ptr;
        _notVisited = ai._notVisited;
        _queue = ai._queue;
        _depth = ai._depth;
      }
    return *this;
  }

  OBMolAtomBFSIter& OBMolAtomBFSIter::operator++()
  {
    if (!_queue.empty())
      {
        _ptr = _queue.front();
        _queue.pop();
      }
    else // are there any disconnected subgraphs?
      {
        int next = _notVisited.FirstBit();
        if (next != _notVisited.EndBit())
          {
            _ptr = _parent->GetAtom(next + 1); // Atom index issue
            if (_ptr != NULL)
              _depth[_ptr->GetIdx()] = 1; // new island
            _notVisited.SetBitOff(next);
          }
        else
          _ptr = NULL;
      }

    if (_ptr)
      {
        vector<OBBond*>::iterator i;
        OBAtom *a;

        for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
          if (_notVisited[a->GetIdx() - 1])
            {
              _queue.push(a);
              _depth[a->GetIdx()] = _depth[_ptr->GetIdx()] + 1;
              _notVisited.SetBitOff(a->GetIdx() - 1);
            }
      }

    return *this;
  }

  OBMolAtomBFSIter OBMolAtomBFSIter::operator++(int)
  {
    OBMolAtomBFSIter tmp(*this);
    operator++();
    return tmp;
  }

  int OBMolAtomBFSIter::CurrentDepth() const
  {
    if (_ptr == NULL)
      return 0;

    return _depth[_ptr->GetIdx()];
  }

  /** \class OBMolBondBFSIter obiter.h <openbabel/obiter.h>

      \since version 2.3

      To facilitate iteration through all bonds in a molecule, without resorting
      to bond indexes (which <strong>will</strong> change in the future), a
      variety of iterator methods are provided.

      This class provides a breadth-first search ordering of bonds. When one
      connected component is exhausted, the iterator will start at another until
      all bonds are visited. No guarantee is made as to the ordering of
      iteration through connected components.

      The iterator maintains an internal queue and list of visited
      atoms. As such it may not be appropriate for memory-constrained
      situations when iterating through large molecules.

      Use of this iterator has been made significantly easier by a series
      of macros in the obiter.h header file:

      \code
      \#define FOR_BONDBFS_OF_MOL(a,m)     for( OBMolBondBFSIter     a(m); a; ++a )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      FOR_BONDBFS_OF_MOL(b, mol)
      {
         // The variable b behaves like OBBond* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*a

      }
      \endcode
  **/

  OBMolBondBFSIter::OBMolBondBFSIter(OBMol *mol, int StartIndex):
    _parent(mol)
  {
    unsigned int numbonds = _parent->NumBonds();
    if (numbonds == 0) {
      _ptr = 0; // mark as invalid
      return;
    }
    _ptr = _parent->GetBond(StartIndex);
    if (!_ptr)
      return;
    
    _notVisited.Resize(numbonds);
    _notVisited.SetRangeOn(0, numbonds - 1);
    
    _notVisited.SetBitOff(_ptr->GetIdx());

    // Set up storage for the depths
    _depth.resize(_parent->NumBonds(), 0);
    _depth[_ptr->GetIdx()] = 1;

    for( OBAtomBondIter b(_ptr->GetBeginAtom()); b; ++b )
      { // Loop thru all bonds attached to the start of the bond
        if (_notVisited[b->GetIdx()]) { // Don't revisit the initial bond
          _queue.push(&(*b));
          _depth[b->GetIdx()] = 2;
          _notVisited.SetBitOff(b->GetIdx());
        }
      }
    for( OBAtomBondIter b(_ptr->GetEndAtom()); b; ++b )
      { // Loop thru all bonds attached to the end of the bond
        if (_notVisited[b->GetIdx()]) { // Don't revisit the initial bond
          _queue.push(&(*b));
          _depth[b->GetIdx()] = 2;
          _notVisited.SetBitOff(b->GetIdx());
        }
      }
  }

  OBMolBondBFSIter::OBMolBondBFSIter(OBMol &mol, int StartIndex):
    _parent(&mol)
  {
    unsigned int numbonds = _parent->NumBonds();
    if (numbonds == 0) {
      _ptr = 0; // mark as invalid
      return;
    }
    _ptr = _parent->GetBond(StartIndex);
    if (!_ptr)
      return;

    _notVisited.Resize(numbonds);
    _notVisited.SetRangeOn(0, numbonds - 1);
    
    _notVisited.SetBitOff(_ptr->GetIdx());

    // Set up storage for the depths
    _depth.resize(numbonds, 0);
    _depth[_ptr->GetIdx()] = 1;

    for( OBAtomBondIter b(_ptr->GetBeginAtom()); b; ++b )
      { // Loop thru all bonds attached to the start of the bond
        if (_notVisited[b->GetIdx()]) { // Don't revisit the initial bond
          _queue.push(&(*b));
          _depth[b->GetIdx()] = 2;
          _notVisited.SetBitOff(b->GetIdx());
        }
      }
    for( OBAtomBondIter b(_ptr->GetEndAtom()); b; ++b )
      { // Loop thru all bonds attached to the end of the bond
        if (_notVisited[b->GetIdx()]) { // Don't revisit the initial bond
          _queue.push(&(*b));
          _depth[b->GetIdx()] = 2;
          _notVisited.SetBitOff(b->GetIdx());
        }
      }
  }

  OBMolBondBFSIter::OBMolBondBFSIter(const OBMolBondBFSIter &ai)
  {
    _parent = ai._parent;
    _ptr = ai._ptr;
    _notVisited = ai._notVisited;
    _queue = ai._queue;
    _depth = ai._depth;
  }

  OBMolBondBFSIter& OBMolBondBFSIter::operator=(const OBMolBondBFSIter &ai)
  {
    if (this != &ai)
      {
        _parent = ai._parent;
        _ptr = ai._ptr;
        _notVisited = ai._notVisited;
        _queue = ai._queue;
        _depth = ai._depth;
      }
    return *this;
  }

  OBMolBondBFSIter& OBMolBondBFSIter::operator++()
  {
    if (!_queue.empty())
    {
      _ptr = _queue.front();
      _queue.pop();
    }
    else // are there any disconnected subgraphs?
    {
      int next = _notVisited.FirstBit();
      if (next != _notVisited.EndBit())
      {
        _ptr = _parent->GetBond(next + 1); // Bond index issue
        if (_ptr != NULL)
          _depth[_ptr->GetIdx()] = 1; // new island
        _notVisited.SetBitOff(next);
      }
      else
        _ptr = NULL;
    }

    if (_ptr) {
      for( OBAtomBondIter b(_ptr->GetBeginAtom()); b; ++b )
        { // Loop thru all bonds attached to the start of the bond
          if (_notVisited[b->GetIdx()]) {
            _queue.push(&(*b));
            _depth[b->GetIdx()] = 2;
            _notVisited.SetBitOff(b->GetIdx());
          }
        }
      for( OBAtomBondIter b(_ptr->GetEndAtom()); b; ++b )
        { // Loop thru all bonds attached to the end of the bond
          if (_notVisited[b->GetIdx()]) {
            _queue.push(&(*b));
            _depth[b->GetIdx()] = 2;
            _notVisited.SetBitOff(b->GetIdx());
          }
        }
      }
    return *this;
  }

  OBMolBondBFSIter OBMolBondBFSIter::operator++(int)
  {
    OBMolBondBFSIter tmp(*this);
    operator++();
    return tmp;
  }

  int OBMolBondBFSIter::CurrentDepth() const
  {
    if (_ptr == NULL)
      return 0;

    return _depth[_ptr->GetIdx()];
  }
///////////////////////////////////////////////////////////////////////

  /** \class OBMolBondIter obiter.h <openbabel/obiter.h>

      To facilitate iteration through all bonds in a molecule, without resorting
      to bond indexes (which may change in the future), a variety of
      iterators are provided.

      This has been made significantly easier by a series of macros in the
      obiter.h header file:

      \code
      \#define FOR_BONDS_OF_MOL(b,m)     for( OBMolBondIter     b(m); b; ++b )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      unsigned int bondOrderSum = 0;
      FOR_BONDS_OF_MOL(b, mol)
      {
         // The variable b behaves like OBBond* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*b
         bondOrderSum +=  b->GetBondOrder();
      }
      \endcode
  **/

  OBMolBondIter::OBMolBondIter(OBMol *mol)
  {
    _parent = mol;
    _ptr = _parent->BeginBond(_i);
  }

  OBMolBondIter::OBMolBondIter(OBMol &mol)
  {
    _parent = &mol;
    _ptr = _parent->BeginBond(_i);
  }

  OBMolBondIter::OBMolBondIter(const OBMolBondIter &bi)
  {
    _parent = bi._parent;
    _ptr = bi._ptr;
    _i = bi._i;
  }

  OBMolBondIter& OBMolBondIter::operator=(const OBMolBondIter &bi)
  {
    if (this != &bi)
      {
        _parent = bi._parent;
        _ptr = bi._ptr;
        _i = bi._i;
      }
    return *this;
  }

  OBMolBondIter& OBMolBondIter::operator++()
  {
    _ptr = _parent->NextBond(_i);
    return *this;
  }

  OBMolBondIter OBMolBondIter::operator++(int)
  {
    OBMolBondIter tmp(*this);
    operator++();
    return tmp;
  }

  /** \class OBAtomAtomIter obiter.h <openbabel/obiter.h>

      To facilitate iteration through all neighbors of an atom, without resorting
      to bond indexes (which may change in the future), a variety of
      iterator classes and methods are provided.

      This has been made significantly easier by a series of macros in the
      obiter.h header file:

      \code
      \#define FOR_NBORS_OF_ATOM(a,p)     for( OBAtomAtomIter     a(p); a; ++a )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      FOR_ATOMS_OF_MOL(a, mol)
      {
         // The variable a behaves like OBAtom* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*a
         FOR_NBORS_OF_ATOM(b, &*a)
         {
            ...
         }
      }
      \endcode
  **/

  OBAtomAtomIter::OBAtomAtomIter(OBAtom *atm)
  {
    _parent = atm;
    _ptr = _parent->BeginNbrAtom(_i);
  }

  OBAtomAtomIter::OBAtomAtomIter(OBAtom &atm)
  {
    _parent = &atm;
    _ptr = _parent->BeginNbrAtom(_i);
  }

  OBAtomAtomIter::OBAtomAtomIter(const OBAtomAtomIter &ai)
  {
    _parent = ai._parent;
    _ptr = ai._ptr;
    _i = ai._i;
  }

  OBAtomAtomIter& OBAtomAtomIter::operator=(const OBAtomAtomIter &ai)
  {
    if (this != &ai)
      {
        _parent = ai._parent;
        _ptr = ai._ptr;
        _i = ai._i;
      }
    return *this;
  }

  OBAtomAtomIter& OBAtomAtomIter::operator++()
  {
    _ptr = _parent->NextNbrAtom(_i);
    return *this;
  }

  OBAtomAtomIter OBAtomAtomIter::operator++(int)
  {
    OBAtomAtomIter tmp(*this);
    operator++();
    return tmp;
  }

  /** \class OBAtomBondIter obiter.h <openbabel/obiter.h>

      To facilitate iteration through all bonds on an atom, without resorting
      to bond indexes (which may change in the future) a variety of
      iterator classes and methods are provided.

      This has been made significantly easier by a series of macros in the
      obiter.h header file:

      \code
      \#define FOR_BONDS_OF_ATOM(b,p)     for( OBAtomBondIter     b(p); b; ++b )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBAtom atom;
      unsigned int tripleBondCount;
      FOR_BONDS_OF_ATOM(b, atom)
      {
         // The variable b behaves like OBBond* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*b
         if (b->GetBondOrder() == 3)
            tripleBondCount++;
      }
      \endcode
  **/

  OBAtomBondIter::OBAtomBondIter(OBAtom *atm)
  {
    _parent = atm;
    _ptr = _parent->BeginBond(_i);
  }

  OBAtomBondIter::OBAtomBondIter(OBAtom &atm)
  {
    _parent = &atm;
    _ptr = _parent->BeginBond(_i);
  }

  OBAtomBondIter::OBAtomBondIter(const OBAtomBondIter &bi)
  {
    _parent = bi._parent;
    _ptr = bi._ptr;
    _i = bi._i;
  }

  OBAtomBondIter& OBAtomBondIter::operator=(const OBAtomBondIter &bi)
  {
    if (this != &bi)
      {
        _parent = bi._parent;
        _ptr = bi._ptr;
        _i = bi._i;
      }
    return *this;
  }

  OBAtomBondIter& OBAtomBondIter::operator++()
  {
    _ptr = _parent->NextBond(_i);
    return *this;
  }

  OBAtomBondIter OBAtomBondIter::operator++(int)
  {
    OBAtomBondIter tmp(*this);
    operator++();
    return tmp;
  }

  /** \class OBResidueIter obiter.h <openbabel/obiter.h>

      To facilitate iteration through all residues in a molecule, without resorting
      to residue indexes (which may change in the future) a variety of
      iterator classes and methods are provided.

      This has been made significantly easier by a series of macros in the
      obiter.h header file:

      \code
      \#define FOR_RESIDUES_OF_MOL(r,m)     for( OBResidueIter     r(m); r; ++r )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      FOR_RESIDUES_OF_MOL(r, mol)
      {
         // The variable r behaves like OBResidue* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*r

         if (r->GetName() == resname && r->GetNum() == rnum)
         {
            // got a match, let's go to work
            ...
         }
      }
      \endcode
  **/

  OBResidueIter::OBResidueIter(OBMol *mol)
  {
    _parent = mol;
    _ptr = _parent->BeginResidue(_i);
  }

  OBResidueIter::OBResidueIter(OBMol &mol)
  {
    _parent = &mol;
    _ptr = _parent->BeginResidue(_i);
  }

  OBResidueIter::OBResidueIter(const OBResidueIter &ri)
  {
    _parent = ri._parent;
    _ptr = ri._ptr;
    _i = ri._i;
  }

  OBResidueIter& OBResidueIter::operator=(const OBResidueIter &ri)
  {
    if (this != &ri)
      {
        _parent = ri._parent;
        _ptr = ri._ptr;
        _i = ri._i;
      }
    return *this;
  }

  OBResidueIter& OBResidueIter::operator++()
  {
    _ptr = _parent->NextResidue(_i);
    return *this;
  }

  OBResidueIter OBResidueIter::operator++(int)
  {
    OBResidueIter tmp(*this);
    operator++();
    return tmp;
  }

  /** \class OBResidueAtomIter obiter.h <openbabel/obiter.h>

      To facilitate iteration through all atoms in a residue, without resorting
      to atom indexes (which may change in the future) a variety of
      iterator classes and methods are provided.

      This has been made significantly easier by a series of macros in the
      obiter.h header file:

      \code
      \#define FOR_ATOMS_OF_RESIDUE(a,r)     for( OBResidueAtomIter     a(r); a; ++a )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      double residueMass = 0.0;
      FOR_RESIDUES_OF_MOL(r, mol)
      {
         // The variable r behaves like OBResidue* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*r

         if (r->GetName() == resname && r->GetNum() == rnum)
         {
            FOR_ATOMS_OF_RESIDUE(a, &*r)
            {
               residueMass += a->GetMass();
            }
         }
      }
      \endcode
  **/

  OBResidueAtomIter::OBResidueAtomIter(OBResidue *res):
    _parent(res), _ptr(_parent->BeginAtom(_i))
  {  }

  OBResidueAtomIter::OBResidueAtomIter(OBResidue &res):
    _parent(&res), _ptr(_parent->BeginAtom(_i))
  {  }

  OBResidueAtomIter::OBResidueAtomIter(const OBResidueAtomIter &ri)
  {
    _parent = ri._parent;
    _ptr    = ri._ptr;
    _i      = ri._i;
  }

  OBResidueAtomIter & OBResidueAtomIter::operator = (const OBResidueAtomIter &ri)
  {
    if (this != &ri)
      {
        _parent = ri._parent;
        _ptr    = ri._ptr;
        _i      = ri._i;
      }

    return (*this);
  }

  OBResidueAtomIter& OBResidueAtomIter::operator++ ()
  {
    _ptr = _parent->NextAtom(_i);
    return (*this);
  }

  OBResidueAtomIter OBResidueAtomIter::operator++ (int)
  {
    OBResidueAtomIter tmp(*this);
    operator++();
    return tmp;
  }

  /** \class OBMolRingIter obiter.h <openbabel/obiter.h>

      \since version 2.1

      To facilitate iteration through all rings in a molecule, without resorting
      to ring indexes (which may change in the future) a variety of
      iterator classes and methods are provided. One word of warning is that
      these iterator methods automatically call OBMol::FindSSSR() which may
      involve a significant performance hit on large molecules.

      Calling iterator classes has been made significantly easier by a series
      of macros in the obiter.h header file:

      \code
      \#define FOR_RINGS_OF_MOL(r,m)     for( OBMolRingIter     r(m); r; ++r )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      FOR_RINGS_OF_MOL(r, mol)
      {
         // The variable r behaves like OBRing* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*r

      }
      \endcode
  **/

  OBMolRingIter::OBMolRingIter(OBMol *mol): _parent(mol)
  {
    if (!_parent->HasSSSRPerceived())
      _parent->FindSSSR();

    _rings = (OBRingData *) _parent->GetData("SSSR");
    if(_rings)
      _ptr = _rings->BeginRing(_i);
  }

  OBMolRingIter::OBMolRingIter(OBMol &mol): _parent(&mol)
  {
    if (!_parent->HasSSSRPerceived())
      _parent->FindSSSR();

    _rings = (OBRingData *) _parent->GetData("SSSR");
    if (_rings)
      _ptr = _rings->BeginRing(_i);
  }

  OBMolRingIter::OBMolRingIter(const OBMolRingIter &ri)
  {
    _parent = ri._parent;
    _ptr = ri._ptr;
    _rings = ri._rings;
    _i = ri._i;
  }

  OBMolRingIter& OBMolRingIter::operator=(const OBMolRingIter &ri)
  {
    if (this != &ri)
      {
        _parent = ri._parent;
        _ptr = ri._ptr;
        _rings = ri._rings;
        _i = ri._i;
      }
    return *this;
  }

  OBMolRingIter& OBMolRingIter::operator++()
  {
    if (_rings)
      _ptr = _rings->NextRing(_i);
    return *this;
  }

  OBMolRingIter OBMolRingIter::operator++(int)
  {
    OBMolRingIter tmp(*this);
    operator++();
    return tmp;
  }

  /** \class OBMolAngleIter obiter.h <openbabel/obiter.h>

      \since version 2.1

      To facilitate iteration through all angles in a molecule, without resorting
      to atom indexes (which <strong>will</strong> change in the future), a
      variety of iterator methods are provided.

      This has been made significantly easier by a series of macros in the
      obiter.h header file:

      \code
      \#define FOR_ANGLES_OF_MOL(a,m)     for( OBMolAngleIter     a(m); a; a++ )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      OpenBabel::OBAtom *a, *b, *c;
      double ang;

      FOR_ANGLES_OF_MOL(angle, mol)
      {
         // The variable a behaves like OBAngle* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*a

         b = _mol.GetAtom((*angle)[0] + 1);
         a = _mol.GetAtom((*angle)[1] + 1);
         c = _mol.GetAtom((*angle)[2] + 1);
         ang = a->GetAngle(b->GetIdx(), c->GetIdx());
      }
      \endcode
  **/

  OBMolAngleIter::OBMolAngleIter(OBMol *mol)
  {
    mol->FindAngles();
    OBAngleData *ad = (OBAngleData *) mol->GetData(OBGenericDataType::AngleData);
    ad->FillAngleArray(_vangle);

    _parent = mol;
    if (!_vangle.empty()) {
      _i = _vangle.begin();
      _angle = *_i;
	} else {
	  _i = _vangle.end();
	}
  }

  OBMolAngleIter::OBMolAngleIter(OBMol &mol)
  {
    mol.FindAngles();
    OBAngleData *ad = (OBAngleData *) mol.GetData(OBGenericDataType::AngleData);
    ad->FillAngleArray(_vangle);

    _parent = &mol;
    if (!_vangle.empty()) {
      _i = _vangle.begin();
      _angle = *_i;
	} else {
      _i = _vangle.end();
	}
  }

  OBMolAngleIter::OBMolAngleIter(const OBMolAngleIter &ai)
  {
    _parent = ai._parent;
    _angle = ai._angle;
    _vangle = ai._vangle;
    _i = ai._i;
  }

  OBMolAngleIter& OBMolAngleIter::operator=(const OBMolAngleIter &ai)
  {
    if (this != &ai)
      {
        _parent = ai._parent;
        _angle = ai._angle;
        _vangle = ai._vangle;
        _i = ai._i;
      }
    return *this;
  }

  OBMolAngleIter& OBMolAngleIter::operator++()
  {
    _i++;

    if (_i != _vangle.end())
      _angle = *_i;

    return *this;
  }

  /** \class OBMolTorsionIter obiter.h <openbabel/obiter.h>

      \since version 2.1

      To facilitate iteration through all torsions in a molecule, without resorting
      to atom indexes (which <strong>will</strong> change in the future), a
      variety of iterator methods are provided.

      This has been made significantly easier by a series of macros in the
      obiter.h header file:

      \code
      \#define FOR_TORSIONS_OF_MOL(t,m)  for( OBMolTorsionIter   t(m); t; t++ )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      OpenBabel::OBAtom *a, *b, *c, *d;
      double tor;

      FOR_TORSIONS_OF_MOL(t, mol)
      {
         // The variable a behaves like OBTorsion* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*t

         a = _mol.GetAtom((*t)[0] + 1); // indices in vector start from 0!!!
         b = _mol.GetAtom((*t)[1] + 1);
         c = _mol.GetAtom((*t)[2] + 1);
         d = _mol.GetAtom((*t)[3] + 1);
         tor = mol.GetTorsion(a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetIdx());
      }
      \endcode
  **/

  OBMolTorsionIter::OBMolTorsionIter(OBMol *mol)
  {
    mol->FindTorsions();
    OBTorsionData *td = (OBTorsionData *) mol->GetData(OBGenericDataType::TorsionData);
    td->FillTorsionArray(_vtorsion);

    _parent = mol;
    if (!_vtorsion.empty()) {
      _i = _vtorsion.begin();
      _torsion = *_i;
	} else {
	  _i = _vtorsion.end();
	}
  }

  OBMolTorsionIter::OBMolTorsionIter(OBMol &mol)
  {
    mol.FindTorsions();
    OBTorsionData *td = (OBTorsionData *) mol.GetData(OBGenericDataType::TorsionData);
    td->FillTorsionArray(_vtorsion);

    _parent = &mol;
    if (!_vtorsion.empty()) {
      _i = _vtorsion.begin();
      _torsion = *_i;
	} else {
	  // Avogadro bug #1972244
	  // always set _i, _i will be compared to _vtorsion.end() in OBMolTorsionIter::bool()
	  _i = _vtorsion.end();
	}
  }

  OBMolTorsionIter::OBMolTorsionIter(const OBMolTorsionIter &ai)
  {
    _parent = ai._parent;
    _torsion = ai._torsion;
    _vtorsion = ai._vtorsion;
    _i = ai._i;
  }

  OBMolTorsionIter& OBMolTorsionIter::operator=(const OBMolTorsionIter &ai)
  {
    if (this != &ai)
      {
        _parent = ai._parent;
        _torsion = ai._torsion;
        _vtorsion = ai._vtorsion;
        _i = ai._i;
      }
      return *this;
  }

  OBMolTorsionIter& OBMolTorsionIter::operator++()
  {
    _i++;

    if (_i != _vtorsion.end())
      _torsion = *_i;

    return *this;
  }

  /** \class OBMolPairIter obiter.h <openbabel/obiter.h>

      \since version 2.1.

      To facilitate iteration through all pairs of atoms in a molecule, without
      resorting to bond indexes (which may change in the future), a variety of
      iterators are provided. These pairs of atoms are separated by 4 atoms
      or more (i.e., these are non-bonded interactions).

      This has been made significantly easier by a series of macros in the
      obiter.h header file:

      \code
      \#define FOR_PAIRS_OF_MOL(p,m)     for( OBMolPairIter     p(m); p; p++ )
      \endcode

      Here is an example:
      \code
      #include <openbabel/obiter.h>
      #include <openbabel/mol.h>

      OpenBabel::OBMol mol;
      OpenBabel::OBAtom *a, *b;
      double rab;

      FOR_PAIRS_OF_MOL(p, mol)
      {
         // The variable b behaves like OBBond* when used with -> and * but
         // but needs to be explicitly converted when appearing as a parameter
         // in a function call - use &*p

         a = mol.GetAtom(p->first);
         b = mol.GetAtom(p->second);
         rab = a->GetDistance(b);
      }
      \endcode
  **/

  OBMolPairIter::OBMolPairIter(OBMol *mol)
  {
    _parent = mol;

    bool foundPair = false;
    OBAtom *a = _parent->BeginAtom(_i);
    if (!a)
      return;
    OBAtom *b = _parent->BeginAtom(_j);
    while (!foundPair) {
      b = _parent->NextAtom(_j);

      if (!b) {
        a = _parent->NextAtom(_i);
        if (!a)
	  return;
        b = _parent->BeginAtom(_j);
      }

      if (a->GetIdx() >= b->GetIdx()) continue;
      if (a->IsConnected(b)) continue;
      if (a->IsOneThree(b)) continue;

      foundPair = true;
    }

    _pair.clear();
    _pair.push_back(a->GetIdx());
    _pair.push_back(b->GetIdx());
  }

  OBMolPairIter::OBMolPairIter(OBMol &mol)
  {
    _parent = &mol;

    bool foundPair = false;
    OBAtom *a = _parent->BeginAtom(_i);
    if (!a)
      return;
    OBAtom *b = _parent->BeginAtom(_j);
    while (!foundPair) {
      b = _parent->NextAtom(_j);

      if (!b) {
        a = _parent->NextAtom(_i);
	if (!a)
          return;
        b = _parent->BeginAtom(_j);
      }

      if (a->GetIdx() >= b->GetIdx()) continue;
      if (a->IsConnected(b)) continue;
      if (a->IsOneThree(b)) continue;

      foundPair = true;
    }

    _pair.clear();
    _pair.push_back(a->GetIdx());
    _pair.push_back(b->GetIdx());
  }

  OBMolPairIter::OBMolPairIter(const OBMolPairIter &ai)
  {
    _parent = ai._parent;
    _pair = ai._pair;
    _i = ai._i;
    _j = ai._j;
  }

  OBMolPairIter& OBMolPairIter::operator=(const OBMolPairIter &ai)
  {
    if (this != &ai) {
      _parent = ai._parent;
      _pair = ai._pair;
      _i = ai._i;
      _j = ai._j;
    }
    return *this;
  }

  OBMolPairIter& OBMolPairIter::operator++()
  {
    _pair.clear();

    bool foundPair = false;
    OBAtom *a, *b;
    a = *_i;
    while (!foundPair) {
      b = _parent->NextAtom(_j);

      if (!b) {
        a = _parent->NextAtom(_i);
	if (!a)
          return *this;
        b = _parent->BeginAtom(_j);
      }

      if (a->GetIdx() >= b->GetIdx()) continue;
      if (a->IsConnected(b)) continue;
      if (a->IsOneThree(b)) continue;


      foundPair = true;
    }

    _pair.push_back(a->GetIdx());
    _pair.push_back(b->GetIdx());

    return *this;
  }

} // namespace OpenBabel

//! \file obiter.cpp
//! \brief STL-style iterators for Open Babel.
