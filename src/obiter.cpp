/**********************************************************************
obiter.cpp - STL-style iterators for Open Babel
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
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
#include "babelconfig.h"

#include <vector>

#include "mol.h"
#include "obiter.h"

using namespace std;

namespace OpenBabel
{

  /** \class OBMolAtomIter

  To facilitate iteration through all atoms in a molecule, without resorting
  to atom indexes (which <strong>will</strong> change in the future), a 
  variety of iterator methods are provided.

  This has been made significantly easier by a series of macros in the 
  obiter.h header file:

  \code
  \#define FOR_ATOMS_OF_MOL(a,m)     for( OBMolAtomIter     a(m); a; a++ )
  \endcode

  Here is an example:
  \code
  #include "obiter.h"
  #include "mol.h"

  OBMol mol;
  double exactMass = 0.0f;
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

  OBMolAtomIter OBMolAtomIter::operator++(int)
  {
    _ptr = _parent->NextAtom(_i);
    return *this;
  }

  /** \class OBMolAtomDFSIter

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
  \#define FOR_DFS_OF_MOL(a,m)     for( OBMolAtomDFSIter     a(m); a; a++ )
  \endcode

  Here is an example:
  \code
  #include "obiter.h"
  #include "mol.h"

  OBMol mol;
  FOR_DFS_OF_MOL(a, mol)
  {
  // The variable a behaves like OBAtom* when used with -> and * but
  // but needs to be explicitly converted when appearing as a parameter
  // in a function call - use &*a

  }
  \endcode
  **/

  OBMolAtomDFSIter::OBMolAtomDFSIter(OBMol *mol):
    _parent(mol), _ptr(_parent->GetFirstAtom())
  {
    _notVisited.Resize(_parent->NumAtoms());
    _notVisited.Negate(); // all on
    _notVisited.SetBitOff(_ptr->GetIdx() - 1);

    vector<OBEdgeBase*>::iterator i;
    OBAtom *a;

    for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
      {
        _stack.push(a);
        _notVisited.SetBitOff(a->GetIdx() - 1);
      }
  }

  OBMolAtomDFSIter::OBMolAtomDFSIter(OBMol &mol):
    _parent(&mol), _ptr(_parent->GetFirstAtom())
  {
    _notVisited.Resize(_parent->NumAtoms());
    _notVisited.Negate(); // all on
    _notVisited.SetBitOff(_ptr->GetIdx() - 1);

    vector<OBEdgeBase*>::iterator i;
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

  OBMolAtomDFSIter OBMolAtomDFSIter::operator++(int)
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
      }

    if (_ptr)
      {
        vector<OBEdgeBase*>::iterator i;
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

  /** \class OBMolAtomBFSIter

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
  \#define FOR_BFS_OF_MOL(a,m)     for( OBMolAtomBFSIter     a(m); a; a++ )
  \endcode

  Here is an example:
  \code
  #include "obiter.h"
  #include "mol.h"

  OBMol mol;
  FOR_BFS_OF_MOL(a, mol)
  {
  // The variable a behaves like OBAtom* when used with -> and * but
  // but needs to be explicitly converted when appearing as a parameter
  // in a function call - use &*a
  
  }
  \endcode
  **/

  OBMolAtomBFSIter::OBMolAtomBFSIter(OBMol *mol):
    _parent(mol), _ptr(_parent->GetFirstAtom())
  {
    _notVisited.Resize(_parent->NumAtoms());
    _notVisited.Negate(); // all on
    _notVisited.SetBitOff(_ptr->GetIdx() - 1);

    vector<OBEdgeBase*>::iterator i;
    OBAtom *a;

    for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
      {
        _queue.push(a);
        _notVisited.SetBitOff(a->GetIdx() - 1);
      }
  }

  OBMolAtomBFSIter::OBMolAtomBFSIter(OBMol &mol):
    _parent(&mol), _ptr(_parent->GetFirstAtom())
  {
    _notVisited.Resize(_parent->NumAtoms());
    _notVisited.Negate(); // all on
    _notVisited.SetBitOff(_ptr->GetIdx() - 1);

    vector<OBEdgeBase*>::iterator i;
    OBAtom *a;

    for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
      {
        _queue.push(a);
        _notVisited.SetBitOff(a->GetIdx() - 1);
      }
  }

  OBMolAtomBFSIter::OBMolAtomBFSIter(const OBMolAtomBFSIter &ai)
  {
    _parent = ai._parent;
    _ptr = ai._ptr;
    _notVisited = ai._notVisited;
    _queue = ai._queue;
  }

  OBMolAtomBFSIter& OBMolAtomBFSIter::operator=(const OBMolAtomBFSIter &ai)
  {
    if (this != &ai)
      {
        _parent = ai._parent;
        _ptr = ai._ptr;
        _notVisited = ai._notVisited;
        _queue = ai._queue;
      }
    return *this;
  }

  OBMolAtomBFSIter OBMolAtomBFSIter::operator++(int)
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
            _ptr = _parent->GetAtom(next + 1);
            _notVisited.SetBitOff(next);
          }
      }

    if (_ptr)
      {
        vector<OBEdgeBase*>::iterator i;
        OBAtom *a;
        
        for (a = _ptr->BeginNbrAtom(i); a; a = _ptr->NextNbrAtom(i))
          if (_notVisited[a->GetIdx() - 1])
            {
              _queue.push(a);
              _notVisited.SetBitOff(a->GetIdx() - 1);
            }
      }

    return *this;
  }

  /** \class OBMolBondIter

  To facilitate iteration through all bonds in a molecule, without resorting
  to bond indexes (which may change in the future), a variety of
  iterators are provided.

  This has been made significantly easier by a series of macros in the 
  obiter.h header file:

  \code
  \#define FOR_BONDS_OF_MOL(b,m)     for( OBMolBondIter     b(m); b; b++ )
  \endcode

  Here is an example:
  \code
  #include "obiter.h"
  #include "mol.h"

  OBMol mol;
  unsigned int bondOrderSum = 0;
  FOR_BONDS_OF_MOL(b, mol)
  {
  // The variable b behaves like OBBond* when used with -> and * but
  // but needs to be explicitly converted when appearing as a parameter
  // in a function call - use &*b
  bondOrderSum +=  b->GetBO();
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

  OBMolBondIter OBMolBondIter::operator++(int)
  {
    _ptr = _parent->NextBond(_i);
    return *this;
  }

  /** \class OBAtomAtomIter

  To facilitate iteration through all neighbors of an atom, without resorting
  to bond indexes (which may change in the future), a variety of
  iterator classes and methods are provided.

  This has been made significantly easier by a series of macros in the 
  obiter.h header file:

  \code
  \#define FOR_NBORS_OF_ATOM(a,p)     for( OBAtomAtomIter     a(p); a; a++ )
  \endcode

  Here is an example:
  \code
  #include "obiter.h"
  #include "mol.h"

  OBMol mol;
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
    _ptr = _ptr->BeginNbrAtom(_i);
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

  OBAtomAtomIter OBAtomAtomIter::operator++(int)
  {
    _ptr = _parent->NextNbrAtom(_i);
    return *this;
  }

  /** \class OBAtomBondIter

  To facilitate iteration through all bonds on an atom, without resorting
  to bond indexes (which may change in the future) a variety of
  iterator classes and methods are provided.

  This has been made significantly easier by a series of macros in the 
  obiter.h header file:

  \code
  \#define FOR_BONDS_OF_ATOM(b,p)     for( OBAtomBondIter     b(p); b; b++ )
  \endcode

  Here is an example:
  \code
  #include "obiter.h"
  #include "mol.h"

  OBAtom atom;
  unsigned int tripleBondCount;
  FOR_BONDS_OF_ATOM(b, atom)
  {
  // The variable b behaves like OBBond* when used with -> and * but
  // but needs to be explicitly converted when appearing as a parameter
  // in a function call - use &*b
  if (b->GetBO() == 3)
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

  OBAtomBondIter OBAtomBondIter::operator++(int)
  {
    _ptr = _parent->NextBond(_i);
    return *this;
  }

  /** \class OBResidueIter

  To facilitate iteration through all residues in a molecule, without resorting
  to residue indexes (which may change in the future) a variety of
  iterator classes and methods are provided.

  This has been made significantly easier by a series of macros in the 
  obiter.h header file:

  \code
  \#define FOR_RESIDUES_OF_MOL(a,m)     for( OBResidueIter     a(m); a; a++ )
  \endcode

  Here is an example:
  \code
  #include "obiter.h"
  #include "mol.h"

  OBMol mol;
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

  OBResidueIter OBResidueIter::operator++(int)
  {
    _ptr = _parent->NextResidue(_i);
    return *this;
  }

  /** \class OBResidueAtomIter

  To facilitate iteration through all atoms in a residue, without resorting
  to atom indexes (which may change in the future) a variety of
  iterator classes and methods are provided.

  This has been made significantly easier by a series of macros in the 
  obiter.h header file:

  \code
  \#define FOR_ATOMS_OF_RESIDUE(a,r)     for( OBResidueAtomIter     a(r); a; a++ )
  \endcode

  Here is an example:
  \code
  #include "obiter.h"
  #include "mol.h"

  OBMol mol;
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

  OBResidueAtomIter::OBResidueAtomIter(OBResidue *res)
  {
    _parent = res;
    _ptr    = _parent->BeginAtom(_i);
  }

  OBResidueAtomIter::OBResidueAtomIter(OBResidue &res)
  {
    _parent = &res;
    _ptr    = _parent->BeginAtom(_i);
  }

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

  OBResidueAtomIter OBResidueAtomIter::operator++ (int)
  {
    _ptr = _parent->NextAtom(_i);
    return (*this);
  }

} // namespace OpenBabel

//! \file obiter.cpp
//! \brief STL-style iterators for Open Babel.
