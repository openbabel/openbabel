/**********************************************************************
obiter.h - STL-style iterators for Open Babel.
 
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

#ifndef OB_OBITER_H
#define OB_OBITER_H

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
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

    OBMolAtomIter() :_parent(NULL), _ptr(NULL) { }
    OBMolAtomIter(OBMol *mol);
    OBMolAtomIter(OBMol &mol);
    OBMolAtomIter(const OBMolAtomIter &ai);

    OBMolAtomIter& operator=(const OBMolAtomIter &ai);
    operator bool() const        { return _ptr != NULL; }
    OBMolAtomIter& operator++();
    OBMolAtomIter  operator++(int);
    OBAtom* operator->() const   { return _ptr;      }
    OBAtom& operator*() const    { return *_ptr;     }
  };

  //! \brief Iterate over all atoms in an OBMol in a depth-first search (DFS)
  class OBAPI OBMolAtomDFSIter {
    OBMol               *_parent;
    OBAtom              *_ptr;
    OBBitVec             _notVisited;
    std::stack<OBAtom *> _stack;
  public:

    OBMolAtomDFSIter() : _parent(NULL), _ptr(NULL) { }
    OBMolAtomDFSIter(OBMol *mol, int StartIndex=1);
    OBMolAtomDFSIter(OBMol &mol, int StartIndex=1);
    OBMolAtomDFSIter(const OBMolAtomDFSIter &ai);

    OBMolAtomDFSIter& operator=(const OBMolAtomDFSIter &ai);
    operator bool() const        { return _ptr != NULL; }
    OBMolAtomDFSIter& operator++();
    OBMolAtomDFSIter  operator++(int);
    OBAtom* operator->() const   { return _ptr;      }
    OBAtom& operator*() const    { return *_ptr;     }
    /// \return NULL if at the last atom in a fragment, else the next atom
    OBAtom* next()
    { 
      if(_stack.empty())
        return NULL; //end of a disconnected fragment
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
  public:

    OBMolAtomBFSIter(): _parent(NULL), _ptr(NULL) { }
    OBMolAtomBFSIter(OBMol *mol);
    OBMolAtomBFSIter(OBMol &mol);
    OBMolAtomBFSIter(const OBMolAtomBFSIter &ai);

    OBMolAtomBFSIter& operator=(const OBMolAtomBFSIter &ai);
    operator bool() const        { return _ptr != NULL; }
    OBMolAtomBFSIter& operator++();
    OBMolAtomBFSIter  operator++(int);
    OBAtom* operator->() const   { return _ptr;      }
    OBAtom& operator*() const    { return *_ptr;     }
  };

  //! \brief Iterate over all bonds in an OBMol
  class OBAPI OBMolBondIter {
    std::vector<OBBond*>::iterator _i;
    OBMol *_parent;
    OBBond *_ptr;
  public:

    OBMolBondIter() : _parent(NULL), _ptr(NULL) {}
    OBMolBondIter(OBMol *mol);
    OBMolBondIter(OBMol &mol);
    OBMolBondIter(const OBMolBondIter &bi);

    OBMolBondIter& operator=(const OBMolBondIter &bi);
    operator bool() const        { return _ptr != NULL; }
    OBMolBondIter& operator++();
    OBMolBondIter  operator++(int);
    OBBond* operator->() const   { return _ptr;      }
    OBBond& operator*() const    { return *_ptr;     }
  };

  //! \brief Iterate over all neighboring atoms to an OBAtom
  class OBAPI OBAtomAtomIter {
    std::vector<OBBond*>::iterator _i;
    OBAtom *_parent;
    OBAtom *_ptr;
  public:

    OBAtomAtomIter() : _parent(NULL), _ptr(NULL) { }
    OBAtomAtomIter(OBAtom *atm);
    OBAtomAtomIter(OBAtom &atm);
    OBAtomAtomIter(const OBAtomAtomIter &ai);

    OBAtomAtomIter& operator=(const OBAtomAtomIter &ai);
    operator bool() const        { return _ptr != NULL; }
    OBAtomAtomIter& operator++();
    OBAtomAtomIter  operator++(int);
    OBAtom* operator->() const   { return _ptr; }
    OBAtom& operator*() const    { return *_ptr;}
  };

  //! \brief Iterate over all bonds on an OBAtom
  class OBAPI OBAtomBondIter {
    std::vector<OBBond*>::iterator _i;
    OBAtom *_parent;
    OBBond *_ptr;
  public:

    OBAtomBondIter(): _parent(NULL), _ptr(NULL) { }
    OBAtomBondIter(OBAtom *atm);
    OBAtomBondIter(OBAtom &atm);
    OBAtomBondIter(const OBAtomBondIter &bi);

    OBAtomBondIter& operator=(const OBAtomBondIter &bi);
    operator bool() const        { return _ptr != NULL; }
    OBAtomBondIter& operator++();
    OBAtomBondIter  operator++(int);
    OBBond* operator->() const   { return _ptr; }
    OBBond& operator*() const    { return *_ptr;}
  };

  //! \brief Iterate over all residues in an OBMol
  class OBAPI OBResidueIter {
    std::vector<OBResidue*>::iterator _i;
    OBResidue *_ptr;
    OBMol *_parent;
  public:

    OBResidueIter() : _ptr(NULL), _parent(NULL) { }
    OBResidueIter(OBMol *mol);
    OBResidueIter(OBMol &mol);
    OBResidueIter(const OBResidueIter &ri);

    OBResidueIter& operator=(const OBResidueIter &ri);
    operator bool() const        { return _ptr != NULL; }
    OBResidueIter& operator++();
    OBResidueIter  operator++(int);
    OBResidue* operator->() const{ return _ptr; }
    OBResidue& operator*() const { return *_ptr;}
  };

  //! \brief Iterate over all atoms in an OBResidue
  class OBAPI OBResidueAtomIter {
    std::vector<OBAtom*>::iterator _i;
    OBResidue *_parent;
    OBAtom    *_ptr;
  public:

    OBResidueAtomIter() : _parent(NULL), _ptr(NULL) { }
    OBResidueAtomIter(OBResidue *res);
    OBResidueAtomIter(OBResidue &res);
    OBResidueAtomIter(const OBResidueAtomIter &ri);

    OBResidueAtomIter &operator = (const OBResidueAtomIter &ri);
    operator bool() const        { return _ptr != NULL; }
    OBResidueAtomIter& operator++ ();
    OBResidueAtomIter  operator++ (int);
    OBAtom *operator->() const   { return _ptr; }
    OBAtom &operator*() const    { return *_ptr;}
  };
  
  //! \brief Iterate over all angles in an OBMol
  class OBAPI OBMolAngleIter {
    OBMol     *_parent;
    std::vector<std::vector<unsigned int> > _vangle;
    std::vector<std::vector<unsigned int> >::iterator _i;
    std::vector<unsigned int> _angle;
  public:

    OBMolAngleIter() :_parent(NULL) { }
    OBMolAngleIter(OBMol *mol);
    OBMolAngleIter(OBMol &mol);
    OBMolAngleIter(const OBMolAngleIter &ai);

    OBMolAngleIter& operator=(const OBMolAngleIter &ai);
    operator bool() const        { return (_i != _vangle.end()); }
    OBMolAngleIter& operator++();
    std::vector<unsigned int> operator*() const    { return _angle;     }
  };

  //! \brief Iterate over all torsions in an OBMol
  class OBAPI OBMolTorsionIter {
    OBMol *_parent;
    std::vector<std::vector<unsigned int> > _vtorsion;
    std::vector<std::vector<unsigned int> >::iterator _i;
    std::vector<unsigned int> _torsion;
  public:

    OBMolTorsionIter() :_parent(NULL) { }
    OBMolTorsionIter(OBMol *mol);
    OBMolTorsionIter(OBMol &mol);
    OBMolTorsionIter(const OBMolTorsionIter &ai);

    OBMolTorsionIter& operator=(const OBMolTorsionIter &ai);
    operator bool() const        { return (_i != _vtorsion.end()); }
    OBMolTorsionIter& operator++();
    std::vector<unsigned int> operator*() const    { return _torsion;     }
  };
  
  //! \brief Iterate over all pairs of atoms (>1-4) in an OBMol
  class OBAPI OBMolPairIter {
    OBMol *_parent;
    std::vector<std::vector<unsigned int> > _vpair;
    std::vector<std::vector<unsigned int> >::iterator _i;
    std::vector<unsigned int> _pair;
 
  public:

    OBMolPairIter() :_parent(NULL) { }
    OBMolPairIter(OBMol *mol);
    OBMolPairIter(OBMol &mol);
    OBMolPairIter(const OBMolPairIter &ai);

    OBMolPairIter& operator=(const OBMolPairIter &ai);
    operator bool() const        { return (_i != _vpair.end()); }
    OBMolPairIter& operator++();
    std::vector<unsigned int> operator*() const    { return _pair;     }
  };


  //! \brief Iterate over all residues in an OBMol
  class OBAPI OBMolRingIter {
    std::vector<OBRing*>::iterator _i;
    OBRing *_ptr;
    OBMol *_parent;
    OBRingData *_rings;
  public:

    OBMolRingIter() : _ptr(NULL), _parent(NULL), _rings(NULL) { }
    OBMolRingIter(OBMol *mol);
    OBMolRingIter(OBMol &mol);
    OBMolRingIter(const OBMolRingIter &ri);

    OBMolRingIter& operator=(const OBMolRingIter &ri);
    operator bool()      const { return _ptr != NULL; }
    OBMolRingIter& operator++();
    OBMolRingIter  operator++(int);
    OBRing* operator->() const { return _ptr; }
    OBRing& operator*()  const { return *_ptr;}
  };

#define FOR_ATOMS_OF_MOL(a,m)     for( OBMolAtomIter     a(m); a; ++a )
#define FOR_BONDS_OF_MOL(b,m)     for( OBMolBondIter     b(m); b; ++b )
#define FOR_NBORS_OF_ATOM(a,p)    for( OBAtomAtomIter    a(p); a; ++a )
#define FOR_BONDS_OF_ATOM(b,p)    for( OBAtomBondIter    b(p); b; ++b )
#define FOR_RESIDUES_OF_MOL(r,m)  for( OBResidueIter     r(m); r; ++r )
#define FOR_ATOMS_OF_RESIDUE(a,r) for( OBResidueAtomIter a(r); a; ++a )
#define FOR_DFS_OF_MOL(a,m)       for( OBMolAtomDFSIter  a(m); a; ++a )
#define FOR_BFS_OF_MOL(a,m)       for( OBMolAtomBFSIter  a(m); a; ++a )
#define FOR_RINGS_OF_MOL(a,m)     for( OBMolRingIter     r(m); r; ++r )
#define FOR_ANGLES_OF_MOL(a,m)    for( OBMolAngleIter    a(m); a; ++a )
#define FOR_TORSIONS_OF_MOL(t,m)  for( OBMolTorsionIter  t(m); t; ++t )
#define FOR_PAIRS_OF_MOL(p,m)     for( OBMolPairIter     p(m); p; ++p )

} // namespace OpenBabel
#endif // OB_OBITER_H

//! \file obiter.h
//! \brief STL-style iterators for Open Babel.
