/**********************************************************************
obiter.h - STL-style iterators for Open Babel.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#include <vector>

namespace OpenBabel
{

  class OBMol;
  class OBAtom;
  class OBBond;
  class OBResidue;

  //! \brief Iterate over all atoms in an OBMol
class OBAPI OBMolAtomIter {
    std::vector<OBNodeBase*>::iterator _i;
    OBMol *_parent;
    OBAtom *_ptr;
public:

    OBMolAtomIter()              { _parent = NULL; _ptr = NULL; }
    OBMolAtomIter(OBMol *mol);
    OBMolAtomIter(OBMol &mol);
    OBMolAtomIter(const OBMolAtomIter &ai);

    OBMolAtomIter& operator=(const OBMolAtomIter &ai);
    operator bool() const        { return _ptr != NULL; }
    OBMolAtomIter operator++(int);
    OBAtom* operator->() const   { return _ptr;      }
    OBAtom& operator*() const    { return *_ptr;     }
};

  //! \brief Iterate over all bonds in an OBMol
class OBAPI OBMolBondIter {
    std::vector<OBEdgeBase*>::iterator _i;
    OBMol *_parent;
    OBBond *_ptr;
public:

    OBMolBondIter()              { _parent = NULL; _ptr = NULL; }
    OBMolBondIter(OBMol *mol);
    OBMolBondIter(OBMol &mol);
    OBMolBondIter(const OBMolBondIter &bi);

    OBMolBondIter& operator=(const OBMolBondIter &bi);
    operator bool() const        { return _ptr != NULL; }
    OBMolBondIter operator++(int);
    OBBond* operator->() const   { return _ptr;      }
    OBBond& operator*() const    { return *_ptr;     }
};

  //! \brief Iterate over all neighboring atoms to an OBAtom
class OBAPI OBAtomAtomIter {
    std::vector<OBEdgeBase*>::iterator _i;
    OBAtom *_parent;
    OBAtom *_ptr;
public:

    OBAtomAtomIter()             { _parent = NULL; _ptr = NULL; }
    OBAtomAtomIter(OBAtom *atm);
    OBAtomAtomIter(OBAtom &atm);
    OBAtomAtomIter(const OBAtomAtomIter &ai);

    OBAtomAtomIter& operator=(const OBAtomAtomIter &ai);
    operator bool() const        { return _ptr != NULL; }
    OBAtomAtomIter operator++(int);
    OBAtom* operator->() const   { return _ptr; }
    OBAtom& operator*() const    { return *_ptr;}
};

  //! \brief Iterate over all bonds on an OBAtom
class OBAPI OBAtomBondIter {
    std::vector<OBEdgeBase*>::iterator _i;
    OBAtom *_parent;
    OBBond *_ptr;
public:

    OBAtomBondIter()             { _parent = NULL; _ptr = NULL; }
    OBAtomBondIter(OBAtom *atm);
    OBAtomBondIter(OBAtom &atm);
    OBAtomBondIter(const OBAtomBondIter &bi);

    OBAtomBondIter& operator=(const OBAtomBondIter &bi);
    operator bool() const        { return _ptr != NULL; }
    OBAtomBondIter operator++(int);
    OBBond* operator->() const   { return _ptr; }
    OBBond& operator*() const    { return *_ptr;}
};

  //! \brief Iterate over all residues in an OBMol
class OBAPI OBResidueIter {
    std::vector<OBResidue*>::iterator _i;
    OBResidue *_ptr;
    OBMol *_parent;
public:

    OBResidueIter()              { _parent = NULL; _ptr = NULL; }
    OBResidueIter(OBMol *mol);
    OBResidueIter(OBMol &mol);
    OBResidueIter(const OBResidueIter &ri);

    OBResidueIter& operator=(const OBResidueIter &ri);
    operator bool() const        { return _ptr != NULL; }
    OBResidueIter operator++(int);
    OBResidue* operator->() const{ return _ptr; }
    OBResidue& operator*() const { return *_ptr;}
};

  //! \brief Iterate over all atoms in an OBResidue
class OBAPI OBResidueAtomIter {
    std::vector<OBAtom*>::iterator _i;
    OBResidue *_parent;
    OBAtom    *_ptr;
public:

    OBResidueAtomIter()          { _parent = NULL; _ptr = NULL; }
    OBResidueAtomIter(OBResidue *res);
    OBResidueAtomIter(OBResidue &res);
    OBResidueAtomIter(const OBResidueAtomIter &ri);

    OBResidueAtomIter &operator = (const OBResidueAtomIter &ri);
    operator bool() const        { return _ptr != NULL; }
    OBResidueAtomIter operator++ (int);
    OBAtom *operator->() const   { return _ptr; }
    OBAtom &operator*() const    { return *_ptr;}
};

#define FOR_ATOMS_OF_MOL(a,m)     for( OBMolAtomIter     a(m); a; a++ )
#define FOR_BONDS_OF_MOL(b,m)     for( OBMolBondIter     b(m); b; b++ )
#define FOR_NBORS_OF_ATOM(a,p)    for( OBAtomAtomIter    a(p); a; a++ )
#define FOR_BONDS_OF_ATOM(b,p)    for( OBAtomBondIter    b(p); b; b++ )
#define FOR_RESIDUES_OF_MOL(r,m)  for( OBResidueIter     r(m); r; r++ )
#define FOR_ATOMS_OF_RESIDUE(a,r) for( OBResidueAtomIter a(r); a; a++ )

} // namespace OpenBabel
#endif // OB_OBITER_H

//! \file obiter.h
//! \brief STL-style iterators for Open Babel.
