/**********************************************************************
obiter.cpp - STL-style iterators for Open Babel
 
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

#include <vector>

#include "mol.h"
#include "obiter.h"

using namespace std;

namespace OpenBabel
{

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
