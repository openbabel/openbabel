/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef _AROM_H
#define _AROM_H

#include "parsmart.h"
#include "data.h"

namespace OpenBabel {

class OBAtomTyper : public OBGlobalDataBase
{
  int                                                      _rc;
  std::vector<std::vector<int> >                           _mlist;
  std::vector<std::pair<OBSmartsPattern*,int> >            _vinthyb;
  std::vector<std::pair<OBSmartsPattern*,int> >            _vimpval;
  std::vector<std::pair<OBSmartsPattern*,std::string> >    _vexttyp;
 public:
  OBAtomTyper();
  ~OBAtomTyper();

  void ParseLine(const char*);
  void AssignHyb(OBMol&);
  void AssignTypes(OBMol&);
  void AssignImplicitValence(OBMol&);
  void CorrectAromaticNitrogens(OBMol&);
};

class OBAromaticTyper : public OBGlobalDataBase
{
  std::vector<bool>             _vpa;   //potentially aromatic atoms
  std::vector<bool>             _visit;
  std::vector<bool>             _root;
  std::vector<std::vector<int> >     _mlist;
  std::vector<OBSmartsPattern*> _vsp;   //smarts of potentially aromatic atoms
  std::vector<std::pair<int,int> >   _verange; //min and max number of electrons
  std::vector<std::pair<int,int> >   _velec;   //num electrons an atom contributes
 public:
  OBAromaticTyper();
  ~OBAromaticTyper();

  void ParseLine(const char*);
  void AssignAromaticFlags(OBMol &);
  void PropagatePotentialAromatic(OBAtom*);
  void ExcludeSmallRing(OBMol &);
  void CheckAromaticity(OBAtom*,int);
  bool TraverseCycle(OBAtom*,OBAtom*,OBBond*,std::pair<int,int>&,int);
};

}

#endif //_AROM_H
