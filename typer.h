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
  int                                            _rc;
  vector<vector<int> >                           _mlist;
  vector<pair<OBSmartsPattern*,int> >            _vinthyb;
  vector<pair<OBSmartsPattern*,int> >            _vimpval;
  vector<pair<OBSmartsPattern*,string> >         _vexttyp;
 public:
  OBAtomTyper();
  ~OBAtomTyper();

  void ParseLine(char*);
  void AssignHyb(OBMol&);
  void AssignTypes(OBMol&);
  void AssignImplicitValence(OBMol&);
  void CorrectAromaticNitrogens(OBMol&);
};

class OBAromaticTyper : public OBGlobalDataBase
{
  vector<bool>             _vpa;   //potentially aromatic atoms
  vector<bool>             _visit;
  vector<bool>             _root;
  vector<vector<int> >     _mlist;
  vector<OBSmartsPattern*> _vsp;   //smarts of potentially aromatic atoms
  vector<pair<int,int> >   _verange; //min and max number of electrons
  vector<pair<int,int> >   _velec;   //num electrons an atom contributes
 public:
  OBAromaticTyper();
  ~OBAromaticTyper();

  void ParseLine(char*);
  void AssignAromaticFlags(OBMol &);
  void PropagatePotentialAromatic(OBAtom*);
  void ExcludeSmallRing(OBMol &);
  void CheckAromaticity(OBAtom*,int);
  bool TraverseCycle(OBAtom*,OBAtom*,OBBond*,pair<int,int>&,int);
};

}

#endif //_AROM_H
