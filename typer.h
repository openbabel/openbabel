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

class OEAtomTyper : public OEGlobalDataBase
{
  int                                            _rc;
  vector<vector<int> >                           _mlist;
  vector<pair<OESmartsPattern*,int> >            _vinthyb;
  vector<pair<OESmartsPattern*,int> >            _vimpval;
  vector<pair<OESmartsPattern*,string> >         _vexttyp;
 public:
  OEAtomTyper();
  ~OEAtomTyper();

  void ParseLine(char*);
  void AssignHyb(OEMol&);
  void AssignTypes(OEMol&);
  void AssignImplicitValence(OEMol&);
  void CorrectAromaticNitrogens(OEMol&);
};

class OEAromaticTyper : public OEGlobalDataBase
{
  vector<bool>             _vpa;   //potentially aromatic atoms
  vector<bool>             _visit;
  vector<bool>             _root;
  vector<vector<int> >     _mlist;
  vector<OESmartsPattern*> _vsp;   //smarts of potentially aromatic atoms
  vector<pair<int,int> >   _verange; //min and max number of electrons
  vector<pair<int,int> >   _velec;   //num electrons an atom contributes
 public:
  OEAromaticTyper();
  ~OEAromaticTyper();

  void ParseLine(char*);
  void AssignAromaticFlags(OEMol &);
  void PropagatePotentialAromatic(OEAtom*);
  void ExcludeSmallRing(OEMol &);
  void CheckAromaticity(OEAtom*,int);
  bool TraverseCycle(OEAtom*,OEAtom*,OEBond*,pair<int,int>&,int);
};

}

#endif //_AROM_H
