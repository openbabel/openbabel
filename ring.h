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

#ifndef __RINGS_H__
#define __RINGS_H__

#include <deque>

namespace OpenBabel {

class OEMol;
class OEAtom; 
class OEBond;

class OERTree
{
  OEAtom  *_atom;
  OERTree *_prv;
 public:
  OERTree(OEAtom*,OERTree*);
  ~OERTree() {}
  int  GetAtomIdx();
  void PathToRoot(vector<OENodeBase*>&);
};

class OERing
{
  OEMol *_parent;
 public:
  //public data members
  vector<int> _path;
  OEBitVec _pathset;
  bool findCenterAndNormal(Vector & center, Vector &norm1, Vector &norm2);

  //constructors
  OERing(){};
  OERing(vector<int>&,int);
	OERing(const OERing &src);
	OERing& operator=(const OERing &src);

  //member functions
  int    Size()     const     {return(_path.size());}
  int    PathSize() const     {return(_path.size());}
  bool   IsMember(OEAtom *a);
	bool	 IsMember(OEBond *b);
  bool   IsAromatic();
  bool   IsInRing(int i)      {return(_pathset.BitIsOn(i));}
  void   SetParent(OEMol *m)  {_parent = m;}
  OEMol *GetParent()          {return(_parent);}
};

bool CompareRingSize(const OERing *,const OERing *);

class OERingSearch
{
  vector<OEBond*> _bonds;
  vector<OERing*> _rlist;
 public:
  OERingSearch(){}
  ~OERingSearch();
  void    SortRings() {sort(_rlist.begin(),_rlist.end(),CompareRingSize);}
  void    RemoveRedundant(int);
  void    AddRingFromClosure(OEMol &,OEBond *,int);
  void    WriteRings();
  bool    SaveUniqueRing(deque<int>&,deque<int>&);
  vector<OERing*>::iterator BeginRings() {return(_rlist.begin());}
  vector<OERing*>::iterator EndRings() {return(_rlist.end());}
};

}

#endif //__RINGS_H__
