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

class OBMol;
class OBAtom; 
class OBBond;

class OBRTree
{
  OBAtom  *_atom;
  OBRTree *_prv;
 public:
  OBRTree(OBAtom*,OBRTree*);
  ~OBRTree() {}
  int  GetAtomIdx();
  void PathToRoot(std::vector<OBNodeBase*>&);
};

class OBRing
{
  OBMol *_parent;
 public:
  //public data members
  std::vector<int> _path;
  OBBitVec _pathset;
  bool findCenterAndNormal(Vector & center, Vector &norm1, Vector &norm2);

  //constructors
  OBRing(){};
  OBRing(std::vector<int>&,int);
	OBRing(const OBRing &src);
	OBRing& operator=(const OBRing &src);

  //member functions
  int    Size()     const     {return(_path.size());}
  int    PathSize() const     {return(_path.size());}
  bool   IsMember(OBAtom *a);
	bool	 IsMember(OBBond *b);
  bool   IsAromatic();
  bool   IsInRing(int i)      {return(_pathset.BitIsOn(i));}
  void   SetParent(OBMol *m)  {_parent = m;}
  OBMol *GetParent()          {return(_parent);}
};

bool CompareRingSize(const OBRing *,const OBRing *);

class OBRingSearch
{
  std::vector<OBBond*> _bonds;
  std::vector<OBRing*> _rlist;
 public:
  OBRingSearch(){}
  ~OBRingSearch();
  void    SortRings() {sort(_rlist.begin(),_rlist.end(),CompareRingSize);}
  void    RemoveRedundant(int);
  void    AddRingFromClosure(OBMol &,OBBond *,int);
  void    WriteRings();
  bool    SaveUniqueRing(std::deque<int>&,std::deque<int>&);
  std::vector<OBRing*>::iterator BeginRings() {return(_rlist.begin());}
  std::vector<OBRing*>::iterator EndRings() {return(_rlist.end());}
};

}

#endif //__RINGS_H__
