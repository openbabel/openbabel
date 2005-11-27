/**********************************************************************
ring.h - Deal with rings, find smallest set of smallest rings (SSSR).
 
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

#ifndef OB_RING_H
#define OB_RING_H

#include <deque>
#include <algorithm>

namespace OpenBabel
{

class OBMol;
class OBAtom;
class OBBond;

//! Internal class for OBRing search algorithms to create a search tree of OBAtom objects
class OBAPI OBRTree
{
  OBAtom  *_atom; //!< Atom represented by this node in the tree
  OBRTree *_prv; //!< Previous (parent) entry in an OBRing tree
public:
  OBRTree(OBAtom*,OBRTree*);
  ~OBRTree()    {}
  
  //! \return the OBAtom::GetIdx() index of the atom in this node
  int  GetAtomIdx();
  //! Recursively find the root of this tree, building up a vector of OBAtom nodes.
  void PathToRoot(std::vector<OBNodeBase*>&);
};

// class introduction in ring.cpp
class OBAPI OBRing
{
    OBMol *_parent;
public:
    //public data members
    std::vector<int> _path;
    OBBitVec _pathset;
    bool findCenterAndNormal(vector3 & center, vector3 &norm1, vector3 &norm2);

    //constructors
    OBRing()    {}
    OBRing(std::vector<int>&,int);
    OBRing(const OBRing &src);
    OBRing& operator=(const OBRing &src);

    //member functions
    int    Size()     const
    {
        return(_path.size());
    }
    int    PathSize() const
    {
        return(_path.size());
    }
    bool   IsMember(OBAtom *a);
    bool	 IsMember(OBBond *b);
    bool   IsAromatic();
    bool   IsInRing(int i)
    {
        return(_pathset.BitIsOn(i));
    }
    void   SetParent(OBMol *m)
    {
        _parent = m;
    }
    OBMol *GetParent()
    {
        return(_parent);
    }
};

bool CompareRingSize(const OBRing *,const OBRing *);


//! Internal class to facilitate OBMol::FindSSSR()
class OBAPI OBRingSearch
{
    std::vector<OBBond*> _bonds;
    std::vector<OBRing*> _rlist;
public:
    OBRingSearch()    {}
    ~OBRingSearch();

    //! Sort ring sizes from smallest to largest
    void    SortRings()
    {
        std::sort(_rlist.begin(),_rlist.end(),CompareRingSize);
    }
    //! Starting with a full ring set - reduce to SSSR set
    void    RemoveRedundant(int);
    void    AddRingFromClosure(OBMol &,OBBond *);
     //! For debugging only, write the rings to std::cout
    void    WriteRings();

    bool    SaveUniqueRing(std::deque<int>&,std::deque<int>&);

    std::vector<OBRing*>::iterator BeginRings()
    {
        return(_rlist.begin());
    }
    std::vector<OBRing*>::iterator EndRings()
    {
        return(_rlist.end());
    }
};

} // end namespace OpenBabel

#endif // OB_RING_H

//! \file ring.h
//! \brief Deal with rings, find smallest set of smallest rings (SSSR).
