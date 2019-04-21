/**********************************************************************
ring.h - Deal with rings, find smallest set of smallest rings (SSSR).

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

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

// TODO: Make this work as a free-standing header
// Currently only used in ring.cpp which imports mol.h beforehand
#include <openbabel/bitvec.h>
#include <openbabel/typer.h>

namespace OpenBabel
{

  class OBMol;
  class OBAtom;
  class OBBond;
  class vector3;

  // class introduction in ring.cpp
  class OBAPI OBRing
  {
    OBMol *_parent; //!< parent molecule for this ring
  public:
    //public data members
    int ring_id;            //!< a unique id to ensure a stable sort in SSSR determination
    std::vector<int> _path; //!< the path of this ring (i.e., the atom indexes)
    OBBitVec _pathset;      //!< the path of this ring as a redundant bit vector

    //! \name Constructors
    //@{
    OBRing()    {}
    //! Initialize a ring from a set of atom indexes @p path and with @p size
    OBRing(std::vector<int>& path, int size);
    OBRing(std::vector<int>& path, OBBitVec set) : _path(path), _pathset(set) {}
    OBRing(const OBRing &src);
    OBRing& operator=(const OBRing &src);
    //@}

    //member functions

    //! \return the size of this ring (i.e., how many atoms in the cycle)
    size_t    Size()     const  {    return(_path.size());  }
    //! \return the size of this ring (i.e., how many atoms in the cycle)
    //! \deprecated Use Size() instead
    size_t    PathSize() const  {    return(_path.size());  }

    //! \return whether this ring is aromatic
    //! If all atoms in this ring are aromatic, the ring will be considered aromatic
    //! \todo This method uses implicit bonding -- bond info is not stored in OBRing
    bool   IsAromatic();

    //! Set the ring type (see OBRingTyper for more)
    void SetType(char *type);
    //! Set the ring type (see OBRingTyper for more)
    void SetType(std::string &type);
    //! \return the ring type
    char *GetType();
    //! \return the index for the root atom. O for furan, S for thiazole, N for pyrrole.
    //! For 6 membered aromatic rings, the first non carbon atom is used for root.
    //! For 5 members rings the O, S or N (BOSum=3, valence=3) will be used for root
    unsigned int GetRootAtom();

    //! \return Whether atom @p a is a member of this ring
    bool   IsMember(OBAtom *a);
    //! \return Whether both atoms in bond @p b are in this ring
    //! \todo This method uses implicit bonding -- bond info is not stored in OBRing
    bool	 IsMember(OBBond *b);
    //! \return Whether @p i as an atom index is in this ring
    bool   IsInRing(int i)
    {
      return(_pathset.BitIsSet(i));
    }

    //! Set the parent of this ring to @p m
    void   SetParent(OBMol *m)  {    _parent = m;    }
    //! \return the parent of this ring, or NULL if none has been defined
    OBMol *GetParent()          {    return(_parent);}

    //! Set the supplied vectors to the @p center of this ring, along with
    //! the @p normal (in both directions).
    //! \param center The center of the ring
    //! \param norm1 The normal of the best-fit plane for this ring
    //! \param norm2 -1 * norm1 (i.e., the opposite direction of norm1)
    //! \return True (success)
    bool findCenterAndNormal(vector3 & center, vector3 &norm1, vector3 &norm2);
  private:
    char _type[30];   //!< ring type
  };

  //! Comparison function for rings, used by OBRingSearch::SortRings()
  //! \return true if a.size() > b.size()
  OBAPI bool CompareRingSize(const OBRing *,const OBRing *);


  /** \class OBRingSearch ring.h <openbabel/ring.h>
      \brief Internal class to facilitate OBMol::FindSSSR()
  **/
  class OBAPI OBRingSearch
  {
    std::vector<OBBond*> _bonds; //!< the internal list of closure bonds (deprecated)
    std::vector<OBRing*> _rlist; //!< the internal list of rings
  public:
    OBRingSearch()    {}
    ~OBRingSearch();

    //! Sort ring sizes from smallest to largest
    void    SortRings()
    {
      std::vector<OBRing*>::iterator j;
      int ring_id; // for each ring, assign a unique id to ensure a stable sort
      
      for (j = _rlist.begin(), ring_id = 0; j != _rlist.end(); ++j, ++ring_id)
        (*j)->ring_id = ring_id;
      std::sort(_rlist.begin(),_rlist.end(),CompareRingSize);
    }
    //! Starting with a full ring set - reduce to SSSR set
    void    RemoveRedundant(int);
    //! Add a new ring from a "closure" bond: See OBBond::IsClosure()
    void    AddRingFromClosure(OBMol &,OBBond *);

    bool    SaveUniqueRing(std::deque<int>&,std::deque<int>&);

    //! For debugging only, write the rings to std::cout
    void    WriteRings();

    //! \name Iterator methods -- see OBMolRingIter for iteration over a molecule
    //@{
    //! \return an iterator pointing to the beginning of the list of rings
    std::vector<OBRing*>::iterator BeginRings()
      {
        return(_rlist.begin());
      }
    //! \return an iterator pointing to the end of the list of rings
    std::vector<OBRing*>::iterator EndRings()
      {
        return(_rlist.end());
      }
    //@}
  };

  /** \class OBRTree ring.h <openbabel/ring.h>
      \brief Internal class for OBRing search algorithms to create a search tree
      of OBAtom objects
  **/
  class OBAPI OBRTree
  {
    OBAtom  *_atom; //!< Atom represented by this node in the tree
    OBRTree *_prv;  //!< Previous (parent) entry in an OBRing tree
  public:
    //! Construct a search tree from a possible parent entry and atom entry
    OBRTree(OBAtom*,OBRTree*);
    ~OBRTree()    {}

    //! \return the OBAtom::GetIdx() index of the atom in this node
    int  GetAtomIdx();
    //! Recursively find the root of this tree, building up a vector of OBAtom nodes.
    void PathToRoot(std::vector<OBAtom*>&);
  };

} // end namespace OpenBabel

#endif // OB_RING_H

//! \file ring.h
//! \brief Deal with rings, find smallest set of smallest rings (SSSR).
