/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2002 by Geoffrey R. Hutchison

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

#ifndef OB_RINGS_H__
#define OB_RINGS_H__

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
/** Ring class
*
* Ring information beyond atom and bond membership is usually
* not necessary, but more information can be had about the rings
* in a molecule. The OBRing class is used to store the information
* from a Smallest Set of Smallest Rings (SSSR) perception of a molecule.
* The OBMol member function GetSSSR() stores the information it perceives
* in a vector<OBRing*> inside the molecule. Perception is only done once
* for a molecule unless the connection table is modified. The following
* code demonstrates how to extract the SSSR information:
\code
   OEMol mol(SDF,SDF);
   cin >> mol;
   vector<OERing*> vr;
   vr = mol.GetSSSR();
\endcode
* OBRing store the atom numbers of the atoms in each of the smallest
* set of smallest rings in both a vector<int> and an OBBitVec.
* An example of how to print out the atom numbers present in all SSSR
* rings is show below.
\code
   vector<OERing*>::iterator i;
   vector<int>::iterator j;
   vector<OERing*> *rlist = (vector<OERing*>*)mol.GetData("RingList");
   for (i = rlist->begin();i != rlist->end();i++)
   {
       for(j = (*i)->_path.begin();j != (*i)->_path.end();j++)
           cout << *j << ` `;
          cout << endl;
   }
\endcode
 * will produce something like the following output for benzene
\code
   1 2 3 4 5 6
\endcode
* Ring information is automatically deleted from an OBMol when it
* goes out of scope or the Clear() member function is called.
*/
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

} // end namespace OpenBabel

#endif // OB_RINGS_H__
