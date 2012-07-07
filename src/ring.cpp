/**********************************************************************
ring.cpp - Deal with rings, find smallest set of smallest rings (SSSR).

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h> // implements some OBMol methods
#include <openbabel/ring.h>

using namespace std;

namespace OpenBabel
{
  OBRingTyper      ringtyper;

  /*! \class OBRing ring.h <openbabel/ring.h>
    \brief Stores information on rings in a molecule from SSSR perception.

    Ring information beyond atom and bond membership is usually not
    necessary, but more information can be had about the rings in a
    molecule. The OBRing class is used to store the information from a
    Smallest Set of Smallest Rings (SSSR) perception of a molecule. The
    OBMol member function OBMol::GetSSSR() stores the information it perceives in
    a vector<OBRing*> inside the molecule. Perception is only done once
    for a molecule unless the connection table is modified. The following
    code demonstrates how to extract the SSSR information:

    \code
    OBMol mol;

    vector<OBRing*> vr;
    vr = mol.GetSSSR();
    \endcode

    OBRings store the atom numbers of the atoms in each of the smallest
    set of smallest rings in both a vector<int> and an OBBitVec.
    An example of how to print out the atom numbers present in all SSSR
    rings is show below:

    \code
    vector<OBRing*>::iterator i;
    vector<int>::iterator j;
    vector<OBRing*> *rlist = (vector<OBRing*>*)mol.GetData("RingList");
    for (i = rlist->begin();i != rlist->end();++i)
    {
       for(j = (*i)->_path.begin(); j != (*i)->_path.end(); ++j)
          cout << *j << ' ';
       cout << endl;
    }
    \endcode

    will produce something like the following output for benzene:

    \code
    1 2 3 4 5 6
    \endcode

    Ring information is automatically deleted from an OBMol when it goes
    out of scope or the OBMol::Clear() member function is called.

    Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#findSmallestSetOfSmallestRings">blue-obelisk:findSmallestSetOfSmallestRings</a>.
  */

  static int DetermineFRJ(OBMol &);
  static void BuildOBRTreeVector(OBAtom*,OBRTree*,vector<OBRTree*>&,OBBitVec&);

  void OBMol::FindSSSR()
  {
    if (HasSSSRPerceived())
      return;
    SetSSSRPerceived();
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::FindSSSR", obAuditMsg);

    // Delete any old data before we start finding new rings
    // The following procedure is slow
    // So if client code is multi-threaded, we don't want to make them wait
    if (HasData("SSSR")) {
      DeleteData("SSSR");
    }

    OBRing *ring;
    vector<OBRing*>::iterator j;

    //get frerejaque taking int account multiple possible spanning graphs
    int frj = DetermineFRJ(*this);
    if (frj)
      {
        vector<OBRing*> vr;
        FindRingAtomsAndBonds();

        OBBond *bond;
        vector<OBBond*> cbonds;
        vector<OBBond*>::iterator k;

        //restrict search for rings around closure bonds
        for (bond = BeginBond(k);bond;bond = NextBond(k))
          if (bond->IsClosure())
            cbonds.push_back(bond);

        if (!cbonds.empty())
          {
            OBRingSearch rs;
            //search for all rings about closures
            vector<OBBond*>::iterator i;

            for (i = cbonds.begin();i != cbonds.end();++i)
              rs.AddRingFromClosure(*this,(OBBond*)*i);

            rs.SortRings();
            rs.RemoveRedundant(frj);
            //store the SSSR set

            for (j = rs.BeginRings();j != rs.EndRings();++j)
              {
                ring = new OBRing ((*j)->_path,NumAtoms()+1);
                ring->SetParent(this);
                vr.push_back(ring);
              }
            //rs.WriteRings();
          }

        OBRingData *rd = new OBRingData();
        rd->SetOrigin(perceived); // to separate from user or file input
        rd->SetAttribute("SSSR");
        rd->SetData(vr);
        SetData(rd);
      }
  }

  std::vector<unsigned int> atomRingToBondRing(OBMol *mol, const std::vector<int> &atoms)
  {
    std::vector<unsigned int> bonds;
    for (unsigned int i = 0; i < atoms.size() - 1; ++i) {
      unsigned int beginIndex = atoms[i];
      unsigned int endIndex = atoms[i+1];
      unsigned int index = mol->GetBond(beginIndex, endIndex)->GetIdx();
      bonds.push_back(index);
    }
    bonds.push_back(mol->GetBond(atoms[0], atoms[atoms.size()-1])->GetIdx());
    return bonds;
  }

  /**
   * This function finds the LSSR containing all relevant cycles. A cycle is
   * relevant if it belongs to at least one minimum cycle basis. Another
   * description is more useful though:
   *
   * A cycle (C) is relevant if:
   * - no smaller cycles C_i, ..., C_k exist such that C = C_1 + ... + C_k
   * - both bonds & atoms are checked
   *
   * This is based on lemma 1 from:
   *
   * P. Vismara, Union of all the minimum cycle bases of a graph, The electronic
   * journal of combinatorics, Vol. 4, 1997
   * http://www.emis.de/journals/EJC/Volume_4/PostScriptfiles/v4i1r9.ps
   */
  void visitRing(OBMol *mol, OBRing *ring, std::vector<OBRing*> &rlist, std::vector<OBRing*> &rignored)
  {
    OBBitVec mask;
    // Make sure mask is the same size as the maximum ring atom/bond index.
    mask.SetBitOn(mol->NumAtoms());
    mask.SetBitOn(mol->NumBonds());

    //
    // Remove larger rings that cover the same atoms as smaller rings.
    //
    mask.Clear();
    for (unsigned int j = 0; j < rlist.size(); ++j)
      // Here we select only smaller rings.
      if (rlist[j]->_path.size() < ring->_path.size())
        mask |= rlist[j]->_pathset;

    mask = mask & ring->_pathset;

    bool containsSmallerAtomRing = (mask == ring->_pathset) ? true : false;

    // Translate ring atom indexes to ring bond indexes.
    std::vector<unsigned int> bonds = atomRingToBondRing(mol, ring->_path);
    OBBitVec bondset;
    for (unsigned int i = 0; i < bonds.size(); ++i)
      bondset.SetBitOn(bonds[i]);

    //
    // Remove larger rings that cover the same bonds as smaller rings.
    //
    mask.Clear();
    for (unsigned int j = 0; j < rlist.size(); ++j) {
      std::vector<unsigned int> otherBonds = atomRingToBondRing(mol, rlist[j]->_path);
      OBBitVec bs;
      for (unsigned int i = 0; i < otherBonds.size(); ++i)
        bs.SetBitOn(otherBonds[i]);

      // Here we select only smaller rings.
      if (otherBonds.size() < bonds.size())
        mask |= bs;
    }

    mask = mask & bondset;

    bool containsSmallerBondRing = (mask == bondset) ? true : false;

    // The ring is part of the LSSR if all it's atoms and bonds are not
    // found in smaller rings.
    if (!containsSmallerAtomRing || !containsSmallerBondRing) {
      rlist.push_back(ring);
    } else {
      rignored.push_back(ring);
    }
  }


  void OBMol::FindLSSR()
  {
    if (HasLSSRPerceived())
      return;
    SetLSSRPerceived();
    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::FindLSSR", obAuditMsg);

    // Delete any old data before we start finding new rings
    // The following procedure is slow
    // So if client code is multi-threaded, we don't want to make them wait
    if (HasData("LSSR")) {
      DeleteData("LSSR");
    }

    OBRing *ring;
    vector<OBRing*>::iterator j;

    //get frerejaque taking int account multiple possible spanning graphs
    int frj = DetermineFRJ(*this);
    if (frj)
      {
        vector<OBRing*> vr;
        FindRingAtomsAndBonds();

        OBBond *bond;
        vector<OBBond*> cbonds;
        vector<OBBond*>::iterator k;

        //restrict search for rings around closure bonds
        for (bond = BeginBond(k);bond;bond = NextBond(k))
          if (bond->IsClosure())
            cbonds.push_back(bond);

        if (!cbonds.empty())
          {
            OBRingSearch rs;
            //search for all rings about closures
            vector<OBBond*>::iterator i;

            for (i = cbonds.begin();i != cbonds.end();++i)
              rs.AddRingFromClosure(*this,(OBBond*)*i);

            rs.SortRings();
            rs.RemoveRedundant(-1); // -1 means LSSR

            //store the LSSR set

            for (j = rs.BeginRings();j != rs.EndRings();++j)
              {
                ring = new OBRing ((*j)->_path,NumAtoms()+1);
                ring->SetParent(this);
                vr.push_back(ring);
              }
            //rs.WriteRings();
          }

        OBRingData *rd = new OBRingData();
        rd->SetOrigin(perceived); // to separate from user or file input
        rd->SetAttribute("LSSR");
        rd->SetData(vr);
        SetData(rd);
      }
  }

  static int DetermineFRJ(OBMol &mol)
  {
    vector<vector<int> >::iterator i;
    vector<vector<int> > cfl;
    //find all continuous graphs in the mol area
    mol.ContigFragList(cfl);

    if (cfl.empty())
      return(0);
    if (cfl.size() == 1)
      return(mol.NumBonds() - mol.NumAtoms() + 1);

    //count up the atoms and bonds belonging to each graph
    OBBond *bond;
    vector<OBBond*>::iterator j;
    int numatoms,numbonds,frj=0;
    OBBitVec frag;
    for (i = cfl.begin();i != cfl.end();++i)
      {
        frag.Clear();
        frag.FromVecInt(*i);
        numatoms = (*i).size();
        numbonds=0;
        for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
          if (frag.BitIsOn(bond->GetBeginAtomIdx()) &&
              frag.BitIsOn(bond->GetEndAtomIdx()))
            numbonds++;
        frj += numbonds - numatoms + 1;
      }

    return(frj);
  }

  void OBRingSearch::RemoveRedundant(int frj)
  {
    int i,j;

    //remove identical rings
    for (i = _rlist.size()-1;i > 0;i--)
      for (j = i-1;j >= 0;j--)
        if ((_rlist[i])->_pathset == (_rlist[j])->_pathset)
          {
            delete _rlist[i];
            _rlist.erase(_rlist.begin()+i);
            break;
          }

    if (_rlist.size() == 0)
      return; // nothing to do

    // handle LSSR
    if (frj < 0) {
      OBMol *mol = _rlist[0]->GetParent();
      std::vector<OBRing*> rlist, rignored;
      for (unsigned int i = 0; i < _rlist.size(); ++i) {
        visitRing(mol, _rlist[i], rlist, rignored);
      }
      for (unsigned int i = 0; i < rignored.size(); ++i)
        delete rignored[i];
      _rlist = rlist;
      return;
    }

    // exit if we already have frj rings
    if (_rlist.size() == (unsigned)frj)
      return;

    //make sure tmp is the same size as the rings
    OBBitVec tmp;
    for (j = 0;j < (signed)_rlist.size();++j)
      tmp = (_rlist[j])->_pathset;

    //remove larger rings that cover the same atoms as smaller rings
    for (i = _rlist.size()-1;i >= 0;i--)
      {
        tmp.Clear();
        for (j = 0;j < (signed)_rlist.size();++j)
          if ((_rlist[j])->_path.size() <= (_rlist[i])->_path.size() && i != j)
            tmp |= (_rlist[j])->_pathset;

        tmp = tmp & (_rlist[i])->_pathset;

        if (tmp == (_rlist[i])->_pathset)
          {
            delete _rlist[i];
            _rlist.erase(_rlist.begin()+i);
          }

        if (_rlist.size() == (unsigned)frj)
          break;
      }
  }


  void OBRingSearch::AddRingFromClosure(OBMol &mol,OBBond *cbond)
  {
    vector<OBRTree*> t1(mol.NumAtoms()+1,(OBRTree*)NULL);
    vector<OBRTree*> t2(mol.NumAtoms()+1,(OBRTree*)NULL);
    OBBitVec bv1,bv2;

    bv1.SetBitOn(cbond->GetEndAtomIdx());
    bv2.SetBitOn(cbond->GetBeginAtomIdx());
    BuildOBRTreeVector(cbond->GetBeginAtom(),NULL,t1,bv1);
    BuildOBRTreeVector(cbond->GetEndAtom(),NULL,t2,bv2);

    bool pathok;
    deque<int> p1,p2;
    vector<OBAtom*> path1,path2;
    vector<OBAtom*>::iterator m,n;
    vector<OBRTree*>::iterator i;

    for (i = t1.begin();i != t1.end();++i)
      if (*i)
        {
          path1.clear();
          (*i)->PathToRoot(path1);

          if (t2[(*i)->GetAtomIdx()])
            {
              pathok = true;
              path2.clear();
              t2[(*i)->GetAtomIdx()]->PathToRoot(path2);

              p1.clear();
              m = path1.begin();
              if (m != path1.end())
                p1.push_back((*m)->GetIdx());
              for (m = path1.begin(),++m;m != path1.end();++m)
                {
                  p1.push_back((*m)->GetIdx());
                  p2.clear();
                  for (n = path2.begin(),++n;n != path2.end();++n)
                    {
                      p2.push_front((*n)->GetIdx());
                      if (*n == *m)//don't traverse across identical atoms
                        {
                          p2.pop_front();
                          if (p1.size()+p2.size() > 2)
                            SaveUniqueRing(p1,p2);
                          pathok = false;
                          break;
                        }
                      if ((*n)->IsConnected(*m) && p1.size()+p2.size() > 2)
                        SaveUniqueRing(p1,p2);
                    }
                  if (!pathok)
                    break;
                }
            }
        }

    //clean up OBRTree vectors
    for (i = t1.begin();i != t1.end();++i)
      if (*i)
        delete *i;

    for (i = t2.begin();i != t2.end();++i)
      if (*i)
        delete *i;

    // set parent for all rings
    for (unsigned int j = 0; j < _rlist.size(); ++j)
      _rlist[j]->SetParent(&mol);
  }

  bool OBRingSearch::SaveUniqueRing(deque<int> &d1,deque<int> &d2)
  {
    vector<int> path;
    OBBitVec bv;
    deque<int>::iterator i;

    for (i = d1.begin();i != d1.end();++i)
      {
        bv.SetBitOn(*i);
        path.push_back(*i);
      }

    for (i = d2.begin();i != d2.end();++i)
      {
        bv.SetBitOn(*i);
        path.push_back(*i);
      }

    vector<OBRing*>::iterator j;
    for (j = _rlist.begin();j != _rlist.end();++j)
      if (bv == (*j)->_pathset)
        return(false);

    OBRing *ring = new OBRing(path, bv);
    _rlist.push_back(ring);

    return(true);
  }

  //! Destructor -- free all rings created from this search
  OBRingSearch::~OBRingSearch()
  {
    vector<OBRing*>::iterator i;
    for (i = _rlist.begin();i != _rlist.end();++i)
      delete *i;
  }

  bool CompareRingSize(const OBRing *a,const OBRing *b)
  {
    return(a->Size() == b->Size() ? a->ring_id < b->ring_id : a->Size() < b->Size()); // ensure stable sort
  }

  void OBRingSearch::WriteRings()
  {
    vector<OBRing*>::iterator i;

    for (i = _rlist.begin();i != _rlist.end();++i)
      cout << (*i)->_pathset << endl;
  }

  static void FindRings(OBMol &mol,vector<int> &path,OBBitVec &avisit,
                        OBBitVec &bvisit, int natom,int depth );

  void OBMol::FindRingAtomsAndBonds()
  {
    if (HasFlag(OB_RINGFLAGS_MOL))
      return;
    SetFlag(OB_RINGFLAGS_MOL);

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::FindRingAtomsAndBonds", obAuditMsg);

    OBBitVec avisit,bvisit;
    avisit.Resize(NumAtoms()+1);
    bvisit.Resize(NumAtoms()+1);
    vector<int> path;
    path.resize(NumAtoms()+1);

    for(unsigned int i=1; i<= NumAtoms(); i++ )
      if(!avisit[i])
        FindRings(*this,path,avisit,bvisit,i,0);
  }

  static void FindRings(OBMol &mol,vector<int> &path,OBBitVec &avisit,
                        OBBitVec &bvisit, int natom,int depth )
  {
    OBAtom *atom;
    OBBond *bond;
    vector<OBBond*>::iterator k;

    // don't return if all atoms are visited
    // (For example, some atoms are in multiple rings!) -GRH

    if (avisit[natom])
      {
        int j = depth-1;
        bond=mol.GetBond(path[j--]);
        bond->SetInRing();
        while( j >= 0 )
          {
            bond=mol.GetBond(path[j--]);
            bond->SetInRing();
            (bond->GetBeginAtom())->SetInRing();
            (bond->GetEndAtom())->SetInRing();
            if(bond->GetBeginAtomIdx()==static_cast<unsigned int>(natom) || bond->
               GetEndAtomIdx()==static_cast<unsigned int>(natom))
              break;
          }
      }
    else
      {
        avisit.SetBitOn(natom);
        atom = mol.GetAtom(natom);
        for(bond = atom->BeginBond(k);bond;bond=atom->NextBond(k))
          if( !bvisit[bond->GetIdx()])
            {
              path[depth] = bond->GetIdx();
              bvisit.SetBitOn(bond->GetIdx());
              FindRings(mol,path,avisit,bvisit,bond->GetNbrAtomIdx(atom),
                        depth+1);
            }
      }
  }

  bool OBRing::IsAromatic()
  {
    OBMol *mol = _parent;
    vector<int>::iterator i;
    for (i = _path.begin();i != _path.end();++i)
      if (!(mol->GetAtom(*i))->IsAromatic())
        return(false);

    return(true);
  }

  void OBRing::SetType(char *type)
  {
    strncpy(_type,type, sizeof(_type) - 1);
    _type[sizeof(_type) - 1] = '\0';
  }

  void OBRing::SetType(std::string &type)
  {
     strncpy(_type,type.c_str(), sizeof(_type) - 1);
    _type[sizeof(_type) - 1] = '\0';
  }

  char* OBRing::GetType()
  {
    OBMol *mol = (OBMol*)GetParent();
    if (mol && !mol->HasRingTypesPerceived())
      ringtyper.AssignTypes(*((OBMol*)GetParent()));

    return(_type);
  }

  unsigned int OBRing::GetRootAtom()
  {
    vector<int>::iterator i;
    OBMol *mol = (OBMol*)GetParent();

    //if (!IsAromatic())
    //  return 0;

    if (Size() == 6)
      for (i = _path.begin();i != _path.end();++i)
        if (!(mol->GetAtom(*i))->IsCarbon())
	  return (*i);

    if (Size() == 5)
      for (i = _path.begin();i != _path.end();++i) {
        OBAtom *atom = mol->GetAtom(*i);
        if (atom->IsSulfur() && (atom->GetValence() == 2))
	  return (*i);
        if (atom->IsOxygen() && (atom->GetValence() == 2))
	  return (*i);
        if (atom->IsNitrogen() && (atom->BOSum() == atom->GetValence()))
	  return (*i);
      }

    return 0;
  }

  bool OBRing::IsMember(OBAtom *a)
  {
    return(_pathset.BitIsOn(a->GetIdx()));
  }

  bool OBRing::IsMember(OBBond *b)
  {
    return((_pathset.BitIsOn(b->GetBeginAtomIdx()))&&(_pathset.BitIsOn(b->GetEndAtomIdx())));
  }

  OBRing::OBRing(vector<int> &path,int size) : _path(path)
  {
    _pathset.FromVecInt(_path);
    _pathset.Resize(size);
    _type[0] = '\0';
  }

  /*!
  **\brief OBRing copy constructor
  **\param src reference to original OBRing object (rhs)
  */
  OBRing::OBRing(const OBRing &src)
    //no base class
    : _path(src._path),	_pathset(src._pathset)	//chain to member classes
  {
    //member data
    _parent = src._parent;  //this may not be what you want, but it's a default
    _type[0] = '\0';
  }

  /*!
  **\brief OBRing assignment operator
  **\param src reference to original OBRing object (rhs)
  **\return reference to modified OBRing object (lhs)
  */
  OBRing& OBRing::operator =(const OBRing &src)
  {
    //on identity, return
    if(this == &src)
      return(*this);

    //no base class

    //memeber classes & data
    _path = src._path;
    _pathset = src._pathset;
    _parent = src._parent; //note, this may not be what you want
    _type[0] = '\0';

    return(*this);
  }

  void BuildOBRTreeVector(OBAtom *atom,OBRTree *prv,vector<OBRTree*> &vt,OBBitVec &bv)
  {
    vt[atom->GetIdx()] = new OBRTree (atom,prv);

    int i;
    OBAtom *nbr;
    OBMol *mol = (OBMol*)atom->GetParent();
    OBBitVec curr,used,next;
    vector<OBBond*>::iterator j;
    curr |= atom->GetIdx();
    used = bv|curr;

#define OB_RTREE_CUTOFF 20

    int level=0;
    for (;;)
      {
        next.Clear();
        for (i = curr.NextBit(0);i != bv.EndBit();i = curr.NextBit(i))
          {
            atom = mol->GetAtom(i);
            for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
              if (!used[nbr->GetIdx()])
                {
                  next |= nbr->GetIdx();
                  used |= nbr->GetIdx();
                  vt[nbr->GetIdx()] = new OBRTree (nbr,vt[atom->GetIdx()]);
                }
          }

        if (next.Empty())
          break;
        curr = next;
        level++;
        if (level > OB_RTREE_CUTOFF)
          break;
      }
#undef OB_RTREE_CUTOFF
  }

  OBRTree::OBRTree(OBAtom *atom,OBRTree *prv)
  {
    _atom = atom;
    _prv = prv;
  }

  //! The supplied path is built up of OBAtom nodes, with the root atom
  //! the last item in the vector.
  void OBRTree::PathToRoot(vector<OBAtom*> &path)
  {
    path.push_back(_atom);
    if (_prv)
      _prv->PathToRoot(path);
  }

  int OBRTree::GetAtomIdx()
  {
    return(_atom->GetIdx());
  }

  bool OBRing::findCenterAndNormal(vector3 & center, vector3 &norm1, vector3 &norm2)
  {
    OBMol *mol= this->_parent;
    int j= 0;
    const int nA= this->_path.size();
    vector3 tmp;

    center.Set(0.0,0.0,0.0);
    norm1.Set(0.0,0.0,0.0);
    norm2.Set(0.0,0.0,0.0);
    for (j = 0; j != nA; ++j)
      {
        center += (mol->GetAtom(_path[j]))->GetVector();
      }
    center/= double(nA);

    for (j = 0; j != nA; ++j)
      {
        vector3 v1= (mol->GetAtom(_path[j]))->GetVector() - center;
        vector3 v2= (mol->GetAtom(_path[j+1==nA?0:j+1]))->GetVector() - center;
        tmp= cross(v1,v2);
        norm1+= tmp;
      }
    norm1/= double(nA);
    norm1.normalize();
    norm2= norm1;
    norm2 *= -1.0;
    return(true);
  }

} // end namespace OpenBabel

//! \file ring.cpp
//! \brief Deal with rings, find smallest set of smallest rings (SSSR).
