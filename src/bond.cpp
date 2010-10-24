/**********************************************************************
bond.cpp - Handle OBBond class

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

#include <openbabel/bond.h>
#include <openbabel/mol.h>
#include <climits>

using namespace std;

namespace OpenBabel
{

  class OBBondPrivate
  {
    public:
      OBBondPrivate() {}
  };

  extern OBAromaticTyper  aromtyper;

  /** \class OBBond bond.h <openbabel/bond.h>
      \brief Bond class

      The OBBond class is straightforward in its data access and
      modification methods. OBBonds store pointers to the atoms on each end
      of the bond. In storing pointers to atoms instead of integer indices,
      the necessity of having to reorder bonds when atoms are shuffled,
      added, or delete is obviated.

      While methods indicate "begin" and "end" atoms in the bond, all methods
      are designed to be independent of atom ordering, with the exception of
      stereochemically aware properties such as IsUp(), IsDown(), IsWedge, or
      IsHash().
  */

  // *******************************
  // *** OBBond member functions ***
  // *******************************

  OBBond::OBBond() /*: d(new OBBondPrivate)*/
  {
    _idx=0;
    _order=0;
    _flags=0;
    _bgn=NULL;
    _end=NULL;
    _vdata.clear();
  }

  OBBond::~OBBond()
  {
    // OBGenericData handled in OBBase parent class.
  }

  /** Mark the main information for a bond
      \param idx The unique bond index for this bond (inside an OBMol)
      \param begin The 'beginning' atom for the bond
      \param end The 'end' atom for the bond
      \param order The bond order (i.e., 1 = single, 2 = double... 5 = aromatic)
      \param flags Any initial property flags
   **/
  void OBBond::Set(int idx,OBAtom *begin,OBAtom *end,int order,int flags)
  {
    SetIdx(idx);
    SetBegin(begin);
    SetEnd(end);
    SetBondOrder(order);
    SetFlag(flags);
  }

  void OBBond::SetBO(int order)
  {
    SetBondOrder(order);
  }

  void OBBond::SetBondOrder(int order)
  {
    _order = (char)order;
    if (order == 5)
      {
        SetAromatic();
        if (_bgn)
          _bgn->SetAromatic();
        if (_end)
          _end->SetAromatic();
      }
    else
      {
        if (order == 1)
          SetKSingle();
        else if (order == 2)
          SetKDouble();
        else if (order == 3)
          SetKTriple();
        else
          UnsetFlag(OB_KSINGLE_BOND | OB_KDOUBLE_BOND | OB_KTRIPLE_BOND);

        UnsetAromatic();
      }
  }

  void OBBond::SetLength(OBAtom *fixed, double length)
  {
    unsigned int i;
    OBMol *mol = (OBMol*)fixed->GetParent();
    vector3 v1,v2,v3,v4,v5;
    vector<int> children;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::SetBondLength", obAuditMsg);

    int a = fixed->GetIdx();
    int b = GetNbrAtom(fixed)->GetIdx();

    if (a == b)
      return; // this would be a problem...

    mol->FindChildren(children,a,b);
    children.push_back(b);

    v1 = GetNbrAtom(fixed)->GetVector();
    v2 = fixed->GetVector();
    v3 = v1 - v2;

    if (IsNearZero(v3.length_2())) { // too small to normalize, move the atoms apart
      obErrorLog.ThrowError(__FUNCTION__,
                            "Atoms are both at the same location, moving out of the way.", obWarning);
      v3.randomUnitVector();
    } else {
      v3.normalize();
    }

    v3 *= length;
    v3 += v2;
    v4 = v3 - v1;

    cerr << "v3: " << v3 << " v4: " << v4 << endl;

    for ( i = 0 ; i < children.size() ; i++ )
      {
        v1 = mol->GetAtom(children[i])->GetVector();
        v1 += v4;
        mol->GetAtom(children[i])->SetVector(v1);
      }
  }

  void OBBond::SetLength(double length)
  {
    OBAtom *atom1 = GetBeginAtom();
    OBAtom *atom2 = GetEndAtom();

    //split the length difference in half, and modify the bond twice
    double firstLength = length + ((GetLength() - length) / 2);

    SetLength(atom1, firstLength);
    SetLength(atom2, length);
  }

  bool OBBond::IsRotor()
  {
    return(_bgn->GetHvyValence() > 1 && _end->GetHvyValence() > 1 &&
           _order == 1 && !IsInRing() && _bgn->GetHyb() != 1 &&
           _end->GetHyb() != 1);
  }

   bool OBBond::IsAmide()
   {
      OBAtom *c,*n;
      c = n = NULL;

      // Look for C-N bond
      if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
      {
         c = (OBAtom*)_bgn;
         n = (OBAtom*)_end;
      }
      if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
      {
         c = (OBAtom*)_end;
         n = (OBAtom*)_bgn;
      }
      if (!c || !n) return(false);
      if (GetBondOrder() != 1) return(false);
      if (n->GetImplicitValence() != 3) return(false);

      // Make sure C is attached to =O
      OBBond *bond;
      vector<OBBond*>::iterator i;
      for (bond = c->BeginBond(i); bond; bond = c->NextBond(i))
      {
         if (bond->IsCarbonyl()) return(true);
      }

      // Return
      return(false);
   }

   bool OBBond::IsPrimaryAmide()
   {
      OBAtom *c,*n;
      c = n = NULL;

      // Look for C-N bond
      if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
      {
         c = (OBAtom*)_bgn;
         n = (OBAtom*)_end;
      }
      if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
      {
         c = (OBAtom*)_end;
         n = (OBAtom*)_bgn;
      }
      if (!c || !n) return(false);
      if (GetBondOrder() != 1) return(false);
      if (n->GetImplicitValence() != 3) return(false);

      // Make sure that N is connected to one non-H
      if (n->GetHvyValence() != 1) return(false);

      // Make sure C is attached to =O
      OBBond *bond;
      vector<OBBond*>::iterator i;
      for (bond = c->BeginBond(i); bond; bond = c->NextBond(i))
      {
         if (bond->IsCarbonyl()) return(true);
      }

      return(false);
   }

   bool OBBond::IsSecondaryAmide()
   {
      OBAtom *c,*n;
      c = n = NULL;

      // Look for C-N bond
      if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
      {
         c = (OBAtom*)_bgn;
         n = (OBAtom*)_end;
      }
      if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
      {
         c = (OBAtom*)_end;
         n = (OBAtom*)_bgn;
      }
      if (!c || !n) return(false);
      if (GetBondOrder() != 1) return(false);
      if (n->GetImplicitValence() != 3) return(false);

      // Make sure that N is connected to two non-H atoms
      if (n->GetHvyValence() != 2) return(false);

      // Make sure C is attached to =O
      OBBond *bond;
      vector<OBBond*>::iterator i;
      for (bond = c->BeginBond(i); bond; bond = c->NextBond(i))
      {
         if (bond->IsCarbonyl()) return(true);
      }

      return(false);
   }

   bool OBBond::IsTertiaryAmide()
   {
      OBAtom *c,*n;
      c = n = NULL;

      // Look for C-N bond
      if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
      {
         c = (OBAtom*)_bgn;
         n = (OBAtom*)_end;
      }
      if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
      {
         c = (OBAtom*)_end;
         n = (OBAtom*)_bgn;
      }
      if (!c || !n) return(false);
      if (GetBondOrder() != 1) return(false);
      if (n->GetImplicitValence() != 3) return(false);

      // Make sure that N is connected to three non-H atoms
      if (n->GetHvyValence() != 3) return(false);

      // Make sure C is attached to =O
      OBBond *bond;
      vector<OBBond*>::iterator i;
      for (bond = c->BeginBond(i); bond; bond = c->NextBond(i))
      {
         if (bond->IsCarbonyl()) return(true);
      }

      return(false);
   }

/*
  bool OBBond::IsAmide()
  {
    OBAtom *a1,*a2;
    a1 = a2 = NULL;

    if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
      {
        a1 = (OBAtom*)_bgn;
        a2 = (OBAtom*)_end;
      }

    if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
      {
        a1 = (OBAtom*)_end;
        a2 = (OBAtom*)_bgn;
      }

    if (!a1 || !a2)
      return(false);
    if (GetBO() != 1)
      return(false);

    OBBond *bond;
    vector<OBBond*>::iterator i;
    for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
      if (bond->IsCarbonyl())
        return(true);

    return(false);
  }

  bool OBBond::IsPrimaryAmide()
  {
    OBAtom *a1,*a2;
    a1 = a2 = NULL;

    if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
      {
        a1 = (OBAtom*)_bgn;
        a2 = (OBAtom*)_end;
      }

    if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
      {
        a1 = (OBAtom*)_end;
        a2 = (OBAtom*)_bgn;
      }

    if (!a1 || !a2)
      return(false);
    if (GetBO() != 1)
      return(false);

    OBBond *bond;
    vector<OBBond*>::iterator i;
    for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
      if (bond->IsCarbonyl())
        if (a2->GetHvyValence() == 2)
          return(true);

    return(false);
  }

  bool OBBond::IsSecondaryAmide()
  {
    OBAtom *a1,*a2;
    a1 = a2 = NULL;

    if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
      {
        a1 = (OBAtom*)_bgn;
        a2 = (OBAtom*)_end;
      }

    if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
      {
        a1 = (OBAtom*)_end;
        a2 = (OBAtom*)_bgn;
      }

    if (!a1 || !a2)
      return(false);
    if (GetBO() != 1)
      return(false);

    OBBond *bond;
    vector<OBBond*>::iterator i;
    for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
      if (bond->IsCarbonyl())
        if (a2->GetHvyValence() == 3)
          return(true);

    return(false);
  }
*/

  bool OBBond::IsEster()
  {
    OBAtom *a1,*a2;
    a1 = a2 = NULL;

    if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 8)
      {
        a1 = (OBAtom*)_bgn;
        a2 = (OBAtom*)_end;
      }

    if (_bgn->GetAtomicNum() == 8 && _end->GetAtomicNum() == 6)
      {
        a1 = (OBAtom*)_end;
        a2 = (OBAtom*)_bgn;
      }

    if (!a1 || !a2)
      return(false);
    if (GetBO() != 1)
      return(false);

    OBBond *bond;
    vector<OBBond*>::iterator i;
    for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
      if (bond->IsCarbonyl())
        return(true);

    return(false);
  }

  bool OBBond::IsCarbonyl()
  {
    if (GetBO() != 2)
      return(false);

    if ((_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 8) ||
        (_bgn->GetAtomicNum() == 8 && _end->GetAtomicNum() == 6))
      return(true);

    return(false);
  }

  bool OBBond::IsSingle()
  {
    if (HasFlag(OB_AROMATIC_BOND))
      return(false);

    if (!((OBMol*)GetParent())->HasAromaticPerceived())
      {
        aromtyper.AssignAromaticFlags(*((OBMol*)GetParent()));
      }

    if ((this->GetBondOrder()==1) && !(HasFlag(OB_AROMATIC_BOND)))
      return(true);

    return(false);
  }

  bool OBBond::IsDouble()
  {
    if	(HasFlag(OB_AROMATIC_BOND))
      return(false);

    if (!((OBMol*)GetParent())->HasAromaticPerceived())
      {
        aromtyper.AssignAromaticFlags(*((OBMol*)GetParent()));
      }

    if ((this->GetBondOrder()==2) && !(HasFlag(OB_AROMATIC_BOND)))
      return(true);

    return(false);
  }

  bool OBBond::IsTriple()
  {
    if	(HasFlag(OB_AROMATIC_BOND))
      return(false);

    if (!((OBMol*)GetParent())->HasAromaticPerceived())
      {
        aromtyper.AssignAromaticFlags(*((OBMol*)GetParent()));
      }

    if ((this->GetBondOrder()==3) && !(HasFlag(OB_AROMATIC_BOND)))
      return(true);

    return(false);
  }

  bool OBBond::IsAromatic() const
  {
    if (((OBBond*)this)->HasFlag(OB_AROMATIC_BOND))
      return(true);

    OBMol *mol = (OBMol*)((OBBond*)this)->GetParent();
    if (!mol->HasAromaticPerceived())
      {
        aromtyper.AssignAromaticFlags(*mol);
        if (((OBBond*)this)->HasFlag(OB_AROMATIC_BOND))
          return(true);
      }

    return(false);
  }

  /*! This method checks if the geometry around this bond looks unsaturated
    by measuring the torsion angles formed by all connected atoms X-start=end-Y
    and checking that they are close to 0 or 180 degrees */
  bool OBBond::IsDoubleBondGeometry()
  {
    double torsion;
    OBAtom *nbrStart,*nbrEnd;
    vector<OBBond*>::iterator i,j;
    // We concentrate on sp2 atoms with valence up to 3 and ignore the rest (like sp1 or S,P)
    // As this is called from PerceiveBondOrders, GetHyb() may still be undefined.
    if (_bgn->GetHyb()==1 || _bgn->GetValence()>3||
        _end->GetHyb()==1 || _end->GetValence()>3)
      return(true);

    for (nbrStart = static_cast<OBAtom*>(_bgn)->BeginNbrAtom(i); nbrStart;
         nbrStart = static_cast<OBAtom*>(_bgn)->NextNbrAtom(i))
      {
        if (nbrStart != _end)
          {
            for (nbrEnd = static_cast<OBAtom*>(_end)->BeginNbrAtom(j);
                 nbrEnd; nbrEnd = static_cast<OBAtom*>(_end)->NextNbrAtom(j))
              {
                if (nbrEnd != _bgn)
                  {
                    torsion=fabs(CalcTorsionAngle(nbrStart->GetVector(),
                                                  static_cast<OBAtom*>(_bgn)->GetVector(),
                                                  static_cast<OBAtom*>(_end)->GetVector(),
                                                  nbrEnd->GetVector()));

                    // >12&&<168 not enough
                    if (torsion > 15.0  && torsion < 160.0)
                      {
                        // Geometry does not match a double bond
                        return(false);
                      }

                  }
              }  // end loop for neighbors of end
          }
      } // end loop for neighbors of start
    return(true);
  }

  void OBBond::SetKSingle()
  {
    _flags &= (~(OB_KSINGLE_BOND|OB_KDOUBLE_BOND|OB_KTRIPLE_BOND));
    _flags |= OB_KSINGLE_BOND;
  }

  void OBBond::SetKDouble()
  {
    _flags &= (~(OB_KSINGLE_BOND|OB_KDOUBLE_BOND|OB_KTRIPLE_BOND));
    _flags |= OB_KDOUBLE_BOND;
  }

  void OBBond::SetKTriple()
  {
    _flags &= (~(OB_KSINGLE_BOND|OB_KDOUBLE_BOND|OB_KTRIPLE_BOND));
    _flags |= OB_KTRIPLE_BOND;
  }

  bool OBBond::IsKSingle()
  {
    if (_flags & OB_KSINGLE_BOND)
      return(true);
    if (!((OBMol*)GetParent())->HasKekulePerceived())
      ((OBMol*)GetParent())->NewPerceiveKekuleBonds();

    return((_flags & OB_KSINGLE_BOND) != 0) ? true : false;
  }

  bool OBBond::IsKDouble()
  {
    if (_flags & OB_KDOUBLE_BOND)
      return(true);
    if (!((OBMol*)GetParent())->HasKekulePerceived())
      ((OBMol*)GetParent())->NewPerceiveKekuleBonds();

    return((_flags & OB_KDOUBLE_BOND) != 0) ? true : false;
  }

  bool OBBond::IsKTriple()
  {
    if (_flags & OB_KTRIPLE_BOND)
      return(true);
    if (!((OBMol*)GetParent())->HasKekulePerceived())
      ((OBMol*)GetParent())->NewPerceiveKekuleBonds();

    return((_flags & OB_KTRIPLE_BOND) != 0) ? true : false;
  }

  bool OBBond::IsInRing() const
  {
    if (((OBBond*)this)->HasFlag(OB_RING_BOND))
      return(true);

    OBMol *mol = (OBMol*)((OBBond*)this)->GetParent();
    if (!mol->HasRingAtomsAndBondsPerceived())
      {
        mol->FindRingAtomsAndBonds();
        if (((OBBond*)this)->HasFlag(OB_RING_BOND))
          return(true);
      }

    return(false);
  }

  // Adapted from OBAtom::IsInRingSize()
  OBRing* OBBond::FindSmallestRing() const
  {
    vector<OBRing*> rlist;
    vector<OBRing*>::iterator i;

    OBMol *mol = (OBMol*)((OBBond*)this)->GetParent();

    rlist = mol->GetSSSR();
    OBRing* result = (OBRing*) NULL;
    size_t min_size = UINT_MAX;
    for (i = rlist.begin();i != rlist.end();++i) {
      if ((*i)->IsMember((OBBond*)this) && (*i)->Size() < min_size) {
        min_size = (*i)->Size();
        result = *i;
      }
    }
    return result;
  }

  bool OBBond::IsClosure()
  {
    OBMol *mol = (OBMol*)GetParent();
    if (!mol)
      return(false);
    if (mol->HasClosureBondsPerceived())
      return(HasFlag(OB_CLOSURE_BOND));

    mol->SetClosureBondsPerceived();

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::PerceiveClosureBonds", obAuditMsg);

    OBBond *bond;
    OBAtom *atom,*nbr;
    OBBitVec uatoms,ubonds;
    vector<OBAtom*> curr,next;
    vector<OBAtom*>::iterator i;
    vector<OBBond*>::iterator j;

    uatoms.Resize(mol->NumAtoms()+1);
    ubonds.Resize(mol->NumAtoms()+1);

    for (;static_cast<unsigned int>(uatoms.CountBits()) < mol->NumAtoms();)
      {
        if (curr.empty())
          for (atom = mol->BeginAtom(i);atom;atom = mol->NextAtom(i))
            if (!uatoms[atom->GetIdx()])
              {
                uatoms |= atom->GetIdx();
                curr.push_back(atom);
                break;
              }

        for (;!curr.empty();)
          {
            for (i = curr.begin();i != curr.end();++i)
              for (nbr = ((OBAtom*)*i)->BeginNbrAtom(j);nbr;nbr = ((OBAtom*)*i)->NextNbrAtom(j))
                if (!uatoms[nbr->GetIdx()])
                  {
                    uatoms |= nbr->GetIdx();
                    ubonds |= (*j)->GetIdx();
                    next.push_back(nbr);
                  }

            curr = next;
            next.clear();
          }
      }

    for (bond = mol->BeginBond(j);bond;bond = mol->NextBond(j))
      if (!ubonds[bond->GetIdx()])
        bond->SetClosure();

    return(HasFlag(OB_CLOSURE_BOND));
  }

  double OBBond::GetEquibLength() const
  {
    double length;
    const OBAtom *begin, *end;
    // CorrectedBondRad will always return a # now
    //  if (!CorrectedBondRad(GetBeginAtom(),rad1)) return(0.0);
    //  if (!CorrectedBondRad(GetEndAtom(),rad2))   return(0.0);

    begin = GetBeginAtom();
    end = GetEndAtom();
    length = etab.CorrectedBondRad(begin->GetAtomicNum(), begin->GetHyb())
      + etab.CorrectedBondRad(end->GetAtomicNum(), end->GetHyb());

    if (IsAromatic())
      length *= 0.93;
    else if (GetBO() == 2)
      length *= 0.91;
    else if (GetBO() == 3)
      length *= 0.87;
    return(length);
  }

  double OBBond::GetLength() const
  {
    double	d2;
    const OBAtom *begin, *end;
    begin = GetBeginAtom();
    end = GetEndAtom();

    d2 = SQUARE(begin->GetX() - end->GetX());
    d2 += SQUARE(begin->GetY() - end->GetY());
    d2 += SQUARE(begin->GetZ() - end->GetZ());

    return(sqrt(d2));
  }

  /*Now in OBBase
  // OBGenericData methods
  bool OBBond::HasData(string &s)
  //returns true if the generic attribute/value pair exists
  {
  if (_vdata.empty())
  return(false);

  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetAttribute() == s)
  return(true);

  return(false);
  }

  bool OBBond::HasData(const char *s)
  //returns true if the generic attribute/value pair exists
  {
  if (_vdata.empty())
  return(false);

  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetAttribute() == s)
  return(true);

  return(false);
  }

  bool OBBond::HasData(unsigned int dt)
  //returns true if the generic attribute/value pair exists
  {
  if (_vdata.empty())
  return(false);

  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetDataType() == dt)
  return(true);

  return(false);
  }

  OBGenericData *OBBond::GetData(string &s)
  //returns the value given an attribute
  {
  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetAttribute() == s)
  return(*i);

  return(NULL);
  }

  OBGenericData *OBBond::GetData(const char *s)
  //returns the value given an attribute
  {
  vector<OBGenericData*>::iterator i;

  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetAttribute() == s)
  return(*i);

  return(NULL);
  }

  OBGenericData *OBBond::GetData(unsigned int dt)
  {
  vector<OBGenericData*>::iterator i;
  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetDataType() == dt)
  return(*i);
  return(NULL);
  }

  void OBBond::DeleteData(unsigned int dt)
  {
  vector<OBGenericData*> vdata;
  vector<OBGenericData*>::iterator i;
  for (i = _vdata.begin();i != _vdata.end();++i)
  if ((*i)->GetDataType() == dt)
  delete *i;
  else
  vdata.push_back(*i);
  _vdata = vdata;
  }

  void OBBond::DeleteData(vector<OBGenericData*> &vg)
  {
  vector<OBGenericData*> vdata;
  vector<OBGenericData*>::iterator i,j;

  bool del;
  for (i = _vdata.begin();i != _vdata.end();++i)
  {
  del = false;
  for (j = vg.begin();j != vg.end();++j)
  if (*i == *j)
  {
  del = true;
  break;
  }
  if (del)
  delete *i;
  else
  vdata.push_back(*i);
  }
  _vdata = vdata;
  }

  void OBBond::DeleteData(OBGenericData *gd)
  {
  vector<OBGenericData*>::iterator i;
  for (i = _vdata.begin();i != _vdata.end();++i)
  if (*i == gd)
  {
  delete *i;
  _vdata.erase(i);
  }

  }
  */

} // end namespace OpenBabel

//! \file bond.cpp
//! \brief Handle OBBond class
