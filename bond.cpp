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

#include "mol.h"
#include "typer.h"

namespace OpenBabel {

extern OEAromaticTyper  aromtyper;

// *******************************
// *** OEBond member functions ***
// *******************************

OEBond::OEBond()
 {
   _idx=0;
   _order=0;
   _flags=0;
   _bgn=NULL;
   _end=NULL;
 }

void OEBond::Set(int idx,OEAtom *begin,OEAtom *end,int order,int flags)
{
  SetIdx(idx);
  SetBegin(begin);
  SetEnd(end);
  SetBO(order);
  SetFlag(flags);
}

void OEBond::SetBO(int order)
{
  _order = (char)order;
  if (order == 5)
    {
      SetAromatic();
      if (_bgn) _bgn->SetAromatic();
      if (_end) _end->SetAromatic();
    }
  else            UnsetAromatic();
}

void OEBond::SetLength(OEAtom *fixed, float length)
{
  unsigned int i;
  OEMol *mol = (OEMol*)fixed->GetParent();
  Vector v1,v2,v3,v4,v5;
  vector<int> children;

  int a = fixed->GetIdx();
  int b = GetNbrAtom(fixed)->GetIdx();

  mol->FindChildren(children,a,b);
  children.push_back(b);

  v1 = GetNbrAtom(fixed)->GetVector();
  v2 = fixed->GetVector();
  v3 = v1 - v2;
  v3.normalize();
  v3 *= length;
  v3 += v2;
  v4 = v3 - v1;

  for ( i = 0 ; i < children.size() ; i++ )
    {
      v1 = mol->GetAtom(children[i])->GetVector();
      v1 += v4;
      mol->GetAtom(children[i])->SetVector(v1);
      /*
      idx = (children[i]-1) * 3;
      c[idx]   += x; 
      c[idx+1] += y;
      c[idx+2] += z;
      */
    }
}
   
bool OEBond::IsRotor()
{
  return(_bgn->GetHvyValence() > 1 && _end->GetHvyValence() > 1 && 
	 _order == 1 && !IsInRing() && _bgn->GetHyb() != 1 &&
	 _end->GetHyb() != 1);
}

bool OEBond::IsAmide()
{
  OEAtom *a1,*a2;
  a1 = a2 = NULL;

  if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
    {a1 = (OEAtom*)_bgn; a2 = (OEAtom*)_end;}

  if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
    {a1 = (OEAtom*)_end; a2 = (OEAtom*)_bgn;}

  if (!a1 || !a2) return(false);
  if (GetBO() != 1) return(false);
  
  OEBond *bond;
  vector<OEEdgeBase*>::iterator i;
  for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
    if (bond->IsCarbonyl())
      return(true);
  
  return(false);
}

bool OEBond::IsPrimaryAmide()
{
  OEAtom *a1,*a2;
  a1 = a2 = NULL;

  if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
    {a1 = (OEAtom*)_bgn; a2 = (OEAtom*)_end;}

  if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
    {a1 = (OEAtom*)_end; a2 = (OEAtom*)_bgn;}

  if (!a1 || !a2) return(false);
  if (GetBO() != 1) return(false);
  
  OEBond *bond;
  vector<OEEdgeBase*>::iterator i;
  for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
    if (bond->IsCarbonyl())
      if (a2->GetHvyValence() == 2)
      return(true);

  
  return(false);
}

bool OEBond::IsSecondaryAmide()
{
  
  
  return(false);
}

bool OEBond::IsEster()
{
  OEAtom *a1,*a2;
  a1 = a2 = NULL;

  if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 8)
    {a1 = (OEAtom*)_bgn; a2 = (OEAtom*)_end;}

  if (_bgn->GetAtomicNum() == 8 && _end->GetAtomicNum() == 6)
    {a1 = (OEAtom*)_end; a2 = (OEAtom*)_bgn;}

  if (!a1 || !a2) return(false);
  if (GetBO() != 1) return(false);
  
  OEBond *bond;
  vector<OEEdgeBase*>::iterator i;
  for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
    if (bond->IsCarbonyl())
      return(true);
  
  
  return(false);
}

bool OEBond::IsCarbonyl()
{
  if (GetBO() != 2) return(false);
  
  if ((_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 8) || 
      (_bgn->GetAtomicNum() == 8 && _end->GetAtomicNum() == 6))
    return(true);

  return(false);
}

bool OEBond::IsSingle()
{
	if (HasFlag(OE_AROMATIC_BOND)) return(false);

	if (!((OEMol*)GetParent())->HasAromaticPerceived())
		{
			aromtyper.AssignAromaticFlags(*((OEMol*)GetParent()));
		}

	if ((this->GetBondOrder()==1) && !(HasFlag(OE_AROMATIC_BOND)))
		return(true);

	return(false);
}

bool OEBond::IsDouble()
{
	if	(HasFlag(OE_AROMATIC_BOND)) return(false);

	if (!((OEMol*)GetParent())->HasAromaticPerceived())
		{
			aromtyper.AssignAromaticFlags(*((OEMol*)GetParent()));
		}
			
	if ((this->GetBondOrder()==2) && !(HasFlag(OE_AROMATIC_BOND)))
		return(true);

	return(false);
}

bool OEBond::IsAromatic() const
{
  if (((OEBond*)this)->HasFlag(OE_AROMATIC_BOND)) return(true);

  OEMol *mol = (OEMol*)((OEBond*)this)->GetParent();
  if (!mol->HasAromaticPerceived())
    {
      aromtyper.AssignAromaticFlags(*mol);
      if (((OEBond*)this)->HasFlag(OE_AROMATIC_BOND)) return(true);
    }

  return(false);
}

void OEBond::SetKSingle()
{
  _flags &= (~(OE_KSINGLE_BOND|OE_KDOUBLE_BOND|OE_KTRIPLE_BOND));
  _flags |= OE_KSINGLE_BOND;
}

void OEBond::SetKDouble()
{
  _flags &= (~(OE_KSINGLE_BOND|OE_KDOUBLE_BOND|OE_KTRIPLE_BOND));
  _flags |= OE_KDOUBLE_BOND;
}

void OEBond::SetKTriple()
{
  _flags &= (~(OE_KSINGLE_BOND|OE_KDOUBLE_BOND|OE_KTRIPLE_BOND));
  _flags |= OE_KTRIPLE_BOND;
}

bool OEBond::IsKSingle()
{
  if (_flags & OE_KSINGLE_BOND) return(true);
  if (!((OEMol*)GetParent())->HasKekulePerceived()) 
	  ((OEMol*)GetParent())->PerceiveKekuleBonds();

  return((_flags & OE_KSINGLE_BOND) != 0) ? true : false;
}

bool OEBond::IsKDouble()
{
  if (_flags & OE_KDOUBLE_BOND) return(true);
  if (!((OEMol*)GetParent())->HasKekulePerceived()) 
	  ((OEMol*)GetParent())->PerceiveKekuleBonds();

  return((_flags & OE_KDOUBLE_BOND) != 0) ? true : false;
}

bool OEBond::IsKTriple()
{
  if (_flags & OE_KTRIPLE_BOND) return(true);
  if (!((OEMol*)GetParent())->HasKekulePerceived()) 
	  ((OEMol*)GetParent())->PerceiveKekuleBonds();

  return((_flags & OE_KTRIPLE_BOND) != 0) ? true : false;
}

bool OEBond::IsInRing() const
{
  if (((OEBond*)this)->HasFlag(OE_RING_BOND)) return(true);
  
  OEMol *mol = (OEMol*)((OEBond*)this)->GetParent();
  if (!mol->HasRingAtomsAndBondsPerceived())
    {
      mol->FindRingAtomsAndBonds();
      if (((OEBond*)this)->HasFlag(OE_RING_BOND)) return(true);
    }

  return(false);
}

bool OEBond::IsClosure()
{
  OEMol *mol = (OEMol*)GetParent();
  if (!mol) return(false);
  if (mol->HasClosureBondsPerceived()) return(HasFlag(OE_CLOSURE_BOND));
  
  mol->SetClosureBondsPerceived();
  
  OEBond *bond;
  OEAtom *atom,*nbr;
  OEBitVec uatoms,ubonds;
  vector<OENodeBase*> curr,next;
  vector<OENodeBase*>::iterator i;
  vector<OEEdgeBase*>::iterator j;

  uatoms.Resize(mol->NumAtoms()+1);
  ubonds.Resize(mol->NumAtoms()+1);

  for (;uatoms.CountBits() < (signed)mol->NumAtoms();)
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
	  for (i = curr.begin();i != curr.end();i++)
	    for (nbr = ((OEAtom*)*i)->BeginNbrAtom(j);nbr;nbr = ((OEAtom*)*i)->NextNbrAtom(j))
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
  
  return(HasFlag(OE_CLOSURE_BOND));
}

float OEBond::GetEquibLength()
{
  float length;
  OEAtom *begin, *end;
  // CorrectedBondRad will always return a # now
  //  if (!CorrectedBondRad(GetBeginAtom(),rad1)) return(0.0);
  //  if (!CorrectedBondRad(GetEndAtom(),rad2))   return(0.0);

  begin = GetBeginAtom();
  end = GetEndAtom();
  length = etab.CorrectedBondRad(begin->GetAtomicNum(), begin->GetHyb())
    + etab.CorrectedBondRad(end->GetAtomicNum(), end->GetHyb());
  
  if (IsAromatic()) length *= 0.93f;
  else if (GetBO() == 2)         length *= 0.91f;
  else if (GetBO() == 3)         length *= 0.87f;
  return(length);
}

float OEBond::GetLength()
{
  float	d2;
  OEAtom *begin, *end;
  begin = GetBeginAtom();
  end = GetEndAtom();
  
  d2 = SQUARE(begin->GetX() - end->GetX());
  d2 += SQUARE(begin->GetY() - end->GetY());
  d2 += SQUARE(begin->GetZ() - end->GetZ());
 
  return(sqrt(d2));
}

}
