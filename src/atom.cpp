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
#include "molchrg.h"
#include "phmodel.h"

#include "math/matrix3x3.h"

using namespace std;

namespace OpenBabel {

extern OBAromaticTyper  aromtyper;
extern OBAtomTyper      atomtyper;
extern OBPhModel        phmodel;

//
// OBAtom member functions
//

OBAtom::OBAtom()
{
  _parent = (OBMol*)NULL;
  Clear();
}

OBAtom::~OBAtom()
{
    if (_residue != NULL)
        _residue->RemoveAtom(this);
}

void OBAtom::Clear()
{
  _c = (float**)NULL;
  _cidx = 0;
  _flags=0;
  _idx = 0;
  _hyb = 0;
  _ele = (char)0;
  _impval = 0;
  _fcharge = 0;
  _type[0] = '\0';
  _pcharge = 0.0;
  _vbond.clear();
  _vbond.reserve(4);
  _residue = (OBResidue*)NULL;
}

OBAtom &OBAtom::operator=(OBAtom &src)
     //copy atom information 
     //bond info is not copied here as ptrs may be invalid
{
  _idx = src.GetIdx();
  _hyb = src.GetHyb();
  _ele = src.GetAtomicNum();
  _fcharge = src.GetFormalCharge();
  strcpy(_type,src.GetType());
  _pcharge = src.GetPartialCharge();
  _v = src.GetVector();
  _flags = src.GetFlag();
  _residue = (OBResidue*)NULL;
  return(*this);
}

bool OBAtom::IsConnected(OBAtom *a1)
{
  vector<OBEdgeBase*>::iterator i;
  OBBond *bond;

  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if (bond->GetBeginAtom() == a1 || bond->GetEndAtom() == a1)
      return(true);

  return(false);
}

bool OBAtom::IsOneThree(OBAtom *a1)
{
  OBAtom *atom1,*atom2;
  OBBond *bond1,*bond2;
  vector<OBEdgeBase*>::iterator i,j;
  atom1 = this;
  atom2 = a1;

  for (bond1 = atom1->BeginBond(i);bond1;bond1 = atom1->NextBond(i))
    for (bond2 = atom2->BeginBond(j);bond2;bond2 = atom2->NextBond(j))
      if (bond1->GetNbrAtom(atom1) == bond2->GetNbrAtom(atom2))
	return(true);

  return(false);
}

bool OBAtom::IsOneFour(OBAtom *a1)
{
  OBAtom *atom1,*atom2;
  OBBond *bond1,*bond2;
  vector<OBEdgeBase*>::iterator i,j;
  atom1 = this;
  atom2 = a1;

  for (bond1 = atom1->BeginBond(i);bond1;bond1 = atom1->NextBond(i))
    for (bond2 = atom2->BeginBond(j);bond2;bond2 = atom2->NextBond(j))
      if ((bond1->GetNbrAtom(atom1))->IsConnected(bond2->GetNbrAtom(atom2)))
	return(true);

  return(false);
}

bool OBAtom::IsAxial()
{
  float tor;
  OBAtom *a,*b,*c;
  vector<OBEdgeBase*>::iterator i,j,k;
  
  for (a = BeginNbrAtom(i);a;a = NextNbrAtom(i))
    if (a->GetHyb() == 3 && a->IsInRing() && !(*i)->IsInRing())
      for (b = a->BeginNbrAtom(j);b;b = a->NextNbrAtom(j))
	if (b != this && b->IsInRing() && b->GetHyb() == 3)
	  for (c = b->BeginNbrAtom(k);c;c = b->NextNbrAtom(k))
	    if (c != a && c->IsInRing())
	      {
			tor = fabs(((OBMol*)GetParent())->GetTorsion(this,a,b,c));
			return(tor > 55.0 && tor < 75.0);
	      }

  return(false);
}


bool OBAtom::HasAlphaBetaUnsat(bool includePandS)
{
  OBAtom *a1,*a2;
  vector<OBEdgeBase*>::iterator i,j;

  for (a1 = BeginNbrAtom(i);a1;a1 = NextNbrAtom(i))
    if (includePandS || (!a1->IsPhosphorus() && !a1->IsSulfur()))
    for (a2 = a1->BeginNbrAtom(j);a2;a2 = a1->NextNbrAtom(j))
      if (a2 != this && ((*j)->GetBO() == 2 || (*j)->GetBO() == 3 || (*j)->GetBO() == 5))
	return(true);

  return(false);
}

bool OBAtom::HasBondOfOrder(unsigned int order)
{
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if (bond->GetBO() == order)
      return(true);

  return(false);
}

int OBAtom::CountBondsOfOrder(unsigned int order)
{
	int count = 0;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if (bond->GetBO() == order)
      count++;

  return(count);
}

bool OBAtom::HasNonSingleBond()
{
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if (bond->GetBO() != 1)
      return(true);

  return(false);
}

bool OBAtom::IsPolarHydrogen()
{
  if (!IsHydrogen()) return(false);
  
  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    {
      atom = bond->GetNbrAtom(this);
      if (atom->GetAtomicNum() == 7) return(true);
      if (atom->GetAtomicNum() == 8) return(true);
      if (atom->GetAtomicNum() == 15) return(true);
      if (atom->GetAtomicNum() == 16) return(true);
    }

  return(false);
}

bool OBAtom::IsNonPolarHydrogen()
{
  if (!IsHydrogen()) return(false);
  
  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    {
      atom = bond->GetNbrAtom(this);
      if (atom->GetAtomicNum() == 6) return(true);
    }

  return(false);
}

vector3 &OBAtom::GetVector()
{
  if (!_c) return(_v);

  _v.Set((*_c)[_cidx],(*_c)[_cidx+1],(*_c)[_cidx+2]);
  return(_v);
}

void OBAtom::SetVector()
{
  obAssert(_c);
  if (_c) _v.Set((*_c)[_cidx],(*_c)[_cidx+1],(*_c)[_cidx+2]);
}

void OBAtom::SetVector(vector3 &v)
{
  if (!_c) _v = v;
  else
    {
      (*_c)[_cidx  ] = v.x();
      (*_c)[_cidx+1] = v.y();
      (*_c)[_cidx+2] = v.z();
    }
}

void OBAtom::SetVector(const float x,const float y,const float z)
{
  if (!_c) _v.Set(x,y,z);
  else
    {
      (*_c)[_cidx  ] = x;
      (*_c)[_cidx+1] = y;
      (*_c)[_cidx+2] = z;
    }
}

OBAtom *OBAtom::GetNextAtom()
{
  OBMol *mol = (OBMol*)GetParent();
  return(((unsigned)GetIdx() == mol->NumAtoms())? NULL : mol->GetAtom(GetIdx()+1));
}

OBResidue *OBAtom::GetResidue()
{
    if (_residue != NULL)
        return _residue;
    else if (!((OBMol*)GetParent())->HasChainsPerceived())
    {
        chainsparser.PerceiveChains(*((OBMol*)GetParent()));
        return _residue;
    }
    else
        return NULL;
}

char *OBAtom::GetType()
{
  OBMol *mol = (OBMol*)GetParent();
  if (mol && !mol->HasAtomTypesPerceived())
    atomtyper.AssignTypes(*((OBMol*)GetParent()));

  if (strlen(_type) == 0) // Somehow we still don't have a type!
    {
      char num[6];
      OBTypeTable tempTable;
      tempTable.SetFromType("ATN");
      tempTable.SetToType("INT");
      snprintf(num, 6, "%d", GetAtomicNum());
      tempTable.Translate(_type, num);
    }

  return(_type);
}

unsigned int OBAtom::GetImplicitValence() const
{
  OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
  if (mol && !mol->HasImplicitValencePerceived())
    atomtyper.AssignImplicitValence(*((OBMol*)((OBAtom*)this)->GetParent()));

  return((unsigned int)_impval);
}

unsigned int OBAtom::GetHyb() const
{
  //hybridization is assigned when atoms are typed
  OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
  if (mol && !mol->HasHybridizationPerceived())
    atomtyper.AssignHyb(*mol);

  return(_hyb);
}


unsigned int OBAtom::GetHvyValence() const
     //returns the number of non-hydrogens connected to an atom
{
  unsigned int count=0;

  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = ((OBAtom*)this)->BeginBond(i); bond; bond = ((OBAtom*)this)->NextBond(i))
    if (!(bond->GetNbrAtom((OBAtom*)this)->IsHydrogen()))
	count++;

  return(count);
}

unsigned int OBAtom::GetHeteroValence() const
     //returns the number of heteroatoms connected to an atom
{
  unsigned int count=0;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  for (bond = ((OBAtom*)this)->BeginBond(i);bond;bond = ((OBAtom*)this)->NextBond(i))
    if (bond->GetNbrAtom((OBAtom*)this)->IsHeteroatom())
      count++;

  return((unsigned int)count);
}

float OBAtom::GetPartialCharge()
{
  if (!GetParent()) return(_pcharge);
  if (!((OBMol*)GetParent())->AutomaticPartialCharge()) return(_pcharge);

  if (!((OBMol*)GetParent())->HasPartialChargesPerceived())
    {
      //seed partial charges are set in the atom typing procedure
      OBAtom *atom;
      OBMol *mol = (OBMol*)GetParent();
      vector<OBNodeBase*>::iterator i;
      for (atom = mol->BeginAtom(i); atom; atom = mol->NextAtom(i))
	atom->SetPartialCharge(0.0f);

      phmodel.AssignSeedPartialCharge(*((OBMol*)GetParent()));
      OBGastChrg gc;
      gc.AssignPartialCharges(*((OBMol*)GetParent()));
    }

  return(_pcharge);
}

bool OBAtom::IsAmideNitrogen()
     //returns true if nitrogen is part of an amide
{
  if (!IsNitrogen()) return(false);

  OBAtom *nbratom,*atom;
  OBBond *abbond,*bond;

  vector<OBEdgeBase*>::iterator i,j;
  atom = this;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    {
      nbratom = bond->GetNbrAtom(atom);
      for (abbond = nbratom->BeginBond(j);abbond;abbond = nbratom->NextBond(j))
	if (abbond->GetBO() == 2 && 
	   (((abbond->GetNbrAtom(nbratom))->GetAtomicNum() == 8) ||
			((abbond->GetNbrAtom(nbratom))->GetAtomicNum() == 16)))
	  return(true);
    }  

  return(false);
}

bool OBAtom::IsAromaticNOxide()
{
	if (!IsNitrogen() || !IsAromatic()) return(false);

	OBAtom *atom;
	vector<OBEdgeBase*>::iterator i;

	for (atom = BeginNbrAtom(i);atom;atom = NextNbrAtom(i))
		if (atom->IsOxygen() && !(*i)->IsInRing() && (*i)->GetBO() == 2)		
			return(true);

	return(false);
}

bool OBAtom::IsCarboxylOxygen()
{
  if (!IsOxygen()) return(false);
  if (GetHvyValence() != 1) return(false);

  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;

  atom = NULL;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if ((bond->GetNbrAtom(this))->IsCarbon())
      {
	atom = bond->GetNbrAtom(this);
	break;
      }
  if (!atom) return(false);
  if (atom->CountFreeOxygens() != 2) return(false);

  //atom is connected to a carbon that has a total 
  //of 2 attached free oxygens
  return(true);
}

bool OBAtom::IsPhosphateOxygen()
{
  if (!IsOxygen()) return(false);
  if (GetHvyValence() != 1) return(false);

  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;

  atom = NULL;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if ((bond->GetNbrAtom(this))->IsPhosphorus())
      {
	atom = bond->GetNbrAtom(this);
	break;
      }
  if (!atom) return(false);
  if (atom->CountFreeOxygens() > 2) return(true);

  //atom is connected to a carbon that has a total 
  //of 2 attached free oxygens
  return(false);
}

bool OBAtom::IsSulfateOxygen()
{
  if (!IsOxygen()) return(false);
  if (GetHvyValence() != 1) return(false);

  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;

  atom = NULL;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if ((bond->GetNbrAtom(this))->IsSulfur())
      {
	atom = bond->GetNbrAtom(this);
	break;
      }
  if (!atom) return(false);
  if (atom->CountFreeOxygens() < 3) return(false);

  //atom is connected to a carbon that has a total 
  //of 2 attached free oxygens
  return(true);
}

bool OBAtom::IsNitroOxygen()
{
  if (!IsOxygen()) return(false);
  if (GetHvyValence() != 1) return(false);

  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;

  atom = NULL;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if ((bond->GetNbrAtom(this))->IsNitrogen())
      {
	atom = bond->GetNbrAtom(this);
	break;
      }
  if (!atom) return(false);
  if (atom->CountFreeOxygens() != 2) return(false);

  //atom is connected to a nitrogen that has a total 
  //of 2 attached free oxygens
  return(true);
}

bool OBAtom::IsHeteroatom()
{
  switch(GetAtomicNum())
    {
    case 7:
    case 8:
    case 9:
    case 15:
    case 16:
    case 17:
    case 35:
    case 53:
      return(true);
    }
  return(false);
}

bool OBAtom::IsAromatic() const
{
  if (((OBAtom*)this)->HasFlag(OB_AROMATIC_ATOM)) return(true);

  OBMol	*mol = (OBMol*)((OBAtom*)this)->GetParent();

  if (!mol->HasAromaticPerceived())
    {
      aromtyper.AssignAromaticFlags(*mol);
      if (((OBAtom*)this)->HasFlag(OB_AROMATIC_ATOM)) return(true);
    }

  return(false);
}

bool OBAtom::IsInRing() const
{
  if (((OBAtom*)this)->HasFlag(OB_RING_ATOM)) return(true);

  OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
  if (!mol->HasRingAtomsAndBondsPerceived())
    {
      mol->FindRingAtomsAndBonds();
      if (((OBAtom*)this)->HasFlag(OB_RING_ATOM)) return(true);
    }

  return(false);
}

bool OBAtom::IsChiral()
{
  if (HasFlag(OB_CHIRAL_ATOM)) return(true);
  
  if (!((OBMol*)GetParent())->HasChiralityPerceived())
    {
      ((OBMol*)GetParent())->FindChiralCenters();
      if (HasFlag(OB_CHIRAL_ATOM)) return(true);
    }

  return(false);
}

bool OBAtom::IsInRingSize(int size) const
{
  vector<OBRing*> rlist;
  vector<OBRing*>::iterator i;

  OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
  if (!mol->HasSSSRPerceived())
    mol->FindSSSR();

 if (!((OBAtom*)this)->HasFlag(OB_RING_ATOM)) return(false);

  rlist = mol->GetSSSR();
  for (i = rlist.begin();i != rlist.end();i++)
    if ((*i)->IsInRing(GetIdx()) && (*i)->PathSize() == size)
      return(true);
  
  return(false);
}

unsigned int OBAtom::MemberOfRingCount() const
{
  vector<OBRing*> rlist;
  vector<OBRing*>::iterator i;
  unsigned int count=0;

  OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();

  if (!mol->HasSSSRPerceived())
    mol->FindSSSR();

  if (!((OBAtom*)this)->IsInRing()) return(0);

  rlist = mol->GetSSSR();

  for (i = rlist.begin();i != rlist.end();i++)
    if ((*i)->IsInRing(GetIdx()))
      count++;
  
  return((unsigned int)count);
}

unsigned int OBAtom::MemberOfRingSize() const
{
  vector<OBRing*> rlist;
  vector<OBRing*>::iterator i;

  OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();

  if (!mol->HasSSSRPerceived())
    mol->FindSSSR();

  if (!((OBAtom*)this)->IsInRing()) return(0);

  rlist = mol->GetSSSR();

  for (i = rlist.begin();i != rlist.end();i++)
    if ((*i)->IsInRing(GetIdx()))
      return((*i)->Size());
  
  return(0);
}

unsigned int OBAtom::CountFreeOxygens() const
{
  unsigned int count = 0;
  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;

  for (bond = ((OBAtom*)this)->BeginBond(i);bond;bond = ((OBAtom*)this)->NextBond(i))
    {
      atom = bond->GetNbrAtom((OBAtom*)this);
      if (atom->IsOxygen() && atom->GetHvyValence() == 1)
	count++;
    }

  return(count);
}

unsigned int OBAtom::BOSum() const
{
  unsigned int bo;
  unsigned int bosum=0;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  
  for (bond = ((OBAtom*)this)->BeginBond(i);bond;bond = ((OBAtom*)this)->NextBond(i))
    {
      bo = bond->GetBO();
      bosum += (bo < 4) ? 2*bo : 3;
    }

  bosum /= 2;
  return(bosum);
}

unsigned int OBAtom::KBOSum() const
{
  OBBond *bond;
  unsigned int bosum = 0;
  vector<OBEdgeBase*>::iterator i;
  
  bosum = GetImplicitValence();

  for (bond = ((OBAtom*)this)->BeginBond(i);bond;bond = ((OBAtom*)this)->NextBond(i))
    {
      if (bond->IsKDouble()) bosum++;
      else if (bond->IsKTriple()) bosum += 2;
    }

  return(bosum);
}

unsigned int OBAtom::ImplicitHydrogenCount() const
     //handles H,C,N,S,O,X
{
  OBMol *mol = (OBMol*)((OBAtom*)this)->GetParent();
  if (mol && !mol->HasImplicitValencePerceived())
    atomtyper.AssignImplicitValence(*((OBMol*)((OBAtom*)this)->GetParent()));

  int impval = _impval - GetHvyValence();
  return((impval>0)?impval:0);
}

unsigned int OBAtom::ExplicitHydrogenCount() const
{
  int numH=0;
  OBAtom *atom;
  vector<OBEdgeBase*>::iterator i;
  for (atom = ((OBAtom*)this)->BeginNbrAtom(i);atom;atom = ((OBAtom*)this)->NextNbrAtom(i))
    if (atom->IsHydrogen())
      numH++;

  return(numH);
}

bool OBAtom::DeleteBond(OBBond *bond)
{
  vector<OBEdgeBase*>::iterator i;
  for (i = _vbond.begin();i != _vbond.end();i++)
    if ((OBBond*)bond == *i)
    {
      _vbond.erase(i);
      return(true);
    }
  return(false);
}

OBBond *OBAtom::BeginBond(vector<OBEdgeBase*>::iterator &i) 
{
	i = _vbond.begin();
	return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
}

OBBond *OBAtom::NextBond(vector<OBEdgeBase*>::iterator &i) 
{
	i++;
	return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
}

OBAtom *OBAtom::BeginNbrAtom(vector<OBEdgeBase*>::iterator &i)
{
	i = _vbond.begin();
  return((i != _vbond.end()) ? ((OBBond*) *i)->GetNbrAtom(this):NULL);
}

OBAtom *OBAtom::NextNbrAtom(vector<OBEdgeBase*>::iterator &i) 
{
  i++; 
  return((i != _vbond.end()) ?  ((OBBond*) *i)->GetNbrAtom(this):NULL);
}


bool OBAtom::GetNewBondVector(vector3 &v,float length)
{
  // ***experimental code***

  OBAtom *atom;
  vector<OBEdgeBase*>::iterator i,j;
  v = VZero;

  if (GetValence() == 0)
    {
      v = VX;
      v *= length;
      v += GetVector();
      return(true);
    }

  if (GetValence() == 1)
    {
      vector3 vtmp,v1,v2;
      atom = BeginNbrAtom(i);
      if (atom)
	vtmp = GetVector() - atom->GetVector();

      if (GetHyb() == 2 || (IsOxygen() && HasAlphaBetaUnsat()))
	{
	  bool quit = false;
	  OBAtom *a1,*a2;
	  v2 = VZero;
	  for (a1 = BeginNbrAtom(i);a1 && !quit;a1 = NextNbrAtom(i))
	    for (a2 = a1->BeginNbrAtom(j);a2 && !quit;a2 = a1->NextNbrAtom(j))
	      if (a1 && a2 && a2 != this)
		{
		  v2 = a1->GetVector() - a2->GetVector();
		  quit = true;
		}

	  if (v2 == VZero)
	    {
	      v1 = cross(vtmp,VX);
	      v2 = cross(vtmp,VY);
	      if (v1.length() < v2.length()) v1 = v2;
	    }
	  else
	      v1 = cross(vtmp,v2);

	  matrix3x3 m;
	  m.RotAboutAxisByAngle(v1,60.0);
	  v = m*vtmp;
	  v.normalize();
	}
      else if (GetHyb() == 3)
	{
	  v1 = cross(vtmp,VX);
	  v2 = cross(vtmp,VY);
	  if (v1.length() < v2.length()) v1 = v2;
	  matrix3x3 m;
	  m.RotAboutAxisByAngle(v1,70.5);
	  v = m*vtmp;
	  v.normalize();
	}
      if (GetHyb() == 1) v = vtmp;

      v *= length;
      v += GetVector();
      return(true);
    }

  if (GetValence() == 2)
    {
      vector3 v1,v2,vtmp,vsum,vnorm;
      atom = BeginNbrAtom(i);if (!atom) return(false);
      v1 = GetVector() - atom->GetVector();
      atom = NextNbrAtom(i); if (!atom) return(false);
      v2 = GetVector() - atom->GetVector();
      v1.normalize();v2.normalize();
      vsum = v1+v2;
      vsum.normalize();

      if (GetHyb() == 2)
	v = vsum;
      else if (GetHyb() == 3)
	{
	  vnorm = cross(v2,v1);
	  vnorm.normalize();

#ifndef ONE_OVER_SQRT3
#define ONE_OVER_SQRT3  0.577350269f
#endif //SQRT_TWO_THIRDS
#ifndef SQRT_TWO_THIRDS
#define SQRT_TWO_THIRDS 0.816496581f
#endif //ONE_OVER_SQRT3

	  vsum *= ONE_OVER_SQRT3;
	  vnorm *= SQRT_TWO_THIRDS;

	  v = vsum + vnorm;
	}

      v *= length;

      v += GetVector();
      return(true);
    }

  if (GetValence() == 3)
    {
      vector3 vtmp,vsum;
      OBAtom *atom;
      vector<OBEdgeBase*>::iterator i;
      for (atom = BeginNbrAtom(i);atom;atom = NextNbrAtom(i))
	{
	  vtmp = GetVector() - atom->GetVector();
	  vtmp.normalize();
	  vtmp /= 3.0;
	  vsum += vtmp;
	}
      vsum.normalize();
      v = vsum;
      v *= length;
      v += GetVector();
      return(true);
    }

  return(true);
}

bool OBAtom::HtoMethyl()
{
  if (!IsHydrogen()) return(false);
  OBMol *mol = (OBMol*)GetParent();

  mol->BeginModify();

  SetAtomicNum(6); SetType("C3"); SetHyb(3);

  OBAtom *atom;
  OBBond *bond;
  vector<OBEdgeBase*>::iterator i;
  atom  = BeginNbrAtom(i);
  bond = (OBBond *)*i;
  if (!atom)
    {
      mol->EndModify();
      return(false);
    }

  float br1,br2;
  br1 = etab.CorrectedBondRad(6,3);
  br2 = etab.CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
  bond->SetLength(atom,br1+br2);

  OBAtom *hatom;
  br2 = etab.CorrectedBondRad(1,0);
  vector3 v;
  
  for (int j = 0;j < 3;j++)
    {
      hatom = mol->NewAtom();
      hatom->SetAtomicNum(1);
      hatom->SetType("H");

      GetNewBondVector(v,br1+br2); 
      hatom->SetVector(v);
      mol->AddBond(GetIdx(),mol->NumAtoms(),1);
    }

  mol->EndModify();
  return(true);
}

static void ApplyRotMatToBond(OBMol &mol,matrix3x3 &m,OBAtom *a1,OBAtom *a2)
{
  vector<int> children;
  mol.FindChildren(children,a1->GetIdx(),a2->GetIdx());
  children.push_back(a2->GetIdx());

  vector3 v;
  vector<int>::iterator i;
  for (i = children.begin();i != children.end();i++)
    {
      v = mol.GetAtom(*i)->GetVector();
      v -= a1->GetVector();
      v *= m;
      v += a1->GetVector();
      mol.GetAtom(*i)->SetVector(v);
    }

}

bool OBAtom::SetHybAndGeom(int hyb)
{
  //if (hyb == GetHyb()) return(true);
  if (GetAtomicNum() == 1) return(false);
  if (hyb == 0 && GetHvyValence() > 1) return(false);
  if (hyb == 1 && GetHvyValence() > 2) return(false);
  if (hyb == 2 && GetHvyValence() > 3) return(false);
  if (hyb == 3 && GetHvyValence() > 4) return(false);

  OBMol *mol = (OBMol*)GetParent();

  OBAtom *nbr;
  vector<OBAtom*> delatm;
  vector<OBEdgeBase*>::iterator i;

  for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
    if (nbr->IsHydrogen())
      delatm.push_back(nbr);
  
  //delete attached hydrogens
  mol->IncrementMod();
  vector<OBAtom*>::iterator j;
  for (j = delatm.begin();j != delatm.end();j++) mol->DeleteAtom(*j);
  mol->DecrementMod();

  float targetAngle;
  if (hyb == 3)     targetAngle = 109.5;
  else if (hyb == 2) targetAngle = 120.0;
  else if (hyb == 1) targetAngle = 180.0;
  else               targetAngle = 0.0;
  
  if (IsInRing())
    targetAngle = 180.0 - (360.0f / MemberOfRingSize());

  //adjust attached acyclic bond lengths
  float br1,br2, length;
  br1 = etab.CorrectedBondRad(GetAtomicNum(),hyb);
  for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
    if (!(*i)->IsInRing())
    {
      br2 = etab.CorrectedBondRad(nbr->GetAtomicNum(),nbr->GetHyb());
      length = br1 + br2;
      if ((*i)->IsAromatic()) length *= 0.93f;
      else if ((*i)->GetBO() == 2)         length *= 0.91f;
      else if ((*i)->GetBO() == 3)         length *= 0.87f;
      ((OBBond*) *i)->SetLength(this, length);
    }

  if (GetValence() > 1)
    {
      float angle;
      matrix3x3 m;
      vector3 v1,v2,v3,v4,n,s;
      OBAtom *r1,*r2,*r3,*a1,*a2,*a3,*a4;
      r1 = r2 = r3 = a1 = a2 = a3 = a4 = NULL;

      //find ring atoms first
      for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
	if ((*i)->IsInRing())
	  if (!r1)     r1 = nbr;
	  else if (!r2) r2 = nbr;
	  else if (!r3) r3 = nbr;

      //find non-ring atoms
      for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
	if (!(*i)->IsInRing())
	  if (!a1)     a1 = nbr;
	  else if (!a2) a2 = nbr;
	  else if (!a3) a3 = nbr;
	  else if (!a4) a4 = nbr;

      //adjust geometries of heavy atoms according to hybridization
      if (hyb == 1)
	{
	  if (a2)
	    {
	      v1 = a1->GetVector()-GetVector();v1.normalize();
	      v2 = a2->GetVector()-GetVector();v2.normalize();
	      n = cross(v1,v2);
	      angle = vectorAngle(v1,v2)-targetAngle;
	      m.RotAboutAxisByAngle(n,-angle);
	      ApplyRotMatToBond(*mol,m,this,a1);
	    }
	}
      else if (hyb == 2)
	{
	  if (r1 && r2 && a1)
	    {
	      v1 = r1->GetVector()-GetVector();v1.normalize();
	      v2 = r2->GetVector()-GetVector();v2.normalize();
	      v3 = a1->GetVector()-GetVector();
	      s = v1+v2; s.normalize(); s *= -1.0f;
	      n = cross(s,v3);
	      angle = vectorAngle(s,v3);
	      m.RotAboutAxisByAngle(n,angle);
	      ApplyRotMatToBond(*mol,m,this,a1);
	    }
	  else
	    {
	      if (a2)
		{
		  v1 = a1->GetVector()-GetVector();v1.normalize();
		  v2 = a2->GetVector()-GetVector();v2.normalize();
		  n = cross(v1,v2);
		  angle = vectorAngle(v1,v2)-targetAngle;
		  m.RotAboutAxisByAngle(n,-angle);
		  ApplyRotMatToBond(*mol,m,this,a1);
		}
	      if (a3)
		{
		  v1 = a1->GetVector()-GetVector();v1.normalize();
		  v2 = a2->GetVector()-GetVector();v2.normalize();
		  v3 = a3->GetVector()-GetVector();
		  s = v1+v2; s.normalize(); s *= -1.0f;
		  n = cross(s,v3);
		  angle = vectorAngle(s,v3);
		  m.RotAboutAxisByAngle(n,angle);
		  ApplyRotMatToBond(*mol,m,this,a3);
		}
	    }
	}
      else if (hyb == 3)
	{
	  if (r1 && r2 && r3 && a1)
	    {
	      v1 = r1->GetVector()-GetVector();v1.normalize();
	      v2 = r2->GetVector()-GetVector();v2.normalize();
	      v3 = r3->GetVector()-GetVector();v3.normalize();
	      v4 = a1->GetVector()-GetVector();
	      s = v1 + v2 + v3; s *= -1.0f; s.normalize();
	      n = cross(s,v4);
	      angle = vectorAngle(s,v4);
	      m.RotAboutAxisByAngle(n,angle);
	      ApplyRotMatToBond(*mol,m,this,a1);
	    }
	  else if (r1 && r2 && !r3 && a1)
	    {
	      v1 = r1->GetVector()-GetVector();v1.normalize();
	      v2 = r2->GetVector()-GetVector();v2.normalize();
	      v3 = a1->GetVector()-GetVector();
	      s = v1+v2; s *= -1.0f; s.normalize(); 
	      n = cross(v1,v2);      n.normalize();
	      s *= 0.577350269f; //1/sqrt(3)
	      n *= 0.816496581f; //sqrt(2/3)
	      s += n; s.normalize();
	      n = cross(s,v3);
	      angle = vectorAngle(s,v3);
	      m.RotAboutAxisByAngle(n,angle);
	      ApplyRotMatToBond(*mol,m,this,a1);

	      if (a2)
		{
		  v1 = r1->GetVector()-GetVector();v1.normalize();
		  v2 = r2->GetVector()-GetVector();v2.normalize();
		  v3 = a1->GetVector()-GetVector();v3.normalize();
		  v4 = a2->GetVector()-GetVector();
		  s = v1 + v2 + v3; s *= -1.0f; s.normalize();
		  n = cross(s,v4);
		  angle = vectorAngle(s,v4);
		  m.RotAboutAxisByAngle(n,angle);
		  ApplyRotMatToBond(*mol,m,this,a2);
		}
	    }
	  else
	    {
	      if (a2)
		{
		  v1 = a1->GetVector()-GetVector();v1.normalize();
		  v2 = a2->GetVector()-GetVector();v2.normalize();
		  n = cross(v1,v2);
		  angle = vectorAngle(v1,v2)-targetAngle;
		  m.RotAboutAxisByAngle(n,-angle);
		  ApplyRotMatToBond(*mol,m,this,a1);
		}
	      if (a3)
		{
		  v1 = a1->GetVector()-GetVector();v1.normalize();
		  v2 = a2->GetVector()-GetVector();v2.normalize();
		  v3 = a3->GetVector()-GetVector();
		  s = v1+v2; s *= -1.0f; s.normalize(); 
		  n = cross(v1,v2);      n.normalize();
		  s *= 0.577350269f; //1/sqrt(3)
		  n *= 0.816496581f; //sqrt(2/3)
		  s += n; s.normalize();
		  n = cross(s,v3);
		  angle = vectorAngle(s,v3);
		  m.RotAboutAxisByAngle(n,angle);
		  ApplyRotMatToBond(*mol,m,this,a3);
		}
	    }
	}
    }

  //add hydrogens back to atom
  int impval=1;
  switch (GetAtomicNum())
    {
    case 6: 
      if (hyb == 3) impval = 4;
      if (hyb == 2) impval = 3;
      if (hyb == 1) impval = 2;
      break;
    case 7: 
      if (hyb == 3) impval = 3;
      if (hyb == 2) impval = 2;
      if (hyb == 1) impval = 1;
      break;
    case 8: 
      if (hyb == 3) impval = 2;
      if (hyb == 2) impval = 2;
      if (hyb == 1) impval = 2;
    case 16: 
      if (hyb == 3) impval = 2;
      if (hyb == 2) impval = 2;
      if (hyb == 1) impval = 0;
      break;
    case 15: 
      if (hyb == 3) impval = 4;
      if (hyb == 2) impval = 3;
      if (hyb == 1) impval = 2;
      break;
    default:
      impval = 1;
    }

  int hcount = impval-GetHvyValence();
  if (hcount)
    {
      int k;
      vector3 v;
      OBAtom *atom;
      float brsum = etab.CorrectedBondRad(1,0)+
	etab.CorrectedBondRad(GetAtomicNum(),GetHyb());
      SetHyb(hyb);

      mol->BeginModify();
      for (k = 0;k < hcount;k++)
	{
	  GetNewBondVector(v,brsum);
	  atom = mol->NewAtom();
	  atom->SetAtomicNum(1);
	  atom->SetType("H");
	  atom->SetVector(v);
	  mol->AddBond(atom->GetIdx(),GetIdx(),1);
	}
      mol->EndModify();
    }


  return(true);
}

OBBond *OBAtom::GetBond(OBAtom *nbr)
{
  OBBond *bond;
  vector<OBEdgeBase *>::iterator i;
  for (bond = BeginBond(i) ; bond ; bond = NextBond(i))
    if (bond->GetNbrAtom(this) == nbr)
      return bond;
  return NULL;
}

}
