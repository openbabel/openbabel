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
#include <deque>

#include "obutil.h"

namespace OpenBabel {

static int DetermineFRJ(OBMol &);
static void BuildOBRTreeVector(OBAtom*,OBRTree*,vector<OBRTree*>&,OBBitVec&);


void OBMol::FindSSSR()
{
  if (HasSSSRPerceived()) return;
  SetSSSRPerceived();

  OBRing *ring;
  vector<OBRing*>::iterator j;

  //get frerejaque taking int account multiple possible spanning graphs
  int frj = DetermineFRJ(*this);
  if (frj)
    {
	  vector<OBRing*> vr;
      FindRingAtomsAndBonds();

      OBBond *bond;
      vector<OBEdgeBase*> cbonds;
      vector<OBEdgeBase*>::iterator k;

      //restrict search for rings around closure bonds
      for (bond = BeginBond(k);bond;bond = NextBond(k))
	if (bond->IsClosure())
	  cbonds.push_back(bond);

      if (!cbonds.empty())
	{
	  OBRingSearch rs;
	  //search for all rings about closures
	  vector<OBEdgeBase*>::iterator i;

	  for (i = cbonds.begin();i != cbonds.end();i++)
	    rs.AddRingFromClosure(*this,(OBBond*)*i,0);
	  
	  rs.SortRings(); //sort ring sizes from smallest to largest
	  rs.RemoveRedundant(frj);  //full ring set - reduce to SSSR set
	  //store the SSSR set

	  for (j = rs.BeginRings();j != rs.EndRings();j++)
	    {
	      ring = new OBRing ((*j)->_path,NumAtoms()+1);
	      ring->SetParent(this);
	      vr.push_back(ring);
	    }
	  //rs.WriteRings(); //for debugging only
	 }

	  if (!HasData(obRingData)) SetData(new OBRingData);
	  OBRingData *rd = (OBRingData*)GetData(obRingData);
	  rd->SetData(vr);
    }
}

static int DetermineFRJ(OBMol &mol)
{
  vector<vector<int> >::iterator i;
  vector<vector<int> > cfl;
  //find all continuous graphs in the mol area
  mol.ContigFragList(cfl);
  
  if (cfl.empty()) return(0);
  if (cfl.size() == 1) return(mol.NumBonds() - mol.NumAtoms() + 1);

  //count up the atoms and bonds belonging to each graph
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;
  int numatoms,numbonds,frj=0;
  OBBitVec frag;
  for (i = cfl.begin();i != cfl.end();i++)
    {
      frag.Clear();
      frag.FromVecInt(*i);
      numatoms = (*i).size();numbonds=0;
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
  OBBitVec tmp;
  register int i,j;

  //remove identical rings
  for (i = _rlist.size()-1;i > 0;i--)
    for (j = i-1;j >= 0;j--)
      if ((_rlist[i])->_pathset == (_rlist[j])->_pathset)
	{
	  delete _rlist[i];
	  _rlist.erase(_rlist.begin()+i);
	  break;
	}

  //make sure tmp is the same size as the rings
  for (j = 0;j < (signed)_rlist.size();j++)
    tmp = (_rlist[j])->_pathset;

  //remove larger rings that cover the same atoms as smaller rings
  for (i = _rlist.size()-1;i >= 0;i--)
  {
    tmp.Clear();
    for (j = 0;j < (signed)_rlist.size();j++)
      if ((_rlist[j])->_path.size() <= (_rlist[i])->_path.size() && i != j)
	  tmp |= (_rlist[j])->_pathset;
    
    tmp = tmp & (_rlist[i])->_pathset;
    
    if (tmp == (_rlist[i])->_pathset)
      {
	delete _rlist[i];
	_rlist.erase(_rlist.begin()+i);
      }

    if (_rlist.size() == (unsigned)frj) break;
  }
}


void OBRingSearch::AddRingFromClosure(OBMol &mol,OBBond *cbond,int level)
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
  vector<OBNodeBase*> path1,path2;
  vector<OBNodeBase*>::iterator m,n;
  vector<OBRTree*>::iterator i;

  for (i = t1.begin();i != t1.end();i++)
    if (*i)
    {
      path1.clear();
      (*i)->PathToRoot(path1);

      if (t2[(*i)->GetAtomIdx()])
	{
	  pathok = true;
	  path2.clear();
	  t2[(*i)->GetAtomIdx()]->PathToRoot(path2);

	  p1.clear();m = path1.begin();
	  if (m != path1.end()) p1.push_back((*m)->GetIdx());
	  for (m = path1.begin(),m++;m != path1.end();m++)
	    {
	      p1.push_back((*m)->GetIdx());
	      p2.clear();
	      for (n = path2.begin(),n++;n != path2.end();n++)
		{
		  p2.push_front((*n)->GetIdx());
		  if (*n == *m)//don't traverse across identical atoms
		    {
		      p2.pop_front();
		      if (p1.size()+p2.size() > 2) SaveUniqueRing(p1,p2);
		      pathok = false;
		      break;  
		    }
		  if ((*n)->IsConnected(*m) && p1.size()+p2.size() > 2)
		    SaveUniqueRing(p1,p2);
		}
	      if (!pathok) break;
	    }
	}
    }

  //clean up OBRTree vectors
  for (i = t1.begin();i != t1.end();i++)
    if (*i) delete *i;
      
  for (i = t2.begin();i != t2.end();i++)
    if (*i) delete *i;
}

bool OBRingSearch::SaveUniqueRing(deque<int> &d1,deque<int> &d2)
{
  vector<int> path;
  OBBitVec bv;
  deque<int>::iterator i;

  for (i = d1.begin();i != d1.end();i++)
    {
      bv.SetBitOn(*i);
      path.push_back(*i);
    }

  for (i = d2.begin();i != d2.end();i++)
    {
      bv.SetBitOn(*i);
      path.push_back(*i);
    }

  vector<OBRing*>::iterator j;
  for (j = _rlist.begin();j != _rlist.end();j++)
    if (bv == (*j)->_pathset)
      return(false);

  OBRing *ring = new OBRing;
  ring->_path = path;
  ring->_pathset = bv;
  _rlist.push_back(ring);

  return(true);
}

OBRingSearch::~OBRingSearch()
{
  vector<OBRing*>::iterator i;
  for (i = _rlist.begin();i != _rlist.end();i++) delete *i;
}

bool CompareRingSize(const OBRing *a,const OBRing *b)
{
  return(a->PathSize() < b->PathSize());
}

void OBRingSearch::WriteRings()
{
  vector<OBRing*>::iterator i;

  for (i = _rlist.begin();i != _rlist.end();i++)
    cout << (*i)->_pathset << endl;
}

static void FindRings(OBMol &mol,vector<int> &path,OBBitVec &avisit,
		      OBBitVec &bvisit, int natom,int depth );

void OBMol::FindRingAtomsAndBonds()
{
  if (HasFlag(OB_RINGFLAGS_MOL)) return;
  SetFlag(OB_RINGFLAGS_MOL);

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
  vector<OBEdgeBase*>::iterator k;

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
	  if(bond->GetBeginAtomIdx()==natom || bond->GetEndAtomIdx()==natom)
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
  for (i = _path.begin();i != _path.end();i++)
    if (!(mol->GetAtom(*i))->IsAromatic())
      return(false);

  return(true);
}

bool OBRing::IsMember(OBAtom *a)
{
  return(_pathset.BitIsOn(a->GetIdx()));
}

bool OBRing::IsMember(OBBond *b)
{
	return((_pathset.BitIsOn(b->GetBeginAtomIdx()))&&(_pathset.BitIsOn(b->GetEndAtomIdx())));
}

OBRing::OBRing(vector<int> &path,int size)
{
  _path = path;
  _pathset.FromVecInt(_path);
  _pathset.Resize(size);
}

/*!
**\brief OBRing copy constructor
**\param src reference to original OBRing object (rhs)
*/
OBRing::OBRing(const OBRing &src)
//no base class
	:	_pathset(src._pathset),	//chain to member classes
		_path(src._path)
{
	//member data
	_parent = src._parent;  //this is messed up, but what can you do?
}

/*!
**\brief OBRing assignemtn operator
**\param src reference to original OBRing object (rhs)
**\return reference to modified OBRing object (lhs)
*/
OBRing& OBRing::operator =(const OBRing &src)
{
	//on identity, return
	if(this == &src)	return(*this);

	//no base class

	//memeber classes & data
	_path = src._path;
	_pathset = src._pathset;
	_parent = src._parent; //note, this may not be what you want

	return(*this);
}
void BuildOBRTreeVector(OBAtom *atom,OBRTree *prv,vector<OBRTree*> &vt,OBBitVec &bv)
{
  vt[atom->GetIdx()] = new OBRTree (atom,prv);

  int i;
  OBAtom *nbr;
  OBMol *mol = (OBMol*)atom->GetParent();
  OBBitVec curr,used,next;
  vector<OBEdgeBase*>::iterator j;
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
      
      if (next.Empty()) break;
      curr = next;
      level++;
      if (level > OB_RTREE_CUTOFF) break;
    }
#undef OB_RTREE_CUTOFF
}

OBRTree::OBRTree(OBAtom *atom,OBRTree *prv)
{
  _atom = atom;
  _prv = prv;
}

void OBRTree::PathToRoot(vector<OBNodeBase*> &path)
{
  path.push_back(_atom);
  if (_prv) _prv->PathToRoot(path);
}

int OBRTree::GetAtomIdx()
{
  return(_atom->GetIdx());
}

bool OBRing::findCenterAndNormal(Vector & center, Vector &norm1, Vector &norm2)
{
    OBMol *mol= this->_parent;
    int j= 0;
    const int nA= this->_path.size();
    Vector tmp;

    center.Set(0.0f,0.0f,0.0f);
    norm1.Set(0.0f,0.0f,0.0f);
    norm2.Set(0.0f,0.0f,0.0f);
    for (j = 0; j != nA; j++){
       center += (mol->GetAtom(_path[j]))->GetVector();
    }
    center/= float(nA);

    for (j = 0; j != nA; j++){
          Vector v1= (mol->GetAtom(_path[j]))->GetVector() - center;
          Vector v2= (mol->GetAtom(_path[j+1==nA?0:j+1]))->GetVector() - center;
          tmp= cross(v1,v2);
          norm1+= tmp;
    }
    norm1/= float(nA);
    norm1.normalize();
    norm2= norm1;
    norm2 *= -1.0f;
    return(true);
}

}
