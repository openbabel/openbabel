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

#include "oeutil.h"

namespace OpenEye {

static int DetermineFRJ(OEMol &);
static void BuildOERTreeVector(OEAtom*,OERTree*,vector<OERTree*>&,OEBitVec&);


void OEMol::FindSSSR()
{
  if (HasSSSRPerceived()) return;
  SetSSSRPerceived();

  OERing *ring;
  vector<OERing*>::iterator j;

  //get frerejaque taking int account multiple possible spanning graphs
  int frj = DetermineFRJ(*this);
  if (frj)
    {
	  vector<OERing*> vr;
      FindRingAtomsAndBonds();

      OEBond *bond;
      vector<OEBond*> cbonds;
      vector<OEBond*>::iterator k;

      //restrict search for rings around closure bonds
      for (bond = BeginBond(k);bond;bond = NextBond(k))
	if (bond->IsClosure())
	  cbonds.push_back(bond);

      if (!cbonds.empty())
	{
	  OERingSearch rs;
	  //search for all rings about closures
	  vector<OEBond*>::iterator i;

	  for (i = cbonds.begin();i != cbonds.end();i++)
	    rs.AddRingFromClosure(*this,*i,0);
	  
	  rs.SortRings(); //sort ring sizes from smallest to largest
	  rs.RemoveRedundant(frj);  //full ring set - reduce to SSSR set
	  //store the SSSR set

	  for (j = rs.BeginRings();j != rs.EndRings();j++)
	    {
	      ring = new OERing ((*j)->_path,NumAtoms()+1);
	      ring->SetParent(this);
	      vr.push_back(ring);
	    }
	  //rs.WriteRings(); //for debugging only
	 }

	  if (!HasData(oeRingData)) SetData(new OERingData);
	  OERingData *rd = (OERingData*)GetData(oeRingData);
	  rd->SetData(vr);
    }
}

static int DetermineFRJ(OEMol &mol)
{
  vector<vector<int> >::iterator i;
  vector<vector<int> > cfl;
  //find all continuous graphs in the mol area
  mol.ContigFragList(cfl);
  
  if (cfl.empty()) return(0);
  if (cfl.size() == 1) return(mol.NumBonds() - mol.NumAtoms() + 1);

  //count up the atoms and bonds belonging to each graph
  OEBond *bond;
  vector<OEBond*>::iterator j;
  int numatoms,numbonds,frj=0;
  OEBitVec frag;
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

void OERingSearch::RemoveRedundant(int frj)
{
  OEBitVec tmp;
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


void OERingSearch::AddRingFromClosure(OEMol &mol,OEBond *cbond,int level)
{
  vector<OERTree*> t1(mol.NumAtoms()+1,(OERTree*)NULL);
  vector<OERTree*> t2(mol.NumAtoms()+1,(OERTree*)NULL);
  OEBitVec bv1,bv2;

  bv1.SetBitOn(cbond->GetEndAtomIdx());
  bv2.SetBitOn(cbond->GetBeginAtomIdx());
  BuildOERTreeVector(cbond->GetBeginAtom(),NULL,t1,bv1);
  BuildOERTreeVector(cbond->GetEndAtom(),NULL,t2,bv2);

  bool pathok;
  deque<int> p1,p2;
  vector<OEAtom*> path1,path2;
  vector<OEAtom*>::iterator m,n;
  vector<OERTree*>::iterator i;

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

  //clean up OERTree vectors
  for (i = t1.begin();i != t1.end();i++)
    if (*i) delete *i;
      
  for (i = t2.begin();i != t2.end();i++)
    if (*i) delete *i;
}

bool OERingSearch::SaveUniqueRing(deque<int> &d1,deque<int> &d2)
{
  vector<int> path;
  OEBitVec bv;
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

  vector<OERing*>::iterator j;
  for (j = _rlist.begin();j != _rlist.end();j++)
    if (bv == (*j)->_pathset)
      return(false);

  OERing *ring = new OERing;
  ring->_path = path;
  ring->_pathset = bv;
  _rlist.push_back(ring);

  return(true);
}

OERingSearch::~OERingSearch()
{
  vector<OERing*>::iterator i;
  for (i = _rlist.begin();i != _rlist.end();i++) delete *i;
}

bool CompareRingSize(const OERing *a,const OERing *b)
{
  return(a->PathSize() < b->PathSize());
}

void OERingSearch::WriteRings()
{
  vector<OERing*>::iterator i;

  for (i = _rlist.begin();i != _rlist.end();i++)
    cout << (*i)->_pathset << endl;
}

static void FindRings(OEMol &mol,vector<int> &path,OEBitVec &avisit,
		      OEBitVec &bvisit, int natom,int depth );

void OEMol::FindRingAtomsAndBonds()
{
  if (HasFlag(OE_RINGFLAGS_MOL)) return;
  SetFlag(OE_RINGFLAGS_MOL);

  OEBitVec avisit,bvisit;
  avisit.Resize(NumAtoms()+1);
  bvisit.Resize(NumAtoms()+1);
  vector<int> path;
  path.resize(NumAtoms()+1);

  for(unsigned int i=1; i<= NumAtoms(); i++ )
    if(!avisit[i])
      FindRings(*this,path,avisit,bvisit,i,0);
}

static void FindRings(OEMol &mol,vector<int> &path,OEBitVec &avisit,
		      OEBitVec &bvisit, int natom,int depth )
{
  OEAtom *atom;
  OEBond *bond;
  vector<OEBond*>::iterator k;

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

bool OERing::IsAromatic()
{
  OEMol *mol = _parent;
  vector<int>::iterator i;
  for (i = _path.begin();i != _path.end();i++)
    if (!(mol->GetAtom(*i))->IsAromatic())
      return(false);

  return(true);
}

bool OERing::IsMember(OEAtom *a)
{
  return(_pathset.BitIsOn(a->GetIdx()));
}

bool OERing::IsMember(OEBond *b)
{
	return((_pathset.BitIsOn(b->GetBeginAtomIdx()))&&(_pathset.BitIsOn(b->GetEndAtomIdx())));
}

OERing::OERing(vector<int> &path,int size)
{
  _path = path;
  _pathset.FromVecInt(_path);
  _pathset.Resize(size);
}

/*!
**\brief OERing copy constructor
**\param src reference to original OERing object (rhs)
*/
OERing::OERing(const OERing &src)
//no base class
	:	_pathset(src._pathset),	//chain to member classes
		_path(src._path)
{
	//member data
	_parent = src._parent;  //this is messed up, but what can you do?
}

/*!
**\brief OERing assignemtn operator
**\param src reference to original OERing object (rhs)
**\return reference to modified OERing object (lhs)
*/
OERing& OERing::operator =(const OERing &src)
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
void BuildOERTreeVector(OEAtom *atom,OERTree *prv,vector<OERTree*> &vt,OEBitVec &bv)
{
  vt[atom->GetIdx()] = new OERTree (atom,prv);

  int i;
  OEAtom *nbr;
  OEMol *mol = (OEMol*)atom->GetParent();
  OEBitVec curr,used,next;
  vector<OEBond*>::iterator j;
  curr |= atom->GetIdx();
  used = bv|curr;

#define OE_RTREE_CUTOFF 20
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
		vt[nbr->GetIdx()] = new OERTree (nbr,vt[atom->GetIdx()]);
	      }
	}
      
      if (next.Empty()) break;
      curr = next;
      level++;
      if (level > OE_RTREE_CUTOFF) break;
    }
#undef OE_RTREE_CUTOFF
}

OERTree::OERTree(OEAtom *atom,OERTree *prv)
{
  _atom = atom;
  _prv = prv;
}

void OERTree::PathToRoot(vector<OEAtom*> &path)
{
  path.push_back(_atom);
  if (_prv) _prv->PathToRoot(path);
}

int OERTree::GetAtomIdx()
{
  return(_atom->GetIdx());
}

bool OERing::findCenterAndNormal(Vector & center, Vector &norm1, Vector &norm2)
{
    OEMol *mol= this->_parent;
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
