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
#include "phmodel.h"
#include "phmodeldata.h"

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

using namespace std;

namespace OpenBabel
{

OBPhModel phmodel;
extern OBAtomTyper atomtyper;

OBPhModel::OBPhModel()
{
  _init = false;
  _dir = DATADIR;
  _envvar = "BABEL_DATADIR";
  _filename = "phmodel.txt";
  _subdir = "data";
  _dataptr = PhModelData;
}

OBPhModel::~OBPhModel()
{
  vector<OBChemTsfm*>::iterator k;
  for (k = _vtsfm.begin();k != _vtsfm.end();k++)     delete *k;
    
  vector<pair<OBSmartsPattern*,vector<float> > >::iterator m;
  for (m = _vschrg.begin();m != _vschrg.end();m++)   delete m->first;
}

void OBPhModel::ParseLine(const char *buffer)
{
  vector<string> vs;
  OBSmartsPattern *sp;

  if (buffer[0] == '#') return;

  if (EQn(buffer,"TRANSFORM",7))
    {
      tokenize(vs,buffer);
      if (vs.empty() || vs.size() < 4) return;

      OBChemTsfm *tsfm = new OBChemTsfm;
      if (!tsfm->Init(vs[1],vs[3])) 
	{
	  delete tsfm; tsfm = NULL;
	  return;
	}

      _vtsfm.push_back(tsfm);
    }
  else if (EQn(buffer,"SEEDCHARGE",10))
    {
      tokenize(vs,buffer);
      if (vs.empty() || vs.size() < 2) return;
      sp = new OBSmartsPattern;
      if (!sp->Init(vs[1]) || ((signed)vs.size()-2) != sp->NumAtoms())
	{
	  delete sp; sp = NULL;
	  return;
	}
      vector<float> vf;
      vector<string>::iterator i;
      for (i = vs.begin()+2;i != vs.end();i++)
	vf.push_back(atof((char*)i->c_str()));

      _vschrg.push_back(pair<OBSmartsPattern*,vector<float> > (sp,vf));
    }
}

void OBPhModel::AssignSeedPartialCharge(OBMol &mol)
{
  if (!_init) Init();

  mol.SetPartialChargesPerceived();
  if (!mol.AutomaticPartialCharge()) return;

  vector<pair<OBSmartsPattern*,vector<float> > >::iterator i;
  for (i = _vschrg.begin();i != _vschrg.end();i++)
    if (i->first->Match(mol))
    {
      _mlist = i->first->GetUMapList();
      unsigned int k;
      vector<vector<int> >::iterator j;

      for (j = _mlist.begin();j != _mlist.end();j++)
	for (k = 0;k < j->size();k++)
	  mol.GetAtom((*j)[k])->SetPartialCharge(i->second[k]);
    }
}

void OBPhModel::CorrectForPH(OBMol &mol)
{
  if (!_init) Init();
  if (mol.IsCorrectedForPH()) return;
  if (!mol.AutomaticFormalCharge()) return;

  mol.SetCorrectedForPH();

  OBAtom *atom;
  vector<OBNodeBase*>::iterator j;
  for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
    atom->SetFormalCharge(0);

  vector<OBChemTsfm*>::iterator i;
  for (i = _vtsfm.begin();i != _vtsfm.end();i++) (*i)->Apply(mol);

  atomtyper.CorrectAromaticNitrogens(mol);
}

bool OBChemTsfm::Init(string &bgn,string &end)
{
  if (!_bgn.Init(bgn)) return(false);
  if (!end.empty())
    if (!_end.Init(end)) return(false);

  //find atoms to be deleted
  int i,j,vb;
  bool found;
  for (i = 0;i < _bgn.NumAtoms();i++)
    if ((vb = _bgn.GetVectorBinding(i)))
      {
	found = false;
	for (j = 0;j < _end.NumAtoms();j++)
	  if (vb == _end.GetVectorBinding(j))
	    {found = true; break;}

	if (!found)
	  _vadel.push_back(i);
      }

  //find elements to be changed
  int ele;
  for (i = 0;i < _bgn.NumAtoms();i++)
    if ((vb = _bgn.GetVectorBinding(i)))
      {
	ele = _bgn.GetAtomicNum(i);
	for (j = 0;j < _end.NumAtoms();j++)
	  if (vb == _end.GetVectorBinding(j))
	    if (ele != _end.GetAtomicNum(j))
	      {
		_vele.push_back(pair<int,int> (i,_end.GetAtomicNum(j))); 
		break;
	      }
      }

  //find charges to modify
  int chrg;
  for (i = 0;i < _bgn.NumAtoms();i++)
	  if ((vb = _bgn.GetVectorBinding(i)))
      {
		  chrg = _bgn.GetFormalCharge(i);
		  for (j = 0;j < _end.NumAtoms();j++)
			  if (vb == _end.GetVectorBinding(j))
				  if (chrg != _end.GetFormalCharge(j))
					  _vchrg.push_back(pair<int,int> (i,_end.GetFormalCharge(j)));
      }

  //find bonds to be modified
  OBEdge *edge1,*edge2;
  OBNode *bgn1,*bgn2,*end1,*end2;
  vector<OBEdgeBase*>::iterator k,m;
  for (edge1 = _bgn.Begin(k);edge1;edge1 = _bgn.Next(k))
  {
	  bgn1 = (OBNode*)edge1->GetBgn();
	  end1 = (OBNode*)edge1->GetEnd();
	  if (!bgn1->GetVectorBinding()) continue;
	  if (!end1->GetVectorBinding()) continue;

	  for (edge2 = _end.Begin(m);edge2;edge2 = _end.Next(m))
	  {
		  bgn2 = (OBNode*)edge2->GetBgn();
		  end2 = (OBNode*)edge2->GetEnd();
		  //check to see edge1 == edge2
		  if ((bgn1->GetVectorBinding() == bgn2->GetVectorBinding() &&
			   end1->GetVectorBinding() == end2->GetVectorBinding()) || 
			  (bgn1->GetVectorBinding() == end2->GetVectorBinding() &&
			   end1->GetVectorBinding() == bgn2->GetVectorBinding()))
		  {
			if (edge1->GetBO() != edge2->GetBO()) //no mod necessary
			_vbond.push_back(pair<pair<int,int>,int> (pair<int,int> (bgn1->GetIdx(),end1->GetIdx()),edge2->GetBO()));
			break;
		  }
	  }
  }

  //make sure there is some kind of transform to do here
  if (_vadel.empty() && _vchrg.empty() && _vbond.empty()) return(false);

  return(true);
}

bool OBChemTsfm::Apply(OBMol &mol)
{
  if (!_bgn.Match(mol)) return(false);

  vector<vector<int> > mlist = _bgn.GetUMapList();

  if (!_vchrg.empty()) //modify charges
    {
      vector<vector<int> >::iterator i;
      vector<pair<int,int> >::iterator j;

      for (i = mlist.begin();i != mlist.end();i++)
	for (j = _vchrg.begin();j != _vchrg.end();j++)
	  if (j->first < (signed)i->size()) //goof proofing
	    mol.GetAtom((*i)[j->first])->SetFormalCharge(j->second);

      mol.UnsetImplicitValencePerceived();
    }

  if (!_vbond.empty()) //modify bond orders
    {
      OBBond *bond;
      vector<vector<int> >::iterator i;
      vector<pair<pair<int,int>,int> >::iterator j;      
      for (i = mlist.begin();i != mlist.end();i++)
	for (j = _vbond.begin();j != _vbond.end();j++)
	  {
	    bond = mol.GetBond((*i)[j->first.first],(*i)[j->first.second]);
	    if (!bond)
	      {
		ThrowError("unable to find bond");
		continue;
	      }

	    bond->SetBO(j->second);
	  }
    }

  if (!_vadel.empty() || !_vele.empty()) //delete atoms and change elements
    {
      vector<int>::iterator j;
      vector<vector<int> >::iterator i;

      if (!_vele.empty())
	  {
		  vector<pair<int,int> >::iterator k;
		  for (i = mlist.begin();i != mlist.end();i++)
			  for (k = _vele.begin();k != _vele.end();k++)
				  mol.GetAtom((*i)[k->first])->SetAtomicNum(k->second);
	  }

      //make sure same atom isn't delete twice
      vector<bool> vda;
      vector<OBNodeBase*> vdel;
      vda.resize(mol.NumAtoms()+1,false);
      for (i = mlist.begin();i != mlist.end();i++)
		  for (j = _vadel.begin();j != _vadel.end();j++)
			  if (!vda[(*i)[*j]])
			  {
				  vda[(*i)[*j]] = true;
				  vdel.push_back(mol.GetAtom((*i)[*j]));
			  }

	  vector<OBNodeBase*>::iterator k;
      for (k = vdel.begin();k != vdel.end();k++)
		  mol.DeleteAtom((OBAtom*)*k); 
    }
  
  return(true);
}

} //namespace OpenBabel
