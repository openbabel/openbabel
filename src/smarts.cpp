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

#include <iostream>
#include <fstream>

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include "base.h"
#include "smarts.h"
#include "parsmart.h"
#include "bitvec.h"

namespace OpenBabel
{

//
//OBSmartsPattern member functions
//

bool OBSmartsPattern::SeedMatch(OBNodeBase *nb)
//routine ONLY used for matching recursive SMARTS
{
	_done = false;
	_single = true;
	_vmatch.clear();
	OBGraphBase *g = nb->GetParent();
	vector<bool> visit(g->NumNodes()); 

	//save visit flags as they are being used further up in the stack
	int i;
	OBNodeBase *node;
	vector<OBNodeBase*>::iterator j;
	for (i=0,node = g->Begin(j);node;node = g->Next(j))
		visit[i++] = node->Visit;
	g->ResetVisitFlags();

	OBNodeBase *seed = GetFirstSeed();

	if (seed->Eval(nb))
	{
		nb->Visit = true;
		seed->SetMatch(nb);
		Match(*g,BgnMatch(),BgnMatch()->second.begin());
	}

	//set visit flags to original values
	for (i=0,node = g->Begin(j);node;node = g->Next(j))
		node->Visit = visit[i++];

	return(FinishedMatch());
}

void OBSmartsPattern::GetMatches(vector<vector<int> > &vm)
{
	vector<OBNodeBase*>::iterator j;
	vector<vector<OBNodeBase*> >::iterator i;
	for (i = _vmatch.begin();i != _vmatch.end();i++)
	{
		vector<int> vtmp;
		vtmp.reserve(i->size());
		for (j = i->begin();j != i->end();j++)
			vtmp.push_back((*j)->GetIdx());
		vm.push_back(vtmp);
	}
}

vector<vector<int> > OBSmartsPattern::GetMapList()
{
	vector<vector<int> > vm;
	vector<OBNodeBase*>::iterator j;
	vector<vector<OBNodeBase*> >::iterator i;
	for (i = _vmatch.begin();i != _vmatch.end();i++)
	{
		vector<int> vtmp;
		vtmp.reserve(i->size());
		for (j = i->begin();j != i->end();j++)
			vtmp.push_back((*j)->GetIdx());
		vm.push_back(vtmp);
	}

	return(vm);
}

vector<vector<int> > OBSmartsPattern::GetUMapList()
{
	vector<vector<int> > vm; 
	vm.clear();

	if (!_vmatch.empty())
	{
		bool ok;
		OBBitVec bv;
		vector<OBBitVec> vbv;
		vector<OBBitVec>::iterator k;
		vector<OBNodeBase*>::iterator j;
		vector<vector<OBNodeBase*> >::iterator i;
		for (i = _vmatch.begin();i != _vmatch.end();i++)
		{
			vector<int> vtmp;
			vtmp.reserve(i->size());
			bv.Clear();
			for (j = i->begin();j != i->end();j++)
			{
				bv |= (*j)->GetIdx();
				vtmp.push_back((*j)->GetIdx());
			}
			ok = true;
			for (k = vbv.begin();k != vbv.end() && ok;k++)
				if (*k == bv)
					ok = false;
			
			if (ok) 
			{
				vm.push_back(vtmp);
				vbv.push_back(bv);
			}
		}
	}

	return(vm);
}

void OBSmartsPattern::PrepForMatch()
{
	_vm.clear();
	SetVisitLock(true);
	ResetVisitFlags();

	OBEdge *edge;
	OBNode *node,*seed,*nbr;
	vector<OBNodeBase*> curr,next;
	vector<OBNodeBase*>::iterator i,j;
	vector<OBEdgeBase*>::iterator k;
	vector<OBEdgeBase*> vedge;

	for (node = Begin(i);node;node = Next(i)) //make sure all atoms are covered
		if (!node->Visit)
		{
			seed = node;
			node->Visit = true;
			curr.push_back(node);
			vedge.clear();

			//traverse graph starting from seed atom and store order of bonds to match
			for (;!curr.empty();)
			{
				next.clear();
				for (j = curr.begin();j != curr.end();j++)
				{
					for (nbr = ((OBNode*)*j)->BeginNbr(k);nbr;nbr = ((OBNode*)*j)->NextNbr(k))
						if (!nbr->Visit)
						{
							if ((*k)->GetBgn() != *j) (*k)->SwapEnds();
							next.push_back(nbr);
							nbr->Visit = true;
							(*k)->Visit = true;
							vedge.push_back(*k);
						}
				}
				curr = next;
			}

			//add closure bonds to list
			for (edge = Begin(k);edge;edge = Next(k))
				if (!edge->Visit && edge->GetBgn()->Visit && edge->GetEnd()->Visit)
				{
					edge->Visit = true;
					edge->SetClosure();
					vedge.push_back(edge);
				}

			//store bonds associated with seed atom
			_vm.push_back(pair<OBNodeBase*,vector<OBEdgeBase*> > (seed,vedge));
		}

	SetVisitLock(false);
}

bool OBSmartsPattern::Init(const char *buf)
{
	OBSmartsParser sp;
	return(sp.Parse(*this,buf));
}

bool OBSmartsPattern::Init(string &s)
{
	OBSmartsParser sp;
	return(sp.Parse(*this,s));
}

OBNode *OBSmartsPattern::NewNode(OBExprBase *xb)
{
	OBNode *node = new OBNode;
	node->SetExpr(xb);
	node->SetIdx(NumNodes());
	node->SetParent(this);
	_vatom.push_back(node);

	return(node);
}

OBEdge* OBSmartsPattern::NewEdge(OBNode *bgn,OBNode *end,OBExprBase *expr)
{
	OBEdge *edge = new OBEdge(bgn,end,expr);
	edge->SetIdx(NumEdges());
	edge->SetParent(this);
	_vbond.push_back(edge);
	if (bgn) bgn->AddEdge(edge);
	if (end) end->AddEdge(edge);
	return(edge);
}

void OBSmartsPattern::ClearMatches()
{
	_vmatch.clear(); //clear out previous matches

	OBNode *node;
	vector<OBNodeBase*>::iterator i;
	for (node = Begin(i);node;node = Next(i))
		node->ResetRecurs();
}

OBNode *OBSmartsPattern::Begin(vector<OBNodeBase*>::iterator &i)
{
	i = _vatom.begin();
	return((i != _vatom.end()) ? (OBNode*)*i : (OBNode*)NULL);
}

OBNode *OBSmartsPattern::Next(vector<OBNodeBase*>::iterator &i)
{
	i++;
	return((i != _vatom.end()) ? (OBNode*)*i : (OBNode*)NULL);
}

OBEdge *OBSmartsPattern::Begin(vector<OBEdgeBase*>::iterator &i) 
{
	i = _vbond.begin();
	return((i != _vbond.end()) ? (OBEdge*)*i : (OBEdge*)NULL);
}

OBEdge *OBSmartsPattern::Next(vector<OBEdgeBase*>::iterator &i)
{
	i++;
	return((i != _vbond.end()) ? (OBEdge*)*i : (OBEdge*)NULL);
}

OBNode *OBNode::BeginNbr(vector<OBEdgeBase*>::iterator &i)
{
	i = _vbond.begin();
	if (i == _vbond.end()) return((OBNode*)NULL);
	return((this == ((OBEdge*)*i)->GetBgn()) ? (OBNode*)((OBEdge*)*i)->GetEnd() : (OBNode*)((OBEdge*)*i)->GetBgn());
}

OBNode *OBNode::NextNbr(vector<OBEdgeBase*>::iterator &i)
{
	i++;
	if (i == _vbond.end()) return((OBNode*)NULL);
	return((this == ((OBEdge*)*i)->GetBgn()) ? (OBNode*)((OBEdge*)*i)->GetEnd() : (OBNode*)((OBEdge*)*i)->GetBgn());
}

OBSmartsPattern::OBSmartsPattern()
{
  _done = false;
  _single = false;
  _smarts = "";
  _vmatch.clear();
  _vm.clear();
}

OBSmartsPattern::OBSmartsPattern(const OBSmartsPattern &src)
{
  string s = src.GetSMARTS();
  Init(s);
  _done = src.FinishedMatch();
  _single = src.SingleMatch();
  src.GetMatches(_vmatch);
}

OBSmartsPattern::~OBSmartsPattern()
{
	vector<OBNodeBase*>::iterator i;
	for (i = _vatom.begin();i != _vatom.end();i++) delete *i;

	vector<OBEdgeBase*>::iterator j;
	for (j = _vbond.begin();j != _vbond.end();j++) delete *j;
}

//
//Expression member functions
//

int OBAndExpr::GetFormalCharge()
{
	int chrg = _lft->GetFormalCharge();
	if (chrg != 0) return(chrg);
	return(_rgt->GetFormalCharge());
}

unsigned int OBAndExpr::GetAtomicNum()
{
	unsigned int ele = _lft->GetAtomicNum();
	if (ele != 0) return(ele);
	return(_rgt->GetAtomicNum());
}

unsigned int OBAndExpr::GetBO()
{
	unsigned int bo = _lft->GetBO();
	if (bo != 0) return(bo);
	return(_rgt->GetBO());
}

bool OBAndExpr::IsAromatic()
{
	if (_lft->IsAromatic()) return(true);
	return(_rgt->IsAromatic());
}

void OBAndExpr::Print(ostream &ofs)
{
	ofs << "((";
	_lft->Print(ofs);
	ofs << ")&(";
	_rgt->Print(ofs);
	ofs << "))";
}

void OBAndExpr::ResetRecurs()
{
	_lft->ResetRecurs();
	_rgt->ResetRecurs();
}

int OBOrExpr::GetFormalCharge()
{
	int chrg = _lft->GetFormalCharge();
	if (chrg != 0) return(chrg);
	return(_rgt->GetFormalCharge());
}
unsigned int OBOrExpr::GetAtomicNum()
{
	unsigned int ele = _lft->GetAtomicNum();
	if (ele != 0) return(ele);
	return(_rgt->GetAtomicNum());
}

unsigned int OBOrExpr::GetBO()
{
	unsigned int bo = _lft->GetBO();
	if (bo != 0) return(bo);
	return(_rgt->GetBO());
}
bool OBOrExpr::IsAromatic() 
{
	return(_lft->IsAromatic() || _rgt->IsAromatic());
}

void OBOrExpr::Print(ostream &ofs)  
{
	ofs << "((";
	_lft->Print(ofs);
	ofs << "),(";
	_rgt->Print(ofs);
	ofs << "))";
}

void OBOrExpr::ResetRecurs()
{
	_lft->ResetRecurs();
	_rgt->ResetRecurs();
}

void OBNotExpr::Print(ostream &ofs)
{
	ofs << "!(";
	_expr->Print(ofs);
	ofs << ")";
}

void OBNotExpr::ResetRecurs()
{
	_expr->ResetRecurs();
}

void OBConstExpr::Print(ostream &ofs)
{
	ofs << "(TRUE)";
}

void OBElementExpr::Print(ostream &ofs)
{
	ofs << "#" << _ele;
}

void OBAromaticExpr::Print(ostream &ofs)  
{
	if (_val) ofs << "Aromatic";
	else      ofs << "Aliphatic";
}

bool OBRingExpr::Eval(OBNodeBase *nb) 
{
  if (_val == -1)     return(nb->IsInRing());
  else if (_val == 0) return(!(nb->IsInRing()));
  else return((signed)nb->MemberOfRingCount() == _val);
}

void OBRingExpr::Print(ostream &ofs)
{
	if (_val == -1) ofs << "IsInRing";
	else            ofs << "Ring=" << _val;
}

bool OBRecursExpr::Eval(OBNodeBase *nb)
{
	if (_vtest.empty()) //first match per molecule
	{
		_vtest.resize(nb->GetParent()->NumNodes()+1,false);
		_vmatch.resize(nb->GetParent()->NumNodes()+1,false);
	}
	else //if atom has already matched return value
		if (_vtest[nb->GetIdx()]) return(_vmatch[nb->GetIdx()]);

	_vtest[nb->GetIdx()] = true;
	_vmatch[nb->GetIdx()] = _sp->SeedMatch(nb);

	return(_vmatch[nb->GetIdx()]);
}

void OBRecursExpr::ResetRecurs()
{
	_vtest.clear();
	_vmatch.clear();
	_sp->ClearMatches();
}

OBRecursExpr::~OBRecursExpr()
{
	if (_sp) delete _sp;
}

}//namespace OpenBabel
