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
//OESmartsPattern member functions
//

bool OESmartsPattern::SeedMatch(OENodeBase *nb)
//routine ONLY used for matching recursive SMARTS
{
	_done = false;
	_single = true;
	_vmatch.clear();
	OEGraphBase *g = nb->GetParent();
	vector<bool> visit(g->NumNodes()); 

	//save visit flags as they are being used further up in the stack
	int i;
	OENodeBase *node;
	vector<OENodeBase*>::iterator j;
	for (i=0,node = g->Begin(j);node;node = g->Next(j))
		visit[i++] = node->Visit;
	g->ResetVisitFlags();

	OENodeBase *seed = GetFirstSeed();

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

void OESmartsPattern::GetMatches(vector<vector<int> > &vm)
{
	vector<OENodeBase*>::iterator j;
	vector<vector<OENodeBase*> >::iterator i;
	for (i = _vmatch.begin();i != _vmatch.end();i++)
	{
		vector<int> vtmp;
		vtmp.reserve(i->size());
		for (j = i->begin();j != i->end();j++)
			vtmp.push_back((*j)->GetIdx());
		vm.push_back(vtmp);
	}
}

vector<vector<int> > OESmartsPattern::GetMapList()
{
	vector<vector<int> > vm;
	vector<OENodeBase*>::iterator j;
	vector<vector<OENodeBase*> >::iterator i;
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

vector<vector<int> > OESmartsPattern::GetUMapList()
{
	vector<vector<int> > vm; 
	vm.clear();

	if (!_vmatch.empty())
	{
		bool ok;
		OEBitVec bv;
		vector<OEBitVec> vbv;
		vector<OEBitVec>::iterator k;
		vector<OENodeBase*>::iterator j;
		vector<vector<OENodeBase*> >::iterator i;
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

void OESmartsPattern::PrepForMatch()
{
	_vm.clear();
	SetVisitLock(true);
	ResetVisitFlags();

	OEEdge *edge;
	OENode *node,*seed,*nbr;
	vector<OENodeBase*> curr,next;
	vector<OENodeBase*>::iterator i,j;
	vector<OEEdgeBase*>::iterator k;
	vector<OEEdgeBase*> vedge;

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
					for (nbr = ((OENode*)*j)->BeginNbr(k);nbr;nbr = ((OENode*)*j)->NextNbr(k))
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
			_vm.push_back(pair<OENodeBase*,vector<OEEdgeBase*> > (seed,vedge));
		}

	SetVisitLock(false);
}

bool OESmartsPattern::Init(const char *buf)
{
	OESmartsParser sp;
	return(sp.Parse(*this,buf));
}

bool OESmartsPattern::Init(string &s)
{
	OESmartsParser sp;
	return(sp.Parse(*this,s));
}

OENode *OESmartsPattern::NewNode(OEExprBase *xb)
{
	OENode *node = new OENode;
	node->SetExpr(xb);
	node->SetIdx(NumNodes());
	node->SetParent(this);
	_vatom.push_back(node);

	return(node);
}

OEEdge* OESmartsPattern::NewEdge(OENode *bgn,OENode *end,OEExprBase *expr)
{
	OEEdge *edge = new OEEdge(bgn,end,expr);
	edge->SetIdx(NumEdges());
	edge->SetParent(this);
	_vbond.push_back(edge);
	if (bgn) bgn->AddEdge(edge);
	if (end) end->AddEdge(edge);
	return(edge);
}

void OESmartsPattern::ClearMatches()
{
	_vmatch.clear(); //clear out previous matches

	OENode *node;
	vector<OENodeBase*>::iterator i;
	for (node = Begin(i);node;node = Next(i))
		node->ResetRecurs();
}

OENode *OESmartsPattern::Begin(vector<OENodeBase*>::iterator &i)
{
	i = _vatom.begin();
	return((i != _vatom.end()) ? (OENode*)*i : (OENode*)NULL);
}

OENode *OESmartsPattern::Next(vector<OENodeBase*>::iterator &i)
{
	i++;
	return((i != _vatom.end()) ? (OENode*)*i : (OENode*)NULL);
}

OEEdge *OESmartsPattern::Begin(vector<OEEdgeBase*>::iterator &i) 
{
	i = _vbond.begin();
	return((i != _vbond.end()) ? (OEEdge*)*i : (OEEdge*)NULL);
}

OEEdge *OESmartsPattern::Next(vector<OEEdgeBase*>::iterator &i)
{
	i++;
	return((i != _vbond.end()) ? (OEEdge*)*i : (OEEdge*)NULL);
}

OENode *OENode::BeginNbr(vector<OEEdgeBase*>::iterator &i)
{
	i = _vbond.begin();
	if (i == _vbond.end()) return((OENode*)NULL);
	return((this == ((OEEdge*)*i)->GetBgn()) ? (OENode*)((OEEdge*)*i)->GetEnd() : (OENode*)((OEEdge*)*i)->GetBgn());
}

OENode *OENode::NextNbr(vector<OEEdgeBase*>::iterator &i)
{
	i++;
	if (i == _vbond.end()) return((OENode*)NULL);
	return((this == ((OEEdge*)*i)->GetBgn()) ? (OENode*)((OEEdge*)*i)->GetEnd() : (OENode*)((OEEdge*)*i)->GetBgn());
}

OESmartsPattern::OESmartsPattern()
{
  _done = false;
  _single = false;
  _smarts = "";
  _vmatch.clear();
  _vm.clear();
}

OESmartsPattern::OESmartsPattern(const OESmartsPattern &src)
{
  string s = src.GetSMARTS();
  Init(s);
  _done = src.FinishedMatch();
  _single = src.SingleMatch();
  src.GetMatches(_vmatch);
}

OESmartsPattern::~OESmartsPattern()
{
	vector<OENodeBase*>::iterator i;
	for (i = _vatom.begin();i != _vatom.end();i++) delete *i;

	vector<OEEdgeBase*>::iterator j;
	for (j = _vbond.begin();j != _vbond.end();j++) delete *j;
}

//
//Expression member functions
//

int OEAndExpr::GetFormalCharge()
{
	int chrg = _lft->GetFormalCharge();
	if (chrg != 0) return(chrg);
	return(_rgt->GetFormalCharge());
}

unsigned int OEAndExpr::GetAtomicNum()
{
	unsigned int ele = _lft->GetAtomicNum();
	if (ele != 0) return(ele);
	return(_rgt->GetAtomicNum());
}

unsigned int OEAndExpr::GetBO()
{
	unsigned int bo = _lft->GetBO();
	if (bo != 0) return(bo);
	return(_rgt->GetBO());
}

bool OEAndExpr::IsAromatic()
{
	if (_lft->IsAromatic()) return(true);
	return(_rgt->IsAromatic());
}

void OEAndExpr::Print(ostream &ofs)
{
	ofs << "((";
	_lft->Print(ofs);
	ofs << ")&(";
	_rgt->Print(ofs);
	ofs << "))";
}

void OEAndExpr::ResetRecurs()
{
	_lft->ResetRecurs();
	_rgt->ResetRecurs();
}

int OEOrExpr::GetFormalCharge()
{
	int chrg = _lft->GetFormalCharge();
	if (chrg != 0) return(chrg);
	return(_rgt->GetFormalCharge());
}
unsigned int OEOrExpr::GetAtomicNum()
{
	unsigned int ele = _lft->GetAtomicNum();
	if (ele != 0) return(ele);
	return(_rgt->GetAtomicNum());
}

unsigned int OEOrExpr::GetBO()
{
	unsigned int bo = _lft->GetBO();
	if (bo != 0) return(bo);
	return(_rgt->GetBO());
}
bool OEOrExpr::IsAromatic() 
{
	return(_lft->IsAromatic() || _rgt->IsAromatic());
}

void OEOrExpr::Print(ostream &ofs)  
{
	ofs << "((";
	_lft->Print(ofs);
	ofs << "),(";
	_rgt->Print(ofs);
	ofs << "))";
}

void OEOrExpr::ResetRecurs()
{
	_lft->ResetRecurs();
	_rgt->ResetRecurs();
}

void OENotExpr::Print(ostream &ofs)
{
	ofs << "!(";
	_expr->Print(ofs);
	ofs << ")";
}

void OENotExpr::ResetRecurs()
{
	_expr->ResetRecurs();
}

void OEConstExpr::Print(ostream &ofs)
{
	ofs << "(TRUE)";
}

void OEElementExpr::Print(ostream &ofs)
{
	ofs << "#" << _ele;
}

void OEAromaticExpr::Print(ostream &ofs)  
{
	if (_val) ofs << "Aromatic";
	else      ofs << "Aliphatic";
}

bool OERingExpr::Eval(OENodeBase *nb) 
{
  if (_val == -1)     return(nb->IsInRing());
  else if (_val == 0) return(!(nb->IsInRing()));
  else return((signed)nb->MemberOfRingCount() == _val);
}

void OERingExpr::Print(ostream &ofs)
{
	if (_val == -1) ofs << "IsInRing";
	else            ofs << "Ring=" << _val;
}

bool OERecursExpr::Eval(OENodeBase *nb)
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

void OERecursExpr::ResetRecurs()
{
	_vtest.clear();
	_vmatch.clear();
	_sp->ClearMatches();
}

OERecursExpr::~OERecursExpr()
{
	if (_sp) delete _sp;
}

}//namespace OpenBabel
