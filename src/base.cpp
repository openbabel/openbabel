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

#include "base.h"

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenBabel {

bool OBGraphBase::Match(OBGraphBase &g,bool singleMatch)
{
	SetFinishedMatch(false);
	SetSingleMatch(singleMatch);
	ClearMatches();
	g.SetVisitLock(true);
	g.ResetVisitFlags();

	OBNodeBase *node;
	OBNodeBase *seed = GetFirstSeed();
	vector<OBNodeBase*>::iterator i;

	for (node = g.Begin(i);node;node = g.Next(i))
		if (!node->Visit && seed->Eval(node))
		{
			node->Visit = true;
			seed->SetMatch(node);
			Match(g,BgnMatch(),BgnMatch()->second.begin());
			seed->SetMatch((OBNodeBase*)NULL);
			node->Visit = false;
			if (SingleMatch() && FinishedMatch()) break;
		}

	g.SetVisitLock(false);

	return(FinishedMatch());
}

bool OBGraphBase::Match(OBGraphBase &g,
							vector<pair<OBNodeBase*,vector<OBEdgeBase*> > >::iterator i,
							vector<OBEdgeBase*>::iterator j)
{
	//bail if only one match has been requested
	if (SingleMatch() && FinishedMatch()) return(true);
	
	//full match completed
	if (i == EndMatch() || (j == i->second.end() && (i+1) == EndMatch())) 
	{
		SetFinishedMatch(true);
		OBNodeBase *node;
		vector<OBNodeBase*> vn;
		vector<OBNodeBase*>::iterator i;
		for (node = Begin(i);node;node = Next(i))
			vn.push_back(node->GetMatch());
		PushBack(vn);
		
		return(true);
	}

	//handle next seed of disconnected pattern
	if (j == i->second.end()) 
	{
		i++;
		OBNodeBase *node;
		OBNodeBase *seed = i->first;
		vector<OBNodeBase*>::iterator k;
		for (node = g.Begin(k);node;node = g.Next(k))
			if (!node->Visit && seed->Eval(node))
			{
				node->Visit = true;
				seed->SetMatch(node);
				Match(g,i,i->second.begin());
				seed->SetMatch((OBNodeBase*)NULL);
				node->Visit = false;
				if (SingleMatch() && FinishedMatch()) break;
			}	
		return(true);
	}

	OBEdgeBase *edge = *j++;
	if (edge->IsClosure()) //check to see if matched atoms are bonded
	{
		if (edge->GetBgn()->GetMatch()->IsConnected(edge->GetEnd()->GetMatch()))
			Match(g,i,j);
	}
	else //bond hasn't been covered yet
	{
		OBNodeBase *nbr;
		OBNodeBase *curr = edge->GetBgn();
		OBNodeBase *next = edge->GetEnd(); 
		OBNodeBase *match = curr->GetMatch();
		vector<OBEdgeBase*>::iterator k;
	
		for (nbr = match->BeginNbr(k);nbr;nbr = match->NextNbr(k))
			if (!nbr->Visit && next->Eval(nbr) && edge->Eval(*k))
			{
				nbr->Visit = true;
				next->SetMatch(nbr);
				Match(g,i,j);
				next->SetMatch(NULL);
				nbr->Visit = false;
			}
	}

	return(false);
}

OBNodeBase *OBGraphBase::Begin(vector<OBNodeBase*>::iterator &i)
{
	i = _vatom.begin();
	return((i != _vatom.end()) ? *i : NULL);
}

OBNodeBase *OBGraphBase::Next(vector<OBNodeBase*>::iterator &i)
{
	i++;
	return((i != _vatom.end()) ? *i : NULL);
}

OBEdgeBase *OBGraphBase::Begin(vector<OBEdgeBase*>::iterator &i) 
{
	i = _vbond.begin();
	return((i != _vbond.end()) ? *i : NULL);
}

OBEdgeBase *OBGraphBase::Next(vector<OBEdgeBase*>::iterator &i)
{
	i++;
	return((i != _vbond.end()) ? *i : NULL);
}

OBNodeBase *OBNodeBase::BeginNbr(vector<OBEdgeBase*>::iterator &i)
{
	i = _vbond.begin();

	if (i == _vbond.end()) return(NULL);
	return((this == (*i)->GetBgn()) ? (*i)->GetEnd() : (*i)->GetBgn());
}

OBNodeBase *OBNodeBase::NextNbr(vector<OBEdgeBase*>::iterator &i)
{
	i++;
	if (i == _vbond.end()) return(NULL);
	return((this == (*i)->GetBgn()) ? (*i)->GetEnd() : (*i)->GetBgn());
}

void OBNodeBase::SetParent(OBGraphBase *p)
{
	_parent = p;
}

void OBEdgeBase::SetParent(OBGraphBase *p)
{
	_parent = p;
}

bool OBNodeBase::IsConnected(OBNodeBase *nb)
{
	vector<OBEdgeBase*>::iterator i;
	for (i = _vbond.begin();i != _vbond.end();i++)
		if (nb == (*i)->GetBgn() || nb == (*i)->GetEnd())
			return(true);

	return(false);
}

void OBGraphBase::ResetVisitFlags()
{
	OBNodeBase *nb;
	vector<OBNodeBase*>::iterator i;
	for (nb = Begin(i);nb;nb = Next(i)) nb->Visit = false;

	OBEdgeBase *eb;
	vector<OBEdgeBase*>::iterator j;
	for (eb = Begin(j);eb;eb = Next(j)) eb->Visit = false;
}

bool OBGraphBase::SetVisitLock(bool v)
{
	if (v && _vlock) return(false);
	_vlock = v;
	return(true);
}

}
