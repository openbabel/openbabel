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

namespace OpenEye
{

bool OEGraphBase::Match(OEGraphBase &g,bool singleMatch)
{
	SetFinishedMatch(false);
	SetSingleMatch(singleMatch);
	ClearMatches();
	g.SetVisitLock(true);
	g.ResetVisitFlags();

	OENodeBase *node;
	OENodeBase *seed = GetFirstSeed();
	vector<OENodeBase*>::iterator i;

	for (node = g.Begin(i);node;node = g.Next(i))
		if (!node->Visit && seed->Eval(node))
		{
			node->Visit = true;
			seed->SetMatch(node);
			Match(g,BgnMatch(),BgnMatch()->second.begin());
			seed->SetMatch((OENodeBase*)NULL);
			node->Visit = false;
			if (SingleMatch() && FinishedMatch()) break;
		}

	g.SetVisitLock(false);

	return(FinishedMatch());
}

bool OEGraphBase::Match(OEGraphBase &g,
							vector<pair<OENodeBase*,vector<OEEdgeBase*> > >::iterator i,
							vector<OEEdgeBase*>::iterator j)
{
	//bail if only one match has been requested
	if (SingleMatch() && FinishedMatch()) return(true);
	
	//full match completed
	if (i == EndMatch() || (j == i->second.end() && (i+1) == EndMatch())) 
	{
		SetFinishedMatch(true);
		OENodeBase *node;
		vector<OENodeBase*> vn;
		vector<OENodeBase*>::iterator i;
		for (node = Begin(i);node;node = Next(i))
			vn.push_back(node->GetMatch());
		PushBack(vn);
		
		return(true);
	}

	//handle next seed of disconnected pattern
	if (j == i->second.end()) 
	{
		i++;
		OENodeBase *node;
		OENodeBase *seed = i->first;
		vector<OENodeBase*>::iterator k;
		for (node = g.Begin(k);node;node = g.Next(k))
			if (!node->Visit && seed->Eval(node))
			{
				node->Visit = true;
				seed->SetMatch(node);
				Match(g,i,i->second.begin());
				seed->SetMatch((OENodeBase*)NULL);
				node->Visit = false;
				if (SingleMatch() && FinishedMatch()) break;
			}	
		return(true);
	}

	OEEdgeBase *edge = *j++;
	if (edge->IsClosure()) //check to see if matched atoms are bonded
	{
		if (edge->GetBgn()->GetMatch()->IsConnected(edge->GetEnd()->GetMatch()))
			Match(g,i,j);
	}
	else //bond hasn't been covered yet
	{
		OENodeBase *nbr;
		OENodeBase *curr = edge->GetBgn();
		OENodeBase *next = edge->GetEnd(); 
		OENodeBase *match = curr->GetMatch();
		vector<OEEdgeBase*>::iterator k;
	
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

OENodeBase *OEGraphBase::Begin(vector<OENodeBase*>::iterator &i)
{
	i = _vatom.begin();
	return((i != _vatom.end()) ? *i : NULL);
}

OENodeBase *OEGraphBase::Next(vector<OENodeBase*>::iterator &i)
{
	i++;
	return((i != _vatom.end()) ? *i : NULL);
}

OEEdgeBase *OEGraphBase::Begin(vector<OEEdgeBase*>::iterator &i) 
{
	i = _vbond.begin();
	return((i != _vbond.end()) ? *i : NULL);
}

OEEdgeBase *OEGraphBase::Next(vector<OEEdgeBase*>::iterator &i)
{
	i++;
	return((i != _vbond.end()) ? *i : NULL);
}

OENodeBase *OENodeBase::BeginNbr(vector<OEEdgeBase*>::iterator &i)
{
	i = _vbond.begin();

	if (i == _vbond.end()) return(NULL);
	return((this == (*i)->GetBgn()) ? (*i)->GetEnd() : (*i)->GetBgn());
}

OENodeBase *OENodeBase::NextNbr(vector<OEEdgeBase*>::iterator &i)
{
	i++;
	if (i == _vbond.end()) return(NULL);
	return((this == (*i)->GetBgn()) ? (*i)->GetEnd() : (*i)->GetBgn());
}

void OENodeBase::SetParent(OEGraphBase *p)
{
	_parent = p;
}

void OEEdgeBase::SetParent(OEGraphBase *p)
{
	_parent = p;
}

bool OENodeBase::IsConnected(OENodeBase *nb)
{
	vector<OEEdgeBase*>::iterator i;
	for (i = _vbond.begin();i != _vbond.end();i++)
		if (nb == (*i)->GetBgn() || nb == (*i)->GetEnd())
			return(true);

	return(false);
}

void OEGraphBase::ResetVisitFlags()
{
	OENodeBase *nb;
	vector<OENodeBase*>::iterator i;
	for (nb = Begin(i);nb;nb = Next(i)) nb->Visit = false;

	OEEdgeBase *eb;
	vector<OEEdgeBase*>::iterator j;
	for (eb = Begin(j);eb;eb = Next(j)) eb->Visit = false;
}

bool OEGraphBase::SetVisitLock(bool v)
{
	if (v && _vlock) return(false);
	_vlock = v;
	return(true);
}

}