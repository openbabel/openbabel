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

#ifndef SMARTS_H
#define SMARTS_H

#ifdef WIN32
#pragma warning (disable: 4786) // warning: long & complicated stl warning
#pragma warning (disable: 4305) // warning: truncation from 'const double' to 'const float'
#endif

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "base.h"

using namespace std;

namespace OpenEye
{

class OEEdge;
class OENode;
class OEGraph;
class OEExprBase;
class OESmartsPattern;

//
// Base expression class from which all other expressions are derived
//

class OEExprBase
{
public:
	OEExprBase()                           {}
	virtual ~OEExprBase()                  {}
	virtual int  GetFormalCharge()         {return(0);}
	virtual void ResetRecurs()             {}
	virtual void Print(ostream&)           {}
	virtual bool Eval(OENodeBase*)         {return(false);}
	virtual bool Eval(OEEdgeBase*)         {return(false);}
	virtual bool IsAromatic()              {return(false);}
	virtual unsigned int GetAtomicNum()    {return(0);}
	virtual unsigned int GetBO()           {return(0);}
};

//
// Logical expressions for both nodes and edges
//

class OEAndExpr : public OEExprBase
{
protected:
	OEExprBase *_lft;
	OEExprBase *_rgt;
public:
	OEAndExpr(OEExprBase *lft,OEExprBase *rgt) {_lft = lft;_rgt = rgt;}
	~OEAndExpr()
	{
		if (_lft) delete _lft;
		if (_rgt) delete _rgt;
	}
	bool Eval(OENodeBase *nb) {return((_lft->Eval(nb) && _rgt->Eval(nb)));}
	bool Eval(OEEdgeBase *eb) {return((_lft->Eval(eb) && _rgt->Eval(eb)));}
	int          GetFormalCharge();
	void         Print(ostream&);
	void         ResetRecurs();
	bool         IsAromatic();
	unsigned int GetAtomicNum();
	unsigned int GetBO();
};

class OEOrExpr : public OEExprBase
{
protected:
	OEExprBase *_lft;
	OEExprBase *_rgt;
public:
	OEOrExpr(OEExprBase *lft,OEExprBase *rgt) {_lft = lft; _rgt = rgt;}
	~OEOrExpr()
	{
		if (_lft) delete _lft;
		if (_rgt) delete _rgt;
	}
	bool Eval(OENodeBase *nb) {return((_lft->Eval(nb) || _rgt->Eval(nb)));}
	bool Eval(OEEdgeBase *eb) {return((_lft->Eval(eb) || _rgt->Eval(eb)));}
	int          GetFormalCharge();
	void         Print(ostream&);
	void         ResetRecurs();
	bool         IsAromatic();
	unsigned int GetAtomicNum();
	unsigned int GetBO();
};

class OENotExpr : public OEExprBase
{
protected:
	OEExprBase *_expr;
public:
	OENotExpr(OEExprBase *expr) {_expr = expr;}
	~OENotExpr()                {if (_expr) delete _expr;}
	void Print(ostream&);
	void ResetRecurs();
	bool Eval(OENodeBase *nb)   {return(!_expr->Eval(nb));}
	bool Eval(OEEdgeBase *eb)   {return(!_expr->Eval(eb));}
	bool IsAromatic()           {return(!_expr->IsAromatic());}
};

class OEConstExpr : public OEExprBase
{
protected:
public:
	OEConstExpr()          {}
	~OEConstExpr()         {}
	void Print(ostream&);
	bool Eval(OENodeBase *nb) {return(nb->GetAtomicNum() != 1);}
	//bool Eval(OENodeBase*) {return(true);}
	bool Eval(OEEdgeBase*) {return(true);}
	bool IsAromatic()      {return(true);}
};

//
// Node expression class definitions
//

class OEElementExpr : public OEExprBase
{
protected:
	unsigned int _ele;
public:
	OEElementExpr(int ele) {_ele = ele;}
	~OEElementExpr() {}
	void Print(ostream &ofs);
	bool Eval(OENodeBase *nb)    {return(_ele == nb->GetAtomicNum());}
	bool IsAromatic()            {return(true);}
	unsigned int GetAtomicNum()  {return(_ele);}
};

class OEAromElemExpr : public OEExprBase
{
protected:
	unsigned int  _ele;
	bool _aro;
public:
	OEAromElemExpr(int ele,bool aro) {_ele = ele; _aro = aro;}
	~OEAromElemExpr() {}
	bool Eval(OENodeBase *nb)
	{
	  return(_ele == nb->GetAtomicNum() && _aro == nb->IsAromatic());
	}
	bool         IsAromatic()   {return(_aro);}
	unsigned int GetAtomicNum() {return(_ele);}
};

class OEAromaticExpr : public OEExprBase //both nodes and edges
{
protected:
	bool _val;
public:
	OEAromaticExpr(bool val) {_val = val;}
	~OEAromaticExpr() {}
	void Print(ostream&);
	bool Eval(OENodeBase *nb) {return(_val == nb->IsAromatic());}
	bool Eval(OEEdgeBase *eb) {return(_val == eb->IsAromatic());}
	bool IsAromatic()         {return(_val);}
};

class OEMassExpr : public OEExprBase
{
protected:
	int _val;
public:
	OEMassExpr(int val) {_val = val;}
	~OEMassExpr(){}
	bool Eval(OENodeBase*) {return(false);}
};

class OEHCountExpr : public OEExprBase
{
protected:
	unsigned int _val;
public:
	OEHCountExpr(int val) {_val = val;}
	~OEHCountExpr() {}
	bool Eval(OENodeBase *nb) 
	  {
            if (nb->ExplicitHydrogenCount() > nb->ImplicitHydrogenCount())
              return(_val == nb->ExplicitHydrogenCount());
            else
              return(_val == nb->ImplicitHydrogenCount());
	  }
};

class OENegChargeExpr : public OEExprBase
{
protected:
	int _val;
public:
	OENegChargeExpr(int val) {_val = -1*abs(val);}
	~OENegChargeExpr()        {}
	bool Eval(OENodeBase *nb) {return(nb->GetFormalCharge() == _val);}
	int  GetFormalCharge()    {return(_val);}
};

class OEPosChargeExpr : public OEExprBase
{
protected:
	int _val;
public:
	OEPosChargeExpr(int val) {_val = abs(val);}
	~OEPosChargeExpr()        {}
	bool Eval(OENodeBase *nb) {return(nb->GetFormalCharge() == _val);}
	int  GetFormalCharge()    {return(_val);}
};

class OEConnectExpr : public OEExprBase
{
protected:
	unsigned int _val;
public:
	OEConnectExpr(int val)    {_val = val;}
	~OEConnectExpr()          {}
	bool Eval(OENodeBase *nb) {return(nb->GetImplicitValence() == _val);}
};

class OEDegreeExpr : public OEExprBase
{
protected:
	unsigned int _val;
public:
	OEDegreeExpr(int val)     {_val = val;}
	~OEDegreeExpr()           {}
	bool Eval(OENodeBase *nb) {return(nb->GetHvyValence() == _val);}
};

class OEImplicitExpr : public OEExprBase
{
protected:
	unsigned int _val;
public:
	OEImplicitExpr(int val)   {_val = val;}
	~OEImplicitExpr()         {}
	bool Eval(OENodeBase *nb) {return(nb->ImplicitHydrogenCount() == _val);}
};

class OERingExpr : public OEExprBase //both nodes and edges
{
protected:
	int _val;
public:
	OERingExpr(int val = -1)   {_val = val;}
	~OERingExpr()              {}
	void Print(ostream&);
	bool Eval(OENodeBase *); 
	bool Eval(OEEdgeBase *eb) {return(eb->IsInRing());}
};

class OESizeExpr : public OEExprBase
{
protected:
	int _val;
public:
	OESizeExpr(int val = 0) {_val = val;}
	~OESizeExpr()           {}
	bool Eval(OENodeBase *nb) 
	{
	  return((_val) ? nb->IsInRingSize(_val) : !nb->IsInRing());
	}
};

class OEValenceExpr : public OEExprBase
{
protected:
	unsigned int _val;
public:
	OEValenceExpr(int val)    {_val = val;}
	~OEValenceExpr()          {}
	bool Eval(OENodeBase *nb) {return(_val == nb->KBOSum());}
};

class OEHybExpr : public OEExprBase
{
protected:
	unsigned int _val;
public:
	OEHybExpr(int val) {_val = val;}
	~OEHybExpr()       {}
	bool Eval(OENodeBase *nb)  {return(_val == nb->GetHyb());}
};

class OERecursExpr : public OEExprBase
{
protected:
	vector<bool>     _vtest;
	vector<bool>     _vmatch;
	OESmartsPattern *_sp;
public:
	OERecursExpr()                    {_sp = NULL;}
	OERecursExpr(OESmartsPattern *sp) {_sp = sp;}
	~OERecursExpr();
	bool Eval(OENodeBase*);
	void ResetRecurs();
};

//
//Edge expression class definitions
//

class OEDefaultEdgeExpr : public OEExprBase
{
protected:
public:
	OEDefaultEdgeExpr()       {}
	~OEDefaultEdgeExpr()      {}
	bool Eval(OEEdgeBase *eb)  {return((eb->GetBO() == 1 || eb->IsAromatic()));}
	unsigned int GetBO() {return(1);}
};

class OESingleExpr : public OEExprBase
{
protected:
public:
	OESingleExpr()  {}
	~OESingleExpr() {}
	void Print(ostream &ofs)  {ofs << "Single";}
	bool Eval(OEEdgeBase *eb) 
	  {
	    return(eb->GetBO() == 1 && !eb->IsAromatic());
	  }
	unsigned int GetBO()      {return(1);}
};

class OEDoubleExpr : public OEExprBase
{
protected:
public:
	OEDoubleExpr()  {}
	~OEDoubleExpr() {}
	void Print(ostream &ofs)  {ofs << "Double";}
	bool Eval(OEEdgeBase *eb) 
	  {
	    return(eb->GetBO() == 2 && !eb->IsAromatic());
	  }
	unsigned int GetBO()      {return(2);}
};

class OETripleExpr : public OEExprBase
{
protected:
public:
	OETripleExpr()  {}
	~OETripleExpr() {}
	bool Eval(OEEdgeBase *eb) {return(eb->GetBO() == 3);}
	unsigned int GetBO()      {return(3);}
};

class OEUpExpr : public OEExprBase
{
protected:
public:
	OEUpExpr() {}
	~OEUpExpr() {}
	bool Eval(OEEdgeBase *eb) {return(false);}
};

class OEDownExpr : public OEExprBase
{
protected:
public:
	OEDownExpr()  {}
	~OEDownExpr() {}
	bool Eval(OEEdgeBase *eb) {return(false);}
};

class OEUpUnspecExpr : public OEExprBase
{
protected:
public:
	OEUpUnspecExpr() {}
	~OEUpUnspecExpr() {}
	bool Eval(OEEdgeBase *eb) {return(false);}
};

class OEDownUnspecExpr : public OEExprBase
{
protected:
public:
	OEDownUnspecExpr()  {}
	~OEDownUnspecExpr() {}
	bool Eval(OEEdgeBase *eb) {return(false);}
};

//
//OESmartsPattern class definition
//

class OEEdge : public OEEdgeBase
{
protected:
	OEExprBase *_expr;
	bool        _closure;
public:
	OEEdge(OENodeBase *bgn,OENodeBase *end,OEExprBase *expr) 
	{
		_bgn = bgn; 
		_end = end; 
		_expr = expr;
		_closure = false;
	}
	void ReplaceExpr(OEExprBase *expr)
	{
		if (_expr) delete _expr;
		_expr = expr;
	}
	void        SetClosure()        {_closure = true;}
	bool        IsClosure()         {return(_closure);}
	bool        Eval(OEEdgeBase *e) {return(_expr->Eval(e));}
	unsigned int GetBO()            {return(_expr->GetBO());}
	OEExprBase *GetExpr()           {return(_expr);}
	OENodeBase *GetBgn()            {return(_bgn);}
	OENodeBase *GetEnd()            {return(_end);}
};

class OENode : public OENodeBase
{
	int        _stereo;
	int        _vb;
	OENodeBase *_match;
	OEExprBase *_expr;
public:
	OENode() {_stereo=0; _match = NULL;}
	virtual ~OENode() {}
	int     GetFormalCharge()         const  {return(_expr->GetFormalCharge());}
	void    SetExpr(OEExprBase *expr)        {_expr = expr;}
	void    SetMatch(OENodeBase *m)          {_match = m;}
	void    SetVectorBinding(int vb)         {_vb = vb;}
	void    ResetRecurs()                    {_expr->ResetRecurs();}
	void    Print(ostream &ofs)              {_expr->Print(ofs);}
	bool	Eval(OENodeBase *n)       const  {return(_expr->Eval(n));}
	bool         IsAromatic()                {return(_expr->IsAromatic());}
	unsigned int GetVectorBinding()          {return(_vb);}
	unsigned int GetAtomicNum()       const  {return(_expr->GetAtomicNum());}
	OENodeBase *GetMatch()                   {return(_match);}
	OENode *BeginNbr(vector<OEEdge*>::iterator &);
	OENode *NextNbr(vector<OEEdge*>::iterator &);
};

class OESmartsPattern : public OEGraphBase
{
protected:
	bool                                            _done;
	bool                                            _single;
	string                                          _smarts;
	vector<vector<OENodeBase*> >                    _vmatch;
	vector<pair<OENodeBase*,vector<OEEdgeBase*> > > _vm;
public:
	OESmartsPattern();
	OESmartsPattern(const OESmartsPattern &);
	~OESmartsPattern();


	//initialization and graph modification methods
	bool    Init(const char *);
	bool    Init(string &);
	OENode *NewNode(OEExprBase*);
	OEEdge *NewEdge(OENode*,OENode*,OEExprBase*);
	
	//iterator methods
	OENode* Begin(vector<OENode*>::iterator &); 
	OENode* Next(vector<OENode*>::iterator &);
	OEEdge* Begin(vector<OEEdge*>::iterator &);
	OEEdge* Next(vector<OEEdge*>::iterator &);

	//routines that perform substructure search
	bool        SingleMatch()                    const{return(_single);}
	bool        FinishedMatch()                  const {return(_done);}
	void        SetSingleMatch(bool state)       {_single = state;}
	void        SetFinishedMatch(bool state)     {_done = state;}
	void        PrepForMatch();
	void        PushBack(vector<OENodeBase*> &v) {_vmatch.push_back(v);}
	void        ClearMatches();
	OENodeBase *GetFirstSeed()                   {return(_vm.begin()->first);}
	bool        SeedMatch(OENodeBase*);
	bool        RestrictedMatch(OEGraphBase&,vector<pair<int,int> >&,bool) 
	{
		cerr << "need to implement OESmartsPattern::RestrictedMatch()" << endl;
		exit(0);
		return(false);
	}
	vector<pair<OENodeBase*,vector<OEEdgeBase*> > >::iterator BgnMatch() {return(_vm.begin());}
	vector<pair<OENodeBase*,vector<OEEdgeBase*> > >::iterator EndMatch() {return(_vm.end());}

	//data mod methods
	void SetSMARTS(const char *buf) {_smarts = buf;}

	//methods to return match data
	void GetMatches(vector<vector<OENodeBase*> > &vm) const {vm = _vmatch;}
	void GetMatches(vector<vector<int> > &vm);
	vector<vector<int> > GetMapList();
	vector<vector<int> > GetUMapList();

	//misc data retrieval methods
	int           GetFormalCharge(int idx)  {return(_vatom[idx]->GetFormalCharge());}
	bool          IsValid()                 {return(!_vatom.empty());}
	bool          Empty()                   {return(_vatom.empty());}
	string        GetSMARTS()               const {return(_smarts);}
	unsigned int  NumMatches()              {return((_vmatch.empty()) ? 0 : _vmatch.size());}
	unsigned int  NumAtoms()                {return(NumNodes());}
	unsigned int  NumBonds()                {return(NumEdges());}
	unsigned int  GetVectorBinding(int idx) {return(((OENode*)_vatom[idx])->GetVectorBinding());}
	unsigned int  GetAtomicNum(int idx)     const {return(_vatom[idx]->GetAtomicNum());}

	//stream output methods
	void WriteMapList(ostream&);
};


} //namespace OpenEye

#endif //SMARTS_H



