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

namespace OpenBabel
{

class OBEdge;
class OBNode;
class OBGraph;
class OBExprBase;
class OBSmartsPattern;

//
// Base expression class from which all other expressions are derived
//

class OBExprBase
{
public:
	OBExprBase()                           {}
	virtual ~OBExprBase()                  {}
	virtual int  GetFormalCharge()         {return(0);}
	virtual void ResetRecurs()             {}
	virtual void Print(std::ostream&)           {}
	virtual bool Eval(OBNodeBase*)         {return(false);}
	virtual bool Eval(OBEdgeBase*)         {return(false);}
	virtual bool IsAromatic()              {return(false);}
	virtual unsigned int GetAtomicNum()    {return(0);}
	virtual unsigned int GetBO()           {return(0);}
};

//
// Logical expressions for both nodes and edges
//

class OBAndExpr : public OBExprBase
{
protected:
	OBExprBase *_lft;
	OBExprBase *_rgt;
public:
	OBAndExpr(OBExprBase *lft,OBExprBase *rgt) {_lft = lft;_rgt = rgt;}
	~OBAndExpr()
	{
		if (_lft) delete _lft;
		if (_rgt) delete _rgt;
	}
	bool Eval(OBNodeBase *nb) {return((_lft->Eval(nb) && _rgt->Eval(nb)));}
	bool Eval(OBEdgeBase *eb) {return((_lft->Eval(eb) && _rgt->Eval(eb)));}
	int          GetFormalCharge();
	void         Print(std::ostream&);
	void         ResetRecurs();
	bool         IsAromatic();
	unsigned int GetAtomicNum();
	unsigned int GetBO();
};

class OBOrExpr : public OBExprBase
{
protected:
	OBExprBase *_lft;
	OBExprBase *_rgt;
public:
	OBOrExpr(OBExprBase *lft,OBExprBase *rgt) {_lft = lft; _rgt = rgt;}
	~OBOrExpr()
	{
		if (_lft) delete _lft;
		if (_rgt) delete _rgt;
	}
	bool Eval(OBNodeBase *nb) {return((_lft->Eval(nb) || _rgt->Eval(nb)));}
	bool Eval(OBEdgeBase *eb) {return((_lft->Eval(eb) || _rgt->Eval(eb)));}
	int          GetFormalCharge();
	void         Print(std::ostream&);
	void         ResetRecurs();
	bool         IsAromatic();
	unsigned int GetAtomicNum();
	unsigned int GetBO();
};

class OBNotExpr : public OBExprBase
{
protected:
	OBExprBase *_expr;
public:
	OBNotExpr(OBExprBase *expr) {_expr = expr;}
	~OBNotExpr()                {if (_expr) delete _expr;}
	void Print(std::ostream&);
	void ResetRecurs();
	bool Eval(OBNodeBase *nb)   {return(!_expr->Eval(nb));}
	bool Eval(OBEdgeBase *eb)   {return(!_expr->Eval(eb));}
	bool IsAromatic()           {return(!_expr->IsAromatic());}
};

class OBConstExpr : public OBExprBase
{
protected:
public:
	OBConstExpr()          {}
	~OBConstExpr()         {}
	void Print(std::ostream&);
	bool Eval(OBNodeBase *nb) {return(nb->GetAtomicNum() != 1);}
	//bool Eval(OBNodeBase*) {return(true);}
	bool Eval(OBEdgeBase*) {return(true);}
	bool IsAromatic()      {return(true);}
};

//
// Node expression class definitions
//

class OBElementExpr : public OBExprBase
{
protected:
	unsigned int _ele;
public:
	OBElementExpr(int ele) {_ele = ele;}
	~OBElementExpr() {}
	void Print(std::ostream &ofs);
	bool Eval(OBNodeBase *nb)    {return(_ele == nb->GetAtomicNum());}
	bool IsAromatic()            {return(true);}
	unsigned int GetAtomicNum()  {return(_ele);}
};

class OBAromElemExpr : public OBExprBase
{
protected:
	unsigned int  _ele;
	bool _aro;
public:
	OBAromElemExpr(int ele,bool aro) {_ele = ele; _aro = aro;}
	~OBAromElemExpr() {}
	bool Eval(OBNodeBase *nb)
	{
	  return(_ele == nb->GetAtomicNum() && _aro == nb->IsAromatic());
	}
	bool         IsAromatic()   {return(_aro);}
	unsigned int GetAtomicNum() {return(_ele);}
};

class OBAromaticExpr : public OBExprBase //both nodes and edges
{
protected:
	bool _val;
public:
	OBAromaticExpr(bool val) {_val = val;}
	~OBAromaticExpr() {}
	void Print(std::ostream&);
	bool Eval(OBNodeBase *nb) {return(_val == nb->IsAromatic());}
	bool Eval(OBEdgeBase *eb) {return(_val == eb->IsAromatic());}
	bool IsAromatic()         {return(_val);}
};

class OBMassExpr : public OBExprBase
{
protected:
	int _val;
public:
	OBMassExpr(int val) {_val = val;}
	~OBMassExpr(){}
	bool Eval(OBNodeBase*) {return(false);}
};

class OBHCountExpr : public OBExprBase
{
protected:
	unsigned int _val;
public:
	OBHCountExpr(int val) {_val = val;}
	~OBHCountExpr() {}
	bool Eval(OBNodeBase *nb) 
	  {
            if (nb->ExplicitHydrogenCount() > nb->ImplicitHydrogenCount())
              return(_val == nb->ExplicitHydrogenCount());
            else
              return(_val == nb->ImplicitHydrogenCount());
	  }
};

class OBNegChargeExpr : public OBExprBase
{
protected:
	int _val;
public:
	OBNegChargeExpr(int val) {_val = -1*abs(val);}
	~OBNegChargeExpr()        {}
	bool Eval(OBNodeBase *nb) {return(nb->GetFormalCharge() == _val);}
	int  GetFormalCharge()    {return(_val);}
};

class OBPosChargeExpr : public OBExprBase
{
protected:
	int _val;
public:
	OBPosChargeExpr(int val) {_val = abs(val);}
	~OBPosChargeExpr()        {}
	bool Eval(OBNodeBase *nb) {return(nb->GetFormalCharge() == _val);}
	int  GetFormalCharge()    {return(_val);}
};

class OBConnectExpr : public OBExprBase
{
protected:
	unsigned int _val;
public:
	OBConnectExpr(int val)    {_val = val;}
	~OBConnectExpr()          {}
	bool Eval(OBNodeBase *nb) {return(nb->GetImplicitValence() == _val);}
};

class OBDegreeExpr : public OBExprBase
{
protected:
	unsigned int _val;
public:
	OBDegreeExpr(int val)     {_val = val;}
	~OBDegreeExpr()           {}
	bool Eval(OBNodeBase *nb) {return(nb->GetHvyValence() == _val);}
};

class OBImplicitExpr : public OBExprBase
{
protected:
	unsigned int _val;
public:
	OBImplicitExpr(int val)   {_val = val;}
	~OBImplicitExpr()         {}
	bool Eval(OBNodeBase *nb) {return(nb->ImplicitHydrogenCount() == _val);}
};

class OBRingExpr : public OBExprBase //both nodes and edges
{
protected:
	int _val;
public:
	OBRingExpr(int val = -1)   {_val = val;}
	~OBRingExpr()              {}
	void Print(std::ostream&);
	bool Eval(OBNodeBase *); 
	bool Eval(OBEdgeBase *eb) {return(eb->IsInRing());}
};

class OBSizeExpr : public OBExprBase
{
protected:
	int _val;
public:
	OBSizeExpr(int val = 0) {_val = val;}
	~OBSizeExpr()           {}
	bool Eval(OBNodeBase *nb) 
	{
	  return((_val) ? nb->IsInRingSize(_val) : !nb->IsInRing());
	}
};

class OBValenceExpr : public OBExprBase
{
protected:
	unsigned int _val;
public:
	OBValenceExpr(int val)    {_val = val;}
	~OBValenceExpr()          {}
	bool Eval(OBNodeBase *nb) {return(_val == nb->KBOSum());}
};

class OBHybExpr : public OBExprBase
{
protected:
	unsigned int _val;
public:
	OBHybExpr(int val) {_val = val;}
	~OBHybExpr()       {}
	bool Eval(OBNodeBase *nb)  {return(_val == nb->GetHyb());}
};

class OBRecursExpr : public OBExprBase
{
protected:
	std::vector<bool>     _vtest;
	std::vector<bool>     _vmatch;
	OBSmartsPattern *_sp;
public:
	OBRecursExpr()                    {_sp = NULL;}
	OBRecursExpr(OBSmartsPattern *sp) {_sp = sp;}
	~OBRecursExpr();
	bool Eval(OBNodeBase*);
	void ResetRecurs();
};

//
//Edge expression class definitions
//

class OBDefaultEdgeExpr : public OBExprBase
{
protected:
public:
	OBDefaultEdgeExpr()       {}
	~OBDefaultEdgeExpr()      {}
	bool Eval(OBEdgeBase *eb)  {return((eb->GetBO() == 1 || eb->IsAromatic()));}
	unsigned int GetBO() {return(1);}
};

class OBSingleExpr : public OBExprBase
{
protected:
public:
	OBSingleExpr()  {}
	~OBSingleExpr() {}
	void Print(std::ostream &ofs)  {ofs << "Single";}
	bool Eval(OBEdgeBase *eb) 
	  {
	    return(eb->GetBO() == 1 && !eb->IsAromatic());
	  }
	unsigned int GetBO()      {return(1);}
};

class OBDoubleExpr : public OBExprBase
{
protected:
public:
	OBDoubleExpr()  {}
	~OBDoubleExpr() {}
	void Print(std::ostream &ofs)  {ofs << "Double";}
	bool Eval(OBEdgeBase *eb) 
	  {
	    return(eb->GetBO() == 2 && !eb->IsAromatic());
	  }
	unsigned int GetBO()      {return(2);}
};

class OBTripleExpr : public OBExprBase
{
protected:
public:
	OBTripleExpr()  {}
	~OBTripleExpr() {}
	bool Eval(OBEdgeBase *eb) {return(eb->GetBO() == 3);}
	unsigned int GetBO()      {return(3);}
};

class OBUpExpr : public OBExprBase
{
protected:
public:
	OBUpExpr() {}
	~OBUpExpr() {}
	bool Eval(OBEdgeBase *eb) {return(false);}
};

class OBDownExpr : public OBExprBase
{
protected:
public:
	OBDownExpr()  {}
	~OBDownExpr() {}
	bool Eval(OBEdgeBase *eb) {return(false);}
};

class OBUpUnspecExpr : public OBExprBase
{
protected:
public:
	OBUpUnspecExpr() {}
	~OBUpUnspecExpr() {}
	bool Eval(OBEdgeBase *eb) {return(false);}
};

class OBDownUnspecExpr : public OBExprBase
{
protected:
public:
	OBDownUnspecExpr()  {}
	~OBDownUnspecExpr() {}
	bool Eval(OBEdgeBase *eb) {return(false);}
};

//
//OBSmartsPattern class definition
//

class OBEdge : public OBEdgeBase
{
protected:
	OBExprBase *_expr;
	bool        _closure;
public:
	OBEdge(OBNodeBase *bgn,OBNodeBase *end,OBExprBase *expr) 
	{
		_bgn = bgn; 
		_end = end; 
		_expr = expr;
		_closure = false;
	}
	~OBEdge() { if (_expr) { delete _expr; _expr = NULL;} }
	void ReplaceExpr(OBExprBase *expr)
	{
		if (_expr) delete _expr;
		_expr = expr;
	}
	void        SetClosure()        {_closure = true;}
	bool        IsClosure()         {return(_closure);}
	bool        Eval(OBEdgeBase *e) {return(_expr->Eval(e));}
	unsigned int GetBO()            {return(_expr->GetBO());}
	OBExprBase *GetExpr()           {return(_expr);}
	OBNodeBase *GetBgn()            {return(_bgn);}
	OBNodeBase *GetEnd()            {return(_end);}
};

class OBNode : public OBNodeBase
{
	int        _stereo;
	int        _vb;
	OBNodeBase *_match;
	OBExprBase *_expr;
public:
	OBNode() {_stereo=0; _match = NULL;}
	virtual ~OBNode();
	int     GetFormalCharge()         const  {return(_expr->GetFormalCharge());}
	void    SetExpr(OBExprBase *expr)        {_expr = expr;}
	void    SetMatch(OBNodeBase *m)          {_match = m;}
	void    SetVectorBinding(int vb)         {_vb = vb;}
	void    ResetRecurs()                    {_expr->ResetRecurs();}
	void    Print(std::ostream &ofs)              {_expr->Print(ofs);}
	bool	Eval(OBNodeBase *n)       const  {return(_expr->Eval(n));}
	bool         IsAromatic()                {return(_expr->IsAromatic());}
	unsigned int GetVectorBinding()          {return(_vb);}
	unsigned int GetAtomicNum()       const  {return(_expr->GetAtomicNum());}
	OBNodeBase *GetMatch()                   {return(_match);}
	OBNode *BeginNbr(std::vector<OBEdgeBase*>::iterator &);
	OBNode *NextNbr(std::vector<OBEdgeBase*>::iterator &);
};

class OBSmartsPattern : public OBGraphBase
{
protected:
	bool                                            _done;
	bool                                            _single;
	std::string                                     _smarts;
	std::vector<std::vector<OBNodeBase*> >          _vmatch;
	std::vector<std::pair<OBNodeBase*,std::vector<OBEdgeBase*> > > _vm;
public:
	OBSmartsPattern();
	OBSmartsPattern(const OBSmartsPattern &);
	~OBSmartsPattern();


	//initialization and graph modification methods
	bool    Init(const char *);
	bool    Init(std::string &);
	OBNode *NewNode(OBExprBase*);
	OBEdge *NewEdge(OBNode*,OBNode*,OBExprBase*);
	
	//iterator methods
	OBNode* Begin(std::vector<OBNodeBase*>::iterator &); 
	OBNode* Next(std::vector<OBNodeBase*>::iterator &);
	OBEdge* Begin(std::vector<OBEdgeBase*>::iterator &);
	OBEdge* Next(std::vector<OBEdgeBase*>::iterator &);

	//routines that perform substructure search
	bool        SingleMatch()                    const{return(_single);}
	bool        FinishedMatch()                  const {return(_done);}
	void        SetSingleMatch(bool state)       {_single = state;}
	void        SetFinishedMatch(bool state)     {_done = state;}
	void        PrepForMatch();
	void        PushBack(std::vector<OBNodeBase*> &v) {_vmatch.push_back(v);}
	void        ClearMatches();
	OBNodeBase *GetFirstSeed()                   {return(_vm.begin()->first);}
	bool        SeedMatch(OBNodeBase*);
	bool        RestrictedMatch(OBGraphBase&,std::vector<std::pair<int,int> >&,bool) 
	{
		std::cerr << "need to implement OBSmartsPattern::RestrictedMatch()" << std::endl;
		exit(0);
		return(false);
	}
	std::vector<std::pair<OBNodeBase*,std::vector<OBEdgeBase*> > >::iterator BgnMatch() {return(_vm.begin());}
	std::vector<std::pair<OBNodeBase*,std::vector<OBEdgeBase*> > >::iterator EndMatch() {return(_vm.end());}

	//data mod methods
	void SetSMARTS(const char *buf) {_smarts = buf;}

	//methods to return match data
	void GetMatches(std::vector<std::vector<OBNodeBase*> > &vm) const {vm = _vmatch;}
	void GetMatches(std::vector<std::vector<int> > &vm);
	std::vector<std::vector<int> > GetMapList();
	std::vector<std::vector<int> > GetUMapList();

	//misc data retrieval methods
	int           GetFormalCharge(int idx)  {return(_vatom[idx]->GetFormalCharge());}
	bool          IsValid()                 {return(!_vatom.empty());}
	bool          Empty()                   {return(_vatom.empty());}
	std::string   GetSMARTS()               const {return(_smarts);}
	unsigned int  NumMatches()              {return((_vmatch.empty()) ? 0 : _vmatch.size());}
	unsigned int  NumAtoms()                {return(NumNodes());}
	unsigned int  NumBonds()                {return(NumEdges());}
	unsigned int  GetVectorBinding(int idx) {return(((OBNode*)_vatom[idx])->GetVectorBinding());}
	unsigned int  GetAtomicNum(int idx)     const {return(_vatom[idx]->GetAtomicNum());}

	//stream output methods
	void WriteMapList(std::ostream&);
};


} //namespace OpenBabel

#endif //SMARTS_H



