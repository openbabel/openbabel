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

#ifndef BASE_H
#define BASE_H

#include <vector>
#include <iostream>

using namespace std;

namespace OpenBabel
{

class OBBase;
class OBNodeBase;
class OBEdgeBase;
class OBGraphBase;

//
//Base Classes
//

class OBBase
{
public:
};

class OBNodeBase : public OBBase
{
protected:
	unsigned short int  _idx;
	OBGraphBase        *_parent;
	vector<OBEdgeBase*> _vbond;
public:
	bool Visit;
	         OBNodeBase()                            {Visit = false;}
	virtual ~OBNodeBase()                            {}
	void                  SetIdx(int idx)            {_idx = idx;}
	void                  AddEdge(OBEdgeBase *b)     {_vbond.push_back(b);}
	void                  SetParent(OBGraphBase*);
	virtual int           GetFormalCharge()          const {((OBNodeBase*)this)->Error(1); return(0);}
	virtual unsigned int  GetIdx()                   const {return(_idx);}
	virtual unsigned int  ExplicitHydrogenCount()    const {((OBNodeBase*)this)->Error(22); return(0);}
	virtual unsigned int  ImplicitHydrogenCount()    const {((OBNodeBase*)this)->Error(22); return(0);}
	virtual unsigned int  GetImplicitValence()       const {((OBNodeBase*)this)->Error(3); return(0);}
	virtual unsigned int  GetValence()               const {return(_vbond.size());}
	virtual unsigned int  GetHvyValence()            const {((OBNodeBase*)this)->Error(4);return(0);}
	virtual unsigned int  KBOSum()                const {((OBNodeBase*)this)->Error(5);return(0);}
	virtual unsigned int  GetHyb()                   const {((OBNodeBase*)this)->Error(6);return(0);}
	virtual unsigned int  MemberOfRingCount()        const {((OBNodeBase*)this)->Error(7);return(0);}
	virtual unsigned int  GetAtomicNum()			 const {((OBNodeBase*)this)->Error(8);return(0);}
	virtual void          SetMatch(OBNodeBase*)      {}
	virtual void          SetAromatic()              {}
	virtual bool          IsInRingSize(int)          const {((OBNodeBase*)this)->Error(9);return(false);}
	virtual bool          IsAromatic()               const {((OBNodeBase*)this)->Error(10);return(false);}
	virtual bool          IsInRing()                 const {((OBNodeBase*)this)->Error(11);return(false);}
	virtual bool          Eval(OBNodeBase*)          const {((OBNodeBase*)this)->Error(12);return(false);}
	virtual bool          IsConnected(OBNodeBase*);
	virtual OBNodeBase   *GetMatch()                 {((OBNodeBase*)this)->Error(13); return(NULL);}
	virtual OBGraphBase  *GetParent()                {return(_parent);}
	OBEdgeBase           *Begin(vector<OBEdgeBase*>::iterator&);
	OBEdgeBase           *Next(vector<OBEdgeBase*>::iterator&);
	OBNodeBase           *BeginNbr(vector<OBEdgeBase*>::iterator&);
	OBNodeBase           *NextNbr(vector<OBEdgeBase*>::iterator&);
	void Error(int f)     {cerr << "atom vf called = " << f << endl;}
};

class OBEdgeBase : public OBBase
{
protected:
	unsigned short int _idx;
	OBNodeBase        *_bgn;
	OBNodeBase        *_end;
	OBGraphBase       *_parent;
public:
	bool Visit;

	         OBEdgeBase() {Visit = false; _bgn = _end = NULL;}
	         OBEdgeBase(OBNodeBase *bgn,OBNodeBase *end) {}
	virtual ~OBEdgeBase()                  {}
	virtual void SetClosure()      {}
	virtual bool IsAromatic()      const {((OBEdgeBase*)this)->Error(1); return(false);}
	virtual bool IsInRing()        const {((OBEdgeBase*)this)->Error(2); return(false);}
	virtual bool IsClosure()       {((OBEdgeBase*)this)->Error(3); return(false);}
	virtual bool Eval(OBEdgeBase*) {((OBEdgeBase*)this)->Error(4); return(false);}
	virtual unsigned int  GetBO()  const {((OBEdgeBase*)this)->Error(5); return(0);}

	virtual OBGraphBase* GetParent() {return(_parent);}
	unsigned int GetIdx()          {return(_idx);}
	void SetIdx(int idx)           {_idx = idx;}
	void SetEnd(OBNodeBase *n)     {_end = n;}
	void SwapEnds()                {swap(_bgn,_end);}
	void SetParent(OBGraphBase*);               
	OBNodeBase *GetBgn()           {return(_bgn);}
	OBNodeBase *GetEnd()           {return(_end);}
	void Error(int f) {cerr << "bond vf err = " << f << endl;}
};

class OBGraphBase : public OBBase
{
protected:
	bool                 _vlock;
	vector<OBNodeBase*>  _vatom;
	vector<OBEdgeBase*>  _vbond;
public:
	OBGraphBase() 
	{
		_vlock = false; 
		_vatom.clear(); 
		_vbond.clear();
	}
	OBGraphBase(const OBGraphBase &src)
	  {
	    
	  }

	virtual     ~OBGraphBase()        {}
	unsigned int NumNodes()           {return(_vatom.empty() ? 0 : _vatom.size());}
	unsigned int NumEdges()           {return(_vbond.empty() ? 0 : _vbond.size());}
	void        ResetVisitFlags();
	bool        SetVisitLock(bool);
	bool        GetVisitLock()        {return(_vlock);}
	OBNodeBase *Begin(vector<OBNodeBase*>::iterator&); 
	OBNodeBase *Next(vector<OBNodeBase*>::iterator&);
	OBEdgeBase *Begin(vector<OBEdgeBase*>::iterator&);
	OBEdgeBase *Next(vector<OBEdgeBase*>::iterator&);

	//substructure search functions
	virtual bool SingleMatch()                  const {return(false);}
	virtual void SetSingleMatch(bool)           {}
	virtual bool FinishedMatch()                const {return(false);}
	virtual void SetFinishedMatch(bool)         {}
	virtual void ClearMatches()                 {}
	virtual void PushBack(vector<OBNodeBase*>&) {}
	virtual void PrepForMatch()                 {}
	virtual vector<pair<OBNodeBase*,vector<OBEdgeBase*> > >::iterator BgnMatch() {return((vector<pair<OBNodeBase*,vector<OBEdgeBase*> > >::iterator) NULL);}
	virtual vector<pair<OBNodeBase*,vector<OBEdgeBase*> > >::iterator EndMatch() {return((vector<pair<OBNodeBase*,vector<OBEdgeBase*> > >::iterator) NULL);}
	virtual OBNodeBase *GetFirstSeed() {return((OBNodeBase*)NULL);}
	bool Match(OBGraphBase &,bool singleMatch=false);
	bool Match(OBGraphBase &,
		       vector<pair<OBNodeBase*,vector<OBEdgeBase*> > >::iterator,
		       vector<OBEdgeBase*>::iterator);
};

} //namespace OpenBabel

#endif //BASE_H
