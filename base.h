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

namespace OpenEye
{

class OEBase;
class OENodeBase;
class OEEdgeBase;
class OEGraphBase;

//
//Base Classes
//

class OEBase
{
public:
};

class OENodeBase : public OEBase
{
protected:
	unsigned short int  _idx;
	OEGraphBase        *_parent;
	vector<OEEdgeBase*> _vbond;
public:
	bool Visit;
	         OENodeBase()                            {Visit = false;}
	virtual ~OENodeBase()                            {}
	void                  SetIdx(int idx)            {_idx = idx;}
	void                  AddEdge(OEEdgeBase *b)     {_vbond.push_back(b);}
	void                  SetParent(OEGraphBase*);
	virtual int           GetFormalCharge()          const {((OENodeBase*)this)->Error(1); return(0);}
	virtual unsigned int  GetIdx()                   const {return(_idx);}
	virtual unsigned int  ExplicitHydrogenCount()    const {((OENodeBase*)this)->Error(22); return(0);}
	virtual unsigned int  ImplicitHydrogenCount()    const {((OENodeBase*)this)->Error(22); return(0);}
	virtual unsigned int  GetImplicitValence()       const {((OENodeBase*)this)->Error(3); return(0);}
	virtual unsigned int  GetValence()               const {return(_vbond.size());}
	virtual unsigned int  GetHvyValence()            const {((OENodeBase*)this)->Error(4);return(0);}
	virtual unsigned int  KBOSum()                const {((OENodeBase*)this)->Error(5);return(0);}
	virtual unsigned int  GetHyb()                   const {((OENodeBase*)this)->Error(6);return(0);}
	virtual unsigned int  MemberOfRingCount()        const {((OENodeBase*)this)->Error(7);return(0);}
	virtual unsigned int  GetAtomicNum()			 const {((OENodeBase*)this)->Error(8);return(0);}
	virtual void          SetMatch(OENodeBase*)      {}
	virtual void          SetAromatic()              {}
	virtual bool          IsInRingSize(int)          const {((OENodeBase*)this)->Error(9);return(false);}
	virtual bool          IsAromatic()               const {((OENodeBase*)this)->Error(10);return(false);}
	virtual bool          IsInRing()                 const {((OENodeBase*)this)->Error(11);return(false);}
	virtual bool          Eval(OENodeBase*)          const {((OENodeBase*)this)->Error(12);return(false);}
	virtual bool          IsConnected(OENodeBase*);
	virtual OENodeBase   *GetMatch()                 {((OENodeBase*)this)->Error(13); return(NULL);}
	virtual OEGraphBase  *GetParent()                {return(_parent);}
	OEEdgeBase           *Begin(vector<OEEdgeBase*>::iterator&);
	OEEdgeBase           *Next(vector<OEEdgeBase*>::iterator&);
	OENodeBase           *BeginNbr(vector<OEEdgeBase*>::iterator&);
	OENodeBase           *NextNbr(vector<OEEdgeBase*>::iterator&);
	void Error(int f)     {cerr << "atom vf called = " << f << endl;}
};

class OEEdgeBase : public OEBase
{
protected:
	unsigned short int _idx;
	OENodeBase        *_bgn;
	OENodeBase        *_end;
	OEGraphBase       *_parent;
public:
	bool Visit;

	         OEEdgeBase() {Visit = false; _bgn = _end = NULL;}
	         OEEdgeBase(OENodeBase *bgn,OENodeBase *end) {}
	virtual ~OEEdgeBase()                  {}
	virtual void SetClosure()      {}
	virtual bool IsAromatic()      const {((OEEdgeBase*)this)->Error(1); return(false);}
	virtual bool IsInRing()        const {((OEEdgeBase*)this)->Error(2); return(false);}
	virtual bool IsClosure()       {((OEEdgeBase*)this)->Error(3); return(false);}
	virtual bool Eval(OEEdgeBase*) {((OEEdgeBase*)this)->Error(4); return(false);}
	virtual unsigned int  GetBO()  const {((OEEdgeBase*)this)->Error(5); return(0);}

	virtual OEGraphBase* GetParent() {return(_parent);}
	unsigned int GetIdx()          {return(_idx);}
	void SetIdx(int idx)           {_idx = idx;}
	void SetEnd(OENodeBase *n)     {_end = n;}
	void SwapEnds()                {swap(_bgn,_end);}
	void SetParent(OEGraphBase*);               
	OENodeBase *GetBgn()           {return(_bgn);}
	OENodeBase *GetEnd()           {return(_end);}
	void Error(int f) {cerr << "bond vf err = " << f << endl;}
};

class OEGraphBase : public OEBase
{
protected:
	bool                 _vlock;
	vector<OENodeBase*>  _vatom;
	vector<OEEdgeBase*>  _vbond;
public:
	OEGraphBase() 
	{
		_vlock = false; 
		_vatom.clear(); 
		_vbond.clear();
	}
	OEGraphBase(const OEGraphBase &src)
	  {
	    
	  }

	virtual     ~OEGraphBase()        {}
	unsigned int NumNodes()           {return(_vatom.empty() ? 0 : _vatom.size());}
	unsigned int NumEdges()           {return(_vbond.empty() ? 0 : _vbond.size());}
	void        ResetVisitFlags();
	bool        SetVisitLock(bool);
	bool        GetVisitLock()        {return(_vlock);}
	OENodeBase *Begin(vector<OENodeBase*>::iterator&); 
	OENodeBase *Next(vector<OENodeBase*>::iterator&);
	OEEdgeBase *Begin(vector<OEEdgeBase*>::iterator&);
	OEEdgeBase *Next(vector<OEEdgeBase*>::iterator&);

	//substructure search functions
	virtual bool SingleMatch()                  const {return(false);}
	virtual void SetSingleMatch(bool)           {}
	virtual bool FinishedMatch()                const {return(false);}
	virtual void SetFinishedMatch(bool)         {}
	virtual void ClearMatches()                 {}
	virtual void PushBack(vector<OENodeBase*>&) {}
	virtual void PrepForMatch()                 {}
	virtual vector<pair<OENodeBase*,vector<OEEdgeBase*> > >::iterator BgnMatch() {return(NULL);}
	virtual vector<pair<OENodeBase*,vector<OEEdgeBase*> > >::iterator EndMatch() {return(NULL);}
	virtual OENodeBase *GetFirstSeed() {return((OENodeBase*)NULL);}
	bool Match(OEGraphBase &,bool singleMatch=false);
	bool Match(OEGraphBase &,
		       vector<pair<OENodeBase*,vector<OEEdgeBase*> > >::iterator,
		       vector<OEEdgeBase*>::iterator);
};

} //namespace OpenEye

#endif //BASE_H
