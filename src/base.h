/**********************************************************************
base.h - Base classes to build a graph

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_BASE_H
#define OB_BASE_H

#include <vector>

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#endif

namespace OpenBabel
{

class OBBase;
class OBNodeBase;
class OBEdgeBase;
class OBGraphBase;

/** \brief Base Class

 The various classes (Atom, Bond, Molecule) inherit from base classes--
 OBBase is just a placeholder class
*/
class OBBase
{
public:
};

/** \brief Node Base Class

The base class for nodes (e.g. atoms) in a graph-theoretical representation.
 */
class OBNodeBase : public OBBase
{
protected:
	//! What is my unique node index? 	GetIdx(), SetIdx(int idx).
        unsigned short int  _idx;
	//! To which graph do I belong?	GetParent(), SetParent(OBGraphBase*).
        OBGraphBase        *_parent;
	/** \brief What edges or bonds do I have? \sa AddEdge(), GetValence(),
	    This node is assumed to be one of the edge endpoints.
	*/
        std::vector<OBEdgeBase*> _vbond;

public:
	//! Used internally by graph traversal algorithms
	bool Visit;

        //! Constructor
        OBNodeBase()                            {Visit = false;}
        //! Destructor
	virtual ~OBNodeBase()                            {}
	virtual unsigned int  GetIdx()                   const {return(_idx);}
        //! Set my unique node number.
	void                  SetIdx(int idx)            {_idx = idx;}
	virtual OBGraphBase  *GetParent()                {return(_parent);}
        //! Set my parent graph.
	void                  SetParent(OBGraphBase*);
        //! Add an edge. Assumes this node is an endpoint of \a b.
	void                  AddEdge(OBEdgeBase *b)     {_vbond.push_back(b);}
	//! How many edges?
	virtual unsigned int  GetValence()               const {return(_vbond.size());}
	//! Iterate over my edges, returning my connected nodes.
	OBNodeBase           *BeginNbr(std::vector<OBEdgeBase*>::iterator&);
	OBNodeBase           *NextNbr(std::vector<OBEdgeBase*>::iterator&);
	//! Iterate over my edges, returning the edges.
	OBEdgeBase           *Begin(std::vector<OBEdgeBase*>::iterator&);
	OBEdgeBase           *Next(std::vector<OBEdgeBase*>::iterator&);
	//! Is \a OBNodeBase* a beginning or end of some edge of mine?
	virtual bool          IsConnected(OBNodeBase*);

	//! Used to signal erroneously calling some of the nonfunctional members
	void Error(int f)     {std::cerr << "atom vf called = " << f << std::endl;}
	/** \defgroup Nonfunctional members of OpenBabel::OBNodeBase */
	//@{
	/**  \name Nonfunctional members
	     Nonfunctional or unused members of OpenBabel::OBNodeBase 
	     are to be added by a derived class (e.g., OpenBabel::OBAtom)
	     @{
	*/
	virtual int           GetFormalCharge()          const {((OBNodeBase*)this)->Error(1); return(0);}
	virtual unsigned int  ExplicitHydrogenCount()    const {((OBNodeBase*)this)->Error(22); return(0);}
	virtual unsigned int  ImplicitHydrogenCount()    const {((OBNodeBase*)this)->Error(22); return(0);}
	virtual unsigned int  GetImplicitValence()       const {((OBNodeBase*)this)->Error(3); return(0);}
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
	virtual OBNodeBase   *GetMatch()                 {((OBNodeBase*)this)->Error(13); return(NULL);}
        //@}
        //@}
};

/** \brief Edge Base Class

The base class for edges (e.g. bonds) in a graph-theoretical representation.
 */
class OBEdgeBase : public OBBase
{
protected:
        /// What is my unique edge index? GetIdx(), SetIdx(int idx)
	unsigned short int _idx;
	OBNodeBase        *_bgn;     //!< I connect one node
	OBNodeBase        *_end;     //!< to another node
	//! To which graph do I belong?	  GetParent(), SetParent(OBGraphBase*).
	OBGraphBase       *_parent;
public:
	bool Visit;

	//! Constructor
	OBEdgeBase() {Visit = false; _bgn = _end = NULL;}
	//! Constructor (populate from two OBNodeBase objects)
	OBEdgeBase(OBNodeBase *bgn,OBNodeBase *end) {}
	//! Destructor
	virtual ~OBEdgeBase()          {}

	virtual OBGraphBase* GetParent() {return(_parent);}
	void SetParent(OBGraphBase*);               
	unsigned int GetIdx()          {return(_idx);}
	void SetIdx(int idx)           {_idx = idx;}
	void SetBgn(OBNodeBase *n)     {_bgn = n;}
	void SetEnd(OBNodeBase *n)     {_end = n;}
	void SwapEnds()                {std::swap(_bgn,_end);}
	OBNodeBase *GetBgn()           {return(_bgn);}
	OBNodeBase *GetEnd()           {return(_end);}

	void Error(int f) {std::cerr << "bond vf err = " << f << std::endl;}
	/** \defgroup Nonfunctional members of OpenBabel::OBEdgeBase */
	//@{
	/**  \name Nonfunctional members
	  Nonfunctional or unused members of OpenBabel::OBEdgeBase
	  are to be added by a derived class (e.g., OpenBabel::OBBond)
	  @{
	*/
	virtual void SetClosure()      {}
	virtual bool IsAromatic()      const {((OBEdgeBase*)this)->Error(1); return(false);}
	virtual bool IsInRing()        const {((OBEdgeBase*)this)->Error(2); return(false);}
	virtual bool IsClosure()       {((OBEdgeBase*)this)->Error(3); return(false);}
	virtual bool Eval(OBEdgeBase*) {((OBEdgeBase*)this)->Error(4); return(false);}
	virtual unsigned int  GetBO()  const {((OBEdgeBase*)this)->Error(5); return(0);}
	//@}
	//@}
};

/*! \brief Graph Base Class

The base class for graphs (e.g. rings, molecules, etc.) in a 
graph-theoretical representation.
 */
class OBGraphBase : public OBBase
{
protected:
	bool                 _vlock;
	std::vector<OBNodeBase*>  _vatom;
	std::vector<OBEdgeBase*>  _vbond;
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
	OBNodeBase *Begin(std::vector<OBNodeBase*>::iterator&); 
	OBNodeBase *Next(std::vector<OBNodeBase*>::iterator&);
	OBEdgeBase *Begin(std::vector<OBEdgeBase*>::iterator&);
	OBEdgeBase *Next(std::vector<OBEdgeBase*>::iterator&);

	//substructure search functions
	virtual bool SingleMatch()                  const {return(false);}
	virtual void SetSingleMatch(bool)           {}
	virtual bool FinishedMatch()                const {return(false);}
	virtual void SetFinishedMatch(bool)         {}
	virtual void ClearMatches()                 {}
	virtual void PushBack(std::vector<OBNodeBase*>&) {}
	virtual void PrepForMatch()                 {}
	virtual std::vector<std::pair<OBNodeBase*,std::vector<OBEdgeBase*> > >::iterator BgnMatch() 
	  {return((std::vector<std::pair<OBNodeBase*,std::vector<OBEdgeBase*> > >::iterator) NULL);}
	virtual std::vector<std::pair<OBNodeBase*,std::vector<OBEdgeBase*> > >::iterator EndMatch()
	  {return((std::vector<std::pair<OBNodeBase*,std::vector<OBEdgeBase*> > >::iterator) NULL);}
	virtual OBNodeBase *GetFirstSeed() {return((OBNodeBase*)NULL);}
	bool Match(OBGraphBase &,bool singleMatch=false);
	bool Match(OBGraphBase &,
		       std::vector<std::pair<OBNodeBase*,std::vector<OBEdgeBase*> > >::iterator,
		       std::vector<OBEdgeBase*>::iterator);
};

} //namespace OpenBabel

#endif // OB_BASE_H
