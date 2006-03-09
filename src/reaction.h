/**********************************************************************
reaction.h - Handle chemical reactions (i.e., lists of reagents and products).

Copyright (C) 2005 by Chris Morley
 
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

#ifndef OB_REACT_H
#define OB_REACT_H

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include <vector>
#include <boost/shared_ptr.hpp>
#include "mol.h"

using boost::shared_ptr;

namespace OpenBabel
{

//!\brief Used to store chemical reactions (i.e., reactants -> products)
//!
//! Reactants and products stored as smart pointers to molecules stored elsewhere.
//!
//! For performing actual reaction transformations (i.e., deleting atoms,
//! changing bonds, etc.) use the OBChemTsfm class.
class OBAPI OBReaction : public OBBase
{
private:
	std::vector<boost::shared_ptr<OBMol> > _reactants;
	std::vector<boost::shared_ptr<OBMol> > _products;
	std::string _title;
	std::string _comment;
	bool _reversible;
public:
	OBReaction() : _reversible(false)
	{}

	int NumReactants() const
	{ return _reactants.size(); }

	int NumProducts()const
	{ return _products.size(); }
	
	void AddReactant(const boost::shared_ptr<OBMol> sp)
	{ _reactants.push_back(sp); }
	
	void AddProduct(const boost::shared_ptr<OBMol> sp)
	{ _products.push_back(sp); }
	
	boost::shared_ptr<OBMol> GetReactant(const unsigned i)
	{
		boost::shared_ptr<OBMol> sp;
		if(i<_reactants.size())
			sp = _reactants[i];
		return sp; //returns empty if out of range
	}
	boost::shared_ptr<OBMol> GetProduct(const unsigned i)
	{
		boost::shared_ptr<OBMol> sp;
		if(i<_products.size())
			sp = _products[i];
		return sp; //returns empty if out of range
	}
	std::string GetTitle()	const { return _title; }
	std::string GetComment()	const { return _comment; }
	void SetTitle(const std::string& title) { _title=title; }
	void SetComment(const std::string& comment) { _comment=comment; }

	bool IsReversible() const       {return _reversible;}
	void SetReversible(bool b=true) {_reversible=b;}

	static const char* ClassDescription()
	{
    return "reactions\n";
	}
};


} //namespace OpenBabel
#endif

//! \file reaction.h
//! \brief Handle chemical reactions (i.e., lists of reagents and products).
