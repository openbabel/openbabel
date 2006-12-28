/**********************************************************************
reaction.h - Handle chemical reactions (i.e., lists of reagents and products).

Copyright (C) 2004 by Chris Morley
 
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
//#include <boost/shared_ptr.hpp>

using namespace std;

namespace OpenBabel
{

class OBBase;

/** \class OBReaction reaction.h <openbabel/reaction.h>
    \brief Used to store chemical reactions (i.e., reactants -> products)

    Reactants and products stored as pointers to molecules stored elsewhere,
    since the molecules may be involved in other reactions.

    For performing actual reaction transformations (i.e., deleting atoms,
    changing bonds, etc.) use the OBChemTsfm class.
**/
class OBAPI OBReaction : public OBBase
{
public:
    vector<OBMol*> reactants;
    vector<OBMol*> products;
    string title;
};

}
#endif

//! \file reaction.h
//! \brief Handle chemical reactions (i.e., lists of reagents and products).
