/**********************************************************************
bondtyper.h - Bond typer to perceive connectivity and bond orders/types.

Copyright (C) 2003-2005 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_BONDTYPER_H
#define OB_BONDTYPER_H

#include <openbabel/parsmart.h>
#include <openbabel/data.h>

namespace OpenBabel
{

// class introduction in bondtyper.cpp
// Used for "perceiving" bonds, e.g. in XYZ or QM files with no bond info.
class OBAPI OBBondTyper : public OBGlobalDataBase
{
    //! SMARTS patterns for functional group typing
    std::vector<std::pair<OBSmartsPattern*, std::vector<int> > >	_fgbonds;
public:
    OBBondTyper();
    ~OBBondTyper();

    //! \name OBBondTyper Database Utilities
    //@{
    void ParseLine(const char*);
    //! \return the size of the database (for error checking)
    size_t GetSize()                 { return _fgbonds.size();}
    //@}

    //! \name Bond Perception Routines
    //@{
    //! Assign bonds to functional groups based on the bond typer database
    void AssignFunctionalGroupBonds(OBMol &mol);
    //@}
};

}

#endif // OB_BONDTYPER_H

//! \file bondtyper.h
//! \brief Bond typer to perceive connectivity and bond orders/types.
