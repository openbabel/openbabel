/**********************************************************************
typer.h - Open Babel atom and aromaticity typer.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison

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

#ifndef OB_TYPER_H
#define OB_TYPER_H

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>

#include <vector>
#include <string>

#include <openbabel/parsmart.h>
#include <openbabel/data.h>

namespace OpenBabel
{

  // Forward declaration
  class OBSmartsPattern;

// class introduction in typer.cpp
class OBAPI OBAtomTyper : public OBGlobalDataBase
{
  std::vector<std::pair<OBSmartsPattern*,int> >            _vinthyb; //!< internal hybridization rules
  std::vector<std::pair<OBSmartsPattern*,std::string> >    _vexttyp; //!< external atom type rules

public:
    OBAtomTyper();
    OBAtomTyper(const OBAtomTyper& rhs) {abort();}
    ~OBAtomTyper();

    //swig is requiring these, but I can't figure out how to make it not, so definte with abort
    const OBAtomTyper& operator=(const OBAtomTyper& rhs) {abort();}

    void ParseLine(const char*);
    //! \return the number of internal hybridization rules
    size_t GetSize()                 { return _vinthyb.size(); }

    //! Assign atomic hybridization (1 = sp, 2 = sp2, 3 = sp3...)
    void AssignHyb(OBMol&);
    //! Assign external atomic types (i.e., EXTTYP lines in atomtyp.txt)
    void AssignTypes(OBMol&);
};

#ifndef THREAD_LOCAL
# define THREAD_LOCAL
#endif
#ifndef OB_EXTERN
#error OB_EXTERN
#endif
//! Global OBAtomTyper for marking internal valence, hybridization,
//!  and atom types (for internal and external use)
THREAD_LOCAL OB_EXTERN OBAtomTyper      atomtyper;

// class introduction in typer.cpp
class OBAPI OBAromaticTyper
{
public:
    OBAromaticTyper() {};
    ~OBAromaticTyper() {};

    //! Assign aromaticity flag to atoms and bonds
    void AssignAromaticFlags(OBMol &);
};

//! Global OBAromaticTyper for detecting aromatic atoms and bonds
THREAD_LOCAL OB_EXTERN OBAromaticTyper  aromtyper;

// class introduction in typer.cpp
class OBAPI OBRingTyper : public OBGlobalDataBase
{
  std::vector<std::pair<OBSmartsPattern*,std::string> >    _ringtyp; //!< ring type rules

public:
    OBRingTyper();
    ~OBRingTyper();

    void ParseLine(const char*);
    //! \return the number of SMARTS patterns
    size_t GetSize()                 { return _ringtyp.size();}

    //! Assign external atomic types (ringtyp.txt)
    void AssignTypes(OBMol&);
};


} //namespace OpenBabel

#endif // OB_TYPER_H

//! \file typer.h
//! \brief Open Babel atom and aromaticity typer.
