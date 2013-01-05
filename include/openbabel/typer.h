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
  //    int                                                      _rc;
  std::vector<std::vector<int> >                           _mlist; //!< match list for atom typing
  std::vector<std::pair<OBSmartsPattern*,int> >            _vinthyb; //!< internal hybridization rules
  std::vector<std::pair<OBSmartsPattern*,int> >            _vimpval; //!< internal implicit valence rules
  std::vector<std::pair<OBSmartsPattern*,std::string> >    _vexttyp; //!< external atom type rules

public:
    OBAtomTyper();
    ~OBAtomTyper();

    void ParseLine(const char*);
    //! \return the number of implicit valence rules
    size_t GetSize()                 { return _vimpval.size();}

    //! Assign atomic hybridization (1 = sp, 2 = sp2, 3 = sp3...)
    void AssignHyb(OBMol&);
    //! Assign external atomic types (i.e., EXTTYP lines in atomtyp.txt)
    void AssignTypes(OBMol&);
    //! Assign implicit valence (i.e., given an atomic type, what is the
    //! expected number of bonds to this atom
    void AssignImplicitValence(OBMol&, bool CanBeLessThanActual=false);
    //! Correct typing, valence, and hybridization for aromatic nitrogen atoms
    void CorrectAromaticNitrogens(OBMol&);
};

// class introduction in typer.cpp
class OBAPI OBAromaticTyper : public OBGlobalDataBase
{
    std::vector<bool>             _vpa;   //!< potentially aromatic atoms
    std::vector<bool>             _visit;
    std::vector<bool>             _root;
    std::vector<std::vector<int> >     _mlist;
    std::vector<OBSmartsPattern*> _vsp;   //!< SMARTS of potentially aromatic atoms
    std::vector<std::pair<int,int> >   _verange; //!< min and max number of electrons
    std::vector<std::pair<int,int> >   _velec;   //!< # electrons an atom contributes
public:
    OBAromaticTyper();
    ~OBAromaticTyper();

    //! \return the number of SMARTS patterns
    size_t GetSize()                 { return _vsp.size();}

    void ParseLine(const char*);
    //! Assign aromaticity flag to atoms and bonds
    void AssignAromaticFlags(OBMol &);
    //! "Anti-alias" potentially aromatic flags around a molecule
    //! (aromatic atoms need to have >= 2 neighboring ring atoms)
    void PropagatePotentialAromatic(OBAtom*);
    // Documentation in typer.cpp
    void SelectRootAtoms(OBMol &, bool avoidInnerRingAtoms = true);
    //! Remove 3-member rings from consideration
    void ExcludeSmallRing(OBMol &);
    //! Check aromaticity starting from the root atom, up to a specified depth
    void CheckAromaticity(OBAtom *root,int searchDepth);
    // Documentation in typer.cpp
    bool TraverseCycle(OBAtom *root, OBAtom *atom, OBBond *prev,
                       std::pair<int,int> &er,int depth);
};

// class introduction in typer.cpp
class OBAPI OBRingTyper : public OBGlobalDataBase
{
  std::vector<std::vector<int> >                           _mlist; //!< match list for atom typing
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
