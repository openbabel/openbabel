/**********************************************************************
data.h - Global data and resource file parsers.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#ifndef OB_DATA_H
#define OB_DATA_H

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#endif

#if HAVE_FSTREAM
#include <fstream>
#elif HAVE_FSTREAM_H
#include <fstream.h>
#endif

#include <vector>
#include <string>

namespace OpenBabel
{

class OBElement;
class OBAtom;
class OBElementTable;

//! \brief Base data table class, handles reading data files
//!
//! Base data table class--reads ASCII data files in various formats
//! -# Checks for the environment variable _envvar (defaults to "BABEL_DATADIR")
//!     - Tries the _subdir directory if defined (def. "data") and then the main directory
//! -# Checks for the directory _dir (def. determined by the build environment)
//!     - Tries the subdirectory corresponding to this version, then the main directory
//! -# Reverts to the compiled-in default data
class OBAPI OBGlobalDataBase
{
protected:
    bool         _init;		//!< has the data been read already
    const char  *_dataptr;	//!< default data table if file is unreadable
    std::string  _filename;	//!< file to search for
    std::string  _dir;		//!< data directory for file if _envvar fails
    std::string  _subdir;		//!< subdirectory (if using environment variable)
    std::string  _envvar;		//!< environment variable to check first
public:
    //! Constructor
    OBGlobalDataBase()
    {
        _init = false;
        _dataptr = (char*)NULL;
    }
    //! Destructor
    virtual ~OBGlobalDataBase()
    {}
    //! Read in the data file, falling back as needed
    void  Init();
    //! Set the directory before calling Init()
    void  SetReadDirectory(char *dir)
    {
        _dir = dir;
    }
    //! Set the environment variable to use before calling Init()
    void  SetEnvironmentVariable(char *var)
    {
        _envvar = var;
    }
    //! Specified by particular table classes (parses an individual data line)
    virtual void ParseLine(const char*)
    {}
}
;

//! \brief Individual element data type
//!
//! Stores a variety of data about an individual element
class OBAPI OBElement
{
    int _num;
    char _symbol[3];
    double _Rcov,_Rbo,_Rvdw,_mass,_elNeg;
    int _maxbonds;
public:
    OBElement()    {}
    OBElement(int num, const char *sym, double rcov, double rbo,
              double rvdw, int maxbo, double mass, double elNeg) :
      _num(num), _Rcov(rcov), _Rbo(rbo), _Rvdw(rvdw), _maxbonds(maxbo),
      _mass(mass), _elNeg(elNeg)
    {
      strncpy(_symbol, sym, 3);
    }

    int GetAtomicNum()         {       return(_num);    }
    char *GetSymbol()          {       return(_symbol); }
    double GetCovalentRad()    {       return(_Rcov);   }
    //! \deprecated Use GetCovalentRad() instead
    double GetBoRad()          {       return(_Rbo);    }
    double GetVdwRad()         {       return(_Rvdw);   }
    double GetMass()           {       return(_mass);   }
    int GetMaxBonds()          {       return(_maxbonds);}
    double GetElectroNeg()     {       return(_elNeg);  }
};

// class introduction in data.cpp
class OBAPI OBElementTable : public OBGlobalDataBase
{
    std::vector<OBElement*> _element;

public:

    OBElementTable(void);
    ~OBElementTable();

    void  ParseLine(const char*);

    //! Returns the number of elements in the periodic table
    int GetNumberOfElements() { return _element.size(); }

    //! \deprecated Does not properly handle 'D' or 'T' hydrogen isotopes
    int   GetAtomicNum(const char *);
    //! Returns the atomic number matching the element symbol passed
    //! or 0 if not defined. For 'D' or 'T' hydrogen isotopes, will return
    //! a value in the second argument
    int   GetAtomicNum(const char *, int &iso);
    //! Returns the element symbol matching the atomic number passed
    char *GetSymbol(int);
    //! Returns the van der Waals radius for this atomic number
    double GetVdwRad(int);
    //! Returns the covalent radius for this atomic number
    double GetCovalentRad(int);
    //! \deprecated -- Use OBElementTable::GetCovalentRad()
    double GetBORad(int);
    //! Returns the average atomic mass for this element.
    //! For exact isotope masses, use OpenBabel::OBIsotopeTable
    double GetMass(int);
    //! Returns a "corrected" bonding radius based on the hybridization.
    //! Scales the covalent radius by 0.95 for sp2 and 0.90 for sp hybrids
    double CorrectedBondRad(int,int = 3); // atomic #, hybridization
    //! Returns a "corrected" vdW radius based on the hybridization.
    //! Scales the van der Waals radius by 0.95 for sp2 and 0.90 for sp hybrids
    double CorrectedVdwRad(int,int = 3); // atomic #, hybridization
    //! Returns the maximum expected number of bonds to this element
    int	GetMaxBonds(int);
    //! Returns the Pauling electronegativity for this element
    double GetElectroNeg(int);
};

// class introduction in data.cpp
class OBAPI OBIsotopeTable : public OBGlobalDataBase
{
    std::vector<std::vector<std::pair <unsigned int, double> > > _isotopes;

public:

    OBIsotopeTable(void);
    ~OBIsotopeTable()
    {}

    void	ParseLine(const char*);
    //! Return the exact masss of the isotope
    //!   (or by default (i.e. "isotope 0") the most abundant isotope)
    double	GetExactMass(const unsigned int atomicNum,
                        const unsigned int isotope = 0);
};

// class introduction in data.cpp
class OBAPI OBTypeTable : public OBGlobalDataBase
{
    int    _linecount;
    int    _ncols,_nrows,_from,_to;
    std::vector<std::string> _colnames;
    std::vector<std::vector<std::string> > _table;

public:

    OBTypeTable(void);
    ~OBTypeTable()
    {}

    void ParseLine(const char*);
    bool SetFromType(char*);
    bool SetToType(char*);
    bool Translate(char*,char*); // to, from
    bool Translate(std::string &,std::string &); // to, from
};

// Used by other code for reading files
#ifdef WIN32
#define FILE_SEP_CHAR "\\"
#else
#define FILE_SEP_CHAR "/"
#endif

}

#endif //DATA_H
