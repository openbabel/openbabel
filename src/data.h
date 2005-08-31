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

#include "babelconfig.h"

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
    std::string  _subdir;	//!< subdirectory (if using environment variable)
    std::string  _envvar;	//!< environment variable to check first

public:
    //! Constructor
    OBGlobalDataBase()
    {
        _init = false;
        _dataptr = (char*)NULL;
    }
    //! Destructor
    virtual ~OBGlobalDataBase()                  {}
    //! Read in the data file, falling back as needed
    void  Init();
    //! Set the directory before calling Init()
    void  SetReadDirectory(char *dir)            { _dir = dir;    }
    //! Set the environment variable to use before calling Init()
    void  SetEnvironmentVariable(char *var)      { _envvar = var; }
    //! Specified by particular table classes (parses an individual data line)
    virtual void ParseLine(const char*)          {}
};

//! \brief Individual element data type
//!
//! Stores a variety of data about an individual element
class OBAPI OBElement
{
    int _num;
    char _symbol[3];
    std::string _name;
    double _Rcov,_Rvdw,_mass,_elNeg,_ionize,_elAffinity;
    double _red, _green, _blue;
    int _maxbonds;
public:
    OBElement()    {}
    OBElement(int num, const char *sym, double rcov, double rvdw,
	      int maxbo, double mass, double elNeg, double ionize,
	      double elAffin, double red, double green, double blue,
	      std::string name) :
      _num(num), _name(name), _Rcov(rcov), _Rvdw(rvdw), _mass(mass), 
      _elNeg(elNeg), _ionize(ionize), _elAffinity(elAffin), 
      _red(red), _green(green), _blue(blue),
      _maxbonds(maxbo)
    {
      strncpy(_symbol, sym, 3);
    }

    //! Returns the atomic number of this element
    int GetAtomicNum()         {       return(_num);    }
    //! Returns the atomic symbol for this element
    char *GetSymbol()          {       return(_symbol); }
    //! Returns the covalent radius of this element
    double GetCovalentRad()    {       return(_Rcov);   }
    //! Returns the van der Waals radius of this element
    double GetVdwRad()         {       return(_Rvdw);   }
    //! \return the standard atomic mass for this element (in amu)
    double GetMass()           {       return(_mass);   }
    //! \return the maximum expected number of bonds to this element
    int GetMaxBonds()          {       return(_maxbonds);}
    //! \return the Pauling electronegativity for this element
    double GetElectroNeg()     {       return(_elNeg);  }
    //! \return the ionization potential (in eV) of this element
    double GetIonization()     {       return(_ionize);  }
    //! \return the electron affinity (in eV) of this element
    double GetElectronAffinity(){      return(_elAffinity);  }
    //! \return the name of this element (in English)
    std::string GetName()      {       return(_name);    }
    //! \return the red component of this element's default visualization color
    double GetRed()            {       return(_red);     }
    //! \return the green component of this element's default color
    double GetGreen()          {       return(_green);   }
    //! \return the blue component of this element's default color
    double GetBlue()           {       return(_blue);    }
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
    int		GetNumberOfElements();

    //! \deprecated Does not properly handle 'D' or 'T' hydrogen isotopes
    int   GetAtomicNum(const char *);
    //! Returns the atomic number matching the element symbol passed
    //! or 0 if not defined. For 'D' or 'T' hydrogen isotopes, will return
    //! a value in the second argument
    int   GetAtomicNum(const char *, int &iso);
    //! \return the element symbol matching the atomic number passed
    char *GetSymbol(int);
    //! \return the van der Waals radius for this atomic number
    double GetVdwRad(int);
    //! \return the covalent radius for this atomic number
    double GetCovalentRad(int);
    //! \return the average atomic mass for this element.
    //! For exact isotope masses, use OpenBabel::OBIsotopeTable
    double GetMass(int);
    //! \return a "corrected" bonding radius based on the hybridization.
    //! Scales the covalent radius by 0.95 for sp2 and 0.90 for sp hybrids
    double CorrectedBondRad(int,int = 3); // atomic #, hybridization
    //! \return a "corrected" vdW radius based on the hybridization.
    //! Scales the van der Waals radius by 0.95 for sp2 and 0.90 for sp hybrids
    double CorrectedVdwRad(int,int = 3); // atomic #, hybridization
    //! \return the maximum expected number of bonds to this element
    int	GetMaxBonds(int);
    //! \return the Pauling electronegativity for this element
    double GetElectroNeg(int);
    //! \return the ionization potential (in eV) for this element
    double GetIonization(int);
    //! \return the electron affinity (in eV) for this element
    double GetElectronAffinity(int);
    //! \return a vector with red, green, blue color values for this element
    std::vector<double> GetRGB(int);
    //! \return the name of this element
    std::string GetName(int);
};

// class introduction in data.cpp
class OBAPI OBIsotopeTable : public OBGlobalDataBase
{
    std::vector<std::vector<std::pair <unsigned int, double> > > _isotopes;

public:

    OBIsotopeTable(void);
    ~OBIsotopeTable()    {}

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
    ~OBTypeTable() {}

    void ParseLine(const char*);

    //! Set the initial atom type to be translated
    bool SetFromType(char*);
    //! Set the destination atom type for translation
    bool SetToType(char*);
    //! Translate atom types
    bool Translate(char *to, char *from); // to, from
    //! Translate atom types
    bool Translate(std::string &to, std::string &from); // to, from

    //! Return the initial atom type to be translated
    std::string GetFromType();
    //! Return the destination atom type for translation
    std::string GetToType();
};

// Used by other code for reading files
#ifdef WIN32
#define FILE_SEP_CHAR "\\"
#else
#define FILE_SEP_CHAR "/"
#endif

} // end namespace OpenBabel

#endif //DATA_H

//! \file data.h
//! \brief Global data and resource file parsers.
