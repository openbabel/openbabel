/**********************************************************************
data.h - Global data and resource file parsers.

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

#include "babelconfig.h"

#ifndef OB_DATA_H
#define OB_DATA_H

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

namespace OpenBabel {

class OBElement;
class OBAtom;
class OBElementTable;

typedef enum { UNDEFINED,
               ALCHEMY, BALLSTICK, BGF, BIOSYM, BMIN, BOX, CACAO,
               CACAOINT, CACHE, CADPAC, CCC, CDX, CHARMM, CHEM3D1,
	       CHEM3D2, CHEMDRAW, CHEMTOOL, CIF, CML, CSR, CSSR, DELPDB, DMOL,
	       DOCK, FDAT, FEATURE, FH, FIX, FRACT, GAMESSIN, GAMESSOUT,
               GAUSSIAN92, GAUSSIAN94, GAUSSIANCART, GAUSSIANOUT,
               GHEMICAL, GROMOS96A, GROMOS96N, GSTAT, HIN, ICON8,
               IDATM, JAGUARIN, JAGUAROUT, M3D, MACCS, MACMOL,
               MICROWORLD, MM2IN, MM2OUT, MM3, MMADS, MMCIF, MMD,
               MOL2, MOLDEN, MOLIN, MOLINVENT, MOPACCART, MOPACINT,
               MOPACOUT, MPQC, MSF, NWCHEMIN, NWCHEMOUT, OEBINARY,
               PCMODEL, PDB, POV, PREP, QCHEMIN, QCHEMOUT, REPORT,
               SCHAKAL, SDF, SHELX, SKC, SMI, SPARTAN, SPARTANMM,
               SPARTANSEMI, TGF, TINKER, TITLE, TURBOMOLE, UNICHEM, VIEWMOL,
               XED, XYZ, ZINDO, CRK2D, CRK3D, PQS
	       // Insert new formats here (at the end)
	       // for backwards compatibility
             } io_type;

//! \brief Base data table class, handles reading data files
//!
//! Base data table class--reads ASCII data files in various formats
//! -# Checks for the environment variable _envvar (defaults to "BABEL_DATADIR")
//!     Tries that directory as well as the _subdir directory of that (def. "data")
//! -# Checks for the directory _dir (def. determined by the build environment)
//! -# Reverts to the compiled-in default data
class OBGlobalDataBase
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
  virtual ~OBGlobalDataBase() {}
  //! Read in the data file, falling back as needed
  void  Init();
  //! Set the directory before calling Init()
  void  SetReadDirectory(char *dir)       {_dir = dir;}
  //! Set the environment variable to use before calling Init()
  void  SetEnvironmentVariable(char *var) {_envvar = var;}
  //! Specified by particular table classes (parses an individual data line)
  virtual void ParseLine(const char*) {}
};

//! \brief Individual element data type
//!
//! Stores a variety of data about an individual element
class OBElement
{
  int _num;
  char _symbol[3];
  double _Rcov,_Rbo,_Rvdw,_mass,_elNeg;
  int _maxbonds;
 public:
  OBElement() {}
  OBElement(int num, const char *sym, double rcov, double rbo, 
	    double rvdw, int maxbo, double mass, double elNeg)
    {
      _num = num;
      strcpy(_symbol,sym);
      _Rcov = rcov;
      _Rbo = rbo;
      _Rvdw = rvdw;
      _maxbonds = maxbo;
      _mass = mass;
      _elNeg = elNeg;
    }
  int GetAtomicNum() {return(_num);}
  char *GetSymbol() {return(_symbol);}
  double GetCovalentRad() {return(_Rcov);}
  double GetBoRad() {return(_Rbo);}
  double GetVdwRad() {return(_Rvdw);}
  double GetMass() {return(_mass);}
  int GetMaxBonds() {return(_maxbonds);}
  double GetElectroNeg() {return(_elNeg);}
};

// class introduction in data.cpp
class OBElementTable : public OBGlobalDataBase
{
  std::vector<OBElement*> _element;

public:

  OBElementTable(void);
  ~OBElementTable();

  int   GetAtomicNum(const char *, unsigned short int iso = 0);
  void  ParseLine(const char*);
  char *GetSymbol(int);
  double GetVdwRad(int);
  double GetCovalentRad(int);
  double GetBORad(int);
  double GetMass(int);
  double CorrectedBondRad(int,int = 3); // atomic #, hybridization
  double CorrectedVdwRad(int,int = 3); // atomic #, hybridization
  int	GetMaxBonds(int);
  double GetElectroNeg(int);
};

// class introduction in data.cpp
class OBIsotopeTable : public OBGlobalDataBase
{
  std::vector<std::vector<std::pair <unsigned int, double> > > _isotopes;

 public:
  
  OBIsotopeTable(void);
  ~OBIsotopeTable() {}

  void	ParseLine(const char*);
  //! Return the exact masss of the isotope
  //!   (or by default (i.e. "isotope 0") the most abundant isotope)
  double	GetExactMass(const unsigned int atomicNum,
			     const unsigned int isotope = 0);
};

// class introduction in data.cpp
class OBTypeTable : public OBGlobalDataBase
{
  int    _linecount;
  int    _ncols,_nrows,_from,_to;
  std::vector<std::string> _colnames;
  std::vector<std::vector<std::string> > _table;

 public:

  OBTypeTable(void);
  ~OBTypeTable() {}

  void ParseLine(const char*);
  bool SetFromType(char*);
  bool SetToType(char*);
  bool Translate(char*,char*); // to, from
  bool Translate(std::string &,std::string &); // to, from
};

// class introduction in data.cpp
class OBExtensionTable : public OBGlobalDataBase
{
  int                     _linecount;
  std::vector<std::vector<std::string> > _table;

 public:

  OBExtensionTable(void);
  ~OBExtensionTable() {}

  bool    CanReadExtension(char *);
  bool    CanWriteExtension(char *);
  bool	  IsReadable(unsigned int);
  bool    IsReadable(io_type);
  bool	  IsWritable(unsigned int);
  bool	  IsWritable(io_type);
  void    ParseLine(const char*);
  void    TypeToExtension(io_type,char*);
  void	  TypeToMIME(io_type,char*);
  void    ExtensionToDescription(char*, char*);

  io_type       GetType(unsigned int);
  io_type       FilenameToType(char *);
  io_type       FilenameToType(std::string &);
  io_type	MIMEToType(char *);
  io_type	MIMEToType(std::string &);
  const char   *GetExtension(unsigned int);
  const char   *GetDescription(unsigned int);
  unsigned int  Count(); 
};

// Used by other code for reading files
#ifdef WIN32
#define FILE_SEP_CHAR "\\"
#else
#define FILE_SEP_CHAR "/"
#endif

}

#endif //DATA_H
