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

#ifndef DATA_H
#define DATA_H

#include <stdio.h>

#ifdef __sgi
#include <iostream.h>
#include <fstream.h>
#else
#include <iostream>
#include <fstream>
#endif

#include <algorithm>
#include <vector>
#include <string>

#ifndef BUFF_SIZE
#define BUFF_SIZE 1024
#endif

#ifdef WIN32
#define FILE_SEP_CHAR "\\"
#else
#define FILE_SEP_CHAR "/"
#endif

namespace OpenBabel {

class OBElement;
class OBElementTable;
typedef enum {UNDEFINED,SDF,MOL2,PDB,DELPDB,SMI,BOX,FIX,
	      OEBINARY,CCC,MMD,ALCHEMY,BALLSTICK,BGF,GHEMICAL,
              XYZ,GAMESSIN,GAMESSOUT,HIN,CACAO,CACAOINT,CACHE, 
	      CHEMDRAW,CML,CSR,CSSR,FEATURE,FENSKEHALL,GROMOS96A,
	      GROMOS96N,QCHEMIN,MPQC,PREP,BIOSYM,CADPAC,CHEM3D1,
	      CHEM3D2,FDAT,GSTAT,DOCK,FRACT,M3D,GAUSSIANZ,
	      GAUSSIANCART,GAUSSIAN92,GAUSSIAN94,MACMOL,MICROWORLD,
	      MM2IN,MM2OUT,MM3,MMADS,MOLIN,MOLINVENT,MOPACCART,
	      MOPACINT,MOPACOUT,PCMODEL,JAGUARIN,JAGUAROUT,
	      REPORT,MSF,SCHAKAL,SHELX,SPARTAN,
	      SPARTANSEMI,SPARTANMM,UNICHEM,XED,BMIN,ICON8,IDATM,
	      MACCS,TINKER,CHARMM,QCHEMOUT,TITLE,DMOL,
	      NWCHEMIN,NWCHEMOUT,RDF,SMIRKS,VIEWMOL } io_type;

class OBGlobalDataBase
{
 protected:
  bool         _init;
  const char  *_dataptr;
  std::string  _filename;
  std::string  _dir;
  std::string  _subdir;
  std::string  _envvar;
 public:
  OBGlobalDataBase()
    {
      _init = false;
      _dataptr = (char*)NULL;
    }
  virtual ~OBGlobalDataBase() {}
  void  Init();
  void  SetReadDirectory(char *dir)       {_dir = dir;}
  void  SetEnvironmentVariable(char *var) {_envvar = var;}
  virtual void ParseLine(const char*) {}
};

class OBElement
{
  int _num;
  char _symbol[3];
  float _Rcov,_Rbo,_Rvdw,_mass,_elNeg;
  int _maxbonds;
 public:
  OBElement() {}
  OBElement(int num,char *sym,float rcov,float rbo,float rvdw,int maxbo,float mass,float elNeg)
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
  float GetCovalentRad() {return(_Rcov);}
  float GetBoRad() {return(_Rbo);}
  float GetVdwRad() {return(_Rvdw);}
  float GetMass() {return(_mass);}
  int GetMaxBonds() {return(_maxbonds);}
  float GetElectroNeg() {return(_elNeg);}
};

class OBElementTable : public OBGlobalDataBase
{
  std::vector<OBElement*> _element;

public:

  OBElementTable(void);
  ~OBElementTable();

  int   GetAtomicNum(const char *);
  void  ParseLine(const char*);
  char *GetSymbol(int);
  float GetVdwRad(int);
  float GetCovalentRad(int);
  float GetBORad(int);
  float GetMass(int);
  float CorrectedBondRad(int,int = 3); // atomic #, hybridization
  float CorrectedVdwRad(int,int = 3); // atomic #, hybridization
  int	GetMaxBonds(int);
  float GetElectroNeg(int);
};

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
  bool	  IsWritable(unsigned int);
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

}

#endif //DATA_H
