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

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include "data.h"
#include "element.h"
#include "types.h"
#include "extable.h"

#ifdef __sgi
#include <strstream.h>
#else
#include <strstream>
#endif

namespace OpenEye {

OEExtensionTable extab;
OEElementTable   etab;
OETypeTable      ttab;

bool tokenize(vector<string>&, char *buf,char *delimstr=" \t\n");
bool tokenize(vector<string> &vcr, string &s,char *delimstr,int limit=-1);
extern void ThrowError(char *);
extern void ThrowError(string&);

OEElementTable::OEElementTable()
{
  _init = false;
  _dir = "";
  _envvar = "OE_DIR";
  _filename = "element.txt";
  _subdir = "data";
  _dataptr = ElementData;
}

OEElementTable::~OEElementTable()
  {
    vector<OEElement*>::iterator i;
    for (i = _element.begin();i != _element.end();i++) delete *i;
  }

void OEElementTable::ParseLine(char *buffer)
{
  int num,maxbonds;
  char symbol[3];
  float Rbo,Rcov,Rvdw,mass;

  if (buffer[0] != '#') // skip comment line (at the top)
    {
      // Ignore RGB columns
      sscanf(buffer,"%d %s %f %f %f %d %f %*f %*f %*f",
	     &num,
	     symbol,
	     &Rcov,
	     &Rbo,
	     &Rvdw,
	     &maxbonds,
	     &mass);
  
	  OEElement *ele = new OEElement(num,symbol,Rcov,Rbo,Rvdw,maxbonds,mass);
	  _element.push_back(ele);
	}
}

char *OEElementTable::GetSymbol(int atomicnum)
{
  if (!_init) Init();

  if (atomicnum < 0 || atomicnum > _element.size())
    return("\0");

  return(_element[atomicnum]->GetSymbol());
}

int OEElementTable::GetMaxBonds(int atomicnum)
{
  if (!_init) Init();

  if (atomicnum < 0 || atomicnum > _element.size())
    return(0);

  return(_element[atomicnum]->GetMaxBonds());
}

float OEElementTable::GetVdwRad(int atomicnum)
{
  if (!_init) Init();

  if (atomicnum < 0 || atomicnum > _element.size())
    return(0.0f);

  return(_element[atomicnum]->GetVdwRad());
}

float OEElementTable::GetBORad(int atomicnum)
{
  if (!_init) Init();

  if (atomicnum < 0 || atomicnum > _element.size())
    return(0.0f);

  return(_element[atomicnum]->GetBoRad());
}

float OEElementTable::CorrectedBondRad(int atomicnum, int hyb)
{
  float rad;
  if (!_init) Init();

  if (atomicnum < 0 || atomicnum > _element.size())
    return(1.0f);

  rad = _element[atomicnum]->GetBoRad();

  if (hyb == 2)      rad *= 0.95f;
  else if (hyb == 1) rad *= 0.90f;

  return(rad);
}

float OEElementTable::CorrectedVdwRad(int atomicnum, int hyb)
{
  float rad;
  if (!_init) Init();

  if (atomicnum < 0 || atomicnum > _element.size())
    return(1.95f);

  rad = _element[atomicnum]->GetVdwRad();

  if (hyb == 2)      rad *= 0.95f;
  else if (hyb == 1) rad *= 0.90f;

  return(rad);
}

float OEElementTable::GetCovalentRad(int atomicnum)
{
  if (!_init) Init();

  if (atomicnum < 0 || atomicnum > _element.size())
    return(0.0f);

  return(_element[atomicnum]->GetCovalentRad());
}

float OEElementTable::GetMass(int atomicnum)
{
  if (!_init) Init();

  if (atomicnum < 0 || atomicnum > _element.size())
    return(0.0f);

  return(_element[atomicnum]->GetMass());
}

int OEElementTable::GetAtomicNum(const char *sym)
{
  if (!_init) Init();

    vector<OEElement*>::iterator i;
    for (i = _element.begin();i != _element.end();i++)
    if (!strcmp(sym,(*i)->GetSymbol()))
        return((*i)->GetAtomicNum());

    return(0);
}

OETypeTable::OETypeTable()
{
  _init = false;
  _dir = "";
  _envvar = "OE_DIR";
  _filename = "types.txt";
  _subdir = "data";
  _dataptr = TypesData;
  _linecount = 0;
  _from = _to = -1;
}

void OETypeTable::ParseLine(char *buffer)
{
  if (_linecount == 0)
    sscanf(buffer,"%d%d",&_ncols,&_nrows);
  else if (_linecount == 1)
    tokenize(_colnames,buffer);
  else
    {
      vector<string> vc;
      tokenize(vc,buffer);
      if (vc.size() == (unsigned)_nrows) _table.push_back(vc);
    }
  _linecount++;
}

bool OETypeTable::SetFromType(char* from)
{
  if (!_init) Init();

    string tmp = from;

    unsigned int i;
    for (i = 0;i < _colnames.size();i++)
    if (tmp == _colnames[i])
    {
        _from = i;
        return(true);
    }

    ThrowError("Requested type column not found");
    
    return(false);
}

bool OETypeTable::SetToType(char* to)
{
  if (!_init) Init();

    string tmp = to;

    unsigned int i;
    for (i = 0;i < _colnames.size();i++)
    if (tmp == _colnames[i])
    {
        _to = i;
        return(true);
    }

    ThrowError("Requested type column not found");
    
    return(false);
}

bool OETypeTable::Translate(char *to, char *from)
{
  if (!_init) Init();

  bool rval;
  string sto,sfrom;
  sfrom = from;
  rval = Translate(sto,sfrom);
  strcpy(to,(char*)sto.c_str());

  return(rval);
}

bool OETypeTable::Translate(string &to,string &from)
{
  if (!_init) Init();

  if (from == "")
    return(false);

  vector<vector<string> >::iterator i;
  for (i = _table.begin();i != _table.end();i++)
    if ((signed)(*i).size() > _from &&  (*i)[_from] == from)
      {
	to = (*i)[_to];
	return(true);
      }
  
  to = from;
  return(false);
}

void Toupper(string &s)
{
  unsigned int i;
  for (i = 0;i < s.size();i++)
    s[i] = toupper(s[i]);
}

void Tolower(string &s)
{
  unsigned int i;
  for (i = 0;i < s.size();i++)
    s[i] = tolower(s[i]);
}

OEExtensionTable::OEExtensionTable()
{
  _init = false;
  _dir = "";
  _envvar = "OE_DIR";
  _filename = "extable.txt";
  _subdir = "data";
  _dataptr = ExtensionTableData;
  _linecount = 0;
}

void OEExtensionTable::ParseLine(char *buffer)
{
  if (_linecount > 0)
    {
      vector <string> vs;
      
      if (buffer[0] != '#') // skip comments
	{
	  tokenize(vs,buffer,"\t\n"); // spaces are a problem
	  if (vs.size() == 5)
	    {
	      Toupper(vs[1]);
	      _table.push_back(vector <string> (vs));
	    }
	}
    }

  _linecount++;
}

io_type TextToType(string typestring)
{
  if (typestring == "MOL2")			return(MOL2);
  else if (typestring == "PDB")			return(PDB);
  else if (typestring == "SDF")			return(SDF);
  else if (typestring == "BOX")			return(BOX);
  else if (typestring == "SMI")			return(SMI);
  else if (typestring == "MMD")			return(MMD);
  else if (typestring == "OEBINARY")		return(OEBINARY);
  else if (typestring == "GHEMICAL")		return(GHEMICAL);
  else if (typestring == "XYZ")			return(XYZ);
  else if (typestring == "GAMESSIN")		return(GAMESSIN);
  else if (typestring == "GAMESSOUT")		return(GAMESSOUT);
  else if (typestring == "HIN")			return(HIN);
  else if (typestring == "CCC")			return(CCC);
  else if (typestring == "BALLSTICK")		return(BALLSTICK);
  else if (typestring == "ALCHEMY")		return(ALCHEMY);
  else if (typestring == "BGF")			return(BGF);
  else if (typestring == "FIX")			return(FIX);
  else if (typestring == "CACAO")		return(CACAO);
  else if (typestring == "CACAOINT")		return(CACAOINT);
  else if (typestring == "CACHE")		return(CACHE);
  else if (typestring == "CHEMDRAW")		return(CHEMDRAW);
  else if (typestring == "CSR")			return(CSR);
  else if (typestring == "CSSR")		return(CSSR);
  else if (typestring == "FEATURE")		return(FEATURE);
  else if (typestring == "FH")			return(FENSKEHALL);
  else if (typestring == "GROMOS96A")		return(GROMOS96A);
  else if (typestring == "GROMOS96N")		return(GROMOS96N);
  else if (typestring == "QCHEMIN")		return(QCHEMIN);
  else if (typestring == "QCHEMOUT")		return(QCHEMOUT);
  else if (typestring == "MPQC")		return(MPQC);
  else if (typestring == "UNICHEM")		return(UNICHEM);
  else if (typestring == "TINKER")		return(TINKER);
  else if (typestring == "PREP")		return(PREP);
  else if (typestring == "BIOSYM")		return(BIOSYM);
  else if (typestring == "CADPAC")		return(CADPAC);
  else if (typestring == "CHEM3D1")		return(CHEM3D1);
  else if (typestring == "CHEM3D2")		return(CHEM3D2);
  else if (typestring == "FDAT")		return(FDAT);
  else if (typestring == "GSTAT")		return(GSTAT);
  else if (typestring == "DOCK")		return(DOCK);
  else if (typestring == "FRACT")		return(FRACT);
  else if (typestring == "M3D")			return(M3D);
  else if (typestring == "GAUSSIANZMAT")	return(GAUSSIANZ);
  else if (typestring == "GAUSSIANCART")	return(GAUSSIANCART);
  else if (typestring == "GAUSSIAN92")		return(GAUSSIAN92);
  else if (typestring == "GAUSSIAN94")		return(GAUSSIAN94);
  else if (typestring == "MACMOL")		return(MACMOL);
  else if (typestring == "MICROWORLD")		return(MICROWORLD);
  else if (typestring == "MM2IN")		return(MM2IN);
  else if (typestring == "MM2OUT")		return(MM2OUT);
  else if (typestring == "MM3")			return(MM3);
  else if (typestring == "MMADS")		return(MMADS);
  else if (typestring == "MOLIN")		return(MOLIN);
  else if (typestring == "MOLINVENT")		return(MOLINVENT);
  else if (typestring == "MOPACCART")		return(MOPACCART);
  else if (typestring == "MOPACINT")		return(MOPACINT);
  else if (typestring == "MOPACOUT")		return(MOPACOUT);
  else if (typestring == "PCMODEL")		return(PCMODEL);
  else if (typestring == "JAGUARIN")		return(JAGUARIN);
  else if (typestring == "JAGUAROUT")		return(JAGUAROUT);
  else if (typestring == "REPORT")		return(REPORT);
  else if (typestring == "MSF")			return(MSF);
  else if (typestring == "SCHAKAL")		return(SCHAKAL);
  else if (typestring == "SHELX")		return(SHELX);
  else if (typestring == "SPARTAN")		return(SPARTAN);
  else if (typestring == "SPARTANSEMI")		return(SPARTANSEMI);
  else if (typestring == "SPARTANMM")		return(SPARTANMM);
  else if (typestring == "XED")			return(XED);
  else if (typestring == "BMIN")		return(BMIN);
  else if (typestring == "ICON8")		return(ICON8);
  else if (typestring == "IDATM")		return(IDATM);
  else if (typestring == "MACCS")		return(MACCS);
  else if (typestring == "CHARMM")		return(CHARMM);
  else if (typestring == "DMOL")		return(DMOL);
  else if (typestring == "NWCHEMIN")		return(NWCHEMIN);
  else if (typestring == "NWCHEMOUT")		return(NWCHEMOUT);
  else if (typestring == "TITLE")		return(TITLE);
  // Add yours here
  else						return(UNDEFINED);

}

io_type OEExtensionTable::FilenameToType(char *filename)
{
  if (!_init) Init();

  vector<vector<string> >::iterator i;

  vector<string> vs;
  tokenize(vs,filename,".\n\t");
  if (vs.empty()) return(UNDEFINED);

  string ext = vs[vs.size()-1];
  Tolower(ext);

  io_type type = UNDEFINED;
  for (i = _table.begin();i != _table.end();i++)
    if ((*i)[0] == ext)
    {
      type = TextToType((*i)[1]);
      break;
    }

  return(type);
}

io_type OEExtensionTable::FilenameToType(string &filename)
{
  if (!_init) Init();

  vector<vector<string> >::iterator i;

  vector<string> vs;
  tokenize(vs,filename,".\n\t");
  if (vs.empty()) return(UNDEFINED);

  string ext = vs[vs.size()-1];
  Tolower(ext);

  io_type type = UNDEFINED;
  for (i = _table.begin();i != _table.end();i++)
    if ((*i)[0] == ext)
    {
      type = TextToType((*i)[1]);
      break;
    }

  return(type);
}

void OEExtensionTable::TypeToExtension(io_type type,char *ext)
{
  if (!_init) Init();

  vector<vector<string> >::iterator i;
  strcpy(ext,"ext");

  for (i = _table.begin();i != _table.end();i++)
    if (type == TextToType((*i)[1]))
    {
      strcpy(ext,(char*)(*i)[1].c_str());
      break;
    }
}

void OEExtensionTable::ExtensionToDescription(char *filename, char *desc)
{
  if (!_init) Init();

  vector<vector<string> >::iterator i;

  vector<string> vs;
  tokenize(vs,filename,".\n\t");
  if (vs.empty()) return;

  string ext = vs[vs.size()-1];

  for (i = _table.begin();i != _table.end();i++)
    if ((*i)[0] == ext)
      {
	strcpy(desc, (char*)(*i)[2].c_str()); break;
      }

  return;
}

bool OEExtensionTable::CanReadExtension(char *filename)
{
  if (!_init) Init();

  vector<vector<string> >::iterator i;

  vector<string> vs;
  tokenize(vs,filename,".\n\t");
  if (vs.empty()) return(false);

  string ext = vs[vs.size()-1];
  bool read = false;
  for (i = _table.begin();i != _table.end();i++)
    if ((*i)[0] == ext && (*i)[3] == "1")
      {
	read = true;
	break;
      }
  return read;
}

bool OEExtensionTable::CanWriteExtension(char *filename)
{
  if (!_init) Init();

  vector<vector<string> >::iterator i;

  vector<string> vs;
  tokenize(vs,filename,".\n\t");
  if (vs.empty()) return(false);

  string ext = vs[vs.size()-1];
  bool write = false;
  for (i = _table.begin();i != _table.end();i++)
    if ((*i)[0] == ext && (*i)[4] == "1")
      {
	write = true;
	break;
      }
  return write;
}

const char *OEExtensionTable::GetExtension(unsigned int n)
{
  if (!_init) Init();

  if (n >= _table.size())
    return NULL;
  else
    {
      ostrstream longDesc; // Need to null-terminate
      longDesc << _table[n][0] << ends;
      return(longDesc.str());
    }
}

const char *OEExtensionTable::GetDescription(unsigned int n)
{
  if (!_init) Init();

  if (n >= _table.size())
    return NULL;
  else
    {
      ostrstream longDesc; // Need to null-terminate
      longDesc << _table[n][2] << ends;
      return(longDesc.str());
    }
}

io_type OEExtensionTable::GetType(unsigned int n)
{
  if (!_init) Init();

  if (n >= _table.size())
    return UNDEFINED;
  else
    {
      // A bit of a pain to discard the const qualifier...
      char *temp = new char[_table[n][0].length()];
      io_type returnVal;

      strcpy(temp, (char*)_table[n][0].c_str());
      returnVal = FilenameToType(temp);
      delete [] temp;
      return returnVal;
    }
}
 
bool OEExtensionTable::IsReadable(unsigned int n)
{
  if (!_init) Init();

  if (n >= _table.size())
    return false;
  else
    return _table[n][3] == "1";
}

bool OEExtensionTable::IsWritable(unsigned int n)
{
  if (!_init) Init();

  if (n >= _table.size())
    return false;
  else
    return _table[n][4] == "1";
}

unsigned int OEExtensionTable::Count()
{
  if (!_init) Init(); 
  return(_table.size());
}

void OEGlobalDataBase::Init()
{
  if (_init) return;
  _init = true;

  char buffer[BUFF_SIZE],subbuffer[BUFF_SIZE];
  ifstream ifs, ifs1, ifs2, ifs3, *ifsP;
  ifs.open((char*)_filename.c_str());
  ifsP= &ifs;

  if (!(*ifsP) && !_dir.empty())
    {
      strcpy(buffer,_dir.c_str());
      strcat(buffer,FILE_SEP_CHAR);
      strcat(buffer,(char*)_filename.c_str());
      ifs1.open(buffer);     // weirdness have to use ifs1 and not ifs
      ifsP = &ifs1;            // set the pointer so we can use it
    }
   
  if (!(*ifsP) && getenv(_envvar.c_str()) != NULL)
  {
    strcpy(buffer,getenv(_envvar.c_str()));
    strcat(buffer,FILE_SEP_CHAR);

    if (!_subdir.empty())
    {
	  strcpy(subbuffer,buffer);
      strcat(subbuffer,_subdir.c_str());
      strcat(subbuffer,FILE_SEP_CHAR);
    }

    strcat(buffer,(char*)_filename.c_str());
	strcat(subbuffer,(char*)_filename.c_str());

	ifs2.open(subbuffer);
	ifsP= &ifs2;
	if (!(*ifsP))
	{
		ifs3.open(buffer);
		ifsP= &ifs3;
	}
  }

  if ((*ifsP))
    {
      for (;ifsP->getline(buffer,BUFF_SIZE);)
		  ParseLine(buffer);
    }
  else
    if (_dataptr)
    {
      char *p1,*p2;
      for (p1 = p2 = _dataptr;*p2 != '\0';p2++)
	if (*p2 == '\n')
	  {
	    *p2 = '\0';
	    ParseLine(p1);
	    p1 = ++p2;
	  }
    }
  else
    {
      string s = "Unable to open data file '"; s += _filename; s += "'";
      ThrowError(s);
    }

  if (ifs)  ifs.close();
  if (ifs1) ifs1.close();
  if (ifs2) ifs2.close();
}

}

