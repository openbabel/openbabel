/**********************************************************************
Copyright (C) 2003 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "bondtyper.h"
#include "bondtyp.h"

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

using namespace std;

namespace OpenBabel {

OBBondTyper  bondtyper;

  /*! \class OBBondTyper
      \brief Assigns bond types for file formats without bond information

  The OBBondTyper class is designed to read in a list of bond typing
  rules and apply them to molecules.
*/
OBBondTyper::OBBondTyper()
{
  _init = false;
  _dir = DATADIR;
  _envvar = "BABEL_DATADIR";
  _filename = "bondtyp.txt";
  _subdir = "data";
  _dataptr = BondTypeData;
}

void OBBondTyper::ParseLine(const char *buffer)
{
  vector<string> vs;
  OBSmartsPattern *sp;

  if (EQn(buffer,"INTHYB",6))
    {
      tokenize(vs,buffer);
      if (vs.empty() || vs.size() < 3) return;
      sp = new OBSmartsPattern;
      if (sp->Init(vs[1]))
	_vinthyb.push_back(pair<OBSmartsPattern*,int> (sp,atoi((char*)vs[2].c_str())));
      else {delete sp; sp = NULL;}
    }
  else if (EQn(buffer,"IMPVAL",6))
    {
      tokenize(vs,buffer);
      if (vs.empty() || vs.size() < 3) return;
      sp = new OBSmartsPattern;
      if (sp->Init(vs[1]))
	_vimpval.push_back(pair<OBSmartsPattern*,int> (sp,atoi((char*)vs[2].c_str())));
      else {delete sp; sp = NULL;}
    }
  else if (EQn(buffer,"EXTTYP",6))
    {
      tokenize(vs,buffer);
      if (vs.empty() || vs.size() < 3) return;
      sp = new OBSmartsPattern;
      if (sp->Init(vs[1]))
	_vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
      else {delete sp; sp = NULL;}
    }
}

OBBondTyper::~OBBondTyper()
{
  vector<pair<OBSmartsPattern*,int> >::iterator i;
  for (i = _vinthyb.begin();i != _vinthyb.end();i++) {delete i->first; i->first = NULL;}
  for (i = _vimpval.begin();i != _vimpval.end();i++) {delete i->first; i->first = NULL;}

  vector<pair<OBSmartsPattern*,string> >::iterator j;
  for (j = _vexttyp.begin();j != _vexttyp.end();j++) {delete j->first; j->first = NULL;}
}

void OBBondTyper::ConnectTheDots(OBMol &mol)
{
}

} //namespace OpenBabel;


