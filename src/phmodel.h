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

#ifndef _PHMODEL_H
#define _PHMODEL_H

#include "parsmart.h"
#include "data.h"

namespace OpenBabel {

class OBChemTsfm
{
  vector<int>                       _vadel;
  vector<pair<int,int> >            _vele;
  vector<pair<int,int> >            _vchrg;
  vector<pair<int,int> >            _vbdel;
  vector<pair<pair<int,int>,int> >  _vbond;
  OBSmartsPattern _bgn,_end;
public:
  OBChemTsfm() {}
  ~OBChemTsfm() {}
  bool Init(string&,string&);
  bool Apply(OBMol&);
};


class OBPhModel : public OBGlobalDataBase
{
  vector<vector<int> >                           _mlist;
  vector<OBChemTsfm*>                            _vtsfm;
  vector<pair<OBSmartsPattern*,vector<float> > > _vschrg;
 public:
  OBPhModel();
  ~OBPhModel();

  void ParseLine(char*);
  void AssignSeedPartialCharge(OBMol&);
  void CorrectForPH(OBMol&);
};



} //namespace OpenBabel

#endif //_PHMODEL_H
