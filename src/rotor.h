/**********************************************************************
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef ROTOR_H
#define ROTOR_H

#include "parsmart.h"
#include "typer.h" 

namespace OpenBabel {

class OBRotor;
class OBRotorList;
class OBRotorRule;
class OBRotorRules;

class OBRotorRule
{
  int              _ref[4];
  float            _delta;
  string           _s;
  OBSmartsPattern* _sp;
  vector<float>    _vals;
 public:
  OBRotorRule(char *buffer,int ref[4],vector<float> &vals,float d)
    {
      _s = buffer;
      _sp = new OBSmartsPattern;
      _sp->Init(buffer);
      memcpy(_ref,ref,sizeof(int)*4);
      _vals = vals;
      _delta = d;
    }
  ~OBRotorRule() {if (_sp) delete _sp;}

  bool             IsValid() {return(_sp->IsValid());}
  void             GetReferenceAtoms(int ref[4])  
    {memcpy(ref,_ref,sizeof(int)*4);}
  void             SetDelta(float d) {_delta = d;}
  float            GetDelta()        {return(_delta);}
  string          &GetSmartsString() {return(_s);}
  vector<float>   &GetTorsionVals() {return(_vals);}
  OBSmartsPattern *GetSmartsPattern() {return(_sp);}
};

class OBRotorRules : public OBGlobalDataBase
{
  bool                 _quiet;
  vector<OBRotorRule*> _vr;
  vector<float>        _sp3sp3;
  vector<float>        _sp3sp2;
  vector<float>        _sp2sp2;
  public:
  OBRotorRules();
  ~OBRotorRules();

  void ParseLine(char*);
  void SetFilename(string &s) {_filename = s;}
  void GetRotorIncrements(OBMol&,OBBond*,int [4],vector<float>&,float &delta);
  void Quiet() {_quiet=true;};
};

class OBRotor
{
  int _idx,_ref[4];
  int *_rotatoms,_size,_numcoords;
  float _delta;
  float _imag,_refang;
  OBBond *_bond;
  vector<int> _torsion;
  OBBitVec _fixedatoms,_evalatoms;
  vector<float> _res;  //torsion resolution
  vector<float> _invmag;
  vector<vector<float> > _sn,_cs,_t;
 public:
  OBRotor();
  ~OBRotor() {if (_rotatoms) delete [] _rotatoms;}
  int     Size()                             {return((_res.empty())?0:_res.size());}
  int     GetIdx() const                     {return(_idx);}
  void    SetNumCoords(int nc)               {_numcoords = nc;}
  void    SetBond(OBBond *bond)              {_bond = bond;}
  void    SetEvalAtoms(OBBitVec &bv)         {_evalatoms = bv;}
  void    SetDihedralAtoms(vector<int> &vi)  {_torsion = vi;}
  void    SetDelta(float d)                  {_delta = d;}
  void    SetDihedralAtoms(int ref[4]);
  void    SetRotAtoms(vector<int>&);
  inline void SetToAngle(float *c,float setang)
    {
      float dx,dy,dz,sn,cs,t,ang,mag;
      ang = setang - CalcTorsion(c);
      if (fabs(ang) < 0.0001) return;
      
      sn = sin(ang); cs = cos(ang);t = 1 - cs;
      dx = c[_torsion[1]]   - c[_torsion[2]];
      dy = c[_torsion[1]+1] - c[_torsion[2]+1];
      dz = c[_torsion[1]+2] - c[_torsion[2]+2];
      mag = sqrt(SQUARE(dx) + SQUARE(dy) + SQUARE(dz));
      Set(c,sn,cs,t,1.0/mag);
    }
  void    SetRotor(float *,int,int prev=-1);
  void    Set(float*,int);
  void    Precompute(float*);
  void    Set(float *c,int ridx,int cidx)
         {Set(c,_sn[cidx][ridx],_cs[cidx][ridx],_t[cidx][ridx],_invmag[cidx]);}
  void    Set(float*,float,float,float,float);
  void    Precalc(vector<float*>&);
  void    SetIdx(int idx) {_idx = idx;}
  void    SetFixedAtoms(OBBitVec &bv)          {_fixedatoms = bv;}
  void    SetTorsionValues(vector<float> &tmp) {_res = tmp;}
  void    RemoveSymTorsionValues(int);
  void    GetDihedralAtoms(int ref[4]) {for (int i=0;i<4;i++)ref[i]=_ref[i];}
  void    *GetRotAtoms() {return(_rotatoms);}
  float   CalcTorsion(float *);
  float   CalcBondLength(float*);
  float   GetDelta() {return(_delta);}
  OBBond *GetBond()  {return(_bond);}
  vector<int> &GetDihedralAtoms()                {return(_torsion);}
  vector<float> &GetResolution()                 {return(_res);}
  vector<float>::iterator BeginTorIncrement()    {return(_res.begin());}
  vector<float>::iterator EndTorIncrement()      {return(_res.end());}
  OBBitVec &GetEvalAtoms() {return(_evalatoms);}
  OBBitVec &GetFixedAtoms() {return(_fixedatoms);}
};

class OBRotorList
{
  bool _quiet,_removesym;
  OBBitVec _fix;
  OBRotorRules _rr;
  vector<int> _dffv;         //distance from fixed
  vector<OBRotor*> _rotor;
  vector<pair<OBSmartsPattern*,pair<int,int> > > _vsym2;
  vector<pair<OBSmartsPattern*,pair<int,int> > > _vsym3;
 public:
  OBRotorList();
  ~OBRotorList();

  int    Size()
    {return((_rotor.empty()) ? 0: _rotor.size());}
  void   Init(string &fname)
    {
      _rr.SetFilename(fname);
      _rr.Init();
    }
  void   SetFixAtoms(OBBitVec &fix)   {_fix = fix;}
  void   ShutTheHellUp()              {_quiet=true;_rr.Quiet();}
  void   NoSymmetryCondoms()          {_removesym=false;}
  void   Clear();
  void   RemoveSymVals(OBMol&);
  void   SetRotAtomsByFix(OBMol&);
  bool   SetRotAtoms(OBMol&);
  bool   Setup(OBMol &);
  bool   FindRotors(OBMol &);
  bool   IdentifyEvalAtoms(OBMol &);
  bool   SetEvalAtoms(OBMol&);
  bool   AssignTorVals(OBMol &);
  bool   IsFixedBond(OBBond*);
  bool   HasFixedAtoms()              {return(!_fix.Empty());}
  OBRotor *BeginRotor(vector<OBRotor*>::iterator &i) 
    {i = _rotor.begin();return((i ==_rotor.end()) ? NULL:*i);}
  OBRotor *NextRotor(vector<OBRotor*>::iterator &i) 
      {i++;return((i ==_rotor.end()) ? NULL:*i);}
  vector<OBRotor*>::iterator BeginRotors() {return(_rotor.begin());}
  vector<OBRotor*>::iterator EndRotors() {return(_rotor.end());}
};


}

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

#endif //ROTOR_H

