/**********************************************************************
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2002 by Geoffrey R. Hutchison

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

#ifndef OB_ROTOR_H
#define OB_ROTOR_H

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
  std::string           _s;
  OBSmartsPattern* _sp;
  std::vector<float>    _vals;
 public:
  OBRotorRule(char *buffer,int ref[4],std::vector<float> &vals,float d)
    {
      _s = buffer;
      _sp = new OBSmartsPattern;
      _sp->Init(buffer);
      memcpy(_ref,ref,sizeof(int)*4);
      _vals = vals;
      _delta = d;
    }
  ~OBRotorRule() {if (_sp) {delete _sp; _sp = NULL;}}

  bool             IsValid() {return(_sp->IsValid());}
  void             GetReferenceAtoms(int ref[4])  
    {memcpy(ref,_ref,sizeof(int)*4);}
  void             SetDelta(float d) {_delta = d;}
  float            GetDelta()        {return(_delta);}
  std::string          &GetSmartsString() {return(_s);}
  std::vector<float>   &GetTorsionVals() {return(_vals);}
  OBSmartsPattern *GetSmartsPattern() {return(_sp);}
};

class OBRotorRules : public OBGlobalDataBase
{
  bool                 _quiet;
  std::vector<OBRotorRule*> _vr;
  std::vector<float>        _sp3sp3;
  std::vector<float>        _sp3sp2;
  std::vector<float>        _sp2sp2;
  public:
  OBRotorRules();
  ~OBRotorRules();

  void ParseLine(const char*);
  void SetFilename(std::string &s) {_filename = s;}
  void GetRotorIncrements(OBMol&,OBBond*,int [4],std::vector<float>&,float &delta);
  void Quiet() {_quiet=true;};
};

class OBRotor
{
  int _idx,_ref[4];
  int *_rotatoms,_size,_numcoords;
  float _delta;
  float _imag,_refang;
  OBBond *_bond;
  std::vector<int> _torsion;
  OBBitVec _fixedatoms,_evalatoms;
  std::vector<float> _res;  //!< torsion resolution
  std::vector<float> _invmag;
  std::vector<std::vector<float> > _sn,_cs,_t;
 public:
  OBRotor();
  ~OBRotor() {if (_rotatoms) delete [] _rotatoms;}
  int     Size()                             {return((_res.empty())?0:_res.size());}
  int     GetIdx() const                     {return(_idx);}
  void    SetNumCoords(int nc)               {_numcoords = nc;}
  void    SetBond(OBBond *bond)              {_bond = bond;}
  void    SetEvalAtoms(OBBitVec &bv)         {_evalatoms = bv;}
  void    SetDihedralAtoms(std::vector<int> &vi)  {_torsion = vi;}
  void    SetDelta(float d)                  {_delta = d;}
  void    SetDihedralAtoms(int ref[4]);
  void    SetRotAtoms(std::vector<int>&);
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
  void    Precalc(std::vector<float*>&);
  void    SetIdx(int idx) {_idx = idx;}
  void    SetFixedAtoms(OBBitVec &bv)          {_fixedatoms = bv;}
  void    SetTorsionValues(std::vector<float> &tmp) {_res = tmp;}
  void    RemoveSymTorsionValues(int);
  void    GetDihedralAtoms(int ref[4]) {for (int i=0;i<4;i++)ref[i]=_ref[i];}
  void    *GetRotAtoms() {return(_rotatoms);}
  float   CalcTorsion(float *);
  float   CalcBondLength(float*);
  float   GetDelta() {return(_delta);}
  OBBond *GetBond()  {return(_bond);}
  std::vector<int> &GetDihedralAtoms()                {return(_torsion);}
  std::vector<float> &GetResolution()                 {return(_res);}
  std::vector<float>::iterator BeginTorIncrement()    {return(_res.begin());}
  std::vector<float>::iterator EndTorIncrement()      {return(_res.end());}
  OBBitVec &GetEvalAtoms() {return(_evalatoms);}
  OBBitVec &GetFixedAtoms() {return(_fixedatoms);}
};

class OBRotorList
{
  bool _quiet,_removesym;
  OBBitVec _fix;
  OBRotorRules _rr;
  std::vector<int> _dffv;         //!< distance from fixed
  std::vector<OBRotor*> _rotor;
  std::vector<std::pair<OBSmartsPattern*,std::pair<int,int> > > _vsym2;
  std::vector<std::pair<OBSmartsPattern*,std::pair<int,int> > > _vsym3;
 public:
  OBRotorList();
  ~OBRotorList();

  int    Size()
    {return((_rotor.empty()) ? 0: _rotor.size());}
  void   Init(std::string &fname)
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
  OBRotor *BeginRotor(std::vector<OBRotor*>::iterator &i) 
    {i = _rotor.begin();return((i ==_rotor.end()) ? NULL:*i);}
  OBRotor *NextRotor(std::vector<OBRotor*>::iterator &i) 
      {i++;return((i ==_rotor.end()) ? NULL:*i);}
  std::vector<OBRotor*>::iterator BeginRotors() {return(_rotor.begin());}
  std::vector<OBRotor*>::iterator EndRotors() {return(_rotor.end());}
};

} // end namespace OpenBabel

#ifndef SQUARE
#define SQUARE(x) ((x)*(x))
#endif

#endif // OB_ROTOR_H

