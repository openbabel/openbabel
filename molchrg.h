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

#ifndef __MOLCHRG_H
#define __MOLCHRG_H

namespace OpenEye {

class GasteigerState;

class OEGastChrg
{
  vector <GasteigerState*> _gsv;
  void InitialPartialCharges(OEMol &);
  bool GasteigerSigmaChi(OEAtom *,float &,float &,float &);
 public:
  OEGastChrg(){}
  ~OEGastChrg();
  bool AssignPartialCharges(OEMol &);
  void GSVResize(int);
};

class GasteigerState 
{
  //helper class for OEGastChrg
 public:
  GasteigerState();
  ~GasteigerState() {}
  void SetValues(float _a,float _b,float _c,float _q) 
    {a = _a;b = _b;c = _c;denom=a+b+c;q = _q;}
  float a, b, c;
  float denom;
  float chi;
  float q;
};

}

#define MX_GASTEIGER_DENOM  20.02
#define MX_GASTEIGER_DAMP   0.5
#define MX_GASTEIGER_ITERS  6


#endif //__MOLCHRG_H
