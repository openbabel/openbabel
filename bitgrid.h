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

#ifndef __BITGRID_H__
#define __BITGRID_H__

#include "mol.h"
#include "grid.h"
#include "bitvec.h"
#include "oeutil.h"
#include "parsmart.h"
#include "patty.h"

namespace OpenEye {

// Define rint for win32 - this is not defined in Math.h for VC6.0
class BitGrid 
{
protected:

  bool  fuzzy;

  float xmin, ymin, zmin;
  float xmax, ymax, zmax;
  float midx, midy, midz;
  int   xdim, ydim, zdim, xydim;

  float spacing, inv_spa;
  int   size;

  OEBitVec grid,lipo,don,acc;

  patty p;
  vector<string> types;

public:
  
  BitGrid(void);
  BitGrid(bool);
  ~BitGrid(void);

  void Init(OEMol &, float);
  void Init(float,float,float,float,float,float,float);
  void Build(OEMol &);
  void Build(OEMol &, OEBitVec &);
  void Build(OEMol &, vector<int> &);
  void SetBits(OEAtom *);

  void Clear(void) { grid.Clear(); lipo.Clear(); don.Clear(); acc.Clear(); }

  OEBitVec &GetBitVec(void)         { return grid; }
  OEBitVec &GetLipoBitVec(void)     { return lipo; }
  OEBitVec &GetDonorBitVec(void)    { return don; }
  OEBitVec &GetAcceptorBitVec(void) { return acc; }
};

}

#endif /* __BITGRID_H__ */
