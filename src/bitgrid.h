/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
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

#ifndef OB_BITGRID_H
#define OB_BITGRID_H

#include "mol.h"
#include "grid.h"
#include "bitvec.h"
#include "obutil.h"
#include "parsmart.h"
#include "patty.h"

namespace OpenBabel {

// Define rint for win32 - this is not defined in Math.h for VC6.0
class BitGrid 
{
protected:

  bool  fuzzy;

  double xmin, ymin, zmin;
  double xmax, ymax, zmax;
  double midx, midy, midz;
  int   xdim, ydim, zdim, xydim;

  double spacing, inv_spa;
  int   size;

  OBBitVec grid,lipo,don,acc;

  patty p;
  std::vector<std::string> types;

public:
  
  BitGrid(void);
  BitGrid(bool);
  ~BitGrid(void);

  void Init(OBMol &, double);
  void Init(double,double,double,double,double,double,double);
  void Build(OBMol &);
  void Build(OBMol &, OBBitVec &);
  void Build(OBMol &, std::vector<int> &);
  void SetBits(OBAtom *);

  void Clear(void) { grid.Clear(); lipo.Clear(); don.Clear(); acc.Clear(); }

  OBBitVec &GetBitVec(void)         { return grid; }
  OBBitVec &GetLipoBitVec(void)     { return lipo; }
  OBBitVec &GetDonorBitVec(void)    { return don; }
  OBBitVec &GetAcceptorBitVec(void) { return acc; }
};

}

#endif // OB_BITGRID_H
