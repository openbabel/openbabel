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

#ifndef _CHIRAL_H
#define _CHIRAL_H

#include "matrix.h"

namespace OpenBabel {

void GraphPotentials(OBMol &mol, std::vector<float> &pot);
void construct_g_matrix(OBMol &mol, std::vector<std::vector<float> > &m);
void construct_c_matrix(OBMol &mol, std::vector<std::vector<float > > &m);
float CalcSignedVolume(OBMol &mol,OBAtom*);
float signed_volume(const Vector &a, const Vector &b, const Vector &c, const Vector &d);
void GetChirality(OBMol &mol, std::vector<int> &chirality);

}

#endif
