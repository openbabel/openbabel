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

void GraphPotentials(OEMol &mol, vector<float> &pot);
void construct_g_matrix(OEMol &mol, vector<vector<float> > &m);
void construct_c_matrix(OEMol &mol,vector<vector<float > > &m);
float CalcSignedVolume(OEMol &mol,OEAtom*);
float signed_volume(const Vector &a, const Vector &b, const Vector &c, const Vector &d);
void GetChirality(OEMol &mol, vector<int> &chirality);
#endif
