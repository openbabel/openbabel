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

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <math.h>

using namespace std;

namespace OpenEye {

class Vector;

/*
 * vector<vector<float> > m : m[row][col] gives appropriate row/col entry
 * float **m                : m[row][col] gives appropriate row/col entry
 * float *m                 : m[(row * numCols) + col] gives appropriate row/col entry
 */

void print_matrix(vector<vector<float> > &m);
void print_matrix_f (float  *m, int rows, int cols);
void print_matrix_ff(float **m, int rows, int cols);
 
bool mult_matrix(vector<vector<float> > &c, vector<vector<float> > &a, vector<vector<float> > &b);
bool mult_matrix_f (float  *c, float  *a, float  *b, int rows, int cols);
bool mult_matrix_ff(float **c, float **a, float **b, int rows, int cols);

bool invert_matrix(vector<vector<float> > &m, float &det);
bool invert_matrix_f (float  *m, float &det, int rows, int cols);
bool invert_matrix_ff(float **m, float &det, int rows, int cols);

bool convert_matrix_f (vector<vector<float> > &src, float  *dst);
bool convert_matrix_ff(vector<vector<float> > &src, float **dst);
bool convert_matrix_f (float  *src, vector<vector<float> > &dst, int rows, int cols);
bool convert_matrix_ff(float **src, vector<vector<float> > &dst, int rows, int cols);
bool convert_matrix_ff_f(float **src, float  *dst, int rows, int cols);
bool convert_matrix_f_ff(float  *src, float **dst, int rows, int cols);

}

#endif //MATRIX_H


