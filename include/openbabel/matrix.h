/**********************************************************************
matrix.h - Operations on arbitrary-sized matrix.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_MATRIX_H
#define OB_MATRIX_H
#include <openbabel/babelconfig.h>
#include <vector>
#include <math.h>

namespace OpenBabel
{

/*
 * vector<vector<double> > m : m[row][col] gives appropriate row/col entry
 * double **m                : m[row][col] gives appropriate row/col entry
 * double *m                 : m[(row * numCols) + col] gives appropriate row/col entry
 */

void print_matrix(std::vector<std::vector<double> > &m);
void print_matrix_f (double  *m, int rows, int cols);
void print_matrix_ff(double **m, int rows, int cols);

bool mult_matrix(std::vector<std::vector<double> > &c, std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &b);
bool mult_matrix_f (double  *c, double  *a, double  *b, int rows, int cols);
bool mult_matrix_ff(double **c, double **a, double **b, int rows, int cols);

bool invert_matrix(std::vector<std::vector<double> > &m, double &det);
bool invert_matrix_f (double  *m, double &det, int rows, int cols);
bool invert_matrix_ff(double **m, double &det, int rows, int cols);

bool convert_matrix_f (std::vector<std::vector<double> > &src, double  *dst);
bool convert_matrix_ff(std::vector<std::vector<double> > &src, double **dst);
bool convert_matrix_f (double  *src, std::vector<std::vector<double> > &dst, int rows, int cols);
bool convert_matrix_ff(double **src, std::vector<std::vector<double> > &dst, int rows, int cols);
bool convert_matrix_ff_f(double **src, double  *dst, int rows, int cols);
bool convert_matrix_f_ff(double  *src, double **dst, int rows, int cols);

}

#endif // OB_MATRIX_H

//! \file matrix.h
//! \brief Operations on arbitrary-sized matrix.
