/**********************************************************************
matrix.cpp - Operations on arbitrary-sized matrix.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison

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
#include <openbabel/babelconfig.h>
#include <openbabel/matrix.h>
#include <cstdio>

using namespace std;

namespace OpenBabel
{

void print_matrix(std::vector<std::vector<double> > &m)
{
    unsigned int i,j;

    for (i = 0; i < m.size(); ++i)
    {
        for (j = 0; j < m[i].size(); ++j)
            printf("%5.2f",m[i][j]);
        printf("\n");
    }
}

void print_matrix_f(double *m, int rows, int cols)
{
    int i,j,idx;

    for (i = 0; i < rows; ++i)
    {
        idx = i * cols;
        for (j = 0; j < cols; ++j)
            printf("%5.2f",m[idx+j]);
        printf("\n");
    }
}

void print_matrix_ff(double **m, int rows, int cols)
{
    int i,j;

    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
            printf("%5.2f",m[i][j]);
        printf("\n");
    }
}

bool mult_matrix(std::vector<std::vector<double> > &c,
		 std::vector<std::vector<double> > &a,
		 std::vector<std::vector<double> > &b)
{
    unsigned int i,j,k;

    if (a.size() != b.size())
        return(false);

    c.resize(a.size());

    for (i = 0; i < a.size(); ++i)
    {
        c[i].resize(b[i].size());
        for (j = 0; j < b[i].size(); ++j)
        {
            c[i][j] = 0.0;
            for (k = 0; k < a[i].size(); ++k)
                c[i][j] = c[i][j] + a[i][k] * b[k][j];
        }
    }

    return(true);
}

bool mult_matrix_f(double *c, double *a, double *b, int rows, int cols)
{
    int i,j,k,idx;

    for ( i = 0 ; i < rows ; i++ )
    {
        idx = i * cols;
        for ( j = 0; j < cols ; j++ )
        {
            c[idx+j] = 0.0;
            for ( k = 0; k < cols ; k++ )
                c[idx+j] = c[idx+j] + a[idx+k] * b[(k*cols)+j];
        }
    }

    return(true);
}

bool mult_matrix_ff(double **c, double **a, double **b, int rows, int cols)
{
    int i,j,k;

    for ( i = 0 ; i < rows ; i++ )
        for ( j = 0; j < cols ; j++ )
        {
            c[i][j] = 0.0;
            for ( k = 0; k < cols ; k++ )
                c[i][j] = c[i][j] + a[i][k] * b[k][j];
        }

    return(true);
}

bool invert_matrix(std::vector<std::vector<double> > &mat, double &det)
{
    int  i, j, k, m, n, row = 0, col = 0;
    double tempo, big, pvt;

    vector<int> pvt_ind;
    vector<vector<int> > index;

    int cols = mat[0].size();
    int rows = mat.size();

    pvt_ind.resize(mat[0].size());

    index.resize(mat.size());
    for (i = 0; (unsigned)i < mat.size(); ++i)
        index[i].resize(2);

    // make sure we have a square matrix
    // #rows == #cols;
    if (cols != rows)
    {
        det = 0.0;
        return(false);
    }

    det = 1.0;

    for (i = 0; i < cols; ++i)
        pvt_ind[i] = rows+1;

    for (i = 0; i < cols; ++i)
    {
        big = 0.0;
        for (j = 0; j < cols; ++j)
        {
            if (pvt_ind[j] != 0)
                for (k = 0; k < cols; ++k)
                {
                    if (fabs(big) < fabs(mat[j][k]))
                    {
                        row = j;
                        col = k;
                        big = mat[j][k];
                    }
                }
        }

        pvt_ind[col]++;
        if (row != col)
        {
            det = -det;
            for (m = 0; m < cols; ++m)
            {
                tempo = mat[row][m];
                mat[row][m] = mat[col][m];
                mat[col][m] = tempo;
            }
        }

        index[i][0] = row;
        index[i][1] = col;
        pvt = mat[col][col];
        det *= pvt;

        mat[col][col] = 1.0;

        for (m = 0; m < cols; ++m)
            mat[col][m] /= pvt;

        for (n = 0; n < cols; ++n)
            if (n != col)
            {
                tempo = mat[n][col];
                mat[n][col] = 0.0;
                for (m = 0; m < cols; ++m)
                    mat[n][m] -= mat[col][m] * tempo;
            }
    }

    for (i = 0; i < cols; ++i)
    {
        m = cols - 1;
        if (index[m][0] != index[m][1])
        {
            row = index[m][0];
            col = index[m][1];
            for (k = 0; k < cols; ++k)
            {
                tempo = mat[k][row];
                mat[k][row] = mat[k][col];
                mat[k][col] = tempo;
            }
        }
    }

    return(true);
}

bool invert_matrix_f(double *mat, double &det, int rows, int cols)
{
    int  i, j, k, m, n, row = 0, col = 0, idx1, idx2;
    double tempo, big, pvt;

    vector<int> pvt_ind;
    vector<vector<int> > index;

    pvt_ind.resize(cols);
    index.resize(rows);

    for (i = 0; i < rows; ++i)
        index[i].resize(2);

    // make sure we have a square matrix
    // #rows == #cols;
    if (cols != rows)
    {
        det = 0.0;
        return(false);
    }

    det = 1.0;

    for (i = 0; i < cols; ++i)
        pvt_ind[i] = rows+1;

    for (i = 0; i < cols; ++i)
    {
        big = 0.0;
        for (j = 0; j < cols; ++j)
        {
            if (pvt_ind[j] != 0)
            {
                idx1 = (j * cols);
                for (k = 0; k < cols; ++k)
                {
                    idx2 = idx1 + k;
                    if (fabs(big) < fabs(mat[idx2]))
                    {
                        row = j;
                        col = k;
                        big = mat[idx2];
                    }
                }
            }
        }

        pvt_ind[col]++;
        if (row != col)
        {
            det  = -det;
            idx1 = row * cols;
            idx2 = col * cols;
            for (m = 0; m < cols; ++m)
            {
                tempo = mat[idx1+m];
                mat[idx1+m] = mat[idx2+m];
                mat[idx2+m] = tempo;
            }
        }

        index[i][0] = row;
        index[i][1] = col;

        idx1 = (col*cols);
        pvt  = mat[idx1+col];
        det *= pvt;

        mat[idx1+col] = 1.0;

        for (m = 0; m < cols; ++m)
            mat[idx1+m] /= pvt;

        for (n = 0; n < cols; ++n)
            if (n != col)
            {
                idx1  = n * cols;
                tempo = mat[idx1 + col];
                mat[idx1 + col] = 0.0;

                idx2 = col * cols;
                for (m = 0; m < cols; ++m)
                    mat[idx1 + m] -= mat[idx2 + m] * tempo;
            }
    }

    for (i = 0; i < cols; ++i)
    {
        m = cols - 1;
        if (index[m][0] != index[m][1])
        {
            row = index[m][0];
            col = index[m][1];
            for (k = 0; k < cols; ++k)
            {
                idx1  = (k * cols);
                idx2  = idx1 + col;
                idx1 += row;

                tempo = mat[idx1];
                mat[idx1] = mat[idx2];
                mat[idx2] = tempo;
            }
        }
    }

    return(true);
}

bool invert_matrix_ff(double **mat, double &det, int rows, int cols)
{
    int  i, j, k, m, n, row = 0, col = 0;
    double tempo, big, pvt;

    vector<int> pvt_ind;
    vector<vector<int> > index;

    pvt_ind.resize(cols);
    index.resize(rows);

    for (i = 0; i < rows; ++i)
        index[i].resize(2);

    // make sure we have a square matrix
    // #rows == #cols;
    if (cols != rows)
    {
        det = 0.0;
        return(false);
    }

    det = 1.0;

    for (i = 0; i < cols; ++i)
        pvt_ind[i] = rows+1;

    for (i = 0; i < cols; ++i)
    {
        big = 0.0;
        for (j = 0; j < cols; ++j)
        {
            if (pvt_ind[j] != 0)
                for (k = 0; k < cols; ++k)
                {
                    if (fabs(big) < fabs(mat[j][k]))
                    {
                        row = j;
                        col = k;
                        big = mat[j][k];
                    }
                }
        }

        pvt_ind[col]++;
        if (row != col)
        {
            det = -det;
            for (m = 0; m < cols; ++m)
            {
                tempo = mat[row][m];
                mat[row][m] = mat[col][m];
                mat[col][m] = tempo;
            }
        }

        index[i][0] = row;
        index[i][1] = col;
        pvt = mat[col][col];
        det *= pvt;

        mat[col][col] = 1.0;

        for (m = 0; m < cols; ++m)
            mat[col][m] /= pvt;

        for (n = 0; n < cols; ++n)
            if (n != col)
            {
                tempo = mat[n][col];
                mat[n][col] = 0.0;
                for (m = 0; m < cols; ++m)
                    mat[n][m] -= mat[col][m] * tempo;
            }
    }

    for (i = 0; i < cols; ++i)
    {
        m = cols - 1;
        if (index[m][0] != index[m][1])
        {
            row = index[m][0];
            col = index[m][1];
            for (k = 0; k < cols; ++k)
            {
                tempo = mat[k][row];
                mat[k][row] = mat[k][col];
                mat[k][col] = tempo;
            }
        }
    }

    return(true);
}

bool convert_matrix_f(std::vector<std::vector<double> > &src, double *dst)
{
  unsigned int i, j, idx = 0;

    for ( i = 0 ; i < src.size() ; i++ )
      for ( j = 0 ; j < src[i].size() ; j++ )
            dst[idx++] = src[i][j];

    return true;
}

bool convert_matrix_ff(std::vector<std::vector<double> > &src, double **dst)
{
    unsigned int i, j;

    for ( i = 0 ; i < src.size() ; i++ )
        for ( j = 0 ; j < src[i].size() ; j++ )
            dst[i][j] = src[i][j];

    return true;
}

bool convert_matrix_f(double *src, std::vector<std::vector<double> > &dst,
		      int rows, int cols)
{
    int i, j, idx;

    dst.resize(rows);
    for ( i = 0 ; i < rows ; i++ )
    {
        idx = i * cols;
        dst[i].resize(cols);
        for ( j = 0 ; j < cols ; j++ )
            dst[i][j] = src[idx+j];
    }

    return true;
}

bool convert_matrix_ff(double **src, std::vector<std::vector<double> > &dst,
		       int rows, int cols)
{
    int i, j;

    dst.resize(rows);
    for ( i = 0 ; i < rows ; i++ )
    {
        dst[i].resize(cols);
        for ( j = 0 ; j < cols ; j++ )
            dst[i][j] = src[i][j];
    }

    return true;
}

bool convert_matrix_f_ff(double *src, double **dst, int rows, int cols)
{
    int i, j, idx;

    for ( i = 0 ; i < rows ; i++ )
    {
        idx = i * cols;
        for ( j = 0 ; j < cols ; j++ )
            dst[i][j] = src[idx+j];
    }

    return true;
}

bool convert_matrix_ff_f(double **src, double *dst, int rows, int cols)
{
    int i, j, idx;

    for ( i = 0 ; i < rows ; i++ )
    {
        idx = i * cols;
        for ( j = 0 ; j < cols ; j++ )
            dst[idx+j] = src[i][j];
    }

    return true;
}

} // end namespace OpenBabel

//! \file matrix.cpp
//! \brief Operations on arbitrary-sized matrix.

