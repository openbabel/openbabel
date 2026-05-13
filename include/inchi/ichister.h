/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.07
 * April 30, 2024
 *
 * MIT License
 *
 * Copyright (c) 2024 IUPAC and InChI Trust
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*
* The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST.
 * Modifications and additions by IUPAC and the InChI Trust.
 * Some portions of code were developed/changed by external contributors
 * (either contractor or volunteer) which are listed in the file
 * 'External-contributors' included in this distribution.
 *
 * info@inchi-trust.org
 *
*/


#ifndef _ICHISTER_H_
#define _ICHISTER_H_

#include "ichicomn.h"

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif
    int bCanAtomBeAStereoCenter( char *elname, S_CHAR charge, S_CHAR radical );
    int bCanInpAtomBeAStereoCenter( inp_ATOM *at, int cur_at, int bPointedEdgeStereo, int bStereoAtZz );
    int bCanAtomHaveAStereoBond( char *elname, S_CHAR charge, S_CHAR radical );
    int bCanAtomBeTerminalAllene( char *elname, S_CHAR charge, S_CHAR radical );
    int bCanAtomBeMiddleAllene( char *elname, S_CHAR charge, S_CHAR radical );
    int bAtomHasValence3( char *elname, S_CHAR charge, S_CHAR radical );

    double dot_prod3(const double a[], const double b[]); 
    void* cross_prod3(const double a[], const double b[], double result[]);

    struct tagCANON_GLOBALS;
    int set_stereo_parity( struct tagCANON_GLOBALS *pCG,
                           inp_ATOM* at,
                           sp_ATOM* at_output,
                           int num_at,
                           int num_removed_H,
                           int *nMaxNumStereoAtoms,
                           int *nMaxNumStereoBonds,
                           INCHI_MODE nMode,
                           int bPointedEdgeStereo,
                           int vABParityUnknown,
                           int bLooseTSACheck,
                           int bStereoAtZz );

    int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt,
                              int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );

#define PES_BIT_POINT_EDGE_STEREO    1
#define PES_BIT_PHOSPHINE_STEREO     2
#define PES_BIT_ARSINE_STEREO        4
#define PES_BIT_FIX_SP3_BUG          8

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif    /* _ICHISTER_H_ */
