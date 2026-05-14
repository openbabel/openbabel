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


#ifndef _INCHIRING_H_
#define _INCHIRING_H_

#include "extr_ct.h"

#define QUEUE_QINT 1
typedef AT_RANK qInt; /* queue optimization: known type */

#if ( QUEUE_QINT == 1 )
#define QINT_TYPE qInt
#else
#define QINT_TYPE void
#endif

typedef struct tagQieue {
    QINT_TYPE *Val;
    int nTotLength;
    int nFirst;  /* element to remove if nLength > 0 */
    int nLength; /* (nFirst + nLength) is next free position */
#if ( QUEUE_QINT != 1 )
    int nSize;
#endif
}QUEUE;

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

    QUEUE *QueueCreate( int nTotLength, int nSize );
    QUEUE *QueueDelete( QUEUE *q );
    int is_bond_in_Nmax_memb_ring( inp_ATOM* atom, int at_no, int neigh_ord, QUEUE *q, AT_RANK *nAtomLevel, S_CHAR *cSource, AT_RANK nMaxRingSize );
    int is_atom_in_3memb_ring( inp_ATOM* atom, int at_no );

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif    /* _INCHIRING_H_ */
