/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02-beta
 * August 23, 2007
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


#ifndef __INCHIRING_H__
#define __INCHIRING_H__
#define QUEUE_QINT 1
typedef AT_RANK qInt; /* queue optimization: known type */

#if( QUEUE_QINT == 1 )
#define QINT_TYPE qInt
#else
#define QINT_TYPE void
#endif

typedef struct tagQieue {
    QINT_TYPE *Val;
    int nTotLength;
    int nFirst;  /* element to remove if nLength > 0 */
    int nLength; /* (nFirst + nLength) is next free position */
#if( QUEUE_QINT != 1 )
    int nSize;
#endif
}QUEUE;

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

QUEUE *QueueCreate( int nTotLength, int nSize );
QUEUE *QueueDelete( QUEUE *q );
int is_bond_in_Nmax_memb_ring( inp_ATOM* atom, int at_no, int neigh_ord, QUEUE *q, AT_RANK *nAtomLevel, S_CHAR *cSource, AT_RANK nMaxRingSize );
int is_atom_in_3memb_ring( inp_ATOM* atom, int at_no );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __INCHIRING_H__ */
