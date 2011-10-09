/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.04
 * September 9, 2011
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST. Modifications and additions by IUPAC 
 * and the InChI Trust.
 *
 * IUPAC/InChI-Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0.
 * Copyright (C) IUPAC and InChI Trust Limited
 * 
 * This library is free software; you can redistribute it and/or modify it under the 
 * terms of the IUPAC/InChI Trust Licence for the International Chemical Identifier 
 * (InChI) Software version 1.0; either version 1.0 of the License, or 
 * (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the IUPAC/InChI Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0 for more details.
 * 
 * You should have received a copy of the IUPAC/InChI Trust Licence for the 
 * International Chemical Identifier (InChI) Software version 1.0 along with 
 * this library; if not, write to:
 * 
 * The InChI Trust
 * c/o FIZ CHEMIE Berlin
 * Franklinstrasse 11
 * 10587 Berlin
 * GERMANY
 * 
 */


#ifndef __INCHIRING_H__
#define __INCHIRING_H__
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


#endif /* __INCHIRING_H__ */
