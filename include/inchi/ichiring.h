/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.01
 * July 21, 2006
 * Developed at NIST
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
