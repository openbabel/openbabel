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

#include "mode.h"
#include "ichiring.h"

#include "bcf_s.h"

/* Local prototypes */
int GetMinRingSize( inp_ATOM* atom, QUEUE *q, AT_RANK *nAtomLevel, S_CHAR *cSource, AT_RANK nMaxRingSize );

/*  add to the queue */
int QueueAdd( QUEUE *q, QINT_TYPE *Val );
/*  read & remove from the queue */
int QueueGet( QUEUE *q, QINT_TYPE *Val );
/*  read from the queue */
int QueueGetAny( QUEUE *q, QINT_TYPE *, int ord );
/*  initialize the queue */
int QueueReinit( QUEUE *q );
/*  current queue length */
int QueueLength( QUEUE *q );
/*  number of used queue internal elements */
int QueueWrittenLength( QUEUE *q );


#if ( QUEUE_QINT == 1 )  /* { */


/****************************************************************************/
QUEUE *QueueCreate( int nTotLength, int nSize )
{
    QUEUE     *q = NULL;
    QINT_TYPE *Val = NULL;
    if (nTotLength < 1 || nSize != ( int )sizeof( QINT_TYPE ) ||
         !( q = (QUEUE     *) inchi_calloc( 1, sizeof( QUEUE ) ) ) ||
         !( Val = (QINT_TYPE *) inchi_calloc( nTotLength, nSize ) ))
    {
        if (q) inchi_free( q );
        return NULL;
    }
    q->Val = Val;
    /* q->nSize      = nSize; */
    q->nTotLength = nTotLength;

    return q;
}


/****************************************************************************/
int QueueAdd( QUEUE *q, QINT_TYPE *Val )
{
    if (q && Val && q->nLength < q->nTotLength)
    {
        q->Val[( q->nFirst + q->nLength ) % q->nTotLength] = *Val;
        q->nLength++;
        return q->nLength;
    }

    return -1;
}


/****************************************************************************/
int QueueGet( QUEUE *q, QINT_TYPE *Val )
{
    if (q && Val && q->nLength > 0)
    {
        *Val = q->Val[q->nFirst];
        /* new: do not allow to overwrite the retrieved value */
        q->nFirst = ( q->nFirst == q->nTotLength - 1 ) ? 0 : q->nFirst + 1;
        q->nLength--;
        /* -- old --
        if ( -- q->nLength ) {
            q->nFirst = (q->nFirst == q->nTotLength - 1)? 0 : q->nFirst + 1;
        }
        */
        return q->nLength;
    }

    return -1;
}


/****************************************************************************/
int QueueGetAny( QUEUE *q, QINT_TYPE *Val, int ord )
{
    if (0 <= ord && ord < q->nTotLength)
    {
        *Val = q->Val[ord];
        return  1; /*  success */
    }
    else
    {
        return -1; /*  error */
    }
}

#else  /* } QUEUE_QINT == 1 { */


/****************************************************************************/
QUEUE *QueueCreate( int nTotLength, int nSize )
{
    QUEUE     *q = NULL;
    QINT_TYPE *Val = NULL;
    if (nTotLength < 1 || nSize < 1 ||
         !( q = (QUEUE     *) inchi_calloc( 1, sizeof( QUEUE ) ) ) ||
         !( Val = (QINT_TYPE *) inchi_calloc( nTotLength, nSize ) ))
    {
        if (q) inchi_free( q );
        return NULL;
    }
    q->Val = Val;
    q->nSize = nSize;
    q->nTotLength = nTotLength;
    return q;
}
int QueueAdd( QUEUE *q, QINT_TYPE *Val )
{
    if (q && Val && q->nLength < q->nTotLength)
    {
        memcpy( (char*) q->Val + ( ( q->nFirst + q->nLength ) % q->nTotLength )*q->nSize, Val, q->nSize );
        q->nLength++;
        return q->nLength;
    }

    return -1;
}


/****************************************************************************/
int QueueGet( QUEUE *q, QINT_TYPE *Val )
{
    if (q && Val && q->nLength > 0)
    {
        memcpy( Val, (char*) q->Val + q->nFirst * q->nSize, q->nSize );
        if (--q->nLength)
        {
            q->nFirst = ( q->nFirst == q->nTotLength - 1 ) ? 0 : q->nFirst + 1;
        }
        return q->nLength;
    }

    return -1;
}


/****************************************************************************/
int QueueGetAny( QUEUE *q, QINT_TYPE *Val, int ord )
{
    if (0 <= ord && ord < q->nTotLength)
    {
        memcpy( Val, (char*) q->Val + ord * q->nSize, q->nSize );
        return  1; /*  success */
    }
    else
    {
        return -1; /*  error */
    }
}

#endif /* } QUEUE_QINT == 1 */


/****************************************************************************/
QUEUE *QueueDelete( QUEUE *q )
{
    if (q)
    {
        if (q->Val) inchi_free( q->Val );
        inchi_free( q );
    }

    return NULL;
}


/****************************************************************************/
int QueueReinit( QUEUE *q )
{
    if (q)
    {
        q->nFirst = 0;
        q->nLength = 0;
        /* memset( q->Val, 0, q->nTotLength*sizeof(q->Val[0])); */ /*  for debug only */
        return q->nTotLength;
    }

    return -1;
}


/****************************************************************************/
int QueueLength( QUEUE *q )
{
    if (q)
    {
        return q->nLength;
    }
    else
    {
        return 0;
    }
}


/****************************************************************************/
int QueueWrittenLength( QUEUE *q )
{
    if (q)
    {
        int len = q->nFirst + q->nLength;
        return ( len > q->nTotLength ) ? q->nTotLength : len;
    }
    else
    {
        return 0;
    }
}


/****************************************************************************
  BFS: Breadth First Search
****************************************************************************/
int GetMinRingSize( inp_ATOM* atom,
                    QUEUE *q,
                    AT_RANK *nAtomLevel,
                    S_CHAR *cSource,
                    AT_RANK nMaxRingSize )
{
    int qLen, i, j;
    AT_RANK nCurLevel, nRingSize, nMinRingSize = MAX_ATOMS + 1;
    qInt at_no, next;
    int  iat_no, inext;

    while ((qLen = QueueLength( q ))) /* djb-rwth: addressing LLVM warning */
    {
        /*  traverse the next level (next outer ring) */
        for (i = 0; i < qLen; i++)
        {
            if (0 <= QueueGet( q, &at_no ))
            {
                iat_no = (int) at_no;
                nCurLevel = nAtomLevel[iat_no] + 1;
                if (2 * nCurLevel > nMaxRingSize + 4)
                {
                    /*  2*nCurLevel = nRingSize + 3 + k, k = 0 or 1  */
                    if (nMinRingSize < MAX_ATOMS + 1)
                    {
                        return ( nMinRingSize >= nMaxRingSize ) ? 0 : nMinRingSize;
                    }
                    return 0; /*  min. ring size > nMaxRingSize */
                }
                for (j = 0; j < atom[iat_no].valence; j++)
                {
                    next = (qInt) atom[iat_no].neighbor[j];
                    inext = (int) next;
                    if (!nAtomLevel[inext])
                    {
                        /*  the at_no neighbor has not been traversed yet. Add it to the queue */
                        if (0 <= QueueAdd( q, &next ))
                        {
                            nAtomLevel[inext] = nCurLevel;
                            cSource[inext] = cSource[iat_no]; /*  keep the path number */
                        }
                        else
                        {
                            return -1; /*  error */
                        }
                    }
                    else
                    {
                        if (nAtomLevel[inext] + 1 >= nCurLevel &&
                                cSource[inext] != cSource[iat_no]
                                /*  && cSource[(int)next] != -1 */
                              )
                        {
                            /*  found a ring closure */
                            /*  debug */
                            if (cSource[inext] == -1)
                            {
                                return -1;  /*  error */
                            }
                            if (( nRingSize = nAtomLevel[inext] + nCurLevel - 2 ) < nMinRingSize)
                            {
                                nMinRingSize = nRingSize;
                            }
                            /* return (nRingSize >= nMaxRingSize)? 0 : nRingSize; */
                        }
                    }
                }
            }
            else
            {
                return -1; /*  error */
            }
        }
    }

    if (nMinRingSize < MAX_ATOMS + 1)
    {
        return ( nMinRingSize >= nMaxRingSize ) ? 0 : nMinRingSize;
    }

    return 0;
}


/****************************************************************************
   Return value:
     0:    nMaxRingSize < 3 or
           min. ring size >= nMaxRingSize or
           not a ring bond (the last is currently impossible: bond is known to belong to a ring system.
     n>0:  min. ring size < nMaxRingSize
     n<0:  error

  Input:
     atom[]
     at_no      number of the 1st atom adjacent to the bond
     neigh_ord  ordering number of the bond in question: at[at_no].bond_type[neigh_ord]
     q          queue structure
     nAtomLevel work array, DFS distance
     cSource    work array, origin mark
****************************************************************************/
int is_bond_in_Nmax_memb_ring( inp_ATOM* atom,
                               int at_no,
                               int neigh_ord,
                               QUEUE *q,
                               AT_RANK *nAtomLevel,
                               S_CHAR *cSource,
                               AT_RANK nMaxRingSize )
{
    int  nMinRingSize = -1, i;
    qInt n;
    int  nTotLen;

    if (nMaxRingSize < 3)
    {
        return 0;
    }

    QueueReinit( q );

    /*  mark the starting atom */
    nAtomLevel[at_no] = 1;
    cSource[at_no] = -1;
    /*  add neighbors */
    for (i = 0; i < atom[at_no].valence; i++)
    {
        n = (qInt) atom[at_no].neighbor[i];
        nAtomLevel[(int) n] = 2;
        cSource[(int) n] = 1 + ( i == neigh_ord );
        QueueAdd( q, &n );
    }

    nMinRingSize = GetMinRingSize( atom, q, nAtomLevel, cSource, nMaxRingSize );
    /*  cleanup */
    nTotLen = QueueWrittenLength( q );
    for (i = 0; i < nTotLen; i++)
    {
        if (0 < QueueGetAny( q, &n, i ))
        {
            nAtomLevel[(int) n] = 0;
            cSource[(int) n] = 0;
        }
    }
    nAtomLevel[at_no] = 0;
    cSource[at_no] = 0;

    /*
    if ( nAtomLevel )
        inchi_free ( nAtomLevel );
    if ( cSource )
        inchi_free ( cSource );
    QueueDelete( q );
    */

    return nMinRingSize;
}


/****************************************************************************/
int is_atom_in_3memb_ring( inp_ATOM* atom, int at_no )
{
    AT_NUMB   neigh_neigh;
    int       i, j, k, val, val_neigh, neigh;

    if (atom[at_no].nNumAtInRingSystem < 3)
    {
        return 0;
    }

    for (i = 0, val = atom[at_no].valence; i < val; i++)
    {
        neigh = (int) atom[at_no].neighbor[i];
        if (atom[at_no].nRingSystem != atom[neigh].nRingSystem)
        {
            continue;
        }
        for (j = 0, val_neigh = atom[neigh].valence; j < val_neigh; j++)
        {
            neigh_neigh = atom[neigh].neighbor[j];
            if ((int) neigh_neigh == at_no)
            {
                continue;
            }
            for (k = 0; k < val; k++)
            {
                if (atom[at_no].neighbor[k] == neigh_neigh)
                {
                    return 1;
                }
            }
        }
    }

    return 0;
}
