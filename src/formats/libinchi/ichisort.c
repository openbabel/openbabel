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

#include <string.h>

#include "mode.h"
#include "ichicomn.h"
#include "ichicant.h"

#include "bcf_s.h"

#if 0
#define RET_MAX 32767
#define CUTOFF 8            /* testing shows that this is good value */
#endif

/* Note: the theoretical number of stack entries required is
no more than 1 + log2(num).  But we switch to insertion
sort for CUTOFF elements or less, so we really only need
1 + log2(num) - log2(CUTOFF) stack entries.  For a CUTOFF
of 8, that means we need no more than 30 stack entries for
32 bit platforms, and 62 for 64-bit platforms. */
#define STKSIZ (8*sizeof(void*) - 2)


/****************************************************************************
 inchi's qsort
****************************************************************************/
void inchi_qsort( void *pParam,
                  void *base,
                  size_t num,
                  size_t width,
                  int( *comp )( const void *, const void *, void * ) )
{
    char *lo, *hi;              /* ends of sub-array currently sorting */
    char *mid;                  /* points to middle of subarray */
    char *loguy, *higuy;        /* traveling pointers for partition step */
    size_t size;                /* size of the sub-array */
    char *lostk[STKSIZ], *histk[STKSIZ];
    int stkptr;                 /* stack for saving sub-array to be processed */

    if (num < 2)
        return;                 /* nothing to do */

    stkptr = 0;                 /* initialize stack */

    lo = (char *) base;
    hi = (char *) base + width * ( num - 1 );        /* initialize limits */

    /* This entry point is for pseudo-recursion calling: setting
    lo and hi and jumping to here is like recursion, but stkptr is
    preserved, locals aren't, so we preserve stuff on the stack */
recurse:

    size = ( hi - lo ) / width + 1;        /* number of el's to sort */

    /* First we pick a partitioning element.  The efficiency of the
    algorithm demands that we find one that is approximately the median
    of the values, but also that we select one fast.  We choose the
    median of the first, middle, and last elements, to avoid bad
    performance in the face of already sorted data, or data that is made
    up of multiple sorted runs appended together.  Testing shows that a
    median-of-three algorithm provides better performance than simply
    picking the middle element for the latter case. */

    mid = lo + ( size / 2 ) * width;      /* find middle element */

    /* Sort the first, middle, last elements into order */
    if (comp( lo, mid, pParam ) > 0)
    {
        inchi_swap( lo, mid, width );
    }
    if (comp( lo, hi, pParam ) > 0)
    {
        inchi_swap( lo, hi, width );
    }
    if (comp( mid, hi, pParam ) > 0)
    {
        inchi_swap( mid, hi, width );
    }

    /* We now wish to partition the array into three pieces, one consisting
    of elements <= partition element, one of elements equal to the
    partition element, and one of elements > than it.  This is done
    below; comments indicate conditions established at every step. */

    loguy = lo;
    higuy = hi;

    /* Note that higuy decreases and loguy increases on every iteration,
    so loop must terminate. */
    for (;;)
    {
        /* lo <= loguy < hi, lo < higuy <= hi,
        A[i] <= A[mid] for lo <= i <= loguy,
        A[i] > A[mid] for higuy <= i < hi,
        A[hi] >= A[mid] */

        /* The doubled loop is to avoid calling comp(mid,mid), since some
        existing comparison funcs don't work when passed the same
        value for both pointers. */

        if (mid > loguy)
        {
            do
            {
                loguy += width;
            }
            while (loguy < mid && comp( loguy, mid, pParam ) <= 0);
        }
        if (mid <= loguy)
        {
            do
            {
                loguy += width;
            }
            while (loguy <= hi && comp( loguy, mid, pParam ) <= 0);
        }

        /* lo < loguy <= hi+1, A[i] <= A[mid] for lo <= i < loguy,
        either loguy > hi or A[loguy] > A[mid] */

        do
        {
            higuy -= width;
        }
        while (higuy > mid && comp( higuy, mid, pParam ) > 0);

         /* lo <= higuy < hi, A[i] > A[mid] for higuy < i < hi,
         either higuy == lo or A[higuy] <= A[mid] */

        if (higuy < loguy)
        {
            break;
        }

        /* if loguy > hi or higuy == lo, then we would have exited, so
        A[loguy] > A[mid], A[higuy] <= A[mid],
        loguy <= hi, higuy > lo */

        inchi_swap( loguy, higuy, width );

        /* If the partition element was moved, follow it.  Only need
        to check for mid == higuy, since before the swap,
        A[loguy] > A[mid] implies loguy != mid. */

        if (mid == higuy)
        {
            mid = loguy;
        }

        /* A[loguy] <= A[mid], A[higuy] > A[mid]; so condition at top
        of loop is re-established */
    }

    /*     A[i] <= A[mid] for lo <= i < loguy,
    A[i] > A[mid] for higuy < i < hi,
    A[hi] >= A[mid]
    higuy < loguy
    implying:
    higuy == loguy-1
    or higuy == hi - 1, loguy == hi + 1, A[hi] == A[mid] */

    /* Find adjacent elements equal to the partition element.  The
    doubled loop is to avoid calling comp(mid,mid), since some
    existing comparison funcs don't work when passed the same value
    for both pointers. */

    higuy += width;
    if (mid < higuy)
    {
        do
        {
            higuy -= width;
        }
        while (higuy > mid && comp( higuy, mid, pParam ) == 0);
    }
    if (mid >= higuy)
    {
        do
        {
            higuy -= width;
        }
        while (higuy > lo && comp( higuy, mid, pParam ) == 0);
    }

    /* OK, now we have the following:
    higuy < loguy
    lo <= higuy <= hi
    A[i]  <= A[mid] for lo <= i <= higuy
    A[i]  == A[mid] for higuy < i < loguy
    A[i]  >  A[mid] for loguy <= i < hi
    A[hi] >= A[mid] */

    /* We've finished the partition, now we want to sort the subarrays
    [lo, higuy] and [loguy, hi].
    We do the smaller one first to minimize stack usage.
    We only sort arrays of length 2 or more.*/

    if (higuy - lo >= hi - loguy)
    {
        if (lo < higuy)
        {
            lostk[stkptr] = lo;
            histk[stkptr] = higuy;
            ++stkptr;
        }                           /* save big recursion for later */

        if (loguy < hi)
        {
            lo = loguy;
            goto recurse;           /* do small recursion */
        }
    }
    else
    {
        if (loguy < hi)
        {
            lostk[stkptr] = loguy;
            histk[stkptr] = hi;
            ++stkptr;               /* save big recursion for later */
        }

        if (lo < higuy)
        {
            hi = higuy;
            goto recurse;           /* do small recursion */
        }
    }

    /* We have sorted the array, except for any pending sorts on the stack.
    Check if there are any, and do them. */

    --stkptr;
    if (stkptr >= 0)
    {
        lo = lostk[stkptr];
        hi = histk[stkptr];
        goto recurse;           /* pop subarray from stack */
    }
    else
    {
        return;                 /* all subarrays done */
    }
}


/****************************************************************************/
void inchi_swap( char *a, char *b, size_t width )
{
    char tmp;
    if (a != b)
    {
        while (width--)
        {
            tmp = *a;
            *a++ = *b;
            *b++ = tmp;
        }
    }
}


/****************************************************************************
 Sort by insertions
****************************************************************************/
int insertions_sort( void *pCG,
                     void *base,
                     size_t num, size_t width,
                     int( *compare )( const void *, const void *, void * ) ) /* djb-rwth: types of variables are sufficient */
{
    char *i, *j, *pk = (char*) base;
    int  num_trans = 0;
    size_t k;
    for (k = 1; k < num; k++, pk += width)
    {
        /*for( i = pk, j = pk + width; j > (char*)base && (*compare)(i,j) > 0; j=i, i -= width )*/
        for (i = j = pk + width;
             j > ( char* )base && ( i -= width, ( *compare )( i, j, pCG ) ) > 0;
             j = i)        /* changed to keep BoundsChecker happy 2007-09-24 DT */
        {
            inchi_swap( i, j, width );
            num_trans++;
        }
    }

    return num_trans;
}


/****************************************************************************
 Sort by insertions
****************************************************************************/
int insertions_sort_AT_NUMBERS( void *pCG,
                                AT_NUMB *base,
                                int num,
                                int( *compare )( const void *e1, const void *e2, void * ) )
{
    AT_NUMB *i, *j, *pk, tmp;
    int  k, num_trans = 0;
    for (k = 1, pk = base; k < num; k++, pk++)
    {
        for (j = ( i = pk ) + 1, tmp = *j; j > base && ( *compare )( i, &tmp, pCG ) > 0; j = i, i--)
        {
            *j = *i;
            num_trans++;
        }
        *j = tmp;
    }

    return num_trans;
}


/****************************************************************************
 Sort neighbors according to ranks in ascending order
****************************************************************************/
void insertions_sort_NeighList_AT_NUMBERS( NEIGH_LIST base, AT_RANK *nRank )
{
    AT_NUMB *i, *j, *pk, tmp;
    AT_RANK rj; /* optimization */
    int k, num = (int) *base++;
    for (k = 1, pk = base; k < num; k++, pk++)
    {
        for (j = ( i = pk ) + 1, rj = nRank[(int) *j]; j > base && nRank[(int) *i] > rj; j = i, i--)
        {
            tmp = *i;
            *i = *j;
            *j = tmp;
        }
    }
}


/****************************************************************************
 Sort neighbors according to ranks in ascending order
****************************************************************************/
int insertions_sort_AT_RANK( AT_RANK *base, int num )
{
    AT_RANK *i, *j, *pk, tmp;
    int  k, num_trans = 0;
    for (k = 1, pk = base; k < num; k++, pk++)
    {
        for (j = ( i = pk ) + 1, tmp = *j; j > base && *i > tmp; j = i, i--)
        {
            *j = *i;
            num_trans++;
        }
        *j = tmp;
    }

    return num_trans;
}


/****************************************************************************
 Sort neighbors according to ranks in ascending order
****************************************************************************/
int insertions_sort_NeighList_AT_NUMBERS3( NEIGH_LIST base, AT_RANK *nRank )
{
    AT_NUMB *i, *j, *pk, tmp;
    AT_RANK rj;
    int k, n, num = (int) *base++;
    for (k = 1, pk = base, n = 0; k < num; k++, pk++)
    {
        for (j = ( i = pk ) + 1, rj = nRank[(int) ( tmp = *j )];
             j > base && nRank[(int) *i] > rj;
             j = i, i--)
        {
            *j = *i;
            n++;
        }
        *j = tmp;
    }

    return n;
}


/****************************************************************************
 Sort neighbors according to symm. ranks (primary key) and canon.
 ranks (secondary key), in descending order
****************************************************************************/
void insertions_sort_NeighListBySymmAndCanonRank( NEIGH_LIST base,
                                                  const AT_RANK *nSymmRank,
                                                  const AT_RANK *nCanonRank )
{
    AT_NUMB *i, *j, *pk, tmp;
    int  diff;
    int k, num = (int) *base++;
    for (k = 1, pk = base; k < num; k++, pk++)
    {
        for (j = ( i = pk ) + 1;
             j > base &&    /*  always j > i */
             ( 0 > ( diff = (int) nSymmRank[(int) *i] - (int) nSymmRank[(int) *j] ) ||
             (!diff && nCanonRank[(int) *i] < nCanonRank[(int) *j]) ); /* djb-rwth: addressing LLVM warning */
             j = i, i--)
        {
            tmp = *i;
            *i = *j;
            *j = tmp;
        }
    }

}


/****************************************************************************
 *
 *  Comparison functions
 *
 ****************************************************************************/


/****************************************************************************/
int CompNeighborsAT_NUMBER( const void* a1, const void* a2, void *p )
{
    CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
#ifdef CT_NEIGH_INCREASE
    return (int) pCG->m_pn_RankForSort[pCG->m_pNeighborsForSort[( int )*(const AT_NUMB*) a1]] -
           (int) pCG->m_pn_RankForSort[pCG->m_pNeighborsForSort[( int )*(const AT_NUMB*) a2]];
#else
    return (int) ( (CANON_GLOBALS *) pCG )->m_pn_RankForSort[pNeighborsForSort[( int )*(const AT_NUMB*) a2]] -
           (int) ( (CANON_GLOBALS *) pCG )->m_pn_RankForSort[pNeighborsForSort[( int )*(const AT_NUMB*) a1]];
#endif
}


/**********************************************************************************/
int comp_AT_RANK( const void* a1, const void* a2, void *p )
{
    return ( int )*(const AT_RANK*) a1 - ( int )*(const AT_RANK*) a2;
}


/**********************************************************************************/
/*  Compare for sorting Ranks only */
int CompRank( const void* a1, const void* a2, void *p )
{
    CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    int ret = (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a1] -
              (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a2];

    return ret;
}


/**********************************************************************************/
int CompRanksOrd( const void* a1, const void* a2, void *p )
{
    int ret;
    CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    ret = (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a1] -
          (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a2];
    if (!ret)
    {
        ret = ( int )*(const AT_RANK*) a1 - ( int )*(const AT_RANK*) a2;
    }

    return ret;
}


/**********************************************************************************/
int CompAtomInvariants2Only( const void* a1, const void* a2, void *p )
{
    CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    const ATOM_INVARIANT2 *pAI1 = pCG->m_pAtomInvariant2ForSort + ( int )*(const AT_RANK*) a1;
    const ATOM_INVARIANT2 *pAI2 = pCG->m_pAtomInvariant2ForSort + ( int )*(const AT_RANK*) a2;
    int i;
    for (i = 0; i < AT_INV_BREAK1; i++)
    {
        if (pAI1->val[i] == pAI2->val[i])
            continue;
        return  (int) pAI1->val[i] - (int) pAI2->val[i];
    }
    if (pAI1->iso_sort_key != pAI2->iso_sort_key)
    {
        return ( pAI1->iso_sort_key > pAI2->iso_sort_key ) ? 1 : -1;
    }
    for (; i < AT_INV_LENGTH; i++)
    {
        if (pAI1->val[i] != pAI2->val[i])
        {
            continue;
        }
        return  (int) pAI1->val[i] - (int) pAI2->val[i];
    }
    if (pAI1->iso_aux_key != pAI2->iso_aux_key)
    {
        return ( pAI1->iso_aux_key > pAI2->iso_aux_key ) ? 1 : -1;
    }

    return 0;
}


/**********************************************************************************/
int CompAtomInvariants2( const void* a1, const void* a2, void *p )
{
    /*  Warning: the following line may be compiler implementation dependent */
    int ret = CompAtomInvariants2Only( a1, a2, p );
    if (!ret)
    {
        ret = ( int )*(const AT_RANK*) a1 - ( int )*(const AT_RANK*) a2;
    }

    return ret;
}


/**********************************************************************************/
/*  Compare two elements lexicographically */
int CompChemElemLex( const void *a1, const void *a2 )
{
    return memcmp( a1, a2, 2 );
}


/****************************************************************************
 Lexicographic compare
****************************************************************************/
int CompareNeighListLex( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank )
{
    int len1 = (int) *pp1++;
    int len2 = (int) *pp2++;
    int len = inchi_min( len1, len2 );
    int diff = 0;
    int ret;
    /* djb-rwth: fixing oss-fuzz issue #25642 */
    while ((len > 0) && !diff) 
    {
        len--;
        diff = (int)nRank[*pp1++] - (int)nRank[*pp2++];
    };

    ret = diff ? diff : (len1 - len2);
    return ret;
}


/****************************************************************************
 Lexicographic compare
****************************************************************************/
int CompareNeighListLexUpToMaxRank( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank, AT_RANK nMaxAtNeighRank )
{
    int len1 = (int) *pp1++;
    int len2 = (int) *pp2++;
    int diff = 0;
    int len;
    while (0 < len1 && nRank[pp1[len1 - 1]] > nMaxAtNeighRank)
    {
        len1--;
    }
    while (0 < len2 && nRank[pp2[len2 - 1]] > nMaxAtNeighRank)
    {
        len2--;
    }
    len = inchi_min( len1, len2 );
    while (len-- > 0 && !( diff = (int) nRank[*pp1++] - (int) nRank[*pp2++] ))
    {
        ;
    }

    return diff ? diff : ( len1 - len2 );
}


/****************************************************************************/
int compare_NeighLists( const NEIGH_LIST *op1, const NEIGH_LIST *op2, void *p )
{
    CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    return CompareNeighListLex( *op1, *op2, pCG->m_pn_RankForSort );
}


/****************************************************************************/
int CompNeighListRanks( const void* a1, const void* a2, void *p )
{
    CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    int ret;
    ret = (int) pCG->m_pn_RankForSort[*( (const AT_RANK*) a1 )] -
          (int) pCG->m_pn_RankForSort[*( (const AT_RANK*) a2 )];
    if (!ret)
    {
        ret = compare_NeighLists( pCG->m_pNeighList_RankForSort + *( (const AT_RANK*) a1 ),
                                     pCG->m_pNeighList_RankForSort + *( (const AT_RANK*) a2 ), p );
    }

    return ret;
}


/****************************************************************************/
int CompNeighLists( const void* a1, const void* a2, void *p )
{
    CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    int ret;
    ret = compare_NeighLists( pCG->m_pNeighList_RankForSort + *( (const AT_RANK*) a1 ),
                              pCG->m_pNeighList_RankForSort + *( (const AT_RANK*) a2 ), p );

    return ret;
}


/****************************************************************************/
int CompNeighListsUpToMaxRank( const void* a1, const void* a2, void *p )
{
    CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
    int ret;
    ret = CompareNeighListLexUpToMaxRank( pCG->m_pNeighList_RankForSort[*( (const AT_RANK*) a1 )],
                                          pCG->m_pNeighList_RankForSort[*( (const AT_RANK*) a2 )],
                                          pCG->m_pn_RankForSort, pCG->m_nMaxAtNeighRankForSort );

    return ret;
}


/****************************************************************************/
int CompNeighListRanksOrd( const void* a1, const void* a2, void *p )
{
    int ret = CompNeighListRanks( a1, a2, p );
    if (!ret)
    {
        ret = ( int )*( (const AT_RANK*) a1 ) - ( int )*( (const AT_RANK*) a2 ); /*  keep original order if identical */
    }

    return ret;
}


/****************************************************************************/
int CompRanksInvOrd( const void* a1, const void* a2, void *p )
{
    return ( int )*(const AT_RANK*) a2 - ( int )*(const AT_RANK*) a1;
}


/****************************************************************************/
int CompNeighborsRanksCountEql( const void* a1, const void* a2, void *p )
{
    CANON_GLOBALS *pCG = (CANON_GLOBALS *) p;
#ifdef CT_NEIGH_INCREASE
    int ret = (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a1] -
              (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a2];
#else
    int ret = (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a2] -
              (int) pCG->m_pn_RankForSort[( int )*(const AT_RANK*) a1];
#endif
    pCG->m_nNumCompNeighborsRanksCountEql += !ret;

    return ret;
}


/****************************************************************************
 *
 * In this neighbor list the (vertex number) = (canonical number) - 1
 * Since LinearCT is sorted so that parents are in ascending order
 * and all neighbors of a parent are smaller than the parent and are
 * in ascending order, the neighbors in the NEIGH_LIST are automatically
 * sorted in ascending order
****************************************************************************/
NEIGH_LIST *CreateNeighListFromLinearCT( AT_NUMB *LinearCT, int nLenCT, int num_atoms )
{
    /* Atom numbers in LinearCT are canonical numbers
     * order: parent[i] > neigh[i][0] < neigh[i][1]...<neigh[i][n] < parent[i+1] > neigh[i+1][0] < ...
     *        parent[i] < parent[i+1]
     */
    int i, j;
    S_CHAR     *valence = NULL;
    NEIGH_LIST *pp = NULL;
    AT_NUMB    *pAtList = NULL;
    AT_RANK     n_vertex, n_neigh;
    int err = 1, num_bonds;
    int length, start;

    if ((int) LinearCT[0] > num_atoms)
    {
        goto exit_function;
    }
    valence = (S_CHAR*)inchi_calloc((long long)num_atoms + 1, sizeof(valence[0]));
    if (!valence) /* djb-rwth: cast operator added */
    {
        goto exit_function;
    }

    for (i = 1, num_bonds = 0, n_vertex = LinearCT[0]; i < nLenCT; i++)
    {
        if (( n_neigh = LinearCT[i] ) < n_vertex)
        {
            valence[n_neigh] ++;
            valence[n_vertex] ++;
            num_bonds += 2;
        }
        else
        {
            if ((int) ( n_vertex = n_neigh ) > num_atoms)
            {
                goto exit_function;
            }
        }
    }
    if ((int) n_vertex != num_atoms)
    {
        goto exit_function;
    }
    length = num_bonds + num_atoms + 1;
    pp = (NEIGH_LIST*)inchi_calloc(((long long)num_atoms + 1), sizeof(NEIGH_LIST));
    pAtList = (AT_NUMB*)inchi_malloc(length * sizeof(AT_NUMB));
    if (pp && pAtList) /* djb-rwth: cast operator added; addressing LLVM warning */
    {
        /*  Create empty connection table */
        for (i = 1, length = 0; i <= num_atoms; i++)
        {
            start = length;
            length += ( valence[i] + 1 );
            pp[i - 1] = pAtList + start;
            pp[i - 1][0] = 0;
        }
        /*  Fill out the CT */
        for (i = 1, n_vertex = LinearCT[0] - 1; i < nLenCT; i++)
        {
            if (( n_neigh = LinearCT[i] - 1 ) < n_vertex)
            {
                /*  Vertex - neighbor connection */
                j = (int) ( ++pp[(int) n_vertex][0] );
                pp[(int) n_vertex][j] = n_neigh;
                /*  neighbor - vertex connection */
                j = (int) ( ++pp[(int) n_neigh][0] );
                pp[(int) n_neigh][j] = n_vertex;
            }
            else
            {
                if ((int) ( n_vertex = n_neigh ) >= num_atoms)
                {
                    goto exit_function;
                }
            }
        }
        err = 0;
    }

exit_function:
    if (valence)
    {
        inchi_free( valence );
    }
    if (err) /* djb-rwth: ignoring LLVM warning */
    {
        if (pAtList)
        {
            inchi_free( pAtList );
        }
        if (pp)
        {
            inchi_free( pp );
            pp = NULL;
        }
    }

    return pp; /* djb-rwth: ignoring LLVM warning: since a pointer is returned, memory should be freed in a function which calls *CreateNeighListFromLinearCT */
}


/****************************************************************************
 NEIGH_LIST pp[] is an array of pointers to the lists of neighboring
 atoms numbers. The first number in each list is a number of neighbors.
 In case of bDoubleBondSquare != 0 neighbors connected by the double bond
 appear 2 times. The first element pp[0] is a pointer to be deallocated
 to free all the lists.
****************************************************************************/
NEIGH_LIST *CreateNeighList( int num_atoms,
                             int num_at_tg,
                             sp_ATOM* at,
                             int bDoubleBondSquare,
                             T_GROUP_INFO *t_group_info )
{
    /*  +1 to add NULL termination */
    NEIGH_LIST *pp = (NEIGH_LIST *) inchi_calloc( ( (long long)num_at_tg + 1 ), sizeof( NEIGH_LIST ) ); /* djb-rwth: cast operator added */
    T_GROUP   *t_group = NULL;
    AT_NUMB   *nEndpointAtomNumber = NULL;
    int        num_t_groups = 0;
    int        nFirstEndpointAtNoPos;

    AT_NUMB   *pAtList = NULL;
    int        length, start, val, i, j;
    if (pp)
    {
        if (num_at_tg > num_atoms)
        {
            t_group = t_group_info->t_group;
            num_t_groups = t_group_info->num_t_groups;
            nEndpointAtomNumber = t_group_info->nEndpointAtomNumber;
        }

        if (!bDoubleBondSquare)
        {
            for (i = 0, length = 0; i < num_atoms; i++)
            {
                length += (int) at[i].valence + ( num_t_groups && at[i].endpoint );
            }
            length += num_atoms;
            for (i = 0; i < num_t_groups; i++)
            {
                length += (int) t_group[i].nNumEndpoints;
            }
            length += num_t_groups;
        }
        else
        {
            for (i = 0, length = 0; i < num_atoms; i++)
            {
                val = (int) at[i].valence;
                for (j = 0; j < val; j++)
                {
                    length += 1 + ( bDoubleBondSquare && BOND_DOUBLE == at[i].bond_type[j] );
                }
                length += ( num_t_groups && at[i].endpoint );
            }
            length += num_atoms;
            for (i = 0; i < num_t_groups; i++)
            {
                length += (int) t_group[i].nNumEndpoints;
            }
            length += num_t_groups;
        }
        length++; /*  +1 to save number of neighbors */
        pAtList = (AT_NUMB*)inchi_malloc(length * sizeof(AT_NUMB));
        if (pAtList) /* djb-rwth: addressing LLVM warning */
        {
            if (!bDoubleBondSquare)
            {
                for (i = 0, length = 0; i < num_atoms; i++)
                {
                    val = at[i].valence;
                    start = length++;
                    for (j = 0; j < val; j++)
                    {
                        pAtList[length++] = at[i].neighbor[j];
                    }
                    /*  add endpoint */
                    if (num_t_groups && at[i].endpoint)
                    {
                        pAtList[length++] = num_atoms + (int) at[i].endpoint - 1;
                    }
                    pAtList[start] = length - start - 1;  /*  number of neighbors before the list of neighbors */
                    pp[i] = pAtList + start;              /*  pointer to the <num.neigh.><list of neigh> */
                }
            }
            else
            {
                for (i = 0, length = 0; i < num_atoms; i++)
                {
                    val = at[i].valence;
                    start = length++;
                    for (j = 0; j < val; j++)
                    {
                        pAtList[length++] = at[i].neighbor[j]; /* djb-rwth: buffer overrun avoided implicitly */
                        if (bDoubleBondSquare && BOND_DOUBLE == at[i].bond_type[j])
                        {
                            pAtList[length++] = at[i].neighbor[j]; /*  a list of neighbor orig. numbers */
                        }
                    }
                    /*  Add endpoint */
                    if (num_t_groups && at[i].endpoint)
                    {
                        pAtList[length++] = num_atoms + (int) at[i].endpoint - 1;
                    }
                    pAtList[start] = length - start - 1;  /*  number of neighbors before the list of neighbors */
                    pp[i] = pAtList + start;              /*  pointer to the <num.neigh.><list of neigh> */
                }
            }

            /*  Add t-groups */
            for (i = 0; i < num_t_groups; i++)
            {
                val = (int) t_group[i].nNumEndpoints;
                start = length++;
                nFirstEndpointAtNoPos = (int) t_group[i].nFirstEndpointAtNoPos;
                for (j = 0; j < val; j++)
                {
                    pAtList[length++] = nEndpointAtomNumber[nFirstEndpointAtNoPos + j];
                }
                pAtList[start] = length - start - 1;  /*  number of neighbors before the list of neighbors */
                pp[num_atoms + i] = pAtList + start;    /*  pointer to the <num.neigh.><list of neigh> */
            }
        }
        else
        {
            inchi_free(pAtList); /* djb-rwth: fixing coverity ID #499598 */
            inchi_free( pp );
            return NULL;
        }
    } /* djb-rwth: ignoring LLVM warning */

    /* djb-rwth: fixing coverity ID #499598 -- pp uses pAtList values */

    return pp; /* djb-rwth: ignoring LLVM warning: since a pointer is returned, memory should be freed in a function which calls *CreateNeighList */
}


/****************************************************************************/
void FreeNeighList( NEIGH_LIST *pp )
{
    if (pp)
    {
        if (pp[0])
        {
            inchi_free( pp[0] );
        }
        inchi_free( pp );
    }

}


/****************************************************************************/
int BreakAllTies( CANON_GLOBALS *pCG,
                  int num_atoms,
                  int num_max,
                  AT_RANK **pRankStack,
                  NEIGH_LIST *NeighList,
                  AT_RANK *nTempRank,
                  CANON_STAT *pCS )
{
    int i, nRet = -1, nNumRanks = 1 /* value does not matter*/;

    AT_RANK *nPrevRank = *pRankStack++;
    AT_RANK *nPrevAtomNumber = *pRankStack++;

    AT_RANK *nNewRank = NULL;
    AT_RANK *nNewAtomNumber = NULL;

    if (!pRankStack[0])
    {
        pRankStack[0] = (AT_RANK *) inchi_malloc( num_max * sizeof( *nNewRank ) );
    }
    if (!pRankStack[1])
    {
        pRankStack[1] = (AT_RANK *) inchi_malloc( num_max * sizeof( *nNewAtomNumber ) );
    }
    if (!pRankStack[0] || !pRankStack[1])
    {
        return CT_OUT_OF_RAM;  /*   <BRKPT> */
    }
    nNewRank = pRankStack[0];
    nNewAtomNumber = pRankStack[1];

    if (nNewRank && nNewAtomNumber)
    {
        memcpy(nNewAtomNumber, nPrevAtomNumber, num_atoms * sizeof(nNewAtomNumber[0]));
        memcpy(nNewRank, nPrevRank, num_atoms * sizeof(nNewRank[0]));
        for (i = 1, nRet = 0; i < num_atoms; i++)
        {
            /*  12-12-2001: replaced Prev... with New... */
            if (nNewRank[(int) nNewAtomNumber[i - 1]] == nNewRank[(int) nNewAtomNumber[i]])
            {
                nNewRank[nNewAtomNumber[i - 1]] = (AT_RANK) i;
                nNumRanks = DifferentiateRanks2( pCG, num_atoms, NeighList,
                                                 nNumRanks, nNewRank, nTempRank,
                                                 nNewAtomNumber,
                                                 &pCS->lNumNeighListIter, 1 );
                pCS->lNumBreakTies++;
                nRet++;
            }
        }
    }

    return nRet;
}


/****************************************************************************
 Int insertions sort
****************************************************************************/
int * iisort( int *list, int num )
{
    int i;
    for (i = 1; i < num; i++)
    {
        int tmp = list[i];
        int j = i - 1;
        while (j >= 0 && list[j] > tmp)
        {
            list[j + 1] = list[j];
            j--;
        }
        list[j + 1] = tmp;
    }

    return list;
}
