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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "mode.h"

#include "incomdef.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichicant.h"
#include "ichicomn.h"

#include "ichicomp.h"

#define MAP_MODE_STD  0 /* Standard approach: switch 2 neighbors */
#define MAP_MODE_C2v  1 /* Check for C2v reflection leading to parity inversion */
#define MAP_MODE_C2   2 /* Check for C2 rotation preserving parities */
#define MAP_MODE_S4   3 /* Check for S4 rotation/reflection leading to parity inversion */
/* important: MAP_MODE_STD < (MAP_MODE_C2v, MAP_MODE_C2) < MAP_MODE_S4 */

/* local prototypes */
void DeAllocateForNonStereoRemoval( AT_RANK **nAtomNumberCanon1, AT_RANK **nAtomNumberCanon2,
                                    NEIGH_LIST **nl, NEIGH_LIST **nl1, NEIGH_LIST **nl2, AT_RANK **nVisited1, AT_RANK **nVisited2 );
int AllocateForNonStereoRemoval( sp_ATOM *at, int num_atoms, const AT_RANK *nSymmRank, AT_RANK *nCanonRank,
                            AT_RANK **nAtomNumberCanon1, AT_RANK **nAtomNumberCanon2,
                            NEIGH_LIST **nl, NEIGH_LIST **nl1, NEIGH_LIST **nl2, AT_RANK **nVisited1, AT_RANK **nVisited2 );
AT_RANK GetMinNewRank(AT_RANK *nAtomRank, AT_RANK *nAtomNumb, AT_RANK nRank1 );
int BreakNeighborsTie(  sp_ATOM *at, int num_atoms, int num_at_tg, int ib, int ia,
                        AT_RANK *neigh_num, int in1, int in2, int mode,
                        AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList,
                        const AT_RANK *nSymmRank, AT_RANK *nCanonRank, NEIGH_LIST *nl1, NEIGH_LIST *nl2, long *lNumIter );
int CheckNextSymmNeighborsAndBonds( sp_ATOM *at, AT_RANK cur1, AT_RANK cur2, AT_RANK n1, AT_RANK n2,
                                    AT_RANK *nAvoidCheckAtom, AT_RANK *nVisited1, AT_RANK *nVisited2,
                                    AT_RANK *nVisitOrd1, AT_RANK *nVisitOrd2, const AT_RANK *nRank1, const AT_RANK *nRank2 );
int CreateCheckSymmPaths( sp_ATOM *at, AT_RANK prev1, AT_RANK cur1, AT_RANK prev2, AT_RANK cur2,
                         AT_RANK *nAvoidCheckAtom, AT_RANK *nVisited1, AT_RANK *nVisited2,
                         AT_RANK *nVisitOrd1, AT_RANK *nVisitOrd2,
                         NEIGH_LIST *nl1, NEIGH_LIST *nl2, const AT_RANK *nRank1, const AT_RANK *nRank2,
                         AT_RANK *nCanonRank, AT_RANK *nLength, int *bParitiesInverted, int mode );
int CalculatedPathsParitiesAreIdentical( sp_ATOM *at, int num_atoms, const AT_RANK *nSymmRank,
                         AT_RANK *nCanonRank, AT_RANK *nAtomNumberCanon, AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2,
                         AT_RANK *nVisited1, AT_RANK *nVisited2,
                         AT_RANK prev_sb_neigh, AT_RANK cur, AT_RANK next1, AT_RANK next2, int nNeighMode,
                         int bParitiesInverted, int mode, CANON_STAT *pCS,
                         int vABParityUnknown);
int RemoveCalculatedNonStereoBondParities( sp_ATOM *at, int num_atoms, int num_at_tg,
                                          AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList,
                                          AT_RANK *nCanonRank, const AT_RANK *nSymmRank,
                                          AT_RANK *nAtomNumberCanon, AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2,
                                          NEIGH_LIST *nl, NEIGH_LIST *nl1, NEIGH_LIST *nl2,
                                          AT_RANK *nVisited1, AT_RANK *nVisited2, CANON_STAT *pCS,
                                          int vABParityUnknown);
int RemoveCalculatedNonStereoCenterParities( sp_ATOM *at, int num_atoms, int num_at_tg,
                                          AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList,
                                          AT_RANK *nCanonRank, const AT_RANK *nSymmRank,
                                          AT_RANK *nAtomNumberCanon, AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2,
                                          NEIGH_LIST *nl, NEIGH_LIST *nl1, NEIGH_LIST *nl2, 
                                          AT_RANK *nVisited1, AT_RANK *nVisited2, CANON_STAT *pCS,
                                          int vABParityUnknown);

int SortNeighLists3( int num_atoms, AT_RANK *nRank, NEIGH_LIST *NeighList, AT_RANK *nAtomNumber );





/**************************************************************************************
 *
 *   Convert sorted equivalence information (nSymmRank) to ranks (nRank)
 *   nSymmRank and nRank may point to the same array
 *
 */
int SortedEquInfoToRanks( const AT_RANK* nSymmRank, AT_RANK* nRank, const AT_RANK* nAtomNumber, int num_atoms, int *bChanged )
{
    AT_RANK        rNew, rOld, nNumDiffRanks;
    int            i, j, nNumChanges = 0;
    for ( i = num_atoms-1, j = (int)nAtomNumber[i],
          rOld = nSymmRank[j], rNew = nRank[j] = (AT_RANK)num_atoms,
          nNumDiffRanks = 1;
             i > 0;
                 i -- ) {

        j = (int)nAtomNumber[i-1];
        
        if ( nSymmRank[j] != rOld ) {
            nNumDiffRanks ++;
            rNew = (AT_RANK)i;
            nNumChanges += (rOld != rNew+1);
            rOld = nSymmRank[j];
        }

        nRank[j] = rNew;
    }
    if ( bChanged ) {
        *bChanged = (0 != nNumChanges);
    }
    return nNumDiffRanks;
}
/**************************************************************************************
 *
 *   Convert sorted ranks (nRank) to sorted equivalence information (nSymmRank)
 *   nSymmRank and nRank may point to the same array
 *
 */
int SortedRanksToEquInfo( AT_RANK* nSymmRank, const AT_RANK* nRank, const AT_RANK* nAtomNumber, int num_atoms )
{
    AT_RANK        rNew, rOld, nNumDiffRanks;
    int            i, j;
    for ( i = 1, j = (int)nAtomNumber[0],
          rOld = nRank[j], rNew = nSymmRank[j] = 1,
          nNumDiffRanks = 1;
          i < num_atoms;
          i ++ ) {
        j = (int)nAtomNumber[i];
        if ( nRank[j] != rOld ) {
            nNumDiffRanks ++;
            rNew = (AT_RANK)(i+1);
            rOld = nRank[j];
        }
        nSymmRank[j] = rNew;
    }
    return nNumDiffRanks;
}

/**************************************************************************************/
void switch_ptrs( AT_RANK **p1, AT_RANK **p2 )
{
    AT_RANK *tmp = *p1;
    *p1 = *p2;
    *p2 = tmp;
}

/**************************************************************************************/
/*  Set ranks from the products vector and previous ranks                             */
/*  nRank[] and nNewRank[] should refer to different arrays for now                   */
/**************************************************************************************/
int SetNewRanksFromNeighLists3( int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank,
                                AT_RANK *nNewRank, AT_RANK *nAtomNumber )
{
    int     i, j, nNumDiffRanks, nNumNewRanks;
    AT_RANK r1, r2;
    /*  -- nAtomNumber[] is already properly set --
    for ( i = 0; i < num_atoms; i++ ) {
        nAtomNumber[i] = (AT_RANK)i;
    }
    */
    /*  set globals for qsort */
    pNeighList_RankForSort = NeighList;
    pn_RankForSort         = nRank;
    nNumDiffRanks          = 0;
    nNumNewRanks           = 0;

    memset(nNewRank, 0, num_atoms*sizeof(nNewRank[0]));
    
    /*  sorting */
    for ( i = 0, r1 = 1; i < num_atoms; r1++ ) {
        if ( r1 == (r2 = nRank[j=(int)nAtomNumber[i]]) ) {
            nNewRank[j] = r2;
            nNumDiffRanks ++;
            i ++;
            continue;
        }
        r1 = r2;
        insertions_sort_AT_NUMBERS( nAtomNumber+i, (int)r2-i, CompNeighLists );
        /*insertions_sort( nAtomNumber+i, r2-i, sizeof( nAtomNumber[0] ), CompNeighLists );*/
        j = r2-1;
        nNewRank[(int)nAtomNumber[j]] = r2;
        nNumDiffRanks ++;
        while( j > i ) {
            if ( CompareNeighListLex( NeighList[(int)nAtomNumber[j-1]],
                                      NeighList[(int)nAtomNumber[j]], nRank ) ) {
                r2 = j;
                nNumDiffRanks ++;
                nNumNewRanks ++;
            }
            j --;
            nNewRank[(int)nAtomNumber[j]] = r2;
        }
        i = r1;
    }
    return nNumNewRanks? -nNumDiffRanks : nNumDiffRanks;
}
/**************************************************************************************/
/*  Set ranks from the products vector and previous ranks                             */
/*  When comparing neigh lists ignore ranks > max_at_no                               */
/*  nRank[] and nNewRank[] should refer to different arrays for now                   */
/**************************************************************************************/
int SetNewRanksFromNeighLists4( int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank,
                                AT_RANK *nNewRank, AT_RANK *nAtomNumber, AT_RANK nMaxAtRank )
{
    int     i, j, nNumDiffRanks, nNumNewRanks;
    AT_RANK r1, r2;
    /*  -- nAtomNumber[] is already properly set --
    for ( i = 0; i < num_atoms; i++ ) {
        nAtomNumber[i] = (AT_RANK)i;
    }
    */
    /*  set globals for CompNeighListsUpToMaxRank */
    pNeighList_RankForSort = NeighList;
    pn_RankForSort         = nRank;
    nNumDiffRanks          = 0;
    nNumNewRanks           = 0;
    nMaxAtNeighRankForSort = nMaxAtRank;

    memset(nNewRank, 0, num_atoms*sizeof(nNewRank[0]));
    
    /*  sorting */
    for ( i = 0, r1 = 1; i < num_atoms; r1++ ) {
        if ( r1 == (r2 = nRank[j=(int)nAtomNumber[i]]) ) {
            /* non-tied rank: singleton */
            nNewRank[j] = r2;
            nNumDiffRanks ++;
            i ++;
            continue;
        }
        /* tied rank r2
           r2-i atoms have rank r2
           next atom after them is in position r2
        */
        r1 = r2;
        insertions_sort_AT_NUMBERS( nAtomNumber+i, (int)r2-i, CompNeighListsUpToMaxRank );
        /*insertions_sort( nAtomNumber+i, r2-i, sizeof( nAtomNumber[0] ),  CompNeighListsUpToMaxRank );*/
        j = r2-1; /* prepare cycle backward, from j to i step -1 */
        nNewRank[(int)nAtomNumber[j]] = r2;
        nNumDiffRanks ++;
        while( j > i ) {
            if ( CompareNeighListLexUpToMaxRank( NeighList[nAtomNumber[j-1]],
                                                 NeighList[nAtomNumber[j]], nRank, nMaxAtRank ) ) {
                r2 = j;
                nNumDiffRanks ++;
                nNumNewRanks ++;
            }
            j --;
            nNewRank[(int)nAtomNumber[j]] = r2;
        }
        i = r1;
    }
    return nNumNewRanks? -nNumDiffRanks : nNumDiffRanks;
}
/**************************************************************************************/
/*  Set ranks from the products vector and previous ranks                             */
/*  nRank[] and nNewRank[] should refer to different arrays for now                   */
/**************************************************************************************/
int SetNewRanksFromNeighLists( int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank,
                             AT_RANK *nAtomNumber, int bUseAltSort, int ( *comp )(const void *, const void *) )
{
    int     i,  nNumDiffRanks;
    AT_RANK nCurrentRank;
    /*  -- nAtomNumber[] is already properly set --
    for ( i = 0; i < num_atoms; i++ ) {
        nAtomNumber[i] = (AT_RANK)i;
    }
    */
    /*  set globals for qsort */
    pNeighList_RankForSort = NeighList;
    pn_RankForSort         = nRank;
    
    /*  sorting */
    if ( bUseAltSort & 1 )
        tsort( nAtomNumber, num_atoms, sizeof( nAtomNumber[0] ), comp /*CompNeighListRanksOrd*/ );
    else
        qsort( nAtomNumber, num_atoms, sizeof( nAtomNumber[0] ), comp /*CompNeighListRanksOrd*/ );

    for ( i=num_atoms-1, nCurrentRank=nNewRank[(int)nAtomNumber[i]] = (AT_RANK)num_atoms, nNumDiffRanks = 1;
          0 < i ;
          i -- ) {
        /*  Note: CompNeighListRanks() in following line implicitly reads nRank pointed by pn_RankForSort */
        if ( CompNeighListRanks( &nAtomNumber[i-1], &nAtomNumber[i] ) ) {
            nNumDiffRanks ++;
            nCurrentRank = (AT_RANK)i;
        }
        nNewRank[(int)nAtomNumber[i - 1]] = nCurrentRank;
    }
    
    return nNumDiffRanks;
}
/**************************************************************************************/
/*   Sort NeighList[] lists of neighbors according to the ranks of the neighbors      */
/**************************************************************************************/
void SortNeighListsBySymmAndCanonRank( int num_atoms, NEIGH_LIST *NeighList, const AT_RANK *nSymmRank, const AT_RANK *nCanonRank )
{
    int i;
    for ( i = 0; i < num_atoms; i ++ ) {
        insertions_sort_NeighListBySymmAndCanonRank( NeighList[i], nSymmRank, nCanonRank );
    }
}
/**************************************************************************************/
int SortNeighLists2( int num_atoms, AT_RANK *nRank, NEIGH_LIST *NeighList, AT_RANK *nAtomNumber )
{
    int k, i;
    AT_RANK nPrevRank = 0;
    /*
     * on entry nRank[nAtomNumber[k]] <= nRank[nAtomNumber[k+1]]  ( k < num_atoms-1 )
     *          nRank[nAtomNumber[k]] >= k+1                      ( k < num_atoms )
     *          nRank[nAtomNumber[k]] == k+1 if this nRank value is not tied OR if
     *                nRank[nAtomNumber[k]] < nRank[nAtomNumber[k+1]] OR if k = num_atoms-1.
     *
     */
    for ( k = 0; k < num_atoms; k ++ ) {
        i = nAtomNumber[k];
        if ( (nRank[i] != k+1 || nRank[i] == nPrevRank) && NeighList[i][0] > 1 ) {
            /*  nRank[i] is tied (duplicated) */
            insertions_sort_NeighList_AT_NUMBERS( NeighList[i], nRank );
        }
        nPrevRank = nRank[i];
    }
    return 0;
}
/**************************************************************************************/
int SortNeighLists3( int num_atoms, AT_RANK *nRank, NEIGH_LIST *NeighList, AT_RANK *nAtomNumber )
{
    int k, i;
    AT_RANK nPrevRank = 0;
    /*
     * on entry nRank[nAtomNumber[k]] <= nRank[nAtomNumber[k+1]]  ( k < num_atoms-1 )
     *          nRank[nAtomNumber[k]] >= k+1                      ( k < num_atoms )
     *          nRank[nAtomNumber[k]] == k+1 if this nRank value is not tied OR if
     *                nRank[nAtomNumber[k]] < nRank[nAtomNumber[k+1]] OR if k = num_atoms-1.
     *
     */
    for ( k = 0; k < num_atoms; k ++ ) {
        i = nAtomNumber[k];
        if ( (nRank[i] != k+1 || nRank[i] == nPrevRank) && NeighList[i][0] > 1 ) {
            /*  nRank[i] is tied (duplicated) */
            insertions_sort_NeighList_AT_NUMBERS3( NeighList[i], nRank );
        }
        nPrevRank = nRank[i];
    }
    return 0;
}
/**************************************************************************************
 *
 *  Differentiate2
 *
 * Note: on entry nAtomNumber[] must contain a valid transposition of num_atoms length
 *       for example, nAtomNumber[i] = i;
 * Note2: this version does not calculate neighbor lists for non-tied ranks
 */
int  DifferentiateRanks2( int num_atoms, NEIGH_LIST *NeighList,
                                 int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                                 AT_RANK *nAtomNumber, long *lNumIter, int bUseAltSort )
{
    /*int nNumPrevRanks;*/

    /*  SortNeighLists2 needs sorted ranks */
    pn_RankForSort = pnCurrRank;
    if ( bUseAltSort & 1 )
        tsort( nAtomNumber, num_atoms, sizeof(nAtomNumber[0]), CompRank /* CompRanksOrd*/ );
    else
        qsort( nAtomNumber, num_atoms, sizeof(nAtomNumber[0]), CompRanksOrd );

    do {
        *lNumIter += 1;
        /*nNumPrevRanks = nNumCurrRanks;*/
        switch_ptrs( &pnCurrRank, &pnPrevRank );
        SortNeighLists2( num_atoms, pnPrevRank, NeighList, nAtomNumber );
        /*  the following call creates pnCurrRank out of pnPrevRank */
        nNumCurrRanks = SetNewRanksFromNeighLists( num_atoms, NeighList, pnPrevRank, pnCurrRank, nAtomNumber,
                                                 1, CompNeighListRanksOrd );
    } while ( /*nNumPrevRanks != nNumCurrRanks ||*/ memcmp( pnPrevRank, pnCurrRank, num_atoms*sizeof(pnCurrRank[0]) ) );

    return nNumCurrRanks;
}
/**************************************************************************************
 *
 *  Differentiate3
 *
 * Note: on entry nAtomNumber[] must contain a valid transposition of num_atoms length
 *       for example, nAtomNumber[i] = i;
 * Note2: this version does not calculate neighbor lists for non-tied ranks
 */
int  DifferentiateRanks3( int num_atoms, NEIGH_LIST *NeighList,
                          int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                          AT_RANK *nAtomNumber, long *lNumIter )
{
/*    
    static long count = 0;
    count ++;
    if ( count == 103 ) {
        int stop=1;
    }
*/

    /*  SortNeighLists3 needs sorted ranks: ranks/atnumbers must have been already sorted */
    do {
        *lNumIter += 1;
        switch_ptrs( &pnCurrRank, &pnPrevRank );
        SortNeighLists3( num_atoms, pnPrevRank, NeighList, nAtomNumber );
        /*  the following call creates pnCurrRank out of pnPrevRank */
        nNumCurrRanks = SetNewRanksFromNeighLists3( num_atoms, NeighList, pnPrevRank,
                                                    pnCurrRank, nAtomNumber);
    } while ( nNumCurrRanks < 0 /* memcmp( pnPrevRank, pnCurrRank, num_atoms*sizeof(pnCurrRank[0]) )*/ );

    return nNumCurrRanks;
}
/**************************************************************************************
 *
 *  Differentiate4: ignore neighbors with rank > num_atoms
 *
 * Note: on entry nAtomNumber[] must contain a valid transposition of num_atoms length
 *       for example, nAtomNumber[i] = i;
 * Note2: this version does not sort neighbor lists for non-tied ranks
 */
int  DifferentiateRanks4( int num_atoms, NEIGH_LIST *NeighList,
                          int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                          AT_RANK *nAtomNumber, AT_RANK nMaxAtRank, long *lNumIter )
{
/*    
    static long count = 0;
    count ++;
    if ( count == 103 ) {
        int stop=1;
    }
*/
    /*  SortNeighLists4 needs sorted ranks: ranks/atnumbers must have been already sorted */
    do {
        *lNumIter += 1;
        switch_ptrs( &pnCurrRank, &pnPrevRank );
        SortNeighLists3( num_atoms, pnPrevRank, NeighList, nAtomNumber );
        /*  the following call creates pnCurrRank out of pnPrevRank */
        nNumCurrRanks = SetNewRanksFromNeighLists4( num_atoms, NeighList, pnPrevRank,
                                                    pnCurrRank, nAtomNumber, nMaxAtRank );
    } while ( nNumCurrRanks < 0 /* memcmp( pnPrevRank, pnCurrRank, num_atoms*sizeof(pnCurrRank[0]) )*/ );

    return nNumCurrRanks;
}
/**************************************************************************************
 *
 *  DifferentiateBasic (sort according to ranks only)
 *
 * Note: on entry nAtomNumber[] must contain a valid transposition of num_atoms length
 *       for example, nAtomNumber[i] = i;
 * Note2: this version does not calculate neighbor lists for non-tied ranks
 */
int  DifferentiateRanksBasic( int num_atoms, NEIGH_LIST *NeighList,
                                 int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                                 AT_RANK *nAtomNumber, long *lNumIter, int bUseAltSort )
{
    int nNumPrevRanks;

    /*  SortNeighLists2 needs sorted ranks */
    pn_RankForSort     = pnCurrRank;
    if ( bUseAltSort & 1 )
        tsort( nAtomNumber, num_atoms, sizeof(nAtomNumber[0]), CompRank );
    else
        qsort( nAtomNumber, num_atoms, sizeof(nAtomNumber[0]), CompRank );

    do {
        *lNumIter += 1;
        nNumPrevRanks = nNumCurrRanks;
        switch_ptrs( &pnCurrRank, &pnPrevRank );
        SortNeighLists2( num_atoms, pnPrevRank, NeighList, nAtomNumber );
        /*  the following call creates pnCurrRank out of pnPrevRank */
        nNumCurrRanks = SetNewRanksFromNeighLists( num_atoms, NeighList, pnPrevRank, pnCurrRank, nAtomNumber, bUseAltSort, CompNeighListRanks );
    } while ( nNumPrevRanks != nNumCurrRanks || memcmp( pnPrevRank, pnCurrRank, num_atoms*sizeof(pnCurrRank[0]) ) );
    return nNumCurrRanks;
}

/**************************************************************************************
 * For the purpose of mapping an atom to an atom:
 * (a) find number of tied ranks
 * (b) if number of tied ranks > 1 then:
 *    1) find the rank for breaking a tie
 *    2) allocate memory for breaking the tie if it has not been allocated
 *    3) find out if atom 1 ("from") has already been mapped
 * Return value:
 *  < 0: error
 *  = 1: has already been mapped, to tie to break
 *  > 1: we need to break a tie
 */
int NumberOfTies( AT_RANK **pRankStack1, AT_RANK **pRankStack2, int length,
                  int at_no1, int at_no2, AT_RANK *nNewRank, int *bAddStack, int *bMapped1 )
{

    AT_RANK *nRank1       = *pRankStack1++;
    AT_RANK *nAtomNumber1 = *pRankStack1++;  /*  ranks for mapping "1", "from" */

    AT_RANK *nRank2       = *pRankStack2++;
    AT_RANK *nAtomNumber2 = *pRankStack2++;  /*  ranks for mapping "2", "to" */

    AT_RANK r, *pTempArray;

    int iMax, i, i1, i2;

    *bAddStack = 0;
    *bMapped1  = 0;
    *nNewRank  = 0;
    r = nRank1[at_no1];
    if ( r != nRank2[at_no2] )
        return CT_MAPCOUNT_ERR; /*  atoms cannot be mapped onto each other: they have different ranks */ /*   <BRKPT> */
    iMax = r - 1;
    /*  find i1 and i2 = numbers of ranks in nRank1[] and nRank2[] equal to r:  */
    for ( i1 = 1; i1 <= iMax && r == nRank1[nAtomNumber1[iMax-i1]]; i1 ++ )
        ;
    for ( i2 = 1; i2 <= iMax && r == nRank2[nAtomNumber2[iMax-i2]]; i2 ++ )
        ;
    if ( i2 != i1 )
        return CT_MAPCOUNT_ERR; /*  program error: must be identical number of equal ranks */ /*   <BRKPT> */
    /*  found i1 equal rank(s); preceding (smaller) non-equal rank is r-i1 */
    /*  To break the tie we have to reduce the rank r to r-i1+1 */

    /************ Note *******************************
     * IF ( i=r-1 && 0 <= i && i < num_atoms AND 
     *      nRank[nAtomNumber1[i]] == r )
     * THEN:
     * nRank[nAtomNumber1[i+1]] >  r; (if i+1 < num_atoms)
     * nRank[nAtomNumber1[i-1]] <= r; (if i > 0)
     *
     * IF r = nRank[i] THEN
     * nRank[nAtomNumber1[r-1]] == r
     * nRank[nAtomNumber1[r-i-1]] <= nRank[nAtomNumber1[r-i]] (for 1 <= i < r )
     */
    if ( i1 > 1 ) {
        /* int bAtFromHasAlreadyBeenMapped = 0; */
        *nNewRank = r - i1 + 1;
        /*  grab an existing or allocate a new array */
        /*  we need 4 arrays: 2 for ranks + 2 for numbers */
        for ( i = 0; i < 4; i ++ ) {
            if ( i < 2 ) {
                pTempArray = *pRankStack1;
                *bMapped1 += (pTempArray && pTempArray[0]);
            } else {
                pTempArray = *pRankStack2;
            }
            if ( !pTempArray && !(pTempArray = (AT_RANK *) inchi_malloc(length)))
                return CT_OUT_OF_RAM;  /*  out of RAM */ /*   <BRKPT> */
            /*  copy "to" contents */
            switch( i ) {
            case 2:
                memcpy( pTempArray, nRank2, length );
                break;
            case 3:
                memcpy( pTempArray, nAtomNumber2, length );
                break;
            }
            if ( i < 2 )
                *pRankStack1 ++ = pTempArray;
            else {
                *pRankStack2 ++ = pTempArray;
            }
        }
        *bAddStack = 2; /*  to break the tie we added 2 more arrays to pRankStack1 and pRankStack2 */
    }
    return i1;
}    


/**************************************************************************************
 *
 *
 *
 *               Stereo Mappings
 *
 *
 *
 **************************************************************************************/

/**************************************************************************************
 * Parity for a half of a stereo bond. If both halfs have the same parity
 * then the bond is "trans" (E,-,1), otherwise it is "cis" (Z,+,2).
 * The advantage of this approach is: The bond parity does not depend on the
 * rank of the atom located on the opposite end of the stereogenic bond.
 * As the result all bond parities of, for example, benzene, can be calculated
 * from equivalence ranks only, without any mappings.
 *
 * Input: at_no1     = number of atom for which the half-bond parity is calculated
 *        i_sb_neigh = ordering number of the stereo bond in at->stereo_bond_neighbor[]
 *
 * Returns: 0=> no parity can be found; 1=> odd parity; 2=> even parity
 *
 */
int HalfStereoBondParity( sp_ATOM *at, int at_no1, int i_sb_neigh, const AT_RANK *nRank )
{
/*
   Suppose neighbors #0,#1,#2 have ranks a, b, c. Remove rank of the neighbor connected
   by the stereogenic bond (NCSB) from the a, b, c list and denote the two left as r[0], r[1],
   in the same order. Let iNCSB be an ordering number (0,1,or 2) of the NCSB.
   Assume the neighbor connected by the stereogenic bond has infinite positive rank.
   Position the half-bond so that the stereogenic bond neighbor is to the right from the atom (see below)
   
   Definition.
   ===========
                   if rank(X) != rank(Y) then Half-bond parity = (rank(X) > rank(Y)), that is,
    Y              
     \             if ( rank(X) < rank(Y) ) then Half-bond parity is Even
      C==NCSB      if ( rank(X) > rank(Y) ) then Half-bond parity is Odd
     /             if ( rank(X) = rank(Y) ) then Half-bond parity cannot be defined
    X
    
    1                          2             1         
     \                          \             \        
      C==NCSB       C==NCSB      C==NCSB       C==NCSB 
     /             /            /                      
    2             1            1                       
                                                       
    Parity = 1    Parity = 1   Parity = 2    Parity = 2
    (Odd)         (Odd)       (Even) or 0   (Even) or 0
    
   Half-bond parity =  (iNCSB + (r[0] > r[1]) + (Atom C geometric parity))%2

   Consider the following cases to prove the formula:

   Case 1: 3 explicit neighbors
   ============================
   If  (1) atom's geometric parity = even (which means neighbors #0, #1, #2 are located clockwise),
   and (2) neighbors other than NCSB have different ranks, then,
   assuming that NCSB always has the largest (infinite) rank (this is consistent with
   the assumption that implicit hydrogens have smallest ranks), we have 3 possibilities:

                             c         a          b       
                              \         \          \      
                               C==a      C==b       C==c  
                              /         /          /      
                             b         c          a       
                                                          
            iNCSB      =      0          1          2    
       Half-bond parity =     b>c        a<c        a>b     (0=even, 1=odd)
                           r[0]>r[1]  r[0]<r[1]  r[0]>r[1]
       Half-bond parity
       for all 3 cases      =    (iNCSB + (r[0] > r[1]))%2

       The following slight modification will work for both odd and even geometric parity:

       Half-bond parity     =    (iNCSB + (r[0] > r[1]) + (Atom C geometric parity))%2

       even parity (0) => atom above the bond has lower rank than the atom below the bond.


   Case 2: 2 explicit neighbors
   ============================
   One implicit hydrogen atom H or hydrogen isotope (implicit rank=0). Assume r[1]=0

                             H         a            Note. The same method 
                              \         \                 works for              
                               C==a      C==b                 
                              /         /             N==a   and   a     
                             b         H             /              \    
                                                    b                N==b
            iNCSB       =      0         1     
       Half-bond parity =     b>0       a<0    
       (r[1]=0, r[0]>0)    r[0]>r[1]  r[0]<r[1] 

       Half-bond parity =  (iNCSB + (r[0] > r[1]) + (Atom C geometric parity))%2

   Case 3: 1 explicit neighbor (NCSB)
   ==================================
   Two implicit hydrogens, (number of neighbors on non-streogenic bonds)==0:

   Atom C geometric parity:  Even               Odd          Note. The same method
                                                                   works for                          
                             D                  H                       
                              \                  \           Even   and   Odd               
                               C==a               C==a                                      
                              /                  /           H               N==a           
                             H                  D             \             /               
                                                               N==a        H     
            iNCSB =           0                0
       Half-bond parity =    (0<0)=0         (0<0)+1 = 1
       (r[1]=0, r[0]=0)    r[1]<r[0]         (r[1]<r[0])+atom_parity

       Half-parity
       for this case  =    (iNCSB + (r[0] > r[1]) + (Atom C geometric parity))%2

*/
    int i, j, k, iNeigh, parity, at1_parity, at_no2;
    AT_RANK r[MAX_NUM_STEREO_BOND_NEIGH];

    if ( at[at_no1].valence > MAX_NUM_STEREO_BOND_NEIGH || ( at1_parity = at[at_no1].parity ) <= 0 ) {
        return 0;
    }
    if ( !PARITY_WELL_DEF( at1_parity ) ) {
        if ( PARITY_KNOWN( at1_parity ) ) {
            return at1_parity;
        }
        return -at1_parity;
    }
    if ( 0 > i_sb_neigh || i_sb_neigh >= MAX_NUM_STEREO_BOND_NEIGH ) {
        return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
    }
    for ( i = 0; i <= i_sb_neigh; i ++ ) {
        if ( !at[at_no1].stereo_bond_neighbor[i] ) {
            return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
        }
    }
    at_no2 = at[at_no1].neighbor[(int)at[at_no1].stereo_bond_ord[i_sb_neigh]];
    memset( r, 0, sizeof( r ) ); 
    for ( i = j = 0, iNeigh = -1; i < at[at_no1].valence; i ++ ) {
        if ( (k = (int)at[at_no1].neighbor[i]) == at_no2 ) {
            iNeigh = i;
        } else {
            r[j++] = nRank[k];
        }
    }
    if ( iNeigh < 0 || iNeigh != at[at_no1].stereo_bond_ord[i_sb_neigh] ) {
        return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
    }
    if ( (j > 0 && !r[0]) || (j > 1 && !r[1]) )
        return 0; /*  undefined ranks */

    if ( (j == 2 && r[0] == r[1]) || iNeigh < 0 ) {
        parity = AB_PARITY_CALC;  /*  cannot calculate bond parity without additional breaking ties. */
    } else {
        parity = 2 - (at[at_no1].parity + iNeigh + (r[1] < r[0])) % 2;
    }
    return parity;
}
/**************************************************************************************/
int parity_of_mapped_half_bond( int from_at, int to_at, int from_neigh, int to_neigh,
                           sp_ATOM *at, EQ_NEIGH *pEN,
                           const AT_RANK *nCanonRankFrom, const AT_RANK *nRankFrom, const AT_RANK *nRankTo )
{
    int     i, j, k, num_neigh;
    int     to_sb_neigh_ord, from_sb_neigh_ord, parity;
    AT_RANK r_to[MAX_NUM_STEREO_BOND_NEIGH], at_no_to[MAX_NUM_STEREO_BOND_NEIGH];
    AT_RANK r_canon_from[MAX_NUM_STEREO_BOND_NEIGH], at_no_from[MAX_NUM_STEREO_BOND_NEIGH];
    AT_RANK r, r_sb_neigh;

    for ( i = 0; i < MAX_NUM_STEREO_BOND_NEIGH; i ++ ) {
        r_to[i] = r_canon_from[i] = 0;
    }
    
    if ( pEN ) {
        memset( pEN, 0, sizeof(*pEN));
    }

    /*  for debug only */
    if ( nRankFrom[from_at] != nRankTo[to_at] ||
         nRankFrom[from_neigh] != nRankTo[to_neigh] ||
         at[to_at].valence != at[from_at].valence ) {
        return 0;  /*  program error: both atoms must be mapped */ /*   <BRKPT> */
    }

    parity = PARITY_VAL(at[to_at].parity);
    num_neigh = at[to_at].valence;
    
    if ( num_neigh > MAX_NUM_STEREO_BOND_NEIGH || num_neigh < MIN_NUM_STEREO_BOND_NEIGH ) {
        /*  2 neighbors are possible in case of stereo bond with implicit H */
        /*  or a stereocenter -CHD- with an implicit H */
        if ( num_neigh == 1 && at[to_at].stereo_bond_neighbor[0] ) {
            /*  1 neighbor can happen in case of a terminal =CHD */
            if ( PARITY_WELL_DEF(parity) )
                return 2 - parity % 2;
            else
            if ( parity )
                return parity;
            else
                return AB_PARITY_UNDF; /*  undefined parity */
        }
        return 0;  /*  program error */ /*   <BRKPT> */
    }
    if ( ATOM_PARITY_KNOWN(parity) ) {
        if ( !ATOM_PARITY_WELL_DEF(parity) )
            return parity;
    } else
    if ( parity ) {
        return 0; /* parity; */
    } else {
        return 0; /* AB_PARITY_UNDF; */ /*  possibly program error: undefined parity */
    }
    /*  locate at[to_at].stereo_bond_neighbor[] ordering numbers */
    for ( i = 0, to_sb_neigh_ord=-1; i < MAX_NUM_STEREO_BONDS && (k=(int)at[to_at].stereo_bond_neighbor[i]); i ++ ) {
        if ( k == to_neigh+1 ) {
            to_sb_neigh_ord = i;
            break;
        }
    }
    if ( to_sb_neigh_ord < 0 ) {
        return 0;  /*  program error: not a stereo bond */ /*   <BRKPT> */
    }
    to_sb_neigh_ord = (int)at[to_at].stereo_bond_ord[to_sb_neigh_ord];
    r_sb_neigh   = nRankTo[(int)at[to_at].neighbor[to_sb_neigh_ord]];
    for ( i = j = 0; i < num_neigh; i ++ ) {
        if ( i != to_sb_neigh_ord ) {
            r_to[j] = nRankTo[(int)(at_no_to[j]=at[to_at].neighbor[i])];
            if ( r_sb_neigh == r_to[j] ) {
                return 0; /*  stereo bond atoms are not fully mapped */
            }
            j ++;
        }
    }
    if ( j+1 != num_neigh ) {
        return 0; /*  program error */ /*   <BRKPT> */
    }
    if ( j == 1 ) {
        /*  only one neighbor; no mapping needed */
        return 2-(parity+1+to_sb_neigh_ord)%2;
    }
    if ( j != 2 ) {
        return 0; /*  program error: j can be only 0, 1, or 2 */ /*   <BRKPT> */
    }
    
    if ( r_to[0] == r_to[1] ) {
        /*  double bond neighbors need to be mapped */
        j = 0;
        from_sb_neigh_ord = -1;
        for ( i = 0; i < num_neigh; i ++ ) {
            k = at[from_at].neighbor[i];
            r = nRankFrom[k];
            if ( r == r_sb_neigh ) {
                from_sb_neigh_ord = i;   /*  we need this value only for error-checking */
            } else
            if ( r == r_to[0] ) {
                r_canon_from[j] = nCanonRankFrom[k];
                at_no_from[j]   = (AT_RANK)k;
                j ++;
            } else {
                return 0; /*  program error: unexpected rank, not fully mapped adjacent to the stereo bond atoms */ /*   <BRKPT> */
            }
        }
        if ( from_sb_neigh_ord < 0 || j != 2 ) {
            return 0; /*  program error: rank of a neighbor not found */ /*   <BRKPT> */
        }
        if ( pEN ) { /*  j == 2 */
            pEN->to_at[0] = at_no_to[0];
            pEN->to_at[1] = at_no_to[1];
            pEN->num_to   = 2;           /*  number of stored in pEN->to_at[] central atom neighbors */
            pEN->rank     = r_to[0];     /*  mapping rank of the tied neighbors */
             /*  i := index of the smaller out of r_canon_from[1] and r_canon_from[0] */
            i = (r_canon_from[1] < r_canon_from[0]);
            pEN->from_at    = at_no_from[i];
            pEN->canon_rank = r_canon_from[i];
        }
        return -((int)r_to[0]);
    }
    /*  double bond neighbors a mapped: r_to[0] != r_to[1] */
    from_sb_neigh_ord = -1;
    for ( i = 0; i < num_neigh; i ++ ) {
        k = at[from_at].neighbor[i];
        r = nRankFrom[k];
        if ( r == r_sb_neigh ) {
            from_sb_neigh_ord = i;  /*  we need this value only for error-checking */
        } else
        if ( r == r_to[0] ) {
            r_canon_from[0] = nCanonRankFrom[k];
            /* at_no_from[0]   = (AT_RANK)k; */
        } else
        if ( r == r_to[1] ) {
            r_canon_from[1] = nCanonRankFrom[k];
            /* at_no_from[1]   = (AT_RANK)k; */
        } else {
            return 0; /*  program error: unexpected rank, not fully mapped adjacent to the stereo bond atoms */ /*   <BRKPT> */
        }
    }
    if ( !r_canon_from[0] || !r_canon_from[1] || from_sb_neigh_ord < 0 ) {
        return 0; /*  program error: neighbor rank not found */ /*   <BRKPT> */
    }
    return 2 - (parity + to_sb_neigh_ord + (r_canon_from[1]<r_canon_from[0]))%2;
}

/**************************************************************************************/
int parity_of_mapped_atom2( int from_at, int to_at, const sp_ATOM *at, EQ_NEIGH *pEN,
                           const AT_RANK *nCanonRankFrom, const AT_RANK *nRankFrom, const AT_RANK *nRankTo )
{
    AT_RANK nNeighRankFrom[4], nNeighNumberFrom[4], nNeighRankTo[4], nNeighNumberTo[4];
    AT_RANK nNeighRankFromCanon[4], nNeighRankToCanon[4];
    int     i, j, k, num_neigh;
    int     r1, r2, r, r_canon_from_min, neigh_canon_from_min, r_canon_from;
    int     num_trans_to, num_trans_from, neigh1, neigh2;


    num_neigh = at[to_at].valence;
    
    if ( pEN ) {
        memset( pEN, 0, sizeof(*pEN));
    }

    /*  for debug only */
    if ( nRankFrom[from_at] != nRankTo[to_at] )
        return 0;  /*  program error */ /*   <BRKPT> */
    if ( num_neigh > MAX_NUM_STEREO_ATOM_NEIGH || num_neigh < 2 ) {
        /*  2 neighbors are possible in case of stereo bond with implicit H */
        /*  or a stereocenter >CHD with two implicit H */
        if ( num_neigh == 1 ) {
            /*  1 neighbor can happen in case of a terminal -CHDT or =CHD */
            if ( at[to_at].parity )
                return at[to_at].parity;
            else
                return AB_PARITY_UNDF; /*  undefined parity */
        }
        return 0;  /*  program error */ /*   <BRKPT> */
    }
    for ( i = 0; i < num_neigh; i ++ ) { /*  initialization of locals */
        nNeighNumberTo[i]      =
        nNeighNumberFrom[i]    = i;
        nNeighRankTo[i]        = nRankTo[(int)at[to_at].neighbor[i]];       /* mapping rank */
        nNeighRankFrom[i]      = nRankFrom[j=(int)at[from_at].neighbor[i]]; /* mapping rank */
        nNeighRankFromCanon[i] = nCanonRankFrom[j];                     /* canonical number */
    }

    pn_RankForSort = nNeighRankFrom;
    nNumCompNeighborsRanksCountEql = 0; /*  sort mapping ranks-from */
    num_trans_from = insertions_sort( nNeighNumberFrom, num_neigh, sizeof(nNeighNumberFrom[0]), CompNeighborsRanksCountEql );

    if ( nNumCompNeighborsRanksCountEql ) {
        /*  At least 2 neighbors have equal mapping ranks (are tied). */
        /*  Find tied from-neighbors with minimal canonical rank (nCanonRankFrom[]) */
        r_canon_from_min = MAX_ATOMS+1; /*  max possible rank + 1 */
        for ( i = 1, r = 0, r1 = nNeighRankFrom[neigh1=nNeighNumberFrom[0]]; i < num_neigh; i ++, r1 = r2, neigh1 = neigh2 ) {
            r2 = nNeighRankFrom[neigh2=nNeighNumberFrom[i]];
            if ( r2 == r1 ) {
                /*  found neighbors with tied ranks */
                if ( r != r2 ) {
                    /*  the 1st pair of neighbor with this rank */
                    r = r2;
                    if ( (r_canon_from=nNeighRankFromCanon[neigh1]) < r_canon_from_min ) {
                        r_canon_from_min     = r_canon_from; /*  min canon rank */
                        neigh_canon_from_min = neigh1;       /*  neighbor number */
                    }
                }
                if ( (r_canon_from=nNeighRankFromCanon[neigh2]) < r_canon_from_min ) {
                    r_canon_from_min     = r_canon_from;
                    neigh_canon_from_min = neigh2;
                }
            }
        }
        if ( r ) {
            /*  neighbors with tied ranks have been found => parity cannot be determined without additional mapping */
            /*  find to-neighbors on which neigh_canon_from_min can be mapped */
            r1 = nNeighRankFrom[neigh_canon_from_min];
            if ( pEN ) {
                for ( i = j = 0; i < num_neigh; i ++ ) {
                    if ( r1 == nNeighRankTo[i] ) {
                        pEN->to_at[j++] = at[to_at].neighbor[i];
                    }
                }
                insertions_sort( pEN->to_at, j, sizeof(pEN->to_at[0]), CompRanksInvOrd );
                pEN->num_to     = j;  /*  number of stored in pEN->to_at[] central atom neighbors */
                pEN->from_at    = at[from_at].neighbor[neigh_canon_from_min]; /*  neighbor with min. canon number */
                pEN->rank       = r1; /*  mapping rank of the tied neighbors */
                pEN->canon_rank = r_canon_from_min;  /*  canon. rank of the pEN->from_at */
            } else {
                /*  debug only */
                for ( i = j = 0; i < num_neigh; i ++ ) {
                    if ( r1 == nNeighRankTo[i] ) {
                        j++;
                    }
                }
            }
            /*  debug only */
            if ( j <= 1 || !r1 || r_canon_from_min > MAX_ATOMS ) {
                return 0; /*  program error */ /*   <BRKPT> */
            }
            return -r; /*  means parity cannot be determined */
        }
        return 0; /* program error */
    }
    /*  All neighbors have different mapping ranks; */
    /*  therefore no additional mapping of the neighbors is necessary */
    if ( !ATOM_PARITY_WELL_DEF(at[to_at].parity) )
        return at[to_at].parity; /*  unknown parity or cannot be determined */

    pn_RankForSort = nNeighRankTo;
    num_trans_to   = insertions_sort( nNeighNumberTo, num_neigh, sizeof(nNeighNumberTo[0]), CompNeighborsRanksCountEql );
    
    /*  Map canonical ranks of neighbors. Mapped on each other "to" and "from" atoms have equal mapping ranks */
    for ( i = 0; i < num_neigh; i ++ ) {
        if ( nNeighRankTo[j=nNeighNumberTo[i]] != nNeighRankFrom[k=nNeighNumberFrom[i]] )
            return 0; /*  program error: mapping ranks not equal, from_at neigborhood cannot be mapped on to_at neighbood. */ /*   <BRKPT> */
        nNeighRankToCanon[j] = nNeighRankFromCanon[k]; /*  potential problem: other atom(s) may have same mapping rank and */
                                                       /*  different canon. rank(s). */
        /*  we may save some memory by eliminating nNeighRankFromCanon[]: */
        /*  nNeighRankToCanon[j] = nCanonRankFrom[at[from_at].neighbor[k]] */
    }

    pn_RankForSort  = nNeighRankToCanon;
    num_trans_to   += insertions_sort( nNeighNumberTo, num_neigh, sizeof(nNeighNumberTo[0]), CompNeighborsRanksCountEql );
#ifndef CT_NEIGH_INCREASE
    num_trans_to   += ((num_neigh*(num_neigh-1))/2)%2;  /*  get correct parity for ascending order of canon. numbers */
#endif

    return 2 - (num_trans_to + at[to_at].parity)%2;
}

/**************************************************************************************
 *
 *   Phase II: map canonicaly numbrered structure onto itself
 *             to obtain a minimal or maximal stereo part of the CT
 *
 **************************************************************************************/

int ClearPreviousMappings( AT_RANK **pRankStack1 )
{
    int i;
    for ( i = 0; pRankStack1[i]; i ++ ) {
        pRankStack1[i][0] = 0;
    }
    return i;

}
/**************************************************************************************/
/*  map one atom ("from") onto another ("to"): untie their mapping ranks if they are tied. */
int map_an_atom2( int num_atoms, int num_max, int at_no1/*from*/, int at_no2/*to*/,
                AT_RANK *nTempRank,
                int nNumMappedRanks, int *pnNewNumMappedRanks,
                CANON_STAT *pCS,
                NEIGH_LIST    *NeighList,
                AT_RANK  **pRankStack1, AT_RANK  **pRankStack2, int *bAddStack )
{
    AT_RANK *nRank1,  *nAtomNumber1;  /*  ranks for mapping "1", "from" */
    AT_RANK *nRank2,  *nAtomNumber2;  /*  ranks for mapping "2", "to" */
    AT_RANK *nNewRank1=NULL,  *nNewAtomNumber1=NULL;  /*  ranks for mapping "1", "from" */
    AT_RANK *nNewRank2=NULL,  *nNewAtomNumber2=NULL;  /*  ranks for mapping "2", "to" */
    int     length = num_max*sizeof(AT_RANK);
    int     nNewNumRanks2, nNewNumRanks1;
    int     i, bAtFromHasAlreadyBeenMapped, nNumTies;
    AT_RANK nNewRank;

    nNumTies = NumberOfTies( pRankStack1, pRankStack2, length, at_no1, at_no2, &nNewRank, bAddStack, &bAtFromHasAlreadyBeenMapped );
    
    if ( RETURNED_ERROR(nNumTies) )
        return nNumTies;  /*  error */

    nRank1       = *pRankStack1++;
    nAtomNumber1 = *pRankStack1++;  /*  ranks for mapping "1", "from" */

    nRank2       = *pRankStack2++;
    nAtomNumber2 = *pRankStack2++;  /*  ranks for mapping "2", "to" */
    
    if ( nNumTies > 1 ) {

        nNewRank1       = *pRankStack1++;
        nNewAtomNumber1 = *pRankStack1++;  /*  ranks for mapping "1", "from" */

        nNewRank2       = *pRankStack2++;
        nNewAtomNumber2 = *pRankStack2++;  /*  ranks for mapping "2", "to" */
        /*  break a tie for "to" */
        memcpy( nNewRank2, nRank2, length );
        memcpy( nNewAtomNumber2, nAtomNumber2, length );
        nNewRank2[at_no2] = nNewRank;
        nNewNumRanks2 = DifferentiateRanks2( num_atoms, NeighList,
                                         nNumMappedRanks, nNewRank2, nTempRank,
                                         nNewAtomNumber2, &pCS->lNumNeighListIter, 1 );
        pCS->lNumBreakTies ++;

        /*  Check whether the old mapping can be reused */
        if ( 2 == bAtFromHasAlreadyBeenMapped && nNewRank == nNewRank1[at_no1] ) {
            for ( i = 0; i < num_atoms; i ++ ) {
                if ( nNewRank1[nNewAtomNumber1[i]] != nNewRank2[nNewAtomNumber2[i]] ) {
                    bAtFromHasAlreadyBeenMapped = 0; /*  It cannot. */
                    break;
                }
            }
        } else {
            bAtFromHasAlreadyBeenMapped = 0;
        }
        if ( 2 != bAtFromHasAlreadyBeenMapped ) {
            /*  break a tie for "from" */
            for ( i = 0; pRankStack1[i]; i ++ ) {
                pRankStack1[i][0] = 0;
            }
            memcpy( nNewRank1, nRank1, length );
            memcpy( nNewAtomNumber1, nAtomNumber1, length );  /* GPF: bad nAtomNumber1 */
            nNewRank1[at_no1] = nNewRank;
            nNewNumRanks1 = DifferentiateRanks2( num_atoms, NeighList,
                                             nNumMappedRanks, nNewRank1, nTempRank,
                                             nNewAtomNumber1, &pCS->lNumNeighListIter, 1 );
            pCS->lNumBreakTies ++;
        } else {
            nNewNumRanks1 = nNewNumRanks2;
        }

        if ( nNewNumRanks1 != nNewNumRanks2 )
            return CT_MAPCOUNT_ERR; /*  program error */ /*   <BRKPT> */
        *pnNewNumMappedRanks = nNewNumRanks2;
        /*  debug only */
        for ( i = 0; i < num_atoms; i ++ ) {
            if ( nNewRank1[nNewAtomNumber1[i]] != nNewRank2[nNewAtomNumber2[i]] ) {
                return CT_MAPCOUNT_ERR; /*  program error */ /*   <BRKPT> */
            }
        }
    } else {
        *pnNewNumMappedRanks = nNumMappedRanks;
    }
    return ( nNewRank1 )? nNewRank1[at_no1] : nRank1[at_no1]; /*  mapping rank value */
}

/**************************************************************************************/
int might_change_other_atom_parity( sp_ATOM *at, int num_atoms, int at_no, AT_RANK *nRank2, AT_RANK *nRank1 )
{
    int     i, j, neighbor_no;
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( nRank2[i] != nRank1[i] ) {
            if ( i != at_no /*&& ATOM_PARITY_WELL_DEF(at[i].parity)*/
                && at[i].bHasStereoOrEquToStereo
                && !(at[i].stereo_atom_parity & KNOWN_PARITIES_EQL )
                && !at[i].stereo_bond_neighbor[0]
                ) {

                return 1; /*  may have changed stereo atoms order */
            }
            for ( j = 0; j < at[i].valence; j ++ ) {
                neighbor_no = at[i].neighbor[j];
                if ( neighbor_no != at_no
                     /*&& ATOM_PARITY_WELL_DEF(at[neighbor_no].parity)*/ 
                     && at[neighbor_no].bHasStereoOrEquToStereo
                     && !(at[neighbor_no].stereo_atom_parity & KNOWN_PARITIES_EQL )
                     && !at[neighbor_no].stereo_bond_neighbor[0]
                   )
                    return 1; /*  may have changed stereo atom parity */
            }
        }
    }
    return 0;
}
/**************************************************************************************/
#if ( REMOVE_CALC_NONSTEREO == 1 ) /* { */
/**************************************************************************************/
void DeAllocateForNonStereoRemoval( AT_RANK **nAtomNumberCanon1, AT_RANK **nAtomNumberCanon2,
                                    NEIGH_LIST **nl, NEIGH_LIST **nl1, NEIGH_LIST **nl2, AT_RANK **nVisited1, AT_RANK **nVisited2 )
{
    if ( *nAtomNumberCanon1 ) {
        inchi_free( *nAtomNumberCanon1 );
        *nAtomNumberCanon1 = NULL;
    }
    if ( *nAtomNumberCanon2 ) {
        inchi_free( *nAtomNumberCanon2 );
        *nAtomNumberCanon2 = NULL;
    }
    if ( *nl ) {
        FreeNeighList( *nl );
        *nl = 0;
    }
    if ( *nl1 ) {
        FreeNeighList( *nl1 );
        *nl1 = 0;
    }
    if ( *nl2 ) {
        FreeNeighList( *nl2 );
        *nl2 = 0;
    }
    if ( *nVisited1 ) {
        inchi_free( *nVisited1 );
        *nVisited1 = NULL;
    }
    if ( *nVisited2 ) {
        inchi_free( *nVisited2 );
        *nVisited2 = NULL;
    }

}
/**************************************************************************************/
int AllocateForNonStereoRemoval( sp_ATOM *at, int num_atoms, const AT_RANK *nSymmRank, AT_RANK *nCanonRank,
                            AT_RANK **nAtomNumberCanon1, AT_RANK **nAtomNumberCanon2,
                            NEIGH_LIST **nl, NEIGH_LIST **nl1, NEIGH_LIST **nl2, AT_RANK **nVisited1, AT_RANK **nVisited2 )
{
    DeAllocateForNonStereoRemoval( nAtomNumberCanon1, nAtomNumberCanon2, nl, nl1, nl2, nVisited1, nVisited2 );
    *nAtomNumberCanon1 = (AT_RANK *) inchi_malloc( num_atoms * sizeof(**nAtomNumberCanon1) );
    *nAtomNumberCanon2 = (AT_RANK *) inchi_malloc( num_atoms * sizeof(**nAtomNumberCanon2) );
    *nl                = CreateNeighList( num_atoms, num_atoms, at, 0, NULL );
    *nl1               = CreateNeighList( num_atoms, num_atoms, at, 0, NULL );
    *nl2               = CreateNeighList( num_atoms, num_atoms, at, 0, NULL );
    *nVisited1         = (AT_RANK *) inchi_malloc( num_atoms * sizeof(**nVisited1) );
    *nVisited2         = (AT_RANK *) inchi_malloc( num_atoms * sizeof(**nVisited2) );

    if ( !*nl || !*nl1 || !*nl2 || !*nVisited1 || !*nVisited2 || !*nAtomNumberCanon1 || !*nAtomNumberCanon2 ) {
        DeAllocateForNonStereoRemoval( nAtomNumberCanon1, nAtomNumberCanon2, nl, nl1, nl2, nVisited1, nVisited2 );
        return 0;
    }
    /*  Sort neighbors according to symm. ranks (primary key) and canon. ranks (secondary key), in descending order */
    SortNeighListsBySymmAndCanonRank( num_atoms, *nl,  nSymmRank, nCanonRank );
    SortNeighListsBySymmAndCanonRank( num_atoms, *nl1, nSymmRank, nCanonRank );
    SortNeighListsBySymmAndCanonRank( num_atoms, *nl2, nSymmRank, nCanonRank );
    return 1;
}
/**************************************************************************************/
AT_RANK GetMinNewRank(AT_RANK *nAtomRank, AT_RANK *nAtomNumb, AT_RANK nRank1 )
{
    int i;
    AT_RANK nRank2;
    for ( i = (int)nRank1-1; 0 <= i && nRank1 == (nRank2 = nAtomRank[(int)nAtomNumb[i]]); i -- )
        ;
    if ( i >= 0 )
        nRank2 ++;
    else
        nRank2 = 1;
    return nRank2;
}
/**************************************************************************************/
int BreakNeighborsTie(  sp_ATOM *at, int num_atoms, int num_at_tg, int ib, int ia,
                        AT_RANK *neigh_num, int in1, int in2, int mode,
                        AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList,
                        const AT_RANK *nSymmRank, AT_RANK *nCanonRank, NEIGH_LIST *nl1, NEIGH_LIST *nl2, long *lNumIter )
{
    AT_RANK nRank1, nRank2;
    int     nNumDiffRanks, nNumDiffRanks1, nNumDiffRanks2, i;
    int n1  = (int)neigh_num[in1];
    int n2  = (int)neigh_num[in2];
    int other_neigh[2], other_neig_ord[2], num_other_neigh;
    /*  asymmetric calculation */

    if ( (mode == MAP_MODE_S4  && in1) || /* for S4 we need only (in1,in2) = (0,1) (0,2) (0,3) pairs of neighbors */
         (mode != MAP_MODE_STD && at[ia].valence != MAX_NUM_STEREO_ATOM_NEIGH) ||
         (mode != MAP_MODE_STD && nSymmRank[n1]  != nSymmRank[n2]) ) {
        return 0;
    }
    /*  1. Create initial ranks from equivalence information stored in nSymmRank */
    memcpy( pRankStack1[0], nSymmRank, num_at_tg * sizeof(pRankStack1[0][0]) );
    pn_RankForSort = pRankStack1[0];
    tsort( pRankStack1[1], num_at_tg, sizeof(pRankStack1[1][0]), CompRanksOrd );
    nNumDiffRanks = SortedEquInfoToRanks( pRankStack1[0]/*inp*/, pRankStack1[0]/*out*/, pRankStack1[1], num_at_tg, NULL );
    
    /* other neighbors */
    num_other_neigh = 0;
    if ( at[ia].valence <= MAX_NUM_STEREO_ATOM_NEIGH && mode ) {
        for ( i = 0; i < at[ia].valence; i ++ ) {
            if ( i != in1 && i != in2 ) {
                other_neigh[num_other_neigh]    = (int)neigh_num[i];
                other_neig_ord[num_other_neigh] = i;
                num_other_neigh ++;
            }
        }
    }
    if ( (mode != MAP_MODE_STD && nSymmRank[other_neigh[0]] != nSymmRank[other_neigh[1]]) ||
         (mode == MAP_MODE_S4  && nSymmRank[n1]             != nSymmRank[other_neigh[1]]) ) {
        return 0;
    }

    /*  2. Fix at[ia] */
    if ( pRankStack1[0][ia] != nSymmRank[ia] ) {
        /*  at[ia] is constitutionally equivalent to some other atom. Fix at[ia]. */
        pRankStack1[0][ia] = nSymmRank[ia];
        nNumDiffRanks = DifferentiateRanksBasic( num_at_tg, NeighList,
                                     nNumDiffRanks, pRankStack1[0], nTempRank,
                                     pRankStack1[1], lNumIter, 1 );
    }
    /*  3. In case of a double bond/cumulene only: */
    /*     fix at[ib] -- the opposite double bond/cumulene atom */
    if ( ib < num_atoms ) {
        /*  find the smallest possible rank */
        nRank1 = pRankStack1[0][ib];
        nRank2 = GetMinNewRank(pRankStack1[0], pRankStack1[1], nRank1 );
        /*  if the rank is smaller than pRankStack1[0][ib] then fix at[ib] */
        if ( nRank2 != nRank1 ) {
            pRankStack1[0][ib] = nRank2;
            nNumDiffRanks = DifferentiateRanksBasic( num_at_tg, NeighList,
                                         nNumDiffRanks, pRankStack1[0], nTempRank,
                                         pRankStack1[1], lNumIter, 1 );
        }
    }
    
    /**************************************************************************************
     * Note: It may (or may not?) make sense to fix "other neighbors":
     *       in case of a stereo center fix neighbors other than n1, n2
     *       in case of a double bond/cumulene fix the opposite atom neighbors
     *       The ranks assigned to the other neighbors in case of their equivalence
     *       should be in the ascending order of their canonical ranks ????
     *       *** For now we do not fix other neighbors ***
     **************************************************************************************/

    /*  4. Check whether the neighbors still have equal ranks */
    if ( pRankStack1[0][n1] != pRankStack1[0][n2] ) {
        return 0; /*  the two neighbors are not constitutionally equivalent */
    }
    /*  5. Find new smallest possible rank for n1 and n2 */
    nRank1 = pRankStack1[0][n1];
    nRank2 = GetMinNewRank(pRankStack1[0], pRankStack1[1], nRank1 );

    /*  6. Copy the results to the 2nd eq. rank arrays */
    memcpy( pRankStack2[0], pRankStack1[0], num_at_tg * sizeof(pRankStack2[0][0]) );
    memcpy( pRankStack2[1], pRankStack1[1], num_at_tg * sizeof(pRankStack2[0][0]) );

    /*  7. Break neighbor tie: map n1(1) <--> n2(2) */
    pRankStack1[0][n1] = nRank2;
    nNumDiffRanks1 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                 nNumDiffRanks, pRankStack1[0], nTempRank,
                                 pRankStack1[1], lNumIter, 1 );
    
    pRankStack2[0][n2] = nRank2;
    nNumDiffRanks2 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                 nNumDiffRanks, pRankStack2[0], nTempRank,
                                 pRankStack2[1], lNumIter, 1 );

    if ( nNumDiffRanks1 != nNumDiffRanks2 ) {
        return -1; /*  <BRKPT> */
    }
    if ( mode == MAP_MODE_C2v || mode == MAP_MODE_C2 ) {
        /* Check for C2v reflection leading to parity inversion (mode=1) or C2 rotation (mode=2) */
        AT_RANK nRank10, nRank20;
        int     nn1, nn2;
        /*
         * C2v & C2: map
         * n1(1) <--> n2(2) -- at this point already done
         * n1(2) <--> n2(1) --> do at i = 0
         *
         * C2v: other neighbors must be unmoved: map
         * other_neigh[0](1) <--> other_neigh[0](2)
         * other_neigh[1](1) <--> other_neigh[1](2)
         *
         * C2:  other neighbors should be mapped on each other
         * other_neigh[0](1) <--> other_neigh[1](2)
         * other_neigh[1](1) <--> other_neigh[0](2)
         */
        for ( i = 0; i <= 2; i ++ ) {
            if ( i == 0 ) {
                /* C2v & C2. Map n2(1) <--> n1(2) */
                nn1     = n2;
                nn2     = n1;
            } else
            if ( mode == MAP_MODE_C2v ) {   /* was '=', pointed by WDI */
                /* i = 1 or 2
                 * C2v. Other neighbors must be unmoved: map
                 * i=1: other_neigh[0](1) <--> other_neigh[0](2)
                 * i=2: other_neigh[1](1) <--> other_neigh[1](2)
                 */
                nn1 = other_neigh[i-1]; /* 0 or 1 */
                nn2 = other_neigh[i-1]; /* 0 or 1 */
            } else
            if ( mode == MAP_MODE_C2 ) {  /* was '=', pointed by WDI */
                /* i = 1 or 2
                 * C2.  Other neighbors should be mapped on each other
                 * i=1: other_neigh[0](1) <--> other_neigh[1](2)
                 * i=2: other_neigh[1](1) <--> other_neigh[0](2)
                 */
                nn1 = other_neigh[i-1]; /* 0 or 1 */
                nn2 = other_neigh[2-i]; /* 1 or 0 */
            } else {
                return -1; /* program error */
            }
            /* map nn1(1) <--> nn2(2) */
            nRank10 = pRankStack1[0][nn1];
            nRank20 = pRankStack2[0][nn2];
            nRank1 = GetMinNewRank(pRankStack1[0], pRankStack1[1], nRank10 );
            nRank2 = GetMinNewRank(pRankStack2[0], pRankStack2[1], nRank20 );
            if ( nRank10 == nRank20 && nRank1 == nRank2 ) {
                if ( nRank10 == nRank1 ) {
                    ;/* atoms are already mapped */
                } else {
                    /* need additional mapping: ranks are not fixed yet */
                    pRankStack1[0][nn1] = nRank1;
                    nNumDiffRanks1 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                                 nNumDiffRanks, pRankStack1[0], nTempRank,
                                                 pRankStack1[1], lNumIter, 1 );
                    pRankStack2[0][nn2] = nRank2;
                    nNumDiffRanks2 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                                 nNumDiffRanks, pRankStack2[0], nTempRank,
                                                 pRankStack2[1], lNumIter, 1 );
                    if ( nNumDiffRanks1 != nNumDiffRanks2 ) {
                        return -1; /*  <BRKPT> */
                    }
                }
            } else {
                return 0;  /* mapping is not possible */
            }
        }
    }
    if ( mode == MAP_MODE_S4 ) {
        /* 
         *  Check for S4 reflection/rotation leading to parity inversion (mode=3)
         *
         * At this point n1(1) <--> n2(2) have been mapped and n1 has index in1 = 0
         * Below indexes in neigh_num[] are in brackets; [i] means neigh_num[i].
         * Numbers (#) in parentheses refer to pRankStack#
         *
         * in2=1: [0](1) <--> [1](2)  mapping has been done; add more mappings:
         *        [1](1) <--> [2](2)  [x]=[2]
         *        [2](1) <--> [3](2)  [y]=[3]
         *        [3](1) <--> [0](2)
         *        this will succeed if C2 axis crosses middle of [0]-[2] and [1]-[3] lines 
         *
         * in2=2: [0](1) <--> [2](2) mapping has been done; add more mappings:
         *        [2](1) <--> [3](2)  [x]=[3]
         *        [3](1) <--> [1](2)  [y]=[1]
         *        [1](1) <--> [0](2)
         *        this will succeed if C2 axis crosses middle of [0]-[3] and [1]-[2] lines
         *
         * in2=3: [0](1) <--> [3](2) mapping has been done; add more mappings:
         *        [3](1) <--> [1](2)  [x]=[1]
         *        [1](1) <--> [2](2)  [y]=[2]
         *        [2](1) <--> [0](2)
         *        this will succeed if C2 axis crosses middle of [0]-[1] and [2]-[3] lines
         *
         * In general:
         *        [in1](1) <--> [in2](2)
         *        [in2](1) <--> [x]  (2)  i=0
         *        [x]  (1) <--> [y]  (2)  i=1
         *        [y]  (1) <--> [in1](2)  i=2
         *
         *    in1=0    always
         *    ===== how to find x, y from in2 ====
         *    in2=1 => x,y = 2, 3  or [x] = other_neigh[0], [y] = other_neigh[1] 
         *    in2=2 => x,y = 3, 1  or [x] = other_neigh[1], [y] = other_neigh[0]
         *    in2=3 => x,y = 1, 2  or [x] = other_neigh[0], [y] = other_neigh[1]
         *    ====================================
         */
        AT_RANK nRank10, nRank20;
        int     nn1, nn2;
        for ( i = 0; i <= 2; i ++ ) {
            switch( i ) {
            case 0:  /* [in2](1) <--> [x](2);  */
                nn1 = n2;                    /* [in2] */
                nn2 = other_neigh[1-in2%2];  /* [x]   */
                break;
            case 1:  /* [x](1) <--> [y](2) */
                nn1 = other_neigh[1-in2%2];  /* [x]   */
                nn2 = other_neigh[  in2%2];  /* [y]   */
                break;
            case 2:
                nn1 = other_neigh[  in2%2];  /* [y]   */
                nn2 = n1;                    /* [in1] */
                break;
            default:
                return -1; /* program error */
            }
            /* map nn1(1) <--> nn2(2) */
            nRank10 = pRankStack1[0][nn1];
            nRank20 = pRankStack2[0][nn2];
            nRank1 = GetMinNewRank(pRankStack1[0], pRankStack1[1], nRank10 );
            nRank2 = GetMinNewRank(pRankStack2[0], pRankStack2[1], nRank20 );
            if ( nRank10 == nRank20 && nRank1 == nRank2 ) {
                if ( nRank10 == nRank1 ) {
                    ;/* atoms are already mapped */
                } else {
                    /* need additional mapping: ranks are not fixed yet */
                    pRankStack1[0][nn1] = nRank1;
                    nNumDiffRanks1 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                                 nNumDiffRanks, pRankStack1[0], nTempRank,
                                                 pRankStack1[1], lNumIter, 1 );
                    pRankStack2[0][nn2] = nRank2;
                    nNumDiffRanks2 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                                 nNumDiffRanks, pRankStack2[0], nTempRank,
                                                 pRankStack2[1], lNumIter, 1 );
                    if ( nNumDiffRanks1 != nNumDiffRanks2 ) {
                        return -1; /*  <BRKPT> */
                    }
                }
            } else {
                return 0;  /* mapping is not possible */
            }
        }
    }



#if ( BREAK_ONE_MORE_SC_TIE == 1 ) /* { */
    /* Check for a very highly symmetrical stereo center 12-06-2002 */
    if ( ib >= num_atoms && at[ia].valence == MAX_NUM_STEREO_ATOM_NEIGH ) {
        int num_eq;
        nRank1 = pRankStack1[0][n2];
        for ( i = 0, num_eq = 0; i < at[ia].valence; i ++ ) {
            num_eq += ( nRank1 == pRankStack1[0][at[ia].neighbor[i]]);
        }
        if ( num_eq == MAX_NUM_STEREO_ATOM_NEIGH-1 ) {
            for ( i = (int)nRank1-1; 0 <= i && nRank1 == (nRank2 = pRankStack1[0][(int)pRankStack1[1][i]]); i -- )
                ;
            if ( i >= 0 )
                nRank2 ++;
            else
                nRank2 = 1;

            /*  7a. Break another neighbor tie */

            nNumDiffRanks = nNumDiffRanks1;

            pRankStack1[0][n2] = nRank2;
            nNumDiffRanks1 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                         nNumDiffRanks, pRankStack1[0], nTempRank,
                                         pRankStack1[1], lNumIter, 1 );
    
            pRankStack2[0][n1] = nRank2;
            nNumDiffRanks2 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                         nNumDiffRanks, pRankStack2[0], nTempRank,
                                         pRankStack2[1], lNumIter, 1 );
        }
    }

    if ( nNumDiffRanks1 != nNumDiffRanks2 ) {
        return -1; /*  <BRKPT> */
    }
#endif /* } BREAK_ONE_MORE_SC_TIE */

#if ( BREAK_ALSO_NEIGH_TIE == 1 )
    /* check whether neighbor's neighbors are tied and untie them */
    if ( at[n1].nRingSystem == at[n2].nRingSystem &&  ib >= num_atoms ) {
        AT_RANK NeighNeighList[MAX_NUM_STEREO_ATOM_NEIGH+1];
        int m, neigh1=-1, neigh2=-1;
        nRank1 = nRank2 = 0;
        /* n1 */
        NeighNeighList[0] = at[n1].valence-1; /* for insertions_sort_NeighListBySymmAndCanonRank() */
        for ( i = 0, m = 1; i < at[n1].valence; i ++ ) {
            int neigh = at[n1].neighbor[i];
            if ( neigh != ia ) {
                NeighNeighList[m ++] = neigh;
            }
        }
        insertions_sort_NeighListBySymmAndCanonRank( NeighNeighList, pRankStack1[0], nCanonRank );
        for ( m = 2; m < at[n1].valence; m ++ ) {
            if ( pRankStack1[0][NeighNeighList[m]] == pRankStack1[0][NeighNeighList[m-1]] ) {
                neigh1 = NeighNeighList[m-1];
                break;
            }
        }
        /* n2 */
        NeighNeighList[0] = at[n2].valence-1; /* for insertions_sort_NeighListBySymmAndCanonRank() */
        for ( i = 0, m = 1; i < at[n2].valence; i ++ ) {
            int neigh = at[n2].neighbor[i];
            if ( neigh != ia ) {
                NeighNeighList[m ++] = neigh;
            }
        }
        insertions_sort_NeighListBySymmAndCanonRank( NeighNeighList, pRankStack2[0], nCanonRank );
        for ( m = 2; m < at[n2].valence; m ++ ) {
            if ( pRankStack2[0][NeighNeighList[m]] == pRankStack2[0][NeighNeighList[m-1]] ) {
#if ( BREAK_ALSO_NEIGH_TIE_ROTATE == 1 )
                neigh2 = NeighNeighList[m];    /* [m] to obtain same axis orientation  around ia<neigh */
#else
                neigh2 = NeighNeighList[m-1];  /* [m-1] to obtain reflection ??? */
#endif
                break;
            }
        }
        if ( neigh1 >= 0 && neigh2 >= 0 && pRankStack1[0][neigh1] == pRankStack2[0][neigh2] ) {
            /* neighbors' neighbors are tied */
            nRank1 = pRankStack1[0][neigh1];
            nRank2 = GetMinNewRank(pRankStack1[0], pRankStack1[1], nRank1 );

            /*  Break neighbor's neighbor tie */

            nNumDiffRanks = nNumDiffRanks1;

            pRankStack1[0][neigh1] = nRank2;
            nNumDiffRanks1 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                         nNumDiffRanks, pRankStack1[0], nTempRank,
                                         pRankStack1[1], lNumIter, 1 );

            pRankStack2[0][neigh2] = nRank2;
            nNumDiffRanks2 = DifferentiateRanksBasic( num_at_tg, NeighList,
                                         nNumDiffRanks, pRankStack2[0], nTempRank,
                                         pRankStack2[1], lNumIter, 1 );
        }
    }
#endif


    /*  for debug only */
    for ( i = 0; i < num_at_tg; i ++ ) {
        if ( pRankStack1[0][(int)pRankStack1[1][i]] != pRankStack2[0][(int)pRankStack2[1][i]] ) {
            return -1;  /*  <BRKPT> */
        }
    }
    /*  Resort lists of  neighbors */
    SortNeighListsBySymmAndCanonRank( num_atoms, nl1,  pRankStack1[0], nCanonRank );
    SortNeighListsBySymmAndCanonRank( num_atoms, nl2,  pRankStack2[0], nCanonRank );

    return nNumDiffRanks1+1;
}

/**************************************************************************************/
int CheckNextSymmNeighborsAndBonds( sp_ATOM *at, AT_RANK cur1, AT_RANK cur2, AT_RANK n1, AT_RANK n2,
                                    AT_RANK *nAvoidCheckAtom, AT_RANK *nVisited1, AT_RANK *nVisited2,
                                    AT_RANK *nVisitOrd1, AT_RANK *nVisitOrd2, const AT_RANK *nRank1, const AT_RANK *nRank2 )
{
    AT_RANK s1, s2;
    int     i1, i2, k1, k2;
    if ( nRank1[n1] != nRank2[n2] ) {
        return -1; /*  parallel traversal in stereo removal failed */ /*   <BRKPT> */
    }
    switch ( !nVisited1[n1] + !nVisited2[n2] ) {
    case 0:
        if ( nVisited1[n1] != n2+1 || nVisited2[n2] != n1+1 ) {
            return -1; /*  0; */ /*  possibly error???: we have come to an alreardy traversed pair and */
                       /*  found that the pair previously has not been traversed synchroneously. */
        }              /*  -- Happens in C60. */
        break;
    case 1:
        return -1; /*  0; */ /*  possibly error: one is zero, another is not a zero. Happens in C60 */

    /*  case 2: */
        /* both are zero, OK. */
    }
    
    if ( nVisitOrd1[n1] != nVisitOrd2[n2] ) {
        return -1; /*  0; */ /*  different DFS trees */
    }
    /*  at[n1] and at[n2] are next to at[cur1] and at[cur2] respectively */
    /*  Even though the bond might have already been checked, check whether */
    /*  it is a stereo bond/cumulene. If it is, check the bond/cumulene parity. */

    /*  Even though the bond or cumulene might have already been checked, check it: this is */
    /*  the only place we can check stereo bonds and cumulenes that are not edges of the DFS tree */
    /*  The code works both for a stereo bond and a stereogenic cumulene. */

    for ( i1 = 0, k1 = 0; i1 < MAX_NUM_STEREO_BONDS &&
                          (s1=at[cur1].stereo_bond_neighbor[i1]) &&
                         !(k1=(at[cur1].neighbor[(int)at[cur1].stereo_bond_ord[i1]] == n1)); i1 ++ )
        ;
    for ( i2 = 0, k2 = 0; i2 < MAX_NUM_STEREO_BONDS &&
                          (s2=at[cur2].stereo_bond_neighbor[i2]) &&
                         !(k2=(at[cur2].neighbor[(int)at[cur2].stereo_bond_ord[i2]] == n2)); i2 ++ )
        ;

    /* -- this does not work in case of cumulenes --
    for ( i1 = 0, k1 = 0; i1 < MAX_NUM_STEREO_BONDS && (s1=at[cur1].stereo_bond_neighbor[i1]) && !(k1=(s1-1 == n1)); i1 ++ )
        ;
    for ( i2 = 0, k2 = 0; i2 < MAX_NUM_STEREO_BONDS && (s2=at[cur2].stereo_bond_neighbor[i2]) && !(k2=(s2-1 == n2)); i2 ++ )
        ;
    */

    if ( k1 != k2 ) {    
        return 0; /*  not an error: a stereo bond and not a stereo bond */
    }
    if ( k1 ) {
        /* here k1 == k2 */
        int bCheckBond1, bCheckBond2;
        s1 --;
        s2 --;

        bCheckBond1 = (cur1 != nAvoidCheckAtom[0] || s1 != nAvoidCheckAtom[1]) &&
                      (cur1 != nAvoidCheckAtom[1] || s1 != nAvoidCheckAtom[0]);
        bCheckBond2 = (cur2 != nAvoidCheckAtom[0] || s2 != nAvoidCheckAtom[1]) &&
                      (cur2 != nAvoidCheckAtom[1] || s2 != nAvoidCheckAtom[0]);
        
        if ( bCheckBond1 != bCheckBond2 )
            return 0;
        
        if ( !bCheckBond1 && !bCheckBond2 ) {
            return 1; /*  do not go any further in this direction */
        }

        if ( at[cur1].stereo_bond_parity[i1] != at[cur2].stereo_bond_parity[i2] ) {
            /*  different values of  at[].stereo_bond_parity: definitely different bonds */
            /*  known parities */
            if ( PARITY_KNOWN(at[cur1].stereo_bond_parity[i1] ) &&
                 PARITY_KNOWN(at[cur2].stereo_bond_parity[i2] )  ) {
                return 0; /*  different currently known stereo bond parities */
            }
#if ( PROPAGATE_ILL_DEF_STEREO != 1 )
            /*  well defined and to be calculated from the ranks */
            if ( !(PARITY_CALCULATE(at[cur1].stereo_bond_parity[i1]) && PARITY_WELL_DEF (at[cur2].stereo_bond_parity[i2]) ||
                   PARITY_WELL_DEF (at[cur1].stereo_bond_parity[i1]) && PARITY_CALCULATE(at[cur2].stereo_bond_parity[i2]) ||
                   PARITY_CALCULATE(at[cur1].stereo_bond_parity[i1]) && PARITY_CALCULATE(at[cur2].stereo_bond_parity[i2]) ) ) {
                /*  do not reject if: "well defined" and "calculate" or "calculate" and "calculate" */
                return 0; 
            }
#endif
        }
   
#if ( PROPAGATE_ILL_DEF_STEREO != 1 )
        if ( (cur1 != cur2 || s1 != s2) && (cur1 != s2 || cur2 != s1) ) {
            /*  two different stereo bonds */
            if ( PARITY_ILL_DEF( at[cur1].stereo_bond_parity[i1] ) ||
                 PARITY_ILL_DEF( at[cur2].stereo_bond_parity[i2] ) ) {
                return 0;
            }
        }
#endif
    }
    return 1; /*  stereo bonds to n1 and n2 have same known parities or are not stereo bonds */
}
/**************************************************************************************/
int CreateCheckSymmPaths( sp_ATOM *at, AT_RANK prev1, AT_RANK cur1, AT_RANK prev2, AT_RANK cur2,
                         AT_RANK *nAvoidCheckAtom, AT_RANK *nVisited1, AT_RANK *nVisited2,
                         AT_RANK *nVisitOrd1, AT_RANK *nVisitOrd2,
                         NEIGH_LIST *nl1, NEIGH_LIST *nl2, const AT_RANK *nRank1, const AT_RANK *nRank2,
                         AT_RANK *nCanonRank, AT_RANK *nLength, int *bParitiesInverted, int mode  )
{
    int k, k1, k2, ret=0, bParitiesInvertedZero=0, *pbParitiesInverted;
    AT_RANK n1, n2;

    nVisited1[cur1] = cur2+1;  /*  symmetrically exchange atom numbers */
    nVisited2[cur2] = cur1+1;

    (*nLength) ++;

    nVisitOrd1[cur1] = *nLength; /*  save DFS visit order */
    nVisitOrd2[cur2] = *nLength;

    /* new version allows all inverted parities */
    if ( PARITY_WELL_DEF(at[cur1].stereo_atom_parity) &&
         PARITY_WELL_DEF(at[cur2].stereo_atom_parity) ) {
        if ( *bParitiesInverted < 0 ) {
            *bParitiesInverted = (at[cur1].stereo_atom_parity + at[cur2].stereo_atom_parity) % 2;
        } else
        if ( *bParitiesInverted != (at[cur1].stereo_atom_parity + at[cur2].stereo_atom_parity) % 2 ) {
            return 0; /*  Different known in advance parities have wrong "inverted" relation */
        }
    } else
    if ( PARITY_KNOWN(at[cur1].stereo_atom_parity) &&
         PARITY_KNOWN(at[cur2].stereo_atom_parity) &&
         at[cur1].stereo_atom_parity != at[cur2].stereo_atom_parity ) {
        return 0;  /*  Different known in advance parities */
    }

    if ( cur1 != cur2 &&
         !at[cur1].stereo_bond_neighbor[0] && !at[cur2].stereo_bond_neighbor[0] &&
         PARITY_KNOWN(at[cur1].parity) != PARITY_KNOWN(at[cur2].parity) ) {
        return 0; /*  one atom is stereogenic, another (presumably equivalent) is not. 9-11-2002 */
    }
#if ( PROPAGATE_ILL_DEF_STEREO != 1 )
    if ( cur1 != cur2 &&
         (PARITY_ILL_DEF(at[cur1].stereo_atom_parity) ||
          PARITY_ILL_DEF(at[cur2].stereo_atom_parity)) 
       ) {
        return 0;  /*  Cannot detect whether the paths are same or different */
    }
#endif

    if ( at[cur1].valence != at[cur2].valence ) {
        return CT_REMOVE_STEREO_ERR; /*  program error */ /*   <BRKPT> */
    }
    if ( at[cur1].valence == 1 ) {
        return 1; /*  so far success */
    }
    
    if ( nl1[(int)cur1][0] != nl2[(int)cur2][0] || nl1[(int)cur1][0] != at[cur1].valence ) {
        return CT_REMOVE_STEREO_ERR; /*  error: different valences */ /*   <BRKPT> */
    }


    for ( k = 1, k1 = 1, k2 = 1; k < at[cur1].valence; k ++, k1 ++, k2 ++ ) {
        if ( (n1 = nl1[(int)cur1][k1]) == prev1 ) {
            n1 = nl1[(int)cur1][++k1]; /*  don't go back */
        }
        if ( (n2 = nl2[(int)cur2][k2]) == prev2 ) {
            n2 = nl2[(int)cur2][++k2]; /*  don't go back */
        }
        
        if ( 0 >= (ret = CheckNextSymmNeighborsAndBonds( at, cur1, cur2, n1, n2, nAvoidCheckAtom,
                                              nVisited1, nVisited2, nVisitOrd1, nVisitOrd2, nRank1, nRank2 ) ) ) {
            return ret; /*  different neighbors or bonds                       */
        }

        if ( !nVisited1[n1] ) { /*  recursion */
            /* allow all inverted parities only inside a single ring system containing the starting point */
            pbParitiesInverted = (at[cur1].nRingSystem == at[n1].nRingSystem)? bParitiesInverted:&bParitiesInvertedZero;
            if ( 0 >= (ret = CreateCheckSymmPaths( at, cur1, n1, cur2, n2, nAvoidCheckAtom,
                                         nVisited1, nVisited2, nVisitOrd1, nVisitOrd2,
                                         nl1, nl2, nRank1, nRank2, nCanonRank, nLength, pbParitiesInverted, mode ) ) ) {
                return ret;
            }
        }
    }
    return 1; /*  Success */

}
/**************************************************************************************/
/*  Compare parities */
#define MAX_OTHER_NEIGH        2
/*  nNeighMode */
#define NEIGH_MODE_RING        1
#define NEIGH_MODE_CHAIN       2

#define CHECKING_STEREOCENTER  1
#define CHECKING_STEREOBOND    2

#define COMP_STEREO_SUCCESS    1
#define NOT_WELL_DEF_UNKN      2
#define NOT_WELL_DEF_UNDF      4

#define PARITY_IMPOSSIBLE    999
/**************************************************************************************
  Note:    the following C2v/S4 stereo center symmetry recognition
           is not included in the final InChI version released in April 2005
           It is disabled in the mode.h (CHECK_C2v_S4_SYMM = 0)
           As the result, the only central atom in S4 or atoms on C2v axis
           may have pasrity (-) even though these atoms are not stereogenic.

  Reason:  Not finished/tested yet
 **************************************************************************************

  In case of stereocenter with 2 pairs of constitutionally identical neighbors :
  
  G(n) > H(m) means group G has n elements; group H has m elements and
                                            group H is a subgroup of G

  Td(24) > D2d(8> > D2(4)
                  > S4(4)  > C2(2) -- Test for S4
                  > C2v(4) > C2(2) -- Test for C2v
                           > Cs(2)
                  
  Td(24) > C3v(6) > C3(3) -- does not have 2 pairs of constitutionally identical neighbors
                  > Cs(2)

  The pair of atoms to check for the existence of a steregenic atom: X, Y

       X   Y
        \ /
         C
        / \
       A   B

  Conditions to check:

  (a) Old #0: Map canonical numbers X1 <--> Y2
      Traverse DFS from X and Y
      If all parities vs. canon. numbers unchanged except that of C
      then C is not stereogenic

  (b) C2v  #1: discover ACB symmetry plain Cv
      o Map canonical numbers X1 <--> Y2, Fix(Ai), Fix(Bi)
      o Make sure that after mapping X1 <--> Y2 the atoms Ai and
        Bi still have equal mapping ranks
      Traverse DFS from X and Y
      In this case canonical numbers will be reflected in plane ACB if it exists.
      o Criterion of the presence of the symmetry plain is:
        --> all stereogenic atoms and allenes parities are inverted
  (c) C2v  #2: discover vertical axis C2
      o Map canonical numbers X1 <--> Y2 and A1 <--> B2
      o Make sure that after mapping X1 <--> Y2 the atoms Ai and
        Bi still have equal mapping ranks
      o Traverse DFS from X1 and Y2
        In this case canonical numbers will be rotated by
                     180 degrees around the vertical axis
         (this may be considered as a superposition of two Cv
          reflections in perpendicular vertical planes)
      o Criterion of the presence of the C2 axis is:
        --> all stereogenic atoms and allenes parities are not changed
  (d) S4  #3: discover axis horizontal S4 axis
      o Map canonical numbers X1 <--> Y2, Y1 <--> A2, A1 <--> B2, B1 <--> X2
      o Traverse DFS from X1 and Y2
        In this case the canonical numbers will be rotated by
                     90 degrees and reflected in a horizontal plane.
        3 attempts corrresponding to transpositions 0132, 0213, 0321
                                 are sufficient (XY=01,02,03)
      o Criterion of the presence of the S4 symmetry axis is:
        --> all stereogenic atoms and allenes parities are inverted

***************************************************************************************/

/**************************************************************************************/
int CalculatedPathsParitiesAreIdentical( sp_ATOM *at, int num_atoms, const AT_RANK *nSymmRank,
                         AT_RANK *nCanonRank, AT_RANK *nAtomNumberCanon,
                         AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2,
                         AT_RANK *nVisited1, AT_RANK *nVisited2,
                         AT_RANK prev_sb_neigh, AT_RANK cur, AT_RANK next1, AT_RANK next2, int nNeighMode,
                         int bParitiesInverted, int mode, CANON_STAT *pCS,
                         int vABParityUnknown)
{
    int i, i01, i02, i11, i12, i21, i22, k, parity, parity1, parity2, parity12, num_other_neigh;
    int nNumEqStereogenic, nCheckingMode, not_well_def_parities;
    AT_RANK other_neigh[MAX_NUM_STEREO_ATOM_NEIGH], neigh, r1, r2;
    int  nNumComparedCenters = 0, nNumComparedBonds = 0, bCurParityInv1=0 /*, bCurParityInv2=0*/;
    int  bCurRotated=0, nNumDiff=0, nNumInv=0;
    int  s1, s2;

    nCheckingMode = ( prev_sb_neigh < num_atoms )? CHECKING_STEREOBOND : CHECKING_STEREOCENTER;
    not_well_def_parities = 0;
    nNumEqStereogenic = 0;

    if ( (nNeighMode != NEIGH_MODE_RING &&
         bParitiesInverted != 0) || abs(bParitiesInverted) != 1 ) {
        bParitiesInverted = 0;
    }

    if ( bParitiesInverted ) {
        for ( i = 0, i11 = i22 = 0; i < num_atoms; i ++ ) {
            /* count number of atoms that have not been visited */
            i11 += !nVisited1[i];
            i22 += !nVisited2[i];
            nAtomNumberCanon1[i] = MAX_ATOMS+1;  /*  mark unchanged */
            nAtomNumberCanon2[i] = MAX_ATOMS+1;  /*  mark unchanged */
        }
        if ( i11 || i22 ) {
            if ( bParitiesInverted == 1 )
                return 0; /* only a part of the structure has been inverted */
            else
                bParitiesInverted = 0;
        }
    } else {
        for ( i = 0; i < num_atoms; i ++ ) {
            nAtomNumberCanon1[i] = MAX_ATOMS+1;  /*  mark unchanged */
            nAtomNumberCanon2[i] = MAX_ATOMS+1;  /*  mark unchanged */
        }
    }
    if ( (bParitiesInverted  > 0 && !(mode == MAP_MODE_C2v || mode == MAP_MODE_S4)) ||
         (bParitiesInverted == 0 && !(mode == MAP_MODE_C2  || mode == MAP_MODE_STD))) {
        return 0;
    }
    /**************************************************************************************
     *    The following discussion assumes that the canonical numbers are
     *    switched for some pairs of constitutionally identical atoms
     *    in such a way that the new numbering is an equivalent to the
     *    nCanonRank[] canonical numbering (the transposition belongs to the
     *    automorphism group of the chemical structure having no stereo).
     *    At this point non-zero elements of nVisited1[] and nVisited2[]
     *    together contain transposition P of the atom numbers.
     *    and P2 respectively of the ordering atom numbers: nVisitedi[k] = Pi(k)+1;
     *    In this implementation:
     *       P1(k)=k for all k
     *       P2(cur)=cur, P2(next1)=next2, P2(next2)=next1 
     *
     *    Below we call one of the numberings "old", another "new".
     *
     *    *IF* the old and the new canonical numberings produce same parities for stereogenic
     *    elements for the same canonical number(s)
     *    (that is, old_parity(canon_number) == new_parity(canon_number)
     *    *except* the currently being tested stereocenter at[cur] or stereobond/cumulene
     *    at[cur]=at[prev_sb_neigh], whose parity MUST be inverted
     *
     *    *THEN* the stereocenter or stereobond/cumulene is not stereogenic with one
     *
     *    *EXCEPTION* If the currently tested stereogenic element is constitutionally
     *    equivalent to two or more other stereogenic elements that have been
     *    permuted then the currently tested one is still stereogenic.
     **************************************************************************************/

     /*
     * 1. replace the assigned in each of the parallel traversals atom numbers
     *    with the canon. ranks corresponding to the atom numbers in the
     *    currently numbered atoms at[].
     *    One of obtained this way canonical numberings (probably nVisited1[])
     *    is same as the nCanonRank[] because usually nVisited1[i] = i+1 or 0
     */
    for ( i = 0; i < num_atoms; i ++ ) {

        if ( nVisited1[i] ) {
            /* canonical number of the atom mapped on atom #i in 'left' path */
            nVisited1[i] = nCanonRank[ (int)nVisited1[i] - 1 ];
            /* reverse: atom # from the mapped canonical rank in 'left' path */
            nAtomNumberCanon1[nVisited1[i] - 1] = i;
        }
        if ( nVisited2[i] ) {
            /* canonical number of the atom mapped on atom #i in 'right' path */
            nVisited2[i] = nCanonRank[ (int)nVisited2[i] - 1 ];
            /* reverse: atom # from the mapped canonical rank in 'right' path */
            nAtomNumberCanon2[nVisited2[i] - 1] = i;
        }
        /* if 'left' and 'right' path do not have atoms in common except the
           starting atom (and in case of stereobond, the end atom) some of
           nVisitedi[i] elements may be zero.
        */
    }
    
    /*
     * if started with a stereobond then check whether its parity has changed.
     * If yes then continue, otherwise parities are different
     *
     * if started with a stereo center then prev_sb_neigh = MAX_ATOMS+1
     *
     * If the transposition of next1 and next2 changes only the parity of the starting stereo atom or stereo bond
     * then the stereo bond or stereo atom is not stereogenic
     *
     * The exception: the stereogenic elememt in question is equivalent
     *    to two or more traversed other stereogenic elememts
     *    (see nNumEqStereogenic below, case similar to trimethylcyclopropane:
     *     3 or more constitutionally equivalent stereogenic elements)
     */
    if ( nCheckingMode == CHECKING_STEREOBOND ) {
        /******************************************************************************
         *
         *  Possibly stereogenic starting bond or cumulene at[cur]-at[prev_sb_neigh]
         *
         *******************************************************************************/
        /*  checking the starting stereo bond */
        if ( nVisited1[prev_sb_neigh] || nVisited2[prev_sb_neigh] ) {
            /*  the bond or cumulene is in the ring and the opposite atom has been visited */
            if ( nVisited1[prev_sb_neigh]  != nVisited2[prev_sb_neigh] ||
                 nCanonRank[prev_sb_neigh] != nVisited2[prev_sb_neigh] ) {
                return 0; /*  error: we came back to the same bond/cumulene and */ /*   <BRKPT> */
                          /*  assigned different canon. ranks to the opposite atom. */
            }
            if ( at[prev_sb_neigh].valence + at[prev_sb_neigh].num_H > 3 )
                return 0; /*  at[prev_sb_neigh] atom can not be adjacent to a stereo bond/cumulene */
                          /*  or does not have 3 attachments (hydrogens are not considered here) */
            for ( i = 0, k = 0; i < MAX_NUM_STEREO_BONDS &&
                                  (neigh=at[prev_sb_neigh].stereo_bond_neighbor[i]) && !(k=(neigh-1 == cur)); i ++ )
                ;
            if ( !k ) {
                return -1; /*  program error: could not locate stereogenic bond mark on the opposite atom */
            }
            k = (int)at[prev_sb_neigh].stereo_bond_ord[i]; /*  seq. number of the double or cumulene bond on at[prev_sb_neigh] */

            for ( i = 0, num_other_neigh = 0; i < at[prev_sb_neigh].valence && num_other_neigh <= MAX_OTHER_NEIGH; i ++ ) {
                if ( i != k ) { /*  do not include the double or cumulene bond */
                    other_neigh[num_other_neigh ++] = at[prev_sb_neigh].neighbor[i];
                }
            }
            if ( num_other_neigh + at[prev_sb_neigh].num_H > MAX_OTHER_NEIGH ) {
                return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
            }
            for ( i = 0; i < num_other_neigh; i ++ ) {
                k = (int)other_neigh[i];
                if ( nVisited1[k] && nVisited1[k] != nCanonRank[k] ) {
                    return 0; /*  parity of the statring stereo bond/cumulene has not changed. */
                }
                if ( nVisited2[k] && nVisited2[k] != nCanonRank[k] ) {
                    return 0; /*  parity of the statring stereo bond/cumulene has not changed. */
                }
            }
        }
    }
    if ( nCheckingMode == CHECKING_STEREOCENTER ) {
        /**************************************************
         *
         *  Possibly stereogenic starting atom at[cur]
         *
         **************************************************/
        /*  checking the starting stereo center */
        for ( i = 0, num_other_neigh = 0; i < at[cur].valence && num_other_neigh <= MAX_OTHER_NEIGH; i ++ ) {
            neigh  = at[cur].neighbor[i];
            if ( neigh != next1 && neigh != next2 ) {
                other_neigh[num_other_neigh ++] = neigh;
            }
        }
        if ( num_other_neigh + at[cur].num_H > MAX_OTHER_NEIGH ) {
            return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
        }
        /*
        if ( bParitiesInverted && at[cur].valence == MAX_NUM_STEREO_ATOM_NEIGH ) {
            if ( nVisited1[other_neigh[0]] == nCanonRank[other_neigh[0]] ||
                 nVisited2[other_neigh[0]] == nCanonRank[other_neigh[0]] ||
                 nVisited1[other_neigh[1]] == nCanonRank[other_neigh[1]] ||
                 nVisited2[other_neigh[1]] == nCanonRank[other_neigh[1]] ) {
                bParitiesInverted = 0;
                bCurRotated = 1;
            }
        }
        */
        /* bParitiesInverted = -1 means no predefined stereocenter has been checked */
        if ( bParitiesInverted && at[cur].valence == MAX_NUM_STEREO_ATOM_NEIGH ) {
            /* special case: 4 canonically eq. neighbors */
            int canon_parity, parity_vis_1, parity_vis_2;
            canon_parity = GetPermutationParity( at+cur, MAX_ATOMS+1, nCanonRank );
            parity_vis_1 = GetPermutationParity( at+cur, MAX_ATOMS+1, nVisited1 );
            parity_vis_2 = GetPermutationParity( at+cur, MAX_ATOMS+1, nVisited2 );
            if ( parity_vis_1 != parity_vis_2 ) {
                return 0;
            }
            if ( bParitiesInverted ==  1 && parity_vis_1 == canon_parity ) {
                return 0; /* not a typical case of inversion during the mapping of D4h stereocenter */
            } else
            if ( bParitiesInverted == -1 ) {
                if ( parity_vis_1 == canon_parity ) {
                    bParitiesInverted = 0;
                } else {
                    bParitiesInverted = 1;
                }
            }
        }
        /* at this point bParitiesInverted >= 0 */
        if ( !bParitiesInverted && !bCurRotated ) {
            for ( i = 0; i < num_other_neigh; i ++ ) {
                k = (int)other_neigh[i];
                if ( nVisited1[k] && nVisited1[k] != nCanonRank[k] ) {
                    return 0; /*  parity of the statring stereo center has not changed. */
                }
                if ( nVisited2[k] && nVisited2[k] != nCanonRank[k] ) {
                    return 0; /*  parity of the statring stereo center has not changed. */
                }
            }
        }
    }
    
    /*****************************************************
     * Check other (non-starting) stereo centers
     ******************************************************/
    for ( i = 0; i < pCS->nLenLinearCTStereoCarb; i ++, nNumComparedCenters += (k > 0) ) {
        r1     = pCS->LinearCTStereoCarb[i].at_num;
        i01    = nAtomNumberCanon[r1-1]; /*  ord. number of the atom that has canon rank r1 */
        
        i11     = nAtomNumberCanon1[r1-1]; /*  = (MAX_ATOMS+1) > num_atoms if the atom has not been traversed */
        i12     = nAtomNumberCanon2[r1-1]; /*  = otherwise < num_atoms */

        s1 = (i11 < num_atoms); /*  1 => the center was traversed on path #1 */
        s2 = (i12 < num_atoms); /*  1 => the center was traversed on path #2 */

        bCurParityInv1 = (bParitiesInverted &&
                          at[cur].nRingSystem == at[i11].nRingSystem &&
                          at[cur].nRingSystem == at[i12].nRingSystem );


        k  = 0;
        
        /*  check whether the two stereo centers (they can be one and the same atom) have been traversed */
        if ( !s1 && !s2 ) {
            continue;  /*  Both stereo centers have not been traversed; check the next pair. */
        }

        if ( nCheckingMode == CHECKING_STEREOCENTER ) {
            /*  check whether the stereocenters are the starting stereocenter */
            switch( (cur == i11) + (cur == i12) ) {
            case 2:
                continue; /*  do not recheck the starting atom */
            case 1:
                return -1; /*  possibly program error */ /*   <BRKPT> */
            /* case 0: */
            /*     break;  */  /*  the stereo centers are not the sarting stereo center */
            }
            if ( cur == i01 ) {
                return -1;  /*  program error: in this case at least one of the i11, i12 must be == cur */ /*   <BRKPT> */
            }
        }
        
        if ( nNeighMode == NEIGH_MODE_RING ) {
            if ( i11 != i12 && !bCurParityInv1 ) {
                return -1; /*  failed: the two stereo atoms have not been traversed synchronously */
            }
            if ( !at[i11].parity || !at[i12].parity ) {
                return 0; /*  another atom does not have parity (it might have been removed) 9-11-2002 */
            }
        }
        if ( nNeighMode == NEIGH_MODE_CHAIN ) {
            if ( s1+s2 != 1 ) {
                return -1; /*  program error: only one out of s1 and s2 must be 1, another must be 0. */
            }
            if ( (s1 && !at[i11].parity) || (s2 && !at[i12].parity) ) {
                return 0; /*  another atom does not have parity (it might have been removed) 9-11-2002 */
            }
        }
        
        parity  = pCS->LinearCTStereoCarb[i].parity;
        if ( (nNeighMode == NEIGH_MODE_RING  && (i11 != i01) && (i12 != i01)) ||
             /*  in NEIGH_MODE_RING case we know that i11 == i12 except bCurParityInv1 == 1 */
             nNeighMode == NEIGH_MODE_CHAIN 
             /*  in NEIGH_MODE_CHAIN case here we always have 2 different atoms */
        ) {
            /****************************************************************
             * Case of two transposed atoms or a circular permutation in D4h
             */
            parity1 = s1? GetStereoCenterParity( at, i11, nVisited1 ) : PARITY_IMPOSSIBLE;
            parity2 = s2? GetStereoCenterParity( at, i12, nVisited2 ) : PARITY_IMPOSSIBLE;
            if ( !ATOM_PARITY_KNOWN(parity1) && !ATOM_PARITY_KNOWN(parity2) ) {
                return -1; /*  should not happen: must have been detected at the time of the traversal */
            }
            if ( s1 && s2 ) {
                if ( bCurParityInv1 ) {
                    int parity1orig = GetStereoCenterParity( at, i11, nCanonRank ); 
                    int parity2orig = GetStereoCenterParity( at, i12, nCanonRank ); 
                    if ( i11 == i12 ||
                         ((parity1 == parity1orig || parity2 == parity2orig || parity1 != parity2) &&
                         ATOM_PARITY_WELL_DEF(parity1)) ||
                         (parity1 != parity2 && (!ATOM_PARITY_WELL_DEF(parity1) ||
                                                !ATOM_PARITY_WELL_DEF(parity2))) )
                        /*return -1; */ /* should be different atoms with inverted parities */
                        nNumDiff ++;
                } else {
                    if ( i11 != i12 || parity1 != parity2 )
                        return -1; /*  program error: must be the same atom */
                }
            }
            parity12 = s1? parity1 : parity2;

            if ( ATOM_PARITY_WELL_DEF(parity) && parity == parity12 ) {
                /*  symmetrical neighbors have well-defined equal parities */
                k ++;
                if ( nCheckingMode == CHECKING_STEREOCENTER && nNeighMode == NEIGH_MODE_RING ) {
                    /*  all 3: cur, i01, i11 are different atoms (here i11==i12) */
                    /*  here nSymmRank[i01]==nSymmRank[i11] due to the parallel traversal */
                    if ( nSymmRank[cur] == nSymmRank[i01] ) {
                        nNumEqStereogenic ++;  /*  all 3 are equ */
                    }
                }
            } else
            if ( ATOM_PARITY_WELL_DEF(parity) && ATOM_PARITY_WELL_DEF(parity12) ) {
                /*  apparently different well-defined parities */
                if ( !bCurParityInv1 ) {
                    nNumInv ++;
                    /* return 0; */
                }
            } else {
#if ( PROPAGATE_ILL_DEF_STEREO == 1 )
                /*  at least one parity is ill-defined. Use parity1 and parity2 to temporarily save bitmaps */
                parity1 = (parity  ==vABParityUnknown /*AB_PARITY_UNKN*/)? NOT_WELL_DEF_UNKN :
                          (parity  ==AB_PARITY_UNDF)? NOT_WELL_DEF_UNDF : 0;
                parity2 = (parity12==vABParityUnknown /*AB_PARITY_UNKN*/)? NOT_WELL_DEF_UNKN :
                          (parity12==AB_PARITY_UNDF)? NOT_WELL_DEF_UNDF : 0;
                if ( parity1 | parity2 ) {
                    not_well_def_parities |= ( parity1 | parity2 );
                    k ++;
                } else {
                    return -1;  /*  program error */ /*   <BRKPT> */
                }
#else
                return 0;
#endif
            }
        } else
        if ( i11 == i01 && i12 == i01 ) {
            /********************************************************************/
            /*  i11 == i12 are same atom as i01, nNeighMode == NEIGH_MODE_RING */
            if ( !s1 || !s2 ) {
                return -1;
            }
            /*  the parity of the new neighbors permutation must be same as the old one */
            /*  this must work for well-defined and ill-defined parities. */
            /*  actual parity (that includes the geometry) is not important here. */
            /*  old permutation */
            parity  = GetPermutationParity( at+i01, MAX_ATOMS+1, nCanonRank );
            /*  new parmutation */
            parity1 = GetPermutationParity( at+i01, MAX_ATOMS+1, nVisited1 );
            parity2 = GetPermutationParity( at+i01, MAX_ATOMS+1, nVisited2 );
            if ( parity != parity1 || parity != parity2 ) {
                return 0;
            }
            k ++;
        } else {
            /* nNeighMode == NEIGH_MODE_RING and only one out of the two (i11 == i01) (i12 == i01) is true */
            return -1; 
        }
        /* nNumComparedCenters += (k > 0); */
    }
    if ( bCurRotated || nNumDiff || nNumInv ) {
        return 0;
    }

     /* !!!! Add here bParitiesInverted == 1 case !!!! */
    /******************************************************/
    /*  Check other (non-starting) stereo bonds/cumulenes */
    /******************************************************/
    for ( i = 0; i < pCS->nLenLinearCTStereoDble; i ++, nNumComparedBonds += (k > 0) ) {
        r1     = pCS->LinearCTStereoDble[i].at_num1;
        r2     = pCS->LinearCTStereoDble[i].at_num2;
        i01    = nAtomNumberCanon[r1-1];  /*  ord. number of the atom that originally has canon rank r1 */
        i02    = nAtomNumberCanon[r2-1];  /*  ord. number of the atom that originally has canon rank r2 */

        i11     = nAtomNumberCanon1[r1-1]; /*  ord. number of the atom that got canon rank r1 during the parallel traversal */
        i12     = nAtomNumberCanon1[r2-1]; /*  ord. number of the atom that got canon rank r2 during the parallel traversal */

        i21     = nAtomNumberCanon2[r1-1];
        i22     = nAtomNumberCanon2[r2-1];


        s1 = (i11 < num_atoms && i12 < num_atoms);
        s2 = (i21 < num_atoms && i22 < num_atoms);

        k  = 0;

        /*  check whether the two stereo bonds/allenes (they can be one and the same) have been traversed */
        if ( !s1 && !s2 ) {
            continue; /*  Both stereo bonds/cumulenes have not been traversed; check the next pair. */
        }

        if ( nCheckingMode == CHECKING_STEREOBOND ) {
            switch ( ((i11 == cur && i12 == prev_sb_neigh) || (i12 == cur && i11 == prev_sb_neigh)) +
                     ((i21 == cur && i22 == prev_sb_neigh) || (i22 == cur && i21 == prev_sb_neigh)) ) {
            case 2:
                continue; /*  do not recheck the starting bond/cumulene */
            case 1:
                return -1; /*  possibly program error  */ /*   <BRKPT> */
            /* case 0: */
            /*     break; */   /*  the stereo centers are not the sarting stereo center */
            }
            if ( (i01 == cur && i02 == prev_sb_neigh) || (i02 == cur && i01 == prev_sb_neigh) ) {
                return -1;  /*  program error: in this case at least one of the i1x, i2x must be == cur */ /*   <BRKPT> */
            }
        }

        if ( nNeighMode == NEIGH_MODE_RING ) {
            if ( (i11 != i21 || i12 != i22) && (i11 != i22 || i12 != i21) ) {
                return -1; /*  failed: the two bonds/cumulenes have not been traversed synchronously */
            }
            if ( 0 > GetStereoNeighborPos( at, i11, i12 ) ) {
                return 0; /*  another bond is not stereo (the stereo might have been removed) 9-11-2002 */
            }

        }
        if ( nNeighMode == NEIGH_MODE_CHAIN ) {
            if ( s1+s2 != 1 ) {
                return -1; /*  program error: only one out of s1 and s2 must be 1, another must be 0. */
            }
            if ( (s1 && 0 > GetStereoNeighborPos( at, i11, i12 )) ||
                 (s2 && 0 > GetStereoNeighborPos( at, i21, i22 )) ) {
                return 0; /*  another bond is not stereo (the stereo might have been removed) 9-11-2002 */
            }
        }

        parity = pCS->LinearCTStereoDble[i].parity;
        /* bMustBeIdentical  = ATOM_PARITY_ILL_DEF(parity); */
        /* nNumEqStereogenic = 0; */

        if ( (nNeighMode == NEIGH_MODE_RING  && (i11 != i01 || i12 != i02) && (i11 != i02 || i12 != i01)) ||
             nNeighMode == NEIGH_MODE_CHAIN                    /*  in NEIGH_MODE_CHAIN case here we always have 2 different atoms */
        ) {
            /*******************************************/
            /*  case of two transposed bonds/cumulenes */
            parity1 = s1? GetStereoBondParity( at, i11, i12, nVisited1 ) : PARITY_IMPOSSIBLE;
            parity2 = s2? GetStereoBondParity( at, i21, i22, nVisited2 ) : PARITY_IMPOSSIBLE;
            if ( !ATOM_PARITY_KNOWN(parity1) && !ATOM_PARITY_KNOWN(parity2) ) {
                return -1; /*  should not happen: must have been detected at the time of traversal */
            }
            if ( s1 && s2 && (((i11 != i21 || i12 != i22) && (i11 != i22 || i12 != i21)) || parity1 != parity2 ) ) {
                return -1; /*  program error: must be the same bond/cumulene */
            }
            parity12 = s1? parity1 : parity2;
            if ( ATOM_PARITY_WELL_DEF(parity) && parity == parity12 ) {
                /*  symmetrical neighbors have well-defined equal parities */
                k ++;
                if ( nCheckingMode == CHECKING_STEREOBOND && nNeighMode == NEIGH_MODE_RING ) {
                    /*  all 3 bonds: cur-prev_sb_neigh, i01-i02, i11-i12 are different */
                    /*  (here <i11,i12>==<i21,i22> compared as unordered pairs) */
                    if ( (nSymmRank[cur] == nSymmRank[i01] && nSymmRank[prev_sb_neigh] == nSymmRank[i02]) ||
                         (nSymmRank[cur] == nSymmRank[i02] && nSymmRank[prev_sb_neigh] == nSymmRank[i01]) ) {
                        nNumEqStereogenic ++;
                    }
                }
            } else
            if ( ATOM_PARITY_WELL_DEF(parity) && ATOM_PARITY_WELL_DEF(parity12) ) {
                /*  apparently different well-defined parities */
                return 0;
            } else {
                /*  at least one parity is ill-defined. Use parity1 and parity2 to temporarily save bitmaps */
#if ( PROPAGATE_ILL_DEF_STEREO == 1 )
                parity1 = (parity  ==vABParityUnknown /*AB_PARITY_UNKN*/)? NOT_WELL_DEF_UNKN :
                          (parity  ==AB_PARITY_UNDF)? NOT_WELL_DEF_UNDF : 0;
                parity2 = (parity12==vABParityUnknown /*AB_PARITY_UNKN*/)? NOT_WELL_DEF_UNKN :
                          (parity12==AB_PARITY_UNDF)? NOT_WELL_DEF_UNDF : 0;
                if ( parity1 | parity2 ) {
                    not_well_def_parities |= ( parity1 | parity2 );
                    k ++;
                } else {
                    return -1;  /*  program error */
                }
#else
                return 0;
#endif
            }
        } else {
            /*****************************************************************************************/
            /*  i11-i12 and i21-i22 are same as i01-i02 bond/cumulene, nNeighMode == NEIGH_MODE_RING */
            AT_NUMB n1, n2;
            int       j;
            if ( !s1 || !s2 ) {
                return -1;
            }
            /*  find neighbors along the stereo bond/cumulene */
            for ( j = 0, n1 = MAX_ATOMS+1; j < MAX_NUM_STEREO_BOND_NEIGH && at[i01].stereo_bond_neighbor[j]; j ++ ) {
                if ( (int)at[i01].stereo_bond_neighbor[j] == i02+1 ) {
                    n1 = at[i01].neighbor[ (int)at[i01].stereo_bond_ord[j] ];
                    break;
                }
            }
            for ( j = 0, n2 = MAX_ATOMS+1; j < MAX_NUM_STEREO_BOND_NEIGH && at[i02].stereo_bond_neighbor[j]; j ++ ) {
                if ( (int)at[i02].stereo_bond_neighbor[j] == i01+1 ) {
                    n2 = at[i02].neighbor[ (int)at[i02].stereo_bond_ord[j] ];
                    break;
                }
            }
            if ( n1 > MAX_ATOMS || n2 > MAX_ATOMS ) {
                return CT_REMOVE_STEREO_ERR;
            }
            /*  the parity of the new neighbors permutation must be same as the old one */
            /*  this must work for well-defined and ill-defined parities. */
            /*  actual parity (that includes the geometry) is not important here. */
            /*  old permutation */
            parity  = GetPermutationParity( at+i01, n1, nCanonRank) + GetPermutationParity( at+i02, n2, nCanonRank);
            /*  new parmutation */
            parity1 = GetPermutationParity( at+i01, n1, nVisited1 ) + GetPermutationParity( at+i02, n2, nVisited1 );
            parity2 = GetPermutationParity( at+i01, n1, nVisited2 ) + GetPermutationParity( at+i02, n2, nVisited2 );
            if ( parity %2 != parity1 % 2 || parity1 % 2 != parity2 % 2 ) {
                return 0;
            }
            k ++;
        }

        /* nNumComparedBonds += ( k > 0 ); */
    }

    if ( nNumEqStereogenic > 0 ) {
        /*  case similar to trimethylcyclopropane: 3 constitutionally equivalent stereogenic elements */
        /*  the transposition does not change the parities */
#if ( bRELEASE_VERSION == 0 )
        pCS->bExtract |= EXTR_2EQL2CENTER_TO_REMOVE_PARITY;
#endif
        return 0;
    }
/* =========================================================================================
    Note
    ====
    At this point the comparison is complete and no difference sufficient to establish
    absence of stereo parity has been found.
    However, non-zero not_well_def_parities means that an ill-defined parity was
    compared to an ill-defined or well-defined parity. This means that the parity
    of the atom or bond being checked cannot be well-defined anymore.
   ========================================================================================*/


    not_well_def_parities |= COMP_STEREO_SUCCESS;

    return not_well_def_parities;

   /*  Add 1 to indicate success. The stereogenic elements might have been */
   /*  removed while checking existence of the previous atom/bond stereo */
   /* return (nNumComparedCenters + nNumComparedBonds + 1);  */
}
/********************************************************************************/
/*  Remove stereo marks from the bonds that are calculated to be non-stereo     */
/*  Such bonds must have 2 constitutionally equivalent attachments              */
/*  (can find two canonical numberings that change only one stereo bond parity) */
int RemoveCalculatedNonStereoBondParities( sp_ATOM *at, int num_atoms, int num_at_tg,
                                          AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList,
                                          AT_RANK *nCanonRank, const AT_RANK *nSymmRank,
                                          AT_RANK *nAtomNumberCanon, AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2,
                                          NEIGH_LIST *nl, NEIGH_LIST *nl1, NEIGH_LIST *nl2, 
                                          AT_RANK *nVisited1, AT_RANK *nVisited2, 
                                          CANON_STAT *pCS,
                                          int vABParityUnknown)
{
    int j, n, m, ret, ret1, ret2, ret_failed=0;
    
    int i1, n1, s2;  /*  n1 must be SIGNED integer */
    AT_RANK nAtomRank1, nAtomRank2, neigh[3], nAvoidCheckAtom[2], opposite_atom, nLength;
    int         nNeighMode = NEIGH_MODE_CHAIN;
    int         nNumEqRingNeigh = 0, bRingNeigh, bSymmNeigh, bParitiesInverted;
    NEIGH_LIST *nl01, *nl02;
    const AT_RANK    *nSymmRank1, *nSymmRank2;

    ret = 0;

second_pass:

    for ( i1 = 0; i1 < num_atoms && !RETURNED_ERROR(ret_failed); i1 ++ ) {
        if ( at[i1].valence != 3 || !at[i1].stereo_bond_neighbor[0] ) {
            continue;
        }
        for ( n1 = 0; n1 < MAX_NUM_STEREO_BONDS && !RETURNED_ERROR(ret_failed) && (s2=at[i1].stereo_bond_neighbor[n1]); n1++ ) {
            if ( !PARITY_CALCULATE(at[i1].stereo_bond_parity[n1]) && PARITY_WELL_DEF(at[i1].stereo_bond_parity[n1]) ) {
                continue;
            }
            opposite_atom = (AT_RANK)(s2-1);
            s2 = at[i1].neighbor[(int)at[i1].stereo_bond_ord[n1]]; /*  different from opposite_atom in case of a cumulene */
            for ( j = 1, n = 0; j <= (int)at[i1].valence; j ++ ) {      
                if ( nl[i1][j] != s2 ) {
                    neigh[n++] = nl[i1][j]; /*  sorting guarantees that canon. rank of neigh[0] is greater or equal */
                }
            }
            if ( n != 2 ) {
                ret = CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                goto exit_function;
            }
            if ( nSymmRank[(int)neigh[0]] != nSymmRank[(int)neigh[1]] ) {
                continue; /*  may happen if another half-bond has not a defined parity */
            }

            bRingNeigh = (at[(int)neigh[0]].nRingSystem == at[(int)neigh[1]].nRingSystem);
            switch ( nNeighMode ) {
            case NEIGH_MODE_CHAIN:
                if ( bRingNeigh ) {
                    nNumEqRingNeigh ++;
                    continue;
                }
                nl01 = nl;
                nl02 = nl;
                nSymmRank1 = nSymmRank;
                nSymmRank2 = nSymmRank;
                break;

            case NEIGH_MODE_RING:
                if ( !bRingNeigh )
                    continue;
                /*  break a tie between the two contitutionally equivalent neighbors, */
                /*  refine the two partitions, sort neighbors lists nl1, nl2 */
                bSymmNeigh = BreakNeighborsTie(  at, num_atoms, num_at_tg, opposite_atom, i1,
                                    neigh, 0, 1, 0,
                                    pRankStack1, pRankStack2, nTempRank, NeighList, nSymmRank, nCanonRank,
                                    nl1, nl2, &pCS->lNumNeighListIter );
                if ( bSymmNeigh <= 0 ) {
                    if ( ret_failed > bSymmNeigh )
                        ret_failed = bSymmNeigh;
                    continue;
                }
                nl01 = nl1;
                nl02 = nl2;
                nSymmRank1 = pRankStack1[0];
                nSymmRank2 = pRankStack2[0];
                break;
            default:
                return CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
            }

            /*  initialize arrays */
            memset( nVisited1, 0, sizeof(nVisited1[0])*num_atoms );
            memset( nVisited2, 0, sizeof(nVisited2[0])*num_atoms );
            memset( nAtomNumberCanon1, 0, sizeof(nAtomNumberCanon1[0])*num_atoms );
            memset( nAtomNumberCanon2, 0, sizeof(nAtomNumberCanon2[0])*num_atoms );
            nLength       = 1;
            nVisited1[i1] = i1+1;   /*  start atoms are the same */
            nVisited2[i1] = i1+1;
            nAtomNumberCanon1[i1] = nLength;
            nAtomNumberCanon2[i1] = nLength;
            nAvoidCheckAtom[0] = i1;
            nAvoidCheckAtom[1] = opposite_atom;
            bParitiesInverted  = (nNeighMode == NEIGH_MODE_RING &&
                                  IS_ALLENE_CHAIN(at[i1].stereo_bond_parity[n1]) &&
                                  PARITY_CALCULATE(at[i1].stereo_bond_parity[n1]) &&
                                  at[i1].nRingSystem == at[opposite_atom].nRingSystem &&
                                  at[opposite_atom].valence==MAX_NUM_STEREO_BONDS)? -1 : 0;
            ret1 = ret2 = 0;
            if ( 0 < (ret1=CreateCheckSymmPaths( at, (AT_RANK)i1, neigh[0], (AT_RANK)i1, neigh[1], nAvoidCheckAtom,
                                       nVisited1, nVisited2, nAtomNumberCanon1, nAtomNumberCanon2,
                                       nl01, nl02, nSymmRank1, nSymmRank2, nCanonRank, &nLength, &bParitiesInverted, 0 ) ) &&
                 0 < (ret2=CalculatedPathsParitiesAreIdentical( at, num_atoms, nSymmRank,
                                       nCanonRank, nAtomNumberCanon, nAtomNumberCanon1, nAtomNumberCanon2,
                                       nVisited1, nVisited2, opposite_atom, (AT_RANK)i1,
                                       neigh[0], neigh[1], nNeighMode, bParitiesInverted, 0, 
                                       pCS, vABParityUnknown ) ) ) {
                if ( ret2 & ( NOT_WELL_DEF_UNKN | NOT_WELL_DEF_UNDF ) ) {
                    /*  possibly change the parity to unknown or undefined */
                    int new_parity = (ret2 & NOT_WELL_DEF_UNKN)? vABParityUnknown /*AB_PARITY_UNKN*/: AB_PARITY_UNDF;
                    if ( (PARITY_ILL_DEF(at[i1].stereo_bond_parity[n1]) && PARITY_VAL(at[i1].stereo_bond_parity[n1]) > new_parity) ||
                         PARITY_CALCULATE(at[i1].stereo_bond_parity[n1]) ) {
                        /*  set new unknown or undefined parity */
                        SetOneStereoBondIllDefParity( at, i1, /* atom number*/ n1 /* stereo bond ord. number*/, new_parity );
                        /*  change in pCS */
                        nAtomRank1 = inchi_max( nCanonRank[i1], nCanonRank[opposite_atom]);
                        nAtomRank2 = inchi_min( nCanonRank[i1], nCanonRank[opposite_atom]);
                        for ( n = 0, m = pCS->nLenLinearCTStereoDble-1; n <= m; n ++ ) {
                            if ( pCS->LinearCTStereoDble[n].at_num1 == nAtomRank1 &&
                                 pCS->LinearCTStereoDble[n].at_num2 == nAtomRank2 ) {
                                pCS->LinearCTStereoDble[n].parity = new_parity;
#if ( bRELEASE_VERSION == 0 )
                                pCS->bExtract |= EXTR_CALC_USED_TO_REMOVE_PARITY;
#endif
                                m = -1;
                                break;
                            }
                        }
                        if ( m >= 0 ) {
                            ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                            goto exit_function;
                        }
                        ret ++;
                    }
                } else {
                    /*  remove the parity */
                    if ( !RemoveOneStereoBond( at, i1, /* atom number*/ n1 /* stereo bond ord. number*/ ) ) {
                        ret = CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                        goto exit_function;
                    }
                    n1 --;  /*  cycle counter may temporarily become negative */
                    /*  Remove from the pCS */
                    nAtomRank1 = inchi_max( nCanonRank[i1], nCanonRank[opposite_atom]);
                    nAtomRank2 = inchi_min( nCanonRank[i1], nCanonRank[opposite_atom]);
                    for ( n = 0, m = pCS->nLenLinearCTStereoDble-1; n <= m; n ++ ) {
                        if ( pCS->LinearCTStereoDble[n].at_num1 == nAtomRank1 &&
                             pCS->LinearCTStereoDble[n].at_num2 == nAtomRank2 ) {
                            if ( n < m ) { /*  remove pCS->LinearCTStereoDble[n] */
                                memmove( pCS->LinearCTStereoDble + n,
                                         pCS->LinearCTStereoDble + n + 1,
                                         (m-n)*sizeof(pCS->LinearCTStereoDble[0]) );
                            }
                            pCS->nLenLinearCTStereoDble --;
#if ( bRELEASE_VERSION == 0 )
                            pCS->bExtract |= EXTR_CALC_USED_TO_REMOVE_PARITY;
#endif
                            m = -1;
                            break;
                        }
                    }
                    if ( m >= 0 ) {
                        ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                        goto exit_function;
                    }
                    ret ++;
                }
            } else {
                if ( !ret_failed ) {
                    ret_failed = (ret1<0)? ret1 : (ret2<0)? ret2 : 0;
                }
                if ( !RETURNED_ERROR(ret_failed) ) {
                    if ( RETURNED_ERROR( ret1 ) )
                        ret_failed = ret1;
                    else
                    if ( RETURNED_ERROR( ret2 ) )
                        ret_failed = ret2;
                }
            }
        }
    }
    if ( nNeighMode == NEIGH_MODE_CHAIN && nNumEqRingNeigh && !RETURNED_ERROR(ret_failed) ) {
        nNeighMode = NEIGH_MODE_RING;
        goto second_pass;
    }

exit_function:

    return RETURNED_ERROR(ret_failed)? ret_failed : ret_failed? -(ret_failed+1) : ret;
}
/****************************************************************************/
/*  Remove stereo marks from the atoms that are calculated to be non-stereo */
/*  (can find two numberings that change only one stereo center parity)     */
int RemoveCalculatedNonStereoCenterParities( sp_ATOM *at, int num_atoms, int num_at_tg,
                                          AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList,
                                          AT_RANK *nCanonRank, const AT_RANK *nSymmRank,
                                          AT_RANK *nAtomNumberCanon, AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2,
                                          NEIGH_LIST *nl, NEIGH_LIST *nl1, NEIGH_LIST *nl2, 
                                          AT_RANK *nVisited1, AT_RANK *nVisited2, 
                                          CANON_STAT *pCS,
                                          int vABParityUnknown)
{
    int j, n, m, ret;
    
    int i, k, ret1, ret2, ret_failed=0, mode, max_mode;
    AT_RANK nAtomRank1, neigh[MAX_NUM_STEREO_ATOM_NEIGH], nAvoidCheckAtom[2], nLength;
    int         nNeighMode = NEIGH_MODE_CHAIN;
    int         nNumEqRingNeigh = 0, bRingNeigh, bSymmNeigh, bParitiesInverted;
    NEIGH_LIST *nl01, *nl02;
    const AT_RANK    *nSymmRank1, *nSymmRank2;
    
    ret = 0;

second_pass:
    for ( i = 0; i < num_atoms && !RETURNED_ERROR(ret_failed); i ++ ) {
        if ( !at[i].parity || at[i].stereo_bond_neighbor[0] ) {
            continue;
        }
        if ( at[i].valence > MAX_NUM_STEREO_ATOM_NEIGH ) {
            continue; /*  error: stereo center cannot have more than 4 neighbors */ /*   <BRKPT> */
        }
        /*  at[i1] is a stereo center */
        if ( !PARITY_CALCULATE(at[i].stereo_atom_parity) && !PARITY_ILL_DEF(at[i].stereo_atom_parity) ) {
            continue;
        }
        /* neighbors sorted according to symm. ranks (primary key) and canon. ranks (secondary key), in descending order */
        /* sorting guarantees that for two constit. equ. neighbors canon. ranks of the first is greater */
        /* !!! previously (but not anymore) the canon. rank of neigh[0] was greater than the others !!! */
        for ( j = 0; j < at[i].valence; j ++ ) {
            neigh[j] = nl[i][j+1]; /*  sorting does NOT guarantee that canon. rank of neigh[0] is greater than others */
        }
        /* 
         *  mode = 0 => Standard approach: switch 2 neighbors
         *         1 => Check for C2v reflection leading to parity inversion
         *         2 => Check for C2 rotation preserving parities
         *         3 => Check for S4 rotation/reflection leading to parity inversion
         */
#if ( CHECK_C2v_S4_SYMM == 1 )
        if ( nNeighMode = NEIGH_MODE_RING && at[i].valence == 4 &&
             nSymmRank[(int)neigh[0]] == nSymmRank[(int)neigh[1]] &&
             nSymmRank[(int)neigh[2]] == nSymmRank[(int)neigh[3]] &&
             !at[i].bCutVertex 
           ) {
            if ( nSymmRank[(int)neigh[1]] == nSymmRank[(int)neigh[2]] ) {
                max_mode = MAP_MODE_S4;
            } else {
                max_mode = inchi_max(MAP_MODE_C2v, MAP_MODE_C2);
            }
        } else {
            max_mode = MAP_MODE_STD;
        }
#else
        max_mode = MAP_MODE_STD;
#endif
        for ( j = 0; j < at[i].valence && at[i].parity && !RETURNED_ERROR(ret_failed); j ++ ) {
            for ( k = j+1; k < at[i].valence && at[i].parity && !RETURNED_ERROR(ret_failed); k ++ ) {
                for ( mode = 0; mode <= max_mode && at[i].parity && !RETURNED_ERROR(ret_failed); mode ++ ) {
                    if ( nSymmRank[(int)neigh[j]] != nSymmRank[(int)neigh[k]] ) {
                        continue; /*  the two neighbors are not constitutionally identical */
                    }
                    bRingNeigh = (at[(int)neigh[j]].nRingSystem == at[(int)neigh[k]].nRingSystem);
                    switch ( nNeighMode ) {
                    case NEIGH_MODE_CHAIN:
                        if ( bRingNeigh ) {
                            nNumEqRingNeigh ++;
                            continue;
                        }
                        nl01 = nl;
                        nl02 = nl;
                        nSymmRank1 = nSymmRank;
                        nSymmRank2 = nSymmRank;
                        break;
                    case NEIGH_MODE_RING:
                        if ( !bRingNeigh )
                            continue;
                        /*  break a tie between the two contitutionally equivalent neighbors, */
                        /*  refine the two partitions, sort neighbors lists nl1, nl2 */
                        bSymmNeigh = BreakNeighborsTie(  at, num_atoms, num_at_tg, MAX_ATOMS+1, i,
                                            neigh, j, k, mode,
                                            pRankStack1, pRankStack2, nTempRank, NeighList, nSymmRank, nCanonRank,
                                            nl1, nl2, &pCS->lNumNeighListIter );
                        if ( bSymmNeigh <= 0 ) {
                            if ( ret_failed > bSymmNeigh )
                                ret_failed = bSymmNeigh;
                            continue;
                        }
                        nl01 = nl1;
                        nl02 = nl2;
                        nSymmRank1 = pRankStack1[0];
                        nSymmRank2 = pRankStack2[0];
                        break;
                    default:
                        return CT_STEREOCOUNT_ERR;  /*  <BRKPT> */
                    }

                    /*  initialize arrays */
                    memset( nVisited1, 0, sizeof(nVisited1[0])*num_atoms );
                    memset( nVisited2, 0, sizeof(nVisited2[0])*num_atoms );
                    memset( nAtomNumberCanon1, 0, sizeof(nAtomNumberCanon1[0])*num_atoms );
                    memset( nAtomNumberCanon2, 0, sizeof(nAtomNumberCanon2[0])*num_atoms );
                    nLength = 1;
                    nVisited1[i] = i+1;   /*  start atom is same */
                    nVisited2[i] = i+1;
                    nAtomNumberCanon1[i] = nLength;
                    nAtomNumberCanon2[i] = nLength;
                    nAvoidCheckAtom[0] = i;
                    nAvoidCheckAtom[1] = MAX_ATOMS+1;
                
                    bParitiesInverted  = (mode==MAP_MODE_C2v || mode==MAP_MODE_S4)? -1 : 0;
                    /*
                    if (nNeighMode==NEIGH_MODE_RING && at[i].valence==MAX_NUM_STEREO_ATOM_NEIGH) {
                        AT_RANK other_neigh[2];
                        int     n;
                        for ( m = n = 0; m < MAX_NUM_STEREO_ATOM_NEIGH; m ++ ) {
                            if ( at[i].neighbor[m] != neigh[j] && at[i].neighbor[m] != neigh[k] )
                                other_neigh[n++] = at[i].neighbor[m];
                        }
                        if ( nSymmRank[(int)other_neigh[0]] == nSymmRank[(int)other_neigh[1]] )
                            bParitiesInverted = -1;
                    }
                    */
                    /* allow matching inverted centers only in case all equivalent neighbors in same ring system */

                    ret2 = 0; /* initilize. 1/8/2002 */
                
                    if ( 0 < (ret1 = CreateCheckSymmPaths( at, (AT_RANK)i, neigh[j], (AT_RANK)i, neigh[k],
                                               nAvoidCheckAtom,
                                               nVisited1, nVisited2, nAtomNumberCanon1, nAtomNumberCanon2,
                                               nl01, nl02, nSymmRank1, nSymmRank2, nCanonRank, &nLength,
                                               &bParitiesInverted, mode ) ) &&
                         0 < (ret2 = CalculatedPathsParitiesAreIdentical( at, num_atoms, nSymmRank,
                                               nCanonRank, nAtomNumberCanon, nAtomNumberCanon1, nAtomNumberCanon2,
                                               nVisited1, nVisited2, (AT_RANK)MAX_ATOMS, (AT_RANK)i,
                                               neigh[j], neigh[k], nNeighMode, 
                                               bParitiesInverted, mode, pCS,
                                               vABParityUnknown) ) ) {
                        if ( ret2 & ( NOT_WELL_DEF_UNKN | NOT_WELL_DEF_UNDF ) ) {
                            /*  possibly change the parity to unknown or undefined */
                            int new_parity = (ret2 & NOT_WELL_DEF_UNKN)? vABParityUnknown /*AB_PARITY_UNKN*/: AB_PARITY_UNDF;
                            if ( (PARITY_ILL_DEF(at[i].stereo_atom_parity) &&
                                 PARITY_VAL(at[i].stereo_atom_parity) > new_parity) ||
                                 PARITY_CALCULATE(at[i].stereo_atom_parity) ) {
                                /*  set new unknown or undefined parity */
                                at[i].stereo_atom_parity = (at[i].stereo_atom_parity ^ PARITY_VAL(at[i].stereo_atom_parity)) | PARITY_VAL(new_parity);
                                at[i].parity = PARITY_VAL(new_parity);
                                /*  Remove from pCS */
                                nAtomRank1 = nCanonRank[i];
                                for ( n = 0, m = pCS->nLenLinearCTStereoCarb-1; n <= m; n ++ ) {
                                    if ( pCS->LinearCTStereoCarb[n].at_num == nAtomRank1 ) {
                                        pCS->LinearCTStereoCarb[n].parity = PARITY_VAL(new_parity);
    #if ( bRELEASE_VERSION == 0 )
                                        pCS->bExtract |= EXTR_CALC_USED_TO_REMOVE_PARITY;
    #endif
                                        m = -1;
                                        break;
                                    }
                                }
                                if ( m >= 0 ) {
                                    ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                                    goto exit_function;
                                }
                                ret ++; /*  number of removed or set unknown/undefined parities */
                            }
                        } else {
                            RemoveOneStereoCenter( at, i /* atom number*/ );
                            /*  Remove from pCS */
                            nAtomRank1 = nCanonRank[i];
                            for ( n = 0, m = pCS->nLenLinearCTStereoCarb-1; n <= m; n ++ ) {
                                if ( pCS->LinearCTStereoCarb[n].at_num == nAtomRank1 ) {
                                    if ( n < m ) { /*  remove pCS->LinearCTStereoDble[n] */
                                        memmove( pCS->LinearCTStereoCarb + n,
                                                 pCS->LinearCTStereoCarb + n + 1,
                                                 (m-n)*sizeof(pCS->LinearCTStereoCarb[0]) );
                                    }
                                    pCS->nLenLinearCTStereoCarb --;
    #if ( bRELEASE_VERSION == 0 )
                                    pCS->bExtract |= EXTR_CALC_USED_TO_REMOVE_PARITY;
    #endif
                                    m = -1;
                                    break;
                                }
                            }
                            if ( m >= 0 ) {
                                ret = CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                                goto exit_function;
                            }
                            ret ++;  /*  number of removed or set unknown/undefined parities */
                        }
                    } else {
                        if ( !ret_failed ) {
                            if ( ret1 < 0 ) {
                                ret_failed = ret1;
                            } else
                            if ( ret2 < 0 ) {
                                ret_failed = ret2;
                            }
                        }
                        if ( !RETURNED_ERROR(ret_failed) ) {
                            if ( RETURNED_ERROR( ret1 ) )
                                ret_failed = ret1;
                            else
                            if ( RETURNED_ERROR( ret2 ) )
                                ret_failed = ret2;
                        }
                    }
                }
            }
        }
    }
    if ( nNeighMode == NEIGH_MODE_CHAIN && nNumEqRingNeigh && !RETURNED_ERROR(ret_failed) ) {
        nNeighMode = NEIGH_MODE_RING;
        goto second_pass;
    }

exit_function:

    return RETURNED_ERROR(ret_failed)? ret_failed : ret_failed? -(ret+1) : ret;
}

/**************************************************************************************/
int RemoveCalculatedNonStereo( sp_ATOM *at, int num_atoms, int num_at_tg,
                              AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList,
                              const AT_RANK *nSymmRank, AT_RANK *nCanonRank, 
                              AT_RANK *nAtomNumberCanon, CANON_STAT *pCS,
                              int vABParityUnknown)
{
    NEIGH_LIST *nl = NULL, *nl1 = NULL, *nl2 = NULL;
    AT_RANK    *nVisited1 = NULL, *nVisited2 = NULL, *nAtomNumberCanon1 = NULL, *nAtomNumberCanon2 = NULL;
    int        nNumRemoved = 0, nTotRemoved = 0, ret = 0, ret1 = 0, ret2 = 0;
    
    if ( !AllocateForNonStereoRemoval( at, num_atoms, nSymmRank, nCanonRank,
                                       &nAtomNumberCanon1, &nAtomNumberCanon2,
                                       &nl, &nl1, &nl2, &nVisited1, &nVisited2 ) ) {
        return CT_OUT_OF_RAM;  /*   <BRKPT> */
    }
    
    do {
        nNumRemoved = 0;
        /*  bonds */
        ret = RemoveCalculatedNonStereoBondParities( at, num_atoms, num_at_tg,
                                              pRankStack1, pRankStack2, nTempRank, NeighList,
                                              nCanonRank, nSymmRank,
                                              nAtomNumberCanon, nAtomNumberCanon1, nAtomNumberCanon2,
                                              nl, nl1, nl2, nVisited1, nVisited2, pCS,
                                              vABParityUnknown);
        if ( RETURNED_ERROR( ret ) ) {
            goto exit_function;
        }
        if ( ret < 0  ) {
            if ( ret < ret1 ) {  /*   <BRKPT> */
                ret1 = ret;          
            }
            ret = - ( ret + 1 ); /*  number of removed */
        }
        nNumRemoved += ret;

        /*  centers */
        ret = RemoveCalculatedNonStereoCenterParities( at, num_atoms, num_at_tg,
                                              pRankStack1, pRankStack2, nTempRank, NeighList,
                                              nCanonRank, nSymmRank,
                                              nAtomNumberCanon, nAtomNumberCanon1, nAtomNumberCanon2,
                                              nl, nl1, nl2, nVisited1, nVisited2, pCS,
                                              vABParityUnknown);
        if ( RETURNED_ERROR( ret ) ) {
            goto exit_function;
        }
        if ( ret < 0  ) {
            if ( ret < ret2 ) {  /*   <BRKPT> */
                ret2 = ret;          
            }
            ret = - ( ret + 1 ); /*  number of removed */
        }
        nNumRemoved += ret;

        nTotRemoved += nNumRemoved;

    } while ( nNumRemoved );

    if ( !RETURNED_ERROR( ret1 ) && !RETURNED_ERROR( ret2 ) ) {
        ret = inchi_min( ret1, ret2 );
        ret = (ret >= 0)? nTotRemoved : -(1+nTotRemoved);
    }

exit_function:
    
    DeAllocateForNonStereoRemoval( &nAtomNumberCanon1, &nAtomNumberCanon2, &nl, &nl1, &nl2, &nVisited1, &nVisited2 );
    
    return ret;
}
#endif /* } REMOVE_CALC_NONSTEREO */
