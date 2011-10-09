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

#define RET_MAX 32767

/**********************************************************************************/
void inchi_swap ( char *a, char *b, size_t width )
{
    char tmp;
    if ( a != b )
        while ( width-- ) {
            tmp = *a;
            *a++ = *b;
            *b++ = tmp;
        }
}
/**********************************************************************************/
/*  Sort by insertions */
int insertions_sort( void *base, size_t num, size_t width, int ( *compare )(const void *e1, const void *e2 ) )
{
  char *i, *j, *pk = (char*)base;
  int  num_trans = 0;
  size_t k;
  for( k=1; k < num; k++, pk += width ) {
     /*for( i = pk, j = pk + width; j > (char*)base && (*compare)(i,j) > 0; j=i, i -= width )*/
     for( i = j = pk + width; j > (char*)base && (i -= width,(*compare)(i,j)) > 0; j=i ) /* changed to keep BoundsChecker happy 2007-09-24 DT */
     {
        inchi_swap( i, j, width );
        num_trans ++;
     }
  }
  return num_trans;
}
/**********************************************************************************/
/*  Sort by insertions */
int insertions_sort_AT_NUMBERS( AT_NUMB *base, int num, int ( *compare )(const void *e1, const void *e2 ) )
{
  AT_NUMB *i, *j, *pk, tmp;
  int  k, num_trans = 0;
  for( k=1, pk = base; k < num; k++, pk ++ ) {
     for( j = (i = pk) + 1, tmp = *j; j > base && (*compare)(i,&tmp) > 0; j=i, i -- ) {
        *j = *i;
        num_trans ++;
     }
     *j = tmp;
  }
  return num_trans;
}
/**********************************************************************************/
/*  Sort neighbors according to ranks in ascending order */
void insertions_sort_NeighList_AT_NUMBERS( NEIGH_LIST base, AT_RANK *nRank )
{
  AT_NUMB *i, *j, *pk, tmp;
  AT_RANK rj; /* optimization */
  int k, num = (int)*base++;
  for( k=1, pk = base; k < num; k++, pk ++ ) {
     for( j = (i = pk) + 1, rj=nRank[(int)*j]; j > base && nRank[(int)*i] > rj; j=i, i -- ) {
         tmp = *i;
         *i = *j;
         *j = tmp;
     }
  }
}
/**********************************************************************************/
/*  Sort neighbors according to ranks in ascending order */
int insertions_sort_AT_RANK( AT_RANK *base, int num )
{
  AT_RANK *i, *j, *pk, tmp;
  int  k, num_trans = 0;
  for( k=1, pk = base; k < num; k++, pk ++ ) {
     for( j = (i = pk) + 1, tmp = *j; j > base && *i > tmp; j=i, i -- ) {
        *j = *i;
        num_trans ++;
     }
     *j = tmp;
  }
  return num_trans;
}
/**********************************************************************************/
/*  Sort neighbors according to ranks in ascending order */
int insertions_sort_NeighList_AT_NUMBERS3( NEIGH_LIST base, AT_RANK *nRank )
{
  AT_NUMB *i, *j, *pk, tmp;
  AT_RANK rj;
  int k, n, num = (int)*base++;
  for( k=1, pk = base, n=0; k < num; k++, pk ++ ) {
     for( j = (i = pk) + 1, rj=nRank[(int)(tmp=*j)]; j > base && nRank[(int)*i] > rj; j=i, i -- ) {
         *j = *i;
         n ++;
     }
     *j = tmp;
  }
  return n;
}

/**********************************************************************************/
/*  Sort neighbors according to symm. ranks (primary key) and canon. ranks (secondary key), in descending order */
void insertions_sort_NeighListBySymmAndCanonRank( NEIGH_LIST base, const AT_RANK *nSymmRank, const AT_RANK *nCanonRank )
{
  AT_NUMB *i, *j, *pk, tmp;
  int  diff;
  int k, num = (int)*base++;
  for( k=1, pk = base; k < num; k++, pk ++ ) {
     for( j = (i = pk) + 1; j > base &&  /*  always j > i */
          ( 0 > (diff = (int)nSymmRank[(int)*i] - (int)nSymmRank[(int)*j]) ||
            !diff && nCanonRank[(int)*i] < nCanonRank[(int)*j]); j=i, i -- ) {
         tmp = *i;
         *i = *j;
         *j = tmp;
     }
  }
}

/*********************************************************************************************
 *
 *  Comparison functions
 *
 *********************************************************************************************/
int CompNeighborsAT_NUMBER( const void* a1, const void* a2)
{
#ifdef CT_NEIGH_INCREASE
    return (int)pn_RankForSort[pNeighborsForSort[(int)*(const AT_NUMB*)a1]] -
           (int)pn_RankForSort[pNeighborsForSort[(int)*(const AT_NUMB*)a2]];
#else
    return (int)pn_RankForSort[pNeighborsForSort[(int)*(const AT_NUMB*)a2]] -
           (int)pn_RankForSort[pNeighborsForSort[(int)*(const AT_NUMB*)a1]];
#endif
}

/**********************************************************************************/
int comp_AT_RANK( const void* a1, const void* a2)
{
    return (int)*(const AT_RANK*)a1 - (int)*(const AT_RANK*)a2;
}

/**********************************************************************************/
/*  Compare for sorting Ranks only */
int CompRank(const void* a1, const void* a2 )
{
    int ret = (int)pn_RankForSort[(int)*(const AT_RANK*)a1] -
              (int)pn_RankForSort[(int)*(const AT_RANK*)a2];
    return ret;
}
/**********************************************************************************/
int CompRanksOrd( const void* a1, const void* a2 )
{
    int ret;
    ret = (int)pn_RankForSort[(int)*(const AT_RANK*)a1] -
          (int)pn_RankForSort[(int)*(const AT_RANK*)a2];
    if ( !ret )
        ret = (int)*(const AT_RANK*)a1 - (int)*(const AT_RANK*)a2;
    return ret;
}
/**********************************************************************************/
int CompAtomInvariants2Only( const void* a1, const void* a2 )
{
    const ATOM_INVARIANT2 *pAI1 = pAtomInvariant2ForSort + (int)*(const AT_RANK*)a1;
    const ATOM_INVARIANT2 *pAI2 = pAtomInvariant2ForSort + (int)*(const AT_RANK*)a2;
    int i;
    for ( i = 0; i < AT_INV_BREAK1; i ++ ) {
        if ( pAI1->val[i] == pAI2->val[i] )
            continue;
        return  (int)pAI1->val[i] - (int)pAI2->val[i];
    }
    if ( pAI1->iso_sort_key != pAI2->iso_sort_key ) {
        return ( pAI1->iso_sort_key > pAI2->iso_sort_key )? 1 : -1;
    }
    for ( ; i < AT_INV_LENGTH; i ++ ) {
        if ( pAI1->val[i] != pAI2->val[i] )
            continue;
        return  (int)pAI1->val[i] - (int)pAI2->val[i];
    }
    if ( pAI1->iso_aux_key != pAI2->iso_aux_key ) {
        return ( pAI1->iso_aux_key > pAI2->iso_aux_key )? 1 : -1;
    }
    return 0;
}
/**********************************************************************************/
int CompAtomInvariants2( const void* a1, const void* a2 )
{
    /*  Warning: the following line may be compiler implementation dependent */
    int ret = CompAtomInvariants2Only( a1, a2 );
    if ( !ret )
        ret = (int)*(const AT_RANK*)a1 - (int)*(const AT_RANK*)a2;
    return ret;
}
/**********************************************************************************/
/*  Compare two elements lexicographically */
int CompChemElemLex( const void *a1, const void *a2 )
{
    return memcmp( a1, a2, 2);
}
/**********************************************************************************/
/*  lexicographic compare */
int CompareNeighListLex( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank)
{
    int len1 = (int)*pp1++;
    int len2 = (int)*pp2++;
    int len  = inchi_min( len1, len2 );
    int diff = 0;
    while ( len -- > 0 && !( diff = (int)nRank[*pp1++] - (int)nRank[*pp2++] ) )
        ;
    return diff? diff : (len1 - len2);

}
/**********************************************************************************/
/*  lexicographic compare */
int CompareNeighListLexUpToMaxRank( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank, AT_RANK nMaxAtNeighRank )
{
    int len1 = (int)*pp1++;
    int len2 = (int)*pp2++;
    int diff = 0;
    int len;
    while( 0 < len1 && nRank[pp1[len1-1]] > nMaxAtNeighRank ) {
        len1 --;
    }
    while( 0 < len2 && nRank[pp2[len2-1]] > nMaxAtNeighRank ) {
        len2 --;
    }
    len  = inchi_min( len1, len2 );
    while ( len -- > 0 && !( diff = (int)nRank[*pp1++] - (int)nRank[*pp2++] ) )
        ;
    return diff? diff : (len1 - len2);

}
/**********************************************************************************/
int compare_NeighLists( const NEIGH_LIST *op1,  const NEIGH_LIST *op2 )
{
    return CompareNeighListLex( *op1, *op2, pn_RankForSort);
}
/**********************************************************************************/
int CompNeighListRanks( const void* a1, const void* a2 )
{
    int ret;
    ret = (int)pn_RankForSort[*((const AT_RANK*)a1)] -
          (int)pn_RankForSort[*((const AT_RANK*)a2)];
    if ( !ret )
        ret = compare_NeighLists( pNeighList_RankForSort + *((const AT_RANK*)a1),
                                  pNeighList_RankForSort + *((const AT_RANK*)a2) );
    return ret;
}
/**********************************************************************************/
int CompNeighLists( const void* a1, const void* a2 )
{
    int ret;
        ret = compare_NeighLists( pNeighList_RankForSort + *((const AT_RANK*)a1),
                                  pNeighList_RankForSort + *((const AT_RANK*)a2) );
    return ret;
}
/**********************************************************************************/
int CompNeighListsUpToMaxRank( const void* a1, const void* a2 )
{
    int ret;
        ret = CompareNeighListLexUpToMaxRank( pNeighList_RankForSort[*((const AT_RANK*)a1)],
                                              pNeighList_RankForSort[*((const AT_RANK*)a2)],
                                              pn_RankForSort, nMaxAtNeighRankForSort );
    return ret;
}
/**********************************************************************************/
int CompNeighListRanksOrd( const void* a1, const void* a2 )
{
    int ret = CompNeighListRanks( a1, a2 );
    if ( !ret )
        ret = (int)*((const AT_RANK*)a1) - (int)*((const AT_RANK*)a2); /*  keep original order if identical */
    return ret;
}
/**********************************************************************************/
int CompRanksInvOrd( const void* a1, const void* a2 )
{
    return (int)*(const AT_RANK*)a2 - (int)*(const AT_RANK*)a1;
}
/**********************************************************************************/
int CompNeighborsRanksCountEql( const void* a1, const void* a2 )
{
#ifdef CT_NEIGH_INCREASE
    int ret = (int)pn_RankForSort[(int)*(const AT_RANK*)a1] -
              (int)pn_RankForSort[(int)*(const AT_RANK*)a2];
#else
    int ret = (int)pn_RankForSort[(int)*(const AT_RANK*)a2] -
              (int)pn_RankForSort[(int)*(const AT_RANK*)a1];
#endif
    nNumCompNeighborsRanksCountEql += !ret;
    return ret;
}
/****************************************************************************************
 *
 * In this neighbor list the (vertex number) = (canonical number) - 1
 * Since LinearCT is sorted so that parents are in ascending order
 * and all neighbors of a parent are smaller than the parent and are
 * in ascending order, the neighbors in the NEIGH_LIST are automatically
 * sorted in ascending order
 */
NEIGH_LIST *CreateNeighListFromLinearCT( AT_NUMB *LinearCT, int nLenCT, int num_atoms )
{
    /* atom numbers in LinearCT are canonical numbers
     * order: parent[i] > neigh[i][0] < neigh[i][1]...<neigh[i][n] < parent[i+1] > neigh[i+1][0] < ...
     *        parent[i] < parent[i+1]
     */
    int i, j;
    S_CHAR     *valence  = NULL;
    NEIGH_LIST *pp = NULL;
    AT_NUMB    *pAtList = NULL;
    AT_RANK     n_vertex, n_neigh;
    int err = 1, num_bonds;
    int length, start;
    if ( (int)LinearCT[0] > num_atoms ) {
        goto exit_function;
    }
    if ( !(valence = (S_CHAR*)inchi_calloc( num_atoms+1, sizeof(valence[0]) ) ) ) {
        goto exit_function;
    }
    for ( i = 1, num_bonds = 0, n_vertex = LinearCT[0]; i < nLenCT; i ++ ) {
        if ( (n_neigh = LinearCT[i]) < n_vertex ) {
            valence[n_neigh] ++;
            valence[n_vertex] ++;
            num_bonds += 2;
        } else
        if ( (int)(n_vertex = n_neigh) > num_atoms ) {
            goto exit_function;
        }
    }
    if ( (int)n_vertex != num_atoms ) {
        goto exit_function;
    }
    length = num_bonds + num_atoms + 1;
    if ( pp = (NEIGH_LIST *) inchi_calloc((num_atoms+1), sizeof(NEIGH_LIST)) ) {
        if ( pAtList = (AT_NUMB *) inchi_malloc( length*sizeof(*pAtList) ) ) {
            /*  create empty connection table */
            for ( i = 1, length = 0; i <= num_atoms; i ++ ) {
                start = length;
                length += (valence[i]+1);
                pp[i-1] = pAtList + start;
                pp[i-1][0] = 0;
            }
            /*  fill out the CT */
            for ( i = 1, n_vertex = LinearCT[0]-1; i < nLenCT; i ++ ) {
                if ( (n_neigh = LinearCT[i]-1) < n_vertex ) {
                    /*  vertex - neighbor connection */
                    j = (int)(++pp[(int)n_vertex][0]);
                    pp[(int)n_vertex][j] = n_neigh;
                    /*  neighbor - vertex connection */
                    j = (int)(++pp[(int)n_neigh][0]);
                    pp[(int)n_neigh][j] = n_vertex;

                } else
                if ( (int)(n_vertex = n_neigh) >= num_atoms ) {
                    goto exit_function;
                }
            }
            err = 0;
        }
    }
exit_function:
    if ( valence ) {
        inchi_free( valence );
    }
    if ( err ) {
        if ( pAtList )
            inchi_free( pAtList );
        if ( pp ) {
            inchi_free( pp );
            pp = NULL;
        }
    }
    return pp;
}

/***********************************************************************************
 * NEIGH_LIST pp[] is an array of pointers to the lists of neighboring atoms numbers
 *  The first number in each list is a number of neighbors.
 *  In case of bDoubleBondSquare != 0 neighbors connected by the double bond appear 2 times
 * The first element pp[0] is a pointer to be deallocated to free all the lists.
 */
NEIGH_LIST *CreateNeighList( int num_atoms, int num_at_tg, sp_ATOM* at,
                             int bDoubleBondSquare, T_GROUP_INFO *t_group_info )
{
    /*  +1 to add NULL termination */
    NEIGH_LIST *pp = (NEIGH_LIST *) inchi_calloc((num_at_tg+1), sizeof(NEIGH_LIST));
    T_GROUP   *t_group             = NULL;
    AT_NUMB   *nEndpointAtomNumber = NULL;
    int        num_t_groups        = 0;
    int        nFirstEndpointAtNoPos;

    AT_NUMB   *pAtList = NULL;
    int        length, start, val, i, j;
    if ( pp ) {
        if ( num_at_tg > num_atoms ) {
            t_group             = t_group_info->t_group;
            num_t_groups        = t_group_info->num_t_groups;
            nEndpointAtomNumber = t_group_info->nEndpointAtomNumber;
        }

        if ( !bDoubleBondSquare ) {
            for ( i = 0, length = 0; i < num_atoms; i ++ ) {
                length += (int)at[i].valence + (num_t_groups && at[i].endpoint);
            }
            length += num_atoms;
            for ( i = 0; i < num_t_groups; i ++ ) {
                length += (int)t_group[i].nNumEndpoints;
            }
            length += num_t_groups;

        } else {
            for ( i = 0, length = 0; i < num_atoms; i ++ ) {
                val = (int)at[i].valence;
                for ( j = 0; j < val; j ++ ) {
                    length += 1 + (bDoubleBondSquare && BOND_DOUBLE == at[i].bond_type[j]);
                }
                length += (num_t_groups && at[i].endpoint);
            }
            length += num_atoms;
            for ( i = 0; i < num_t_groups; i ++ ) {
                length += (int)t_group[i].nNumEndpoints;
            }
            length += num_t_groups;
        }
        length ++; /*  +1 to save number of neighbors */
        if ( pAtList = (AT_NUMB *) inchi_malloc( length*sizeof(*pAtList) ) ) {
            if ( !bDoubleBondSquare ) {
                for ( i = 0, length = 0; i < num_atoms; i ++ ) {
                    val = at[i].valence;
                    start = length ++;
                    for ( j = 0; j < val; j ++ ) {
                        pAtList[length ++] = at[i].neighbor[j];
                    }
                    /*  add endpoint */
                    if (num_t_groups && at[i].endpoint) {
                        pAtList[length ++] = num_atoms + (int)at[i].endpoint - 1;
                    }
                    pAtList[start] = length - start - 1;  /*  number of neighbors before the list of neighbors */
                    pp[i] = pAtList + start;              /*  pointer to the <num.neigh.><list of neigh> */
                }

            } else {
                for ( i = 0, length = 0; i < num_atoms; i ++ ) {
                    val = at[i].valence;
                    start = length ++;
                    for ( j = 0; j < val; j ++ ) {
                        pAtList[length ++] = at[i].neighbor[j];
                        if ( bDoubleBondSquare && BOND_DOUBLE == at[i].bond_type[j] ) {
                            pAtList[length ++] = at[i].neighbor[j]; /*  a list of neighbor orig. numbers */
                        }
                    }
                    /*  add endpoint */
                    if (num_t_groups && at[i].endpoint) {
                        pAtList[length ++] = num_atoms + (int)at[i].endpoint - 1;
                    }
                    pAtList[start] = length - start - 1;  /*  number of neighbors before the list of neighbors */
                    pp[i] = pAtList + start;              /*  pointer to the <num.neigh.><list of neigh> */
                }
            }

            /*  add t-groups */
            for ( i = 0; i < num_t_groups; i ++ ) {
                val = (int)t_group[i].nNumEndpoints;
                start = length ++;
                nFirstEndpointAtNoPos = (int)t_group[i].nFirstEndpointAtNoPos;
                for ( j = 0; j < val; j ++ ) {
                    pAtList[length ++] = nEndpointAtomNumber[nFirstEndpointAtNoPos+j];
                }
                pAtList[start] = length - start - 1;  /*  number of neighbors before the list of neighbors */
                pp[num_atoms+i] = pAtList + start;    /*  pointer to the <num.neigh.><list of neigh> */
            }
        } else {
            inchi_free ( pp );
            return NULL;
        }
    }
    return pp;
}
/**********************************************************************************/
void FreeNeighList( NEIGH_LIST *pp )
{
    if ( pp ) {
        if ( pp[0] ) {
            inchi_free( pp[0] );
        }
        inchi_free( pp );
    }
}

/**********************************************************************************/
int BreakAllTies( int num_atoms, int num_max, AT_RANK **pRankStack,
                     NEIGH_LIST *NeighList, AT_RANK *nTempRank, CANON_STAT *pCS)
{
    int i, nRet = -1, nNumRanks=1 /* value does not matter*/;
    
    AT_RANK *nPrevRank       = *pRankStack ++;
    AT_RANK *nPrevAtomNumber = *pRankStack ++;

    AT_RANK *nNewRank        = NULL;
    AT_RANK *nNewAtomNumber  = NULL;
        
    if ( !pRankStack[0] ) {
        pRankStack[0] = (AT_RANK *) inchi_malloc(num_max*sizeof(*nNewRank));
    }
    if ( !pRankStack[1] ) {
        pRankStack[1] = (AT_RANK *) inchi_malloc(num_max*sizeof(*nNewAtomNumber));
    }
    if ( !pRankStack[0] || !pRankStack[1] )
        return CT_OUT_OF_RAM;  /*   <BRKPT> */
    nNewRank       = pRankStack[0];
    nNewAtomNumber = pRankStack[1];

    if ( nNewRank && nNewAtomNumber ) {
        memcpy( nNewAtomNumber, nPrevAtomNumber, num_atoms*sizeof(nNewAtomNumber[0]));
        memcpy( nNewRank, nPrevRank, num_atoms*sizeof(nNewRank[0]));
     
        for ( i = 1, nRet=0; i < num_atoms; i ++ ) { /*  12-12-2001: replaced Prev... with New... */
            if ( nNewRank[(int)nNewAtomNumber[i-1]] == nNewRank[(int)nNewAtomNumber[i]] ) {
                nNewRank[nNewAtomNumber[i-1]] = (AT_RANK)i;
                nNumRanks = DifferentiateRanks2( num_atoms, NeighList,
                                                 nNumRanks, nNewRank, nTempRank,
                                                 nNewAtomNumber, &pCS->lNumNeighListIter, 1 );
                pCS->lNumBreakTies ++;
                nRet ++;
            }
        }
    }
    return nRet;
}
