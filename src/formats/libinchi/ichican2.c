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
#include <ctype.h>
#include <limits.h>

/* #define CHECK_WIN32_VC_HEAP */
#include "mode.h"

#include "ichi.h"
#include "util.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "inpdef.h"
#include "ichinorm.h"
#include "ichicant.h"
#include "ichicano.h"
#include "ichicomn.h"

#include "ichicomp.h"

#define MAX_CELLS    1024
#define MAX_NODES    1024
#define MAX_SET_SIZE 2048 /*16384*/
#define MAX_LAYERS   7

#define INFINITY       0x3FFF
#define EMPTY_CT       0
#define EMPTY_H_NUMBER (INFINITY-1)
#define BASE_H_NUMBER  ((INFINITY-1)/2)
#define EMPTY_ISO_SORT_KEY LONG_MAX

#define SEPARATE_CANON_CALLS 0

/* #define INCHI_CANON_USE_HASH */
#define INCHI_CANON_MIN

/****************************************************************/
#ifdef INCHI_CANON_USE_HASH
typedef unsigned long  U_INT_32;
typedef unsigned char  U_INT_08;
typedef U_INT_32       CtHash;
CtHash hash_mark_bit;
#endif

/* -- moved to ichi_bns.h --
typedef U_SHORT  bitWord;
#define BIT_WORD_MASK  ((bitWord)~0)
*/

static bitWord *bBit = NULL;
static int    num_bit = 0;     
/*bitWord      mark_bit; */    /* highest bit in AT_NUMB */
/*bitWord      mask_bit; */    /* ~mark_bit */

AT_NUMB       rank_mark_bit;
AT_NUMB       rank_mask_bit;


typedef AT_NUMB    Node;
typedef NEIGH_LIST Graph;
/*
typedef struct tagGraph {
    int dummy;
} Graph;
*/
typedef struct tagUnorderedPartition {
    /* AT_NUMB *next; */ /* links */
    AT_NUMB *equ2; /* mcr */
} UnorderedPartition;

typedef struct tagCell {
    int       first; /* index of the first cell element in Partition::AtNumber[] */
    int       next;   /* next after the last index */
    int       prev;   /* position of the previously returned cell element */
} Cell;

#ifdef NEVER /* moved to ichi_bns.h */
typedef struct tagNodeSet {
    bitWord **bitword;
    int num_set; /* number of sets */
    int len_set; /* number of bitWords in each set */
} NodeSet;
#endif

typedef struct tagTransposition {
    AT_NUMB *nAtNumb;
} Transposition;


typedef struct tagCTable {
    AT_RANK  *Ctbl;     /* connection table */
    /* Format-atoms:   atom_rank[k] neigh_rank[k][1]...neigh_rank[k][n]
                       1) atom_rank[k1] < atom_rank[k2]  <=> k1 < k2
              where    2) atom_rank[k] > neigh_rank[k][i], i=1..n
                       3) neigh_rank[k][i] < neigh_rank[k][j]  <=>  i < j

       Format-tgroup:  tgroup_rank[k] endpoint_rank[k][1]...endpoint_rank[k][n]
              where    1) tgroup_rank[k1] < tgroup_rank[k2] <=> k1 < k2
                       2) endpoint_rank[k][i] < endpoint_rank[k][j] <=> i < j
                           
              Note:    tgroup_rank[k] > endpoint_rank[k][j] for all j by construction
    */

    int       lenCt;        /* used length */
    int       nLenCTAtOnly; /* to split Ctnl comparison in case of bDigraph != 0 */
    int       maxlenCt;     /* allocated length of Ctbl */
    int       maxPos;       /* allocated length of nextCtblPos */
    int       maxVert;      /* max number of vertices to separate atoms from taut groups */
    int       lenPos;       /* first unused element of nextCtblPos */
    AT_RANK  *nextAtRank;   /* rank (k value) after the last node of the Ctbl portion*/
    AT_NUMB  *nextCtblPos;  /* first unused element of Ctbl */

    /* hydrogen atoms fixed in tautomeric representation:
       compare before diff sign inversion: (+) <=> Ct1->() > Ct2->() */
    NUM_H          *NumH;      
    int             lenNumH;    /* used length */
    int             maxlenNumH; /*  n + T_NUM_NO_ISOTOPIC*(n_tg-n) + 1 */

    /* hydrogen atoms fixed in non-tautomeric representation only:
       compare before diff sign inversion: (+) <=> Ct1->() > Ct2->() */
    NUM_H           *NumHfixed;          
    /*int              lenNumHfixed;    */   /* used length */   
    /*int              maxlenNumHfixed; */   /* max length = n+1  */

    /* isotopic atoms (without tautomeric H) and isotopic tautomeric groups */
    /* note: AT_ISO_SORT_KEY and T_GROUP_ISOWT are identical types: long    */
    AT_ISO_SORT_KEY *iso_sort_key;        
    int              len_iso_sort_key;    /* used length */
    int              maxlen_iso_sort_key; /* max length = n_tg+1 */

    S_CHAR          *iso_exchg_atnos;
    int              len_iso_exchg_atnos;
    int              maxlen_iso_exchg_atnos;

    /* isotopic hydrogen atoms fixed in non-tautomeric representation only */
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    AT_ISO_SORT_KEY *iso_sort_key_Hfixed;     
    int              len_iso_sort_key_Hfixed;    /* used length */    
    int              maxlen_iso_sort_key_Hfixed; /* max length = n+1  */
#endif
    
#ifdef INCHI_CANON_USE_HASH
    CtHash   *hash;
#endif
} ConTable;
/**************************************************************/
typedef struct tagkLeast {
    int k;
    int i;
} kLeast;
/*************local prototypes **********************************/
int CanonGraph( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  );

void CtPartFill( Graph *G, CANON_DATA *pCD, Partition *p,
                 ConTable *Ct, int k, int n, int n_tg );
void CtPartClear( ConTable *Ct, int k );
void CtPartInfinity( ConTable *Ct, S_CHAR *cmp, int k );
int CtPartCompare( ConTable *Ct1, ConTable *Ct2, S_CHAR *cmp,
                   kLeast *kLeastForLayer, int k, int bOnlyCommon, int bSplitTautCompare );
int CtFullCompare( ConTable *Ct1, ConTable *Ct2, int bOnlyCommon, int bSplitTautCompare );
void CtPartCopy( ConTable *Ct1 /* to */, ConTable *Ct2 /* from */, int k );
void CtFullCopy( ConTable *Ct1, ConTable *Ct2 );

int CtFullCompareLayers( kLeast *kLeastForLayer );
int CtCompareLayersGetFirstDiff( kLeast *kLeast_rho, int nOneAdditionalLayer,
                                 int *L_rho, int *I_rho, int *k_rho );
int CtPartCompareLayers( kLeast *kLeast_rho, int L_rho_fix_prev, int nOneAdditionalLayer );
void UpdateCompareLayers( kLeast kLeastForLayer[], int hzz );
int GetOneAdditionalLayer( CANON_DATA *pCD, ConTable *pzb_rho_fix );

void CleanNumH( NUM_H *NumH, int len );
int CleanCt( AT_RANK *Ct, int len );
void CleanIsoSortKeys( AT_ISO_SORT_KEY * isk, int len );
void MergeCleanIsoSortKeys( AT_ISO_SORT_KEY * isk1, AT_ISO_SORT_KEY * isk2, int len );

int  UnorderedPartitionJoin( UnorderedPartition *p1, UnorderedPartition *p2, int n );
Node GetUnorderedPartitionMcrNode( UnorderedPartition *p1, Node v );
int  nJoin2Mcrs2( AT_RANK *nEqArray, AT_RANK n1, AT_RANK n2 );
AT_RANK nGetMcr2( AT_RANK *nEqArray, AT_RANK n );
int  AllNodesAreInSet(  NodeSet *cur_nodes, int lcur_nodes, NodeSet *set, int lset );
void NodeSetFromVertices( NodeSet *cur_nodes, int l, Node *v, int num_v);
void CellMakeEmpty( Cell *baseW, int k );
Node CellGetMinNode( Partition *p, Cell *W, Node v, CANON_DATA *pCD );
int  CellGetNumberOfNodes( Partition *p, Cell *W );
int  CellIntersectWithSet( Partition *p, Cell *W, NodeSet *Mcr, int l );

int  PartitionColorVertex( Graph *G, Partition *p, Node v, int n, int n_tg, int n_max, int bDigraph, int nNumPrevRanks );
void PartitionCopy( Partition *To, Partition *From, int n );
int  PartitionSatisfiesLemma_2_25( Partition *p, int n );
void PartitionGetTransposition( Partition *pFrom, Partition *pTo, int n, Transposition *gamma );
void PartitionGetMcrAndFixSet( Partition *p, NodeSet *Mcr, NodeSet *Fix, int n, int l );
int  PartitionGetFirstCell( Partition *p, Cell *baseW, int k, int n );
int  PartitionIsDiscrete( Partition *p, int n);
void PartitionFree( Partition *p );
int  PartitionCreate( Partition *p, int n);
void UnorderedPartitionMakeDiscrete( UnorderedPartition *p, int n);
void UnorderedPartitionFree( UnorderedPartition *p );
int  UnorderedPartitionCreate( UnorderedPartition *p, int n );
void CTableFree( ConTable *Ct );
int  CTableCreate( ConTable *Ct, int n, CANON_DATA *pCD );
void TranspositionFree( Transposition *p );
int  TranspositionCreate( Transposition *p, int n );
void TranspositionGetMcrAndFixSetAndUnorderedPartition( Transposition *gamma, NodeSet *McrSet, NodeSet *FixSet, int n, int l, UnorderedPartition *p );

void insertions_sort_NeighList_AT_NUMBERS2( NEIGH_LIST base, AT_RANK *nRank, AT_RANK max_rj );
int  WriteGraph( Graph *G, int n, int gnum, char *fname, char *fmode );

int  SetInitialRanks2( int num_atoms, ATOM_INVARIANT2* pAtomInvariant2, AT_RANK *nNewRank, AT_RANK *nAtomNumber );
void FillOutAtomInvariant2( sp_ATOM* at, int num_atoms, int num_at_tg, ATOM_INVARIANT2* pAtomInvariant,
                           int bIgnoreIsotopic, int bHydrogensInRanks, int bHydrogensFixedInRanks,
                           int bDigraph, int bTautGroupsOnly, T_GROUP_INFO *t_group_info );

int  GetCanonRanking2( int num_atoms, int num_at_tg, int num_max, int bDigraph, sp_ATOM* at,
                     AT_RANK **pRankStack,  int nNumPrevRanks,
                     AT_RANK *nSymmRank,  AT_RANK *nCanonRank,
                     NEIGH_LIST *NeighList, AT_RANK *nTempRank,
                     CANON_STAT* pCS );


#if ( SEPARATE_CANON_CALLS == 1 )
/* for profiling purposes */

int CanonGraph01( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph02( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph03( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph04( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph05( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph06( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph07( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph08( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph09( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph10( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph11( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
int CanonGraph12( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{
    return
    CanonGraph(     n,     n_tg,     n_max,     bDigraph,        G,           pi  ,
                         nSymmRank,           nCanonRank,          nAtomNumberCanon,
                            pCD,               pCC,
                           pp_zb_rho_inp,            pp_zb_rho_out  );
}
#else

#define CanonGraph01 CanonGraph
#define CanonGraph02 CanonGraph
#define CanonGraph03 CanonGraph
#define CanonGraph04 CanonGraph
#define CanonGraph05 CanonGraph
#define CanonGraph06 CanonGraph
#define CanonGraph07 CanonGraph
#define CanonGraph08 CanonGraph
#define CanonGraph09 CanonGraph
#define CanonGraph10 CanonGraph
#define CanonGraph11 CanonGraph
#define CanonGraph12 CanonGraph

#endif

#ifdef INCHI_CANON_USE_HASH        
/****************************************************************/
static call_fill_crc32_data = 1;
static U_INT_32 crc32_data[256];

void fill_crc32_data()
{
  U_INT_32 c;
  int n, k;
 
  for (n = 0; n < 256; n++)
  {
    c = (U_INT_32)n;
    for (k = 0; k < 8; k++)
      c = c & 1 ? 0xEDB88320L ^ (c >> 1) : c >> 1;
    crc32_data[n] = c;
  }
  call_fill_crc32_data = 0;
}
/****************************************************************/
unsigned long add2crc32( unsigned long crc32, AT_NUMB n )
{
    U_INT_08 chr;
    if (call_fill_crc32_data) {
      fill_crc32_data();
    }
    chr = n % 128;
    crc32 = crc32_data[((int)crc32 ^ (int)chr) & 0xff] ^ (crc32 >> 8);
    chr = n / 128;
    crc32 = crc32_data[((int)crc32 ^ (int)chr) & 0xff] ^ (crc32 >> 8);
    return crc32;
}
#endif
/****************************************************************/
int TranspositionCreate( Transposition *p, int n )
{
    p->nAtNumb = (AT_NUMB*)inchi_calloc( n, sizeof(p->nAtNumb[0]) );
    if ( p->nAtNumb ) {
        return 1;
    }
    return 0;
}
/****************************************************************/
void TranspositionFree( Transposition *p )
{
    if ( p && p->nAtNumb ) {
        inchi_free( p->nAtNumb );
        p->nAtNumb = NULL;
    }
}
/****************************************************************/
int NodeSetCreate( NodeSet *pSet, int n, int L )
{
    int i, len;

    len = (n+ num_bit - 1)/num_bit;

    pSet->bitword = (bitWord**)inchi_calloc(L, sizeof(pSet->bitword[0]));

    if ( !pSet->bitword ) {
        return 0;
    }
    pSet->bitword[0] = (bitWord*)inchi_calloc(len*L, sizeof(pSet->bitword[0][0]));
    if ( !pSet->bitword[0] ) {
        /* cleanup */
        inchi_free( pSet->bitword );
        pSet->bitword = NULL;
        return 0; /* failed */
    }
    for ( i = 1; i < L; i ++ ) {
        pSet->bitword[i] = pSet->bitword[i-1]+len;
    }
    pSet->len_set = len;
    pSet->num_set = L;
    return 1;
}
/****************************************************************/
void NodeSetFree( NodeSet *pSet )
{
    if ( pSet && pSet->bitword ) {
        if ( pSet->bitword[0] ) {
            inchi_free( pSet->bitword[0] );
        }
        inchi_free( pSet->bitword );
        pSet->bitword = NULL;
    }
}
/****************************************************************/
int CTableCreate( ConTable *Ct, int n, CANON_DATA *pCD )
{
    int maxlenCt        = pCD->nMaxLenLinearCT + 1; /* add one element for CtPartInfinity() */
    int maxlenNumH      = pCD->NumH?       (pCD->maxlenNumH + 1)      : 0;
    int maxlenNumHfixed = pCD->NumHfixed?  (pCD->maxlenNumHfixed + 1) : 0;
    int maxlenIso       = pCD->maxlen_iso_sort_key? (pCD->maxlen_iso_sort_key+1) : 0;
    int maxlenIsoExchg  = pCD->iso_exchg_atnos? (pCD->maxlen_iso_exchg_atnos+1) : 0;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    int maxlenIsoHfixed = pCD->maxlen_iso_sort_key_Hfixed? (pCD->maxlen_iso_sort_key_Hfixed+1):0;
#endif
    
    memset( Ct, 0, sizeof(Ct[0]) );
    
    Ct->maxVert = n;
    
    n ++;

    Ct->Ctbl           = (AT_RANK*) inchi_calloc(maxlenCt, sizeof(Ct->Ctbl[0]) );
    Ct->nextCtblPos    = (AT_NUMB*) inchi_calloc(n, sizeof(Ct->nextCtblPos[0]) );
    Ct->nextAtRank     = (AT_RANK*) inchi_calloc(n, sizeof(Ct->nextAtRank[0]) );
    if ( maxlenNumH ) {
        Ct->NumH            = (NUM_H *)  inchi_calloc(maxlenNumH, sizeof(Ct->NumH[0]));
    }
    if ( maxlenNumHfixed ) {
        Ct->NumHfixed = (NUM_H *)  inchi_calloc(maxlenNumHfixed, sizeof(Ct->NumH[0]));
    }
    if ( maxlenIso ) {
        Ct->iso_sort_key = (AT_ISO_SORT_KEY *)inchi_calloc(maxlenIso, sizeof(Ct->iso_sort_key[0]));
    }
    if ( maxlenIsoExchg ) {
        Ct->iso_exchg_atnos = (S_CHAR *)inchi_calloc( maxlenIsoExchg, sizeof(Ct->iso_exchg_atnos[0]));
    }
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    if ( maxlenIsoHfixed ) {
        Ct->iso_sort_key_Hfixed = (AT_ISO_SORT_KEY *)inchi_calloc(maxlenIsoHfixed, sizeof(Ct->iso_sort_key_Hfixed[0]));
    }
#endif
#ifdef INCHI_CANON_USE_HASH
    Ct->hash           = (CtHash*) inchi_calloc(n, sizeof(Ct->hash[0]) );
#endif
    Ct->lenCt          = 0;
    Ct->nLenCTAtOnly   = pCD->nLenCTAtOnly;
    Ct->maxlenCt       = maxlenCt;
    Ct->lenNumH        = 0;
    Ct->maxlenNumH     = maxlenNumH;
    Ct->len_iso_sort_key           = 0;
    Ct->maxlen_iso_sort_key        = maxlenIso;
    Ct->len_iso_exchg_atnos        = 0;
    Ct->maxlen_iso_exchg_atnos     = maxlenIso;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    Ct->len_iso_sort_key_Hfixed    = 0;
    Ct->maxlen_iso_sort_key_Hfixed = maxlenIsoHfixed;
#endif
    Ct->maxPos         = n;
    Ct->lenPos         = 0;
    Ct->nextAtRank[0]  = 0;
    Ct->nextCtblPos[0] = 0;
    if ( Ct->Ctbl && Ct->nextCtblPos &&
         (!maxlenNumH || Ct->NumH)   &&
         (!maxlenNumHfixed || Ct->NumHfixed ) ) {
        return 1;
    }
    return 0;
}
/****************************************************************/
void CTableFree( ConTable *Ct )
{
    if ( Ct ) {
        if ( Ct->Ctbl )
            inchi_free( Ct->Ctbl );
        if ( Ct->nextCtblPos )
            inchi_free( Ct->nextCtblPos );
        if ( Ct->nextAtRank )
            inchi_free( Ct->nextAtRank );
        if ( Ct->NumH )
            inchi_free( Ct->NumH );
        if ( Ct->NumHfixed )
            inchi_free( Ct->NumHfixed );
        if ( Ct->iso_sort_key )
            inchi_free( Ct->iso_sort_key );
        if ( Ct->iso_exchg_atnos )
            inchi_free( Ct->iso_exchg_atnos );
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
        if ( Ct->iso_sort_key_Hfixed )
            inchi_free( Ct->iso_sort_key_Hfixed );
#endif
#ifdef INCHI_CANON_USE_HASH
        if ( Ct->hash )
            inchi_free( Ct->hash );
#endif
        memset( Ct, 0, sizeof( Ct[0] ) );
    }
}
/****************************************************************/
int UnorderedPartitionCreate( UnorderedPartition *p, int n )
{
    p->equ2 = (AT_NUMB*)inchi_calloc( n, sizeof(p->equ2[0]));
    /* p->next = (AT_NUMB*)inchi_calloc( n, sizeof(p->next[0])); */
    if ( p->equ2 /*&& p->next*/ )
        return 1;
    return 0;
}
/****************************************************************/
void UnorderedPartitionFree( UnorderedPartition *p )
{
    if (p->equ2) inchi_free(p->equ2);
    /* if (p->next) inchi_free(p->next); */
    p->equ2 = NULL;
    /* p->next = NULL; */
}
/****************************************************************/
void UnorderedPartitionMakeDiscrete( UnorderedPartition *p, int n)
{
    int i;
    for ( i = 0; i < n; i ++ ) {
        p->equ2[i] = (AT_NUMB)i;
        /* p->next[i] = INFINITY; */
    }
    INCHI_HEAPCHK
}
/****************************************************************/
int PartitionCreate( Partition *p, int n)
{
    p->AtNumber = (AT_NUMB*)inchi_calloc( n, sizeof(p->AtNumber[0]));
    p->Rank     = (AT_RANK*)inchi_calloc( n, sizeof(p->Rank[0]));
    if ( p->AtNumber && p->Rank ) {
        return 1;
    }
    return 0;
}
/****************************************************************/
void PartitionFree( Partition *p )
{
    if ( p ) {
        if ( p->AtNumber ) {
            inchi_free( p->AtNumber );
            p->AtNumber = NULL;
        }
        if ( p->Rank ) {
            inchi_free( p->Rank );
            p->Rank = NULL;
        }
    }
}
/****************************************************************/
int PartitionIsDiscrete( Partition *p, int n)
{
    int i;
    AT_RANK r;
    for ( i = 0, r = 1; i < n; i ++, r ++ ) {
        if ( r != (rank_mask_bit & p->Rank[p->AtNumber[i]]) ) {
            INCHI_HEAPCHK
            return 0;
        }
    }
    INCHI_HEAPCHK
    return 1;
}
/****************************************************************/
int PartitionGetFirstCell( Partition *p, Cell *baseW, int k, int n )
{
    int i;
    AT_RANK r;
    Cell *W = baseW+k-1;

    i = (k > 1)? baseW[k-2].first+1 : 0;
    if ( i < n ) {
        /* bypass single vertex cells */
        for ( r = (AT_RANK)(i+1); i < n && r == (rank_mask_bit & p->Rank[(int)p->AtNumber[i]]); i ++, r++ )
            ;
    }
    if ( i < n ) {
        W->first = i;
        for ( r = (rank_mask_bit & p->Rank[(int)p->AtNumber[i]]), i++ ;
                i < n && r == (rank_mask_bit & p->Rank[(int)p->AtNumber[i]]);
                   i ++ )
            ;
        W->next = i;
        INCHI_HEAPCHK
        return (W->next - W->first);
    }
    W->first = INFINITY;
    W->next  = 0;
    INCHI_HEAPCHK
    return 0;
}
/****************************************************************/
void CellMakeEmpty( Cell *baseW, int k )
{
    k --;
    baseW[k].first = INFINITY;
    baseW[k].next  = 0;
    baseW[k].prev  = -1;
    INCHI_HEAPCHK
}
/****************************************************************/
void NodeSetFromVertices( NodeSet *cur_nodes, int l, Node *v, int num_v)
{
    bitWord *Bits = cur_nodes->bitword[l-1];
    int      len  = cur_nodes->len_set*sizeof(bitWord);
    int      i, j;

    memset( Bits, 0, len );

    for ( i = 0; i < num_v; i ++ ) {
        j = (int)v[i]-1;
        Bits[ j / num_bit ] |= bBit[ j % num_bit ];
    }
    INCHI_HEAPCHK
}
/****************************************************************/
int AllNodesAreInSet(  NodeSet *cur_nodes, int lcur_nodes, NodeSet *set, int lset )
{
    int i;
    int n = cur_nodes->len_set;
    bitWord *BitsNode = cur_nodes->bitword[lcur_nodes-1];
    bitWord *BitsSet  = set->bitword[lset-1];
    /* find any BitsNode[i] bit not in BitsSet[i] */
    for ( i = 0; i < n; i ++ ) {
        if ( BitsNode[i] & ~BitsSet[i] ) {
            INCHI_HEAPCHK
            return 0;
        }
    }
    INCHI_HEAPCHK
    return 1;
}
/****************************************************************/
void PartitionGetMcrAndFixSet( Partition *p, NodeSet *Mcr, NodeSet *Fix, int n, int l )
{
    int i, j1, j2;
    AT_RANK r, r1;
    bitWord *McrBits = Mcr->bitword[l-1];
    bitWord *FixBits = Fix->bitword[l-1];
    int     len      = Mcr->len_set*sizeof(bitWord);

    memset( McrBits, 0, len );
    memset( FixBits, 0, len );
    for ( i = 0, r = 1; i < n; i ++, r ++ ) {
        if ( r == (r1=(rank_mask_bit&p->Rank[j1=(int)p->AtNumber[i]])) ) {
            FixBits[j1 / num_bit] |= bBit[j1 % num_bit];
            McrBits[j1 / num_bit] |= bBit[j1 % num_bit];
        } else {
            for ( r = r1; i+1 < n && r == (rank_mask_bit&p->Rank[j2=(int)p->AtNumber[i+1]]); i ++ ) {
                if ( j1 > j2 ) {
                    j1 = j2;
                }
            }
            McrBits[j1 / num_bit] |= bBit[j1 % num_bit];
        }
    }
    INCHI_HEAPCHK
}
/************* used in ichi_bns.c ********************************/
void NodeSetFromRadEndpoints( NodeSet *cur_nodes, int k, /*Node *v*/ Vertex RadEndpoints[], int num_v)
{
    bitWord *Bits = cur_nodes->bitword[k];
    int      len  = cur_nodes->len_set*sizeof(bitWord);
    int      i, j;

    memset( Bits, 0, len );

    for ( i = 1; i < num_v; i += 2 ) {
        j = (int)RadEndpoints[i];
        Bits[ j / num_bit ] |= bBit[ j % num_bit ];
    }
}
/************* used in ichi_bns.c ********************************/
void RemoveFromNodeSet( NodeSet *cur_nodes, int k, Vertex v[], int num_v)
{
    if ( cur_nodes->bitword ) {
        bitWord *Bits = cur_nodes->bitword[k];
        /*int      len  = cur_nodes->len_set*sizeof(bitWord);*/
        int      i, j;

        for ( i = 0; i < num_v; i ++ ) {
            j = (int) v[i];
            Bits[ j / num_bit ] &= ~bBit[ j % num_bit ];
        }
    }
}
/************* used in ichi_bns.c ********************************/
int DoNodeSetsIntersect( NodeSet *cur_nodes, int k1, int k2)
{
    if ( cur_nodes->bitword ) {
        bitWord *Bits1 = cur_nodes->bitword[k1];
        bitWord *Bits2 = cur_nodes->bitword[k2];
        int      len  = cur_nodes->len_set;
        int      i;

        for ( i = 0; i < len; i ++ ) {
            if ( Bits1[i] & Bits2[i] )
                return 1;
        }
    }
    return 0;
}
/************* used in ichi_bns.c ********************************/
int IsNodeSetEmpty( NodeSet *cur_nodes, int k)
{
    if ( cur_nodes->bitword ) {
        bitWord *Bits = cur_nodes->bitword[k];
        int      len  = cur_nodes->len_set;
        int      i;

        for ( i = 0; i < len; i ++ ) {
            if ( Bits[i] )
                return 0;
        }
    }
    return 1;
}
/************* used in ichi_bns.c ********************************/
void AddNodeSet2ToNodeSet1( NodeSet *cur_nodes, int k1, int k2)
{
    if ( cur_nodes->bitword ) {
        bitWord *Bits1 = cur_nodes->bitword[k1];
        bitWord *Bits2 = cur_nodes->bitword[k2];
        int      len  = cur_nodes->len_set;
        int      i;

        for ( i = 0; i < len; i ++ ) {
            Bits1[i] |= Bits2[i];
        }
    }
}
/************* used in ichi_bns.c ********************************/
int AddNodesToRadEndpoints( NodeSet *cur_nodes, int k, Vertex RadEndpoints[], Vertex vRad, int nStart, int nLen )
{
    int n = nStart;
    if ( cur_nodes->bitword ) {
        bitWord *Bits = cur_nodes->bitword[k];
        int      len  = cur_nodes->len_set;
        int      i, j;
        Vertex   v;

        for ( i = 0, v = 0; i < len; i ++ ) {
            if ( Bits[i] ) {
                for ( j = 0; j < num_bit; j ++, v ++ ) {
                    if ( Bits[i] & bBit[j] ) {
                        if ( n >= nLen ) {
                            return -1; /* overflow */
                        }
                        RadEndpoints[n ++] = vRad;
                        RadEndpoints[n ++] = v;
                    }
                }
            } else {
                v += num_bit;
            }
        }
    }
    return n;
}
/****************************************************************/
void PartitionGetTransposition( Partition *pFrom, Partition *pTo, int n, Transposition *gamma )
{
    int i;
    for ( i = 0; i < n; i ++ ) {
        gamma->nAtNumb[(int)pFrom->AtNumber[i]] =pTo->AtNumber[i];
    }
    INCHI_HEAPCHK
}
/**************************************************************************************/
/*  Get minimal set (class) representative and partially compress the partitioning */
/*  mcr = minimal class representative. */
AT_RANK nGetMcr2( AT_RANK *nEqArray, AT_RANK n )
{
    AT_RANK n1, n2, mcr; /*  recursive version is much shorter. */
    INCHI_HEAPCHK
    n1=nEqArray[(int)n];
    if ( n == n1 ) {
        return n;
    }
    /*  1st pass: find mcr */
    while ( n1 != (n2=nEqArray[(int)n1])) {
        n1 = n2;
    }
    /*  2nd pass: copy mcr to each element of the set starting from nEqArray[n] */
    mcr = n1;
    n1  = n;
    while ( /*n1*/ mcr != (n2=nEqArray[(int)n1]) ) {
        nEqArray[(int)n1]=mcr;
        n1 = n2;
    }
    INCHI_HEAPCHK
    return ( mcr );
}
/**************************************************************************************/
/*  Join 2 sets (classes) that have members n1 and n2 */
int nJoin2Mcrs2( AT_RANK *nEqArray, AT_RANK n1, AT_RANK n2 )
{
    n1 = nGetMcr2( nEqArray, n1 );
    n2 = nGetMcr2( nEqArray, n2 );
    if ( n1 < n2 ) {
        nEqArray[n2] = n1;
        INCHI_HEAPCHK
        return 1; /*  a change has been made */
    }
    if ( n2 < n1 ) {
        nEqArray[n1] = n2;
        INCHI_HEAPCHK
        return 1; /*  a change has been made */
    }
    INCHI_HEAPCHK
    return 0; /*  no changes */
}
/****************************************************************/
Node GetUnorderedPartitionMcrNode( UnorderedPartition *p1, Node v )
{
    Node ret = (Node)(1+ nGetMcr2( p1->equ2, (AT_RANK)(v-1) ));
    INCHI_HEAPCHK
    return ret;
}
/****************************************************************/
/* change p2 to (p2 v p1)  */
int UnorderedPartitionJoin( UnorderedPartition *p1, UnorderedPartition *p2, int n )
{
    int i, j;
    int nNumChanges = 0;
    for ( i = 0; i < n; i ++ ) {
        if ( (j=(int)p1->equ2[i]) == i || p2->equ2[(int)i] == p2->equ2[(int)j] ) {
            continue;
        }
        nNumChanges += nJoin2Mcrs2(p2->equ2, (AT_NUMB)i, (AT_NUMB)j );
    }
    INCHI_HEAPCHK
    return nNumChanges;
}
/****************************************************************/
int PartitionSatisfiesLemma_2_25( Partition *p, int n )
{
    int nPartitionSize = 0;
    int nNumNonTrivialCells = 0;
    AT_RANK r;
    int i, num;
    for ( i = num = 0, r=1; i < n; i ++, r++ ) {
        if ( (rank_mask_bit & p->Rank[(int)p->AtNumber[i]]) == r ) {
            nPartitionSize ++;
            if ( num ) {
                /* num+1 = cell size > 1 */
                nNumNonTrivialCells ++;
                num = 0;
            }
        } else {
            num ++;
        }
    }
    /* check Lemma_2_25 conditions */
    if ( n <= nPartitionSize+4 ||
         n == nPartitionSize + nNumNonTrivialCells ||
         n == nPartitionSize + nNumNonTrivialCells + 1 ) {
        return 1;
    }
    return 0;
}
/****************************************************************/
void PartitionCopy( Partition *To, Partition *From, int n )
{
    int i;
    memcpy( To->AtNumber, From->AtNumber, n*sizeof(To->AtNumber[0]));
    memcpy( To->Rank, From->Rank, n*sizeof(To->AtNumber[0]));
    for ( i = 0; i < n; i ++ ) {
        To->Rank[i] &= rank_mask_bit;
    }
    INCHI_HEAPCHK
}
/****************************************************************/
/* makes new equitable partition (p+1) out of p; first reduce the rank of vertex v */
int PartitionColorVertex( Graph *G, Partition *p, Node v, int n, int n_tg, int n_max, int bDigraph, int nNumPrevRanks )
{
    int     nNumNewRanks, i, j;
    long    lNumNeighListIter = 0;
    AT_RANK rv, r;
    AT_NUMB s, sv;
    for ( i = 1; i <= 2; i ++ ) {
        if ( !p[i].AtNumber ) {
            p[i].AtNumber = (AT_NUMB *) inchi_malloc(n_max*sizeof(p[0].AtNumber[0]));
        }
        if ( !p[i].Rank ) {
            p[i].Rank = (AT_RANK *) inchi_malloc(n_max*sizeof(p[0].Rank[0]));
        }
        if ( !p[i].AtNumber || !p[i].Rank ) {
            INCHI_HEAPCHK
            return CT_OUT_OF_RAM;
        }
    }
    PartitionCopy( p+1, p, n_tg );
    sv  = v-1;          /* atom number we are looking for */
    if ( sv >= (AT_NUMB) n_tg ) {
        INCHI_HEAPCHK
        return CT_CANON_ERR; /* !!! severe program error: sv not found !!! */
    }
    rv = p[1].Rank[(int)sv];  /* rank of this atom */
    /* second, locate sv among all vertices that have same rank as v */
    s = n_max + 1; /* always greater than sv; this initialization is needed only to keep the compiler happy */
    for ( j = (int)rv-1;  0 <= j && rv == (r = p[1].Rank[(int)(s=p[1].AtNumber[j])]) && s != sv; j -- )
        ;
    if ( s != sv ) {
        INCHI_HEAPCHK
        return CT_CANON_ERR; /* !!! severe program error: sv not found !!! */
    }
    /* shift preceding atom numbers to the right to fill the gap after removing sv */
    r = rv-1; /* initialization only to keep compiler happy */
    for ( i = j--;  0 <= j && rv == (r = p[1].Rank[(int)(s=p[1].AtNumber[j])]); i = j, j -- ) {
        p[1].AtNumber[i] = s;
    }
    r = (i > 0)? (r+1):1;  /* new reduced rank = (next lower rank)+1 or 1 */
    /* insert sv and adjust its rank */
    p[1].AtNumber[i] = sv;
    p[1].Rank[(int)sv] = r;


    /* make equitable partition */
    if ( bDigraph ) {

        /*
        nNumNewRanks = DifferentiateRanks2( n_tg, G,
                                         nNumPrevRanks+1, p[1].Rank, p[2].Rank,
                                         p[1].AtNumber, &lNumNeighListIter, 1 );
        */
        nNumNewRanks = DifferentiateRanks4( n_tg, G, 
                                         nNumPrevRanks+1, p[1].Rank, p[2].Rank /* temp array */,
                                         p[1].AtNumber,  (AT_RANK)n, &lNumNeighListIter );
    
    
    } else {
        /*
        nNumNewRanks = DifferentiateRanks2( n_tg, G,
                                         nNumPrevRanks+1, p[1].Rank, p[2].Rank,
                                         p[1].AtNumber, &lNumNeighListIter, 1 );
        */
        nNumNewRanks = DifferentiateRanks3( n_tg, G,
                                         nNumPrevRanks+1, p[1].Rank, p[2].Rank /* temp array */,
                                         p[1].AtNumber, &lNumNeighListIter );
    }
    INCHI_HEAPCHK

    return nNumNewRanks;
}
typedef struct tagNodeValues {
    NUM_H            NumH;
    AT_ISO_SORT_KEY  iso_sort_key;
    NUM_H            NumHfixed;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    AT_ISO_SORT_KEY  iso_sort_key_Hfixed;
#endif
    AT_NUMB          nAtNumber;
} NV;

/****************************************************************/
/* return min node > vPrev or INFINITY if not found */
/* Input: v = previous atom number + 1 or 0 on first call*/
Node CellGetMinNode( Partition *p, Cell *W, Node v, CANON_DATA *pCD )
{
    AT_NUMB i;
    AT_NUMB uCurAtNumb, uMinAtNumb = INFINITY;
    /* in case of emty cell:  (W->first=INFINITY) > (W->next=0); returns INFINITY */
    if ( W->first > W->next ) {
        return INFINITY;
    }
#if ( USE_AUX_RANKING == 1 )
    if ( pCD && pCD->nAuxRank )
    {
        AT_RANK uMinAuxRank, uCurAuxRank;
        int     nCurAtNumb;
#if ( USE_AUX_RANKING_ALL == 1 )
        AT_RANK uInpAuxRank;
        int     nInpAtNumb, nMinAtNumb;
#endif
        for ( i = W->first; i < W->next; i ++ ) {
            uCurAtNumb = p->AtNumber[(int)i];
            if ( !(p->Rank[(int)uCurAtNumb] & rank_mark_bit) ) {
                break; /* found the first unmarked yet node */
            }
        }
        if ( i == W->next ) {
            return INFINITY;
        }
#if ( USE_AUX_RANKING_ALL == 1 )
        /*==== vertex ordering definition ===
         * vertex v1 < v2 <=> (AuxRank(v1)==AuxRank(v2) && AtNumb(v1) < AtNumb(v2)) || (AuxRank(v1) < AuxRank(v2))
         * vertex v1 > v2 <=> (AuxRank(v1)==AuxRank(v2) && AtNumb(v1) > AtNumb(v2)) || (AuxRank(v1) > AuxRank(v2))
         * vertex v1 = v2 <=> (AuxRank(v1)==AuxRank(v2) && AtNumb(v1) == AtNumb(v2))
         */
        
        /* set initial vMin so that vMin > any vertex */
        uMinAuxRank = INFINITY;
        nMinAtNumb  = INFINITY;
        /* set vInp */
        if ( v ) {
            nInpAtNumb  = (int)v - 1; /* previous vertex */
            uInpAuxRank = pCD->nAuxRank[nInpAtNumb];
        } else {
            nInpAtNumb  = -1; /* less than any vertex */
            uInpAuxRank =  0;
        }
        /* find vMin = min { vCur : (vCur > vInp) && (vCur in W) } */ 
        for ( ; i < W->next; i ++ ) {
            nCurAtNumb = (int)p->AtNumber[(int)i];
            if ( !(p->Rank[nCurAtNumb] & rank_mark_bit) ) {
                /* vertex nCurAtNumb is not marked, find whether it fits the conditions */
                uCurAuxRank = pCD->nAuxRank[nCurAtNumb];
                if ( ((uCurAuxRank == uInpAuxRank) && (nCurAtNumb > nInpAtNumb)) ||
                      (uCurAuxRank > uInpAuxRank) ) {
                    /* here vCur > vInp */
                    if ( uCurAuxRank == uMinAuxRank && nCurAtNumb < nMinAtNumb ) {
                        /* vCur < vMin (1) */
                        nMinAtNumb = nCurAtNumb;
                    } else
                    if ( uCurAuxRank < uMinAuxRank ) {
                        /* vCur < vMin (2) */
                        uMinAuxRank = uCurAuxRank;
                        nMinAtNumb = nCurAtNumb;
                    }
                }
            }
        }
        uMinAtNumb = (nMinAtNumb==INFINITY)? INFINITY : (AT_NUMB)nMinAtNumb;
#else
        if ( v ) {
            nCurAtNumb = (int)v-1;
            /* any valid found node must have nAuxRank == uMinAuxRank */
            uMinAuxRank   = pCD->nAuxRank[nCurAtNumb];
        } else {
            /* any valid found node must have minimal uMinAuxRank from pCD->nAuxRank[] */
            uMinAuxRank   = INFINITY; /* undefined */
        }

        for ( ; i < W->next; i ++ ) {
            uCurAtNumb = p->AtNumber[(int)i];
            nCurAtNumb = (int)uCurAtNumb;
            if ( uCurAtNumb >= v && !(p->Rank[nCurAtNumb] & rank_mark_bit) ) {
                uCurAuxRank = pCD->nAuxRank[nCurAtNumb];
                if ( v ) {
                    /* get next node */
                    /* find node with smallest uCurAtNumb among nodes with aux. ranks equal to uMinAuxRank */
                    if ( uCurAuxRank == uMinAuxRank && uCurAtNumb < uMinAtNumb ) {
                        uMinAtNumb = uCurAtNumb;
                    }
                } else {
                    /* get first node */
                    /* find node with smallest smallest uCurAtNumb among nodes with smallest aux. ranks */
                    if ( uMinAuxRank > uCurAuxRank ) {
                        uMinAuxRank = uCurAuxRank;
                        uMinAtNumb = uCurAtNumb;
                    } else
                    if ( uMinAuxRank == uCurAuxRank && uCurAtNumb < uMinAtNumb ) {
                        uMinAtNumb = uCurAtNumb;
                    }
                }
            }
        }
#endif
    } else
#endif /* } USE_AUX_RANKING */
    {
        for ( i = W->first; i < W->next; i ++ ) {
            uCurAtNumb = p->AtNumber[(int)i];
            if ( uCurAtNumb >= v && !(p->Rank[(int)uCurAtNumb] & rank_mark_bit) && uCurAtNumb < uMinAtNumb ) {
                uMinAtNumb = uCurAtNumb;
            }
        }
    }
    if ( uMinAtNumb != INFINITY ) uMinAtNumb ++;
    INCHI_HEAPCHK
    return uMinAtNumb;
}
/****************************************************************/
int CellGetNumberOfNodes( Partition *p, Cell *W )
{
    int first = W->first;
    int next  = W->next;
    int i, num;
    for ( i = first, num = 0; i < next; i ++ ) {
        if ( !( rank_mark_bit & p->Rank[(int)p->AtNumber[i]] ) ) {
            num++;
        }
    }
    INCHI_HEAPCHK
    return num;
}
/****************************************************************/
int CellIntersectWithSet( Partition *p, Cell *W, NodeSet *Mcr, int l )
{
    bitWord *McrBits = Mcr->bitword[l-1];
    int first = W->first;
    int next  = W->next;
    int i, j, k;
    if ( first >= next ) { /* for testing only */
        return 0;
    }
    for ( i = first, k = 0; i < next; i ++ ) {
        j = (int)p->AtNumber[i];
        if ( !(McrBits[ j / num_bit ] & bBit[ j % num_bit ]) ) { /* BC: reading uninit memory ???-not examined yet */
            k += !(p->Rank[j] & rank_mark_bit); /* for testing only */
            p->Rank[j] |= rank_mark_bit;
        }
    }
    INCHI_HEAPCHK
    return k;
}
/****************************************************************/
void CtPartClear( ConTable *Ct, int k )
{
    int start;
    int len;
    /* connection table */
    start = k>1? Ct->nextCtblPos[k-1] : 0;
    len   = Ct->lenCt - start;
    if ( len > 0 ) {
        memset( Ct->Ctbl + start, 0, (Ct->lenCt - start)*sizeof(Ct->Ctbl[0]) );
    }
    Ct->lenCt = start;
    Ct->lenPos = k;

    INCHI_HEAPCHK
}
/**********************************************************************************/
/*  Sort neighbors according to ranks in ascending order */
void insertions_sort_NeighList_AT_NUMBERS2( NEIGH_LIST base, AT_RANK *nRank, AT_RANK max_rj )
{
  AT_NUMB *i, *j, *pk, tmp, rj;
  int k, num = (int)*base++;
  for( k=1, pk = base; k < num; k++, pk ++ ) {
     i = pk;
     j = i + 1;
     rj = (rank_mask_bit & nRank[(int)*j]);
     if ( rj < max_rj ) {
         while ( j > base && rj < (rank_mask_bit & nRank[(int)*i])) {
             tmp = *i;
             *i = *j;
             *j = tmp;
             j = i --;
         }
     }
  }
  INCHI_HEAPCHK
}
/****************************************************************/
/* may need previous Lambda */
void CtPartFill( Graph *G, CANON_DATA *pCD, Partition *p,
                 ConTable *Ct, int k, int n, int n_tg )
 /*  k = (new index in Ct->nextAtRank[] and Ct->nextCtblPos[]) + 1 */
{
    int     startCtbl;
    int     startAtOrd;
    AT_RANK r, rj, nn, j, rj_prev;
    int     i, m;
#ifdef INCHI_CANON_USE_HASH
    CtHash  hash = 0;
#endif
        static int count; /* for debug only */
        count ++;


    INCHI_HEAPCHK
    
    k --;
    if ( k ) {
        startCtbl  = Ct->nextCtblPos[k-1];
        startAtOrd = Ct->nextAtRank[k-1]-1;  /* here  p->Rank[p->AtNumber[r-1]] = r */
    } else {
        startCtbl  = 0;
        startAtOrd = 0;
    }
    /******* well-defined (by fixed ranks) part of the connection table ************/
    r = (rank_mask_bit & p->Rank[(int)p->AtNumber[startAtOrd]]);
    for ( i = startAtOrd; i < n_tg && r == (rank_mask_bit&p->Rank[m=(int)p->AtNumber[i]]); i++, r ++ ) {
        Ct->Ctbl[startCtbl++] = r;
        insertions_sort_NeighList_AT_NUMBERS2( G[m], p->Rank, r );
        nn = G[m][0]; /* number of neighbors */
        rj_prev = 0; /* debug only */
#ifdef INCHI_CANON_USE_HASH        
        hash = add2crc32( hash, (AT_NUMB)(r + n) );
#endif
        for ( j = 1; j <= nn && (rj=(rank_mask_bit&p->Rank[(int)G[m][j]])) < r; j ++ ) {
            Ct->Ctbl[startCtbl++] = rj;
#ifdef INCHI_CANON_USE_HASH        
            hash = add2crc32( hash, rj );
#endif
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
            /* debug only */
            if ( rj < rj_prev ) {
                int stop = 1;   /* <BRKPT> */
            }
#endif
            rj_prev = rj;
        }
    }
    
    INCHI_HEAPCHK

    /****************** well-defined part of base hydrogen atoms *******************/
    if ( pCD->NumH && Ct->NumH ) {
        nn = inchi_min(n, i);
        for ( j = startAtOrd; j < nn; j ++ ) { /* atoms */
            Ct->NumH[j] = pCD->NumH[p->AtNumber[j]];
        }
        for ( ; j < i; j ++ ) {  /* t-groups */
            int data_pos = n + T_NUM_NO_ISOTOPIC * ((int)p->AtNumber[j] - n);
            for ( m = 0; m < T_NUM_NO_ISOTOPIC; m ++ ) {
                Ct->NumH[nn ++] = pCD->NumH[data_pos ++];
            }
        }
        Ct->lenNumH = nn;
    } else {
        Ct->lenNumH = 0;
    }

    INCHI_HEAPCHK

    /****************** well-defined part of fixed hydrogen atoms *******************/
    if ( pCD->NumHfixed && Ct->NumHfixed ) {
        nn = inchi_min(n, i);
        for ( j = startAtOrd; j < nn; j ++ ) {
            Ct->NumHfixed[j] = pCD->NumHfixed[p->AtNumber[j]];

    INCHI_HEAPCHK

        }
        /* Ct->lenNumHfixed = nn; */
    } else {
        ;/* Ct->lenNumHfixed = 0; */
    }

    INCHI_HEAPCHK

    /****************** well-defined part of isotopic keys ***************************/
    if ( pCD->iso_sort_key && Ct->iso_sort_key ) {
        for ( j = startAtOrd; j < i; j ++ ) {
            Ct->iso_sort_key[j] = pCD->iso_sort_key[p->AtNumber[j]];
        }
        Ct->len_iso_sort_key = i;
    } else {
        Ct->len_iso_sort_key = 0;
    }

    INCHI_HEAPCHK

    /****************** well-defined part of isotopic iso_exchg_atnos ***************************/
    if ( pCD->iso_exchg_atnos && Ct->iso_exchg_atnos ) {
        for ( j = startAtOrd; j < i; j ++ ) {
            Ct->iso_exchg_atnos[j] = pCD->iso_exchg_atnos[p->AtNumber[j]];
        }
        Ct->len_iso_exchg_atnos = i;
    } else {
        Ct->len_iso_exchg_atnos = 0;
    }

    INCHI_HEAPCHK
    /******** well-defined part of isotopic keys for fixed hydrogen atoms ************/
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    if ( pCD->iso_sort_key_Hfixed && Ct->iso_sort_key_Hfixed ) {
        nn = inchi_min(n, i);
        for ( j = startAtOrd; j < nn; j ++ ) {
            Ct->iso_sort_key_Hfixed[j] = pCD->iso_sort_key_Hfixed[p->AtNumber[j]];
        }
        Ct->len_iso_sort_key_Hfixed = nn;
    } else {
        Ct->len_iso_sort_key_Hfixed = 0;
    }
#endif

    INCHI_HEAPCHK

    Ct->lenCt          = startCtbl; /* not aways increases */
    Ct->nextCtblPos[k] = startCtbl;
    Ct->nextAtRank[k]  = r;
    Ct->lenPos = k+1;
    /* the rest of the CTable */
#ifdef INCHI_CANON_USE_HASH        
    while ( i < n ) {
        r = (rank_mask_bit&p->Rank[m=(int)p->AtNumber[i]]);
        hash = add2crc32( hash, (AT_NUMB)(r + n) );
        r++;
        insertions_sort_NeighList_AT_NUMBERS2( G[m], p->Rank, r );
        nn = G[m][0];
        rj_prev = 0; /* debug only */
        for ( j = 1; j <= nn && (rj=(rank_mask_bit&p->Rank[(int)G[m][j]])) < r; j ++ ) {
            hash = add2crc32( hash, rj );
        }
        i ++;
    }
    Ct->hash[k] = hash;
#endif        

    INCHI_HEAPCHK
}
/****************************************************************/
void CtPartInfinity( ConTable *Ct, S_CHAR *cmp, int k )
{
    int     startCtbl;
    /*int     startAtOrd;*/
    k --;
    if ( k ) {
        startCtbl  = Ct->nextCtblPos[k-1];
        /*startAtOrd = Ct->nextAtRank[k-1]-1;*/  /* here  p->Rank[p->AtNumber[r-1]] = r */
        if ( cmp ) {
            memset( cmp, 0, k*sizeof(cmp[0]) );
        }
    } else {
        startCtbl  = 0;
        /*startAtOrd = 0;*/
    }
    if ( !startCtbl || Ct->Ctbl[startCtbl-1] != EMPTY_CT ) {
        Ct->Ctbl[startCtbl] = EMPTY_CT;
    }
    INCHI_HEAPCHK
}
/****************************************************************/
/* Return value:
  -1 <=> *Lambda1 < *Lambda2
   0 <=> *Lambda1 = *Lambda2
  +1 <=> *Lambda1 > *Lambda2

  Input: k+1 = value of level at which the comparison is executed
               (that is, in the calling program k(caller) = k+1)

  Stars (*) below mark the differences:

  bSplitTautCompare != 0 => directed graph; compare:
                   non-tautomeric part of CT in layer 0; (*)
                   non-tautomeric H          in layer 1; (*)
                   tautomeric part of CT & H in layer 2; (*)
                   fixed H                   in layer 3;
                   isotopic atoms, non-taut
                               H & t-groups  in layer 4;
                   fixed isotopic H          in layer 5; <- move to layer 4

  bSplitTautCompare == 0 => undirected graph; compare:
                   full CT                   in Layer 0; (*)
                   taut and non-taut H       in Layer 1; (*)
                   * nothing *               in layer 2; (*)
                   fixed H                   in layer 3;
                   isotopic atoms, non-taut
                               H & t-groups  in layer 4;
                   fixed isotopic H          in layer 5; <- move to layer 4

*/
int CtPartCompare( ConTable *Ct1, ConTable *Ct2, S_CHAR *cmp,
                   kLeast *kLeastForLayer, int k, int bOnlyCommon, int bSplitTautCompare )
{
    int     startCt1, endCt1, startCt2, endCt2; /*endCt,*/
    int     startAt1, endAt1, startAt2, endAt2; /*endCt,*/
    int     midCt /* end of atoms only Ct */, midNumH=0 /* end of atoms only NumH */, maxVert;
    int     diff, i, k1, k2, lenNumH, len_iso_sort_key, /*mid_iso_sort_key,*/ midAt;
    int     nLayer = 0;

    k --;
    i = -1;
    /* set kLeastForLayer[nLayer].k = (k+1) or -(k+1)
           kLeastForLayer[nLayer].i = iDiff
        if all the conditions are met:
        1) kLeastForLayer[nLayer].k = 0
        2) diff==0 for all layers < nLayer

        sign:
        if the final diff < 0 then kLeastForLayer[nLayer].k = -(k+1) else
        if the final diff > 0 then kLeastForLayer[nLayer].k = +(k+1)

        k+1 instead of k takes into account k--; statememt above)

        meaning:
        ========
        abs(kLeastForLayer[nLayer].k) is the greatest level k at which
        difference at layer nLayer are zeroes of hidden by differences in smaller nLayer.
        
        "Hidden by difference in smaller level" means that nLayer of comparison
        has not been reached because the difference was discovered at a previous layer.

  
        Lambda vs zf_zeta comparison
        =============================================
        accept only diff == 0

        Lambda vs pzb_rho and pzb_rho_fix comparison
        =============================================
        Maintain kLeastForLayer[] and kLeastForLayerFix[]

        The algorithm provides that pzb_rho(m-1) < pzb_rho(m) <= pzb_rho_fix

        Definition: pzb_rho(m-1) < pzb_rho(m) means that
        -----------------------------------------------
        pzb_rho(m-1)[nLayerCurr] == pzb_rho(m)[nLayerCurr] for nLayerCurr = 0..nLayerDiff-1
        pzb_rho(m-1)[nLayerDiff] <  pzb_rho(m)[nLayerDiff]

        Definition: pzb_rho(m-1)[nLayerDiff] <  pzb_rho(m)[nLayerDiff] means that
        -------------------------------------------------------------------------
        pzb_rho(m-1)[nLayerDiff][i]     == pzb_rho(m)[nLayerDiff][i] for i=0..iDdiff-1
        pzb_rho(m-1)[nLayerDiff][iDdiff] < pzb_rho(m)[nLayerDiff][iDdiff]

        This defines nLayerDiff(pzb1, pzb2) where pszb1 = pzb_rho(a), pzb2=pzb_rho(b) (a<b) or pzb_rho_fix
               and   iDdiff    (pzb1, pzb2).
        In case pzb_rho(m)[nLayerCurr] == pzb_rho_fix[nLayerCurr] for all non-NULL nLayerCurr in pzb_rho_fix,
           nLayerDiff(pzb_rho(m), pzb_rho_fix) = the first layer in pzb_rho(m) not present in pzb_rho_fix
           iDdiff    (pzb_rho(m), pzb_rho_fix) = -1
        Case when such a layer does not exist means program error

        Suppose L_rho = nLayerDiff(Lambda, pzb_rho(m))
                L_fix = nLayerDiff(Lambda, pzb_rho_fix)
                I_rho = iDdiff    (Lambda, pzb_rho(m))
                I_fix = iDdiff    (Lambda, pzb_rho_fix)
                kLeastForLayer determined from Lambda vs pzb_rho(m) comparison
        Then:

        1. Comparison Lambda vs pzb_rho_fix before reaching discrete partition
        ----------------------------------------------------------------------
        a)    0 < abs(kLeastForLayerFix[L_fix].k) <= k-1 (* in this case I_fix >= 0 *)  &&
              ((L_fix < L_rho) || (L_fix == L_rho && I_fix < I_rho))
              =>
              qzb_rho_fix = kLeastForLayerFix[L_fix].k if prevoiusly qzb_rho_fix == 0

        b)    otherwise do not change qzb_rho_fix, except the following:

        c)    Special case L_rho == L_fix && I_rho == I_fix. Let L=L_rho, I = I_rho.

              Compare 3 valirs: Lambda[L][I], pzb_rho(m)[L][I], pzb_rho_fix[L][I]
              The algorithm provides pzb_rho(m)[L][I] < pzb_rho_fix[L][I]
              (pzb_rho(m)[L][I]==pzb_rho_fix[L][I] <=> pzb_rho(m)[L][I]==pzb_rho_fix[L][I] 
               is impossible by construction)
              There are 3 possibilities:
              c1) Lambda[L][I]     < pzb_rho(m)[L][I]  < pzb_rho_fix[L][I] <=>
                  kLeastForLayer[L].k  < 0 && kLeastForLayerFix[L].k < 0
                  => qzb_rho := kLeastForLayer[L].k, reject too small Lambda
              c2) pzb_rho(m)[L][I] < Lambda[L][I]      < pzb_rho_fix[L][I]
                  kLeastForLayer[L].k  > 0 && kLeastForLayerFix[L].k < 0
                  => qzb_rho := kLeastForLayer[L].k, accept Lambda, rho:=nu
              c3) pzb_rho(m)[L][I] < pzb_rho_fix[L][I] < Lambda[L][I]
                  kLeastForLayer[L].k  > 0 && kLeastForLayerFix[L].k > 0
                  => qzb_rho_fix := kLeastForLayerFix[L].k, reject too big Lambda

              Case
                  kLeastForLayer[L].k  < 0 && kLeastForLayerFix[L].k > 0 is impossible
                  because it means
                  pzb_rho_fix < Lambda < pzb_rho(m) <=> pzb_rho_fix < pzb_rho(m)


            Case (c3) occurs in case of (a)
            Case (c1) 

        2. Comparison Lambda vs pzb_rho before reaching discrete partition
        ----------------------------------------------------------------------
        a) (L_rho < L_fix) || (L_rho == L_fix && I_rho < I_fix)  =>

           Lambda differs from pzb_rho(m) in the part of pzb_rho(m) that will never change
           qzb_rho = kLeastForLayer[L_rho].k; reject Labmda or accept pzb_rho(m+1):=Labmda

        b) (L_rho == L_fix && I_rho > I_fix) && kLeastForLayer[L_rho].k < 0
           Lambda < pzb_rho(m), therefore
           qzb_rho = kLeastForLayer[L_rho].k; reject Labmda

        c) (L_rho > L_fix) =>
           qzb_rho := 0 because more significant difference may be discovered
           in layer < L_rho later. The final comparison may be needed at the
           level of discrete partition.


    */

    if ( cmp ) {
        for ( i = 0; i <= k && !cmp[i]; i++ )
            ;
        if ( i < k ) {
            cmp[k] = cmp[i];
            return (int)cmp[i];
        }
    }
    k1 = Ct1->lenPos-1;
    k2 = Ct2->lenPos-1;

#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    if ( k > k1 || k > k2 ) {
        int stop = 1;
    }
#endif
    diff = 0;

    if ( k ) {
        startCt1  = Ct1->nextCtblPos[k-1];
        startCt2  = Ct2->nextCtblPos[k-1];
        startAt1  = Ct1->nextAtRank[k-1]-1;
        startAt2  = Ct2->nextAtRank[k-1]-1;
    } else {
        startCt1 = startCt2 = 0;
        startAt1 = startAt2 = 0;
    }

    endCt1 = Ct1->nextCtblPos[k];
    endCt2 = Ct2->nextCtblPos[k];
    endAt1 = (int)Ct1->nextAtRank[k]-1;
    endAt2 = (int)Ct2->nextAtRank[k]-1;

    maxVert = inchi_min(Ct1->maxVert, Ct2->maxVert);

#ifdef INCHI_CANON_USE_HASH        
    if ( !diff ) {
        if ( Ct1->hash[k] > Ct2->hash[k] )
            diff = 1;
        else
        if ( Ct1->hash[k] < Ct2->hash[k] )
            diff = -1;
    }
    if ( diff ) {
        goto done;
    }
#endif
    
    /************************** lengths **************************************************/
    if ( (diff = -(startCt1 - startCt2)) ) {
        /* comparing two INFINITY terminations */
        if ( bOnlyCommon &&
             startCt1 >= Ct1->nLenCTAtOnly && startCt2 >= Ct2->nLenCTAtOnly &&
             Ct1->Ctbl[startCt1] == EMPTY_CT && Ct2->Ctbl[startCt2] == EMPTY_CT ) {
            return 0;
        }
        if ( bOnlyCommon ) {
            startCt1 = startCt2 = inchi_min(startCt1, startCt2);
            startAt1 = startAt2 = inchi_min(startAt1, startAt2);
            if ( Ct1->lenCt == Ct2->lenCt ) {
                endCt1 = endCt2 = inchi_max(endCt1, endCt2);
                endAt1 = endAt2 = inchi_max(endAt1, endAt2);
            }
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
            else {
                int stop = 1;
            }
#endif
        } else
        /* comparing (taut tail) vs INFINITY termination -- ??? */
        if ( startCt1 > startCt2 &&
             Ct1->maxVert > Ct2->maxVert &&
             startAt2 == Ct2->maxVert ) {
            return 0;
        } else {
            goto done;
        }
    }
    
    lenNumH          = Ct1->lenNumH;
    len_iso_sort_key = Ct1->len_iso_sort_key;
    
    if ( (diff = -(endCt1 - endCt2)) ) { /* negative sign reproduces results for NSC=28393 */
        if ( bOnlyCommon ) {
            endCt1 = endCt2 = inchi_min(endCt1, endCt2);
            endAt1 = endAt2 = inchi_min(endAt1, endAt2);
            lenNumH = inchi_min(Ct1->lenNumH, Ct2->lenNumH);
            len_iso_sort_key = inchi_min(Ct1->len_iso_sort_key, Ct1->len_iso_sort_key);
        } else
        /* take care of case when comparing tautomeric vs non-tautomeric:
           since (taut)->maxVert > (non-taut)->maxVert, --???
           (taut)->maxlenCt  > (non-taut)->maxlenCt     --!!!
           compare up to min out of the two, ignoring INFINITY in the last position */
        if ( endCt1 > endCt2 && Ct1->maxlenCt > Ct2->maxlenCt ) {
            if ( endAt2 == Ct2->maxVert + 1 ) {
                /* remove INFINITY termination of the shorter CT */
                /* should never happen */
                endAt2 --;
                len_iso_sort_key = lenNumH = endAt1 = endAt2;
                endCt2 --;
                endCt1 = endCt2;
                diff = 0;
            } else
            if ( endAt2 == Ct2->maxVert ) {
                /* remove INFINITY termination of CT */
                len_iso_sort_key = lenNumH = endAt1 = endAt2;
                endCt1 = endCt2;
                diff = 0;
            } else {
                goto done;
            }
        } else {
            goto done;
        }
    }
    
    if ( bSplitTautCompare ) {
        midCt = inchi_min(Ct1->nLenCTAtOnly, Ct2->nLenCTAtOnly);
        if ( midCt > endCt1 ) {
            midCt = endCt1;
        }
        midAt = inchi_min(maxVert, endAt1); 
    } else {
        midCt = endCt1;
        midAt = endAt1;
    }

    /*endCt   = min(endCt1, endCt2);*/
    /*************************************************************************/
    /************ layer 0: connection table without tautomeric groups ********/
    /*************************************************************************/
    for ( i = startCt1; i < midCt && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i ++ )
    /*for ( i = startCt1; i < endCt && !(diff = (int)Ct1->Ctbl[i] - (int)Ct2->Ctbl[i]); i ++ )*/
        ;
    if ( i < midCt ) {
        diff = (int)Ct1->Ctbl[i] - (int)Ct2->Ctbl[i];
        goto done;
    }
    /*************************************************************************/
    /******** layer 1 NumH: H atoms without tautomeric H *********************/
    /*************************************************************************/
    nLayer ++;
    /*============= check limits for consistency  ==========*/
    if ( (diff = -(startAt1 - startAt2)) ) {
        goto done;   /* should not happen */
    }
    if ( (diff = -(endAt1 - endAt2)) ) {
        goto done;   /* should not happen */
    }
    /*============= comparison =============================*/
    if ( Ct1->NumH && Ct2->NumH ) {
        if ( endAt1 < maxVert ) {
            midNumH = lenNumH = endAt1;
        } else
        if ( bSplitTautCompare ) {
            midNumH = maxVert;
        } else {
            midNumH = lenNumH;
        }
        /* lenNumH = (endAt2 >= maxVert)? lenNumH : endAt2; */
        /* endAt1 = (endAt2 == n)? lenNumH : endAt2; */
        
        for ( i = startAt1; i < midNumH && Ct1->NumH[i] == Ct2->NumH[i]; i ++ )
            ;
        if ( i < midNumH ) {
            diff = (int)Ct1->NumH[i] - (int)Ct2->NumH[i];
            goto done;
        }
    }
    /*************************************************************************/
    /************** layer 2: tautomeric part of CT and tautomeric H **********/
    /*************************************************************************/
    nLayer ++;
    for ( i = midCt; i < endCt1 && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i ++ )
        ; /* compare tautomeric groups part of CT */
    if ( i < endCt1 ) {
        diff = (int)Ct1->Ctbl[i] - (int)Ct2->Ctbl[i];
        goto done;
    }
    if ( Ct1->NumH && Ct2->NumH ) {
        for ( i = midNumH; i < lenNumH && Ct1->NumH[i] == Ct2->NumH[i]; i ++ )
            ; /* compare tautomeric H */
        if ( i < lenNumH ) {
            diff = (int)Ct1->NumH[i] - (int)Ct2->NumH[i];
            i += endCt1 - midCt;
            goto done;
        }
    }
    /*************************************************************************/
    /************** layer 3: Fixed H atoms ***********************************/
    /*************************************************************************/
    nLayer ++;
    if ( Ct1->NumHfixed && Ct2->NumHfixed ) {
        for ( i = startAt1; i < midAt && Ct1->NumHfixed[i] == Ct2->NumHfixed[i]; i ++ )
            ;
        if ( i < midAt ) {
            diff = (int)Ct1->NumHfixed[i] - (int)Ct2->NumHfixed[i];
            goto done;
        }
    }
    /*************************************************************************/
    /************** layer 4: isotopic atoms H, incl. tautomeric **************/
    /*************************************************************************/
    nLayer ++;
    if ( Ct1->iso_sort_key && Ct2->iso_sort_key ) {
        for ( i = startAt1; i < endAt1 && Ct1->iso_sort_key[i] == Ct2->iso_sort_key[i]; i ++ )
            ;
        if ( i < endAt1 ) {
            diff = Ct1->iso_sort_key[i] > Ct2->iso_sort_key[i]? 1:-1;
            goto done;
        }
    }
    if ( Ct1->iso_exchg_atnos && Ct2->len_iso_exchg_atnos ) {
        for ( i = startAt1; i < endAt1 && Ct1->iso_exchg_atnos[i] == Ct2->iso_exchg_atnos[i]; i ++ )
            ;
        if ( i < endAt1 ) {
            diff = Ct1->iso_exchg_atnos[i] > Ct2->iso_exchg_atnos[i]? 1:-1;
            goto done;
        }
    }
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    /*************************************************************************/
    /************** layer 6: Fixed isotopic H atoms **************************/
    /*************************************************************************/
    nLayer ++;
    if ( Ct1->iso_sort_key_Hfixed && Ct2->iso_sort_key_Hfixed ) {
        for ( i = startAt1; i < midAt && Ct1->iso_sort_key_Hfixed[i] == Ct2->iso_sort_key_Hfixed[i]; i ++ )
            ;
        if ( i < midAt ) {
            diff = Ct1->iso_sort_key_Hfixed[i] > Ct2->iso_sort_key_Hfixed[i]? 1:-1;
            goto done;
        }
    }
#endif

done:
#ifdef INCHI_CANON_MIN
    diff = -diff;
#endif

    if ( diff ) {
        diff = (diff > 0)? (nLayer+1) : -(nLayer+1); /* return the discovered difference layer number >= 1 */
        if ( kLeastForLayer ) {
#if ( bRELEASE_VERSION != 1 )
            if ( abs(kLeastForLayer[nLayer].k) > k+1 ) { /* for debug only */
                int stop = 1; /* <BRKPT> */
            }
#endif
            if ( !kLeastForLayer[nLayer].k ) {
                kLeastForLayer[nLayer].k = (diff > 0)? (k+1) : -(k+1);
                kLeastForLayer[nLayer].i = i;
            }
            if ( nLayer /* && !bOnlyCommon */) {
                diff = 0;
            }
        }
    }
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    else {
        int stop = 1;  /* for debug only */
    }
#endif
    if ( cmp ) {
        cmp[k] = (diff > 0)? 1 : (diff < 0)? -1 : 0;
    }
    return diff;
}
/**************************************************************************************************************/
int CtFullCompare( ConTable *Ct1, ConTable *Ct2, int bOnlyCommon, int bSplitTautCompare )
{
    int     startCt1, endCt1, startCt2, endCt2; /*endCt,*/
    int     startAt1, endAt1, startAt2, endAt2; /*endCt,*/
    int     midCt   /* end of atoms only in Ctbl */,
            midNumH = 0 /* end of atoms only NumH */,
            midAt   /* end of atoms only */;
    int     diff, i, k1, k2, lenNumH1, lenNumH2, lenNumH, maxVert /* min num atoms */;
    int     len_iso_sort_key1, len_iso_sort_key2, len_iso_sort_key /*, mid_iso_sort_key*/;
    int     nLayer = 0;

    k1 = Ct1->lenPos-1;
    k2 = Ct2->lenPos-1;

    diff = 0;

    startCt1 = startCt2 = 0;
    startAt1 = startAt2 = 0;

    endCt1 = Ct1->nextCtblPos[k1];
    endCt2 = Ct2->nextCtblPos[k2];
    endAt1 = (int)Ct1->nextAtRank[k1]-1;
    endAt2 = (int)Ct2->nextAtRank[k2]-1;

    maxVert = inchi_min(Ct1->maxVert, Ct2->maxVert);

    if ( bOnlyCommon ) {
        endCt1 = inchi_min(endCt1, endCt2);
        endCt1 = endCt2 = inchi_min(endCt1, Ct1->lenCt);
        endAt1 = endAt2 = inchi_min(endAt1, endAt2);
        if ( Ct1->Ctbl[endCt1] == EMPTY_CT || Ct1->Ctbl[endCt1] == 0 ||
             Ct2->Ctbl[endCt1] == EMPTY_CT || Ct2->Ctbl[endCt1] == 0 ) {
            endCt1 = endCt2 = endCt1-1;
        }
        lenNumH  =
        lenNumH1 =
        lenNumH2 = inchi_min(Ct1->lenNumH, Ct2->lenNumH);
        len_iso_sort_key  =
        len_iso_sort_key1 =
        len_iso_sort_key2 = inchi_min(Ct1->len_iso_sort_key, Ct1->len_iso_sort_key);
    } else {
        if ( Ct1->Ctbl[endCt1-1] == EMPTY_CT ) {
            endCt1 --;
        }
        if ( Ct2->Ctbl[endCt2-1] == EMPTY_CT ) {
            endCt2 --;
        }
        lenNumH1          = Ct1->lenNumH;
        lenNumH2          = Ct2->lenNumH;
        lenNumH           = inchi_min(lenNumH1, lenNumH2);
        len_iso_sort_key1 = Ct1->len_iso_sort_key;
        len_iso_sort_key2 = Ct2->len_iso_sort_key;
        len_iso_sort_key  = inchi_min(len_iso_sort_key1, len_iso_sort_key2);
    }

    if ( (diff = -(endCt1 - endCt2)) ) { /* negative sign reproduces results for NSC=28393 */
        goto done;
    }
    
    if ( bSplitTautCompare ) {
        midCt = inchi_min(Ct1->nLenCTAtOnly, Ct2->nLenCTAtOnly);
        if ( midCt > endCt1 ) {
            midCt = endCt1;
        }
        midAt = inchi_min(maxVert, endAt1);
    } else {
        midCt = endCt1;
        midAt = endAt1;
    }

    /*************************************************************************/
    /************ layer 0: connection table without tautomeric groups ********/
    /*************************************************************************/
    for ( i = startCt1; i < midCt && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i ++ )
        ;
    if ( i < midCt ) {
        diff = (int)Ct1->Ctbl[i] - (int)Ct2->Ctbl[i];
        goto done;
    }
    /*************************************************************************/
    /************* layer 1: H atoms without tautomeric H *********************/
    /*************************************************************************/
    nLayer ++;
    if ( Ct1->NumH && Ct2->NumH ) {
        if ( (diff = -(lenNumH1 - lenNumH2)) ) { /* negative sign reproduces results for NSC=28393 */
            goto done;
        }
        if ( endAt1 < maxVert ) {
            midNumH = lenNumH1 = endAt1;
        } else
        if ( bSplitTautCompare ) {
            midNumH = maxVert;
        } else {
            midNumH = lenNumH1;
        }
        for ( i = startAt1; i < midNumH && Ct1->NumH[i] == Ct2->NumH[i]; i ++ )
            ;
        if ( i < midNumH ) {
            diff = (int)Ct1->NumH[i] - (int)Ct2->NumH[i];
            goto done;
        }
    }
    /*************************************************************************/
    /************** layer 2: tautomeric part of CT and tautomeric H **********/
    /*************************************************************************/
    nLayer ++;
    for ( i = midCt; i < endCt1 && Ct1->Ctbl[i] == Ct2->Ctbl[i]; i ++ )
        ; /* compare tautomeric groups part of CT */
    if ( i < endCt1 ) {
        diff = (int)Ct1->Ctbl[i] - (int)Ct2->Ctbl[i];
        goto done;
    }
    if ( Ct1->NumH && Ct2->NumH ) {
        for ( i = midNumH; i < lenNumH1 && Ct1->NumH[i] == Ct2->NumH[i]; i ++ )
            ; /* compare tautomeric H */
        if ( i < lenNumH1 ) {
            diff = (int)Ct1->NumH[i] - (int)Ct2->NumH[i];
            goto done;
        }
    }
    /*************************************************************************/
    /************** layer 3: Fixed H atoms ***********************************/
    /*************************************************************************/
    nLayer ++;
    if ( Ct1->NumHfixed && Ct2->NumHfixed ) {
        for ( i = startAt1; i < endAt1 && Ct1->NumHfixed[i] == Ct2->NumHfixed[i]; i ++ )
            ;
        if ( i < endAt1 ) {
            diff = (int)Ct1->NumHfixed[i] - (int)Ct2->NumHfixed[i];
            goto done;
        }
    }
    /*************************************************************************/
    /************** layer 4: isotopic atoms, H and isotopic taut H ***********/
    /*************************************************************************/
    nLayer ++;
    if ( Ct1->iso_sort_key && Ct2->iso_sort_key ) {
        if ( (diff = -(len_iso_sort_key1 - len_iso_sort_key2)) ) { /* negative sign reproduces results for NSC=28393 */
            goto done;
        }
        for ( i = startAt1; i < endAt1 && Ct1->iso_sort_key[i] == Ct2->iso_sort_key[i]; i ++ )
            ;
        if ( i < endAt1 ) {
            diff = Ct1->iso_sort_key[i] > Ct2->iso_sort_key[i]? 1:-1;
            goto done;
        }
    }
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    /*************************************************************************/
    /************** layer 6: Fixed isotopic H atoms **************************/
    /*************************************************************************/
    nLayer ++;
    if ( Ct1->iso_sort_key_Hfixed && Ct2->iso_sort_key_Hfixed ) {
        for ( i = startAt1; i < midAt && Ct1->iso_sort_key_Hfixed[i] == Ct2->iso_sort_key_Hfixed[i]; i ++ )
            ;
        if ( i < midAt ) {
            diff = Ct1->iso_sort_key_Hfixed[i] > Ct2->iso_sort_key_Hfixed[i]? 1:-1;
            goto done;
        }
    }
#endif

done:
#ifdef INCHI_CANON_MIN
    diff = -diff;
#endif

    if ( diff ) {
        diff = (diff > 0)? (nLayer+1) : -(nLayer+1); /* return the discovered difference layer number >= 1 */
    }
    return diff;
}
/****************************************************************/
int CtFullCompareLayers( kLeast *kLeastForLayer )
{
    int iLayer;
    /* check for the rejection condition: Lambda > zb_rho_fix */
    for ( iLayer = 0; iLayer < MAX_LAYERS; iLayer ++ ) {
        if ( kLeastForLayer[iLayer].k ) {
            return (kLeastForLayer[iLayer].k > 0)? (iLayer+1) : -(iLayer+1);
        }
    }
    return 0;
}
/**************************************************************/
int CtCompareLayersGetFirstDiff( kLeast *kLeast_rho, int nOneAdditionalLayer,
                                 int *L_rho, int *I_rho, int *k_rho )
{
    int iLayer;
    if ( kLeast_rho ) {
        for ( iLayer = 0; iLayer < MAX_LAYERS; iLayer ++ ) {
            if ( kLeast_rho[iLayer].k ) {
                *L_rho = iLayer;
                *I_rho = kLeast_rho[iLayer].i;
                *k_rho = kLeast_rho[iLayer].k;
                break;
            }
        }
        if ( iLayer == MAX_LAYERS ) {
            if ( nOneAdditionalLayer ) {
                *L_rho = nOneAdditionalLayer;  /* ??? subtract 1 ??? */
                *I_rho = -1;
                *k_rho =  0;
                return 0; /* difference may be in the first additional layer */
            } else {
                *L_rho = INFINITY;
                *I_rho = -1;
                *k_rho =  0;
                return 0; /* no difference found */
            }
        } else {
            return 1; /* difference in a real layer */
        }
    } else {
        return -1; /* no input, should not happen */
    }
}
/**************************************************************/
int CtPartCompareLayers( kLeast *kLeast_rho, int L_rho_fix_prev, int nOneAdditionalLayer )
{
    int L_rho, I_rho, k_rho;
    if ( 0 < CtCompareLayersGetFirstDiff( kLeast_rho, nOneAdditionalLayer, &L_rho, &I_rho, &k_rho ) &&
         /* differences has been found in a real layer or all real layers are identical */
         L_rho <= L_rho_fix_prev ) {
         /* in this layer pzb_rho == pzb_rho_fix or in the previous real layer */
        return k_rho > 0? (L_rho+1) : -(L_rho+1);
    }
    return 0;
}
/****************************************************************/
void UpdateCompareLayers( kLeast kLeastForLayer[], int hzz )
{
    int i;
    if ( kLeastForLayer ) {
        for ( i = 0; i < MAX_LAYERS; i ++ ) {
            if ( abs(kLeastForLayer[i].k) >= hzz ) {
                kLeastForLayer[i].k = 0;
                kLeastForLayer[i].i = 0;
            }
        }                                          
    }
}
    
/****************************************************************/
void CtPartCopy( ConTable *Ct1 /* to */, ConTable *Ct2 /* from */, int k )
{
    int     startCt1, startCt2, endCt1, endCt2;
    int     len2, len2H, len2Hfixed, len2iso_sort_key, len2iso_exchg_atnos, i;
    int     startAt1, endAt1, startAt2, endAt2; /*endCt,*/
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    int     len2iso_sort_key_Hfixed;
#endif
    k --;
    if ( k ) {
        startCt1  = Ct1->nextCtblPos[k-1];
        startCt2  = Ct2->nextCtblPos[k-1];
        startAt1  = Ct1->nextAtRank[k-1]-1;
        startAt2  = Ct2->nextAtRank[k-1]-1;
    } else {
        startCt1 = startCt2 = 0;
        startAt1 = startAt2 = 0;
    }

    endCt1 = Ct1->nextCtblPos[k];
    endCt2 = Ct2->nextCtblPos[k];
    endAt1 = (int)Ct1->nextAtRank[k]-1;
    endAt2 = (int)Ct2->nextAtRank[k]-1;

    len2   = endCt2-startCt2;
    /* len    = min(len1, len2); */
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    if ( startCt1 != startCt2 || startAt1 != startAt2  ) {
        int stop = 1;
    }
#endif

    /* copy connection table: Ctbl */
    for ( i = 0; i < len2; i ++ ) {
        Ct1->Ctbl[startCt1+i] = Ct2->Ctbl[startCt2+i];
    }
    /* copy number of H: NumH */
    len2H = 0;
    if ( Ct1->NumH && Ct2->NumH ) {
        len2H   = endAt2-startAt2;
        if ( endAt2 > Ct2->maxVert ) {
            len2H = Ct2->lenNumH - startAt2;
        }
        for ( i = 0; i < len2H; i ++ ) {
            Ct1->NumH[startAt1+i] = Ct2->NumH[startAt2+i];
        }
    }
    /* copy number of fixed H */
    len2Hfixed = 0;
    if ( Ct1->NumHfixed && Ct2->NumHfixed ) {
        len2Hfixed   = endAt2-startAt2;
        for ( i = 0; i < len2Hfixed; i ++ ) {
            Ct1->NumHfixed[startAt1+i] = Ct2->NumHfixed[startAt2+i];
        }
    }
    /* copy isotopic keys */
    len2iso_sort_key = 0;
    if ( Ct1->iso_sort_key && Ct2->iso_sort_key ) {
        len2iso_sort_key = endAt2-startAt2;
        for ( i = 0; i < len2iso_sort_key; i ++ ) {
            Ct1->iso_sort_key[startAt1+i] = Ct2->iso_sort_key[startAt2+i];
        }
    }
    len2iso_exchg_atnos = 0;
    if ( Ct1->iso_exchg_atnos && Ct2->iso_exchg_atnos ) {
        len2iso_exchg_atnos = endAt2-startAt2;
        for ( i = 0; i < len2iso_exchg_atnos; i ++ ) {
            Ct1->iso_exchg_atnos[startAt1+i] = Ct2->iso_exchg_atnos[startAt2+i];
        }
    }
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    len2iso_sort_key_Hfixed = 0;
    if ( Ct1->iso_sort_key_Hfixed && Ct2->iso_sort_key_Hfixed ) {
        len2iso_sort_key_Hfixed = endAt2-startAt2;
        for ( i = 0; i < len2iso_sort_key; i ++ ) {
            Ct1->iso_sort_key_Hfixed[startAt1+i] = Ct2->iso_sort_key_Hfixed[startAt2+i];
        }
    }
#endif
    Ct1->lenCt          = startCt1 + len2;
    Ct1->nextCtblPos[k] = startCt1 + len2;
    Ct1->nextAtRank[k]  = Ct2->nextAtRank[k];
    if ( len2H ) {
        Ct1->lenNumH        = startAt1 + len2H;
    }
    /*
    if ( len2Hfixed ) {
        Ct1->lenNumHfixed   = startAt1 + len2Hfixed;
    }
    */
    if ( len2iso_sort_key ) {
        Ct1->len_iso_sort_key = startAt1 + len2iso_sort_key;
    }
    if ( len2iso_exchg_atnos ) {
        Ct1->len_iso_exchg_atnos = startAt1 + len2iso_exchg_atnos;
    }

#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    if ( len2iso_sort_key_Hfixed ) {
        Ct1->len_iso_sort_key_Hfixed = startAt1 + len2iso_sort_key_Hfixed;
    }
#endif
#ifdef INCHI_CANON_USE_HASH        
    Ct1->hash[k]        = Ct2->hash[k];
#endif
    Ct1->lenPos         = k+1;
    INCHI_HEAPCHK
}
/****************************************************************/
void CtFullCopy( ConTable *Ct1, ConTable *Ct2 )
{
    /* Ct1 does not have INFINITY termination */
    int k;
    for ( k = 0; k < Ct2->lenPos; k ++ ) {
        CtPartCopy( Ct1 /* to */, Ct2 /* from */, k+1 );
    }
}
/****************************************************************/
void TranspositionGetMcrAndFixSetAndUnorderedPartition( Transposition *gamma, NodeSet *McrSet, NodeSet *FixSet, int n, int l, UnorderedPartition *p )
{
    int i, j, k, mcr, num;
    AT_RANK next;
    bitWord *McrBits = McrSet->bitword[l-1];
    bitWord *FixBits = FixSet->bitword[l-1];
    int     len      = McrSet->len_set*sizeof(bitWord);

    memset( McrBits, 0, len );
    memset( FixBits, 0, len );
    for ( i = 0; i < n; i ++ ) {
        p->equ2[i] = INFINITY; /* for debug only */
    }

    for ( i = 0; i < n; i ++ ) {
        j = (int)(next = gamma->nAtNumb[i]);
        if ( j == i ) {
            FixBits[ i / num_bit ] |= bBit[ i % num_bit ];
            McrBits[ i / num_bit ] |= bBit[ i % num_bit ];
            /* p->next[i] = INFINITY; */ /* no link to same orbit points */
            p->equ2[i] = next;  /* fixed point */
        } else
        if ( !(rank_mark_bit & next) ) {
            gamma->nAtNumb[i] |= rank_mark_bit;
            mcr = inchi_min(j, i);
            num = 0;
            /* mark all nodes in the cycle to ignore later; find mcr */
            while( !(rank_mark_bit & (next = gamma->nAtNumb[j])) ) {
                gamma->nAtNumb[j] |= rank_mark_bit;
                j = (int)next;
                if ( mcr > j ) {
                    mcr = j;
                }
                num ++;
            }
            McrBits[ mcr / num_bit ] |= bBit[mcr % num_bit]; /* save mcr */
            /* fill out the unordered partition, the mcr first, other in the cycle after that */
            p->equ2[mcr] = mcr;
            for ( k = mcr; mcr != (j = (int)(rank_mask_bit & gamma->nAtNumb[k])); k = j ) {
                p->equ2[j] = mcr;
            }
        }
    }
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    /* for debug only */
    for ( i = 0; i < n; i ++ ) {
        if ( p->equ2[i] >= n ) {
            int stop = 1;
        }
    }
#endif
    /* remove the marks */
    for ( i = 0; i < n; i ++ ) {
        gamma->nAtNumb[i] &= rank_mask_bit;
    }
    INCHI_HEAPCHK
}
/****************************************************************/
int SetBitCreate( void )
{
    bitWord  b1, b2;
    AT_NUMB n1, n2;
#ifdef INCHI_CANON_USE_HASH        
    CtHash  h1, h2;
#endif
    int    i;

    if ( bBit ) {
        INCHI_HEAPCHK
        return 0; /* already created */
    }

    b1 = 1;
    num_bit  = 1;
    for ( b1=1, num_bit=1; b1 < (b2 = (bitWord)((b1 << 1)& BIT_WORD_MASK)); b1 = b2, num_bit ++ )
        ;
    bBit = (bitWord*)inchi_calloc( num_bit, sizeof(bitWord));
    if ( !bBit ) {
        INCHI_HEAPCHK
        return -1; /* failed */
    }
    for ( i = 0, b1=1; i < num_bit; i++, b1 <<= 1 ) {
        bBit[i] = b1;
    }

    for ( n1 = 1; n1 < (n2 = (AT_RANK)((n1 << 1)& AT_RANK_MASK) ); n1 = n2 )
        ;
    rank_mark_bit = n1;
    rank_mask_bit = ~n1;

#ifdef INCHI_CANON_USE_HASH        
    for ( h1 = 1; h1 < (h2 = (h1 << 1)); h1 = h2 )
        ;
    hash_mark_bit = h1;
#endif
    INCHI_HEAPCHK
    return 1;
}
/****************************************************************/
int SetBitFree( void )
{
    if ( bBit ) {
        inchi_free( bBit );
        bBit = NULL;
        INCHI_HEAPCHK
        return 1; /* success */
    }
    INCHI_HEAPCHK
    return 0; /* already destroyed */
}
#ifdef NEVER  /* { how to renumber a graph */
/*********************************************************************/
void RenumberTheGraph( int n, NEIGH_LIST *NeighList, AT_NUMB *old2new,
                       AT_NUMB *new2old, S_CHAR *mark, int bDone )
{
    int        i, k, j;
    NEIGH_LIST nl;
    
    /* renumber neighbors */
    for ( i = 0; i < n; i ++ ) {
        for ( j = 1; j <= NeighList[i][0]; j ++ ) {
            NeighList[i][j] = old2new[NeighList[i][j]];
        }
    }
    /* rearrange NeighList in situ using new2old[] */
    for ( k = 0; k < n; k ++ ) {
        if ( mark[k] & bDone )
            continue;
        if ( k == ( j = new2old[k] ) ) {
            mark[k] |= bDone;
            continue;
        }
        /* transposition cycle */
        i = k;
        nl = NeighList[k];
        do {
            NeighList[i] = NeighList[j];
            mark[i] |= bDone;
            i = j;
        } while ( k != ( j = new2old[i] ) );
        NeighList[i] = nl;
        mark[i] |= bDone;
    }
#ifdef NEVER
    /* rearrange NeighList in situ using old2new[] */
    s = 0;
    for ( k = 0; k < n; k ++ ) {
        if ( mark[k] & bDone )
            continue;
        if ( k == ( j = old2new[k] ) ) {
            mark[k] |= bDone;
            continue;
        }
        /* transposition cycle */
        i = k;
        /* NeighList[j] goes to ith position */
        nl2[s] = NeighList[j];
        s ^= 1;
        do {
            nl2[s] = NeighList[i];
            NeighList[i] = nl2[s ^= 1];
            mark[i] |= bDone;
            i = j;
            j = old2new[i];
        } while ( k != ( j = old2new[i] ) );
        NeighList[j] = nl2[s ^= 1];
        mark[j] |= bDone;
    }
#endif
}
/***************************************************************************************/
void RearrangeAtRankArray ( int n, AT_RANK *nRank, AT_NUMB *new2old, S_CHAR *mark, int bDone )
{
    int i, k, j;
    AT_RANK r;    
    /* rearrange the array in situ using new2old[] */
    for ( k = 0; k < n; k ++ ) {
        if ( mark[k] & bDone )
            continue;
        if ( k == ( j = new2old[k] ) ) {
            mark[k] |= bDone;
            continue;
        }
        /* transposition cycle */
        i = k;
        r = nRank[k];
        do {
            nRank[i] = nRank[j];
            mark[i] |= bDone;
            i = j;
        } while ( k != ( j = new2old[i] ) );
        nRank[i] = r;
        mark[i] |= bDone;
    }

}
/***************************************************************************************/
void RenumberAtNumbArray( int n, AT_NUMB *nAtNumb, AT_NUMB *old2new )
{
    int i;
    for ( i = 0; i < n; i ++ ) {
        nAtNumb[i] = old2new[nAtNumb[i]];
    }
}
/****************************************************************/
int  GetCanonRanking2( int num_atoms, int num_at_tg, int num_max, int bDigraph, sp_ATOM* at,
                     AT_RANK **pRankStack,  int nNumPrevRanks,
                     AT_RANK *nSymmRank,  AT_RANK *nCanonRank,
                     NEIGH_LIST *NeighList, AT_RANK *nTempRank,
                     CANON_STAT* pCS )
{
    void *pzb_rho=NULL;
    int ret, cmp1=0, cmp2=0;
    
    int i, j, k, n;
    AT_NUMB *old2new = NULL, *new2old = NULL, m, r1, r2;
    S_CHAR  *mark = NULL;

    int      nMaxLenCt        =  pCS->nMaxLenLinearCT;
    AT_RANK *pCt              =  pCS->LinearCT;
    int      nLenCt           =  pCS->nLenLinearCT;
    AT_RANK *pCt0             =  pCS->LinearCT;
    int      nLenCt0          =  pCS->nLenLinearCT;


    CANON_DATA CanonData;
    CANON_DATA *pCD = &CanonData;
    CANON_COUNTS CanonCounts;
    CANON_COUNTS *pCC = &CanonCounts;
    memset (pCD, 0, sizeof(pCD[0]));
    memset (pCC, 0, sizeof(pCC[0]));
    /* pointers */
    pCD->LinearCT           = pCS->LinearCT;
    /* variables - unchanged */
    pCD->ulTimeOutTime      = pCS->ulTimeOutTime;
    pCD->nMaxLenLinearCT    = pCS->nMaxLenLinearCT;
    /* return values & input/output */
    pCD->nLenLinearCT       = pCS->nLenLinearCT;

    pCC->lNumBreakTies      = pCS->lNumBreakTies;
    pCC->lNumDecreasedCT    = pCS->lNumDecreasedCT;
    pCC->lNumRejectedCT     = pCS->lNumRejectedCT;
    pCC->lNumEqualCT        = pCS->lNumEqualCT;
    pCC->lNumTotCT          = pCS->lNumTotCT;

    ret = CanonGraph( num_atoms, num_at_tg, num_max, bDigraph, NeighList, (Partition *)pRankStack,
                       nSymmRank,  nCanonRank, pCS->nPrevAtomNumber, pCD, pCC, NULL, &pzb_rho );

    pCS->nLenLinearCT       = pCD->nLenLinearCT;
    pCS->lNumBreakTies      = pCC->lNumBreakTies;
    pCS->lNumDecreasedCT    = pCC->lNumDecreasedCT;
    pCS->lNumRejectedCT     = pCC->lNumRejectedCT;
    pCS->lNumEqualCT        = pCC->lNumEqualCT;
    pCS->lNumTotCT          = pCC->lNumTotCT;

   
    /* save the connection table for comparison with the 2nd one */
    pCt0 = (AT_RANK*)inchi_calloc(nMaxLenCt, sizeof(pCt0[0]));
    memcpy(pCt0, pCS->LinearCT, nMaxLenCt*sizeof(pCt0[0]));
    nLenCt0 =  pCS->nLenLinearCT;

    /**********************************************************/
    /* rearrange numbering to make canon. numbering the first */
    /**********************************************************/
    n = num_at_tg;
    /* 1. get transpositions */
    old2new = (AT_NUMB*) inchi_calloc( n, sizeof(old2new[0]) ); 
    new2old = (AT_NUMB*) inchi_calloc( n, sizeof(new2old[0]) );
    mark    = (S_CHAR *) inchi_calloc( n, sizeof(mark[0]) );
    for ( i = 0; i < n; i ++ ) {
        /* forward transposition: at[i] -> at[old2new[i]] position */
        old2new[i] = m = nCanonRank[i]-1;
        /* forward transposition: at[new2old[i]] -> at[i] position */
        new2old[m] = i;
    }
    /* rearrange input data according to the new numbering */
    RenumberTheGraph( n, NeighList, old2new, new2old, mark, 1 );
    RearrangeAtRankArray ( n, pRankStack[0], new2old, mark, 2 );
    RenumberAtNumbArray ( n, pRankStack[1], old2new );
    /* make sure the atom numbers are sorted */
    for ( i = k = 0, r1 = pRankStack[0][pRankStack[1][i]]; i < n; r1 = r2) {
        for ( j = i++; i < n && r1 == (r2 = pRankStack[0][pRankStack[1][i]]); i ++ )
            ;
        if ( i - j > 1 ) {
            k += insertions_sort_AT_RANK( pRankStack[1]+j, i-j );
        }
    }

    ret = CanonGraph( num_atoms, num_at_tg, num_max, bDigraph, NeighList, (Partition *)pRankStack,
                       nSymmRank,  nCanonRank, pCS->nPrevAtomNumber, pCD, pCC, &pzb_rho, NULL );

    pCS->nLenLinearCT       = pCD->nLenLinearCT;
    pCS->lNumBreakTies      = pCC->lNumBreakTies;
    pCS->lNumDecreasedCT    = pCC->lNumDecreasedCT;
    pCS->lNumRejectedCT     = pCC->lNumRejectedCT;
    pCS->lNumEqualCT        = pCC->lNumEqualCT;
    pCS->lNumTotCT          = pCC->lNumTotCT;


    /* compare the connection tables */
    cmp1 = nLenCt0 - pCS->nLenLinearCT;
    cmp2 = memcmp( pCt0, pCS->LinearCT, pCS->nLenLinearCT*sizeof(pCt0[0]));
#ifdef _DEBUG            
    if ( cmp1 || cmp2 ) {
        int stop = 1;
    }
#endif
    /**********************************************************/
    /* rearrange numbering back to the original numbering     */
    /**********************************************************/
    /* restore the input data to its original numbering */
    RenumberTheGraph( n, NeighList, new2old, old2new, mark, 4 );
    RearrangeAtRankArray ( n, pRankStack[0], old2new, mark, 8 );
    RenumberAtNumbArray ( n, pRankStack[1], new2old );
    /* rearrange the output data to the original numbering */
    RearrangeAtRankArray ( n, nCanonRank, old2new, mark, 16 );
    RenumberAtNumbArray ( n, pCS->nPrevAtomNumber, new2old );
    RearrangeAtRankArray ( n, nSymmRank, old2new, mark, 32 );

    /* free memory */
    CTableFree( pzb_rho );
    if ( pzb_rho ) {
        inchi_free( pzb_rho );
    }
    inchi_free( old2new );
    inchi_free( new2old );
    inchi_free( mark );
    inchi_free( pCt0 );

    return ret;
}
#endif  /* } */

#define QZFIX_OK(X) ((X)<=0)

int GetOneAdditionalLayer( CANON_DATA *pCD, ConTable *pzb_rho_fix )
{
    int nLastLayer = -1, nNumLast = 0, nLayer = 0;

    if ( !pCD || !pzb_rho_fix ) {
        return 0;
    }

    nLayer ++; /* 1 */
    if ( pCD->NumH && !pzb_rho_fix->NumH ) {
        nLastLayer = nLayer;
        nNumLast ++;
    }
    nLayer ++; /* 2 */
    if ( pCD->nLenCTAtOnly < pCD->nLenLinearCT && pzb_rho_fix->nLenCTAtOnly == pzb_rho_fix->lenCt ) {
        nLastLayer = nLayer;
        nNumLast ++;
    }
    nLayer ++; /* 3 */
    if ( pCD->NumHfixed && !pzb_rho_fix->NumHfixed ) {
        nLastLayer = nLayer;
        nNumLast ++;
    }
    nLayer ++; /* 4 */
    if ( pCD->iso_sort_key && !pzb_rho_fix->iso_sort_key ) {
        nLastLayer = nLayer;
        nNumLast ++;
    }
    /*
    nLayer ++; // 5
    if ( pCD->nLenCTAtOnly < pCD->nLenLinearCT && pCD->iso_sort_key &&
        (pzb_rho_fix->nLenCTAtOnly == pzb_rho_fix->lenCt || !pzb_rho_fix->iso_sort_key ) ) {
        nLastLayer = nLayer;
        nNumLast ++;
    }
    */
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    nLayer ++; /* 6 */
    if ( pCD->iso_sort_key_Hfixed && !pzb_rho_fix->iso_sort_key_Hfixed ) {
        nLastLayer = nLayer;
        nNumLast ++;
    }
#endif
    if ( 1 == nNumLast ) {
        return nLastLayer;
    }
    return 0;
}


/*#define QZFIX_OK(X) (!(X))*/

/********************** CanonGraph *************************************
 *   A naive implementation of graph canonical numbering algorithm     *
 *   from "Practical Graph Isomorphism" by Brendan D. McKay,           *
 *   Congressus Numerantium, Vol. 30 (1981), pp. 45 - 87.              *
 *   Note: Several typos fixed, added chem. struct. specifics          *
 ***********************************************************************/

/* on entry: pi[0] is equitable */
/*******************************************************10/21/2003******
 * Later add optimization: if Aut(G) <= Aut(G0) due to some additional *
 * layer of coloring applied to G and the following is known about G0: *                                                    *
 *                                                                     *
 * 0) canonical numbering of G should be same as that of G0            *
 *                                                                     *
 * 1) canonical numbering as v= v0[n] {vertex number v from            *
 *                                     G0 canonical number n)          *
 * 2) orbits of Aut(G0) as UnorderedPartition theta0                   *
 *                                                                     *
 * then when choosing next v[i] for refining the partition consider    *
 * only vertices from the Aut(G0) orbit of V(i), that is, only such    *
 * v[i] that:                                                          *
 *                                                                     *
 *       GetUnorderedPartitionMcrNode( &theta0, v[i] ) ==              *
 *       GetUnorderedPartitionMcrNode( &theta0, v0[i] )                *
 ***********************************************************************/
int CanonGraph( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[],
                AT_RANK *nSymmRank,  AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon,
                CANON_DATA *pCD, CANON_COUNTS *pCC,
                ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out  )
{   
    /* bDigraph != 0
       means consider edges from atoms to t-groups
       as directed, that is, do not include
       t-group ranks in comparing neighbors
       when refining partition
    */

    /* Always set
       lab = true
       dig = true
    */

    /* in the comments:
        m = |zeta|
        r = |rho|

        m < n or r < n in case pi[k] in P (i.e. satisfies Lemma 2.25)

        Just after passing point B:
        ===========================
        K = k-1
        wi = v[i], i = 1..K
        Gamma(0) = Gamma = Aut(G)pi
        Gamma(i) = Gamma(w1,w2,...,wi) pointwise stabilizer for i=1..K
        zeta is a terminal node =>
          the coarsest equitable partition that fixes w1,...,wK is discrete =>
            Gamma(K)=1
        At point A only:
            index = |Gamma(k-1)|/|Gamma(k)|
        At points A and B:
            size  = |Gamma(k-1)|
            theta = theta(Gamma(k-1));
            Gamma(k-1) = <Y>, where Y is the set of all automprhisms output up
              to the present stage (in Step 10 = L10 )
            |Y| <= n - |theta|
    */

    AT_RANK    *pCt              =  pCD->LinearCT;
    /*int         nMaxLenCt        =  pCD->nMaxLenLinearCT;*/
    int        *nLenCt           = &pCD->nLenLinearCT;
    CANON_DATA *pCD1             = pCD;

    int i, k, k2, index, l, ok, ret=0, res;
    int t_Lemma;   /* hh: if pi[k-1] satisfies Lemma 2-25 then 
                          t_Lemma = min{i| i=1..k && pi[i-1] satisfies Lemma 2-25}*/
                   /* otherwise t_Lemma = k --> here this is always the case */
    int t_eq_zeta; /* ht: min{i|i=1..m && all terminal modes descended from or equal
                      to zeta(i) have been shown to be equivalent}. */
    int h_zeta;    /* h: the longest common ancestor of zeta and nu is nu(h_zeta) */
    int h_rho;     /* hb: the longest common ancestor of rho and nu is nu(h_rho) */
    int hz_rho;    /* hzb: max{i|i=1..min(k,r) && Lambda(G,pi,nu(i)) == Lambda(G,pi,rho(i))} */
    int hz_zeta;   /* hzf: max{i|i=1..min(k,m) && Lambda(G,pi,nu(i)) == Lambda(G,pi,zeta(i))} */
    int qzb_rho;   /* Ct(Lambda[k]) - Ct(rho[k]) */
    double size;     /* |Aut(G)| */
    int  nNumLayers = (NULL != pCD->NumH) + (NULL != pCD->NumHfixed) +
                      /* (bDigraph && pCD->nLenLinearCT > pCD->nLenCTAtOnly)*/ /* ??? tautomeric */
                      (NULL != pCD->iso_sort_key)
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
                      + (NULL != pCD->iso_sort_key_Hfixed)
#endif
                      ;
    int  dig        = (bDigraph || nNumLayers);
    int  bSplitTautCompare = (bDigraph || nNumLayers); /* compare taut. H and tgroups connections after H */
                   /* digraph: 1=>do not use Lemma 2.25, 0 => use */
    int  lab = 1;  /* label: 1=>find canonical numbering;
                      0=>do not find canonical numbering, do not use rho */
    int  r;        /* |rho| */
    int  bZetaEqRho = lab;
    int  bZetaIsomorph;

    long lNumEqlZeta;

    
    const int L = MAX_SET_SIZE;
    UnorderedPartition theta, theta_from_gamma;
    Cell *W;  /* W[i] is the first non-trivial cell of pi[i+1] */ 
    Node *v;  /* v[i] is in W[i] to create T(G,pi,nu[i+1]) */
    Node tvc, tvh;
    S_CHAR *e, *qzb=NULL;   /* qzb = NULL always */
    /* current node CT */
    ConTable Lambda;
    /* first leaf CT */
    ConTable zf_zeta; /* Ct for zeta,  the first discovered terminal node */
    /* best leaf/node CT: find the greatest pzb_rho possibly subject to pzb_rho[k] <= pzb_rho_fix[k] condition */
    ConTable *pzb_rho = NULL;  /* Ct for rho,   the best discovered terminal node */
    /* fixed input CT: for all k pzb_rho[k] <= pzb_rho_fix[k]; at the end pzb_rho == pzb_rho_fix */
    ConTable *pzb_rho_fix = (pp_zb_rho_inp && *pp_zb_rho_inp)? *pp_zb_rho_inp:NULL;
    
    NodeSet Omega; /* MAX_SET_SIZE */
    NodeSet Phi;   /* MAX_SET_SIZE */
    NodeSet cur_nodes;   /* 1 each */
    Transposition gamma;
    Partition zeta;      /* the first discovered terminal node */
    Partition rho;       /* the best discovered terminal node */
    int nNumFoundGenerators=0;
    int qzb_rho_fix = 0;
    int hzb_rho_fix = 0;
    int bRhoIsDiscrete = 1;
    kLeast kLeast_rho[MAX_LAYERS];
    kLeast kLeast_rho_fix[MAX_LAYERS];
    int nOneAdditionalLayer;
    int pzb_rho_fix_reached = 0;
    int L_rho_fix_prev = 0, I_rho_fix_prev=-1, k_rho_fix_prev=0;

    /* Note: Layered comparison should be consistent, especially in layer numbers.
             Layered comparison is implemented in:
                     CtFullCompare()
                     CtPartCompare()
                     GetOneAdditionalLayer()

             The partial comparison results in kLeast[] are used in
                     CtFullCompareLayers()
                     CtPartCompareLayers()
                     CtCompareLayersGetFirstDiff()
                     UpdateCompareLayers()
     */
    nOneAdditionalLayer = GetOneAdditionalLayer( pCD1, pzb_rho_fix );


    /* next 2 lines for debug only */
    /* num_g++;  */
    /* WriteGraph( G, n_tg, num_g, "V:\\IChI_v10\\Gordon-Graphs\\hard\\k06g08v312-alt.dre", "a+" ); */
    
    /* memory allocation */

    if ( 0 > SetBitCreate() ) {
        return -1;
    }
    if ( pzb_rho_fix && pzb_rho_fix->nLenCTAtOnly != pCD->nLenCTAtOnly ) {
        /* consistency check */
        return -2;
    }
    ok = 1;


    ok &= UnorderedPartitionCreate( &theta, n_tg );
    ok &= UnorderedPartitionCreate( &theta_from_gamma, n_tg );

    ok &= (NULL != (W   = (Cell*)inchi_calloc(n_tg, sizeof(W[0]))));
    ok &= (NULL != (v   = (Node*)inchi_calloc(n_tg, sizeof(v[0]))));
    ok &= (NULL != (e   = (S_CHAR*)inchi_calloc(n_tg, sizeof(e[0]))));

/*
    ok &= (NULL != (v   = (Node*)inchi_calloc(n_tg, sizeof(W[0]))));
    ok &= (NULL != (e   = (S_CHAR*)inchi_calloc(n_tg, sizeof(W[0]))));
*/

    /*    ok &= (NULL != (qzb = (S_CHAR*)calloc(n_tg, sizeof(W[0])))); */
    ok &= CTableCreate( &Lambda, n, pCD );
    ok &= CTableCreate( &zf_zeta, n, pCD );
    ok &= ( (pzb_rho = (ConTable *)inchi_calloc( 1, sizeof( *pzb_rho ) ) ) &&
            CTableCreate( pzb_rho, n, pCD ) );

    ok &= NodeSetCreate( &Omega, n_tg, MAX_SET_SIZE );
    ok &= NodeSetCreate( &Phi, n_tg, MAX_SET_SIZE );
    ok &= NodeSetCreate( &cur_nodes, n_tg, 1 );

    ok &= PartitionCreate( &zeta, n_tg);
    ok &= PartitionCreate( &rho, n_tg);
    ok &= TranspositionCreate( &gamma, n_tg );

    INCHI_HEAPCHK

    
/*L1:*/
    k = 1;
    size = 1.0;
    h_zeta = hz_rho = index = l = 0;

    if ( !ok ) {
        goto exit_function; /* initialization failed */
    }

    UnorderedPartitionMakeDiscrete(&theta, n_tg);
    t_Lemma = 2;

    pCC->lNumBreakTies   = 0;
    pCC->lNumDecreasedCT = 0;
    pCC->lNumRejectedCT  = 0;
    pCC->lNumEqualCT     = 1;
    pCC->lNumTotCT       = 0;
    lNumEqlZeta          = 1;
    
    hzb_rho_fix          = 1;

    memset( kLeast_rho, 0, sizeof(kLeast_rho) );
    memset( kLeast_rho_fix, 0, sizeof(kLeast_rho_fix) );

    if ( PartitionIsDiscrete( &pi[k-1], n_tg ) ) {
        /* added the following 3 lines to the original to create Ct */
        PartitionCopy( &rho, &pi[k-1], n_tg );
        CtPartFill( G, pCD, &pi[k-1], pzb_rho, 1, n, n_tg );
        CtPartInfinity( pzb_rho, qzb, 2 );
        pCC->lNumTotCT ++;
        r = k;
        /* goto L18; */
        goto exit_function;
    }
    if ( !dig && PartitionSatisfiesLemma_2_25( &pi[0], n ) )
        t_Lemma = 1;
    /*
    PartitionGetFirstCell( &pi[k-1], &W[k-1], k, n );
    v[k-1] = CellGetMinNode( &pi[k-1], &W[k-1], 0, pCD1 );
    CtPartClear( &Lambda, 1 ); 
    e[k-1] = 0;
    */
    CtPartClear( &Lambda, 1 ); 
    INCHI_HEAPCHK

/* L2: reach the first leaf and save it in zeta and rho */
    while( k ) {
        /* the two next lines intentionally switched */
        /* Create equitable partition in pi[k]  */
        PartitionGetFirstCell( &pi[k-1], W, k, n );
        v[k-1] = CellGetMinNode( &pi[k-1], &W[k-1], 0, pCD1 );
        e[k-1] = 0;
        if ( dig || !PartitionSatisfiesLemma_2_25(&pi[k-1], n) )
            t_Lemma = k+1;
        /* e[k-1] = 0; */
        { Node vv = v[k-1];
        if ( 0 > (ret=PartitionColorVertex( G, &pi[k-1], vv /*v[k-1]*/, n, n_tg, n_max, bDigraph, 0 )) ) {
            goto exit_error;
        }}
        pCC->lNumBreakTies ++;
        k ++;
        CtPartFill( G, pCD, &pi[k-1], &Lambda, k-1, n, n_tg );
        /* return -1; *//* debug only */
        /* if(h_zeta==0)goto L5; L5: */
        /* the first terminal node has not been reached yet */
        /* search for the predefined numbering */
        if ( pzb_rho_fix && QZFIX_OK(qzb_rho_fix) ) {
            qzb_rho_fix = CtPartCompare( &Lambda, pzb_rho_fix, qzb, kLeast_rho_fix, k-1, 1, bSplitTautCompare );
            if ( QZFIX_OK(qzb_rho_fix) ) {
                hzb_rho_fix = k;
            }
        }

        if ( lab && QZFIX_OK(qzb_rho_fix) )  /* DCh */
            CtPartCopy( pzb_rho, &Lambda, k-1 );
        CtPartCopy( &zf_zeta, &Lambda, k-1 );
        /*goto L4; L4:*/
        if ( PartitionIsDiscrete( &pi[k-1], n ) ) {
            break;  /* goto L7; */
        }
        /* goto L2; */
    }
    pCC->lNumTotCT ++;
    /* L7; L7: */
    /* if ( h_zeta == 0 ) goto L18; L18:*/
    h_zeta = t_eq_zeta = hz_zeta = k;
    CtPartInfinity( &zf_zeta, NULL, k );
    /******************** <<<===== B **************************/
    PartitionCopy( &zeta, &pi[k-1], n_tg );
    if ( lab ) {
        if ( pzb_rho_fix ) {
            if ( 0 == qzb_rho_fix ) {
                qzb_rho_fix = CtFullCompare( &Lambda, pzb_rho_fix, 1, bSplitTautCompare );
                if ( qzb_rho_fix > 0 ) {
                    hzb_rho_fix = 1;
                }
            }
            if ( hzb_rho_fix > 1 ) {
                PartitionCopy( &rho, &pi[hzb_rho_fix-1], n_tg );
                /*CtPartInfinity( pzb_rho, qzb, k );*/
            }
            hz_rho = h_rho = hzb_rho_fix;
            bRhoIsDiscrete = (hzb_rho_fix == k);
            if ( bRhoIsDiscrete ) {
                CtPartInfinity( pzb_rho, qzb, k );
                pzb_rho_fix_reached = !qzb_rho_fix;
                CtCompareLayersGetFirstDiff( kLeast_rho_fix, nOneAdditionalLayer,
                                 &L_rho_fix_prev, &I_rho_fix_prev, &k_rho_fix_prev );
            }
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
            else {
                int stop = 1;
            }
#endif
        } else {
            PartitionCopy( &rho, &pi[k-1], n_tg );
            hz_rho = h_rho = k;
            CtPartInfinity( pzb_rho, qzb, k );
        }
        qzb_rho = 0;
    }
    r = k;
    v[k-1] = INFINITY;     /* DCh */
    CellMakeEmpty( W, k ); /* DCh */
    k --;
    goto L13;


L2:
    /* the two next lines intentionally switched */
    /* Create equitable partition in pi[k]  */
    if ( 0 > (ret=PartitionColorVertex( G, &pi[k-1], v[k-1], n, n_tg, n_max, bDigraph, 0 )) ) {
        goto exit_error;
    }
    pCC->lNumBreakTies ++;
    k ++;
    CtPartFill( G, pCD, &pi[k-1], &Lambda, k-1, n, n_tg );
    e[k-1] = 0;         /* moved  */
    v[k-1] = INFINITY;  /* added by DCh. */
    CellMakeEmpty( W, k ); /* DCh */

    if ( hz_zeta == k-1 && 0 == CtPartCompare( &Lambda, &zf_zeta, NULL, NULL, k-1, 0, bSplitTautCompare ) ) {
        hz_zeta = k; /* max{k|Lambda(G,pi,nu(k))==Lambda(G,pi,zeta) }  */
    } /* added */
    
    /* -- old code ---
    if ( pzb_rho_fix && QZFIX_OK(qzb_rho_fix) ) {
        qzb_rho_fix = CtPartCompare( &Lambda, pzb_rho_fix, qzb, kLeast_rho_fix, k-1, 1, bSplitTautCompare );
        if ( QZFIX_OK(qzb_rho_fix) ) {
            hzb_rho_fix = k;
        } else {
            pCC->lNumRejectedCT ++;
        }
    }
    */

    /* --- new code ---*/
    if ( pzb_rho_fix && !qzb_rho_fix ) {
        qzb_rho_fix = CtPartCompare( &Lambda, pzb_rho_fix, qzb, kLeast_rho_fix, k-1, 1, bSplitTautCompare );
        if ( !qzb_rho_fix && bRhoIsDiscrete ) {
            qzb_rho_fix = CtPartCompareLayers( kLeast_rho_fix, L_rho_fix_prev, nOneAdditionalLayer );

#if ( FIX_ChCh_CONSTIT_CANON_BUG == 1 )
            if ( qzb_rho_fix ) {
                int L_rho_fix_diff = abs(qzb_rho_fix)-1;
                if ( L_rho_fix_diff < L_rho_fix_prev ||
                     ((L_rho_fix_diff == L_rho_fix_prev) &&
                     (kLeast_rho_fix[L_rho_fix_diff].i < I_rho_fix_prev)) ) {
                    qzb_rho_fix = L_rho_fix_diff+1; /* positive difference will be rejected */
                }
            }
#endif

#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
            if ( qzb_rho_fix ) {
                int stop = 1; /* debug only */
            }
#endif
        }
        if ( !QZFIX_OK(qzb_rho_fix) ) {
            pCC->lNumRejectedCT ++;
        }
    }
    if ( pzb_rho_fix && QZFIX_OK(qzb_rho_fix) ) {
        hzb_rho_fix = k;
    }
    /* if (!lab) goto L3; */
    if ( lab && QZFIX_OK(qzb_rho_fix) ) {
        /* once the difference has been found it is meaningful as long as k increments */
        /* cur_qzb2 = CtPartCompare( &Lambda, pzb_rho, qzb, k-1 ); */ /* rho compare */
        if ( hz_rho == k-1 && !qzb_rho && bRhoIsDiscrete ) {
            int qzb_rho_temp = 0;
            qzb_rho = CtPartCompare( &Lambda, pzb_rho, qzb, kLeast_rho, k-1, 0, bSplitTautCompare );
            /* old code */
            if ( !qzb_rho && pzb_rho_fix_reached &&
                  nOneAdditionalLayer && 0 > kLeast_rho[nOneAdditionalLayer].k ) {
                qzb_rho_temp = -(nOneAdditionalLayer+1);
                /* qzb_rho = -(nOneAdditionalLayer+1); *//* early rejection */
            }
            /* new code */
            if ( !qzb_rho && bRhoIsDiscrete ) {
                qzb_rho = CtPartCompareLayers( kLeast_rho, L_rho_fix_prev, 0 );
#if ( FIX_ChCh_CONSTIT_CANON_BUG == 1 )
                if ( qzb_rho ) {
                    int L_rho_diff = abs(qzb_rho)-1;
                    if ( L_rho_diff < L_rho_fix_prev ||
                         ((L_rho_diff == L_rho_fix_prev) &&
                         (kLeast_rho[L_rho_diff].i < I_rho_fix_prev)) ) {
                        qzb_rho = -(L_rho_diff+1); /* negative difference will be rejected */
                    }
                }
#endif
            }
            /* compare old results to new */
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
            if ( qzb_rho_temp && qzb_rho_temp != qzb_rho ) {
                int stop = 1; /* <BRKPT> */
            }
#endif
            if ( !qzb_rho ) {
                hz_rho = k;
            } else
            if ( qzb_rho < 0 ) {
                pCC->lNumRejectedCT ++;
            }
        }
        if ( qzb_rho > 0 || (!qzb_rho && !bRhoIsDiscrete) ) {
            /* found better rho */
            if ( !nNumLayers ) {
                CtPartCopy( pzb_rho, &Lambda, k-1 );
            }
        }
    }
/*L3:*/
    /*if ( hz_rho == k || (lab && qzb_rho >= 0 ) )*/
    /*if ( hz_zeta == k || hz_rho == k || (lab && qzb_rho >= 0 ) ) goto L4; else goto L6;*/
    if ( hz_zeta == k || hz_rho == k || (lab && qzb_rho >= 0 && QZFIX_OK(qzb_rho_fix) ) ) {
        /*L4: check for possible isomorphism or found a better rho */
        if ( PartitionIsDiscrete( &pi[k-1], n ) ) {
            pCC->lNumTotCT ++;
            goto L7;
        }
        PartitionGetFirstCell( &pi[k-1], W, k, n );
        v[k-1] = CellGetMinNode( &pi[k-1], &W[k-1], 0, pCD1 );
        if ( !dig && PartitionSatisfiesLemma_2_25(&pi[k-1], n) ) {
            ; /* found additional isomprphism */
        } else {
            t_Lemma = k+1;
        }
        e[k-1] = 0;  /* created new cell W[k-1] */
        goto L2;
    }
L6:
    /* a better rho or no good node was found at this level; return to smaller k */
    k2 = k;
    k = inchi_min(t_Lemma-1, inchi_max(t_eq_zeta-1, hz_rho));
    if ( k2 == t_Lemma )
        goto L13;
    /* store isomorphism found from Lemma 2.25. should be dig=0 !!! */
    if ( dig ) {
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
        int stop = 1;
#endif
        goto L13;
    } 
    l = inchi_min(l+1, L);
    PartitionGetMcrAndFixSet( &pi[t_Lemma-1], &Omega, &Phi, n_tg, l );
    goto L12;
L7:
    /* from L4: pi[k-1] is discrete */
    if ( h_zeta == 0 ) {
        /*goto L18;*/  /* error. the first T(nu) leaf was found */
        ret = CT_CANON_ERR;
        goto exit_error;
    }
    if ( k != hz_zeta )
        goto L8;
    /*  here zeta^gamma == nu */
    /*  if ( G^gamma == G ) goto L10; */
    if ( 0 == (res=CtFullCompare( &Lambda, &zf_zeta, 0, bSplitTautCompare )) ) {
        PartitionGetTransposition( &zeta, &pi[k-1], n_tg, &gamma );
        bZetaIsomorph = 1; /* for testing only */
        lNumEqlZeta ++;
        goto L10;
    } else
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    {
        int stop = 1;
    }
#endif
    /* !!! we should never come here !!! */
    if ( !nNumLayers ) {
        ret = -2;
        goto exit_error;
    }

L8: /* here nu is discrete: check rho for being a bettere leaf or isomorphism */
    /*if ( !lab || qzb_rho < 0 || !QZFIX_OK(qzb_rho_fix) )*/
    if ( !lab || ((qzb_rho < 0) && ( !pzb_rho_fix || qzb_rho_fix > 0 )) )
        goto L6;
    if ( pzb_rho_fix && kLeast_rho_fix && 0 == qzb_rho_fix ) {
        /* check for the rejection condition: Lambda > zb_rho_fix */
        if ( kLeast_rho_fix ) {
            int qzb_rho_fix_alt;
            qzb_rho_fix     = CtFullCompareLayers( kLeast_rho_fix );
            /* for debug only */
            qzb_rho_fix_alt =  CtFullCompare( &Lambda, pzb_rho_fix, 1, bSplitTautCompare );
            if ( qzb_rho_fix != qzb_rho_fix_alt ) {
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
                int stop = 1;
#endif
                qzb_rho_fix = qzb_rho_fix_alt;
            }
            /* end debug */
        } else {
            qzb_rho_fix = CtFullCompare( &Lambda, pzb_rho_fix, 1, bSplitTautCompare );
        }
        if ( !pzb_rho_fix_reached ) {
            pzb_rho_fix_reached = !qzb_rho_fix;
        }
        if ( 0 < qzb_rho_fix ) {
            /* Lambda > pzb_rho_fix, ignore this node */
            /* hzb_rho_fix = min( hzb_rho_fix, hz_rho ); */ /* ??? */
            qzb_rho_fix = 0;
            goto L6;
        }
        qzb_rho_fix = 0;
    }

    if ( qzb_rho < 0 )
        goto L6;
    if ( qzb_rho > 0 || !bRhoIsDiscrete )
        goto L9; /* note: p67 says k > PartitionSize( &rho, n ) */
    if ( k < r ) {
        goto L9; /* cannot understand it... */
    }

    /* !!! we should never come here if G(nu) != G(rho): CtPartCompare must be enough !!! */
    
    /* if ( G(nu) > G(rho) ) goto L9; */
    if ( kLeast_rho ) {
        int cur_qzb_alt;
        qzb_rho =     CtFullCompareLayers( kLeast_rho );
        /* for debug only */
        cur_qzb_alt = CtFullCompare( &Lambda, pzb_rho, 0, bSplitTautCompare );
        if ( qzb_rho != cur_qzb_alt ) {
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
            int stop = 1;
#endif
            qzb_rho = cur_qzb_alt;
        }
        /* end debug */
    } else {
        qzb_rho = CtFullCompare( &Lambda, pzb_rho, 0, bSplitTautCompare );
    }
    /* qzb_rho difference can be due to layers 1..MAX_LAYERS-1 only */
    if ( 0 < qzb_rho ) {
        /* CtFullCompare( &Lambda, pzb_rho, 0, bSplitTautCompare ); */
        qzb_rho = 0;
        goto L9;
    }
    /* if ( G(nu) < G(rho) ) goto L6; */
    if ( 0 > qzb_rho ) {
        qzb_rho = 0;
        goto L6;
    }
    /* nu^gamma == rho */
    if ( r != k ) {  /* if() is for debug only */
        r = k;
    }
    PartitionGetTransposition( &pi[k-1], &rho, n_tg, &gamma );
    bZetaIsomorph = 0; /* DCh */
    pCC->lNumEqualCT ++;
    goto L10;
L9:
    /* rho := nu; */
    PartitionCopy( &rho, &pi[k-1], n_tg );
    if ( nNumLayers ) {
        CtFullCopy( pzb_rho, &Lambda );
    }
    bZetaEqRho = 0;
    qzb_rho = 0;
    CtCompareLayersGetFirstDiff( kLeast_rho_fix, nOneAdditionalLayer,
                                 &L_rho_fix_prev, &I_rho_fix_prev, &k_rho_fix_prev );
    memset( kLeast_rho, 0, sizeof(kLeast_rho) );
    h_rho = hz_rho = k;
    CtPartInfinity( pzb_rho, qzb, k );
    pCC->lNumDecreasedCT ++;
    pCC->lNumEqualCT = 1;
    bRhoIsDiscrete = 1;
    goto L6;

L10: /* discrete pi[k-1] && G^gamma == G */

    pCC->lNumEqualCT += bZetaEqRho || !(bZetaIsomorph || qzb_rho);
    l = inchi_min(l+1, L);
    /* Omega[l] := mcr(gamma);
       Phi[l]   := fix(gamma);
    */
    TranspositionGetMcrAndFixSetAndUnorderedPartition( &gamma, &Omega, &Phi, n_tg, l, &theta_from_gamma );

    /*
    if ( theta(gamma) <= theta ) goto L11;
    theta := theta v theta(gamma);
    UnorderedPartitionJoin() returns 0 if theta_from_gamma is finer than theta,
    which means no changes in theta: theta_from_gamma ^ theta == theta.
    */
    if ( !UnorderedPartitionJoin( &theta_from_gamma, &theta, n_tg ) )
        goto L11; /* no new isomorphism found */
    /*  Output gamma (it is the Aut(G) generator) -- omitted -- */
    nNumFoundGenerators ++;
    /* if ( tvc in mcr(theta) ) goto L11; */
    if ( tvc == GetUnorderedPartitionMcrNode( &theta, tvc ) )
        goto L11;
    k = h_zeta;
    goto L13;
L11:
    k = lab? h_rho : h_zeta; /***Changed*** originally was k = h_rho; */
L12:
    /* if ( e[k-1] == 1 ) */
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    if ( e[k-1] == 1 && v[k-1] == INFINITY ) {
        int stop = 1;          /* <BRKPT> testing only */
    }
#endif
    if ( e[k-1] == 1 && v[k-1] != INFINITY ) { /* INFINITY for testing only */
        CellIntersectWithSet( &pi[k-1], &W[k-1], &Omega, l );
    }
L13:
    
    if ( (UserAction && USER_ACTION_QUIT == (*UserAction)()) ||
         (ConsoleQuit && (*ConsoleQuit)()) ) {
        ret = CT_USER_QUIT_ERR;
        goto exit_error;
    }
    if ( bInchiTimeIsOver(pCD->ulTimeOutTime) ) {
        ret = CT_TIMEOUT_ERR;
        goto exit_error;
    }
    
    if ( k == 0 )
        goto exit_function; /* stop */

    if ( lab && k < h_rho ) {  /***Added***/
        h_rho = k;
    }
    if ( k >  h_zeta ) {
        if ( v[k-1] == INFINITY ) {/*** Added by DCh for testing only ****/
            k --;
            goto L13;
        }
        goto L17;
    }
    if ( k == h_zeta )
        goto L14;
    h_zeta = k;
    tvc = tvh = CellGetMinNode( &pi[k-1], &W[k-1], 0, pCD1 );
L14:
    /* if v[k] and tvh are in the same cell of theta then index ++ */
    if ( GetUnorderedPartitionMcrNode( &theta, v[k-1] ) ==
         GetUnorderedPartitionMcrNode( &theta, tvh ) ) {
        index ++;
    }
    v[k-1] = CellGetMinNode(  &pi[k-1], &W[k-1], v[k-1], pCD1 );

    if ( v[k-1] == INFINITY )
        goto L16;
    if ( v[k-1] != GetUnorderedPartitionMcrNode( &theta, v[k-1] ) )
        goto L14;
L15:
    t_Lemma = inchi_min(t_Lemma, k+1);
    hz_zeta = inchi_min(hz_zeta, k);
/*    
    if ( lab && hz_rho >= k ) {
        hz_rho = k;
        qzb_rho = 0;
    }
*/
    if ( lab ) {
        if ( hz_rho >= k /*-1*/ )
            qzb_rho = 0;
        if ( hz_rho > k )
            hz_rho = k;
        UpdateCompareLayers( kLeast_rho, hz_rho );
    }
    if ( pzb_rho_fix ) {
        if ( hzb_rho_fix >= k /*-1*/ )
            qzb_rho_fix = 0;
        if ( hzb_rho_fix > k )
            hzb_rho_fix = k;
        UpdateCompareLayers( kLeast_rho_fix, hzb_rho_fix );
    }

    goto L2;
L16:
    if ( t_eq_zeta == k+1 && index == CellGetNumberOfNodes( &pi[k-1], &W[k-1] ) )
        t_eq_zeta = k;
    size *= (double)index;
    /******************** <<<===== A **************************/
    /* passed K times after passing point A. At these passes
       k = K, K-1, ..., 1 in this order
    */
    index = 0;
    k --;
    goto L13;
L17:
    /* if ( e[k-1] == 0 ) */
#if ( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    if ( e[k-1] == 0 && v[k-1] == INFINITY ) { /* testing only */
        int stop = 1;  /* <BRKPT> */
    }
#endif
    /*
    if ( e[k] == 0 set W[k] = Intersection(W[k], Omega[i]) for each i = 1..l,
         such that {v[1]..v[k-1]} in Phi[i]
    */
    if ( e[k-1] == 0 && v[k-1] != INFINITY ) /* Added v[k-1]!=... DCh */
    {
        NodeSetFromVertices( &cur_nodes, 1, v, k-1 );
        for ( i = 1; i <= l; i ++ ) {
            if ( AllNodesAreInSet( &cur_nodes, 1, &Phi, i ) ) {
                CellIntersectWithSet( &pi[k-1], &W[k-1],  &Omega, i );
            }
        }
    }
    e[k-1] = 1;
    v[k-1] = CellGetMinNode(  &pi[k-1], &W[k-1], v[k-1], pCD1 );
    if ( v[k-1] != INFINITY )
        goto L15;
    k --;
    goto L13;
/* L18: see above */

exit_function:
    /* CtPartFill( G, pCD, &rho, pzb_rho, 1, n, n_tg ); */
    if ( !bRhoIsDiscrete ) {
        ret = CT_CANON_ERR;
        goto exit_error;
    }
    if ( pzb_rho_fix ) {
        qzb_rho_fix = CtFullCompare( pzb_rho_fix, pzb_rho, 1, bSplitTautCompare );
        if ( qzb_rho_fix ) {
            ret = CT_CANON_ERR;
            goto exit_error;
        }
    }
    /* SymmRank */
    memset( nSymmRank, 0, n_tg * sizeof(nSymmRank[0]) );
    for ( i = 0; i < n_tg; i ++ ) {
        k = rho.AtNumber[i];
        k2 = (int)GetUnorderedPartitionMcrNode( &theta, (AT_NUMB)(k+1) ) - 1;
        if ( !nSymmRank[k2] || nSymmRank[k2] > rho.Rank[k] ) {
            nSymmRank[k2] = rho.Rank[k];
        }
    }
    for ( i = 0; i < n_tg; i ++ ) {
        k = rho.AtNumber[i];
        k2 = (int)GetUnorderedPartitionMcrNode( &theta, (AT_NUMB)(k+1) ) - 1;
        nSymmRank[k] = nSymmRank[k2];
    }
    /* CanonRank, nAtomNumberCanon */
    memcpy( nCanonRank, rho.Rank, n_tg * sizeof(nCanonRank[0]) );
    memcpy( nAtomNumberCanon, rho.AtNumber, n_tg * sizeof(nAtomNumberCanon[0]) );
    /* LinearCT */
    *nLenCt = pzb_rho->lenCt-1;
    if ( pCt ) {
        memcpy( pCt, pzb_rho->Ctbl, *nLenCt*sizeof(pCt[0]) );
    }
    pCC->lNumTotCT = pCC->lNumDecreasedCT + pCC->lNumRejectedCT + pCC->lNumEqualCT;
    pCC->dGroupSize             = size;
    pCC->lNumGenerators         = nNumFoundGenerators;
    pCC->lNumStoredIsomorphisms = l;
    /* Note: check nNumFoundGenerators */

    if ( pp_zb_rho_out && !*pp_zb_rho_out ) {
        *pp_zb_rho_out = pzb_rho;
        pzb_rho = NULL;
    }
    
exit_error:
    INCHI_HEAPCHK

    UnorderedPartitionFree( &theta );
    UnorderedPartitionFree( &theta_from_gamma );
    if ( W )   inchi_free( W );
    if ( v )   inchi_free( v );
    if ( e )   inchi_free( e );
    if ( qzb ) inchi_free( qzb );
    CTableFree( &Lambda );
    CTableFree( &zf_zeta );
    if ( pzb_rho ) {
        CTableFree( pzb_rho );
        inchi_free( pzb_rho );
        pzb_rho = NULL;
    }

/*    CTableFree( &zf_zeta2 ); */


    NodeSetFree( &Omega );
    NodeSetFree( &Phi );
  /* NodeSetFree( &mcr_theta, n, 1 ); */
    NodeSetFree( &cur_nodes );

    PartitionFree( &zeta );
/*    PartitionFree( &zeta2 ); */
    PartitionFree( &rho );
    TranspositionFree( &gamma );


    return ret;
}
/**********************************************************************************************
 * SetInitialRanks2: Set initial ranks in nRank according to pAtomInvariant[] values
 *                  Make sure enough prines have been generated.
 **********************************************************************************************/
/* Upon exit: */
/* nAtomNumber[i]: number (from 0) of an atom in the ith (from 0) position of the sorted order */
/* nNewRank[i]:    initial rank of the atom[i] based on atom invariants; from 1 to num_atoms */
/* Return value:   Number of different ranks */
int SetInitialRanks2( int num_atoms, ATOM_INVARIANT2* pAtomInvariant2, AT_RANK *nNewRank, AT_RANK *nAtomNumber )
{
    int i, nNumDiffRanks;
    AT_RANK nCurrentRank;

    for ( i = 0; i < num_atoms; i++ )
        nAtomNumber[i] = (AT_RANK)i;
    
    /* global for qsort */
    pAtomInvariant2ForSort = pAtomInvariant2;

    qsort( nAtomNumber, num_atoms, sizeof(nAtomNumber[0]), CompAtomInvariants2 );

    /* nNewRank[i]: non-decreading order; do not increment nCurrentRank */
    /*           if consecutive sorted atom invariants are identical */

    for ( i=num_atoms-1, nCurrentRank=nNewRank[nAtomNumber[i]] = (AT_RANK)num_atoms, nNumDiffRanks = 1; 0 < i ; i -- ) {
        /* Note: CompAtomInvariants2Only() in following line implicitly reads pAtomInvariant2 pointed by pAtomInvariant2ForSort */
        if ( CompAtomInvariants2Only( &nAtomNumber[i-1], &nAtomNumber[i] ) ) {
            nNumDiffRanks ++;
            nCurrentRank = (AT_RANK)i;
        }
        nNewRank[nAtomNumber[i - 1]] = nCurrentRank;
    }
    
    
    return nNumDiffRanks;
}

/****************************************************************************/
void FillOutAtomInvariant2( sp_ATOM* at, int num_atoms, int num_at_tg, ATOM_INVARIANT2* pAtomInvariant,
                           int bIgnoreIsotopic, int bHydrogensInRanks, int bHydrogensFixedInRanks,
                           int bDigraph, int bTautGroupsOnly, T_GROUP_INFO *t_group_info )
{
    int i, k, j, i_t_group;
    /* tautomers */
    T_GROUP          *t_group=NULL;
    int               num_t_groups     = 0;
    int               num_tautomer_iso = 0;
#define ELEM_NAME_LEN  2
    char ChemElements[ELEM_NAME_LEN*NUM_CHEM_ELEMENTS+ELEM_NAME_LEN];
    char CurElement[ELEM_NAME_LEN + ELEM_NAME_LEN], *pCurElem;
    int  nNumChemElements  = 0;
    int  nNumHydrogenAtoms = 0;
    int  nNumCarbonAtoms   = 0;
    memset( ChemElements, 0, sizeof(ChemElements) );
    memset( CurElement,   0, sizeof(CurElement)   );
    nNumChemElements = 0;

    if ( num_at_tg > num_atoms && t_group_info ) {
        t_group          = t_group_info->t_group;
        num_t_groups     = t_group_info->num_t_groups;
        num_tautomer_iso = t_group_info->bIgnoreIsotopic? 0 : T_NUM_ISOTOPIC;
    }

    if ( !bTautGroupsOnly ) {

        for ( i = 0; i < num_atoms; i ++ ) {
            if ( !strcmp( at[i].elname, "C" ) ) {
                nNumCarbonAtoms ++;
            } else
            if ( !strcmp( at[i].elname, "H" ) ||
                 !strcmp( at[i].elname, "D" ) ||
                 !strcmp( at[i].elname, "T" ) ) {
                nNumHydrogenAtoms ++;
            } else {
                CurElement[0] = at[i].elname[0];
                CurElement[1] = at[i].elname[1]? at[i].elname[1] : ' ';
                if ( ! (pCurElem = strstr( ChemElements, CurElement ) ) ) {
                    strcat( ChemElements, CurElement );
                    nNumChemElements ++;
                }
            }
        }
        if ( nNumChemElements > 1 ) {
            qsort( ChemElements, nNumChemElements, ELEM_NAME_LEN, CompChemElemLex );
        }
        if ( nNumCarbonAtoms ) {
            if ( nNumChemElements ) {
                memmove( ChemElements + ELEM_NAME_LEN, ChemElements, ELEM_NAME_LEN*nNumChemElements );
            }
            ChemElements[0] = 'C';
            ChemElements[1] = ' ';
            nNumChemElements ++;
        }
        if ( nNumHydrogenAtoms ) {
            ChemElements[ ELEM_NAME_LEN*nNumChemElements   ] = 'H';
            ChemElements[ ELEM_NAME_LEN*nNumChemElements+1 ] = ' ';
            nNumChemElements ++;
        }


        /* general */
        for ( i = 0; i < num_atoms; i ++ ) {
            memset( &pAtomInvariant[i], 0, sizeof(pAtomInvariant[0]) );
            CurElement[0] = at[i].elname[0];
            CurElement[1] = at[i].elname[1]? at[i].elname[1] : ' ';
            pCurElem = strstr( ChemElements, CurElement );
            if ( pCurElem ) {
                j = (pCurElem - ChemElements)/ELEM_NAME_LEN + 1;
            } else {
                j = nNumChemElements; /* must be D or T */
            }
            /* at[i].hill_type = (U_CHAR) j; */
            pAtomInvariant[i].val[AT_INV_HILL_ORDER] = j;

            pAtomInvariant[i].val[AT_INV_NUM_CONNECTIONS] = at[i].valence;
            if ( bHydrogensInRanks ) {
                pAtomInvariant[i].val[AT_INV_NUM_H] = ((t_group && at[i].endpoint>0)? 0 : at[i].num_H);
            }
            if ( bHydrogensFixedInRanks ) {
                pAtomInvariant[i].val[AT_INV_NUM_H_FIX] = ((t_group && at[i].endpoint>0)? at[i].num_H : 0);
            }
            if ( !bDigraph &&  t_group && (i_t_group = (int)at[i].endpoint-1) >= 0 && i_t_group < num_t_groups ) {
                pAtomInvariant[i].val[AT_INV_NUM_TG_ENDPOINTS]  = t_group[i_t_group].nNumEndpoints;
                for ( j = 0; j < T_NUM_NO_ISOTOPIC; j ++ ) {
                    pAtomInvariant[i].val[AT_INV_TG_NUMBERS+j] = t_group[i_t_group].num[j];
                }
                for ( j = 0; j < num_tautomer_iso; j ++ ) {
                    pAtomInvariant[i].val[AT_INV_TAUT_ISO+j] = t_group[i_t_group].num[j + T_NUM_NO_ISOTOPIC];
                }
            }
            pAtomInvariant[i].iso_sort_key = bIgnoreIsotopic? 0 : at[i].iso_sort_key;
        }
    } else {
        /* fill tautomeric groups only */
        memset ( pAtomInvariant, 0, num_at_tg*sizeof(pAtomInvariant[0]) );
    }
    /**************************************/
    /*          tautomeric groups         */
    /**************************************/
    for ( i = num_atoms; i < num_at_tg; i ++ ) {

        k = i - num_atoms;
        memset( &pAtomInvariant[i], 0, sizeof(pAtomInvariant[0]) );
        if ( !t_group )
            continue;
        /* make sure ranks of t-groups are larger than that of any atom */
         /* greater than for any real atom */
        pAtomInvariant[i].val[AT_INV_HILL_ORDER] = bTautGroupsOnly? num_at_tg : nNumChemElements+1;
        /* greater than for any real atom */
        pAtomInvariant[i].val[AT_INV_NUM_CONNECTIONS] = MAXVAL+1;
        if ( k < num_t_groups ) {
            pAtomInvariant[i].val[AT_INV_NUM_TG_ENDPOINTS] = t_group[k].nNumEndpoints;
            for ( j = 0; j < T_NUM_NO_ISOTOPIC; j ++ ) {
                pAtomInvariant[i].val[AT_INV_TAUT_ISO+j] = t_group[k].num[j];
            }
            for ( j = 0; j < num_tautomer_iso; j ++ ) {
                pAtomInvariant[i].val[AT_INV_TAUT_ISO+j] = t_group[k].num[j + T_NUM_NO_ISOTOPIC];
            }
        }
    }
}
/*****************************************************************************/
void CleanNumH( NUM_H *NumH, int len )
{
    int i;
    if ( NumH ) {
        for ( i = 0; i < len; i ++ ) {
            if ( NumH[i] == EMPTY_H_NUMBER ) {
                NumH[i] = 0;
            } else {
                NumH[i] -= BASE_H_NUMBER;
            }
        }
    }
}
/*****************************************************************************/
int CleanCt( AT_RANK *Ct, int len )
{
    if ( Ct && Ct[len] == EMPTY_CT ) {
        Ct[len] = 0;
        return 1;
    }
    return 0;
}
/*****************************************************************************/
void CleanIsoSortKeys( AT_ISO_SORT_KEY * isk, int len )
{
    int i;
    if ( isk ) {
        for ( i = 0; i < len; i ++ ) {
            if ( isk[i] == EMPTY_ISO_SORT_KEY ) {
                isk[i] = 0;
            }
        }
    }
}
/*****************************************************************************/
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
void MergeCleanIsoSortKeys( AT_ISO_SORT_KEY * isk1, AT_ISO_SORT_KEY * isk2, int len )
{
    int i;
    AT_ISO_SORT_KEY k1, k2;
    if ( isk1 && isk2 ) {
        for ( i = 0; i < len; i ++ ) {
            k1 = (isk1[i] == EMPTY_ISO_SORT_KEY)? 0 : isk1[i];
            k2 = (isk2[i] == EMPTY_ISO_SORT_KEY)? 0 : isk2[i];
            isk1[i] = k1 | k2;
        }
    } else
    if ( isk1 ) {
        CleanIsoSortKeys( isk1, len );
    }
}
#endif
#define FREE_CONTABLE(X) if (X) {CTableFree(X);inchi_free(X);}
#define FREE_ARRAY(X) if (X) inchi_free(X);
/*****************************************************************************/
void DeAllocBCN( BCN *pBCN )
{
    int    i, k;
    FTCN  *ftcn;
    if ( !pBCN )
        return;
    if ( pBCN->pRankStack ) {
        for ( i = 0; i < pBCN->nMaxLenRankStack; i ++ ) {
            FREE_ARRAY( pBCN->pRankStack[i] )
        }
        FREE_ARRAY( pBCN->pRankStack )
    }
    for ( k = 0; k < TAUT_NUM; k ++ ) {
        ftcn = pBCN->ftcn + k;
        FreeNeighList( ftcn->NeighList );
        FREE_ARRAY( ftcn->LinearCt )
        PartitionFree( &ftcn->PartitionCt );
        FREE_ARRAY( ftcn->nSymmRankCt )
        FREE_ARRAY( ftcn->nNumHOrig )
        FREE_ARRAY( ftcn->nNumH )
        FREE_ARRAY( ftcn->nNumHOrigFixH )
        FREE_ARRAY( ftcn->nNumHFixH )
        PartitionFree( &ftcn->PartitionCtIso );
        FREE_ARRAY( ftcn->nSymmRankCtIso )
        FREE_ARRAY( ftcn->iso_sort_keys )
        FREE_ARRAY( ftcn->iso_sort_keysOrig )
        FREE_ARRAY( ftcn->iso_exchg_atnos )
        FREE_ARRAY( ftcn->iso_exchg_atnosOrig )
    }
}
#undef FREE_CONTABLE
#undef FREE_ARRAY

/*****************************************************************************/
#if ( bRELEASE_VERSION == 0 && FIND_CANON_NE_EQUITABLE == 1 )
/* debug: find whether canonical equivalence is different from equitable partition */
int bCanonIsFinerThanEquitablePartition( int num_atoms, sp_ATOM* at, AT_RANK *nSymmRank )
{
    AT_RANK *nRank                    = NULL;
    AT_RANK *nAtomNumber              = NULL;
    AT_RANK *nTempRank                = NULL;
    AT_RANK nCurSymm, nCurRank;
    ATOM_INVARIANT2 *pAtomInvariant   = NULL;
    NEIGH_LIST      *NeighList        = NULL;
    int              nNumCurrRanks, i, is, ir, j;
    long             lCount;
    int              bIsNotSame  = 0;
    if ( at && nSymmRank ) {
        if ( !(nRank          = (AT_RANK*)inchi_calloc( num_atoms, sizeof(nRank[0]))) ||
             !(nAtomNumber    = (AT_RANK*)inchi_calloc( num_atoms, sizeof(nAtomNumber[0]))) ||
             !(nTempRank      = (AT_RANK*)inchi_calloc( num_atoms, sizeof(nTempRank[0]))) ||
             !(pAtomInvariant = (ATOM_INVARIANT2 *)inchi_calloc( num_atoms, sizeof(pAtomInvariant[0])))
            ) {
            goto exit_err;
        }
        if ( !(NeighList = CreateNeighList( num_atoms, num_atoms, at, 0, NULL )) ) {
            goto exit_err;
        }

        FillOutAtomInvariant2( at, num_atoms, num_atoms, pAtomInvariant, 1 /*bIgnoreIsotopic*/,
                               1 /*bHydrogensInRanks*/, 1 /*bHydrogensFixedInRanks*/, 0 /*bTaut=bDigraph*/,
                               0 /* bTautGroupsOnly */, NULL /*t_group_info*/ );
        /* initial partitioning of a hydrogenless skeleton: create equitable partition (assign initial ranks) */
        nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariant, nRank, nAtomNumber );

        lCount = 0;
        /* make equitable partition in pBCN->pRankStack[0,1] */
        nNumCurrRanks = DifferentiateRanks2( num_atoms, NeighList,
                                            nNumCurrRanks, nRank,
                                            nTempRank, nAtomNumber, &lCount, 0 /* 0 means use qsort */ );
        /* at this point the equitable partition is in nRank; the order of atoms is in nAtomNumber*/
        /* compare */
        nCurSymm = nCurRank = 0;
        for ( i = 0; i < num_atoms; i ++ ) {
            j = (int)nAtomNumber[i];
            if ( nCurSymm != nSymmRank[j] ) {
                nCurSymm = nSymmRank[j];
                is = i;
            }
            if ( nCurRank != nRank[j] ) {
                nCurRank = nRank[j];
                ir = i;
            }
            if ( is != ir ) {
                bIsNotSame = 1;
                break;
            }
        }
    }
exit_err:
    if ( nRank )
        inchi_free( nRank );
    if ( nAtomNumber )
        inchi_free( nAtomNumber );
    if ( nTempRank )
        inchi_free( nTempRank );
    if ( pAtomInvariant )
        inchi_free( pAtomInvariant );
    if ( NeighList ) 
        FreeNeighList( NeighList );
    return bIsNotSame;
}
#endif
/*****************************************************************************/
int GetBaseCanonRanking( int num_atoms, int num_at_tg, sp_ATOM* at[],
                         T_GROUP_INFO *t_group_info, ATOM_SIZES s[], BCN *pBCN, 
                         struct tagInchiTime *ulTimeOutTime, int bFixIsoFixedH )
{
    int ret = 0;
    int iBase;                   /* base structure index, always valid; = TAUT_YES except special fully non-taut mode */
    int iOther;                  /* other than basic structure index, usually non-taut; may be = iBase */
    int bReqNonTaut;             /* 1 => requested non-tautomeric results */
    int bReqTaut;                /* 1 => requested tautomeric results and the base structure is tautomeric */
    int bChanged;
    sp_ATOM *at_base        = NULL;
    sp_ATOM *at_other       = NULL;
    int bTautIgnoreIsotopic = 0;
    /*int bIgnoreIsotopic     = 0;*/
    int nNumCurrRanks       = 0;
    int nMaxLenRankStack    = 0;
    int num_max             = num_at_tg;
    long lCount;
    /* local allocations */
    ATOM_INVARIANT2 *pAtomInvariant    = NULL;
    NEIGH_LIST     *NeighList[TAUT_NUM];
    ConTable *Ct_Temp                 = NULL;

    /* initial partition for canonicalization */
    AT_RANK *nRank                    = NULL;
    AT_NUMB *nAtomNumber              = NULL;

    /* canonicalization output */
    
    ConTable *Ct_NoH                  = NULL;
    AT_RANK *nCanonRankNoH            = NULL;
    AT_NUMB *nAtomNumberCanonNoH      = NULL;
    AT_RANK *nSymmRankNoH             = NULL;
    
    ConTable *Ct_NoTautH              = NULL;
    AT_RANK *nSymmRankNoTautH         = NULL;
    AT_RANK *nCanonRankNoTautH        = NULL;
    AT_NUMB *nAtomNumberCanonNoTautH  = NULL;
    NUM_H   *numHNoTautH              = NULL;
    int      lenNumHNoTautH;
    int      maxlenNumHNoTautH;
    
    ConTable *Ct_Base                 = NULL;
    AT_RANK *nSymmRankBase            = NULL;
    AT_RANK *nCanonRankBase           = NULL;
    AT_NUMB *nAtomNumberCanonBase     = NULL;
    NUM_H   *numH                     = NULL;
    int      lenNumH;
    int      maxlenNumH               = 0;

#if ( USE_AUX_RANKING == 1 )
    AT_RANK *nRankAux                 = NULL;
    AT_NUMB *nAtomNumberAux           = NULL;
    ATOM_INVARIANT2 *pAtomInvariantAux= NULL;
#endif    


    ConTable *Ct_FixH                 = NULL;
    AT_RANK *nSymmRankFixH            = NULL;
    AT_RANK *nCanonRankFixH           = NULL;
    AT_NUMB *nAtomNumberCanonFixH     = NULL;
    NUM_H   *NumHfixed                = NULL;
    int      maxlenNumHfixed;

    /* isotopic canonicalization */
    
    ConTable *Ct_NoTautHIso               = NULL;
    AT_RANK *nSymmRankNoTautHIso          = NULL;
    AT_RANK *nCanonRankNoTautHIso         = NULL;
    AT_NUMB *nAtomNumberCanonNoTautHIso   = NULL;
    AT_ISO_SORT_KEY *iso_sort_key_NoTautH = NULL;
    int              maxlen_iso_sort_key_NoTautH;
    int              len_iso_sort_key_NoTautH;
    int num_iso_NoTautH, num_iso_NoAuxBase;

    ConTable *Ct_BaseIso                  = NULL;
    AT_RANK *nSymmRankBaseIso             = NULL;
    AT_RANK *nCanonRankBaseIso            = NULL;
    AT_NUMB *nAtomNumberCanonBaseIso      = NULL;
    
    AT_ISO_SORT_KEY *iso_sort_keyBase     = NULL;
    int              maxlen_iso_sort_keyBase;
    int              len_iso_sort_keyBase;

    int              bUseIsoAuxBase[TAUT_NUM];
    S_CHAR          *iso_exchg_atnos      = NULL;
    int              len_iso_exchg_atnos;
    int              maxlen_iso_exchg_atnos;
    int num_iso_Base;

    AT_ISO_SORT_KEY  iso_sort_key; 

    ConTable *Ct_FixHIso                  = NULL;
    AT_RANK *nSymmRankFixHIso             = NULL;
    AT_RANK *nCanonRankFixHIso            = NULL;
    AT_NUMB *nAtomNumberCanonFixHIso      = NULL;

#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    AT_ISO_SORT_KEY  iso_sort_key2; 
    AT_ISO_SORT_KEY *iso_sort_key_Hfixed  = NULL;
    int              maxlen_iso_sort_key_Hfixed;
    int              len_iso_sort_key_Hfixed;
    int num_iso_Hfixed;
#endif

    AT_RANK *nTempRank               = NULL;

    CANON_DATA    pCD[3]; /* = &CanonData; */
    CANON_COUNTS  CanonCounts;
    CANON_COUNTS *pCC = &CanonCounts;

    int i, j, k, m;
    int nCanonFlags[2];

    /*^^^ */
    int iflag;

    memset (pCD, 0, sizeof(pCD));
    memset (pCC, 0, sizeof(pCC[0]));
    memset ( bUseIsoAuxBase, 0, sizeof(bUseIsoAuxBase) );
    memset ( nCanonFlags, 0, sizeof(nCanonFlags) );
    NeighList[TAUT_NON] = NULL;
    NeighList[TAUT_YES] = NULL;

    /* select base structure, find whether it is tautomeric or not */
    if ( at[TAUT_YES] && s[TAUT_YES].nLenCT &&
         t_group_info && ((s[TAUT_YES].nLenLinearCTTautomer > 0 && /* ordinary tautomerism */
                          t_group_info->t_group && t_group_info->num_t_groups > 0) ||
                          /* protons have been moved */
                          (t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT) ||
                          /* tautomerism due to possible isotopic proton exchange */
                          (t_group_info->nNumIsotopicEndpoints > 1 &&
                          (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE|TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))) ) ) {
        /* tautomeric: (1) has tautomeric atoms OR
                       (2) H-atoms have been rearranged due to proton addition/removal OR
                       (3) Found isotopic H-atoms on tautomeric or hetero atoms
         */
        iBase    = TAUT_YES;
        bReqTaut = 1;
        bUseIsoAuxBase[iBase] = (s[iBase].nLenIsotopicEndpoints > 1) &&
                                (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE|TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE));
        if ( at[TAUT_NON] && s[TAUT_NON].nLenCT ) {
            iOther      = TAUT_NON; /* tautomeric and non-tautomeric */
            bReqNonTaut = 1;
        } else {
            iOther      = iBase; /* tautomeric only */
            bReqNonTaut = 0;
        }
    } else
    if ( at[TAUT_NON] && s[TAUT_NON].nLenCT ) {
        /* force pure non-tautomeric processing; happens for testing only */
        iBase       = TAUT_NON;
        bReqTaut    = 0;
        iOther      = iBase;
        bReqNonTaut = 1;
        num_at_tg   = num_atoms;
    } else
    if ( at[TAUT_YES] && s[TAUT_YES].nLenCT ) {
        /* although the user requested tautomeric processing, tautomerism has not been found */
        /* however, the results should be saved in the TAUT_YES elements of the arrays */
        iBase       = TAUT_YES;
        bReqTaut    = 0;
        bUseIsoAuxBase[iBase] = (s[iBase].nLenIsotopicEndpoints > 1);
        iOther      = iBase;
        bReqNonTaut = 1;
        num_at_tg   = num_atoms;
    } else {
        ret = CT_UNKNOWN_ERR;
        goto exit_error;
    }
    if ( bReqTaut ) {
        /* save "process isotopic" mark; temporarily set it to NO */
        bTautIgnoreIsotopic = t_group_info->bIgnoreIsotopic;
        t_group_info->bIgnoreIsotopic = 1;
    }
    lenNumH                  = num_atoms;

    /* isotopic canonicalization */
    num_iso_NoTautH             = 0;
    len_iso_sort_key_NoTautH    = 0;
    maxlen_iso_sort_key_NoTautH = 0;
    num_iso_Base                = 0;
    len_iso_sort_keyBase        = 0;
    maxlen_iso_sort_keyBase     = 0;
    len_iso_exchg_atnos         = 0;
    maxlen_iso_exchg_atnos      = 0;
    len_iso_exchg_atnos         = 0;
    maxlen_iso_exchg_atnos      = 0;

#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    num_iso_Hfixed              =
    len_iso_sort_key_Hfixed     =
    maxlen_iso_sort_key_Hfixed  = 0;
#endif

    /* prepare initial data */
    at_base  = at[iBase];
    at_other = at[iOther];
    pAtomInvariant      = (ATOM_INVARIANT2 *)inchi_calloc( num_max,     sizeof(pAtomInvariant[0]) );
    nSymmRankNoH        = (AT_RANK *)        inchi_calloc( num_max,     sizeof(nSymmRankNoH[0]       ) );
    nCanonRankNoH       = (AT_RANK *)        inchi_calloc( num_max,     sizeof(nCanonRankNoH[0]      ) );
    nAtomNumberCanonNoH = (AT_NUMB *)        inchi_calloc( num_max,     sizeof(nAtomNumberCanonNoH[0]) );
    nRank               = (AT_RANK *)        inchi_calloc( num_max,     sizeof(nRank[0]      ) );
    nAtomNumber         = (AT_NUMB *)        inchi_calloc( num_max,     sizeof(nAtomNumber[0]) );
    nTempRank           = (AT_RANK *)        inchi_calloc( num_max,     sizeof(nTempRank[0]  ) );

    if ( !pAtomInvariant ||
         !nSymmRankNoH   || !nCanonRankNoH    || !nAtomNumberCanonNoH ||
         !nRank          || !nAtomNumber      || !nTempRank             ) {
        goto exit_error_alloc;
    }
#if ( USE_AUX_RANKING == 1 )
    nRankAux            = (AT_RANK *)  inchi_calloc( num_max, sizeof(nRankAux[0]            ) );
    nAtomNumberAux      = (AT_NUMB *)  inchi_calloc( num_max, sizeof(nAtomNumberAux[0]      ) );
    pAtomInvariantAux   = (ATOM_INVARIANT2 *) inchi_malloc( num_max * sizeof(pAtomInvariantAux[0]) );
    if ( !nRankAux || !nAtomNumberAux || !pAtomInvariantAux ) {
        goto exit_error_alloc;
    }
#endif    

    if ( bReqTaut ) {
        if ( !(NeighList[TAUT_YES] = CreateNeighList( num_atoms, num_at_tg, at_base, 0, t_group_info )) )
            goto exit_error_alloc;
        /* needed for the hydrogenless structure */
        if ( !(NeighList[TAUT_NON] = CreateNeighList( num_atoms, num_atoms, at_base, 0, NULL )) )
            goto exit_error_alloc;
    } else {
        if ( !(NeighList[TAUT_NON] = CreateNeighList( num_atoms, num_atoms, at_base, 0, NULL )) )
            goto exit_error_alloc;
        NeighList[TAUT_YES] = NULL;
        INCHI_HEAPCHK
    }

    /* avoid memory leaks in case of error */
    /*
    pBCN->ftcn[TAUT_NON].NeighList          = NeighList[TAUT_NON];
    pBCN->ftcn[TAUT_YES].NeighList          = NeighList[TAUT_YES];
    */
    pBCN->nMaxLenRankStack                  = 0;
    pBCN->num_max                           = num_max;        /* allocated nRank[] arrays lengths in pRankStack */
    pBCN->num_at_tg                         = num_at_tg;  /* all of the following arrays have this length */
    pBCN->num_atoms                         = num_atoms;
    pBCN->ulTimeOutTime                     = ulTimeOutTime;

    /* initial partitioning of a hydrogenless skeleton: fill out the inveriant */
    FillOutAtomInvariant2( at_base, num_atoms, num_atoms, pAtomInvariant, 1 /*bIgnoreIsotopic*/,
                           0 /*bHydrogensInRanks*/, 0 /*bHydrogensFixedInRanks*/, 0 /*bTaut=bDigraph*/,
                           0 /* bTautGroupsOnly */, NULL /*t_group_info*/ );
    /* initial partitioning of a hydrogenless skeleton: create equitable partition (assign initial ranks) */
    nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariant, nRank, nAtomNumber );

    lCount = 0;
    /* make equitable partition in pBCN->pRankStack[0,1] */
    nNumCurrRanks = DifferentiateRanks2( num_atoms, NeighList[TAUT_NON],
                                        nNumCurrRanks, nRank,
                                        nTempRank, nAtomNumber, &lCount, 0 /* 0 means use qsort */ );

    /* allocate partition stack */
    nMaxLenRankStack = 2*(num_at_tg-nNumCurrRanks) + 8;  /* was 2*(...) + 6 */
    pBCN->pRankStack = (AT_RANK **) inchi_calloc( nMaxLenRankStack, sizeof(pBCN->pRankStack[0]) );
    if ( !pBCN->pRankStack ) {
        pBCN->nMaxLenRankStack = 0; /* avoid memory leaks in case of error */
        goto exit_error_alloc;
    }
    pBCN->nMaxLenRankStack = nMaxLenRankStack; /* avoid memory leaks in case of error */
    /* init partition stack */
    pBCN->pRankStack[0] = nRank;
    pBCN->pRankStack[1] = nAtomNumber;
    
    /********************************************************************************************/
    /* get NoH/no taut groups  canonical numbering, connection table, and equivalence partition */
    /********************************************************************************************/

    /* pointers */
    pCD[iOther].LinearCT                   = NULL;
    pCD[iOther].NumH                       = NULL;
    pCD[iOther].NumHfixed                  = NULL;
    pCD[iOther].iso_sort_key               = NULL;
    pCD[iOther].iso_exchg_atnos            = NULL;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
    pCD[iOther].iso_sort_key_Hfixed        = NULL;
#endif                                     
    /* variables - unchanged */            
    pCD[iOther].ulTimeOutTime              = pBCN->ulTimeOutTime;
    pCD[iOther].nMaxLenLinearCT            = s[iOther].nLenCTAtOnly + 1;
    /* return values & input/output */     
    pCD[iOther].nLenLinearCT               = s[iOther].nLenCTAtOnly;
    pCD[iOther].nLenCTAtOnly               = s[iOther].nLenCTAtOnly;
    pCD[iOther].lenNumH                    = 0;
    pCD[iOther].lenNumHfixed               = 0;
    pCD[iOther].len_iso_sort_key           = 0;
    pCD[iOther].maxlen_iso_sort_key        = 0;
    pCD[iOther].len_iso_exchg_atnos        = 0;
    pCD[iOther].maxlen_iso_exchg_atnos     = 0;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
    pCD[iOther].len_iso_sort_key_Hfixed    = 0;
    pCD[iOther].maxlen_iso_sort_key_Hfixed = 0;
#endif
    ret = CanonGraph01( num_atoms, num_atoms, num_max, 0, NeighList[TAUT_NON], (Partition *)pBCN->pRankStack,
                      nSymmRankNoH,  nCanonRankNoH, nAtomNumberCanonNoH, pCD+iOther, pCC, NULL, &Ct_NoH );
    if ( ret < 0 ) {
        goto exit_error;
    }
    /* update initial partitioning */
    nNumCurrRanks = FixCanonEquivalenceInfo( num_atoms, nSymmRankNoH, nRank, nTempRank, nAtomNumber, &bChanged );
    /* repartition if necessary */
    if ( bChanged & 3 ) {
        if ( Ct_NoH ) {
            CTableFree( Ct_NoH );
            inchi_free( Ct_NoH );
            Ct_NoH = NULL;
        }
        pCD[iOther].nCanonFlags |= CANON_FLAG_NO_H_RECANON;

        ret = CanonGraph02( num_atoms, num_atoms, num_max, 0, NeighList[TAUT_NON], (Partition *)pBCN->pRankStack,
                          nSymmRankNoH,  nCanonRankNoH, nAtomNumberCanonNoH, pCD+iOther, pCC, NULL, &Ct_NoH );
        if ( ret < 0 ) {
            goto exit_error;
        }
    }
    /********************************************************************************/
    /* get NoTautH canonical numbering, connection table, and equivalence partition */
    /********************************************************************************/
    maxlenNumHNoTautH        = num_atoms + 1;
    nSymmRankNoTautH         = (AT_RANK *)  inchi_calloc( num_max, sizeof(nSymmRankNoTautH[0]       ) );
    nCanonRankNoTautH        = (AT_RANK *)  inchi_calloc( num_max, sizeof(nCanonRankNoTautH[0]      ) );
    nAtomNumberCanonNoTautH  = (AT_NUMB *)  inchi_calloc( num_max, sizeof(nAtomNumberCanonNoTautH[0]) );
    numHNoTautH              = (NUM_H *)    inchi_calloc( maxlenNumHNoTautH, sizeof(numHNoTautH[0]) );
    if ( !numHNoTautH || !nSymmRankNoTautH || !nCanonRankNoTautH || !nAtomNumberCanonNoTautH ) {
        goto exit_error_alloc;
    }
    /* find number of H atoms attached to not-a-tautomeric-endpoint atoms */
    for ( i = 0; i < num_atoms; i ++ ) {
        numHNoTautH[i] = (!at_base[i].endpoint && at_base[i].num_H)? at_base[i].num_H+BASE_H_NUMBER : EMPTY_H_NUMBER;
    }
    /* pointers */
    pCD[iOther].LinearCT                   = NULL;
    pCD[iOther].NumH                       = numHNoTautH;
    pCD[iOther].NumHfixed                  = NULL;
    pCD[iOther].iso_sort_key               = NULL;
    pCD[iOther].iso_exchg_atnos            = NULL;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
    pCD[iOther].iso_sort_key_Hfixed        = NULL;
#endif
    /* variables - unchanged */             
    pCD[iOther].ulTimeOutTime              = pBCN->ulTimeOutTime;
    pCD[iOther].nMaxLenLinearCT            = s[iOther].nLenCTAtOnly + 1;
    pCD[iOther].maxlenNumH                 = maxlenNumHNoTautH;
    /* return values & input/output */     
    pCD[iOther].nLenLinearCT               = s[iOther].nLenCTAtOnly;
    pCD[iOther].nLenCTAtOnly               = s[iOther].nLenCTAtOnly;
    pCD[iOther].lenNumH                    = lenNumHNoTautH = num_atoms;
    pCD[iOther].lenNumHfixed               = 0;
    pCD[iOther].len_iso_sort_key           = 0;
    pCD[iOther].maxlen_iso_sort_key        = 0;
    pCD[iOther].len_iso_exchg_atnos        = 0;
    pCD[iOther].maxlen_iso_exchg_atnos     = 0;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
    pCD[iOther].len_iso_sort_key_Hfixed    = 0;
    pCD[iOther].maxlen_iso_sort_key_Hfixed = 0;
#endif
    pCD[iOther].nAuxRank                   = NULL;

    /* check whether we need NoTautH cononicalization */
    memset( nTempRank, 0, num_max * sizeof(nTempRank[0]) );
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( nTempRank[nSymmRankNoH[i]-1] < i ) {
            nTempRank[nSymmRankNoH[i]-1] = i; /* greatest class representative */
        }
    }
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( numHNoTautH[i] != numHNoTautH[nTempRank[nSymmRankNoH[i]-1]] ) {
            pCD[iOther].nCanonFlags |= CANON_FLAG_NO_TAUT_H_DIFF;
            break; /* atoms so far found to be equivalent have different number of H; the canonicalization is needed */
        }
    }
    /* i = 0; *//* debug: force to call the canonicalization */
    if ( i < num_atoms ) {
        /* needs canonicalization */
        /* get aux canonical ranking of the structure with attached H */
#if ( USE_AUX_RANKING == 1 )
        /* refine no-H partition according to not-a-taut-H distribution */
        memset( pAtomInvariantAux, 0, num_max * sizeof(pAtomInvariantAux[0]) );
        for ( i = 0; i < num_atoms; i ++ ) {
            pAtomInvariantAux[i].val[0] = nSymmRankNoH[i];
            pAtomInvariantAux[i].val[1] = numHNoTautH[i]; /* additional differentiation: not-a-taut-H distribution */
        }
        /* initial partitioning */
        nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariantAux, nRankAux, nAtomNumberAux );
        /* make equitable partition */
        nNumCurrRanks = DifferentiateRanks2( num_atoms, NeighList[TAUT_NON],
                                            nNumCurrRanks, nRankAux,
                                            nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means use qsort */ );
        /* to accelerate do not call CanonGraph() to find really equivalent atoms */
        pCD[iOther].nAuxRank = nRankAux;
#endif    

        ret = CanonGraph03( num_atoms, num_atoms, num_max, 1 /* digraph?? was 0 */, NeighList[TAUT_NON], (Partition *)pBCN->pRankStack,
                          nSymmRankNoTautH,  nCanonRankNoTautH, nAtomNumberCanonNoTautH, pCD+iOther, pCC, &Ct_NoH, &Ct_NoTautH );
        if ( ret < 0 ) {
            goto exit_error;
        }
        /* in case of non-tautomeric structure the final results are in:

                   nSymmRankNoTautH
                   nCanonRankNoTautH
                   nAtomNumberCanonNoTautH
                   Ct_NoTautH
                   numHNoTautH (original H positions)
        */
    } else {
        /* copy the results of the previous (no H) canonicalization */
        /* in this case numHNoTautH[] is not needed for the next canonicalization(s) */
        if ( (Ct_Temp = (ConTable *)inchi_calloc( 1, sizeof( *Ct_Temp ) ) ) &&
             CTableCreate( Ct_Temp, num_atoms, pCD+iOther) ) {
            CtFullCopy( Ct_Temp, Ct_NoH );
            /* since Ct_NoH does not have Ct_NoH->NumH we have to fill out Ct_Temp->NumH separately */
            for ( i = 0; i < num_atoms; i ++ ) {
                Ct_Temp->NumH[nCanonRankNoH[i]-1] = numHNoTautH[i];
                /*Ct_Temp->NumH[i] = numHNoTautH[nAtomNumberCanonNoH[i]]; -- alternative */
            }
            Ct_Temp->lenNumH = num_atoms;
        } else {
            goto exit_error_alloc;
        }
        Ct_NoTautH = Ct_Temp;
        Ct_Temp    = NULL;
        memcpy( nSymmRankNoTautH, nSymmRankNoH, num_atoms*sizeof(nSymmRankNoTautH[0]) );
        memcpy( nCanonRankNoTautH, nCanonRankNoH, num_atoms*sizeof(nCanonRankNoTautH[0]) );
        memcpy( nAtomNumberCanonNoTautH, nAtomNumberCanonNoH, num_atoms*sizeof(nAtomNumberCanonNoTautH[0]) );
    }
    /* in case of non-tautomeric component this is the final result */
    /* i = CtFullCompare( Ct_NoTautH, Ct_Temp, num_atoms, 0, 0 );*/

    /*******************************************************************************************/
    /* If only Isotopic atoms and isotopic H, tautomerism has not been found:                  */
    /* get isotopic canonical numbering, connection table, and equivalence partition           */
    /*******************************************************************************************/
    
    if ( s[iOther].num_isotopic_atoms && !s[iOther].bIgnoreIsotopic && !bReqTaut && bReqNonTaut ) {
        
        maxlen_iso_sort_key_NoTautH = num_atoms+1;
        nSymmRankNoTautHIso        = (AT_RANK *)  inchi_calloc( num_max, sizeof(nSymmRankNoTautHIso[0]       ) );
        nCanonRankNoTautHIso       = (AT_RANK *)  inchi_calloc( num_max, sizeof(nCanonRankNoTautHIso[0]       ) );
        nAtomNumberCanonNoTautHIso = (AT_NUMB *)  inchi_calloc( num_max, sizeof(nAtomNumberCanonNoTautHIso[0]) );
        iso_sort_key_NoTautH       = (AT_ISO_SORT_KEY *) inchi_calloc( maxlen_iso_sort_key_NoTautH, sizeof(iso_sort_key_NoTautH[0]) );

        if ( !nSymmRankNoTautHIso || !nCanonRankNoTautHIso || !nAtomNumberCanonNoTautHIso || !iso_sort_key_NoTautH ) {
            goto exit_error_alloc;
        }

        /* fill out isotopic non-tautomeric keys */
        num_iso_NoTautH = 0;
        for ( i = 0; i < num_atoms; i ++ ) {
            if ( at_base[i].endpoint ) {
                /* should not happen */
                iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, 0, 0, 0);
            } else {
                iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2]);
            }
            if ( iso_sort_key ) {
                iso_sort_key_NoTautH[i] = iso_sort_key;
                num_iso_NoTautH ++;
            } else {
                iso_sort_key_NoTautH[i] = EMPTY_ISO_SORT_KEY;
            }
        }
        /* pointers */
        pCD[iOther].LinearCT                   = NULL; /* LinearCT; */
        pCD[iOther].NumH                       = numHNoTautH;
        pCD[iOther].NumHfixed                  = NULL;
        pCD[iOther].iso_sort_key               = iso_sort_key_NoTautH;
        pCD[iOther].iso_exchg_atnos            = NULL;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
        pCD[iOther].iso_sort_key_Hfixed        = NULL;
#endif
        /* variables - unchanged */             
        pCD[iOther].ulTimeOutTime              = pBCN->ulTimeOutTime;
        pCD[iOther].nMaxLenLinearCT            = s[iOther].nLenCTAtOnly + 1;
        pCD[iOther].maxlenNumH                 = maxlenNumHNoTautH;
        /* return values & input/output */    
        pCD[iOther].nLenLinearCT               = s[iOther].nLenCTAtOnly;
        pCD[iOther].nLenCTAtOnly               = s[iOther].nLenCTAtOnly;
        pCD[iOther].lenNumH                    = lenNumHNoTautH /*= num_atoms*/;
        pCD[iOther].lenNumHfixed               = 0;
        pCD[iOther].len_iso_sort_key           = len_iso_sort_key_NoTautH = num_atoms;
        pCD[iOther].maxlen_iso_sort_key        = maxlen_iso_sort_key_NoTautH;
        pCD[iOther].len_iso_exchg_atnos        = 0;
        pCD[iOther].maxlen_iso_exchg_atnos     = 0;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
        pCD[iOther].len_iso_sort_key_Hfixed    = 0;
        pCD[iOther].maxlen_iso_sort_key_Hfixed = 0;
#endif
        pCD[iOther].nAuxRank                   = NULL;

        if ( num_iso_NoTautH ) {
            /* check whether we need NoTautH cononicalization */
            memset( nTempRank, 0, num_max * sizeof(nTempRank[0]) );
            for ( i = 0; i < num_atoms; i ++ ) {
                if ( nTempRank[nSymmRankNoTautH[i]-1] < i ) {
                    nTempRank[nSymmRankNoTautH[i]-1] = i; /* greatest class representative */
                }
            }
            for ( i = 0; i < num_atoms; i ++ ) {
                if ( iso_sort_key_NoTautH[i] != iso_sort_key_NoTautH[nTempRank[nSymmRankNoTautH[i]-1]] ) {
                    pCD[iOther].nCanonFlags |= CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF;
                    break; /* atoms so far found to be equivalent differ in isotopes; the canonicalization is needed */
                }
            }
        } else {
            i = num_atoms;
        }
        /* i = 0; *//* debug: force to call the canonicalization */
        if ( i < num_atoms ) {
            /* we need canonicalization */
            /* get aux canonical ranking of the structure with isotopic non-tautomeric H */

#if ( USE_AUX_RANKING == 1 )
            /* refine no-taut-H partition according to non-taut H isotopic distribution */
            memset( pAtomInvariantAux, 0, num_max * sizeof(pAtomInvariantAux[0]) );
            for ( i = 0; i < num_atoms; i ++ ) {
                pAtomInvariantAux[i].val[0] = nSymmRankNoTautH[i];
                pAtomInvariantAux[i].iso_sort_key = iso_sort_key_NoTautH[i]; /* additional differentiation */
            }
            /* initial ranks for non-taut H isotopic distribution */
            nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariantAux, nRankAux, nAtomNumberAux );
            /* make equitable */
            nNumCurrRanks = DifferentiateRanks2( num_atoms, NeighList[TAUT_NON],
                                                nNumCurrRanks, nRankAux,
                                                nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means use qsort */ );
            /* to accelerate do not call CanonGraph() to find really equivalent atoms */
            pCD[iOther].nAuxRank = nRankAux;
#endif    

            ret = CanonGraph04( num_atoms, num_atoms, num_max, 1 /* digraph?? was 0 */, NeighList[TAUT_NON], (Partition *)pBCN->pRankStack,
                              nSymmRankNoTautHIso,  nCanonRankNoTautHIso, nAtomNumberCanonNoTautHIso, pCD+iOther, pCC, &Ct_NoTautH, &Ct_NoTautHIso );
            if ( ret < 0 ) {
                goto exit_error;
            }
            /* in case of non-tautomeric structure the final results are in:

                       nSymmRankNoTautHIso
                       nCanonRankNoTautHIso
                       nAtomNumberCanonNoTautHIso
                       Ct_NoTautHIso
                       iso_sort_key_NoTautH (original isotopic atom positions)
            */
        } else {
            /* copy the results of the previous (no taut H) canonicalization */
            /* in this case numHNoTautH[] is not needed for the next canonicalization(s) */
            if ( (Ct_Temp = (ConTable *)inchi_calloc( 1, sizeof( *Ct_Temp ) ) ) &&
                 CTableCreate( Ct_Temp, num_atoms, pCD+iOther) ) {
                CtFullCopy( Ct_Temp, Ct_NoTautH );
                /* since Ct_NoTautH does not have Ct_NoTautH->iso_sort_key we have to fill out Ct_Temp->iso_sort_key separately */
                for ( i = 0; i < num_atoms; i ++ ) {
                    Ct_Temp->iso_sort_key[nCanonRankNoTautH[i]-1] = iso_sort_key_NoTautH[i];
                }
                Ct_Temp->len_iso_sort_key = num_atoms;
            } else {
                goto exit_error_alloc;
            }
            Ct_NoTautHIso = Ct_Temp;
            Ct_Temp    = NULL;
            memcpy( nSymmRankNoTautHIso,  nSymmRankNoTautH, num_atoms*sizeof(nSymmRankNoTautHIso[0]) );
            memcpy( nCanonRankNoTautHIso, nCanonRankNoTautH, num_atoms*sizeof(nCanonRankNoTautHIso[0]) );
            memcpy( nAtomNumberCanonNoTautHIso, nAtomNumberCanonNoTautH, num_atoms*sizeof(nAtomNumberCanonNoTautHIso[0]) );
        }
        /* in case of non-tautomeric component this is the final result */
        /* i = CtFullCompare( Ct_NoTautHIso, Ct_Temp, num_atoms, 0, 0 );*/
    }


    if ( bReqTaut ) {
        /*****************************************************************************/
        /* Tautomeric Structure Canonicalizaton:                                     */
        /* get base canonical numbering, connection table, and equivalence partition */
        /*****************************************************************************/
        /* find H atoms attached to non-tautomeric-endpoints and to tautomeric endpoints */
        maxlenNumH            = num_atoms + T_NUM_NO_ISOTOPIC*(num_at_tg-num_atoms) + 1; /* including negative charges */
        nSymmRankBase         = (AT_RANK *)  inchi_calloc( num_max, sizeof(nSymmRankBase[0]       ) );
        nCanonRankBase        = (AT_RANK *)  inchi_calloc( num_max, sizeof(nCanonRankBase[0]      ) );
        nAtomNumberCanonBase  = (AT_NUMB *)  inchi_calloc( num_max, sizeof(nAtomNumberCanonBase[0]) );
        numH                  = (NUM_H *)    inchi_calloc( maxlenNumH, sizeof(numH[0]) );
        if ( !numH || !nSymmRankBase || !nCanonRankBase || !nAtomNumberCanonBase ) {
            goto exit_error_alloc;
        }
        /* non-tautomeric H counts */
        for ( i = 0; i < num_atoms; i ++ ) {
            numH[i] = (!at_base[i].endpoint && at_base[i].num_H)? at_base[i].num_H+BASE_H_NUMBER : EMPTY_H_NUMBER;
        }
        /* tautomeric H and negative charge counts */
        for ( i = k = num_atoms; i < num_at_tg; i ++ ) {
            m = i-num_atoms;
            for ( j = 0; j < T_NUM_NO_ISOTOPIC; j ++ ) {
                /* non-zeroes for j=1 are negative charge counts; T_NUM_NO_ISOTOPIC=2 entry per t-group */
                numH[k ++] = t_group_info->t_group[m].num[j]? t_group_info->t_group[m].num[j]+BASE_H_NUMBER : EMPTY_H_NUMBER;
            }
        }
        /* pointers */
        pCD[iBase].LinearCT                   = NULL;
        pCD[iBase].NumH                       = numH; /* num_atoms non-tautomeric H; num_tg pairs of H and (-) in t-groups */
        pCD[iBase].NumHfixed                  = NULL;
        pCD[iBase].iso_sort_key               = NULL;
        pCD[iBase].iso_exchg_atnos            = NULL;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
        pCD[iBase].iso_sort_key_Hfixed        = NULL;
#endif
        /* variables - unchanged */              
        pCD[iBase].ulTimeOutTime              = pBCN->ulTimeOutTime;
        pCD[iBase].nMaxLenLinearCT            = s[iBase].nLenCT + 1;
        pCD[iBase].maxlenNumH                 = maxlenNumH;
        /* return values & input/output */    
        pCD[iBase].nLenLinearCT               = s[iBase].nLenCT;
        pCD[iBase].nLenCTAtOnly               = s[iBase].nLenCTAtOnly;
        pCD[iBase].lenNumH                    = lenNumH = k;
        pCD[iBase].lenNumHfixed               = 0;
        pCD[iBase].len_iso_sort_key           = 0;
        pCD[iBase].maxlen_iso_sort_key        = 0;
        pCD[iBase].len_iso_exchg_atnos        = 0;
        pCD[iBase].maxlen_iso_exchg_atnos     = 0;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
        pCD[iBase].len_iso_sort_key_Hfixed    = 0;
        pCD[iBase].maxlen_iso_sort_key_Hfixed = 0;
#endif
        pCD[iBase].nAuxRank                   = NULL;

        /* make sure the initial partition is equitable (at this point t-groups do not have ranks yet) */
        FillOutAtomInvariant2( at_base, num_atoms, num_at_tg, pAtomInvariant, 1 /*bIgnoreIsotopic*/,
                       0 /*bHydrogensInRanks*/, 0 /*bHydrogensFixedInRanks*/, 1 /*bTaut=bDigraph*/,
                       1 /* bTautGroupsOnly */, t_group_info );
        for ( i = 0; i < num_atoms; i ++ ) {
            pAtomInvariant[i].val[0] = pBCN->pRankStack[0][i];
        }
        /* initial ranks for t-group(s) only */
        nNumCurrRanks = SetInitialRanks2( num_at_tg, pAtomInvariant, nRank, nAtomNumber );
        /* make equitable, call digraph procedure;
           pBCN->pRankStack[0] is nRank, pBCN->pRankStack[1] is nAtomNumber
           This should only split ranks of tautomeric groups */
        nNumCurrRanks = DifferentiateRanks4( num_at_tg, NeighList[TAUT_YES], 
                                         nNumCurrRanks, pBCN->pRankStack[0], nTempRank /* temp array */,
                                         pBCN->pRankStack[1],  (AT_RANK)num_atoms, &lCount );
#if ( USE_AUX_RANKING == 1 )
        /* refine no-H partition according to non-taut H distribution */
        memset( pAtomInvariantAux, 0, num_max * sizeof(pAtomInvariantAux[0]) );
        for ( i = 0; i < num_atoms; i ++ ) {
            pAtomInvariantAux[i].val[0] = nSymmRankNoTautH[i];
            pAtomInvariantAux[i].val[1] = numH[i]; /* additional differentiation */
        }
        for ( j = i; i < num_at_tg; i ++ ) {
            pAtomInvariantAux[i].val[0] = nRank[i];
        }

        /* initial ranks for t-group(s) */
        nNumCurrRanks = SetInitialRanks2( num_at_tg, pAtomInvariantAux, nRankAux, nAtomNumberAux );
        /* make equitable, call digraph procedure */
        nNumCurrRanks = DifferentiateRanks4( num_at_tg, NeighList[TAUT_YES], 
                                         nNumCurrRanks, nRankAux, nTempRank /* temp array */,
                                         nAtomNumberAux,  (AT_RANK)num_atoms, &lCount );
        /* to accelerate do not call CanonGraph() to find really equivalent atoms */
        pCD[iBase].nAuxRank = nRankAux;
#endif    


        ret = CanonGraph05( num_atoms, num_at_tg, num_max, 1 /* digraph*/, NeighList[TAUT_YES], (Partition *)pBCN->pRankStack,
                          nSymmRankBase,  nCanonRankBase, nAtomNumberCanonBase, pCD+iBase, pCC, &Ct_NoTautH, &Ct_Base );
        if ( ret < 0 ) {
            goto exit_error;
        }
    
        /* tautomeric isotopic structure */
        /**************************************************************************************/
        /* Isotopic atoms and isotopic H atoms and isotopic tautomeric groups                 */
        /* get isotopic canonical numbering, connection table, and equivalence partition      */
        /**************************************************************************************/
        if ( (s[iBase].num_isotopic_atoms        && !s[iBase].bIgnoreIsotopic) ||
             (s[iBase].bHasIsotopicTautGroups    && !bTautIgnoreIsotopic) || 
             (bUseIsoAuxBase[iBase]              && !bTautIgnoreIsotopic) ) {
                               
            t_group_info->bIgnoreIsotopic = bTautIgnoreIsotopic;

            nSymmRankBaseIso        = (AT_RANK *)  inchi_calloc( num_max, sizeof(nSymmRankBaseIso[0]       ) );
            nCanonRankBaseIso       = (AT_RANK *)  inchi_calloc( num_max, sizeof(nCanonRankBaseIso[0]      ) );
            nAtomNumberCanonBaseIso = (AT_NUMB *)  inchi_calloc( num_max, sizeof(nAtomNumberCanonBaseIso[0]) );
            if ( bUseIsoAuxBase[iBase] ) {
                maxlen_iso_exchg_atnos = num_max+1;
                iso_exchg_atnos     = (S_CHAR  *)  inchi_calloc( maxlen_iso_exchg_atnos, sizeof(iso_exchg_atnos[0]) );
            } 
            maxlen_iso_sort_keyBase = num_max+1; /* num_at_tg+1;*/
            iso_sort_keyBase        = (AT_ISO_SORT_KEY *) inchi_calloc( maxlen_iso_sort_keyBase, sizeof(iso_sort_keyBase[0]) );
            if ( !nSymmRankBaseIso || !nCanonRankBaseIso || !nAtomNumberCanonBaseIso ||
                 !iso_sort_keyBase ||
                 (maxlen_iso_exchg_atnos  && !iso_exchg_atnos) ) {
                goto exit_error_alloc;
            }
            /* atoms */
            num_iso_NoTautH   = 0;
            num_iso_NoAuxBase = 0;
            if ( iso_exchg_atnos ) {
                len_iso_exchg_atnos = num_at_tg;
            }
            for ( i = 0; i < num_atoms; i ++ ) {
                if ( at_base[i].endpoint || (iso_exchg_atnos && (at_base[i].cFlags & AT_FLAG_ISO_H_POINT)) ) {
                    /* tautomeric or may have exchangeable isotopic H */
                    iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, 0, 0, 0);
                    if ( iso_exchg_atnos ) {
                        num_iso_NoAuxBase += !at_base[i].endpoint; /* these non-taut atom may exchange isotopic H as tautomeric atoms do */
                    }
                } else {
                    /* non-mobile H */
                    iso_sort_key = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2]);
                    if ( iso_exchg_atnos ) {
                        iso_exchg_atnos[i] = 1; /* atom cannot have exchangable isotopic H atom(s) */
                    }
                }
                if ( iso_sort_key ) {
                    num_iso_NoTautH ++;
                    iso_sort_keyBase[i] = iso_sort_key;
                } else {
                    iso_sort_keyBase[i] = EMPTY_ISO_SORT_KEY;
                }
            }
            /* check marking and count of non-taut atoms that may exchange isotopic H -- debug only */
            if ( iso_exchg_atnos ) {
                if ( num_iso_NoAuxBase != t_group_info->nIsotopicEndpointAtomNumber[0] ) {
                    ret = CT_ISOCOUNT_ERR;
                    goto exit_error;
                }
                for ( i = 1; i <= num_iso_NoAuxBase; i ++ ) {
                    j = t_group_info->nIsotopicEndpointAtomNumber[i];
                    if ( at_base[j].endpoint || !(at_base[j].cFlags & AT_FLAG_ISO_H_POINT) ) {
                        ret = CT_ISOCOUNT_ERR;
                        goto exit_error;
                    }
                }
            }
            /* t-groups */
            num_iso_Base = 0;
            if ( iso_exchg_atnos ) {
                for ( i = num_atoms; i < num_at_tg; i ++ ) {
                    iso_sort_keyBase[i] = EMPTY_ISO_SORT_KEY; /* new mode: do not provide info about isotopic tautomeric H */
                }
            } else {
                for ( i = num_atoms; i < num_at_tg; i ++ ) { /* should not happen anymore */
                    m = i-num_atoms;
                    if ( (iso_sort_key = t_group_info->t_group[m].iWeight) ) {
                        /* old approach: each t-group has its own isotopic "weight" */
                        num_iso_Base ++;
                        iso_sort_keyBase[i] = iso_sort_key;
                    } else {
                        iso_sort_keyBase[i] = EMPTY_ISO_SORT_KEY;
                    }
                }
            }
            if ( !num_iso_NoAuxBase && iso_exchg_atnos ) {
                /* all atoms that may exchange isotopic H are either tautomeric or not present */
                inchi_free( iso_exchg_atnos );
                iso_exchg_atnos = NULL;
                len_iso_exchg_atnos    = 0;
                maxlen_iso_exchg_atnos = 0;
            }
            if ( !num_iso_NoTautH && !num_iso_Base && iso_sort_keyBase ) {
                /* no isotopic atoms present */
                inchi_free( iso_sort_keyBase );
                iso_sort_keyBase = NULL;
                maxlen_iso_sort_keyBase = 0;
            } else {
                len_iso_sort_keyBase = num_at_tg;
            }
            if ( !iso_exchg_atnos && !iso_sort_keyBase ) {
                /* no isotopic part at all or only tautomeric groups */
                inchi_free( nSymmRankBaseIso );        nSymmRankBaseIso        = NULL;
                inchi_free( nCanonRankBaseIso );       nCanonRankBaseIso       = NULL;
                inchi_free( nAtomNumberCanonBaseIso ); nAtomNumberCanonBaseIso = NULL;
            } else {
                /* proceed with tautomeric isotopic canonicalization */
                /* pointers */
                pCD[iBase].LinearCT                   = NULL;
                pCD[iBase].NumH                       = numH; /* num_atoms non-tautomeric H; num_tg pairs of H and (-) in t-groups */
                pCD[iBase].NumHfixed                  = NULL;
                pCD[iBase].iso_sort_key               = iso_sort_keyBase;
                pCD[iBase].iso_exchg_atnos            = iso_exchg_atnos;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
                pCD[iBase].iso_sort_key_Hfixed        = NULL;
#endif
                /* variables - unchanged */              
                pCD[iBase].ulTimeOutTime              = pBCN->ulTimeOutTime;
                pCD[iBase].nMaxLenLinearCT            = s[iBase].nLenCT + 1;
                pCD[iBase].maxlenNumH                 = maxlenNumH;
                /* return values & input/output */    
                pCD[iBase].nLenLinearCT               = s[iBase].nLenCT;
                pCD[iBase].nLenCTAtOnly               = s[iBase].nLenCTAtOnly;
                pCD[iBase].lenNumH                    = lenNumH /* = k */;
                pCD[iBase].lenNumHfixed               = 0;
                pCD[iBase].len_iso_sort_key           = len_iso_sort_keyBase;
                pCD[iBase].maxlen_iso_sort_key        = maxlen_iso_sort_keyBase;
                pCD[iBase].len_iso_exchg_atnos        = len_iso_exchg_atnos;
                pCD[iBase].maxlen_iso_exchg_atnos     = maxlen_iso_exchg_atnos;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
                pCD[iBase].len_iso_sort_key_Hfixed    = 0;
                pCD[iBase].maxlen_iso_sort_key_Hfixed = 0;
#endif
                pCD[iBase].nAuxRank                   = NULL;
        
                if ( num_iso_NoTautH || num_iso_Base || num_iso_NoAuxBase ) {
                    /* check whether we need actual canonicalization */
                    memset( nTempRank, 0, num_max * sizeof(nTempRank[0]) );
                    for ( i = 0; i < num_at_tg; i ++ ) {
                        if ( nTempRank[nSymmRankBase[i]-1] < i ) {
                            nTempRank[nSymmRankBase[i]-1] = i; /* greatest class representative */
                        }
                    }
                    for ( i = 0; i < num_at_tg; i ++ ) {
                        if ( (iso_sort_keyBase? (iso_sort_keyBase[i] != iso_sort_keyBase[nTempRank[nSymmRankBase[i]-1]]):0) ||
                             (iso_exchg_atnos? (iso_exchg_atnos[i] != iso_exchg_atnos[nTempRank[nSymmRankBase[i]-1]]):0)) {
                            pCD[iBase].nCanonFlags |= CANON_FLAG_ISO_TAUT_DIFF;
                            break; /* atoms so far found to be equivalent have different number of H; the canonicalization is needed */
                        }
                    }
                } else {
                    i = num_at_tg; /* should not happen */
                }
                /* i = 0; *//* debug: force to call the canonicalization */
                if ( i < num_at_tg ) {
                    /* we need canonicalization */
                    /* get aux canonical ranking of the structure with isotopic non-tautomeric H */

    #if ( USE_AUX_RANKING == 1 )
                    /* refine no-taut-H partition according to non-taut H + t-groups isotopic distribution */
                    memset( pAtomInvariantAux, 0, num_max * sizeof(pAtomInvariantAux[0]) );
                    for ( i = 0; i < num_at_tg; i ++ ) {
                        pAtomInvariantAux[i].val[0] = nSymmRankBase[i];
                        pAtomInvariantAux[i].iso_sort_key = iso_sort_keyBase? iso_sort_keyBase[i] : 0; /* additional differentiation */
                        pAtomInvariantAux[i].iso_aux_key  = iso_exchg_atnos? iso_exchg_atnos[i] : 0;
                    }
                    /* initial ranks for non-taut H isotopic distribution */
                    nNumCurrRanks = SetInitialRanks2( num_at_tg, pAtomInvariantAux, nRankAux, nAtomNumberAux );
                    /* make equitable, not a digraph procedure */
                    nNumCurrRanks = DifferentiateRanks2( num_at_tg, NeighList[TAUT_YES],
                                                        nNumCurrRanks, nRankAux,
                                                        nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means first use qsort */ );
                    /* to accelerate do not call CanonGraph() to find really equivalent atoms */
                    pCD[iBase].nAuxRank = nRankAux;
    #endif    


                    ret = CanonGraph06( num_atoms, num_at_tg, num_max, 1 /* digraph */, NeighList[TAUT_YES], (Partition *)pBCN->pRankStack,
                                      nSymmRankBaseIso,  nCanonRankBaseIso, nAtomNumberCanonBaseIso, pCD+iBase, pCC, &Ct_Base, &Ct_BaseIso );
                    if ( ret < 0 ) {
                        goto exit_error;
                    }
                    /* in case of a tautomeric structure the final results are in:

                               nSymmRankBaseIso
                               nCanonRankBaseIso
                               nAtomNumberCanonBaseIso
                               Ct_BaseIso
                               iso_sort_keyBase (original isotopic atom & t-group positions)
                               Ct_BaseIso->iso_exchg_atnos: 0=>can exchange isotopic H, including tautomeric atoms
                               iso_exchg_atnos            : same, in order of t_group_info->nIsotopicEndpointAtomNumber[]
                    */
                } else {
                    /* copy the results of the previous (no taut H) canonicalization */
                    /* in this case numHNoTautH[] is not needed for the next canonicalization(s) */
                    if ( (Ct_Temp = (ConTable *)inchi_calloc( 1, sizeof( *Ct_Temp ) ) ) &&
                         CTableCreate( Ct_Temp, num_atoms, pCD+iBase) ) {
                        CtFullCopy( Ct_Temp, Ct_Base );
                        /* since Ct_Base does not have Ct_Base->iso_sort_key we
                           have to fill out Ct_Temp->iso_sort_key separately */
                        if ( iso_sort_keyBase ) {
                            for ( i = 0; i < num_at_tg; i ++ ) {
                                Ct_Temp->iso_sort_key[nCanonRankBase[i]-1] = iso_sort_keyBase[i];
                            }
                            Ct_Temp->len_iso_sort_key = num_at_tg;
                        } else {
                            Ct_Temp->len_iso_sort_key = 0;
                        }
                        if ( iso_exchg_atnos ) {
                            for ( i = 0; i < num_atoms; i ++ ) {
                                Ct_Temp->iso_exchg_atnos[nCanonRankBase[i]-1] = iso_exchg_atnos[i];
                            }
                            Ct_Temp->len_iso_exchg_atnos = num_at_tg;
                        } else {
                            Ct_Temp->len_iso_exchg_atnos = 0;
                        }
                    } else {
                        goto exit_error_alloc;
                    }
                    Ct_BaseIso = Ct_Temp;
                    Ct_Temp    = NULL;
                    memcpy( nSymmRankBaseIso,  nSymmRankBase, num_at_tg*sizeof(nSymmRankBaseIso[0]) );
                    memcpy( nCanonRankBaseIso, nCanonRankBase, num_at_tg*sizeof(nCanonRankBaseIso[0]) );
                    memcpy( nAtomNumberCanonBaseIso, nAtomNumberCanonBase, num_at_tg*sizeof(nAtomNumberCanonBaseIso[0]) );
                }
                /* in case of non-tautomeric component this is the final result */
                /* i = CtFullCompare( Ct_BaseIso, Ct_Temp, num_at_tg, 0, 0 );*/

                t_group_info->bIgnoreIsotopic = 1;
            }
        }
    }

    /**********************************************************************************/
    /* get "fixed H" canonical numbering, connection table, and equivalence partition */
    /**********************************************************************************/

    if ( bReqTaut && bReqNonTaut ) {
        maxlenNumHfixed       = num_atoms + 1;
        nSymmRankFixH         = (AT_RANK *)  inchi_calloc( num_max, sizeof(nSymmRankFixH[0]       ) );
        nCanonRankFixH        = (AT_RANK *)  inchi_calloc( num_max, sizeof(nCanonRankFixH[0]      ) );
        nAtomNumberCanonFixH  = (AT_NUMB *)  inchi_calloc( num_max, sizeof(nAtomNumberCanonFixH[0]) );
        NumHfixed             = (NUM_H *)    inchi_calloc( maxlenNumHfixed, sizeof(NumHfixed[0]) );
        if ( !NumHfixed || !nSymmRankFixH || !nCanonRankFixH || !nAtomNumberCanonFixH ) {
            goto exit_error_alloc;
        }
        for ( i = 0; i < num_atoms; i ++ ) {
            /* fixed and non-tautomeric H different in taut and non-taut structures */
            if ( at_base[i].endpoint ) {
                NumHfixed[i] = at_other[i].num_H? at_other[i].num_H+BASE_H_NUMBER : EMPTY_H_NUMBER;
            } else
            if ( at_other[i].num_H != at_base[i].num_H ) {
                NumHfixed[i] = (NUM_H)at_other[i].num_H - (NUM_H)at_base[i].num_H + BASE_H_NUMBER;
            } else {
                NumHfixed[i] = EMPTY_H_NUMBER;
            }
        }
        /* pointers */
        pCD[iOther].LinearCT                   = NULL; /* LinearCT; */
        pCD[iOther].NumH                       = numHNoTautH;
        pCD[iOther].NumHfixed                  = NumHfixed;/* variables - unchanged */
        pCD[iOther].iso_sort_key               = NULL;
        pCD[iOther].iso_exchg_atnos            = NULL;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
        pCD[iOther].iso_sort_key_Hfixed        = NULL;
#endif
        pCD[iOther].ulTimeOutTime              = pBCN->ulTimeOutTime;
        pCD[iOther].nMaxLenLinearCT            = s[iOther].nLenCTAtOnly + 1;
        pCD[iOther].maxlenNumH                 = maxlenNumHNoTautH;
        pCD[iOther].maxlenNumHfixed            = maxlenNumHfixed;
        /* return values & input/output */     
        pCD[iOther].nLenLinearCT               = s[iOther].nLenCTAtOnly;
        pCD[iOther].nLenCTAtOnly               = s[iOther].nLenCTAtOnly;
        pCD[iOther].lenNumH                    = lenNumHNoTautH = num_atoms;
        pCD[iOther].lenNumHfixed               = num_atoms;
        pCD[iOther].len_iso_sort_key           = 0;
        pCD[iOther].maxlen_iso_sort_key        = 0;
        pCD[iOther].len_iso_exchg_atnos        = 0;
        pCD[iOther].maxlen_iso_exchg_atnos     = 0;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
        pCD[iOther].len_iso_sort_key_Hfixed    = 0;
        pCD[iOther].maxlen_iso_sort_key_Hfixed = 0;
#endif
        pCD[iOther].nAuxRank                   = NULL;

#if ( USE_AUX_RANKING == 1 )
        if ( !nRankAux )
            nRankAux            = (AT_RANK *)  inchi_calloc( num_max, sizeof(nRankAux[0]            ) );
        if ( !nAtomNumberAux )
            nAtomNumberAux      = (AT_NUMB *)  inchi_calloc( num_max, sizeof(nAtomNumberAux[0]      ) );
        if ( !pAtomInvariantAux )
            pAtomInvariantAux   = (ATOM_INVARIANT2 *) inchi_malloc( num_max * sizeof(pAtomInvariantAux[0]) );

        if ( !nRankAux || !nAtomNumberAux ||
             !pAtomInvariantAux ) {
            goto exit_error_alloc;
        }
        /* refine no-H partition according to non-taut H distribution */
        memset( pAtomInvariantAux, 0, num_max * sizeof(pAtomInvariantAux[0]) );
        for ( i = 0; i < num_atoms; i ++ ) {
            pAtomInvariantAux[i].val[0] = nSymmRankBase[i];
            pAtomInvariantAux[i].val[1] = NumHfixed[i]; /* additional differentiation */
        }

        /* initial ranks for t-group(s) */
        nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariantAux, nRankAux, nAtomNumberAux );
        /* make equitable, digraph procedure */
        nNumCurrRanks = DifferentiateRanks2( num_atoms, NeighList[TAUT_NON],
                                            nNumCurrRanks, nRankAux,
                                            nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means use qsort */ );
        /* to accelerate do not call CanonGraph() to find really equivalent atoms */
        pCD[iOther].nAuxRank = nRankAux;
#endif    
    
        ret = CanonGraph07( num_atoms, num_atoms, num_max, 0, NeighList[TAUT_NON], (Partition *)pBCN->pRankStack,
                          nSymmRankFixH,  nCanonRankFixH, nAtomNumberCanonFixH, pCD+iOther, pCC, &Ct_NoTautH, &Ct_FixH );
        if ( ret < 0 ) {
            goto exit_error;
        }

        /*******************************************************************************************/
        /* get "fixed H" isotopic canonical numbering, connection table, and equivalence partition */
        /*******************************************************************************************/
        iflag = (s[iBase].num_isotopic_atoms && !s[iBase].bIgnoreIsotopic) ||
             (s[iBase].bHasIsotopicTautGroups && !bTautIgnoreIsotopic);
        if (bFixIsoFixedH) /* #if ( FIX_ISO_FIXEDH_BUG == 1 )  */
             /* fix bug when iso H was removed as a proton and fixed-H isotopic layer is missing -  2008-09-24 DT*/
             iflag = iflag || (s[iOther].num_isotopic_atoms && !s[iOther].bIgnoreIsotopic);
        if (iflag) {

#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
            maxlen_iso_sort_key_Hfixed  =
#endif
            maxlen_iso_sort_key_NoTautH = num_atoms+1;
            nSymmRankFixHIso        = (AT_RANK *)  inchi_calloc( num_max, sizeof(nSymmRankFixHIso[0]       ) );
            nCanonRankFixHIso       = (AT_RANK *)  inchi_calloc( num_max, sizeof(nCanonRankFixHIso[0]      ) );
            nAtomNumberCanonFixHIso = (AT_NUMB *)  inchi_calloc( num_max, sizeof(nAtomNumberCanonFixHIso[0]) );
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
            iso_sort_key_Hfixed     = (AT_ISO_SORT_KEY *) inchi_calloc( maxlen_iso_sort_key_Hfixed, sizeof(iso_sort_key_Hfixed[0]) );
#endif
            iso_sort_key_NoTautH    = (AT_ISO_SORT_KEY *) inchi_calloc( maxlen_iso_sort_key_NoTautH, sizeof(iso_sort_key_NoTautH[0]) );

            if ( !nSymmRankFixHIso || !nCanonRankFixHIso || !nAtomNumberCanonFixHIso ||
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
                 !iso_sort_key_Hfixed ||
#endif
                 !iso_sort_key_NoTautH ) {
                goto exit_error_alloc;
            }

#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
            /* fill out isotopic non-tautomeric keys */
            for ( i = 0; i < num_atoms; i ++ ) {
                if ( at_base[i].endpoint ) {
                    iso_sort_key  = make_iso_sort_key( at_base[i].iso_atw_diff, 0, 0, 0);
                    iso_sort_key2 = make_iso_sort_key( 0, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2]);
                } else {
                    iso_sort_key  = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2]);
                    iso_sort_key2 = 0;
                }
                if ( iso_sort_key ) {
                    iso_sort_key_NoTautH[i] = iso_sort_key;
                    num_iso_NoTautH ++;
                } else {
                    iso_sort_key_NoTautH[i] = EMPTY_ISO_SORT_KEY;
                }
                if ( iso_sort_key2 ) {
                    num_iso_Hfixed ++;
                    iso_sort_key_Hfixed[i] = iso_sort_key2;
                } else {
                    iso_sort_key_Hfixed[i] = EMPTY_ISO_SORT_KEY;
                }
            }
#else
            /* fill out isotopic non-tautomeric keys */
            for ( i = 0; i < num_atoms; i ++ ) {

                if (bFixIsoFixedH) /* #if ( FIX_ISO_FIXEDH_BUG == 1 )  */
                {
                    /* fix bug when iso H was removed as a proton and fixed-H isotopic layer is missing -  2008-09-24 DT*/
                    if ( at_other ) 
                    {
                        iso_sort_key  = make_iso_sort_key( at_other[i].iso_atw_diff, at_other[i].num_iso_H[0], at_other[i].num_iso_H[1], at_other[i].num_iso_H[2]);
                    } 
                    else 
                    {
                        iso_sort_key  = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2]);
                    }
                }
                else
                    iso_sort_key  = make_iso_sort_key( at_base[i].iso_atw_diff, at_base[i].num_iso_H[0], at_base[i].num_iso_H[1], at_base[i].num_iso_H[2]);




                if ( iso_sort_key ) {
                    iso_sort_key_NoTautH[i] = iso_sort_key;
                    num_iso_NoTautH ++;
                } else {
                    iso_sort_key_NoTautH[i] = EMPTY_ISO_SORT_KEY;
                }
            }
#endif
            /* pointers */
            pCD[iOther].LinearCT                   = NULL; /* LinearCT; */
            pCD[iOther].NumH                       = numHNoTautH;
            pCD[iOther].NumHfixed                  = NumHfixed;/* variables - unchanged */
            pCD[iOther].iso_sort_key               = iso_sort_key_NoTautH;
            pCD[iOther].iso_exchg_atnos            = NULL;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
            pCD[iOther].iso_sort_key_Hfixed        = iso_sort_key_Hfixed;
#endif
            pCD[iOther].ulTimeOutTime              = pBCN->ulTimeOutTime;
            pCD[iOther].nMaxLenLinearCT            = s[iOther].nLenCTAtOnly + 1;
            pCD[iOther].maxlenNumH                 = maxlenNumHNoTautH;
            pCD[iOther].maxlenNumHfixed            = maxlenNumHfixed;
            /* return values & input/output */     
            pCD[iOther].nLenLinearCT               = s[iOther].nLenCTAtOnly;
            pCD[iOther].nLenCTAtOnly               = s[iOther].nLenCTAtOnly;
            pCD[iOther].lenNumH                    = lenNumHNoTautH = num_atoms;
            pCD[iOther].lenNumHfixed               = num_atoms;
            pCD[iOther].len_iso_sort_key           = len_iso_sort_key_NoTautH = num_atoms;
            pCD[iOther].maxlen_iso_sort_key        = maxlen_iso_sort_key_NoTautH;
            pCD[iOther].len_iso_exchg_atnos        = 0;
            pCD[iOther].maxlen_iso_exchg_atnos     = 0;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
            pCD[iOther].len_iso_sort_key_Hfixed    = len_iso_sort_key_Hfixed  = num_atoms;
            pCD[iOther].maxlen_iso_sort_key_Hfixed = maxlen_iso_sort_key_Hfixed;
#endif
            pCD[iOther].nAuxRank                   = NULL;

#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
            if ( num_iso_Hfixed || num_iso_NoTautH )
#else
            if ( num_iso_NoTautH )
#endif
            {
                /* check whether we need NoTautH cononicalization */
                memset( nTempRank, 0, num_max * sizeof(nTempRank[0]) );
                for ( i = 0; i < num_atoms; i ++ ) {
                    if ( nTempRank[nSymmRankFixH[i]-1] < i ) {
                        nTempRank[nSymmRankFixH[i]-1] = i; /* greatest class representative */
                    }
                }
                for ( i = 0; i < num_atoms; i ++ ) {
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )       
                    if ( iso_sort_key_Hfixed[i] != iso_sort_key_Hfixed[nTempRank[nSymmRankFixH[i]-1]] )
                        break;
#endif
                    if ( iso_sort_key_NoTautH[i] != iso_sort_key_NoTautH[nTempRank[nSymmRankFixH[i]-1]])
                        break; /* atoms so far found to be equivalent have different isotopic shifts; the canonicalization is needed */
                }
            } else {
                i = num_atoms; /* should not happen */
            }
            /* i = 0; *//* debug: force to call the canonicalization */
            if ( i < num_atoms ) {
                pCD[iOther].nCanonFlags |= CANON_FLAG_ISO_FIXED_H_DIFF;
                /* we need canonicalization */
                /* get aux canonical ranking of the structure with isotopic non-tautomeric H */

#if ( USE_AUX_RANKING == 1 )
                /* refine fixed-taut-H partition according to the isotopic distribution */
                memset( pAtomInvariantAux, 0, num_max * sizeof(pAtomInvariantAux[0]) );
                for ( i = 0; i < num_atoms; i ++ ) {
                    pAtomInvariantAux[i].val[0] = nSymmRankFixH[i];
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
                    iso_sort_key = 0;
                    if (iso_sort_key_NoTautH[i]!=EMPTY_ISO_SORT_KEY)
                        iso_sort_key |= iso_sort_key_NoTautH[i];
                    if (iso_sort_key_Hfixed[i] !=EMPTY_ISO_SORT_KEY)
                        iso_sort_key |= iso_sort_key_Hfixed[i];
                    if ( !iso_sort_key )
                        iso_sort_key = EMPTY_ISO_SORT_KEY;
#else
                    iso_sort_key = iso_sort_key_NoTautH[i];
#endif
                    pAtomInvariantAux[i].iso_sort_key = iso_sort_key; /* additional differentiation */
                }

                /* initial ranks for non-taut H isotopic distribution */
                nNumCurrRanks = SetInitialRanks2( num_atoms, pAtomInvariantAux, nRankAux, nAtomNumberAux );
                /* make equitable, digraph procedure */
                nNumCurrRanks = DifferentiateRanks2( num_atoms, NeighList[TAUT_NON],
                                                    nNumCurrRanks, nRankAux,
                                                    nTempRank, nAtomNumberAux, &lCount, 0 /* 0 means use qsort */ );
                /* to accelerate do not call CanonGraph() to find really equivalent atoms */
                pCD[iOther].nAuxRank = nRankAux;
#endif    


                ret = CanonGraph08( num_atoms, num_atoms, num_max, 1 /* digraph?? was 0 */, NeighList[TAUT_NON], (Partition *)pBCN->pRankStack,
                                  nSymmRankFixHIso,  nCanonRankFixHIso, nAtomNumberCanonFixHIso, pCD+iOther, pCC, &Ct_FixH, &Ct_FixHIso );
                if ( ret < 0 ) {
                    goto exit_error;
                }
                /* in case of non-tautomeric structure the final results are in:

                           nSymmRankFixHIso
                           nCanonRankFixHIso
                           nAtomNumberCanonFixHIso
                           Ct_FixHIso
                           iso_sort_keyBase     ([0..num_atoms] original isotopic atom positions)
                           iso_sort_key_Hfixed  (original fixed tautomeric H distribution)
                */
            } else {
                /* copy the results of the previous (no taut H) canonicalization */
                /* in this case numHNoTautH[] is not needed for the next canonicalization(s) */
                if ( (Ct_Temp = (ConTable *)inchi_calloc( 1, sizeof( *Ct_Temp ) ) ) &&
                     CTableCreate( Ct_Temp, num_atoms, pCD+iOther) ) {
                    CtFullCopy( Ct_Temp, Ct_FixH );
                    /* since Ct_FixH does not have Ct_FixH->iso_sort_key and Ct_FixH->iso_sort_key_Hfixed we
                       have to fill out Ct_Temp->iso_sort_key and Ct_Temp->iso_sort_key_Hfixed separately */
                    for ( i = 0; i < num_atoms; i ++ ) {
                        Ct_Temp->iso_sort_key[nCanonRankFixH[i]-1]        = iso_sort_key_NoTautH[i];
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
                        Ct_Temp->iso_sort_key_Hfixed[nCanonRankFixH[i]-1] = iso_sort_key_Hfixed[i];
#endif
                    }
                    Ct_Temp->len_iso_sort_key        = num_atoms;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
                    Ct_Temp->len_iso_sort_key_Hfixed = num_atoms;
#endif
                    /*Ct_Temp->lenNumH = num_atoms;*/
                } else {
                    goto exit_error_alloc;
                }
                Ct_FixHIso = Ct_Temp;
                Ct_Temp    = NULL;
                memcpy( nSymmRankFixHIso,  nSymmRankFixH, num_atoms*sizeof(nSymmRankFixHIso[0]) );
                memcpy( nCanonRankFixHIso, nCanonRankFixH, num_atoms*sizeof(nCanonRankFixHIso[0]) );
                memcpy( nAtomNumberCanonFixHIso, nAtomNumberCanonFixH, num_atoms*sizeof(nAtomNumberCanonFixHIso[0]) );
            }
            /* in case of non-tautomeric component this is the final result */
            /* i = CtFullCompare( Ct_NoTautHIso, Ct_Temp, num_atoms, 0, 0 );*/
        }
    } /* "fixed H" canonical numbering */

    /* consistency check: compare canonical connection tables, H-atoms, isotopic H & taut groups */
    ret = 0;
    ret |= (Ct_NoH->lenCt != Ct_NoTautH->lenCt) || memcmp(Ct_NoH->Ctbl, Ct_NoTautH->Ctbl, Ct_NoH->lenCt * sizeof(Ct_NoH->Ctbl[0]));
    if ( bReqTaut ) {
    if ( Ct_FixH ) {
    ret |= (Ct_NoTautH->lenCt != Ct_FixH->lenCt) || memcmp(Ct_NoTautH->Ctbl, Ct_FixH->Ctbl, Ct_NoTautH->lenCt * sizeof(Ct_NoTautH->Ctbl[0]));
    ret |= (Ct_NoTautH->lenNumH != Ct_FixH->lenNumH) || memcmp( Ct_NoTautH->NumH, Ct_FixH->NumH, Ct_NoTautH->lenNumH*sizeof(Ct_Base->NumH[0]));
    }
    ret |= (Ct_NoTautH->lenCt > Ct_Base->lenCt) || memcmp(Ct_NoTautH->Ctbl, Ct_Base->Ctbl, Ct_NoTautH->lenCt * sizeof(Ct_NoTautH->Ctbl[0]));
    ret |= (Ct_NoTautH->lenNumH > Ct_Base->lenNumH) || memcmp( Ct_NoTautH->NumH, Ct_Base->NumH, Ct_NoTautH->lenNumH*sizeof(Ct_Base->NumH[0]));
    }

    /* isotopic canonicalization */
    if ( Ct_NoTautHIso ) {
    ret |= (Ct_NoH->lenCt != Ct_NoTautHIso->lenCt) || memcmp(Ct_NoH->Ctbl, Ct_NoTautHIso->Ctbl, Ct_NoH->lenCt * sizeof(Ct_NoH->Ctbl[0]));
    ret |= (Ct_NoTautH->lenNumH != Ct_NoTautHIso->lenNumH) || memcmp( Ct_NoTautH->NumH, Ct_NoTautHIso->NumH, Ct_NoTautH->lenNumH*sizeof(Ct_Base->NumH[0]));
    } else
    if ( Ct_BaseIso ) {
    ret |= (Ct_BaseIso->lenCt != Ct_Base->lenCt) || memcmp(Ct_BaseIso->Ctbl, Ct_Base->Ctbl, Ct_BaseIso->lenCt * sizeof(Ct_BaseIso->Ctbl[0]));
    ret |= (Ct_BaseIso->lenNumH != Ct_Base->lenNumH) || memcmp( Ct_BaseIso->NumH, Ct_Base->NumH, Ct_BaseIso->lenNumH*sizeof(Ct_BaseIso->NumH[0]));
    if ( Ct_FixHIso ) {
    ret |= (Ct_FixHIso->lenCt > Ct_BaseIso->lenCt) || memcmp(Ct_FixHIso->Ctbl, Ct_BaseIso->Ctbl, Ct_FixHIso->lenCt * sizeof(Ct_FixHIso->Ctbl[0]));
    ret |= (Ct_FixHIso->lenNumH > Ct_BaseIso->lenNumH) || memcmp( Ct_FixHIso->NumH, Ct_BaseIso->NumH, Ct_FixHIso->lenNumH*sizeof(Ct_BaseIso->NumH[0]));
    }
    }

    if ( ret ) {
        goto exit_error;
    }

    if ( bReqTaut ) {
        /* restore save "process isotopic" mark; temporarily set it to NO */
        t_group_info->bIgnoreIsotopic = bTautIgnoreIsotopic;
    }


    /* output the canonicalization results */
    pBCN->num_max          = num_max;
    pBCN->num_at_tg        = num_at_tg;
    pBCN->num_atoms        = num_atoms;

    pBCN->ftcn[TAUT_NON].NeighList          = NeighList[TAUT_NON]; NeighList[TAUT_NON] = NULL;
    pBCN->ftcn[TAUT_YES].NeighList          = NeighList[TAUT_YES]; NeighList[TAUT_YES] = NULL;

    if ( bReqTaut ) {  /* tautomeric results */
        /* base tautomeric structure, iBase = TAUT_YES */
    
        pBCN->ftcn[TAUT_YES].num_at_tg          = num_at_tg;
        pBCN->ftcn[TAUT_YES].num_atoms          = num_atoms;
                                               
        pBCN->ftcn[TAUT_YES].LinearCt           = Ct_Base->Ctbl;            Ct_Base->Ctbl        = NULL;
        pBCN->ftcn[TAUT_YES].nLenLinearCtAtOnly = s[iBase].nLenCTAtOnly;
        pBCN->ftcn[TAUT_YES].nMaxLenLinearCt    = s[iBase].nLenCT+1;
        pBCN->ftcn[TAUT_YES].nLenLinearCt       = s[iBase].nLenCT;
    
        pBCN->ftcn[TAUT_YES].PartitionCt.Rank     = nCanonRankBase;        nCanonRankBase       = NULL;
        pBCN->ftcn[TAUT_YES].PartitionCt.AtNumber = nAtomNumberCanonBase;  nAtomNumberCanonBase = NULL;                     
        pBCN->ftcn[TAUT_YES].nSymmRankCt          = nSymmRankBase;         nSymmRankBase        = NULL;

        pBCN->ftcn[TAUT_YES].nNumHOrig          = numH;                 numH          = NULL;
        pBCN->ftcn[TAUT_YES].nNumH              = Ct_Base->NumH;        Ct_Base->NumH = NULL;
        pBCN->ftcn[TAUT_YES].nLenNumH           = inchi_min(maxlenNumH, Ct_Base->maxlenNumH);

        /* fixed H structure: exists only if the structure is tautomeric */
        pBCN->ftcn[TAUT_YES].nNumHOrigFixH      = NULL;
        pBCN->ftcn[TAUT_YES].nNumHFixH          = NULL;
        pBCN->ftcn[TAUT_YES].nLenNumHFixH       = 0;
        pBCN->ftcn[TAUT_YES].nCanonFlags       |= pCD[iBase].nCanonFlags;
        
        CleanNumH( pBCN->ftcn[TAUT_YES].nNumHOrig, pBCN->ftcn[TAUT_YES].nLenNumH );
        CleanNumH( pBCN->ftcn[TAUT_YES].nNumH,     pBCN->ftcn[TAUT_YES].nLenNumH );
        CleanCt  ( pBCN->ftcn[TAUT_YES].LinearCt,  pBCN->ftcn[TAUT_YES].nLenLinearCt );

        /* isotopic canonicalization */
        if ( Ct_BaseIso ) {
            pBCN->ftcn[TAUT_YES].PartitionCtIso.Rank     = nCanonRankBaseIso;           nCanonRankBaseIso           = NULL;
            pBCN->ftcn[TAUT_YES].PartitionCtIso.AtNumber = nAtomNumberCanonBaseIso;     nAtomNumberCanonBaseIso     = NULL;
            pBCN->ftcn[TAUT_YES].nSymmRankCtIso          = nSymmRankBaseIso;            nSymmRankBaseIso            = NULL;
            pBCN->ftcn[TAUT_YES].iso_sort_keys           = Ct_BaseIso->iso_sort_key;    Ct_BaseIso->iso_sort_key    = NULL;
            pBCN->ftcn[TAUT_YES].iso_sort_keysOrig       = iso_sort_keyBase;            iso_sort_keyBase            = NULL;
            pBCN->ftcn[TAUT_YES].len_iso_sort_keys       = len_iso_sort_keyBase;
            pBCN->ftcn[TAUT_YES].iso_exchg_atnos         = Ct_BaseIso->iso_exchg_atnos; Ct_BaseIso->iso_exchg_atnos = NULL;
            pBCN->ftcn[TAUT_YES].iso_exchg_atnosOrig     = iso_exchg_atnos;             iso_exchg_atnos             = NULL;

            CleanIsoSortKeys( pBCN->ftcn[TAUT_YES].iso_sort_keys, pBCN->ftcn[TAUT_YES].len_iso_sort_keys );
            CleanIsoSortKeys( pBCN->ftcn[TAUT_YES].iso_sort_keysOrig, pBCN->ftcn[TAUT_YES].len_iso_sort_keys );
        }

    } /* tautomeric results */

    if ( bReqNonTaut  ) { /* non-tautomeric results */
        /* TAUT_NON if tautomeric + non-tautomeric or special non-taut request
           TAUT_YES if the structure happened to be non-tautomeric while user requested tautomeric processing
           In both cases the correct index is iOther. TAUT_NON replaced with iOther 4-2-2004 */
    
        if ( !bReqTaut ) {
            /* rearrange the results for a non-tautomeric structure */
            nSymmRankFixH           = nSymmRankNoTautH;           nSymmRankNoTautH        = NULL;
            nCanonRankFixH          = nCanonRankNoTautH;          nCanonRankNoTautH       = NULL;
            nAtomNumberCanonFixH    = nAtomNumberCanonNoTautH;    nAtomNumberCanonNoTautH = NULL;
            Ct_FixH                 = Ct_NoTautH;                 Ct_NoTautH              = NULL;
            /* isotopic canonicalization */
            nSymmRankFixHIso        = nSymmRankNoTautHIso;        nSymmRankNoTautHIso        = NULL;
            nCanonRankFixHIso       = nCanonRankNoTautHIso;       nCanonRankNoTautHIso       = NULL;
            nAtomNumberCanonFixHIso = nAtomNumberCanonNoTautHIso; nAtomNumberCanonNoTautHIso = NULL;
            Ct_FixHIso              = Ct_NoTautHIso;              Ct_NoTautHIso              = NULL;

            if ( iOther == TAUT_YES && pBCN->ftcn[TAUT_NON].NeighList && !pBCN->ftcn[TAUT_YES].NeighList ) {
                /* here only non-taut results go to pBCN->ftcn[TAUT_YES]
                   Since non-taut NeighList is always in pBCN->ftcn[TAUT_NON].NeighList, move it to
                   pBCN->ftcn[TAUT_YES].NeighList. 2004-04-02.
                */
                pBCN->ftcn[TAUT_YES].NeighList = pBCN->ftcn[TAUT_NON].NeighList;
                pBCN->ftcn[TAUT_NON].NeighList = NULL;
            }
        }
        pBCN->ftcn[iOther].num_at_tg          = num_atoms;
        pBCN->ftcn[iOther].num_atoms          = num_atoms;
    
        pBCN->ftcn[iOther].LinearCt           = Ct_FixH->Ctbl;       Ct_FixH->Ctbl        = NULL;
        pBCN->ftcn[iOther].nLenLinearCtAtOnly = s[iOther].nLenCTAtOnly;
        pBCN->ftcn[iOther].nMaxLenLinearCt    = s[iOther].nLenCTAtOnly+1;
        pBCN->ftcn[iOther].nLenLinearCt       = s[iOther].nLenCTAtOnly;
    
        pBCN->ftcn[iOther].PartitionCt.Rank      = nCanonRankFixH;        nCanonRankFixH       = NULL;
        pBCN->ftcn[iOther].PartitionCt.AtNumber  = nAtomNumberCanonFixH;  nAtomNumberCanonFixH = NULL;                     
        pBCN->ftcn[iOther].nSymmRankCt           = nSymmRankFixH;         nSymmRankFixH        = NULL;

        pBCN->ftcn[iOther].nNumHOrig          = numHNoTautH;              numHNoTautH          = NULL;
        pBCN->ftcn[iOther].nNumH              = Ct_FixH->NumH;            Ct_FixH->NumH        = NULL;
        pBCN->ftcn[iOther].nLenNumH           = inchi_min(maxlenNumHNoTautH,Ct_FixH->maxlenNumH);

        /* fixed H structure: exists only if the structure is tautomeric */
        pBCN->ftcn[iOther].nNumHOrigFixH      = NumHfixed;            NumHfixed          = NULL;
        pBCN->ftcn[iOther].nNumHFixH          = Ct_FixH->NumHfixed;   Ct_FixH->NumHfixed = NULL;
        pBCN->ftcn[iOther].nLenNumHFixH       = num_atoms;
        pBCN->ftcn[iOther].nCanonFlags       |= pCD[iOther].nCanonFlags;

        /* original H */
        CleanNumH( pBCN->ftcn[iOther].nNumHOrig,     pBCN->ftcn[iOther].nLenNumH );
        CleanNumH( pBCN->ftcn[iOther].nNumHOrigFixH, pBCN->ftcn[iOther].nLenNumH );
        /* canonical H positions */
        CleanNumH( pBCN->ftcn[iOther].nNumH,     pBCN->ftcn[iOther].nLenNumH );
        CleanNumH( pBCN->ftcn[iOther].nNumHFixH, pBCN->ftcn[iOther].nLenNumH );
        /* connection table */
        CleanCt( pBCN->ftcn[iOther].LinearCt, pBCN->ftcn[iOther].nLenLinearCt );

       /* isotopic canonicalization */
        if ( Ct_FixHIso ) {
            pBCN->ftcn[iOther].PartitionCtIso.Rank     = nCanonRankFixHIso;        nCanonRankFixHIso        = NULL;
            pBCN->ftcn[iOther].PartitionCtIso.AtNumber = nAtomNumberCanonFixHIso;  nAtomNumberCanonFixHIso  = NULL;
            pBCN->ftcn[iOther].nSymmRankCtIso          = nSymmRankFixHIso;         nSymmRankFixHIso         = NULL;
            pBCN->ftcn[iOther].iso_sort_keys           = Ct_FixHIso->iso_sort_key; Ct_FixHIso->iso_sort_key = NULL;
            pBCN->ftcn[iOther].iso_sort_keysOrig       = iso_sort_key_NoTautH;     iso_sort_key_NoTautH     = NULL;
            pBCN->ftcn[iOther].len_iso_sort_keys       = len_iso_sort_key_NoTautH;
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
            MergeCleanIsoSortKeys( pBCN->ftcn[iOther].iso_sort_keys, Ct_FixHIso->iso_sort_key_Hfixed, pBCN->ftcn[iOther].len_iso_sort_keys );
            MergeCleanIsoSortKeys( pBCN->ftcn[iOther].iso_sort_keysOrig, iso_sort_key_Hfixed, pBCN->ftcn[iOther].len_iso_sort_keys );
#else
            CleanIsoSortKeys( pBCN->ftcn[iOther].iso_sort_keys, pBCN->ftcn[iOther].len_iso_sort_keys );
            CleanIsoSortKeys( pBCN->ftcn[iOther].iso_sort_keysOrig, pBCN->ftcn[iOther].len_iso_sort_keys );
#endif
        }

    }  /* non-tautomeric results */
    goto exit_function;

exit_error_alloc:
    ret = CT_OUT_OF_RAM;
    goto exit_function;

exit_error:
    if ( !RETURNED_ERROR(ret) ) {
        ret = CT_CANON_ERR;
    }
    goto exit_function;

exit_function:    

#define FREE_CONTABLE(X) if (X) {CTableFree(X);inchi_free(X);}
#define FREE_ARRAY(X) if (X) inchi_free(X);

    FreeNeighList( NeighList[TAUT_NON] );
    FreeNeighList( NeighList[TAUT_YES] );

    FREE_CONTABLE( Ct_NoH )
    FREE_CONTABLE( Ct_NoTautH )
    FREE_CONTABLE( Ct_Base )
    FREE_CONTABLE( Ct_FixH )
    FREE_CONTABLE( Ct_Temp )
    /* isotopic canonicalization */
    FREE_CONTABLE( Ct_NoTautHIso )
    FREE_CONTABLE( Ct_BaseIso )
    FREE_CONTABLE( Ct_FixHIso )
    
    /* free the first two pointers from pBCN->pRankStack */
    FREE_ARRAY( nRank )
    FREE_ARRAY( nAtomNumber )
    
    if ( pBCN->pRankStack ) {
        pBCN->pRankStack[0] =
        pBCN->pRankStack[1] = NULL;
    }

#if ( USE_AUX_RANKING == 1 )
    FREE_ARRAY( nRankAux            )
    FREE_ARRAY( nAtomNumberAux      )
    FREE_ARRAY( pAtomInvariantAux   )
#endif

    FREE_ARRAY( pAtomInvariant )
    
    FREE_ARRAY( nCanonRankNoH )
    FREE_ARRAY( nAtomNumberCanonNoH )
    FREE_ARRAY( nSymmRankNoH )
    
    FREE_ARRAY( nSymmRankNoTautH )
    FREE_ARRAY( nCanonRankNoTautH )
    FREE_ARRAY( nAtomNumberCanonNoTautH )
    FREE_ARRAY( numHNoTautH )

    FREE_ARRAY( nSymmRankBase )
    FREE_ARRAY( nCanonRankBase )
    FREE_ARRAY( nAtomNumberCanonBase )
    FREE_ARRAY( numH )

    FREE_ARRAY( nSymmRankFixH )
    FREE_ARRAY( nCanonRankFixH )
    FREE_ARRAY( nAtomNumberCanonFixH )
    FREE_ARRAY( NumHfixed )

    /* isotopic canonicalization */

    FREE_ARRAY( nSymmRankNoTautHIso )
    FREE_ARRAY( nCanonRankNoTautHIso )
    FREE_ARRAY( nAtomNumberCanonNoTautHIso )
    FREE_ARRAY( iso_sort_key_NoTautH )

    FREE_ARRAY( nSymmRankBaseIso )
    FREE_ARRAY( nCanonRankBaseIso )
    FREE_ARRAY( nAtomNumberCanonBaseIso )
    FREE_ARRAY( iso_sort_keyBase )
    FREE_ARRAY( iso_exchg_atnos )
    
    FREE_ARRAY( nSymmRankFixHIso )
    FREE_ARRAY( nCanonRankFixHIso )
    FREE_ARRAY( nAtomNumberCanonFixHIso )
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    FREE_ARRAY( iso_sort_key_Hfixed )
#endif
    
    FREE_ARRAY( nTempRank )

#undef FREE_CONTABLE
#undef FREE_ARRAY

    return ret;
}
