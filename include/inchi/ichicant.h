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


#ifndef _INCHICANT_H_
#define _INCHICANT_H_


#include "ichisize.h"
#include "ichinorm.h"


/*
    Canonicalization definitions
*/

#ifndef INCHI_US_SHORT_DEF
typedef signed short S_SHORT;
typedef unsigned short U_SHORT;
#define INCHI_US_SHORT_DEF
#endif

/*typedef unsigned long  INCHI_MODE;*/

typedef union tagSplitLong {
    unsigned long  ul;
    U_SHORT        us[2];
}SU_LONG;

#define _HI 1            /* Intel platform */
#define _LO 0

#define NEIGH_LIST_LEN 4
#define U_LONG_LEN  2

#ifndef defined_NEIGH_LIST
typedef AT_RANK  *NEIGH_LIST;
#define defined_NEIGH_LIST
#endif

typedef struct tagEQUIV_INFO {
    int nNumSets;
    int *nCutVertexAtom; /* cut-vertex atom for the set of equivalent atoms */
    int *nFirstInSet;    /* first of equivalent atoms in the connected to the cut-vertex atom parts of the structure */
    int *nNumInSet;      /* number of the equivalent atoms connected to the cut-vertex atom */
    int *nAtomNo;        /* eqivalent atom number */
    int *nAddToRank;     /* number to add to the rank to normalize */
} EQUIV_INFO;

#define MOL_PART_MASK  (~0x0U ^ 0x07U)

typedef struct tagAtData_dch {
    char element[3];
    int valence;
}AT_DATA;


#define MAXVAL 20       /* maximum valence */

#define ATOM_EL_LEN 6

typedef struct tagAtomInvariantBytes {
    S_CHAR cNotExactlyHillOrderNumber;
    S_CHAR cNumberOfConnections;
    /* S_CHAR cNumberOfNonHydrogenBonds; */
    S_CHAR cAtomicNumber;
#if ( HYDROGENS_IN_INIT_RANKS == 1 )
    S_CHAR cNumberOfAttachedHydrogens;
#endif
} ATOM_INVARIANT_BYTES;

typedef struct tagAtomInvariant {
    /* non-isotopic part */
#if ( USE_DISTANCES_FOR_RANKING == 1 )
    AT_RANK         nDistanceFromTerminal;
#endif
    ATOM_INVARIANT_BYTES b;
    AT_RANK         cNum_tautomer;        /* 0 or for tautomer endpoint: number of endpoints in the group */
    AT_RANK         cNum_tautomer_num[T_NUM_NO_ISOTOPIC]; /* 0 or numbers from t_gtroup */
    /* isotopic part */
    AT_ISO_SORT_KEY iso_sort_key;
    AT_RANK         cNum_tautomer_iso[T_NUM_ISOTOPIC]; /* 0 or numbers from t_group */
} ATOM_INVARIANT;
/**********************************/
typedef enum tagAtInvariantIndexes {
    AT_INV_HILL_ORDER,
    AT_INV_NUM_CONNECTIONS,
    AT_INV_NUM_H,
    /* for endpoint + undirected graph, otherwise 0 */
    AT_INV_NUM_TG_ENDPOINTS,
    AT_INV_TG_NUMBERS,       /* num H, num (-) */
    AT_INV_NUM_H_FIX = AT_INV_TG_NUMBERS + T_NUM_NO_ISOTOPIC,
    AT_INV_BREAK1,
    /* here compare iso sort key */
    AT_INV_TAUT_ISO = AT_INV_BREAK1,
    AT_INV_LENGTH = AT_INV_TAUT_ISO + T_NUM_ISOTOPIC
} AT_INV_INDEXES;

typedef struct tagAtomInvariant2 {
    AT_NUMB         val[AT_INV_LENGTH];
    AT_ISO_SORT_KEY iso_sort_key;
    S_CHAR          iso_aux_key;
} ATOM_INVARIANT2;

/******************* Partition **********************************/
typedef struct tagPartition {
    AT_RANK *Rank;
    AT_NUMB *AtNumber;
} Partition;

/********************* BFCN *************************************/
typedef struct tagFixHOrTautCanonNumbering {

    int             num_at_tg;  /* = num_atoms for non-taut */
    int             num_atoms;
    int             nCanonFlags;
    NEIGH_LIST     *NeighList;  /* length = num_at_tg */
    /****************************/
    /*     base structure       */
    /****************************/
    AT_RANK        *LinearCt;      /* connection table atoms (+taut. groups, directed graph)*/
    int             nLenLinearCtAtOnly;
    int             nLenLinearCt;
    int             nMaxLenLinearCt;

    Partition       PartitionCt;  /* canonical numbering */
    AT_RANK        *nSymmRankCt;  /* orbits */

    /* orig. fixed by tautomerism H positions */
    NUM_H          *nNumHOrig;  /* original  H atoms positions + taut. info, excluding tautomeric H */
    NUM_H          *nNumH;      /* canonical H atoms positions + taut. info, excluding tautomeric H */
    int             nLenNumH;   /* length = num_atoms + 2*num_taut_groups */

    /* fixed H: original positions of tautomeric H; exists obly for tautomeric structures */
    NUM_H          *nNumHOrigFixH;  /* original fixed positions of tautomeric H */
    NUM_H          *nNumHFixH;      /* canonical fixed positions of tautomeric H */
    int             nLenNumHFixH;   /* length = num_atoms */

    /*******************************************************************************/
    /* the following exists only if isotopic and isotopic results requested        */
    /*******************************************************************************/
    Partition       PartitionCtIso;     /* canonical numbering of isotopic base structure, defined later */
    AT_RANK        *nSymmRankCtIso;     /* orbits of isotopic structure */
    AT_ISO_SORT_KEY *iso_sort_keys;     /* original isotopic sort keys for atoms and taut groups */
    AT_ISO_SORT_KEY *iso_sort_keysOrig; /* canonical isotopic sort keys for atoms and taut groups */
    int              len_iso_sort_keys;
    S_CHAR          *iso_exchg_atnos;     /* canonical: 0=> tautomeric or may have isotopic H exchanged */
    S_CHAR          *iso_exchg_atnosOrig; /* original: 0=> tautomeric or may have isotopic H exchanged */
} FTCN;

/******************** BCN *************************************/
typedef struct tagBaseCanonNumbering {

    AT_RANK            **pRankStack;
    int                  nMaxLenRankStack;
    int                  num_max;        /* allocated nRank[] arrays lengths in pRankStack */
    int                  num_at_tg;  /* all of the following arrays have this length */
    int                  num_atoms;
    struct tagInchiTime *ulTimeOutTime;
    FTCN                 ftcn[TAUT_NUM];
} BCN;

/***********************************
 *
 *  CANON_STAT
 */
typedef struct tagCanonStat {
    /*  statistics */
    long                 lNumBreakTies;
    long                 lNumNeighListIter;
    long                 lNumTotCT;
    long                 lNumDecreasedCT;
    long                 lNumRejectedCT;
    long                 lNumEqualCT;
    struct tagInchiTime *ulTimeOutTime;
    long                 lTotalTime;

    /* control */
    int                  bFirstCT;
    int                  bKeepSymmRank;
    int                  bStereoIsBetter;

    int nCanonFlags;

    /* data : */

    AT_NUMB          *LinearCT;        /* connection table only */
    AT_ISOTOPIC      *LinearCTIsotopic;
    AT_ISO_TGROUP    *LinearCTIsotopicTautomer;
    AT_STEREO_DBLE   *LinearCTStereoDble;
    AT_STEREO_CARB   *LinearCTStereoCarb;
    AT_STEREO_DBLE   *LinearCTStereoDbleInv;
    AT_STEREO_CARB   *LinearCTStereoCarbInv;
    AT_STEREO_DBLE   *LinearCTIsotopicStereoDble;
    AT_STEREO_CARB   *LinearCTIsotopicStereoCarb;
    AT_STEREO_DBLE   *LinearCTIsotopicStereoDbleInv;
    AT_STEREO_CARB   *LinearCTIsotopicStereoCarbInv;
    AT_TAUTOMER      *LinearCTTautomer;  /*  minimal */

/* second copies of line notation arrays */

    AT_NUMB          *LinearCT2;   /* to save non-isotopic CT */

    int               nLenLinearCTStereoDble;
    int               nLenLinearCTStereoDbleInv;
    int               nMaxLenLinearCTStereoDble;  /* new */

    int               bCmpStereo;         /* 0 => no stereo to invert;
                                             1 => StereoCtInv < StereoCt;
                                             2 => StereoCtInv = StereoCt;
                                             3 => StereoCtInv > StereoCt;
                                           */
    int               nLenLinearCTStereoCarb;
    int               nLenLinearCTStereoCarbInv;
    int               nMaxLenLinearCTStereoCarb;  /* new */

    int               nLenLinearCTIsotopic;
    int               nMaxLenLinearCTIsotopic;

    int               nLenLinearCTIsotopicTautomer;
    int               nMaxLenLinearCTIsotopicTautomer;

    int               nLenLinearCT;         /* connection table only  */
    int               nLenLinearCT2;        /* connection table only, non-isotopic result  */
    int               nLenLinearCTAtOnly;   /* connection table only without tautomeric pseudoatoms  */
    int               nLenLinearCTAtOnly2;  /* connection table only, non-isotopic result without tautomeric pseudoatoms  */
    int               nMaxLenLinearCT;      /* connection table only  */

    int               nLenLinearCTTautomer;
    int               nMaxLenLinearCTTautomer;

    int               bCmpIsotopicStereo; /* 0 => no stereo to invert;
                                             1 => StereoCtInv < StereoCt;
                                             2 => StereoCtInv = StereoCt;
                                             3 => StereoCtInv > StereoCt;
                                           */
    int               nLenLinearCTIsotopicStereoDble;
    int               nLenLinearCTIsotopicStereoDbleInv;
    int               nMaxLenLinearCTIsotopicStereoDble;

    int               nLenLinearCTIsotopicStereoCarb; /*  new */
    int               nLenLinearCTIsotopicStereoCarbInv; /*  new */
    int               nMaxLenLinearCTIsotopicStereoCarb;
    S_CHAR           *bRankUsedForStereo;  /* canon. rank used for stereo mapping */
    S_CHAR           *bAtomUsedForStereo;  /* 0 if not a stereo atom or during a canon. rank being mapped on this atom; */
                                           /* STEREO_AT_MARK if an unpapped stereogenic atom */
                                           /* or a number of stereogenic bonds adjacent to an atom */

    AT_RANK          *nPrevAtomNumber;

    AT_RANK          *nCanonOrd;       /* atom numbers in order of increasing canon. ranks  */
    AT_RANK          *nSymmRank;       /* symmetry numbers in order of atoms  */
    AT_RANK          *nCanonOrdTaut;   /* t-group numbers numbers in order of increasing canon. ranks  */
    AT_RANK          *nSymmRankTaut;   /* t-group symmetry numbers in order of t-groups  */

    AT_RANK          *nCanonOrdStereo;     /* atom numbers in order of increasing canon. ranks */
    AT_RANK          *nCanonOrdStereoInv;     /* atom numbers in order of increasing canon. ranks */
    AT_RANK          *nCanonOrdStereoTaut; /* t-group numbers in order of increasing canon. ranks */

    AT_RANK          *nSymmRankIsotopic;
    AT_RANK          *nCanonOrdIsotopic;        /* atom numbers in order of increasing canon. ranks */
    AT_RANK          *nSymmRankIsotopicTaut;    /* !!! */
    AT_RANK          *nCanonOrdIsotopicTaut;    /*/ t-group numbers in order of increasing canon. ranks */

    AT_RANK          *nCanonOrdIsotopicStereo;
    AT_RANK          *nCanonOrdIsotopicStereoInv;
    AT_RANK          *nCanonOrdIsotopicStereoTaut;    /*  !!! */

                      /* actual lengths if successfully calculated */

    int               nLenCanonOrd;               /* Superceded by any of the following > 0 */
    int               nLenCanonOrdTaut;           /* !!! Superceded by any of the following > 0 */
    int               nLenCanonOrdIsotopic;
    int               nLenCanonOrdIsotopicTaut;   /* !!! */
    int               nLenCanonOrdStereo;
    int               nLenCanonOrdStereoTaut;     /* !!! */
    int               nLenCanonOrdIsotopicStereo;
    int               nLenCanonOrdIsotopicStereoTaut; /* !!! */

                      /*  other */

    int               bHasIsotopicInTautomerGroups;
    T_GROUP_INFO     *t_group_info;
    int               bIgnoreIsotopic;
    int               bDoubleBondSquare; /* 0 or 2 */
    INCHI_MODE         nMode;
#if ( bRELEASE_VERSION == 0 )
    int               bExtract;          /* for debug only */
#endif
    NEIGH_LIST       *NeighList;
    BCN              *pBCN;
    S_CHAR    *nNum_H;      /* number of terminal hydrogen atoms on each atom except tautomeric [num_atoms], in order of canonical numbers */
    S_CHAR    *nNum_H_fixed;/* number of terminal hydrogen atoms on tautomeric atoms (for non-atautomeric representation) [num_atoms] */
    S_CHAR    *nExchgIsoH;
} CANON_STAT;




typedef struct tagCANON_GLOBALS
{
    const NEIGH_LIST      *m_pNeighList_RankForSort;
    const ATOM_INVARIANT2 *m_pAtomInvariant2ForSort;
    const AT_NUMB         *m_pNeighborsForSort;
    const AT_RANK         *m_pn_RankForSort;
    AT_RANK m_nMaxAtNeighRankForSort;
    int m_nNumCompNeighborsRanksCountEql;
    bitWord *m_bBit;
    int m_bBitInitialized;
    int m_num_bit;
} CANON_GLOBALS;

int  SetBitCreate( struct tagCANON_GLOBALS *pCG );

void inchi_qsort( void *pParam, void *base, size_t num, size_t width, int( *comp )( const void *, const void *, void * ) );





/**************************************************/
typedef struct tagCanonData {

    /* same names/types as in ConTable; here the order is from original numbering */

    AT_NUMB *LinearCT;  /* output ?? */

    int      nMaxLenLinearCT;
    int      nLenLinearCT;
    int      nLenCTAtOnly;
    int      nCanonFlags;
    /* hydrogen atoms fixed in tautomeric representation:
       compare before diff sign inversion: (+) <=> Ct1->() > Ct2->() */
    NUM_H          *NumH;
    int             lenNumH;    /* used length */
    int             maxlenNumH; /*  n + T_NUM_NO_ISOTOPIC*(n_tg-n) + 1 */

    /* hydrogen atoms fixed in non-tautomeric representation only:
       compare before diff sign inversion: (+) <=> Ct1->() > Ct2->() */
    NUM_H           *NumHfixed;
    int              lenNumHfixed;       /* used length */
    int              maxlenNumHfixed;    /* max length = n+1  */

    /* isotopic atoms (without tautomeric H) and isotopic tautomeric groups */
    /* note: AT_ISO_SORT_KEY and T_GROUP_ISOWT are identical types: long    */
    AT_ISO_SORT_KEY *iso_sort_key;
    int              len_iso_sort_key;    /* used length */
    int              maxlen_iso_sort_key; /* max length = n_tg+1 */
    S_CHAR          *iso_exchg_atnos;
    int              len_iso_exchg_atnos;    /* used length */
    int              maxlen_iso_exchg_atnos;

    /* isotopic hydrogen atoms fixed in non-tautomeric representation only */
#if ( USE_ISO_SORT_KEY_HFIXED == 1 )
    AT_ISO_SORT_KEY *iso_sort_key_Hfixed;
    int              len_iso_sort_key_Hfixed;    /* used length */
    int              maxlen_iso_sort_key_Hfixed; /* max length = n+1  */
#endif
    /* auxiliary ranking */

    AT_RANK  *nAuxRank;

    struct tagInchiTime *ulTimeOutTime;  /* timeout */
} CANON_DATA;
/**************************************************/

typedef struct tagCanonCounts {
    long     lNumBreakTies;
    long     lNumDecreasedCT;
    long     lNumRejectedCT;
    long     lNumEqualCT;
    long     lNumTotCT;
    double    dGroupSize;
    long     lNumGenerators;
    long     lNumStoredIsomorphisms;
} CANON_COUNTS;
/***********************************************
 tree structure: one segment

   canon. rank
   at.no          orig. atom numbers on which the canon. rank has been successfully mapped
   ...
   at.no          except the last at.no: it is not known if it has been mapped until all atoms are mapped
   num.at+1       number of atoms in this segment
*/

typedef struct tagCurTree {
    AT_NUMB   *tree;
    int       max_len;  /* allocated length of tree in sizeof(tree[0]) units */
    int       cur_len;  /* currently used length */
    int       incr_len; /* reallocation increment */
} CUR_TREE;


#endif /* _INCHICANT_H_ */
