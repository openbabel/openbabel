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


#ifndef __INCHICOMN_H__
#define __INCHICOMN_H__

#include "ichierr.h"


/********************************
 *
 * Globals for sorting
 */

extern const NEIGH_LIST      *pNeighList_RankForSort; 
extern const ATOM_INVARIANT2 *pAtomInvariant2ForSort;
extern const AT_NUMB         *pNeighborsForSort;
extern const AT_RANK         *pn_RankForSort;

extern AT_RANK         nMaxAtNeighRankForSort;

extern int             nNumCompNeighborsRanksCountEql;


#define tsort insertions_sort

typedef struct tagEquNeigh {
    int     num_to;   /* number of neighbors with equal mapping ranks; one of them has min. canon. number */
    AT_RANK to_at[4]; /* to_atom neighbors #s with equal mapping ranks */
    AT_RANK from_at;  /* from_at neighbor # which has min. canon. number and can be mapped on any of the above to_at[] */
    AT_RANK rank;     /* equal mapping rank value */
    AT_RANK canon_rank;  /* min. canon. number */
} EQ_NEIGH;

#ifndef INCHI_ALL_CPP
#ifdef  __cplusplus
extern "C" {
#endif
#endif

/*******************************************************************/
/* ichiisot.c */
int unpack_iso_sort_key( AT_ISO_SORT_KEY iso_sort_key, S_CHAR *num_1H, S_CHAR *num_2H, S_CHAR *num_3H, S_CHAR *iso_atw_diff );
AT_ISO_SORT_KEY make_iso_sort_key( int iso_atw_diff, int num_1H, int num_2H, int num_3H);
int set_atom_iso_sort_keys(  int num_at, sp_ATOM *at, T_GROUP_INFO* t_group_info, int *bHasIsotopicInTautomerGroups );
/***********************************************************************/
/* ichisort.c */

void insertions_sort_NeighList_AT_NUMBERS( NEIGH_LIST base, AT_RANK *nRank );
int insertions_sort_NeighList_AT_NUMBERS3( NEIGH_LIST base, AT_RANK *nRank );
int insertions_sort_AT_RANK( AT_RANK *base, int num );
void insertions_sort_NeighListBySymmAndCanonRank( NEIGH_LIST base, const AT_RANK *nSymmRank, const AT_RANK *nCanonRank );
int CompareNeighListLex( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank);
int CompareNeighListLexUpToMaxRank( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank, AT_RANK nMaxAtNeighRank );
int compare_NeighLists( const NEIGH_LIST *op1,  const NEIGH_LIST *op2 );
int CompNeighborsAT_NUMBER( const void* a1, const void* a2);
int comp_AT_RANK( const void* a1, const void* a2);
int CompRank(const void* a1, const void* a2 );
int CompRanksOrd( const void* a1, const void* a2 );

int CompAtomInvariants2Only( const void* a1, const void* a2 );
int CompAtomInvariants2( const void* a1, const void* a2 );

int CompNeighListRanks( const void* a1, const void* a2 );
int CompNeighListRanksOrd( const void* a1, const void* a2 );
int CompNeighLists( const void* a1, const void* a2 );
int CompNeighListsUpToMaxRank( const void* a1, const void* a2 );
int CompNeighborsRanksCountEql( const void* a1, const void* a2 );
int CompRanksInvOrd( const void* a1, const void* a2 );
int CompChemElemLex( const void *a1, const void *a2 );



NEIGH_LIST *CreateNeighList( int num_atoms, int num_at_tg, sp_ATOM* at, int bDoubleBondSquare, T_GROUP_INFO *t_group_info );
NEIGH_LIST *CreateNeighListFromLinearCT( AT_NUMB *LinearCT, int nLenCT, int num_atoms );

void FreeNeighList( NEIGH_LIST *pp );
int BreakAllTies( int num_atoms, int num_max, AT_RANK **pRankStack,
                     NEIGH_LIST *NeighList, AT_RANK *nTempRank, CANON_STAT *pCS);

/******************************************************************************/
/* ichimap.c */
void switch_ptrs( AT_RANK **p1, AT_RANK **p2 );

int SortedEquInfoToRanks( const AT_RANK* nSymmRank, AT_RANK* nRank, const AT_RANK* nAtomNumber, int num_atoms, int *bChanged );
int SortedRanksToEquInfo( AT_RANK* nSymmRank, const AT_RANK* nRank, const AT_RANK* nAtomNumber, int num_atoms );

int SetNewRanksFromNeighLists( int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank,
                             AT_RANK *nAtomNumber, int bUseAltSort, int ( *comp )(const void *, const void *) );
int SetNewRanksFromNeighLists3( int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank,
                             AT_RANK *nAtomNumber );
int SetNewRanksFromNeighLists4( int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank,
                             AT_RANK *nAtomNumber, AT_RANK nMaxAtRank );
void SortNeighListsBySymmAndCanonRank( int num_atoms, NEIGH_LIST *NeighList, const AT_RANK *nSymmRank, const AT_RANK *nCanonRank );
int SortNeighLists2( int num_atoms, AT_RANK *nRank, NEIGH_LIST *NeighList, AT_RANK *nAtomNumber );
int  DifferentiateRanks2( int num_atoms, NEIGH_LIST *NeighList,
                                 int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                                 AT_RANK *nAtomNumber, long *lNumIter, int bUseAltSort );

int  DifferentiateRanks3( int num_atoms, NEIGH_LIST *NeighList,
                                 int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                                 AT_RANK *nAtomNumber, long *lNumIter );
int  DifferentiateRanks4( int num_atoms, NEIGH_LIST *NeighList,
                          int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                          AT_RANK *nAtomNumber, AT_RANK nMaxAtRank, long *lNumIter );
int  DifferentiateRanksBasic( int num_atoms, NEIGH_LIST *NeighList,
                                 int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                                 AT_RANK *nAtomNumber, long *lNumIter, int bUseAltSort );

int parity_of_mapped_atom2( int from_at, int to_at, const sp_ATOM *at, EQ_NEIGH *pEN,
                           const AT_RANK *nCanonRankFrom, const AT_RANK *nRankFrom, const AT_RANK *nRankTo );

int parity_of_mapped_half_bond( int from_at, int to_at, int from_neigh, int to_neigh,
                           sp_ATOM *at, EQ_NEIGH *pEN,
                           const AT_RANK *nCanonRankFrom, const AT_RANK *nRankFrom, const AT_RANK *nRankTo );

int HalfStereoBondParity( sp_ATOM *at, int at_no1, int i_sb_neigh, const AT_RANK *nRank );

int NumberOfTies( AT_RANK **pRankStack1, AT_RANK **pRankStack2, int length, int at_no1, int at_no2, AT_RANK *nNewRank, int *bAddStack, int *bMapped1 );

int map_an_atom2( int num_atoms, int num_max, int at_no1/*from*/, int at_no2/*to*/,
                AT_RANK *nTempRank,
                int nNumMappedRanks, int *pnNewNumMappedRanks,
                CANON_STAT *pCS,
                NEIGH_LIST    *NeighList,
                AT_RANK  **pRankStack1, AT_RANK  **pRankStack2, int *bAddStack );
int ClearPreviousMappings( AT_RANK **pRankStack1 );

int SetOneStereoBondIllDefParity( sp_ATOM *at, int jc, /* atom number*/ int k /* stereo bond ord. number*/, int new_parity );
int RemoveOneStereoBond( sp_ATOM *at, int jc, /* atom number*/ int k /* stereo bond number*/ );
int RemoveOneStereoCenter( sp_ATOM *at, int jc /* atom number*/ );
int RemoveCalculatedNonStereo( sp_ATOM *at, int num_atoms, int num_at_tg,
                              AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList,
                              const AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_RANK *nAtomNumberCanon, CANON_STAT *pCS );


int might_change_other_atom_parity( sp_ATOM *at, int num_atoms, int at_no, AT_RANK *nRank2, AT_RANK *nRank1 );

int map_stereo_bonds4 ( 
                sp_ATOM *at, int num_atoms, int num_at_tg, int num_max, int bAllene,
                const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, /*  non-stereo canon ranking */
                AT_RANK *nCanonRankTo, /* output canonical numbering*/
                const AT_RANK *nSymmRank,      AT_RANK **pRankStack1/*from*/,  AT_RANK **pRankStack2/*to*/,
                AT_RANK *nTempRank,      int nNumMappedRanksInput,
                AT_RANK *nSymmStereo,    NEIGH_LIST *NeighList,
                CANON_STAT *pCS,         CUR_TREE *cur_tree,  int nNumMappedBonds );

int map_stereo_atoms4 ( 
                sp_ATOM *at, int num_atoms, int num_at_tg, int num_max,
                const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, AT_RANK *nCanonRankTo, /*  canonical numbering to be mapped */
                const AT_RANK *nSymmRank,      AT_RANK **pRankStack1/*from*/, AT_RANK **pRankStack2/*to*/,
                AT_RANK *nTempRank,      int nNumMappedRanksInput,
                AT_RANK *nSymmStereo,    NEIGH_LIST *NeighList,
                CANON_STAT *pCS,         CUR_TREE *cur_tree, int nNumMappedAtoms );


int CurTreeAlloc( CUR_TREE *cur_tree, int num_atoms );
int CurTreeReAlloc( CUR_TREE *cur_tree );
void CurTreeFree( CUR_TREE *cur_tree );
int CurTreeAddRank( CUR_TREE *cur_tree, AT_NUMB rank );
int CurTreeRemoveLastRank( CUR_TREE *cur_tree );
int CurTreeReplaceLastRank( CUR_TREE *cur_tree, AT_NUMB rank );
int CurTreeFindTheRankPos( CUR_TREE *cur_tree, AT_NUMB rank );
int CurTreeGetPos( CUR_TREE *cur_tree );
int CurTreeSetPos( CUR_TREE *cur_tree, int len );
int CurTreeAddAtom( CUR_TREE *cur_tree, int at_no );
int CurTreeRemoveLastAtom( CUR_TREE *cur_tree );
int CurTreeIsLastRank( CUR_TREE *cur_tree, AT_NUMB rank );
int CurTreeIsLastAtomEqu( CUR_TREE *cur_tree, int at_no, AT_NUMB *nSymmStereo );
int CurTreeRemoveIfLastAtom( CUR_TREE *cur_tree, int at_no );
int CurTreeRemoveLastRankIfNoAtoms( CUR_TREE *cur_tree );
void CurTreeKeepLastAtomsOnly( CUR_TREE *cur_tree, int tpos, int shift );


void SetUseAtomForStereo( S_CHAR *bAtomUsedForStereo, sp_ATOM *at, int num_atoms );

int     nJoin2Mcrs( AT_RANK *nEqArray, AT_RANK n1, AT_RANK n2 );
AT_RANK nGetMcr( AT_RANK *nEqArray, AT_RANK n );
int     bUniqueAtNbrFromMappingRank( AT_RANK **pRankStack, AT_RANK nAtRank, AT_NUMB *nAtNumber );

int Next_SB_At_CanonRanks2( AT_RANK *canon_rank1, AT_RANK *canon_rank2, /*  canonical numbers */
                          AT_RANK *canon_rank1_min, AT_RANK *canon_rank2_min,
                          int *bFirstTime, S_CHAR *bAtomUsedForStereo,
                          const ppAT_RANK pRankStack1, const ppAT_RANK pRankStack2,
                          const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom,
                          const sp_ATOM *at, int num_atoms, int bAllene );

int Next_SC_At_CanonRank2( AT_RANK *canon_rank1,    /*  1st call input: largest canon number mapped so far or 0 */
                                                    /*  output: suggested canon. rank > than input if success */
                          AT_RANK *canon_rank1_min, /*  1st call:0 next calls: first tried canon. number */
                          int *bFirstTime,          /*  1 at the time of the 1st call  */
                          S_CHAR *bAtomUsedForStereo, /*  STEREO_AT_MARK if the atom has not been mapped yet */
                    const ppAT_RANK pRankStack1,    /*  mapping ranks/sort order of atoms with canon. numbers (from) */
                    const ppAT_RANK pRankStack2,    /*  mapping ranks/sort order of atoms with stereo (to) */
                    const AT_RANK *nAtomNumberCanonFrom, /*  sorted order of the canon. numbers */
                    int num_atoms );

int NextStereoParity2Test( int *stereo_bond_parity, int *sb_parity_calc,
                         int nNumBest, int nNumWorse, int nNumUnkn, int nNumUndf, int nNumCalc);

int All_SB_Same(  AT_RANK canon_rank1, AT_RANK canon_rank2, /*  canonical numbers */
                  const ppAT_RANK pRankStack1, const ppAT_RANK pRankStack2,
                  const AT_RANK *nAtomNumberCanonFrom,
                  sp_ATOM *at );

int All_SC_Same(  AT_RANK canon_rank1, /*  canonical number */
                  const ppAT_RANK pRankStack1, const ppAT_RANK pRankStack2,
                  const AT_RANK *nAtomNumberCanonFrom,
                  const sp_ATOM *at );


int CompareLinCtStereoDoubleToValues( AT_STEREO_DBLE *LinearCTStereoDble,
                              AT_RANK at_rank_canon1, AT_RANK at_rank_canon2, U_CHAR bond_parity );

int CompareLinCtStereoAtomToValues( AT_STEREO_CARB *LinearCTStereoCarb,
                              AT_RANK at_rank_canon1, U_CHAR parity );
int CompareLinCtStereoDble ( AT_STEREO_DBLE *LinearCTStereoDble1, int nLenLinearCTStereoDble1,
                             AT_STEREO_DBLE *LinearCTStereoDble2, int nLenLinearCTStereoDble2 );
int CompareLinCtStereoCarb ( AT_STEREO_CARB *LinearCTStereoCarb1, int nLenLinearCTStereoCarb1,
                             AT_STEREO_CARB *LinearCTStereoCarb2, int nLenLinearCTStereoCarb2 );
int CompareLinCtStereo ( AT_STEREO_DBLE *LinearCTStereoDble1, int nLenLinearCTStereoDble1,
                         AT_STEREO_CARB *LinearCTStereoCarb1, int nLenLinearCTStereoCarb1,
                         AT_STEREO_DBLE *LinearCTStereoDble2, int nLenLinearCTStereoDble2,
                         AT_STEREO_CARB *LinearCTStereoCarb2, int nLenLinearCTStereoCarb2 );
/***************************************************************************/
/* ichicans.c */
int UnmarkNonStereo( sp_ATOM *at, int num_atoms, const AT_RANK *nRank, const AT_RANK *nAtomNumber, int bIsotopic );
int FillSingleStereoDescriptors(sp_ATOM *at, int i, int num_trans, const AT_RANK *nRank
                          , AT_STEREO_CARB *LinearCTStereoCarb, int *nStereoCarbLen, int nMaxStereoCarbLen
                          , AT_STEREO_DBLE *LinearCTStereoDble, int *nStereoDbleLen, int nMaxStereoDbleLen
                          , int bAllene );
void SwitchAtomStereoAndIsotopicStereo( sp_ATOM *at, int num_atoms, int *bSwitched );
void SetCtToIsotopicStereo( CANON_STAT *pCS, CANON_STAT *pCS2 );
void SetCtToNonIsotopicStereo( CANON_STAT *pCS, CANON_STAT *pCS2 );
int FillAllStereoDescriptors( sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon, CANON_STAT *pCS );
int FillOutStereoParities( sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon,
                           const AT_RANK *nRank, const AT_RANK *nAtomNumber, CANON_STAT *pCS, int bIsotopic );
int InvertStereo( sp_ATOM *at, int num_at_tg,
                  AT_RANK *nCanonRank, AT_RANK *nAtomNumberCanon,
                  CANON_STAT *pCS, int bInvertLinearCTStereo );
int find_atoms_with_parity( sp_ATOM *at, S_CHAR *visited, int from_atom, int cur_atom );
int GetStereoNeighborPos( sp_ATOM *at, int iAt1, int iAt2 );
int GetStereoBondParity(sp_ATOM *at, int i, int n, AT_RANK *nRank );
int GetStereoCenterParity(sp_ATOM *at, int i, AT_RANK *nRank );
int GetPermutationParity( sp_ATOM *at, AT_RANK nAvoidNeighbor, AT_RANK *nCanonRank );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __INCHICOMN_H__ */
