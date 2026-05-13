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


#ifndef _INCHICOMN_H_
#define _INCHICOMN_H_


#include "ichierr.h"
#include "ichicano.h"

/*
    Globals for sorting
 */

#define tsort insertions_sort

typedef struct tagEquNeigh {
    int     num_to;   /* number of neighbors with equal mapping ranks; one of them has min. canon. number */
    AT_RANK to_at[4]; /* to_atom neighbors #s with equal mapping ranks */
    AT_RANK from_at;  /* from_at neighbor # which has min. canon. number and can be mapped on any of the above to_at[] */
    AT_RANK rank;     /* equal mapping rank value */
    AT_RANK canon_rank;  /* min. canon. number */
} EQ_NEIGH;

struct tagINCHI_CLOCK;
#ifndef COMPILE_ALL_CPP
#ifdef  __cplusplus
extern "C" {
#endif
#endif


    struct tagINCHI_CLOCK;



    /*******************************************************************/
    /* ichiisot.c */
    int unpack_iso_sort_key( AT_ISO_SORT_KEY iso_sort_key, S_CHAR *num_1H, S_CHAR *num_2H, S_CHAR *num_3H, S_CHAR *iso_atw_diff );
    AT_ISO_SORT_KEY make_iso_sort_key( int iso_atw_diff, int num_1H, int num_2H, int num_3H );
    int set_atom_iso_sort_keys( int num_at, sp_ATOM *at, T_GROUP_INFO* t_group_info, int *bHasIsotopicInTautomerGroups );
    /***********************************************************************/
    /* ichisort.c */

    void insertions_sort_NeighList_AT_NUMBERS( NEIGH_LIST base, AT_RANK *nRank );
    int insertions_sort_NeighList_AT_NUMBERS3( NEIGH_LIST base, AT_RANK *nRank );
    int insertions_sort_AT_RANK( AT_RANK *base, int num );
    void insertions_sort_NeighListBySymmAndCanonRank( NEIGH_LIST base, const AT_RANK *nSymmRank, const AT_RANK *nCanonRank );
    int CompareNeighListLex( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank );
    int CompareNeighListLexUpToMaxRank( NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank, AT_RANK nMaxAtNeighRank );
    int compare_NeighLists( const NEIGH_LIST *op1, const NEIGH_LIST *op2, void *p );
    int CompNeighborsAT_NUMBER( const void* a1, const void* a2, void *p );
    int comp_AT_RANK( const void* a1, const void* a2, void * );
    int CompRank( const void* a1, const void* a2, void *p );
    int CompRanksOrd( const void* a1, const void* a2, void *p );

    int CompAtomInvariants2Only( const void* a1, const void* a2, void *p );
    int CompAtomInvariants2( const void* a1, const void* a2, void *p );

    int CompNeighListRanks( const void* a1, const void* a2, void *p );
    int CompNeighListRanksOrd( const void* a1, const void* a2, void *p );
    int CompNeighLists( const void* a1, const void* a2, void *p );
    int CompNeighListsUpToMaxRank( const void* a1, const void* a2, void *p );
    int CompNeighborsRanksCountEql( const void* a1, const void* a2, void *p );
    int CompRanksInvOrd( const void* a1, const void* a2, void * );
    int CompChemElemLex( const void *a1, const void *a2 );
    int* iisort( int *list, int num );


    NEIGH_LIST *CreateNeighList( int num_atoms, int num_at_tg, sp_ATOM* at, int bDoubleBondSquare, T_GROUP_INFO *t_group_info );
    NEIGH_LIST *CreateNeighListFromLinearCT( AT_NUMB *LinearCT, int nLenCT, int num_atoms );

    void FreeNeighList( NEIGH_LIST *pp );
    int BreakAllTies( CANON_GLOBALS *pCG, int num_atoms, int num_max, AT_RANK **pRankStack,
                         NEIGH_LIST *NeighList, AT_RANK *nTempRank, CANON_STAT *pCS );

    /******************************************************************************/
    /* ichimap.c */
    void switch_ptrs( AT_RANK **p1, AT_RANK **p2 );

    int SortedEquInfoToRanks( const AT_RANK* nSymmRank, AT_RANK* nRank, const AT_RANK* nAtomNumber, int num_atoms, int *bChanged );
    int SortedRanksToEquInfo( AT_RANK* nSymmRank, const AT_RANK* nRank, const AT_RANK* nAtomNumber, int num_atoms );

    int SetNewRanksFromNeighLists( CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank,
                                 AT_RANK *nAtomNumber, int bUseAltSort, int( *comp )( const void *, const void *, void * ) );
    int SetNewRanksFromNeighLists3( CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank,
                                 AT_RANK *nAtomNumber );
    int SetNewRanksFromNeighLists4( CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank,
                                 AT_RANK *nAtomNumber, AT_RANK nMaxAtRank );
    void SortNeighListsBySymmAndCanonRank( int num_atoms, NEIGH_LIST *NeighList, const AT_RANK *nSymmRank, const AT_RANK *nCanonRank );
    int SortNeighLists2( int num_atoms, AT_RANK *nRank, NEIGH_LIST *NeighList, AT_RANK *nAtomNumber );
    int  DifferentiateRanks2( CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList,
                                     int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                                     AT_RANK *nAtomNumber, long *lNumIter, int bUseAltSort );

    int  DifferentiateRanks3( CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList,
                                     int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                                     AT_RANK *nAtomNumber, long *lNumIter );
    int  DifferentiateRanks4( CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList,
                              int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                              AT_RANK *nAtomNumber, AT_RANK nMaxAtRank, long *lNumIter );
    int  DifferentiateRanksBasic( CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList,
                                     int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank,
                                     AT_RANK *nAtomNumber, long *lNumIter, int bUseAltSort );

    int parity_of_mapped_atom2( CANON_GLOBALS *pCG, int from_at, int to_at, const sp_ATOM *at, EQ_NEIGH *pEN,
                               const AT_RANK *nCanonRankFrom, const AT_RANK *nRankFrom, const AT_RANK *nRankTo );

    int parity_of_mapped_half_bond( int from_at, int to_at, int from_neigh, int to_neigh,
                               sp_ATOM *at, EQ_NEIGH *pEN,
                               const AT_RANK *nCanonRankFrom, const AT_RANK *nRankFrom, const AT_RANK *nRankTo );

    int HalfStereoBondParity( sp_ATOM *at, int at_no1, int i_sb_neigh, const AT_RANK *nRank );

    int NumberOfTies( AT_RANK **pRankStack1, AT_RANK **pRankStack2, int length, int at_no1, int at_no2, AT_RANK *nNewRank, int *bAddStack, int *bMapped1 );

    int map_an_atom2( CANON_GLOBALS *pCG, int num_atoms, int num_max, int at_no1/*from*/, int at_no2/*to*/,
                    AT_RANK *nTempRank,
                    int nNumMappedRanks, int *pnNewNumMappedRanks,
                    CANON_STAT *pCS,
                    NEIGH_LIST    *NeighList,
                    AT_RANK  **pRankStack1, AT_RANK  **pRankStack2, int *bAddStack );
    int ClearPreviousMappings( AT_RANK **pRankStack1 );

    int SetOneStereoBondIllDefParity( sp_ATOM *at, int jc, /* atom number*/ int k /* stereo bond ord. number*/, int new_parity );
    int RemoveOneStereoBond( sp_ATOM *at, int jc, /* atom number*/ int k /* stereo bond number*/ );
    int RemoveOneStereoCenter( sp_ATOM *at, int jc /* atom number*/ );
    int RemoveCalculatedNonStereo( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg,
                                  AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList,
                                  const AT_RANK *nSymmRank, AT_RANK *nCanonRank,
                                  AT_RANK *nAtomNumberCanon, CANON_STAT *pCS,
                                  int vABParityUnknown );


    int might_change_other_atom_parity( sp_ATOM *at, int num_atoms, int at_no, AT_RANK *nRank2, AT_RANK *nRank1 );

    int map_stereo_bonds4(
                    struct tagINCHI_CLOCK *ic, CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg, int num_max, int bAllene,
                    const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, /*  non-stereo canon ranking */
                    AT_RANK *nCanonRankTo, /* output canonical numbering*/
                    const AT_RANK *nSymmRank, AT_RANK **pRankStack1/*from*/, AT_RANK **pRankStack2/*to*/,
                    AT_RANK *nTempRank, int nNumMappedRanksInput,
                    AT_RANK *nSymmStereo, NEIGH_LIST *NeighList,
                    CANON_STAT *pCS, CUR_TREE *cur_tree, int nNumMappedBonds,
                    int vABParityUnknown );

    int map_stereo_atoms4( struct tagINCHI_CLOCK *ic, CANON_GLOBALS *pCG,
                    sp_ATOM *at, int num_atoms, int num_at_tg, int num_max,
                    const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, AT_RANK *nCanonRankTo, /*  canonical numbering to be mapped */
                    const AT_RANK *nSymmRank, AT_RANK **pRankStack1/*from*/, AT_RANK **pRankStack2/*to*/,
                    AT_RANK *nTempRank, int nNumMappedRanksInput,
                    AT_RANK *nSymmStereo, NEIGH_LIST *NeighList,
                    CANON_STAT *pCS, CUR_TREE *cur_tree, int nNumMappedAtoms,
                    int vABParityUnknown );

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
                             int nNumBest, int nNumWorse, int nNumUnkn, int nNumUndf, int nNumCalc,
                             int vABParityUnknown );

    int All_SB_Same( AT_RANK canon_rank1, AT_RANK canon_rank2, /*  canonical numbers */
                      const ppAT_RANK pRankStack1, const ppAT_RANK pRankStack2,
                      const AT_RANK *nAtomNumberCanonFrom,
                      sp_ATOM *at );

    int All_SC_Same( AT_RANK canon_rank1, /*  canonical number */
                      const ppAT_RANK pRankStack1, const ppAT_RANK pRankStack2,
                      const AT_RANK *nAtomNumberCanonFrom,
                      const sp_ATOM *at );


    int CompareLinCtStereoDoubleToValues( AT_STEREO_DBLE *LinearCTStereoDble,
                                  AT_RANK at_rank_canon1, AT_RANK at_rank_canon2, U_CHAR bond_parity );

    int CompareLinCtStereoAtomToValues( AT_STEREO_CARB *LinearCTStereoCarb,
                                  AT_RANK at_rank_canon1, U_CHAR parity );
    int CompareLinCtStereoDble( AT_STEREO_DBLE *LinearCTStereoDble1, int nLenLinearCTStereoDble1,
                                 AT_STEREO_DBLE *LinearCTStereoDble2, int nLenLinearCTStereoDble2 );
    int CompareLinCtStereoCarb( AT_STEREO_CARB *LinearCTStereoCarb1, int nLenLinearCTStereoCarb1,
                                 AT_STEREO_CARB *LinearCTStereoCarb2, int nLenLinearCTStereoCarb2 );
    int CompareLinCtStereo( AT_STEREO_DBLE *LinearCTStereoDble1, int nLenLinearCTStereoDble1,
                             AT_STEREO_CARB *LinearCTStereoCarb1, int nLenLinearCTStereoCarb1,
                             AT_STEREO_DBLE *LinearCTStereoDble2, int nLenLinearCTStereoDble2,
                             AT_STEREO_CARB *LinearCTStereoCarb2, int nLenLinearCTStereoCarb2 );
    /***************************************************************************/
    /* ichicans.c */
    int UnmarkNonStereo( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nRank, const AT_RANK *nAtomNumber, int bIsotopic );
    int FillSingleStereoDescriptors( CANON_GLOBALS *pCG, sp_ATOM *at, int i, int num_trans, const AT_RANK *nRank
                              , AT_STEREO_CARB *LinearCTStereoCarb, int *nStereoCarbLen, int nMaxStereoCarbLen
                              , AT_STEREO_DBLE *LinearCTStereoDble, int *nStereoDbleLen, int nMaxStereoDbleLen
                              , int bAllene );
    void SwitchAtomStereoAndIsotopicStereo( sp_ATOM *at, int num_atoms, int *bSwitched );
    void SetCtToIsotopicStereo( CANON_STAT *pCS, CANON_STAT *pCS2 );
    void SetCtToNonIsotopicStereo( CANON_STAT *pCS, CANON_STAT *pCS2 );
    int FillAllStereoDescriptors( CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon, CANON_STAT *pCS );
    int FillOutStereoParities( sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon,
                               const AT_RANK *nRank, const AT_RANK *nAtomNumber, CANON_STAT *pCS, CANON_GLOBALS *pCG, int bIsotopic );
    int InvertStereo( sp_ATOM *at, int num_at_tg,
                      AT_RANK *nCanonRank, AT_RANK *nAtomNumberCanon,
                      CANON_STAT *pCS, int bInvertLinearCTStereo );
    int find_atoms_with_parity( sp_ATOM *at, S_CHAR *visited, int from_atom, int cur_atom );
    int GetStereoNeighborPos( sp_ATOM *at, int iAt1, int iAt2 );
    int GetStereoBondParity( sp_ATOM *at, int i, int n, AT_RANK *nRank );
    int GetStereoCenterParity( CANON_GLOBALS *pCG, sp_ATOM *at, int i, AT_RANK *nRank );
    int GetPermutationParity( CANON_GLOBALS *pCG, sp_ATOM *at, AT_RANK nAvoidNeighbor, AT_RANK *nCanonRank );

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif    /* _INCHICOMN_H_ */
