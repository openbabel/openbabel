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


#ifndef __INCHICANO_H__
#define __INCHICANO_H__

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


int GetCanonLengths( int num_at,  sp_ATOM* at, ATOM_SIZES *s, T_GROUP_INFO *t_group_info );

int AllocateCS( CANON_STAT *pCS, int num_at, int num_at_tg, int nLenCT, int nLenCTAtOnly,
                     int nLenLinearCTStereoDble, int nLenLinearCTIsotopicStereoDble,
                     int nLenLinearCTStereoCarb, int nLenLinearCTIsotopicStereoCarb,
                     int nLenLinearCTTautomer, int nLenLinearCTIsotopicTautomer,
                     int nLenIsotopic, INCHI_MODE nMode, BCN *pBCN );

int DeAllocateCS( CANON_STAT *pCS );

void DeAllocBCN( BCN *pBCN );

int Canon_INChI(  int num_atoms, int num_at_tg, sp_ATOM* at, CANON_STAT* pCS, INCHI_MODE nMode, int bTautFtcn);
int GetBaseCanonRanking( int num_atoms, int num_at_tg, sp_ATOM* at[],
                         T_GROUP_INFO *t_group_info, ATOM_SIZES s[], BCN *pBCN, struct tagInchiTime *ulTimeOutTime );
int bCanonIsFinerThanEquitablePartition( int num_atoms, sp_ATOM* at, AT_RANK *nSymmRank );
int UpdateFullLinearCT( int num_atoms, int num_at_tg, sp_ATOM* at, AT_RANK *nRank, AT_RANK *nAtomNumber,
                        CANON_STAT* pCS, int bFirstTime );

int FixCanonEquivalenceInfo( int num_at_tg, AT_RANK *nSymmRank, AT_RANK *nCurrRank,
                             AT_RANK *nTempRank, AT_NUMB *nAtomNumber, int *bChanged);

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __INCHICANO_H__ */
