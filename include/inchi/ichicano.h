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


#ifndef __INCHICANO_H__
#define __INCHICANO_H__

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


int GetCanonLengths(int num_at,  sp_ATOM* at, ATOM_SIZES *s, T_GROUP_INFO *t_group_info );

int AllocateCS(CANON_STAT *pCS, int num_at, int num_at_tg, int nLenCT, int nLenCTAtOnly,
               int nLenLinearCTStereoDble, int nLenLinearCTIsotopicStereoDble,
               int nLenLinearCTStereoCarb, int nLenLinearCTIsotopicStereoCarb,
               int nLenLinearCTTautomer, int nLenLinearCTIsotopicTautomer,
               int nLenIsotopic, INCHI_MODE nMode, BCN *pBCN );

int DeAllocateCS(CANON_STAT *pCS );

void DeAllocBCN(BCN *pBCN );

int Canon_INChI(int num_atoms, int num_at_tg, sp_ATOM* at, 
                CANON_STAT* pCS, INCHI_MODE nMode, int bTautFtcn);
int GetBaseCanonRanking(int num_atoms, int num_at_tg, sp_ATOM* at[],
                        T_GROUP_INFO *t_group_info, ATOM_SIZES s[], BCN *pBCN, 
                        struct tagInchiTime *ulTimeOutTime, int bFixIsoFixedH );
int bCanonIsFinerThanEquitablePartition( int num_atoms, sp_ATOM* at, AT_RANK *nSymmRank );
int UpdateFullLinearCT(int num_atoms, int num_at_tg, sp_ATOM* at, AT_RANK *nRank, AT_RANK *nAtomNumber,
                       CANON_STAT* pCS, int bFirstTime );

int FixCanonEquivalenceInfo(int num_at_tg, AT_RANK *nSymmRank, AT_RANK *nCurrRank,
                            AT_RANK *nTempRank, AT_NUMB *nAtomNumber, int *bChanged);

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __INCHICANO_H__ */
