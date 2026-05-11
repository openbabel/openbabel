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


#ifndef _INCHICANO_H_
#define _INCHICANO_H_


#include "ichicant.h"

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


    int GetCanonLengths( int num_at,
                         sp_ATOM* at,
                         ATOM_SIZES *s,
                         T_GROUP_INFO *t_group_info );
    int AllocateCS( CANON_STAT *pCS,
                    int num_at,
                    int num_at_tg,
                    int nLenCT,
                    int nLenCTAtOnly,
                    int nLenLinearCTStereoDble,
                    int nLenLinearCTIsotopicStereoDble,
                    int nLenLinearCTStereoCarb,
                    int nLenLinearCTIsotopicStereoCarb,
                    int nLenLinearCTTautomer,
                    int nLenLinearCTIsotopicTautomer,
                    int nLenIsotopic,
                    INCHI_MODE nMode,
                    BCN *pBCN );
    int DeAllocateCS( CANON_STAT *pCS );
    void DeAllocBCN( BCN *pBCN );

    struct tagINCHI_CLOCK;

    int Canon_INChI( struct tagINCHI_CLOCK *ic,
                     int num_atoms,
                     int num_at_tg,
                     sp_ATOM* at,
                     CANON_STAT* pCS,
                     CANON_GLOBALS *pCG,
                     INCHI_MODE nMode,
                     int bTautFtcn );
    int GetBaseCanonRanking( struct tagINCHI_CLOCK *ic,
                             int num_atoms,
                             int num_at_tg,
                             sp_ATOM* at[],
                             T_GROUP_INFO *t_group_info,
                             ATOM_SIZES s[],
                             BCN *pBCN,
                             struct tagInchiTime *ulTimeOutTime,
                             CANON_GLOBALS *pCG,
                             int bFixIsoFixedH,
                             int LargeMolecules );
    int bCanonIsFinerThanEquitablePartition( int num_atoms,
                                             sp_ATOM* at,
                                             AT_RANK *nSymmRank );
    int UpdateFullLinearCT( int num_atoms,
                            int num_at_tg,
                            sp_ATOM* at,
                            AT_RANK *nRank,
                            AT_RANK *nAtomNumber,
                            CANON_STAT* pCS,
                            CANON_GLOBALS *pCG,
                            int bFirstTime );
    int FixCanonEquivalenceInfo( CANON_GLOBALS *pCG,
                                 int num_at_tg,
                                 AT_RANK *nSymmRank,
                                 AT_RANK *nCurrRank,
                                 AT_RANK *nTempRank,
                                 AT_NUMB *nAtomNumber,
                                 int *bChanged );
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* _INCHICANO_H_ */
