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

#include "bcf_s.h"

/****************************************************************************
  Check if all equivalent to cr1 possibly stereogenic atoms:
  1) have KNOWN parity, and
  2) their parities are same
****************************************************************************/
int All_SC_Same( AT_RANK canon_rank1, /*  canonical number */
                  const ppAT_RANK pRankStack1,
                  const ppAT_RANK pRankStack2,
                  const AT_RANK *nAtomNumberCanonFrom,
                  const sp_ATOM *at )
{
    int     n1 = (int) nAtomNumberCanonFrom[(int) canon_rank1 - 1];
    AT_RANK r1 = pRankStack1[0][n1];
    int     iMax1 = (int) r1;
    int     i1, s1;
    int     bFound = 0, stereo_atom_parity = -1;

    /*  find one stereo atom such that canon_rank1 can be mapped on it */
    for (i1 = 1; i1 <= iMax1 && r1 == pRankStack2[0][s1 = (int) pRankStack2[1][iMax1 - i1]]; i1++)
    {
        if (at[s1].stereo_bond_neighbor[0])
        {
            bFound = 0; /* at[s1] is not sp3-stereogenic: it belongs to a stereobond */
            break;
        }
        else
            if (i1 == 1)
            {
                stereo_atom_parity = PARITY_VAL( at[s1].stereo_atom_parity );
                if (!ATOM_PARITY_KNOWN( stereo_atom_parity ))
                {
                    bFound = 0;  /* at[s1] does not have a KNOWN parity */
                    break;
                }
            }
            else
                if (stereo_atom_parity != PARITY_VAL( at[s1].stereo_atom_parity ))
                {
                    bFound = 0; /* two equivalent atoms have different parities */
                    break;
                }
        bFound++;
    }

    return bFound;
}


/****************************************************************************
  get next available (not mapped yet) rank for a stereo center atom
****************************************************************************/
int Next_SC_At_CanonRank2( AT_RANK *canon_rank1,        /*  1st call input: largest canon number mapped so far or 0 */
                                                        /*  output: suggested canon. rank > than input if success */
                           AT_RANK *canon_rank1_min,    /*  1st call:0 next calls: first tried canon. number */
                           int *bFirstTime,             /*  1 at the time of the 1st call  */
                           S_CHAR *bAtomUsedForStereo,  /*  STEREO_AT_MARK if the atom has not been mapped yet */
                           const ppAT_RANK pRankStack1, /*  mapping ranks/sort order of atoms with canon. numbers (from) */
                           const ppAT_RANK pRankStack2, /*  mapping ranks/sort order of atoms with stereo (to) */
                           const AT_RANK *nAtomNumberCanonFrom, /*  sorted order of the canon. numbers */
                           int num_atoms )
{
    AT_RANK canon_rank1_inp = *canon_rank1;
    AT_RANK cr1;  /*  canonical rank (canonical number) */
    AT_RANK r1;   /*  mapping rank */
    int     n1;   /*  ord. number of an atom with the canon. number */
    int     s1;   /*  ord. number of an atom with stereo */
    int     i1, bFound = 0;
    int     iMax1;

    if (canon_rank1_inp < *canon_rank1_min)
    {
        canon_rank1_inp = *canon_rank1_min;
    }
    else
    {
        if (canon_rank1_inp < 1)
        {
            canon_rank1_inp = 1;
        }
        else
        {
            canon_rank1_inp++; /*  next canonical rank */
        }
    }
    cr1 = canon_rank1_inp;

    while ((int) cr1 <= num_atoms)
    {
        n1 = (int) nAtomNumberCanonFrom[(int) cr1 - 1]; /*  atom1 (which has canon. rank cr1) ord. number */
        iMax1 = (int) ( r1 = pRankStack1[0][n1] ); /*  mapping rank of atom1 */
        /*  find atoms "to" to which the canon. number can be mapped; they have mapping rank r1, number s1 */
        for (i1 = 1; i1 <= iMax1 && r1 == pRankStack2[0][s1 = (int) pRankStack2[1][iMax1 - i1]]; i1++)
        {
            /*  looking for a stereo center atom that has mapping rank r1 */
            if (bAtomUsedForStereo[s1] == STEREO_AT_MARK)
            {
                /*  found a sterogenic atom that has not been mapped yet */
                bFound = 1;
                break;
            }
        }
        if (bFound)
        {
            /*  one sterogenic not mapped yet atom "to" has been found */
            if (*bFirstTime)
            {
                *canon_rank1_min = cr1;
                *bFirstTime = 0;
            }
            break;
        }
        else
        {
             /*  a not mapped yet stereogenic atom has not found */
             /*  for the mapping rank r1 defined by the canonical rank cr1; try next cr1 */
            cr1++;
        }
    }
    if (bFound)
    {
        /*  success */
        *canon_rank1 = cr1;
        return 1;
    }

    return 0;
}


/****************************************************************************/
int CompareLinCtStereoDble( AT_STEREO_DBLE *LinearCTStereoDble1,
                            int nLenLinearCTStereoDble1,
                            AT_STEREO_DBLE *LinearCTStereoDble2,
                            int nLenLinearCTStereoDble2 )
{
    int i, num, ret = 0;

    /* compare double bonds */
    if (LinearCTStereoDble1 && LinearCTStereoDble2)
    {
        num = inchi_min( nLenLinearCTStereoDble1, nLenLinearCTStereoDble2 );
        for (i = 0; i < num; i++)
        {
            if ((ret = (int) LinearCTStereoDble1[i].at_num1 - (int) LinearCTStereoDble2[i].at_num1)) /* djb-rwth: addressing LLVM warning */
                break;
            if ((ret = (int) LinearCTStereoDble1[i].at_num2 - (int) LinearCTStereoDble2[i].at_num2)) /* djb-rwth: addressing LLVM warning */
                break;
            if ((ret = (int) LinearCTStereoDble1[i].parity - (int) LinearCTStereoDble2[i].parity)) /* djb-rwth: addressing LLVM warning */
                break;
        }
        if (!ret)
        {
            ret = nLenLinearCTStereoDble1 - nLenLinearCTStereoDble2;
        }
    }
    else
    {
        if (LinearCTStereoDble1 && nLenLinearCTStereoDble1 > 0)
        {
            ret = 1;
        }
        else
        {
            if (LinearCTStereoDble2 && nLenLinearCTStereoDble2 > 0)
            {
                ret = -1;
            }
        }
    }

    return ret;
}


/****************************************************************************/
int CompareLinCtStereoCarb( AT_STEREO_CARB *LinearCTStereoCarb1,
                            int nLenLinearCTStereoCarb1,
                            AT_STEREO_CARB *LinearCTStereoCarb2,
                            int nLenLinearCTStereoCarb2 )
{
    int i, num, ret = 0;

    /* compare stereocenters */
    if (LinearCTStereoCarb1 && LinearCTStereoCarb2)
    {
        num = inchi_min( nLenLinearCTStereoCarb1, nLenLinearCTStereoCarb2 );
        for (i = 0; i < num; i++)
        {
            if ((ret = (int) LinearCTStereoCarb1[i].at_num - (int) LinearCTStereoCarb2[i].at_num)) /* djb-rwth: addressing LLVM warning */
                break;
            if ((ret = (int) LinearCTStereoCarb1[i].parity - (int) LinearCTStereoCarb2[i].parity)) /* djb-rwth: addressing LLVM warning */
                break;
        }
        if (!ret)
        {
            ret = nLenLinearCTStereoCarb1 - nLenLinearCTStereoCarb2;
        }
    }
    else
        if (LinearCTStereoCarb1 && nLenLinearCTStereoCarb1 > 0)
        {
            ret = 1;
        }
        else
            if (LinearCTStereoCarb2 && nLenLinearCTStereoCarb2 > 0)
            {
                ret = -1;
            }

    return ret;
}


/****************************************************************************/
int CompareLinCtStereo( AT_STEREO_DBLE *LinearCTStereoDble1,
                        int nLenLinearCTStereoDble1,
                        AT_STEREO_CARB *LinearCTStereoCarb1,
                        int nLenLinearCTStereoCarb1,
                        AT_STEREO_DBLE *LinearCTStereoDble2,
                        int nLenLinearCTStereoDble2,
                        AT_STEREO_CARB *LinearCTStereoCarb2,
                        int nLenLinearCTStereoCarb2 )
{
    int ret;

    /* compare double bonds */
    ret = CompareLinCtStereoDble( LinearCTStereoDble1, nLenLinearCTStereoDble1,
                                   LinearCTStereoDble2, nLenLinearCTStereoDble2 );
    if (!ret)
    {
        ret = CompareLinCtStereoCarb( LinearCTStereoCarb1, nLenLinearCTStereoCarb1,

                                      LinearCTStereoCarb2, nLenLinearCTStereoCarb2 );
    }
    return ret;
}


/****************************************************************************/
int CompareLinCtStereoAtomToValues( AT_STEREO_CARB *LinearCTStereoCarb,
                                    AT_RANK at_rank_canon1,
                                    U_CHAR parity )
{
    if (LinearCTStereoCarb->at_num CT_GREATER_THAN at_rank_canon1)
    {
        return 1;
    }
    if (LinearCTStereoCarb->at_num != at_rank_canon1)
    {
        return -1;
    }
    if (LinearCTStereoCarb->parity CT_GREATER_THAN parity)
    {
        return 1;
    }
    if (LinearCTStereoCarb->parity != parity)
    {
        return -1;
    }

    return 0;
}


/****************************************************************************
  Find atom number from the mapping rank and return 1, or
  if the mapping rank is tied and the atom number is not unique then return 0
****************************************************************************/
int bUniqueAtNbrFromMappingRank( AT_RANK **pRankStack, AT_RANK nAtRank, AT_NUMB *nAtNumber )
{
    int       r = (int) nAtRank - 1;
    AT_NUMB   i = pRankStack[1][r];
    if (nAtRank == pRankStack[0][(int) i] &&
        ( !r || nAtRank != pRankStack[0][pRankStack[1][r - 1]] )
       )
    {
        *nAtNumber = i;
        return 1;
    }

    return 0;
}


/****************************************************************************
  Get minimal set (class) representative and partially compress the partitioning
  mcr = minimal class representative.
****************************************************************************/
AT_RANK nGetMcr( AT_RANK *nEqArray, AT_RANK n )
{
    AT_RANK n1, n2, mcr; /*  recursive version is much shorter. */

    n1 = nEqArray[(int) n];
    if (n == n1)
    {
        return n;
    }
    /*  1st pass: find mcr */
    while (n1 != ( n2 = nEqArray[(int) n1] ))
    {
        n1 = n2;
    }
    /*  2nd pass: copy mcr to each element of the set starting from nEqArray[n] */
    mcr = n1;
    n1 = n;
    while ( /*n1*/ mcr != ( n2 = nEqArray[(int) n1] ))
    {
        nEqArray[(int) n1] = mcr;
        n1 = n2;
    }

    return ( mcr );
}


/****************************************************************************
  Join 2 sets (classes) that have members n1 and n2
****************************************************************************/
int nJoin2Mcrs( AT_RANK *nEqArray, AT_RANK n1, AT_RANK n2 )
{
    n1 = nGetMcr( nEqArray, n1 );
    n2 = nGetMcr( nEqArray, n2 );
    if (n1 < n2)
    {
        nEqArray[n2] = n1;
        return 1; /*  a change has been made */
    }
    if (n2 < n1)
    {
        nEqArray[n1] = n2;
        return 1; /*  a change has been made */
    }

    return 0; /*  no changes */
}

/*********************************************************************************
 *  For all pairs of atoms that are:                                             *
 *  (a) connected by a possibly stereogenic bond                                 *
 *  (b) "equivalent" at this point to canon_rank1-canon_rank2 :                  *
 *  Check if they:                                                               *
 *  1) are connected by a stereo bond or cumulene bonds of the same length       *
 *  2) have KNOWN parity, and                                                    *
 *  3) their parities are same                                                   *
 *********************************************************************************/
int All_SB_Same( AT_RANK canon_rank1,
                 AT_RANK canon_rank2, /*  canonical numbers */
                 const ppAT_RANK pRankStack1,
                 const ppAT_RANK pRankStack2,
                 const AT_RANK *nAtomNumberCanonFrom,
                 sp_ATOM *at )
{
    int     n1 = (int) nAtomNumberCanonFrom[(int) canon_rank1 - 1]; /* at1 has canon_rank1 */
    int     n2 = (int) nAtomNumberCanonFrom[(int) canon_rank2 - 1]; /* at2 has canon_rank2 */
    AT_RANK r1 = pRankStack1[0][n1]; /* at1 mapping rank */
    AT_RANK r2 = pRankStack1[0][n2]; /* at2 mapping rank */
    AT_RANK rNeigh1, rNeigh2;
    int     iMax1 = (int) r1;
    /* int     iMax2 = (int)r2; */
    int     i1, i2, s1 = 0, s2 = 0, k1 = 0, k2, m, k, num_equal;
    int     bNotFound = 1, cumulene_len, stereo_bond_parity;

    /*  at the first atom that possibly may have canon_rank1 find one stereo bond such that */
    /*  canon_rank1-canon_rank2 possibly may be mapped on it */
    for (i1 = 1; i1 <= iMax1 && r1 == pRankStack2[0][s1 = (int) pRankStack2[1][iMax1 - i1]]; i1++)
    {
        /* at[n1] may be possible to map on at[s1] */
        for (k1 = 0, s2 = 0, bNotFound = 1;
              k1 < MAX_NUM_STEREO_BONDS && ( s2 = (int) at[s1].stereo_bond_neighbor[k1] ) &&
              ( bNotFound = ( r2 != pRankStack2[0][--s2] ) ); k1++)
            ; /* continue until the 1st at[s2] (to which at[n2] may be mapped) have been found */
        if (!bNotFound)
        {
            break; /* stop at 1st found */
        }
    }
    if (bNotFound)
    {
        return -1; /*  error: no mapping exists */
    }
    for (k2 = 0, m = 0; k2 < MAX_NUM_STEREO_BONDS && ( m = (int) at[s2].stereo_bond_neighbor[k2] ) && m - 1 != s1; k2++)
        ;
    if (m - 1 != s1)
    {
        return -1; /*  program error: stereo bond in opposite direction not found */
    }
    stereo_bond_parity = at[s1].stereo_bond_parity[k1];
    if (!PARITY_KNOWN( stereo_bond_parity ))
    {
        return 0;
    }
    cumulene_len = BOND_CHAIN_LEN( stereo_bond_parity );
    rNeigh1 = pRankStack2[0][(int) at[s1].neighbor[(int) at[s1].stereo_bond_ord[k1]]];
    rNeigh2 = pRankStack2[0][(int) at[s2].neighbor[(int) at[s2].stereo_bond_ord[k2]]];

    num_equal = 0;
    /*  Search among ALL neighbors because sometimes a stereo bond may be mapped on a non-stereo bond. */
    /*  If is so then return 0: not all mappings are stereo-equivalent */
    for (s1 = 1; s1 <= iMax1 && r1 == pRankStack2[0][i1 = (int) pRankStack2[1][iMax1 - s1]]; s1++)
    {
        for (k = 0; k < at[i1].valence; k++)
        {
            n1 = at[i1].neighbor[k];
            if (rNeigh1 != pRankStack2[0][n1])
            {
                continue; /*  wrong neighbor */
            }
            if (cumulene_len)
            {
                int prev, next, len, j;
                for (prev = i1, len = 0, next = n1; len < cumulene_len; len++)
                {
                    if (at[next].valence == 2 && !at[next].num_H)
                    {
                        j = ( (int) at[next].neighbor[0] == prev );
                        prev = next;
                        next = at[next].neighbor[j];
                    }
                    else
                    {
                        break; /*  cannot continue */
                    }
                }
                if (len != cumulene_len ||
                     r2 != pRankStack2[0][next] ||
                     rNeigh2 != pRankStack2[0][prev])
                {
                    /*  cumulene chain not found */
                    continue;
                }
                i2 = next;
            }
            else
            {
                i2 = n1;
            }
            /*  find if a stereogenic bond between at[i1]-at[i2] exists */
            for (k1 = 0; k1 < MAX_NUM_STEREO_BONDS &&
                ( m = (int) at[i1].stereo_bond_neighbor[k1] ) && m - 1 != i2; k1++)
                ;
            if (m - 1 != i2)
            {
                return 0;
            }
            for (k2 = 0; k2 < MAX_NUM_STEREO_BONDS &&
                ( m = (int) at[i2].stereo_bond_neighbor[k2] ) && m - 1 != i1; k2++)
                ;
            if (m - 1 != i1)
            {
                return 0;
            }
            if (at[i1].stereo_bond_parity[k1] != at[i2].stereo_bond_parity[k2])
            {
                return -1; /*  program error */
            }
            if (stereo_bond_parity != at[i1].stereo_bond_parity[k1])
            {
                return 0;
            }
            num_equal++;
        }
    }

    return num_equal;
}


/****************************************************************************
  get min. ranks for the stereo bond atoms
****************************************************************************/
int Next_SB_At_CanonRanks2( AT_RANK *canon_rank1,
                            AT_RANK *canon_rank2, /*  canonical numbers */
                            AT_RANK *canon_rank1_min,
                            AT_RANK *canon_rank2_min,
                            int *bFirstTime,
                            S_CHAR *bAtomUsedForStereo,
                            const ppAT_RANK pRankStack1,
                            const ppAT_RANK pRankStack2,
                            const AT_RANK *nCanonRankFrom,
                            const AT_RANK *nAtomNumberCanonFrom,
                            const sp_ATOM *at,
                            int num_atoms,
                            int bAllene )
{
    AT_RANK canon_rank1_inp = *canon_rank1;
    AT_RANK canon_rank2_inp = *canon_rank2;
    AT_RANK cr1, cr2; /*  canonical ranks (canonical numbers) */
    AT_RANK r1, r2;   /*  mapping ranks */
    int     n1, n2;   /*  ord. numbers of atoms with stereo */
    int     s1, s2;   /*  ord. numbers of atoms with canon. numbers */
    int     i1, i2, k, m;
    int     iMax1, iMax2;

    if (canon_rank1_inp < *canon_rank1_min ||
         (canon_rank1_inp == *canon_rank1_min &&
         canon_rank2_inp < *canon_rank2_min)) /* djb-rwth: addressing LLVM warning */
    {

        canon_rank1_inp = *canon_rank1_min;
        canon_rank2_inp = *canon_rank2_min;
    }
    else
        if (canon_rank1_inp < 2)
        {
            canon_rank1_inp = 2;
            canon_rank2_inp = 0;
        }
    cr1 = canon_rank1_inp;
    cr2 = num_atoms; /* initialize. 1/8/2002 */
    while ((int) cr1 <= num_atoms)
    {
        cr2 = cr1;
        n1 = (int) nAtomNumberCanonFrom[(int) cr1 - 1]; /*  atom1=at[n1] (which has canon. rank) ord. number */
        iMax1 = (int) ( r1 = pRankStack1[0][n1] ); /*  mapping rank of atom1 */
        for (i1 = 1; i1 <= iMax1 && r1 == pRankStack2[0][s1 = (int) pRankStack2[1][iMax1 - i1]]; i1++)
        {
            /*  looking for a stereo bond atom that has mapping rank r1 */
            /*  found at[s1] such that rank cr1 can be mapped on at[s1] because cr1 and s1 have equal */
            /*  mapping rank = r1. Check at[s1] stereo bonds */
            if (bAtomUsedForStereo[s1] && bAtomUsedForStereo[s1] < STEREO_AT_MARK)
            {
                for (k = 0; k < MAX_NUM_STEREO_BONDS && ( s2 = (int) at[s1].stereo_bond_neighbor[k] ); k++) /* djb-rwth: removing redundant code */
                {
                    /*  stereo bond at[s1]-at[s2] has been found */
                    if (bAtomUsedForStereo[--s2])
                    {
                        /*  stereo bonds have not been mapped. however, this check is not needed */
                        int cumulene_len = BOND_CHAIN_LEN( at[s1].stereo_bond_parity[k] );
                        if ((cumulene_len % 2 && !bAllene) || /* 09-26-2003 */
                             (!( cumulene_len % 2 ) && bAllene)) /* djb-rwth: addressing LLVM warning */
                        { /* 08-17-2003 Fix05 */
                            continue;
                        }
                        iMax2 = (int) ( r2 = pRankStack2[0][s2] ); /*  mapping rank of atom2 */
                        /*  Go back to canonical ranks and find an atom that has mapping rank r2 */
                        /*  and is connected to the atom with canonical rank cr1 (possibly by cumulene chain) */
                        /*  These cr1-cr2 canon. ranks possibly can be mapped on at[s1]-at[s2] stereo bond */
                        for (i2 = 1; i2 <= iMax2 && r2 == pRankStack1[0][n2 = (int) pRankStack1[1][iMax2 - i2]]; i2++)
                        {
                            if (cumulene_len)
                            {
                                int prev, next, len, j;
                                for (m = 0; m < at[n1].valence; m++)
                                {
                                    for (prev = n1, len = 0, next = (int) at[n1].neighbor[m]; len < cumulene_len; len++)
                                    {
                                        if (at[next].valence == 2 && !at[next].num_H)
                                        {
                                            j = ( (int) at[next].neighbor[0] == prev );
                                            prev = next;
                                            next = at[next].neighbor[j];
                                        }
                                        else
                                        {
                                            break; /*  cannot continue */
                                        }
                                    }
                                    if (len == cumulene_len && n2 == next)
                                    {
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                for (m = 0; m < at[n1].valence && n2 != (int) at[n1].neighbor[m]; m++)
                                    ;
                            }
                            if (m < at[n1].valence &&
                                 nCanonRankFrom[n2] < cr2 &&
                                 nCanonRankFrom[n2] > canon_rank2_inp)
                            {

                                cr2 = nCanonRankFrom[n2]; /*  found a candidate for cr2 */
                            }
                        }
                    }
                }
            }
        }
        if (cr2 >= cr1)
        {
            /*  not found for this r1 */
            cr1++;
            canon_rank2_inp = 0;
        }
        else
        {
             /* found cr2 < cr1 */
            if (*bFirstTime)
            {
                *canon_rank1_min = cr1;
                *canon_rank2_min = cr2;
                *bFirstTime = 0;
            }
            break;
        }
    }
    if (cr1 > cr2 && cr1 <= num_atoms)
    {
        /*  success */
        *canon_rank1 = cr1;
        *canon_rank2 = cr2;
        return 1;
    }

    return 0;
}


/****************************************************************************/
int NextStereoParity2Test( int *stereo_bond_parity,
                           int *sb_parity_calc,
                           int nNumBest,
                           int nNumWorse,
                           int nNumUnkn,
                           int nNumUndf,
                           int nNumCalc,
                           int vABParityUnknown )
{
    /* sequence of (stereo_bond_parity, sb_parity_calc) pairs:

          (BEST_PARITY, BEST_PARITY)  <calc>
                      |
          (BEST_PARITY, WORSE_PARITY) <known>
                      |
          (WORSE_PARITY, WORSE_PARITY) <calc>                (BEST_PARITY, 0) <known>
                       \___________________________________________/
                                              |
                                       (WORSE_PARITY, 0)   <known>
                                              |
                                       (AB_PARITY_UNKN, 0) <known>
                                              |
                                       (AB_PARITY_UNDF, 0) <known>
                                              |
                                       <next pair of ranks>
      Meaning:
      stereo_bond_parity is the parity we are looking for
      stereo_bond_parity==sb_parity_calc  => parity to be calculated from canonical numbers
      stereo_bond_parity!=sb_parity_calc  => parity is already known
     */
get_next_parity:
    switch (*stereo_bond_parity)
    {
        case BEST_PARITY:
            switch (*sb_parity_calc)
            {
                case 0:                                 /*  BEST_PARITY(known) : (BEST_PARITY, 0) -> */
                    *stereo_bond_parity = WORSE_PARITY;  /*  WORSE_PARITY(known): (WORSE_PARITY, 0) */
                    if (!nNumWorse)
                    {
                        goto get_next_parity;
                    }
                    break;
                case BEST_PARITY:                       /*  BEST_PARITY(calc) : (BEST_PARITY, BEST_PARITY) -> */
                    *sb_parity_calc = WORSE_PARITY;      /*  BEST_PARITY(known): (BEST_PARITY, WORSE_PARITY) */
                    if (!nNumBest)
                    {
                        goto get_next_parity;
                    }
                    break;
                case WORSE_PARITY:                      /*  BEST_PARITY(known): (BEST_PARITY, WORSE_PARITY)-> */
                    *stereo_bond_parity = WORSE_PARITY;  /*  WORSE_PARITY(calc): (WORSE_PARITY,WORSE_PARITY) */
                    if (!nNumCalc)
                    { /* added 12-17-2003 */
                        goto get_next_parity;
                    }
                    break;
            }
            break;
        case WORSE_PARITY:
            switch (*sb_parity_calc)
            {
                case 0:                                 /*  WORSE_PARITY(known)  : (WORSE_PARITY, 0) -> */
                    *stereo_bond_parity = vABParityUnknown /* AB_PARITY_UNKN */;/*  AB_PARITY_UNKN(known): (AB_PARITY_UNKN, 0) */
                    if (!nNumUnkn)
                    {
                        goto get_next_parity;
                    }
                    break;
                case BEST_PARITY:                       /*  error */
                    return CT_STEREOCOUNT_ERR;          /*   <BRKPT> */
                case WORSE_PARITY:                      /*  WORSE_PARITY(calc) : (WORSE_PARITY,WORSE_PARITY)-> */
                    *sb_parity_calc = 0;                 /*  WORSE_PARITY(known): (WORSE_PARITY, 0) */
                    if (!nNumWorse)
                    {
                        goto get_next_parity;
                    }
                    break;
            }
            break;

        case AB_PARITY_UNKN:                        /* AB_PARITY_UNKN(known): (AB_PARITY_UNKN, 0) -> */
            if (*sb_parity_calc)                 /*  error */
            {
                return CT_STEREOCOUNT_ERR;          /*   <BRKPT> */
            }
            *stereo_bond_parity = AB_PARITY_UNDF;    /* AB_PARITY_UNDF(known): (AB_PARITY_UNDF, 0) */
            if (!nNumUndf)
            {
                return 1; /*goto next_canon_ranks;*/
            }
            break;

        case AB_PARITY_UNDF:                        /*  AB_PARITY_UNDF(known): (AB_PARITY_UNDF, 0) -> */
            if (*sb_parity_calc)
            {                /*  error */
                return CT_STEREOCOUNT_ERR;          /*   <BRKPT> */
            }
            return 1; /*goto next_canon_ranks;*/     /*  next canon ranks */
    }
    return 0;
}


/****************************************************************************/
int CompareLinCtStereoDoubleToValues( AT_STEREO_DBLE *LinearCTStereoDble,
                                      AT_RANK at_rank_canon1,
                                      AT_RANK at_rank_canon2,
                                      U_CHAR bond_parity )
{
    if (LinearCTStereoDble->at_num1 CT_GREATER_THAN at_rank_canon1)
    {
        return 1;
    }
    if (LinearCTStereoDble->at_num1 != at_rank_canon1)
    {
        return -1;
    }
    if (LinearCTStereoDble->at_num2 CT_GREATER_THAN at_rank_canon2)
    {
        return 1;
    }
    if (LinearCTStereoDble->at_num2 != at_rank_canon2)
    {
        return -1;
    }
    if (LinearCTStereoDble->parity CT_GREATER_THAN bond_parity)
    {
        return 1;
    }
    if (LinearCTStereoDble->parity != bond_parity)
    {
        return -1;
    }

    return 0;
}


/****************************************************************************
  Set for at[i]:
   0                  if atom has no parity
   STEREO_AT_MARK=8   if atom has stereo parity and has no stereo bonds
   num_stereo_bonds   number of stereogenic bonds adjacent to the atom <= 3
****************************************************************************/
void SetUseAtomForStereo( S_CHAR *bAtomUsedForStereo, sp_ATOM *at, int num_atoms )
{
    int i, k;
    memset( bAtomUsedForStereo, 0, sizeof( bAtomUsedForStereo[0] )*num_atoms ); /* djb-rwth: memset_s C11/Annex K variant? */
    for (i = 0; i < num_atoms; i++)
    {
        if (at[i].parity)
        {
            for (k = 0; k < MAX_NUM_STEREO_BONDS && at[i].stereo_bond_neighbor[k]; k++)
            {
                ;
            }
            bAtomUsedForStereo[i] = k ? k : STEREO_AT_MARK;
        }
    }
}


 /****************************************************************************/
int CurTreeAlloc( CUR_TREE *cur_tree, int num_atoms )
{
    if (cur_tree)
    {
        if (cur_tree->tree && cur_tree->max_len > 0 && !( cur_tree->max_len % num_atoms ))
        {
            /*  do not reallocate */
            cur_tree->cur_len = 0;
            cur_tree->incr_len = num_atoms;
            memset( cur_tree->tree, 0, cur_tree->max_len * sizeof( cur_tree->tree[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            return 0; /*  ok */
        }
        inchi_free( cur_tree->tree );
        memset( cur_tree, 0, sizeof( *cur_tree ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        if ((cur_tree->tree = (AT_NUMB *) inchi_calloc( num_atoms, sizeof( cur_tree->tree[0] ) ))) /* djb-rwth: addressing LLVM warning */
        {
            cur_tree->incr_len =
                cur_tree->max_len = num_atoms;
            return 0; /*  ok */
        }
    }

    return -1; /*  error */ /*   <BRKPT> */
}


/****************************************************************************/
int CurTreeReAlloc( CUR_TREE *cur_tree )
{
    if (cur_tree)
    {
        if (cur_tree->tree && cur_tree->max_len > 0 && cur_tree->incr_len > 0)
        {
            void *p = cur_tree->tree;
            if ((cur_tree->tree = (AT_NUMB *) inchi_calloc( (long long)cur_tree->max_len + (long long)cur_tree->incr_len, sizeof( cur_tree->tree[0] ) ))) /* djb-rwth: cast operators added; addressing LLVM warning */
            {
                memcpy(cur_tree->tree, p, cur_tree->cur_len * sizeof(cur_tree->tree[0]));
                inchi_free( p );
                cur_tree->max_len += cur_tree->incr_len;
                return 0; /*  ok */
            }
        }
    }

    return -1; /*  error */ /*   <BRKPT> */
}


/****************************************************************************/void CurTreeFree( CUR_TREE *cur_tree )
{
    if (cur_tree)
    {
        inchi_free( cur_tree->tree );
        memset( cur_tree, 0, sizeof( *cur_tree ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    }
}


/****************************************************************************/int CurTreeAddRank( CUR_TREE *cur_tree, AT_NUMB rank )
{
    if (cur_tree)
    {
        if (cur_tree->cur_len + 2 > cur_tree->max_len)
        {
            if (CurTreeReAlloc( cur_tree ))
            {
                return -1; /*  error */ /*   <BRKPT> */
            }
        }
        cur_tree->tree[cur_tree->cur_len++] = rank;
        cur_tree->tree[cur_tree->cur_len++] = 1;
        return 0;
    }

    return -1;  /*  error  */ /*   <BRKPT> */
}


/****************************************************************************/
int CurTreeIsLastRank( CUR_TREE *cur_tree, AT_NUMB rank )
{
    if (cur_tree && cur_tree->cur_len > 0)
    {
        int rank_pos;
        rank_pos = cur_tree->cur_len - 1;
        rank_pos -= cur_tree->tree[rank_pos];
        if (rank_pos >= 0)
        {
            return ( rank == cur_tree->tree[rank_pos] );
        }
    }

    return 0;  /*  not found */
}


/****************************************************************************/
int CurTreeRemoveLastRankIfNoAtoms( CUR_TREE *cur_tree )
{
    if (cur_tree && cur_tree->tree && cur_tree->cur_len >= 2)
    {
        if (1 == cur_tree->tree[cur_tree->cur_len - 1])
        {
            return CurTreeRemoveLastRank( cur_tree ); /*  0=> success, -1=>failed */
        }
        return 1; /*  cannot remove */
    }
    return -1; /*  error */ /*   <BRKPT> */
}


/****************************************************************************/
int CurTreeAddAtom( CUR_TREE *cur_tree, int at_no )
{
    if (cur_tree)
    {
        if (cur_tree->cur_len + 1 > cur_tree->max_len)
        {
            if (CurTreeReAlloc( cur_tree ))
            {
                return -1; /*  error */ /*   <BRKPT> */
            }
        }
        if (cur_tree->cur_len > 0)
        {
            AT_NUMB new_len = cur_tree->tree[--cur_tree->cur_len] + 1;
            cur_tree->tree[cur_tree->cur_len++] = (AT_NUMB) at_no;
            cur_tree->tree[cur_tree->cur_len++] = new_len;
            return 0;
        }
    }

    return -1;
}


/****************************************************************************/
void CurTreeKeepLastAtomsOnly( CUR_TREE *cur_tree, int tpos, int shift )
{   /*  on first entry: shift = 1; other values may occur in subsequent recursion */
    /*  cur_tree[cur_tree->cur_len - shift] is the length of a segment */
    /*  action: remove all atoms except the last from all segments
                that have length value positon to the right from tpos */
    int cur_length_pos;
    if (cur_tree && cur_tree->tree && ( cur_length_pos = cur_tree->cur_len - shift ) > tpos)
    {
        if (cur_tree->tree[cur_length_pos] > 2)
        {
            /*  current segment contains more than 1 atom. Leave in the segment: rank, the last atom, length value */
            /*  subtract (old segment length)-(new segment length) from the tree length  */
            /*  actual segment length including segment length value = (cur_tree->tree[cur_length_pos]+1) */
            cur_tree->cur_len -= (int) cur_tree->tree[cur_length_pos] - 2;
            memmove(cur_tree->tree + cur_length_pos - cur_tree->tree[cur_length_pos] + 1, /*  1st atom pos */
                cur_tree->tree + cur_length_pos - 1,  /*  last atom in the current segment position */
                ((long long)shift + 1) * sizeof(cur_tree->tree[0])); /* djb-rwth: cast operator added */
            /*  (current segment length) distance from the last tree element has not changed */
            cur_tree->tree[cur_tree->cur_len - shift] = 2;
            /*  add 3 to move to the previous segment length position */
            shift += 3; /*  lenghth = 3 accounts for 3 currently present. segment items:
                            (1) the last atom, (2) rank, (3) length value */
        }
        else
        {
            shift += (int) cur_tree->tree[cur_length_pos] + 1; /*  cur_tree->cur_len - (previous segment length position) */
        }
        CurTreeKeepLastAtomsOnly( cur_tree, tpos, shift );
    }
}


/****************************************************************************/
int CurTreeRemoveIfLastAtom( CUR_TREE *cur_tree, int at_no )
{
    if (cur_tree && cur_tree->tree && cur_tree->cur_len > 2)
    {
        AT_NUMB len = cur_tree->tree[cur_tree->cur_len - 1];
        if (len >= 2 && (int) cur_tree->tree[cur_tree->cur_len - 2] == at_no)
        {
            cur_tree->tree[--cur_tree->cur_len - 1] = len - 1;
            return 0;
        }
        return 1; /*  not found */
    }

    return -1; /*  error */ /*   <BRKPT> */
}


/****************************************************************************/
int CurTreeGetPos( CUR_TREE *cur_tree )
{
    if (cur_tree)
    {
        return cur_tree->cur_len;
    }

    return -1;
}


/****************************************************************************/
int CurTreeSetPos( CUR_TREE *cur_tree, int len )
{
    if (cur_tree)
    {
        cur_tree->cur_len = len;
        return 0;
    }

    return -1;
}


/****************************************************************************/
int CurTreeRemoveLastRank( CUR_TREE *cur_tree )
{
    if (cur_tree && cur_tree->cur_len > 0)
    {
        cur_tree->cur_len -= cur_tree->tree[cur_tree->cur_len - 1] + 1;
        if (cur_tree->cur_len >= 0)
        {
            return 0;
        }
    }

    return -1;
}


/****************************************************************************
  Find if the atom is equivalent to already successfully tried current atoms
****************************************************************************/
int CurTreeIsLastAtomEqu( CUR_TREE *cur_tree, int at_no, AT_NUMB *nSymmStereo )
{
    if (cur_tree && cur_tree->tree && nSymmStereo && cur_tree->cur_len > 1)
    {
        AT_NUMB nEq = nSymmStereo[at_no];
        int end = cur_tree->cur_len - 1;
        int len = cur_tree->tree[end] - 1;
        for (; len > 0; len--)
        {
            if (nSymmStereo[(int) cur_tree->tree[end - len]] == nEq)
                return 1;
        }
        return 0;
    }

    return -1; /*  error */ /*   <BRKPT> */
}


#ifdef NEVER /* not used */
/****************************************************************************/
int CurTreeRemoveLastAtom( CUR_TREE *cur_tree )
{
    if (cur_tree && cur_tree->tree && cur_tree->cur_len > 2)
    {
        AT_NUMB len = cur_tree->tree[--cur_tree->cur_len];
        if (len >= 2)
        {
            cur_tree->tree[cur_tree->cur_len - 1] = len - 1;
            return 0;
        }
    }
    return -1;
}


/****************************************************************************/
int CurTreeReplaceLastRank( CUR_TREE *cur_tree, AT_NUMB rank )
{
    if (!CurTreeRemoveLastRank( cur_tree ))
    {
        return CurTreeAddRank( cur_tree, rank );
    }
    return -1;
}


/****************************************************************************
  returns cur_tree->cur_len for the block containing the rank
****************************************************************************/
int CurTreeFindTheRankPos( CUR_TREE *cur_tree, AT_NUMB rank )
{
    int i, k;
    if (cur_tree && cur_tree->tree && ( i = cur_tree->cur_len ) > 0)
    {
        while (0 <= ( k = i - (int) cur_tree->tree[i - 1] - 1 ))
        {
            if (cur_tree->tree[k] == rank)
            {
                return i;
            }
            i = k;
        }
    }

    return -1; /*  error */ /*   <BRKPT> */
}
#endif
