/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.03
 * May 9, 2010
 *
 * Originally developed at NIST
 * Modifications and additions by IUPAC and the InChI Trust
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-2.1.php
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




#define SB_DEPTH 6
/************************************************
   map_stereo_bonds4 and map_stereo_atoms4 use
   the following members of CANON_STAT *pCS:

        pCS->bKeepSymmRank  // ??? almost unused, replaced with nSymmStereo != NULL ???
        pCS->bFirstCT
        pCS->bStereoIsBetter
        pCS->lNumNeighListIter
        pCS->lNumBreakTies
        pCS->lNumRejectedCT
        pCS->lNumTotCT
        pCS->lNumEqualCT
        pCS->lNumDecreasedCT
        pCS->bExtract (bRELEASE_VERSION == 0)
        pCS->ulTimeOutTime

        pCS->bRankUsedForStereo
        pCS->bAtomUsedForStereo

        pCS->LinearCTStereoDble
        pCS->LinearCTStereoCarb
        pCS->nLenLinearCTStereoCarb
        pCS->nLenLinearCTStereoDble

        pCS->nPrevAtomNumber
 ************************************************/
/********************************************************************************/
int map_stereo_bonds4 ( 
                sp_ATOM *at, int num_atoms, int num_at_tg, int num_max, int bAllene,
                const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, /*  non-stereo canon ranking */
                AT_RANK *nCanonRankTo, /* output canonical stereo numbering*/
                const AT_RANK *nSymmRank, AT_RANK   **pRankStack1/*from*/,  AT_RANK **pRankStack2/*to*/,
                AT_RANK *nTempRank,       int         nNumMappedRanksInput,
                AT_RANK *nSymmStereo,     NEIGH_LIST *NeighList,
                CANON_STAT *pCS,          CUR_TREE   *cur_tree,  int nNumMappedBonds,
                int vABParityUnknown)
{
    int nTotSuccess = 0; /* 1=>full mapping has been completed;
                          * 2=>obtained a better stereo;
                          * 4=>restart (stereo bond or atom removed from the stereo CT)
                          */
    int tpos1;
    AT_STEREO_DBLE prevBond;
    tpos1 = CurTreeGetPos( cur_tree );

total_restart:

    if ( !nNumMappedBonds ) {
        
        memset( pCS->bRankUsedForStereo, 0, sizeof( pCS->bRankUsedForStereo[0] )*num_atoms );
        SetUseAtomForStereo( pCS->bAtomUsedForStereo, at, num_atoms );

        if ( pCS->bFirstCT && nSymmStereo && !pCS->bKeepSymmRank ) {
            int i;
            for ( i = 0; i < num_at_tg; i ++ ) {
                                    /*  nSymmStereo[i] = min. {k | at[k] stereo eq. to at[i]} */
                nSymmStereo[i] = i; /*  for union-join to keep track of stereo-equivalent atoms */
            }
        }
    }


    if ( nNumMappedBonds < pCS->nLenLinearCTStereoDble ) {

        int at_rank1, at_rank2, bStereoIsBetterWasSetHere;
        /* AT_RANK *nRankFrom=*pRankStack1++,  AT_RANK *nAtomNumberFrom=pRankStack1++; */
        /* AT_RANK *nRankTo  =*pRankStack2++,  AT_RANK *nAtomNumberTo  =pRankStack2++; */
        AT_RANK canon_min1, canon_min2;
        int bFirstCanonRank;
        int i, j, j1, j2, at_from1, at_from2, at_to1, at_to2, iMax, c;
        int nStackPtr[SB_DEPTH], nNumMappedRanks[SB_DEPTH], LastMappedTo1;
        int istk, istk2, istk3, bAddStack, nNumAtTo1Success;
        int ret1, ret2, parity1, parity2;

        AT_RANK at_rank_canon1; /*  = pCS->LinearCTStereoDble[nNumMappedBonds].at_num1; */ /*  canonical numbers of atoms */
        AT_RANK at_rank_canon2; /*  = pCS->LinearCTStereoDble[nNumMappedBonds].at_num2; */ /*  adjacent to the stereogenic bond */
        int nNumChoices, nNumUnkn, nNumUndf, nNumBest, nNumWorse, nNumCalc, sb_parity_calc;
        int stereo_bond_parity, prev_stereo_bond_parity, pass, bAllParitiesIdentical, bAllParitiesIdentical2;
        AT_STEREO_DBLE prevBond2;

        prevBond = pCS->LinearCTStereoDble[nNumMappedBonds];
        bFirstCanonRank=1;
        canon_min1=canon_min2=0;
/*        
        // find candidates for atom_from1, atom_to1; they must have identical mapping ranks
        at_rank1=pRankStack1[0][at_from1=nAtomNumberCanonFrom[(int)at_rank_canon1 - 1]]; // rank "from" for mapping
        at_rank2=pRankStack1[0][at_from2=nAtomNumberCanonFrom[(int)at_rank_canon2 - 1]]; // rank "from" for mapping
*/
        if ( nNumMappedBonds ) {
            at_rank_canon1 = pCS->LinearCTStereoDble[nNumMappedBonds-1].at_num1;
            at_rank_canon2 = pCS->LinearCTStereoDble[nNumMappedBonds-1].at_num2;
        } else {
            at_rank_canon1 = 0;
            at_rank_canon2 = 0;
        }
        goto bypass_next_canon_ranks_check;

next_canon_ranks:

        /*  Save time: avoid calling Next_SB_At_CanonRanks2() */
        if ( !pCS->bStereoIsBetter /* ??? && !pCS->bFirstCT ???*/ &&
              at_rank_canon1 >  pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 ||
              at_rank_canon1 == pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 &&
              at_rank_canon2 >= pCS->LinearCTStereoDble[nNumMappedBonds].at_num2  ) {

            if ( !nTotSuccess ) {
                pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond;
            }
            CurTreeSetPos( cur_tree, tpos1 );
            return nTotSuccess;
        }

bypass_next_canon_ranks_check:

        CurTreeSetPos( cur_tree, tpos1 );

        /*  find next available canon. numbers for a stereogenic bond pair of atoms */
        /*  process allenes AFTER all double bonds and odd-number-of-double-bonds cumulenes */
        if ( !(ret1 = Next_SB_At_CanonRanks2( &at_rank_canon1, &at_rank_canon2, /*  canonical numbers */
                                              &canon_min1, &canon_min2,
                                              &bFirstCanonRank, pCS->bAtomUsedForStereo,
                                              pRankStack1, pRankStack2,
                                              nCanonRankFrom, nAtomNumberCanonFrom,
                                              at, num_atoms, bAllene ) ) ) {
            /* failed to find next stereo bond to assign parity */
            if ( !bAllene && bFirstCanonRank ) {
                /* all stereobond have been processed; try to find allene to continue */
                AT_RANK at_rank_canon1_Allene = 0, canon_min1_Allene = 0;
                AT_RANK at_rank_canon2_Allene = 0, canon_min2_Allene = 0;
                if ( ret1 = Next_SB_At_CanonRanks2( &at_rank_canon1_Allene, &at_rank_canon2_Allene,
                                              &canon_min1_Allene, &canon_min2_Allene,
                                              &bFirstCanonRank, pCS->bAtomUsedForStereo,
                                              pRankStack1, pRankStack2,
                                              nCanonRankFrom, nAtomNumberCanonFrom,
                                              at, num_atoms, 1 ) ) {
                    at_rank_canon1 = at_rank_canon1_Allene;
                    at_rank_canon2 = at_rank_canon2_Allene;
                    canon_min1     = canon_min1_Allene;
                    canon_min2     = canon_min2_Allene;
                    bAllene        = 1; /* switch to allenes */
                } 
            }
        }
        
        if ( !ret1 || !pCS->bStereoIsBetter &&
             (at_rank_canon1 >  pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 ||
              at_rank_canon1 == pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 &&
              at_rank_canon2 >  pCS->LinearCTStereoDble[nNumMappedBonds].at_num2 ) ) {
            /* new ranks provide greater pCS->LinearCTStereoDble[nNumMappedBonds] and therefore rejected */
            if ( !nTotSuccess ) {
                pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond; /* restore stereo bond CT for the current bond */
            }
            return nTotSuccess;
        }
        /* current stereo bond initialization */
        nNumChoices = 0;
        nNumUnkn  = 0;
        nNumUndf  = 0;
        nNumBest  = 0;
        nNumWorse = 0;
        nNumCalc  = 0;
        pass=0;
        prev_stereo_bond_parity = 0;

        at_rank1=pRankStack1[0][at_from1=nAtomNumberCanonFrom[(int)at_rank_canon1 - 1]]; /* atom 1 rank "from" for mapping */
        at_rank2=pRankStack1[0][at_from2=nAtomNumberCanonFrom[(int)at_rank_canon2 - 1]]; /* atom 2 rank "from" for mapping */
        /* we are going to map bond (at[at_from1], at[at_from2]) and
           canonical ranks of its atoms (at_rank_canon1, at_rank_canon2)
           onto a stereogenic bond (at[at_to1], at[at_to2])
         */
        iMax = at_rank1-1;
        /*  test correctness: sorted pRankStack2[0][] and pRankStack1[0][] should have same ranks for both atoms */
        if ( at_rank1 != pRankStack2[0][pRankStack2[1][at_rank1-1]] ||
             at_rank2 != pRankStack2[0][pRankStack2[1][at_rank2-1]] ) {
            /* program error: "from" and "to" mapping ranks are not equal */
            return CT_STEREOCOUNT_ERR; /*   <BRKPT> */
        }
        /* -- do not check stereo features of "from" atoms:
           -- in case of "bond/charge isomerism" they may be missing.
        if ( !at[at_from1].stereo_bond_neighbor[0] ||
             !at[at_from2].stereo_bond_neighbor[0] )
            return CT_STEREOCOUNT_ERR; // program error
        */

        /*  find out if we have a choice in mapping: check all possible pairs (at_to1, at_to2)
            such that at_from1 is possibly constitutionally equivalent to at_to1, at_from2 to at_to2 */
        for ( j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1=pRankStack2[1][iMax-j1]]; j1 ++ ) {
            if ( !at[at_to1].stereo_bond_neighbor[0] )
                continue; /*  at_to1 does not belong to a stereo bond */
            for( j2 = 0; j2 < MAX_NUM_STEREO_BONDS &&
                         (at_to2  =at[at_to1].stereo_bond_neighbor[j2]); j2 ++ ) {
                at_to2 --;
                if ( pRankStack1[0][at_from2] != pRankStack2[0][at_to2] )
                    continue; /*  at_from2 cannot be mapped on at_to2 */
                stereo_bond_parity = PARITY_VAL(at[at_to1].stereo_bond_parity[j2]);
                i = 0;
                switch(stereo_bond_parity) {

                case AB_PARITY_UNDF: nNumUndf  ++; 
                                     break; /*  4 */
                case AB_PARITY_UNKN: nNumUnkn  ++; 
                                     break; /*  3 (occurs if forced different to UNDF)*/

                case BEST_PARITY:    nNumBest  ++; break; /*  1 */
                case WORSE_PARITY:   nNumWorse ++; break; /*  2 */
                case AB_PARITY_CALC: nNumCalc  ++; break; /*  6 */
                case AB_PARITY_NONE: i ++;         break; /*  0 */
                }
                nNumChoices += !i;
            }
        }
        if ( nNumChoices != nNumCalc + nNumUndf + nNumUnkn + nNumBest + nNumWorse ) {
            return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
        }
        if ( !nNumChoices ) {
            goto next_canon_ranks;
        }
        /*  Determine the first parity to search */
        sb_parity_calc = ( nNumCalc > 0 )? BEST_PARITY : 0;

        /*  ==============================================================
            Search sequence:           sb_parity_calc    stereo_bond_parity
            ==============================================================
            BEST_PARITY   (calc)       BEST_PARITY     BEST_PARITY
            BEST_PARITY   (known)      BEST_PARITY     WORSE_PARITY  or 0
            WORSE_PARITY  (calc)       WORSE_PARITY    WORSE_PARITY
            WORSE_PARITY  (known)      WORSE_PARITY    0
            AB_PARITY_UNKN(known)      AB_PARITY_UNKN  0
            AB_PARITY_UNDF(known)      AB_PARITY_UNDF  0

            if (sb_parity_calc==stereo_bond_parity) then "calc" else "known"
         */

repeat_all:

        if ( !nNumMappedBonds )
            pCS->bStereoIsBetter = 0;  /*  the first stereo feature in the canonical CT; moved here 7-13-2002 */

        if ( !pass ++ ) {
            /*  select the smallest (best) parity to search */
            if ( sb_parity_calc ) {
                stereo_bond_parity = BEST_PARITY;
            } else {
                stereo_bond_parity = nNumBest?   BEST_PARITY :
                                     nNumWorse?  WORSE_PARITY :
                                     nNumUnkn?   AB_PARITY_UNKN :
                                     nNumUndf?   AB_PARITY_UNDF : AB_PARITY_NONE;
            }
        } else {
            /* second pass: since the first pass failed, search for a worse result */
            prev_stereo_bond_parity = stereo_bond_parity;
            i = NextStereoParity2Test( &stereo_bond_parity, &sb_parity_calc,
                                     nNumBest, nNumWorse, nNumUnkn, nNumUndf, nNumCalc, vABParityUnknown);
            switch ( i ) {
            case 0:
                break; /* obtained next parity to test */
            case 1:
                goto next_canon_ranks;
            default:
                return i; /* program error */
            }
        }
        if ( stereo_bond_parity == AB_PARITY_NONE ) {
            /*  error? */
            return CT_STEREOCOUNT_ERR;                   /*   <BRKPT> */
        }
        /*  check if the new requested parity is good (small) enough */
        if ( !pCS->bStereoIsBetter ) {
            c = CompareLinCtStereoDoubleToValues( nTotSuccess? pCS->LinearCTStereoDble+nNumMappedBonds : &prevBond,
                              at_rank_canon1, at_rank_canon2, (U_CHAR)stereo_bond_parity );
            if ( c < 0 ) {
                if ( !nTotSuccess ) {
                    pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond;
                }
                CurTreeSetPos( cur_tree, tpos1 );
                return nTotSuccess;
            }
        }

        bAllParitiesIdentical     =  0;
        bAllParitiesIdentical2    =  0;
        LastMappedTo1             = -1;
        bStereoIsBetterWasSetHere =  0;
        istk = istk2 = istk3      =  0;

        if ( !nNumMappedBonds && prev_stereo_bond_parity != stereo_bond_parity )
            pCS->bStereoIsBetter = 0;  /*  the first stereo feature in the canonical CT; moved here 5-24-2002 */

        if ( prev_stereo_bond_parity != stereo_bond_parity ) {
            CurTreeSetPos( cur_tree, tpos1 );  /*  start over */
        }

        /* Mapping: here at_rank1 = nRankTo, at_to1 = nAtomNumberTo */
        for ( j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1=pRankStack2[1][iMax-j1]]; j1 ++ ) {
            nNumAtTo1Success = 0;
            if ( !at[at_to1].stereo_bond_neighbor[0] )
                continue; /*  at_to1 does not belong to a stereo bond */
            if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 )  &&
                 1 == CurTreeIsLastAtomEqu( cur_tree, at_to1, nSymmStereo ) ) {
                /* at_to1 is known to be stereogenically equivalent to another atom tried with at_rank_canon1 */
                continue;
            }
            bAllParitiesIdentical2 = 0;
            for( j2 = 0; j2 < MAX_NUM_STEREO_BONDS && (at_to2  =at[at_to1].stereo_bond_neighbor[j2]); j2 ++ ) {
                EQ_NEIGH  EN1[2], EN2[2];
                int bond_parity, num1, num2;
                AT_RANK at_rank_canon_n1, at_rank_canon_n2;

                at_to2 --;
                if ( pRankStack1[0][at_from2] != pRankStack2[0][at_to2] )
                    continue; /*  at_from2 cannot be mapped on at_to2 even without mapping at_from1 to at_to1 */
                
                /*  check whether the bond parity corresponds to the requested bond parity */
                if ( PARITY_KNOWN(at[at_to1].stereo_bond_parity[j2]) ) {
                    if ( stereo_bond_parity == sb_parity_calc ) {
                        continue;  /*  requested parity to be calculated, found known parity */
                    }
                    if ( stereo_bond_parity != PARITY_VAL(at[at_to1].stereo_bond_parity[j2]) ) {
                        continue;  /*  parity differs from the requested parity */
                    }
                } else
                if ( PARITY_CALCULATE( at[at_to1].stereo_bond_parity[j2]) ) {
                    if ( stereo_bond_parity != sb_parity_calc ) {
                        continue;  /*  requested known parity, found parity to be calculated */
                    }
                } else {
                    return CT_STEREOCOUNT_ERR;  /*  unknown parity type */ /*   <BRKPT> */
                }
                /*  initialize stack pointer nStackPtr[istk] for "hand-made" recursion */
                /*  stacks are pRankStack1[], pRankStack2[], nNumMappedRanks[] */
                istk                = 0;
                nStackPtr[0]        = 0;
                nNumMappedRanks[0]  = nNumMappedRanksInput;
                bAddStack           = 0;
                bAllParitiesIdentical  = ((at[at_to1].stereo_bond_parity[j2] & KNOWN_PARITIES_EQL )) &&
                                          PARITY_KNOWN(at[at_to1].stereo_bond_parity[j2]);
                
                if ( !bAllParitiesIdentical && !nNumCalc &&
                     (!nNumUndf + !nNumUnkn + !nNumBest + !nNumWorse )==3) {
                    /* only one kind of bond parity is present; check whether all parities are really same */
                    bAllParitiesIdentical = All_SB_Same(  at_rank_canon1, at_rank_canon2, /*  canonical numbers */
                                                          pRankStack1, pRankStack2,
                                                          nAtomNumberCanonFrom, at );
                    if ( bAllParitiesIdentical < 0 ) {
                        return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
                    }
                }
                
                /*****************************************************************
                 * do the mapping only if parities are not same
                 */
                if ( !bAllParitiesIdentical ) {
                    /*  map atom 1 or reuse previous mapping */
                    if ( LastMappedTo1 != at_to1 ) {
                        /*  avoid repetitve mapping to the same first at_to1 using LastMappedTo1 variable */
                        /*  map atom 1 */
                        ret1 = map_an_atom2( num_at_tg, num_max, at_from1, at_to1,
                                            nTempRank, nNumMappedRanks[istk], &nNumMappedRanks[istk+1], pCS,
                                            NeighList, pRankStack1+nStackPtr[istk], pRankStack2+nStackPtr[istk],
                                            &bAddStack );
                        if ( RETURNED_ERROR(ret1) ) {
                            return ret1; /*  error */
                        }
                        nStackPtr[istk+1] = nStackPtr[istk]+bAddStack;
                        LastMappedTo1 = at_to1;
                        if ( bAddStack ) {
                            if ( tpos1 == CurTreeGetPos( cur_tree ) ||
                                 0 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
                                CurTreeAddRank( cur_tree, at_rank_canon1 );
                            }
                            CurTreeAddAtom( cur_tree, at_to1 );
                        }
                    }
                    istk ++; /*  = 1 */
                    /*  check if we can map atom 2 */
                    if ( pRankStack1[nStackPtr[istk]][at_from2] != pRankStack2[nStackPtr[istk]][at_to2] ) {
                        /*
                         * This may happen when:
                         * A) Charge/bond isomerism, for example cyclopentadiene(-), or
                         * B) possibly stereogenic bond in an alternating ring has heighbors
                         * in 2 symmetrically attached rings.
                         * Such an alternating bond cannot be mapped on possibly stereogenic bond
                         * that has neighbors belonging to 1 of the symmetrically attached rings only.
                         * For example:
                         *   A---B---C---D  If all atoms are Carbons then B, C, F, G are constitutionally
                         *  ||  ||  ||  ||  equivalent. However, bonds B-C, F-G are not equivalent to
                         *  ||  ||  ||  ||  B-F and C-G and cannot be mapped on them.
                         *   E---F---G---H  If at_from1=B, at_from2=F, at_to1=B, then at_from2 cannot be mapped on at_to2=C
                         *                  If at_from1=B, at_from2=F, at_to1=C, then at_from2 cannot be mapped on at_to2=B
                         *                  etc.
                         */
                        if ( sb_parity_calc != stereo_bond_parity) {
                            /* can be passed only once for each bond */
                            nNumChoices --;
                            nNumUndf -= (stereo_bond_parity == AB_PARITY_UNDF);
                            nNumUnkn -= (stereo_bond_parity == AB_PARITY_UNKN);
                            nNumBest -= (stereo_bond_parity == BEST_PARITY);
                            nNumWorse-= (stereo_bond_parity == WORSE_PARITY);
                            /* nNumCalc  = nNumChoices - (nNumUndf + nNumUnkn + nNumBest + nNumWorse); */
                        } else
                        if ( sb_parity_calc == BEST_PARITY ) {
                            /* can be passed 2 times: for BEST_PARITY and WORSE_PARITY in this order */
                            nNumChoices --; /*  do not repeate for WORSE_PARITY */
                            nNumCalc    --;
                        }
                        continue;  /*  Happens for ID=80036,80253,91354,95532,101532,103788 */
                    }
                    if ( nStackPtr[istk] > nStackPtr[istk-1] ) { 
                        bAllParitiesIdentical2 = All_SB_Same(  at_rank_canon1, at_rank_canon2,
                                                  pRankStack1+nStackPtr[istk], pRankStack2+nStackPtr[istk],
                                                  nAtomNumberCanonFrom, at );
                        if ( bAllParitiesIdentical2 < 0 ) {
                            return CT_STEREOBOND_ERROR;  /*   <BRKPT> */
                        }
                    } else {
                        bAllParitiesIdentical2 = 0;
                    }
                    if ( bAllParitiesIdentical2 ) {
                        /*  do no mapping when all equivalent bonds have same parity */
                        /*  stereo_bond_parity = PARITY_VAL(at[at_to1].stereo_bond_parity[j2]); */
                        ClearPreviousMappings( pRankStack1+nStackPtr[istk]+2 );
                    } else {
                        if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                             1 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ) &&
                             1 == CurTreeIsLastAtomEqu( cur_tree, at_to2, nSymmStereo ) ) {
                                continue;
                        }
                        /*  map atom 2 */
                        ret2 = map_an_atom2( num_at_tg, num_max, at_from2, at_to2,
                                            nTempRank, nNumMappedRanks[istk], &nNumMappedRanks[istk+1], pCS,
                                            NeighList, pRankStack1+nStackPtr[istk], pRankStack2+nStackPtr[istk],
                                            &bAddStack );
                        if ( RETURNED_ERROR(ret2) ) {
                            return ret2; /*  program error */
                        }
                        nStackPtr[istk+1] = nStackPtr[istk]+bAddStack;
                        istk ++; /*  = 2 */
                        if ( bAddStack ) {
                            if ( tpos1 == CurTreeGetPos( cur_tree ) ||
                                 0 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ) ) {
                                CurTreeAddRank( cur_tree, at_rank_canon2 );
                            }
                            CurTreeAddAtom( cur_tree, at_to2 );
                        }
                    }
                } else {
                    /*  do no mapping when all equivalent bonds have same parity */
                    /*  stereo_bond_parity = PARITY_VAL(at[at_to1].stereo_bond_parity[j2]); */
                    ClearPreviousMappings( pRankStack1+2 );
                }

                /*  we have a precalculated (known) bond parity */
                
                
                /************************************************************
                 *
                 *   Known Bond Parity case: do not map stereo bond neighbors
                 */
                if ( stereo_bond_parity != sb_parity_calc ) /*  parity is known */
                {
                    /*  accept bond parity and do not map the neighbors */
                    bond_parity = stereo_bond_parity;
                    /*  same code as under " make a decision to accept current mapping" comment below */
                    /*  with one exception: istk instead of istk3 */
                    c = CompareLinCtStereoDoubleToValues( pCS->LinearCTStereoDble+nNumMappedBonds,
                                              at_rank_canon1, at_rank_canon2, (U_CHAR)bond_parity );
                    if ( c < 0 && !pCS->bStereoIsBetter ) {
                        
                        /*  reject */
                        
                        pCS->lNumRejectedCT ++;
                        /*  remove failed atom2 from the tree */
                        if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                             1 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ) ) {
                            CurTreeRemoveIfLastAtom( cur_tree, at_to2 );
                            CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                        }
                        continue;  /*  to next at_to2; Reject this at_to2: not a minimal CT. */
                    }  else {

                        /*  accept */

                        if ( c > 0 && !pCS->bStereoIsBetter ) {
                            /*  bond entry is less than the previusly found */
                            pCS->bStereoIsBetter = bStereoIsBetterWasSetHere = 1;
                            prevBond2 = pCS->LinearCTStereoDble[nNumMappedBonds];
                        }
                        pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 = at_rank_canon1;
                        pCS->LinearCTStereoDble[nNumMappedBonds].at_num2 = at_rank_canon2;
                        pCS->LinearCTStereoDble[nNumMappedBonds].parity  = bond_parity;
                        /*  recursive call */
                        pCS->bRankUsedForStereo[at_from1] ++;
                        pCS->bRankUsedForStereo[at_from2] ++;
                        if ( !bAllParitiesIdentical ) {
                            pCS->bAtomUsedForStereo[at_to1] --;
                            pCS->bAtomUsedForStereo[at_to2] --;
                        }
                        ret2 = map_stereo_bonds4 ( at, num_atoms, num_at_tg, num_max, bAllene, nCanonRankFrom, nAtomNumberCanonFrom, nCanonRankTo,
                                                  nSymmRank, pRankStack1+nStackPtr[istk],  pRankStack2+nStackPtr[istk],
                                                  nTempRank, nNumMappedRanks[istk], nSymmStereo, NeighList,
                                                  pCS, cur_tree, nNumMappedBonds+1 ,
                                                  vABParityUnknown);
                        if ( !bAllParitiesIdentical ) {
                            pCS->bAtomUsedForStereo[at_to1] ++;
                            pCS->bAtomUsedForStereo[at_to2] ++;
                        }
                        pCS->bRankUsedForStereo[at_from1] --;
                        pCS->bRankUsedForStereo[at_from2] --;
                        if ( ret2 == 4 ) {
                            if ( nNumMappedBonds ) {
                                return ret2;
                            } else {
                                pCS->bFirstCT = 1;
                                goto total_restart;
                            }
                        }

                        if ( RETURNED_ERROR(ret2) ) {
                            if ( ret2 == CT_TIMEOUT_ERR )
                                return ret2;
                            else
                                return ret2; /*  program error */
                        }
                        if ( ret2 > 0 ) {
                            nTotSuccess |= 1;
                            nNumAtTo1Success ++;
                            if ( bStereoIsBetterWasSetHere || (ret2 & 2) ) {
                                CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                                nTotSuccess |= 2; /*  Obtained a smaller CT */
                            }
                        } else {
                            if ( bStereoIsBetterWasSetHere ) { /*  rollback */
                                pCS->bStereoIsBetter = 0;
                                pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond2;
                            }
                            /*  remove failed atom2 from the tree */
                            if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ) ) {
                                CurTreeRemoveIfLastAtom( cur_tree, at_to2 );
                                CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                            }
                            /*
                            if ( 1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
                                CurTreeRemoveLastAtom( cur_tree, at_to1 );
                                CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                            }
                            */
                        }
                        bStereoIsBetterWasSetHere = 0;
                    }
                    if ( bAllParitiesIdentical || bAllParitiesIdentical2 ) {
                        break; /* j2 cycle, at_to2 (no need to repeat) */
                    }
                    continue; /*  to next at_to2 */
                }
                /***************************************************************************
                 *
                 *   Unknown Bond Parity case: may need to map stereo bond neighbors
                 *
                 ****************************************************************************
                 * Ranks are not known in advance
                 * check if at_from1/at_to1 half-bond has neighbors with equal mapping ranks
                 */

                parity1 = parity_of_mapped_half_bond( at_from1, at_to1, at_from2, at_to2, at, &EN1[0],
                           nCanonRankFrom, pRankStack1[nStackPtr[istk]], pRankStack2[nStackPtr[istk]] );
                /* old approach -- before E/Z parities 
                parity1 = parity_of_mapped_atom2( at_from1, at_to1, at, &EN1[0],
                               nCanonRankFrom, pRankStack1[nStackPtr[istk]], pRankStack2[nStackPtr[istk]] );
                 */
                /*  the following commented out statement is not needed here. */
                /*  parity2 = parity_of_mapped_atom2( at_from2, at_to2, at, &EN2[0], 
                                                     nCanonRankFrom, pRankStack1[nStackPtr[istk]],
                                                     pRankStack2[nStackPtr[istk]] );
                 */
                if ( !parity1 ) {
                    return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                }
                num1 = parity1 > 0? 1:2; /*  parity < 0 means additional mapping is needed to set parity */
                /*  --- try all possible mappings of the stereo bond ending atoms' neighbors --- */
                at_rank_canon_n1 = 0;
                at_rank_canon_n2 = 0;
                for ( i = 0; i < num1; i ++ ) {
                    int at_from_n1, at_to_n1, at_no_n1_num_success = 0;
                    istk2 = istk;
                    if ( num1 == 2 ) {
                        at_rank_canon_n1 = nCanonRankFrom[EN1[0].from_at];
                        /*  an additional neighbor mapping is necessary; */
                        /*  we need to map only one at_from1 neighbor to make all neighbors have different ranks */

                        at_from_n1 = EN1[0].from_at;
                        at_to_n1   = EN1[0].to_at[i];

                        if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                             1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n1 ) &&
                             1 == CurTreeIsLastAtomEqu( cur_tree, at_to_n1, nSymmStereo ) )
                            continue;
                        /*
                        if ( nSymmStereo && !pCS->bFirstCT ) {
                            if ( i && nSymmStereo[at_to_n1] == nSymmStereo[(int)EN1[0].to_at[0]] ) {
                                continue; // do not test stereo equivalent atoms except the first one
                            }
                        }
                        */
                        /*  neighbors are tied. Untie them by breaking a tie on ONE of them. */
                        ret1 = map_an_atom2( num_at_tg, num_max, at_from_n1, at_to_n1,
                                            nTempRank, nNumMappedRanks[istk2], &nNumMappedRanks[istk2+1], pCS,
                                            NeighList, pRankStack1+nStackPtr[istk2], pRankStack2+nStackPtr[istk2],
                                            &bAddStack );
                        if ( RETURNED_ERROR(ret1) ) {
                            return ret1; /*  program error */ /*   <BRKPT> */
                        }
                        nStackPtr[istk2+1] = nStackPtr[istk2] + bAddStack;
                        istk2 ++;  /*  <= 3 */
                        /*  debug */
                        if ( istk2 >= SB_DEPTH ) {
                            return CT_OVERFLOW; /*  program error */ /*   <BRKPT> */
                        }
                        if ( bAddStack ) {
                            if ( tpos1 == CurTreeGetPos( cur_tree ) ||
                                 0 == CurTreeIsLastRank( cur_tree, at_rank_canon_n1 ) ) {
                                CurTreeAddRank( cur_tree, at_rank_canon_n1 );
                            }
                            CurTreeAddAtom( cur_tree, at_to_n1 );
                        }


                        /*  now that all at_from1 neighbors have been mapped the parity must be defined */
                        parity1 = parity_of_mapped_half_bond( at_from1, at_to1, at_from2, at_to2, at, &EN1[1],
                           nCanonRankFrom, pRankStack1[nStackPtr[istk2]], pRankStack2[nStackPtr[istk2]] );
                        if ( parity1 <= 0 )
                            return CT_STEREOCOUNT_ERR;  /*  program error */ /*   <BRKPT> */
                    } else {
                        nNumMappedRanks[istk2+1] = nNumMappedRanks[istk2];
                        nStackPtr[istk2+1] = nStackPtr[istk2];
                        istk2 ++;  /*  <= 3 */
                    }

                    /*  check if at_from2/at_to2 half-bond has neighbors with equal mapping ranks */
                    parity2 = parity_of_mapped_half_bond( at_from2, at_to2, at_from1, at_to1, at, &EN2[0],
                           nCanonRankFrom, pRankStack1[nStackPtr[istk2]], pRankStack2[nStackPtr[istk2]] );
                    if ( !parity2 ) {
                        return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                    }
                    num2 = parity2 > 0? 1:2;
                    at_rank_canon_n2 = 0;
                    for ( j = 0; j < num2; j ++ ) {
                        int at_from_n2, at_to_n2;
                        istk3 = istk2;
                        if ( num2 == 2 ) {
                            at_rank_canon_n2 = nCanonRankFrom[EN2[0].from_at];
                            /*  we need to map only one at_from2 neighbor to make its neighbors have different ranks */
                            at_from_n2 = EN2[0].from_at;
                            at_to_n2   = EN2[0].to_at[j];

                            if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ) &&
                                 1 == CurTreeIsLastAtomEqu( cur_tree, at_to_n2, nSymmStereo ) )
                                continue;
                            
                            /*
                            if ( nSymmStereo && !pCS->bFirstCT ) {
                                if ( j && nSymmStereo[at_to_n2] == nSymmStereo[(int)EN2[0].to_at[0]] ) {
                                    continue; // do not test stereo equivalent atoms except the first one
                                }
                            }
                            */
                            /*  neighbors are tied. Untie them by breaking a tie on ONE of them. */
                            ret1 = map_an_atom2( num_at_tg, num_max, at_from_n2, at_to_n2,
                                                nTempRank, nNumMappedRanks[istk3], &nNumMappedRanks[istk3+1], pCS,
                                                NeighList, pRankStack1+nStackPtr[istk3],
                                                pRankStack2+nStackPtr[istk3],
                                                &bAddStack );
                            if ( RETURNED_ERROR(ret1) ) {
                                return ret1; /*  program error */
                            }
                            nStackPtr[istk3+1] = nStackPtr[istk3]+bAddStack;
                            istk3 ++;  /*  <= 4 */

                            if ( bAddStack ) {
                                if ( tpos1 == CurTreeGetPos( cur_tree ) ||
                                     0 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ) ) {
                                    CurTreeAddRank( cur_tree, at_rank_canon_n2 );
                                }
                                CurTreeAddAtom( cur_tree, at_to_n2 );
                            }

                            parity2 = parity_of_mapped_half_bond( at_from2, at_to2, at_from1, at_to1, at, &EN2[1],
                                     nCanonRankFrom, pRankStack1[nStackPtr[istk3]], pRankStack2[nStackPtr[istk3]] );
                            if ( parity2 <= 0 ) {
                                return CT_STEREOCOUNT_ERR;  /*  program error */ /*   <BRKPT> */
                            }
                        } else {
                            /*  no additional mapping is needed to set atom's parity */
                            nNumMappedRanks[istk3+1] = nNumMappedRanks[istk3];
                            nStackPtr[istk3+1] = nStackPtr[istk3];
                            istk3 ++;  /*  <= 4 */
                        }
                        
                        /*******************************************************************
                         * at this point the stereo bond is fully mapped to find its parity
                         *******************************************************************/

                        if ( parity1 <= 0 || parity2 <= 0 ) {
                            return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                        }
                        
                        /*  find current bond parity  AB_PARITY_ODD */
                        if ( ATOM_PARITY_WELL_DEF(parity1) && ATOM_PARITY_WELL_DEF(parity2) ) {
                            bond_parity = 2 - (parity1 + parity2)%2;
                        } else {
                            bond_parity = inchi_max(parity1, parity2);
                        }
                        if ( ATOM_PARITY_WELL_DEF(bond_parity) && at[at_to1].stereo_bond_z_prod[j2] < 0 )
                            bond_parity = 2 - (bond_parity+1)%2; /*  invert the bond parity */


                        /********************************************************
                         * make a decision whether to accept the current mapping
                         */
                        c = CompareLinCtStereoDoubleToValues( pCS->LinearCTStereoDble+nNumMappedBonds,
                                                  at_rank_canon1, at_rank_canon2, (U_CHAR)bond_parity );
                        if ( sb_parity_calc != bond_parity ||
                             c < 0 && !pCS->bStereoIsBetter ) {
                            /*  reject */
                            pCS->lNumRejectedCT ++;
                            /*  remove failed atom2 from the tree */
                            if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ) ) {
                                CurTreeRemoveIfLastAtom( cur_tree, at_to_n2 );
                                CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                            }
                            continue;  /*  Reject: not a minimal CT. */

                        }  else {

                            /*  try to accept */

                            if ( c > 0 && !pCS->bStereoIsBetter ) {
                                /*  bond_parity is less than the previusly found */
                                pCS->bStereoIsBetter = bStereoIsBetterWasSetHere = 1;
                                prevBond2 = pCS->LinearCTStereoDble[nNumMappedBonds];
                            }
                            /*  accept */
                            pCS->LinearCTStereoDble[nNumMappedBonds].at_num1 = at_rank_canon1;
                            pCS->LinearCTStereoDble[nNumMappedBonds].at_num2 = at_rank_canon2;
                            pCS->LinearCTStereoDble[nNumMappedBonds].parity  = bond_parity;
                            /*  recursive call */
                            pCS->bRankUsedForStereo[at_from1] ++;
                            pCS->bRankUsedForStereo[at_from2] ++;
                            pCS->bAtomUsedForStereo[at_to1] --;
                            pCS->bAtomUsedForStereo[at_to2] --;
                            ret2 = map_stereo_bonds4 ( at, num_atoms, num_at_tg, num_max, bAllene, nCanonRankFrom, nAtomNumberCanonFrom, nCanonRankTo,
                                                      nSymmRank, pRankStack1+nStackPtr[istk3],  pRankStack2+nStackPtr[istk3],
                                                      nTempRank, nNumMappedRanks[istk3], nSymmStereo, NeighList,
                                                      pCS, cur_tree, nNumMappedBonds+1 ,
                                                      vABParityUnknown);
                            pCS->bRankUsedForStereo[at_from1] --;
                            pCS->bRankUsedForStereo[at_from2] --;
                            pCS->bAtomUsedForStereo[at_to1] ++;
                            pCS->bAtomUsedForStereo[at_to2] ++;
                            if ( ret2 == 4 ) {
                                if ( nNumMappedBonds ) {
                                    return ret2;
                                } else {
                                    pCS->bFirstCT = 1;
                                    goto total_restart;
                                }
                            }
                            if ( RETURNED_ERROR(ret2) ) {
                                if ( ret2 == CT_TIMEOUT_ERR )
                                    return ret2;
                                else
                                    return ret2; /*  program error */
                            }
                            if ( ret2 > 0 ) {
                                nTotSuccess |= 1;
                                nNumAtTo1Success ++;
                                if ( bStereoIsBetterWasSetHere || (ret2 & 2) ) {
                                    CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                                    nTotSuccess |= 2; /*  Obtained a smaller CT */
                                }
                                at_no_n1_num_success ++;
                            } else {
                                if ( bStereoIsBetterWasSetHere ) {  /*  rollback */
                                    pCS->bStereoIsBetter = 0;
                                    pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond2;
                                }
                                if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                                     1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ) ) {
                                    CurTreeRemoveIfLastAtom( cur_tree, at_to_n2 );
                                    CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                                }
                            }
                            bStereoIsBetterWasSetHere = 0;
                        }
                    } /*  end choices in mapping neighbors of the 2nd half-bond */
                    if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                         1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n2 ) ) {
                         CurTreeRemoveLastRank( cur_tree );
                    }
                    /* added 2006-07-20 */
                    if ( !at_no_n1_num_success && tpos1 < CurTreeGetPos( cur_tree ) &&
                        1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n1 ) ) {
                         CurTreeRemoveIfLastAtom( cur_tree, at_to_n1 );
                    }
                         
                } /*  end choices in mapping neighbors of the 1st half-bond */
                if (  tpos1 < CurTreeGetPos( cur_tree ) &&
                      1 == CurTreeIsLastRank( cur_tree, at_rank_canon_n1 ) ) {
                     CurTreeRemoveLastRank( cur_tree );
                }
            } /*  end of choices in mapping at_from2 */
            if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon2 ) ) {
                 CurTreeRemoveLastRank( cur_tree );
            }
            if ( !nNumAtTo1Success ) {
                if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                     1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
                    CurTreeRemoveIfLastAtom( cur_tree, at_to1 );
                    CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                }
            }
            if ( bAllParitiesIdentical /*&& !nSymmStereo*/ ) {
                break;
            }
        } /*  end of choices in mapping at_from1 */
        
        if ( tpos1 < CurTreeGetPos( cur_tree ) &&
             1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
            CurTreeRemoveLastRank( cur_tree );
        } else
        /*  CurTree consistecy check (debug only) */
        if ( tpos1 != CurTreeGetPos( cur_tree ) ) {
            return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
        }

        if ( !nTotSuccess || stereo_bond_parity == sb_parity_calc ) {
            goto repeat_all; /*  repeat with next parity if no success or with the same parity, now known */
        }

        /*  Previously the control flow never came here... */
        if ( !nTotSuccess ) {
            pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond;
            CurTreeSetPos( cur_tree, tpos1 );
            /*  Occurs when atoms are not really equvalent ( -O= without positive charge in "aromatic" ring) */
            return 0; /* Happens for ID=92439,100318,100319 when EXCL_ALL_AROM_BOND_PARITY=0 and
                       * nNumChoices=0.
                       * Results from impossible previous mapping of symmetric relatively
                       * to a central ring aromatic circles while central ring is not symmetrical due to
                       * alternate bonds (in the central ring number of pi-electrons, atoms and bonds
                       * are symmetrical).
                       * Does not happen when alternate bonds of the central ring
                       * are treated as aromatic by attaching a (+) charge to the oxygen.
                       */
        }
    } else

    {
        int ret;

        if ( !nNumMappedBonds ) {
            pCS->bStereoIsBetter = 0;  /*  the first stereo feature in the canonical CT has not been processed yet */
        }

        if ( nNumMappedBonds < pCS->nLenLinearCTStereoDble ) {
            prevBond = pCS->LinearCTStereoDble[nNumMappedBonds];
        }

        /*  all stereo bonds have been mapped; now start processing stereo atoms... */
        ret = map_stereo_atoms4 ( at, num_atoms, num_at_tg, num_max, nCanonRankFrom, nAtomNumberCanonFrom, nCanonRankTo,
                        nSymmRank, pRankStack1,  pRankStack2, nTempRank, nNumMappedRanksInput,
                        nSymmStereo,  NeighList, pCS, cur_tree,  0 , vABParityUnknown);
        if ( ret == 4 ) {
            if ( nNumMappedBonds ) {
                return ret;
            } else {
                pCS->bFirstCT = 1;
                goto total_restart;
            }
        }
        if ( RETURNED_ERROR(ret) ) {
            if ( ret == CT_TIMEOUT_ERR )
                return ret;
            else
                return ret; /*  program error */
        }
        if ( ret > 0 ) {
            nTotSuccess |= 1;
            if ( ret & 2 ) {
                CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                nTotSuccess |= 2; /*  Obtained a smaller CT */
            }
        }
    }
    if ( !nTotSuccess && pCS->nLenLinearCTStereoDble &&
         nNumMappedBonds < pCS->nLenLinearCTStereoDble ) {
        pCS->LinearCTStereoDble[nNumMappedBonds] = prevBond;
    }
    return nTotSuccess;  /*  ok */
}







/****************************************************************************************
 */
int map_stereo_atoms4 ( 
                sp_ATOM *at, int num_atoms, int num_at_tg, int num_max,
                const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, AT_RANK *nCanonRankTo, /*  canonical numbering to be mapped */
                const AT_RANK *nSymmRank,      AT_RANK **pRankStack1/*from*/, AT_RANK **pRankStack2/*to*/,
                AT_RANK *nTempRank,      int nNumMappedRanksInput,
                AT_RANK *nSymmStereo,    NEIGH_LIST *NeighList,
                CANON_STAT *pCS,         CUR_TREE *cur_tree, int nNumMappedAtoms ,
                int vABParityUnknown )
{
/*
 *   Do not check whether "from" atoms have any stereo features.
 */
    int            nTotSuccess = 0;
    AT_STEREO_CARB prevAtom;
    int            tpos1;

    tpos1 = CurTreeGetPos( cur_tree );

    if ( nNumMappedAtoms < pCS->nLenLinearCTStereoCarb ) {
        /* AT_RANK *nRankFrom=*pRankStack1++,  AT_RANK *nAtomNumberFrom=pRankStack1++; */
        /* AT_RANK *nRankTo  =*pRankStack2++,  AT_RANK *nAtomNumberTo  =pRankStack2++; */
        int j1, at_from1, at_to1, /*at_from2, at_to2,*/ iMax, lvl, bStereoIsBetterWasSetHere;
        int istk, istk2, bAddStack, nNumAtTo1Success, c, bFirstTime=1, bAllParitiesIdentical;
        EQ_NEIGH EN[5], *pEN;
        int nStackPtr[5], nMappedRanks[5], j[5], *nSP, *nMR, bLastLvlFailed;
        
        AT_RANK at_rank_canon1, cr[5], at_to[5];
        AT_RANK canon_rank1_min = 0;
        int at_rank1; /*  rank for mapping */
        int nNumChoices, nNumUnkn, nNumUndf, nNumWorse, nNumBest, nNumCalc;
        int stereo_center_parity, prev_stereo_center_parity, sb_parity_calc, pass;
        AT_STEREO_CARB prevAtom2;

        prevAtom = pCS->LinearCTStereoCarb[nNumMappedAtoms]; /*  save to restore in case of failure */
        at_rank_canon1 = nNumMappedAtoms? pCS->LinearCTStereoCarb[nNumMappedAtoms-1].at_num:0;

        goto bypass_next_canon_rank_check;

next_canon_rank:

        if ( !pCS->bStereoIsBetter /*??? && !pCS->bFirstCT ???*/ &&
              at_rank_canon1 >= pCS->LinearCTStereoCarb[nNumMappedAtoms].at_num) {
            /*  cannot find next available canonical number */
            if ( !nTotSuccess ) {
                pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom; /*  restore because of failure */
            }
            CurTreeSetPos( cur_tree, tpos1 );
            return nTotSuccess;
        }

bypass_next_canon_rank_check:

        CurTreeSetPos( cur_tree, tpos1 );

        /*  find next available canon. number for a stereogenic atom */
        if ( !Next_SC_At_CanonRank2( &at_rank_canon1, &canon_rank1_min, &bFirstTime,
                          pCS->bAtomUsedForStereo, pRankStack1, pRankStack2,
                          nAtomNumberCanonFrom, num_atoms ) ||
              !pCS->bStereoIsBetter && 
              at_rank_canon1 > pCS->LinearCTStereoCarb[nNumMappedAtoms].at_num) {
            /*  cannot find next available canonical number */
            if ( !nTotSuccess ) {
                pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom; /*  restore because of failure */
            }
            return nTotSuccess;
        }

        nNumChoices = 0;
        nNumUnkn    = 0;
        nNumUndf    = 0;
        nNumBest    = 0;
        nNumWorse   = 0;
        nNumCalc    = 0;
        pass        = 0;
        prev_stereo_center_parity = 0;

        /*  get mapping rank for the canon. number */
        at_rank1 = pRankStack1[0][at_from1=(int)nAtomNumberCanonFrom[at_rank_canon1 - 1]];
        iMax = at_rank1-1;
        /*  for debug only */
        if ( at_rank1 != pRankStack2[0][pRankStack2[1][at_rank1-1]] )
            return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */

        /*  count special parities of the not mapped yet "to" atoms */
        for ( j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1  =pRankStack2[1][iMax-j1]]; j1 ++ ) {
            if ( !at[at_to1].stereo_bond_neighbor[0] && pCS->bAtomUsedForStereo[at_to1] == STEREO_AT_MARK ) {
                int no_choice = 0;
                stereo_center_parity = PARITY_VAL(at[at_to1].stereo_atom_parity);
                switch(stereo_center_parity) {

                case AB_PARITY_UNDF: nNumUndf  ++; break; /*  4 */

                case AB_PARITY_UNKN: nNumUnkn  ++; 
                                     break; /*  3 */

                case BEST_PARITY:    nNumBest  ++; break; /*  1 */
                case WORSE_PARITY:   nNumWorse ++; break; /*  2 */
                case AB_PARITY_CALC: nNumCalc  ++; break;
                case AB_PARITY_NONE: no_choice ++; break; /*  0 */
                }
                nNumChoices += !no_choice;
            }
        }
        if ( nNumChoices != nNumCalc + nNumUndf + nNumUnkn + nNumBest + nNumWorse ) {
            return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
        }
        if ( !nNumChoices ) {
            goto next_canon_rank;
        }
        /*  Determine the first parity to search */
        sb_parity_calc = ( nNumCalc > 0 )? BEST_PARITY : 0;

        /*  ==============================================================
            Search sequence:           sb_parity_calc    stereo_center_parity
            ==============================================================
            BEST_PARITY   (calc)       BEST_PARITY     BEST_PARITY
            BEST_PARITY   (known)      BEST_PARITY     WORSE_PARITY  or 0
            WORSE_PARITY  (calc)       WORSE_PARITY    WORSE_PARITY
            WORSE_PARITY  (known)      WORSE_PARITY    0
            AB_PARITY_UNKN(known)      AB_PARITY_UNKN  0
            AB_PARITY_UNDF(known)      AB_PARITY_UNDF  0

            if (sb_parity_calc==stereo_center_parity) then "calc" else "known"
         */

repeat_all:

        if ( !pass ++ ) {
            /*  select the smallest parity to search */
            if ( sb_parity_calc ) {
                stereo_center_parity = BEST_PARITY;
            } else {
                stereo_center_parity = nNumBest?   BEST_PARITY :
                                       nNumWorse?  WORSE_PARITY :
                                       nNumUnkn?   AB_PARITY_UNKN :
                                       nNumUndf?   AB_PARITY_UNDF : AB_PARITY_NONE;
            }
        } else {
            prev_stereo_center_parity = stereo_center_parity;
            j1 = NextStereoParity2Test( &stereo_center_parity, &sb_parity_calc,
                                     nNumBest, nNumWorse, nNumUnkn, nNumUndf, nNumCalc, 
                                     vABParityUnknown );
            switch ( j1 ) {
            case 0:
                break; /* obtained next parity to test */
            case 1:
                goto next_canon_rank;
            default:
                return j1; /* program error */
            }
        }
        if ( stereo_center_parity == AB_PARITY_NONE ) {
            /*  error? */
            return CT_STEREOCOUNT_ERR;                  /*   <BRKPT> */
        }
        /*  check if the new requested parity is small enough */
        if ( !pCS->bStereoIsBetter ) {
            c = CompareLinCtStereoAtomToValues( nTotSuccess? pCS->LinearCTStereoCarb+nNumMappedAtoms : &prevAtom,
                              at_rank_canon1, (U_CHAR)stereo_center_parity );
            if ( c < 0 ) {
                if ( !nTotSuccess ) {
                    pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom;
                }
                CurTreeSetPos( cur_tree, tpos1 );
                return nTotSuccess;
            }
        }


        bAllParitiesIdentical     = 0;
        bStereoIsBetterWasSetHere = 0;
        istk = istk2              = 0;
        CurTreeSetPos( cur_tree, tpos1 );  /*  start over */
        /*
        if ( prev_stereo_center_parity != stereo_center_parity ) {
            CurTreeSetPos( cur_tree, tpos1 );
        }
        */                                      /*  nRankTo                 nAtomNumberTo */
        for ( j1 = 0; j1 <= iMax && at_rank1 == pRankStack2[0][at_to1  =pRankStack2[1][iMax-j1]]; j1 ++ ) {
            int ret, ret1, ret2, parity1;
            nNumAtTo1Success = 0;
            /*
            if ( !(at[at_to1].stereo_atom_parity && !at[at_to1].stereo_bond_neighbor[0] &&
                   pCS->bAtomUsedForStereo[at_to1] == STEREO_AT_MARK ) )
            */
            if ( !at[at_to1].stereo_atom_parity || at[at_to1].stereo_bond_neighbor[0] ||
                  pCS->bAtomUsedForStereo[at_to1] != STEREO_AT_MARK ) /* simplify 12-17-2003 */
                continue;  
                           /* Do not map on non-stereogenic atom constitutionally
                            * equivalent to a steregenic atom. Here
                            * at[at_to1] is not a sterereo center;  |       |
                            * bonds tautomerism is a usual cause.  -P(+)-CH=P-
                            * For example, consider a fragment:     |       |
                            * The two atoms P may be constitutionally
                            * equivalent, P(+) may be seen as a stereocenter
                            * while another P has a double bond (Now such a P(V) IS a stereocenter).
                           */
            /*  check whether the stereocenter parity corresponds to the requested stereocenter parity */
            if ( PARITY_KNOWN(at[at_to1].stereo_atom_parity) ) {
                if ( stereo_center_parity == sb_parity_calc ) {
                    continue;  /*  requested parity to be calculated, found known parity */
                }
                if ( stereo_center_parity != PARITY_VAL(at[at_to1].stereo_atom_parity) ) {
                    continue;  /*  parity differs from the requested parity */
                }
            } else
            if ( PARITY_CALCULATE( at[at_to1].stereo_atom_parity) ) {
                if ( stereo_center_parity != sb_parity_calc ) {
                    continue;  /*  requested known parity, found patity to be calculated */
                }
            } else {
                return CT_STEREOCOUNT_ERR;  /*  unknown parity type */
            }
            
            bAllParitiesIdentical = (( at[at_to1].stereo_atom_parity & KNOWN_PARITIES_EQL ) &&
                                     PARITY_KNOWN(at[at_to1].stereo_atom_parity));

            if ( !bAllParitiesIdentical && !nNumCalc &&
                 (!nNumUndf + !nNumUnkn + !nNumBest + !nNumWorse)==3   ) {
                /* only one kind of stereocenter parity is present; check whether all parities are really same */
                bAllParitiesIdentical = All_SC_Same(  at_rank_canon1, /*  canonical number */
                                                      pRankStack1, pRankStack2,
                                                      nAtomNumberCanonFrom, at );
                if ( bAllParitiesIdentical < 0 ) {
                    return CT_STEREOCOUNT_ERR;
                }
            }
            if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                 1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 )  &&
                 1 == CurTreeIsLastAtomEqu( cur_tree, at_to1, nSymmStereo ) )
                continue;

            /*  initialize stack pointer nStackPtr[istk] for "hand-made" recursion */
            /*  stacks are pRankStack1[], pRankStack2[], nNumMappedRanks[] */
            istk               = 0;
            nStackPtr[istk]    = 0;
            nMappedRanks[istk] = nNumMappedRanksInput;
            bAddStack          = 0;
            /*  if all equivalent atoms have same known parity, do not map any of them here */
            if ( !bAllParitiesIdentical ) {
                /*  map the central atom */
                /*  this mapping is always possible */
                ret1 = map_an_atom2( num_at_tg, num_max, at_from1, at_to1,
                                    nTempRank, nMappedRanks[istk], &nMappedRanks[istk+1], pCS,
                                    NeighList, pRankStack1+nStackPtr[istk], pRankStack2+nStackPtr[istk],
                                    &bAddStack );
                if ( RETURNED_ERROR(ret1) ) {
                    return ret1; /*  error */
                }
                nStackPtr[istk+1] = nStackPtr[istk] + bAddStack;
                istk ++;
            } else {
                ClearPreviousMappings( pRankStack1+2 ); /*  precaution */
            }
            
            /*********************************************************************************
             *
             *   Unknown Stereocenter Parity case: possibly need to map stereo center neighbors
             */
            if ( stereo_center_parity == sb_parity_calc )
            {
                /*  find out the parity */
                parity1 = parity_of_mapped_atom2( at_from1, at_to1, at, &EN[istk],
                                                 nCanonRankFrom, pRankStack1[nStackPtr[istk]],
                                                 pRankStack2[nStackPtr[istk]] );
                /*  if parity is well-defined then returned EN[istk].num_to=0 */
                if ( !parity1 ) {
                    return CT_STEREOCOUNT_ERR; /*  program error */ /*   <BRKPT> */
                }
                if ( !EN[istk].num_to && parity1 != sb_parity_calc ) {
                    continue; /*  looking for the parity value = sb_parity_calc */
                }

            } else {
                /*  Known parity */
                parity1 = stereo_center_parity;
                EN[istk].num_to = 0;
            }

            /***********************************************************************
             * no need to map the neighbors: parity is known or has been calculated
             */
            if ( stereo_center_parity == sb_parity_calc && !EN[istk].num_to ||
                 /*  now well-defined, but unknown in advance atom parity  OR   */
                 stereo_center_parity != sb_parity_calc )
                 /*  known in advance parity = stereo_center_parity */
            {
                /*  do not need to map the neighbors */
                c = CompareLinCtStereoAtomToValues( pCS->LinearCTStereoCarb+nNumMappedAtoms,
                                            at_rank_canon1, (U_CHAR)parity1 );
                if ( c < 0 && !pCS->bStereoIsBetter ) {
                    /*  reject */
                    pCS->lNumRejectedCT ++;
                    continue;  /*  Reject: not a minimal CT. Should not happen */
                }  else  {
                    /*  accept */

                    if ( bAddStack ) {
                        if ( tpos1 == CurTreeGetPos( cur_tree ) ||
                             0 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
                            CurTreeAddRank( cur_tree, at_rank_canon1 );
                        }
                        CurTreeAddAtom( cur_tree, at_to1 );
                    }

                    if ( c > 0 && !pCS->bStereoIsBetter ) {
                        /*  stereo center entry is less than the previusly found */
                        pCS->bStereoIsBetter = bStereoIsBetterWasSetHere = 1;
                        prevAtom2 = pCS->LinearCTStereoCarb[nNumMappedAtoms];
                    }
                    pCS->LinearCTStereoCarb[nNumMappedAtoms].parity = parity1;
                    pCS->LinearCTStereoCarb[nNumMappedAtoms].at_num = at_rank_canon1;
                    pCS->bRankUsedForStereo[at_from1] = 3;
#if( FIX_ChCh_STEREO_CANON_BUG == 1 )
                    if ( !bAllParitiesIdentical )
#endif
                        pCS->bAtomUsedForStereo[at_to1] -= STEREO_AT_MARK;

                    ret = map_stereo_atoms4 ( at, num_atoms, num_at_tg, num_max, nCanonRankFrom, nAtomNumberCanonFrom, nCanonRankTo,
                                       nSymmRank, pRankStack1+nStackPtr[istk],  pRankStack2+nStackPtr[istk],
                                       nTempRank,  nMappedRanks[istk],  nSymmStereo,  NeighList,
                                       pCS, cur_tree, nNumMappedAtoms+1 ,
                                       vABParityUnknown);
                    pCS->bRankUsedForStereo[at_from1] = 0;
#if( FIX_ChCh_STEREO_CANON_BUG == 1 )
                    if ( !bAllParitiesIdentical )
#endif
                        pCS->bAtomUsedForStereo[at_to1] += STEREO_AT_MARK;
                    if ( ret == 4 ) {
                        return ret;
                    }
                    if ( RETURNED_ERROR(ret) ) {
                        if ( ret == CT_TIMEOUT_ERR )
                            return ret;
                        else
                            return ret; /*  program error */
                    }
                    if ( ret > 0 ) {
                        nTotSuccess |= 1;
                        nNumAtTo1Success ++;
                        if ( bStereoIsBetterWasSetHere || (ret & 2) ) {
                            CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                            nTotSuccess |= 2; /*  Obtained a smaller CT */
                        }
                    } else {
                        if ( bStereoIsBetterWasSetHere ) {
                            pCS->bStereoIsBetter = 0;
                            pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom2;
                        }
                        /*  remove failed atom1 from the tree */
                        if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                             1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
                            CurTreeRemoveIfLastAtom( cur_tree, at_to1 );
                            CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                        }

                    }
                    bStereoIsBetterWasSetHere = 0;
                }
                /*
                if ( (at[at_to1].stereo_atom_parity & KNOWN_PARITIES_EQL ) &&
                     ATOM_PARITY_KNOWN(stereo_center_parity) && !nSymmStereo ) { // ??? add && !nSymmStereo ???
                    break; // do not repeat for the same kind of stereo atom with the parity known in advance
                }
                */
                if ( bAllParitiesIdentical ) {
                    break;  /*  do not repeat for the same kind of stereo atom with the parity known in advance */
                }
                continue;
                
            }
            
            /***************************************************
             *
             * Need to map the neighbors
             */
            if ( stereo_center_parity != sb_parity_calc ) {
                return CT_STEREOCOUNT_ERR;  /*  program error */ /*   <BRKPT> */
            }
            /* -- has already been calculated --
            parity1 = parity_of_mapped_atom2( at_from1, at_to1, at, &EN[istk],
                                             nCanonRankFrom, pRankStack1[nStackPtr[istk]], pRankStack2[nStackPtr[istk]] );
            */
            if ( !parity1 ) {
                return CT_STEREOCOUNT_ERR; /*  1/25/2002 */ /*   <BRKPT> */
            }
            
            if ( bAddStack ) {
                if ( tpos1 == CurTreeGetPos( cur_tree ) ||
                     0 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
                    CurTreeAddRank( cur_tree, at_rank_canon1 );
                }
                CurTreeAddAtom( cur_tree, at_to1 );
            }
            /******************************************************
             * Need to fix the neighbors to define the atom parity
             ******************************************************/
            /*  a recursion replaced with the hand-made stack */
            lvl = 0;              /*  the "recursion" depth level */
            nSP = &nStackPtr[istk];
            nMR = &nMappedRanks[istk];
            pEN = &EN[istk];
            bLastLvlFailed = 0;

            /*  entering "recursion" depth level lvl */
next_lvl:            
            if ( pEN[lvl].num_to ) {
                /* Found tied neighbors. Try all transpositions of the tied neighbors.
                 * j is a number of the "to" tied neighbor in the pEN[lvl].to_at[*] to
                 * which the pEN[lvl].from_at "from" neighbor's canonical number is mapped
                 */
                j[lvl] = 0;
next_j:         
                cr[lvl]    = nCanonRankFrom[pEN[lvl].from_at];
                at_to[lvl] = pEN[lvl].to_at[j[lvl]];
                if ( j[lvl] &&  tpos1 < CurTreeGetPos( cur_tree ) &&
                     1 == CurTreeIsLastRank( cur_tree, cr[lvl] ) &&
                     1 == CurTreeIsLastAtomEqu( cur_tree, at_to[lvl], nSymmStereo ) ) {
                    lvl ++;
                    bLastLvlFailed = 0;
                    goto backup; /*  do not test stereo equivalent atoms except the first one */
                }

                ret2 = map_an_atom2( num_at_tg, num_max,
                                    pEN[lvl].from_at,        /* from */
                                    pEN[lvl].to_at[j[lvl]],  /* to */
                                    nTempRank, nMR[lvl], &nMR[lvl+1], pCS,
                                    NeighList, pRankStack1+nSP[lvl], pRankStack2+nSP[lvl],
                                    &bAddStack );
                if ( RETURNED_ERROR(ret2) ) {
                    return ret2; /*  program error */
                }

                /*  next recursion depth level */
                if ( bAddStack ) {
                    if ( tpos1 == CurTreeGetPos( cur_tree ) ||
                         0 == CurTreeIsLastRank( cur_tree, cr[lvl] ) ) {
                        CurTreeAddRank( cur_tree, cr[lvl] );
                    }
                    CurTreeAddAtom( cur_tree, at_to[lvl] );
                }
                nSP[lvl+1] = nSP[lvl] + bAddStack;
                lvl ++; /*  upon increment lvl = number of additionally mapped neighbors
                         *  (entering next recursion level) */
                /*  check if the mapping has defined the parity */
                parity1 = parity_of_mapped_atom2( at_from1, at_to1, at, &pEN[lvl],
                                             nCanonRankFrom, pRankStack1[nSP[lvl]], pRankStack2[nSP[lvl]] );
                if ( !parity1 ) {
                    return CT_STEREOCOUNT_ERR; /*  1/25/2002 */ /*   <BRKPT> */
                }
                if ( parity1 < 0 ) {
                    goto next_lvl; /*  we need at least one more mapping to define the parity */
                }
                
                /**********************************************************
                 *
                 *  Check the parity
                 *
                 **********************************************************
                 *  make a decision whether to accept the current mapping */

                c = CompareLinCtStereoAtomToValues( pCS->LinearCTStereoCarb+nNumMappedAtoms,
                                            at_rank_canon1, (U_CHAR)parity1 );
                if ( sb_parity_calc != parity1 ||
                     c < 0 && !pCS->bStereoIsBetter ) {
                    pCS->lNumRejectedCT ++;
                    bLastLvlFailed = 1;
                }  else
                /*  the parity has been defined (all neighbors have untied ranks) */
                /*  if ( bAcceptAllParities || parity1 == BEST_PARITY ) */
                {
                    /*********************************************************************
                     *
                     * Process the parity here. We are at the top of the recursion stack.
                     *
                     *********************************************************************/
                    /*  try to accept current neighbors mapping */
                    if ( c > 0  && !pCS->bStereoIsBetter ) {
                        pCS->bStereoIsBetter = bStereoIsBetterWasSetHere = 1;
                        prevAtom2 = pCS->LinearCTStereoCarb[nNumMappedAtoms];
                    }
                    pCS->LinearCTStereoCarb[nNumMappedAtoms].parity = parity1;
                    pCS->LinearCTStereoCarb[nNumMappedAtoms].at_num = at_rank_canon1;
                    pCS->bRankUsedForStereo[at_from1] = 3;
                    pCS->bAtomUsedForStereo[at_to1] -= STEREO_AT_MARK;

                    ret = map_stereo_atoms4 ( at, num_atoms, num_at_tg, num_max, nCanonRankFrom, nAtomNumberCanonFrom, nCanonRankTo,
                                       nSymmRank, pRankStack1+nSP[lvl],  pRankStack2+nSP[lvl],
                                       nTempRank,  nMR[lvl],  nSymmStereo,  NeighList,
                                       pCS, cur_tree, nNumMappedAtoms+1 ,
                                       vABParityUnknown );
                    pCS->bRankUsedForStereo[at_from1] = 0;
                    pCS->bAtomUsedForStereo[at_to1] += STEREO_AT_MARK;
                    if ( ret == 4 ) {
                        return ret;
                    }
                    if ( RETURNED_ERROR(ret) ) {
                        if ( ret == CT_TIMEOUT_ERR )
                            return ret;
                        else
                            return ret; /*  program error */
                    }
                    if ( ret > 0 ) {
                        nTotSuccess |= 1;
                        nNumAtTo1Success ++;
                        if ( bStereoIsBetterWasSetHere || (ret & 2) ) {
                            CurTreeKeepLastAtomsOnly( cur_tree, tpos1, 1 );  /*  start over */
                            nTotSuccess |= 2; /*  Obtained a smaller CT */
                        }
                    } else {
                        if ( bStereoIsBetterWasSetHere ) {
                            pCS->bStereoIsBetter = 0;
                            pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom2;
                        }
                        bLastLvlFailed = 1;
                    }
                    bStereoIsBetterWasSetHere = 0;

                    /*  avoid redundant repetitions: */
                    /*  check if neighbors mappings have altered another stereo center parity */
                    if ( !nSymmStereo && !might_change_other_atom_parity( at, num_atoms, at_to1,
                                       pRankStack2[nSP[lvl]] /* ranks after neigbors mapping */,
                                       pRankStack2[nStackPtr[istk]] /* ranks before the mapping neighbors */) ) {
                        goto done;
                    }
                }
                /*  Continue the cycle. Go to the previous "recursion" level */
backup:
                while (lvl -- > 0 ) {

                    j[lvl] ++; /*  next neighbor at this level */
                    if ( j[lvl] < pEN[lvl].num_to ) {
                        if ( bLastLvlFailed ) {
                            if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                                 1 == CurTreeIsLastRank( cur_tree, cr[lvl] ) ) {
                                CurTreeRemoveIfLastAtom( cur_tree, at_to[lvl] );
                                CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                            }
                            bLastLvlFailed = 0;
                        }
                        /*  Done with this level. Go back one level */
                        goto next_j;
                    }
                    /*  remove failed atom from the tree */
                    if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                         1 == CurTreeIsLastRank( cur_tree, cr[lvl] ) ) {
                        CurTreeRemoveLastRank( cur_tree );
                    }
                }
                goto done;
            } else {
                cr[lvl] = 0;
            }

done:;      /*  at this point lvl=0. */
            if ( !nNumAtTo1Success ) {
                if ( tpos1 < CurTreeGetPos( cur_tree ) &&
                     1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
                    CurTreeRemoveIfLastAtom( cur_tree, at_to1 );
                    CurTreeRemoveLastRankIfNoAtoms( cur_tree );
                }
            }
        } /*  end of stereo atom mapping cycle */

        if ( tpos1 < CurTreeGetPos( cur_tree ) &&
             1 == CurTreeIsLastRank( cur_tree, at_rank_canon1 ) ) {
            CurTreeRemoveLastRank( cur_tree );
        } else
        /*  CurTree consistency check (debug only) */
        if ( tpos1 != CurTreeGetPos( cur_tree ) ) {
            return CT_STEREOCOUNT_ERR;  /*   <BRKPT> */
        }

        if ( !nTotSuccess || stereo_center_parity == sb_parity_calc ) {
            goto repeat_all; /*  repeat with next parity if no success or with the same parity, now known */
        }

    } else {

        /****************************************************
         *
         *  All stereogenic atoms and bonds have been mapped 
         *
         ****************************************************/

        if ( UserAction && USER_ACTION_QUIT == (*UserAction)()  ||
             ConsoleQuit && (*ConsoleQuit)() ) {
            return CT_USER_QUIT_ERR;
        }

        if ( pCS->bStereoIsBetter || pCS->bFirstCT ) {
            /* All stereo atoms have been mapped. Current stereo name is better than all previous.
             * Create new numbering for the new CT
             * break all remaining "from" ties
             */
            int i1, ret;
            AT_RANK rc, n1, n2;
            ret=BreakAllTies( num_at_tg, num_max, pRankStack1, NeighList, nTempRank, pCS);
            if ( RETURNED_ERROR( ret ) ) {
                return ret;
            }
            /*  break all remaining "from" ties */
            ret=BreakAllTies( num_at_tg, num_max, pRankStack2, NeighList, nTempRank, pCS);
            if ( RETURNED_ERROR( ret ) ) {
                return ret;
            }
            /*  move stack pointers to the "nAtomNumber[*]" after all ties are broken */
            pRankStack1 += 2;
            pRankStack2 += 2;
            /* Now final mapping ranks of "to" atom (*pRankStack2)[i] and "from" atom (*pRankStack1)[i]
             * are equal and all ranks are different, that is, we have a full mapping
             * Copy so far best canonical numbering from "from" to "to".
             */
            memset( pCS->nPrevAtomNumber, 0, num_at_tg*sizeof(pCS->nPrevAtomNumber[0]) );
            for ( i1 = 0; i1 < num_at_tg; i1 ++ ) {
                n1 = pRankStack1[1][i1];
                rc = nCanonRankFrom[n1]; /*  new canon. rank */
                n2 = pRankStack2[1][i1];                  /*  orig. atom number */
                nCanonRankTo[n2] = rc;                    /*  assign new canon. number to the atom */
                /*  use this array to find stereo-equivalent atoms */
                pCS->nPrevAtomNumber[rc-1] = n2; /*  ord. number of the atom having canon. rank = rc */
                nSymmStereo[i1] = i1;            /*  restart search for stereo equivalent atoms */
                /* check mapping correctness */
                if ( pRankStack1[0][n1] != pRankStack2[0][n2] ||
                     nSymmRank[n1]      != nSymmRank[n2] ) {
                    return CT_STEREO_CANON_ERR; /* stereo mapping error */
                }
            }
            /*  statistics */
            pCS->lNumTotCT ++;
            pCS->lNumEqualCT = 1;
            pCS->lNumDecreasedCT ++;
            pCS->bStereoIsBetter = 0; /*  prepare to start over */
            nTotSuccess = 1;
            pCS->bFirstCT = 0;
#if( REMOVE_CALC_NONSTEREO == 1 ) /* { */
            if ( !(pCS->nMode & CMODE_REDNDNT_STEREO ) ) {
                i1 = RemoveCalculatedNonStereo( at, num_atoms, num_at_tg,
                                  pRankStack1, pRankStack2, nTempRank, NeighList,
                                  nSymmRank, nCanonRankTo, pCS->nPrevAtomNumber, pCS,
                                  vABParityUnknown);
                if ( RETURNED_ERROR( i1 ) ) {
                    return i1;
                }
                if ( i1 < 0 ) {
#if( bRELEASE_VERSION == 0 )
                    pCS->bExtract |= EXTR_REMOVE_PARITY_WARNING;
#endif
                    i1 = -(1+i1);
                }
                if ( i1 > 0 ) {
                    return 4; /*  total restart: due to newly found stereo equivalence */
                              /*  the length of the stereo CT has changed */
                }
            }
#endif /* } REMOVE_CALC_NONSTEREO */
            pRankStack1 -= 2;
            pRankStack2 -= 2;
        } else {
            /*  current stereo name is same as previous. We do not need a full mapping. */
            if ( nSymmStereo ) {
                int num_changes = 0;
                AT_RANK r, n1, n2, r_max, cr;
                r_max = (AT_RANK)num_at_tg;
                for ( r = 1; r <= r_max; r ++ ) {
                    if ( bUniqueAtNbrFromMappingRank( pRankStack1, r, &n1 ) ) {
                        if ( bUniqueAtNbrFromMappingRank( pRankStack2, r, &n2 ) ) {
                            /*  atoms at[n1], at[n2] have identical untied mapping rank r */
                            cr = nCanonRankFrom[(int)n1]-1; /*  (new at[n2] canonical rank)-1 */
                            /*  pCS->nPrevAtomNumber[(int)cr] = */
                            /*    previous ordering number of an atom with the canon. rank = cr+1; */
                            /*    make this atom equivalent to atom at[n2]: */
                            num_changes += nJoin2Mcrs( nSymmStereo, pCS->nPrevAtomNumber[(int)cr], n2 );
                        } else {
                            return CT_MAPCOUNT_ERR; /*  mapping ranks must be either both tied or untied. */ /*   <BRKPT> */
                        }
                    }
                }
                if ( num_changes ) { /*  compress trees to stars */
                    for ( r = r_max-1; r; r -- ) {
                        nGetMcr( nSymmStereo, r );
                    }
                }
            }
            /*  statistics */
            pCS->lNumEqualCT ++;
            pCS->lNumTotCT ++;
            nTotSuccess = 1;
        }
        if ( bInchiTimeIsOver( pCS->ulTimeOutTime ) ) {
            return CT_TIMEOUT_ERR;
        }
    }
    if ( !nTotSuccess && nNumMappedAtoms < pCS->nLenLinearCTStereoCarb ) {
        pCS->LinearCTStereoCarb[nNumMappedAtoms] = prevAtom;
        CurTreeSetPos( cur_tree, tpos1 );
    }
    return nTotSuccess;  /*  return to the previous level of the recursion. */
}
