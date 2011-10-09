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

/*#define CHECK_WIN32_VC_HEAP*/
#include "mode.h"

#if ( READ_INCHI_STRING == 1 )

#include "ichi.h"
#include "ichitime.h"

#include "inpdef.h"
#include "ichimain.h"
#include "ichierr.h"
#include "incomdef.h"
#include "ichiring.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichinorm.h"
#include "util.h"

#include "ichicomp.h"
#include "ichister.h"

#include "ichi_bns.h"

#include "strutil.h"

#include "ichirvrs.h"


#define INC_ADD_EDGE 64
/***********************************************************************************************/
int FixRestoredStructureStereo( INCHI_MODE cmpInChI, ICR *icr, INCHI_MODE cmpInChI2, ICR *icr2,
                        ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
                        StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA,
                        ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ppt_group_info, inp_ATOM **ppat_norm,
                        inp_ATOM **ppat_prep, INChI *pInChI[], long num_inp,
                        int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask, int forbidden_stereo_edge_mask)
{
    /*--------- process extra or missing Fixed-H on non-tautomeric atoms ------*/
    /* at2 should be the most recently restored atom, Fixed-H */
    int i, j, k, delta, tot_succes, max_success, cur_success, ret = 0;
    int err, iOrigInChI, iRevrInChI;
    int j12, v1, v2, e, vRad;
    BNS_VERTEX *pv1, *pv2, *pvRad;
    BNS_EDGE   *pe, *peRad;
    EDGE_LIST AllChargeEdges, CurrEdges, NFlowerEdges, OtherNFlowerEdges, FixedStereoEdges, AllRadList;
    EDGE_LIST TautMinusEdges[2]; /* 0 -> O & O(+), 1=> N & N(+) */

    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
    INChI_Stereo *pStereoInChI, *pStereo2InChI, *pStereoRevrs, *pStereo2Revrs; 
    
    /* Stereo */

    /* currently being processed layer */
    pStereoInChI = (pInChI[0]->StereoIsotopic &&
               pInChI[0]->StereoIsotopic->nNumberOfStereoBonds +
               pInChI[0]->StereoIsotopic->nNumberOfStereoCenters)?
               pInChI[0]->StereoIsotopic : pInChI[0]->Stereo;

    /* mobile-H layer in case of Fixed-H */
    pStereo2InChI = (pStruct->bMobileH == TAUT_YES || !pInChI[1] ||
               !pInChI[1]->nNumberOfAtoms || pInChI[1]->bDeleted)?
                        NULL:
               (pInChI[1]->StereoIsotopic &&
               pInChI[1]->StereoIsotopic->nNumberOfStereoBonds +
               pInChI[1]->StereoIsotopic->nNumberOfStereoCenters)?
                        pInChI[1]->StereoIsotopic :
                        pInChI[1]->Stereo;

    /* currently being processed layer */
    pStereoRevrs = (pStruct->pOneINChI[0]->StereoIsotopic &&
               pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
               pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters)?
               pStruct->pOneINChI[0]->StereoIsotopic : pStruct->pOneINChI[0]->Stereo;

    /* mobile-H layer in case of Fixed-H */
    pStereo2Revrs = (pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
               !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted)?
                        NULL:
               (pStruct->pOneINChI[1]->StereoIsotopic &&
               pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
               pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters)?
                        pStruct->pOneINChI[1]->StereoIsotopic :
                        pStruct->pOneINChI[1]->Stereo;
    
    INCHI_HEAPCHK

    AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &CurrEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &NFlowerEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &FixedStereoEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &AllRadList, EDGE_LIST_CLEAR );

    AllocEdgeList( TautMinusEdges+0, EDGE_LIST_CLEAR );
    AllocEdgeList( TautMinusEdges+1, EDGE_LIST_CLEAR );

    cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
    if ( cmpInChI & IDIF_PROBLEM ) {
        ret = RI_ERR_PROGR; /* severe restore problem */
        goto exit_function;
    }
    if ( err ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }

    cmpInChI2 = 0;

    if ( pStruct->bMobileH == TAUT_NON ) {
        /* these indexes are used to compare Mobile-H InChI */
        iOrigInChI = (pInChI[1] && pInChI[1]->nNumberOfAtoms && !pInChI[1]->bDeleted)? 1 : 0;
        iRevrInChI = (pStruct->pOneINChI[1] &&pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted)? 1 : 0;
    } else {
        iOrigInChI = 0;
        iRevrInChI = 0;
    }

    memset ( icr2, 0, sizeof(*icr2) );
    if ( iRevrInChI || iOrigInChI ) {
        /* additional mobile-H compare in case of Fixed-H */
        cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
        if ( cmpInChI & IDIF_PROBLEM ) {
            ret = RI_ERR_PROGR; /* severe restore problem */
            goto exit_function;
        }
        if ( err ) {
            ret = RI_ERR_ALLOC;
            goto exit_function;
        }
    }


    if ( !(cmpInChI & IDIFF_SB) && !(cmpInChI2 & IDIFF_SB) ) {
        goto exit_function;
    }
    /* need to temporarily remove fixing of stereogenic bonds */
    for ( i = 0; i < pStruct->num_atoms; i ++ ) {
        pv1 = pBNS->vert + i;
        for ( j = 0; j < at2[i].valence; j ++ ) {
            pe = pBNS->edge + (e=pv1->iedge[j]);
            if ( j == pe->neighbor1 ) {
                /* do not store same bond 2 times */
                if ( (pe->forbidden & forbidden_stereo_edge_mask) &&
                     (ret = AddToEdgeList( &FixedStereoEdges, e, INC_ADD_EDGE ) ) ) {
                    goto exit_function;
                }
            }
        }
    }


    tot_succes  = 0;
    cur_success = 0;
    if ( (cmpInChI & IDIF_SB_MISS) && (!cmpInChI2 || (cmpInChI2 & IDIF_SB_MISS)) &&
         0 < (max_success = pBNS->tot_st_cap - pBNS->tot_st_flow) ) {
        /*----------------------------------------------------*/
        /* case 01: extra stereogenic bond, radical present   */
        /* X=N-O*  => X=N=O and eliminate radical             */
        /*----------------------------------------------------*/
        int aN;
        BNS_VERTEX *pvO, *pvN;
        BNS_EDGE   *peNO;

        RemoveForbiddenEdgeMask( pBNS, &FixedStereoEdges, forbidden_stereo_edge_mask  );

        for ( i = 0; i < icr->num_sb_in2_only && cur_success < max_success; i ++ ) {
            j12 = icr->sb_in2_only[i];
            pv1 = pBNS->vert + (v1 = pStereoInChI->nBondAtom1[j12]-1);
            pv2 = pBNS->vert + (v2 = pStereoInChI->nBondAtom2[j12]-1);
            for ( k = 0; k < at2[v1].valence; k ++ ) {
                pe = pBNS->edge + (e = pv1->iedge[k]);
                if ( v2 == (pe->neighbor12 ^ v1) )
                    break; /* the edge has been found */
            }
            if ( k == at2[v1].valence ) {
                ret = RI_ERR_SYNTAX;
                goto exit_function;
            }
            /* check v1 */
            pv1->st_edge.cap --;
            pv1->st_edge.flow --;
            pv2->st_edge.flow --;
            pe->flow --; /* new radical on v2 */
            vRad = NO_VERTEX;
            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
            pv1->st_edge.cap ++;
            pv1->st_edge.flow ++;
            pv2->st_edge.flow ++;
            pe->flow ++; /* remove new radical on v2 */

            if ( ret == 1 /*&& !nDeltaH*/ && !nDeltaCharge && (v2 == vPathStart || v2 == vPathEnd) ) {
                vRad = (v2 == vPathStart)? vPathEnd : vPathStart;
            } else {
                pv2->st_edge.cap --;
                pv2->st_edge.flow --;
                pv1->st_edge.flow --;
                pe->flow --; /* new radical on v1 */
                vRad = NO_VERTEX;
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                pv2->st_edge.cap ++;
                pv2->st_edge.flow ++;
                pv1->st_edge.flow ++;
                pe->flow ++; /* remove new radical on v1 */
                if ( ret == 1 /*&& !nDeltaH*/ && !nDeltaCharge && (v1 == vPathStart || v1 == vPathEnd) ) {
                    vRad = (v1 == vPathStart)? vPathEnd : vPathStart;
                }
            }
            if ( vRad == NO_VERTEX ) {
                continue; /* radical did not affect this bond */
            }
            pvRad = pBNS->vert + vRad;
            /* detect =N-O*  */
            if ( pVA[vRad].cNumValenceElectrons == 6 && at2[vRad].valence == 1 &&
                 (peRad = pBNS->edge + pvRad->iedge[0])->flow == 0 &&
                 pVA[aN = peRad->neighbor12 ^ vRad].cNumValenceElectrons == 5 &&
                 at2[aN].valence == 2 ) {
                /*------------------------------------------------------------ 
                  Fix Metal disconnection/normalization inconsistency :
                                         disconnected  restored
                  R=N(+)-M     R=N--M     R=N  + M     R=N   + M                         
                    |       ->   ||    ->   ||     ->    |
                    O(-)         O          O            O* <- radical                        

                  The correct     R=N    + M(+)      
                  disconnection     |
                  would be this:    O(-)        
                --------------------------------------------------------------*/
                pvN  = pBNS->vert + aN;
                pvO  = pvRad;
                peNO = peRad;

                /* N-O*  => N=O */
                peNO->flow ++;
                pvO->st_edge.flow ++;
                pvN->st_edge.cap  ++;
                pvN->st_edge.flow ++;
                pBNS->tot_st_cap  += 1;
                pBNS->tot_st_flow += 2;
                cur_success ++;
            } else {
                /* all other radicals that affect stereo */
                delta = pvRad->st_edge.cap - pvRad->st_edge.flow;
                pvRad->st_edge.cap -= delta;
                pBNS->tot_st_cap   -= delta;
            }
        }
/*exit_case_01:*/
        SetForbiddenEdgeMask( pBNS, &FixedStereoEdges, forbidden_stereo_edge_mask  );
        if ( cur_success ) {
            tot_succes += cur_success;
            /* recalculate InChI from the structure */
            if ( 0 > (ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                            ppt_group_info, ppat_norm, ppat_prep ) ) ) {
                goto exit_function;
            }
            if ( ret = FillOutExtraFixedHDataRestr( pStruct ) ) {
                goto exit_function;
            }
            /*
            if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            */
            cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
            if ( cmpInChI & IDIF_PROBLEM ) {
                ret = RI_ERR_PROGR; /* severe restore problem */
                goto exit_function;
            }
            if ( err ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            cmpInChI2 = 0;
            memset ( icr2, 0, sizeof(*icr2) );
            if ( iRevrInChI || iOrigInChI ) {
                /* additional mobile-H compare in case of Fixed-H */
                cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                if ( cmpInChI & IDIF_PROBLEM ) {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if ( err ) {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
            }

            pStereoRevrs = (pStruct->pOneINChI[0]->StereoIsotopic &&
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters)?
                       pStruct->pOneINChI[0]->StereoIsotopic : pStruct->pOneINChI[0]->Stereo;

            
            pStereo2Revrs = (pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                       !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted)?
                                NULL:
                       (pStruct->pOneINChI[1]->StereoIsotopic &&
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters)?
                                pStruct->pOneINChI[1]->StereoIsotopic :
                                pStruct->pOneINChI[1]->Stereo;

        }
    }

    cur_success = 0;
    if ( !(cmpInChI & IDIF_SB_MISS) && (cmpInChI2 & IDIF_SB_MISS) &&
         icr2->num_sb_in2_only &&
         0 < (max_success = pBNS->tot_st_cap - pBNS->tot_st_flow) ) {
        /*----------------------------------------------------*/
        /* case 02: missing stereogenic bond in Mobile-H only */
        /* X=N-O*  => X=N=O and eliminate radical             */
        /*----------------------------------------------------*/
        int retC, ret2C, retS, ret2S;
        INCHI_MODE cmpInChI_Prev, cmpInChI2_Prev;
        ICR  icr_Prev, icr2_Prev;

        /* blind attepmt */
        icr_Prev       = *icr;
        icr2_Prev      = *icr2;
        cmpInChI_Prev  = cmpInChI;
        cmpInChI2_Prev = cmpInChI2;
        for ( i = AllRadList.num_edges = 0; i < pStruct->num_atoms; i ++ ) {
            if ( pBNS->vert[i].st_edge.cap - pBNS->vert[i].st_edge.flow == 1 &&
                 (ret = AddToEdgeList( &AllRadList, i, INC_ADD_EDGE ) ) ) {
                goto exit_function;
            }
        }
        for ( i = 0; i < AllRadList.num_edges; i ++ ) {
            j = AllRadList.pnEdges[i];
            pBNS->vert[j].st_edge.cap -= 1;
            pBNS->tot_st_cap -= 1;
        }
        /*-------------------------------------------------*/
        /* re-create InChI and see whether it looks better */
        /*-------------------------------------------------*/
        if ( 0 > (ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                        ppt_group_info, ppat_norm, ppat_prep ) ) ) {
            goto exit_function;
        }
        if ( ret = FillOutExtraFixedHDataRestr( pStruct ) ) {
            goto exit_function;
        }
        /*
        if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
            goto exit_function;
        }
        */
        cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
        if ( cmpInChI & IDIF_PROBLEM ) {
            ret = RI_ERR_PROGR; /* severe restore problem */
            goto exit_function;
        }
        if ( err ) {
            ret = RI_ERR_ALLOC;
            goto exit_function;
        }
        cmpInChI2 = 0;
        memset ( icr2, 0, sizeof(*icr2) );
        if ( iRevrInChI || iOrigInChI ) {
            /* additional mobile-H compare in case of Fixed-H */
            cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
            if ( cmpInChI & IDIF_PROBLEM ) {
                ret = RI_ERR_PROGR; /* severe restore problem */
                goto exit_function;
            }
            if ( err ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
        }
        retC  = CompareIcr( icr, &icr_Prev, NULL, NULL, IDIFF_CONSTIT );
        retS  = CompareIcr( icr, &icr_Prev, NULL, NULL, IDIFF_STEREO );
        ret2C = CompareIcr( icr2, &icr2_Prev, NULL, NULL, IDIFF_CONSTIT );
        ret2S = CompareIcr( icr2, &icr2_Prev, NULL, NULL, IDIFF_STEREO );

        if ( 0 >= retC &&
             0 >= retS &&
             0 >= ret2C &&
             0 > ret2S ) {
            ; /* accept */
        } else {
            /* reject */
            for ( i = 0; i < AllRadList.num_edges; i ++ ) {
                j = AllRadList.pnEdges[i];
                pBNS->vert[j].st_edge.cap += 1;
                pBNS->tot_st_cap += 1;
            }

            /*-------------------------------------------------*/
            /* re-create InChI-- return to previous state      */
            /*-------------------------------------------------*/
            if ( 0 > (ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                            ppt_group_info, ppat_norm, ppat_prep ) ) ) {
                goto exit_function;
            }
            if ( ret = FillOutExtraFixedHDataRestr( pStruct ) ) {
                goto exit_function;
            }
            /*
            if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            */
            cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
            if ( cmpInChI & IDIF_PROBLEM ) {
                ret = RI_ERR_PROGR; /* severe restore problem */
                goto exit_function;
            }
            if ( err ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            cmpInChI2 = 0;
            memset ( icr2, 0, sizeof(*icr2) );
            if ( iRevrInChI || iOrigInChI ) {
                /* additional mobile-H compare in case of Fixed-H */
                cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                if ( cmpInChI & IDIF_PROBLEM ) {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if ( err ) {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
            }
            pStereoRevrs = (pStruct->pOneINChI[0]->StereoIsotopic &&
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters)?
                       pStruct->pOneINChI[0]->StereoIsotopic : pStruct->pOneINChI[0]->Stereo;

            
            pStereo2Revrs = (pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                       !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted)?
                                NULL:
                       (pStruct->pOneINChI[1]->StereoIsotopic &&
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters)?
                                pStruct->pOneINChI[1]->StereoIsotopic :
                                pStruct->pOneINChI[1]->Stereo;
        }
/*exit_case_02:;*/
    }

    cur_success = 0;
    if ( pStruct->bMobileH == TAUT_NON && (cmpInChI & IDIF_SB_EXTRA_UNDF) &&
         pStruct->endpoint ) {
        /*------------------------------------------------------*/
        /* case 03: extra stereogenic bond in Fixed-H  only     */
        /* in Mobile-H this bond is not stereogenic.            */
        /* Since this bond parity is not known, it is UNDEFINED */
        /*------------------------------------------------------*/
        int bDone, num_endpoints;

        TautMinusEdges[0].num_edges = 0;
        TautMinusEdges[1].num_edges = 0;
        AllChargeEdges.num_edges = 0;
        /* in1 => in restored structure; in2 => in original InChI */
        for ( i = 0; i < icr->num_sb_undef_in1_only; i ++ ) {
            j12 = icr->sb_undef_in1_only[i];
            pv1 = pBNS->vert + (v1 = pStereoRevrs->nBondAtom1[j12]-1);
            pv2 = pBNS->vert + (v2 = pStereoRevrs->nBondAtom2[j12]-1);

            if ( pStereo2Revrs ) {
                /* reject if it is extra in Mobile-H also */
                if ( icr2->num_sb_undef_in1_only ) {
                    for ( j = 0; j < icr2->num_sb_undef_in1_only; j ++ ) {
                        k = icr2->sb_undef_in1_only[j];
                        if ( v1 == pStereo2Revrs->nBondAtom1[k] &&
                             v2 == pStereo2Revrs->nBondAtom2[k] ) {
                            break;
                        }
                    }
                    if ( j < icr->num_sb_in1_only ) {
                        continue; /* extra stereobond in Mobile H also */
                    }
                }
            }
            /* reject if it is a stereobond in Mobile-H also */
            if ( pStereo2InChI && pStereo2InChI->nNumberOfStereoBonds ) {
                for ( j = 0; j < pStereo2InChI->nNumberOfStereoBonds; j ++ ) {
                    if ( v1 == pStereo2InChI->nBondAtom1[j] &&
                         v2 == pStereo2InChI->nBondAtom1[j]  ) {
                        break;
                    }
                }
                if ( j < pStereo2InChI->nNumberOfStereoBonds ) {
                    continue; /* ignore this extra stereo bond: it is in Mobile-H */
                }
            }
            /* find the edge between v1 and v2 */
            for ( k = 0; k < at2[v1].valence; k ++ ) {
                pe = pBNS->edge + (e = pv1->iedge[k]);
                if ( v2 == (pe->neighbor12 ^ v1) )
                    break; /* the edge has been found */
            }
            if ( k == at2[v1].valence ) {
                ret = RI_ERR_SYNTAX;
                goto exit_function;
            }
            /* Fix all charges except negative charges on tautomeric endpoints */
            if ( !AllChargeEdges.num_edges && !TautMinusEdges[0].num_edges && !TautMinusEdges[1].num_edges ) {
                for ( j = 0; j < pStruct->num_atoms; j ++ ) {
                    if ( (k=pVA[j].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                        if ( !pStruct->endpoint[j] ) {
                            if (ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                                goto exit_function;
                            }
                        } else
                        if ( pVA[j].cNumValenceElectrons == 6 ) {
                            /* O */
                            if (ret = AddToEdgeList( TautMinusEdges+0, k, INC_ADD_EDGE ) ) {
                                goto exit_function;
                            }
                        } else {
                            /* N */
                            if (ret = AddToEdgeList( TautMinusEdges+1, k, INC_ADD_EDGE ) ) {
                                goto exit_function;
                            }
                        }
                    }
                    if ( (k=pVA[j].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                        if ( ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                        /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                        if ( pVA[j].cNumValenceElectrons == 5 && !pVA[j].cMetal && /* N, P, As */
                             NO_VERTEX != (k = GetChargeFlowerUpperEdge( pBNS, pVA, k ))) {

                            if ( !pBNS->edge[j].forbidden && pBNS->edge[k].flow ) {
                                if ( ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                                    goto exit_function;
                                }
                            }
                        }
                    }
                }
            }
            if ( !pe->flow )
                continue;
            /* fix all charges except tautomeric; first allow only O, then only N, finally both N and O */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
            for ( k = 1, bDone = 0; k < 4 && !bDone; k ++ ) {
                /* fix tautomeric charges */
                num_endpoints = (TautMinusEdges+0)->num_edges + (TautMinusEdges+1)->num_edges;
                if ( k == 2 ) {
                    /* fix charges on O */
                    SetForbiddenEdgeMask( pBNS, TautMinusEdges+0, forbidden_edge_mask  );
                    num_endpoints -= (TautMinusEdges+0)->num_edges;
                }
                if ( k == 1 ) {
                    SetForbiddenEdgeMask( pBNS, TautMinusEdges+1, forbidden_edge_mask  );
                    num_endpoints -= (TautMinusEdges+1)->num_edges;
                }
                if ( num_endpoints >= 2 ) {
                    delta = 1;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask; /* fix stereobond */
                    pe->flow -= delta;                    /* decrement stereobond order */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2*delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 0 ) {
                        /* Negative charge has been moved, no change in number of charges */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if ( ret > 0 ) {
                            (*pnNumRunBNS) ++;
                            cur_success ++; /* 01 */
                            bDone = 1;
                        }
                    } else {
                        pe->forbidden &= ~forbidden_edge_mask;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2*delta;
                    }
                }
                /* unfix tautomeric charges */
                if ( k == 2 )
                    RemoveForbiddenEdgeMask( pBNS, TautMinusEdges+0, forbidden_edge_mask  );
                if ( k == 1 )
                    RemoveForbiddenEdgeMask( pBNS, TautMinusEdges+1, forbidden_edge_mask  );
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
/*exit_case_03:*/
        if ( cur_success ) {
            tot_succes += cur_success;
            /* recalculate InChI from the structure */
            if ( 0 > (ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                            ppt_group_info, ppat_norm, ppat_prep ) ) ) {
                goto exit_function;
            }
            if ( ret = FillOutExtraFixedHDataRestr( pStruct ) ) {
                goto exit_function;
            }
            /*
            if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            */
            cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
            if ( cmpInChI & IDIF_PROBLEM ) {
                ret = RI_ERR_PROGR; /* severe restore problem */
                goto exit_function;
            }
            if ( err ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            cmpInChI2 = 0;
            memset ( icr2, 0, sizeof(*icr2) );
            if ( iRevrInChI || iOrigInChI ) {
                /* additional mobile-H compare in case of Fixed-H */
                cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                if ( cmpInChI & IDIF_PROBLEM ) {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if ( err ) {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
            }
            pStereoRevrs = (pStruct->pOneINChI[0]->StereoIsotopic &&
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters)?
                       pStruct->pOneINChI[0]->StereoIsotopic : pStruct->pOneINChI[0]->Stereo;

            
            pStereo2Revrs = (pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                       !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted)?
                                NULL:
                       (pStruct->pOneINChI[1]->StereoIsotopic &&
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters)?
                                pStruct->pOneINChI[1]->StereoIsotopic :
                                pStruct->pOneINChI[1]->Stereo;

        }
    }

    cur_success = 0;
    if ( (cmpInChI & IDIF_SB_EXTRA_UNDF) ) {
        /*------------------------------------------------------*/
        /* case 04: extra stereogenic bond                      */
        /* Since this bond parity is not known, it is UNDEFINED */
        /*------------------------------------------------------*/
        int bDone, num_endpoints;

        TautMinusEdges[0].num_edges = 0;
        TautMinusEdges[1].num_edges = 0;
        AllChargeEdges.num_edges = 0;
        /* in1 => in restored structure; in2 => in original InChI */
        for ( i = 0; i < icr->num_sb_undef_in1_only; i ++ ) {
            j12 = icr->sb_undef_in1_only[i];
            pv1 = pBNS->vert + (v1 = pStereoRevrs->nBondAtom1[j12]-1);
            pv2 = pBNS->vert + (v2 = pStereoRevrs->nBondAtom2[j12]-1);

            /* find the edge between v1 and v2 */
            for ( k = 0; k < at2[v1].valence; k ++ ) {
                pe = pBNS->edge + (e = pv1->iedge[k]);
                if ( v2 == (pe->neighbor12 ^ v1) )
                    break; /* the edge has been found */
            }
            if ( k == at2[v1].valence ) {
                ret = RI_ERR_SYNTAX;
                goto exit_function;
            }
            if ( pStereo2Revrs ) {
                /* reject if it is not extra in Mobile-H also */
                if ( icr2->num_sb_undef_in1_only ) {
                    for ( j = 0; j < icr2->num_sb_undef_in1_only; j ++ ) {
                        k = icr2->sb_undef_in1_only[j];
                        if ( v1 == pStereo2Revrs->nBondAtom1[k] &&
                             v2 == pStereo2Revrs->nBondAtom2[k] ) {
                            break;
                        }
                    }
                    if ( j == icr->num_sb_in1_only ) {
                        continue; /* extra stereobond only in Fixed-H, not in Mobile H also */
                    }
                }
            }
            
            /* Fix all charges except negative charges on tautomeric endpoints */
            if ( !AllChargeEdges.num_edges && !TautMinusEdges[0].num_edges && !TautMinusEdges[1].num_edges ) {
                for ( j = 0; j < pStruct->num_atoms; j ++ ) {
                    if ( (k=pVA[j].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                        if (ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                    }
                    if ( (k=pVA[j].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                        int bMayBeUnfixed = !at2[j].num_H && !(pStruct->endpoint && pStruct->endpoint[j]);
                        if ( bMayBeUnfixed && pVA[j].cNumValenceElectrons == 6 || 
                             pVA[j].cNumValenceElectrons == 5 && pVA[j].cPeriodicRowNumber > 1 ) {
                            /* O & P */
                            if (ret = AddToEdgeList( TautMinusEdges+0, k, INC_ADD_EDGE ) ) {
                                goto exit_function;
                            }
                        } else
                        if ( bMayBeUnfixed &&
                             pVA[j].cNumValenceElectrons == 5 && pVA[j].cPeriodicRowNumber == 1 ) {
                            /* N */
                            if (ret = AddToEdgeList( TautMinusEdges+1, k, INC_ADD_EDGE ) ) {
                                goto exit_function;
                            }
                        } else {
                            if ( ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                                goto exit_function;
                            }
                        }
                        /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                        if ( pVA[j].cNumValenceElectrons == 5 && !pVA[j].cMetal && /* N, P, As */
                             NO_VERTEX != (k = GetChargeFlowerUpperEdge( pBNS, pVA, k ))) {
                            if ( !pBNS->edge[j].forbidden && pBNS->edge[k].flow ) {
                                if ( ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                                    goto exit_function;
                                }
                            }
                        }
                    }
                }
            }
            if ( !pe->flow )
                continue;
            /* fix all charges except tautomeric; first allow only O, then only N, finally both N and O */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
            for ( k = 1, bDone = 0; k < 4 && !bDone; k ++ ) {
                /* fix positive charges on heteroatoms */
                num_endpoints = (TautMinusEdges+0)->num_edges + (TautMinusEdges+1)->num_edges;
                if ( k == 2 ) {
                    /* fix charges on O */
                    SetForbiddenEdgeMask( pBNS, TautMinusEdges+0, forbidden_edge_mask  );
                    num_endpoints -= (TautMinusEdges+0)->num_edges;
                }
                if ( k == 1 ) {
                    /* fix charges on N */
                    SetForbiddenEdgeMask( pBNS, TautMinusEdges+1, forbidden_edge_mask  );
                    num_endpoints -= (TautMinusEdges+1)->num_edges;
                }
                if ( num_endpoints >= 2 ) {
                    delta = 1;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask; /* fix stereobond */
                    pe->flow -= delta;                    /* decrement stereobond order */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2*delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 0 ) {
                        /* Negative charge has been moved, no change in number of charges */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if ( ret > 0 ) {
                            (*pnNumRunBNS) ++;
                            cur_success ++; /* 01 */
                            bDone = 1;
                        }
                    } else {
                        pe->forbidden &= ~forbidden_edge_mask;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2*delta;
                    }
                }
                /* unfix tautomeric charges */
                if ( k == 2 )
                    RemoveForbiddenEdgeMask( pBNS, TautMinusEdges+0, forbidden_edge_mask  );
                if ( k == 1 )
                    RemoveForbiddenEdgeMask( pBNS, TautMinusEdges+1, forbidden_edge_mask  );
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
/*exit_case_04:*/
        if ( cur_success ) {
            tot_succes += cur_success;
            /* recalculate InChI from the structure */
            if ( 0 > (ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                            ppt_group_info, ppat_norm, ppat_prep ) ) ) {
                goto exit_function;
            }
            if ( ret = FillOutExtraFixedHDataRestr( pStruct ) ) {
                goto exit_function;
            }
            /*
            if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            */
            cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
            if ( cmpInChI & IDIF_PROBLEM ) {
                ret = RI_ERR_PROGR; /* severe restore problem */
                goto exit_function;
            }
            if ( err ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            cmpInChI2 = 0;
            memset ( icr2, 0, sizeof(*icr2) );
            if ( iRevrInChI || iOrigInChI ) {
                /* additional mobile-H compare in case of Fixed-H */
                cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                if ( cmpInChI & IDIF_PROBLEM ) {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if ( err ) {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
            }
            pStereoRevrs = (pStruct->pOneINChI[0]->StereoIsotopic &&
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters)?
                       pStruct->pOneINChI[0]->StereoIsotopic : pStruct->pOneINChI[0]->Stereo;

            
            pStereo2Revrs = (pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                       !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted)?
                                NULL:
                       (pStruct->pOneINChI[1]->StereoIsotopic &&
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters)?
                                pStruct->pOneINChI[1]->StereoIsotopic :
                                pStruct->pOneINChI[1]->Stereo;

        }
    }

    cur_success = 0;
    if ( pStruct->bMobileH == TAUT_YES &&
         (cmpInChI & IDIF_SB_EXTRA_UNDF &&
         !pStruct->ti.num_t_groups)
         /*pStruct->bMobileH == TAUT_NON && (cmpInChI2 & IDIF_SB_EXTRA_UNDF)*/) {
        /*----------------------------------------------------------*/
        /* case 05: extra stereogenic bond on =NH2(+), (B, Mobile-H)*/
        /*                   H             H                        */
        /* original: N(+)=-N<   ->  N--==N/                         */
        /* (A)               H                                      */
        /*                             double bond is marked as     */
        /*                             not stereogenic due to       */
        /*                             its change during proton     */
        /*                             removal => No Stereo bond    */
        /*                             (=NH may be tautomeric)      */
        /*                                                          */
        /*                   H             H                        */
        /* original: N=-N(+)<   ->  N--==N/                         */
        /* (B)               H                                      */
        /*                             double bond was not          */
        /*                             changed during proton        */
        /* In Fixed-H this bond        removal => Undef Stereo      */
        /* may not be stereogenic      (=NH is not tautomeric)      */
        /* (a) due to (+) movement                                  */
        /* (b) due to symmetry (2H), even if isotopic               */
        /*                                                          */
        /* Fixed-H: move (+) to or from NH2 for Undef or No stereo  */
        /*          respectively                                    */
        /* Mobile-H: Add H(+) to =NH and move the charge to =N-     */
        /*           to eliminate Undef stereo                      */
        /*          Move charge from N to -NH2 to create            */
        /*           Undef Stereo                                   */
        /* Since this bond parity is not known, it is UNDEFINED     */
        /*                                                          */
        /* Solution: Add H(+) to =NH and move charge to -N=         */
        /*                                                          */
        /*----------------------------------------------------------*/
        int aN, aC, i1, i2, vPlusMinus;
        AllChargeEdges.num_edges = 0;
        /* in1 => in restored structure; in2 => in original InChI */
        for ( i = 0; i < icr->num_sb_undef_in1_only; i ++ ) {
            j12 = icr->sb_undef_in1_only[i];
            pv1 = pBNS->vert + (v1 = pStereoRevrs->nBondAtom1[j12]-1);
            pv2 = pBNS->vert + (v2 = pStereoRevrs->nBondAtom2[j12]-1);
            /* indicators of -NH: */
            i1 = at2[v1].valence == 1 && at2[v1].num_H == 1 && !at2[v1].endpoint &&
                 pVA[v1].cNumValenceElectrons == 5 && pVA[v1].cPeriodicRowNumber == 1;
            i2 = at2[v2].valence == 1 && at2[v2].num_H == 1 && !at2[v2].endpoint &&
                 pVA[v2].cNumValenceElectrons == 5 && pVA[v2].cPeriodicRowNumber == 1;
            if ( !i1 && !i2  || i1 && i2 ) {
                continue;
            }
            /* find the edge between v1 and v2 */
            for ( k = 0; k < at2[v1].valence; k ++ ) {
                pe = pBNS->edge + (e = pv1->iedge[k]);
                if ( v2 == (pe->neighbor12 ^ v1) )
                    break; /* the edge has been found */
            }
            if ( k == at2[v1].valence ) {
                ret = RI_ERR_SYNTAX;
                goto exit_function;
            }
            if ( pe->flow != 1 ) {
                continue; /* already charged */
            }
            aN = i1? v1 : v2; /* -NH atom */
            aC = i1? v2 : v1; /* neighbor */
            /* Replace =NH with -NH2
               Create such a charge on some -N< that may be moved to NH2 to remove H(+):
               transformation:
               from:  HN=C-=-N=(+vert)-Y=(+super)-(+/-)
               to:   2HN-C*-=-N=(+vert)-Y=(+super)-(+/-)*
               Run BNS to obtain:
                     2HN-C=-=N(+)-(+vert)=Y-(+super)=(+/-)
            */
            vPlusMinus = GetPlusMinusVertex( pBNS, pTCGroups, 1, 0 );
            if ( NO_VERTEX == vPlusMinus ) {
                break; /* cannot do anything */
            }
            /* increase edges to -Y-(+/-)-Y- capacities */
            delta = 1;
            for ( i1 = 0; i1 < pBNS->vert[vPlusMinus].num_adj_edges; i1 ++ ) {
                i2 = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i1]].neighbor12 ^ vPlusMinus;
                for ( k = 0; k < pBNS->vert[i2].num_adj_edges; k ++ ) {
                    pBNS->edge[pBNS->vert[i2].iedge[k]].cap += delta;
                }
            }
            /* Fix all charges except (+) on -N< */
            if ( !AllChargeEdges.num_edges ) {
                for ( j = 0; j < pStruct->num_atoms; j ++ ) {
                    if ( (k=pVA[j].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                        if (ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                    }
                    if ( (k=pVA[j].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                        if ( pVA[j].cNumValenceElectrons == 5 && pVA[j].cPeriodicRowNumber == 1 &&
                             !at2[j].num_H && at2[j].valence == 3 &&
                             !(at2[j].endpoint || pStruct->endpoint && pStruct->endpoint[j]) ) {
                                 ; /* do not fix -N< or =N(+)< */
                        } else {
                            /* all others */
                            if (ret = AddToEdgeList( TautMinusEdges+0, k, INC_ADD_EDGE ) ) {
                                goto exit_function;
                            }
                        } 
                        /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                        if ( pVA[j].cNumValenceElectrons == 5 && !pVA[j].cMetal && /* N, P, As */
                             NO_VERTEX != (k = GetChargeFlowerUpperEdge( pBNS, pVA, k ))) {
                            if ( !pBNS->edge[j].forbidden && pBNS->edge[k].flow ) {
                                if ( ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                                    goto exit_function;
                                }
                            }
                        }
                    }
                }
            }
            /* Make bond to =NH single, add radical to aC */
            pe->flow                    -= delta; /* make single bond */
            pBNS->vert[aN].st_edge.flow -= delta;
            pBNS->vert[aN].st_edge.cap  -= delta; /* avoid radical on N */
            pBNS->vert[aC].st_edge.flow -= delta; /* create radical on C */
            pBNS->vert[vPlusMinus].st_edge.cap += delta; /* create radical on (+/-) */
            pBNS->tot_st_flow           -=  2*delta;
            /* fix C-NH bond */
            if ( ret = AddToEdgeList( &AllChargeEdges, e, INC_ADD_EDGE ) ) {
                goto exit_function;
            }
            /* pBNS->tot_st_cap is unchanged */
            /* find all aC edges except pe to fix them */
            /* 2. Check whether it would work and do if it would */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );/* fix aC edges */
            pe->cap ++;
            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && (vPathEnd == vPlusMinus && vPathStart == aC ||
                              vPathEnd == aC && vPathStart == vPlusMinus) && nDeltaCharge == 1 ) {
                /* Negative charge has been moved, no change in number of charges */
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                if ( ret > 0 ) {
                    (*pnNumRunBNS) ++;
                    /* 3. Add  H to -NH and register increaded charge */
                    pStruct->at[aN].num_H ++;
                    pTCGroups->total_charge ++;
                    cur_success ++; /* 01 */
                }
            } else {
                pe->flow                    += delta; /* make single bond */
                pBNS->vert[aN].st_edge.flow += delta;
                pBNS->vert[aN].st_edge.cap  += delta; /* avoid radical on N */
                pBNS->vert[aC].st_edge.flow += delta; /* create radical on C */
                pBNS->vert[vPlusMinus].st_edge.cap -= delta; /* create radical on (+/-) */
                pBNS->tot_st_flow           +=  2*delta;
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );/* fix aC edges */
                AllChargeEdges.num_edges --; /* remove pe from the list */
                CurrEdges.num_edges = 0;
                continue; /* should not happen */
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );/* fix aC edges */
            AllChargeEdges.num_edges --; /* remove pe from the list */
            CurrEdges.num_edges = 0;
        }
/*exit_case_05:*/
        if ( cur_success ) {
            tot_succes += cur_success;
            /* recalculate InChI from the structure */
            if ( 0 > (ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                            ppt_group_info, ppat_norm, ppat_prep ) ) ) {
                goto exit_function;
            }
            if ( ret = FillOutExtraFixedHDataRestr( pStruct ) ) {
                goto exit_function;
            }
            /*
            if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            */
            cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr, &err );
            if ( cmpInChI & IDIF_PROBLEM ) {
                ret = RI_ERR_PROGR; /* severe restore problem */
                goto exit_function;
            }
            if ( err ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            cmpInChI2 = 0;
            memset ( icr2, 0, sizeof(*icr2) );
            if ( iRevrInChI || iOrigInChI ) {
                /* additional mobile-H compare in case of Fixed-H */
                cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                if ( cmpInChI & IDIF_PROBLEM ) {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if ( err ) {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
            }
            pStereoRevrs = (pStruct->pOneINChI[0]->StereoIsotopic &&
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters)?
                       pStruct->pOneINChI[0]->StereoIsotopic : pStruct->pOneINChI[0]->Stereo;

            
            pStereo2Revrs = (pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                       !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted)?
                                NULL:
                       (pStruct->pOneINChI[1]->StereoIsotopic &&
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters)?
                                pStruct->pOneINChI[1]->StereoIsotopic :
                                pStruct->pOneINChI[1]->Stereo;

        }
    }

    cur_success = 0;
    if ( pStruct->bMobileH == TAUT_NON && pStereo2Revrs /* added check 2006-04-05 */ &&
         (cmpInChI2 & IDIF_SB_EXTRA_UNDF &&
         !pStruct->ti.num_t_groups)
         /*pStruct->bMobileH == TAUT_NON && (cmpInChI2 & IDIF_SB_EXTRA_UNDF)*/) {
        /*----------------------------------------------------------*/
        /* case 06: extra stereogenic bond on =NH2(+), (B, Fixed-H) */
        /*                   H                H        ===========  */
        /* original: N(+)=-N<   ->  N--==N(+)<                      */
        /* (A)               H                H                     */
        /*                             double bond in Mobile-H      */
        /*                             layer has Undef stereo       */
        /*                                                          */
        /*                                                          */
        /* Fixed-H: move (+) to or from NH2 for Undef or No stereo  */
        /*          respectively                                    */
        /* Mobile-H: Add H(+) to =NH and move the charge to =N-     */
        /*           to eliminate Undef stereo                      */
        /*          Move charge from N to -NH2 to create            */
        /*           Undef Stereo                                   */
        /* Since this bond parity is not known, it is UNDEFINED     */
        /*                                                          */
        /* Solution: Move (+) from -NH2(+) to othe -N<              */
        /*                                                          */
        /*----------------------------------------------------------*/
        int aN, aC, i1, i2, ePlus;
        BNS_EDGE   *pePlus;
        AllChargeEdges.num_edges = 0;
        /* in1 => in restored structure; in2 => in original InChI */
        for ( i = 0; i < icr2->num_sb_undef_in1_only; i ++ ) {
            j12 = icr2->sb_undef_in1_only[i];
            pv1 = pBNS->vert + (v1 = pStereo2Revrs->nBondAtom1[j12]-1);
            pv2 = pBNS->vert + (v2 = pStereo2Revrs->nBondAtom2[j12]-1);
            /* indicators of -NH: */
            i1 = at2[v1].valence == 1 && at2[v1].num_H == 2 && !at2[v1].endpoint &&
                 pVA[v1].cNumValenceElectrons == 5 && pVA[v1].cPeriodicRowNumber == 1;
            i2 = at2[v2].valence == 1 && at2[v2].num_H == 2 && !at2[v2].endpoint &&
                 pVA[v2].cNumValenceElectrons == 5 && pVA[v2].cPeriodicRowNumber == 1;
            if ( !i1 && !i2  || i1 && i2 ) {
                continue;
            }
            /* find the edge between v1 and v2 */
            for ( k = 0; k < at2[v1].valence; k ++ ) {
                pe = pBNS->edge + (e = pv1->iedge[k]);
                if ( v2 == (pe->neighbor12 ^ v1) )
                    break; /* the edge has been found */
            }
            if ( k == at2[v1].valence ) {
                ret = RI_ERR_SYNTAX;
                goto exit_function;
            }
            if ( pe->flow != 1 ) {
                continue; /* already charged */
            }
            aN = i1? v1 : v2; /* -NH atom */
            aC = i1? v2 : v1; /* neighbor */
            if ( 0 > (ePlus = pVA[aN].nCPlusGroupEdge-1) ||
                 (pePlus = pBNS->edge + ePlus)->flow ||  /* must be (+) charged */
                 pePlus->forbidden ) {
                continue;
            }
            /* Move (+) from =NH2(+) to some other -N<
            */
            /* Fix all charges except (+) on -N< */
            if ( !AllChargeEdges.num_edges ) {
                for ( j = 0; j < pStruct->num_atoms; j ++ ) {
                    if ( (k=pVA[j].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                        if (ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                    }
                    if ( (k=pVA[j].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[k].forbidden ) {
                        if ( pVA[j].cNumValenceElectrons == 5 && pVA[j].cPeriodicRowNumber == 1 &&
                             !at2[j].num_H && at2[j].valence == 3 &&
                             !(at2[j].endpoint || pStruct->endpoint && pStruct->endpoint[j]) ) {
                                 ; /* do not fix -N< or =N(+)< */
                        } else {
                            /* all others */
                            if (ret = AddToEdgeList( TautMinusEdges+0, k, INC_ADD_EDGE ) ) {
                                goto exit_function;
                            }
                        } 
                        /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
                        if ( pVA[j].cNumValenceElectrons == 5 && !pVA[j].cMetal && /* N, P, As */
                             NO_VERTEX != (k = GetChargeFlowerUpperEdge( pBNS, pVA, k ))) {
                            if ( !pBNS->edge[j].forbidden && pBNS->edge[k].flow ) {
                                if ( ret = AddToEdgeList( &AllChargeEdges, k, INC_ADD_EDGE ) ) {
                                    goto exit_function;
                                }
                            }
                        }
                    }
                }
            }
            /* pePlus edge is already fixed; unfix it */
            /* To decrement (+) on =NH2(+) decrement its double bond order */
            delta = 1;
            if ( !pe->flow )
                continue;
            pv1 = pBNS->vert + (v1 = pe->neighbor1);
            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

            delta = 1;
            pe->flow -= delta;
            pv1->st_edge.flow -= delta;
            pv2->st_edge.flow -= delta;
            pBNS->tot_st_flow -= 2*delta;

            pe->forbidden     |= forbidden_edge_mask;
            pePlus->forbidden &= ~forbidden_edge_mask;

            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                              vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 0 ) {
                /* (+)charge was just moved, no change in number of charges */
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                if ( ret > 0 ) {
                    (*pnNumRunBNS) ++;
                    cur_success ++; /* 01 */
                }
            } else {
                pe->flow += delta; /* roll back */
                pv1->st_edge.flow += delta;
                pv2->st_edge.flow += delta;
                pBNS->tot_st_flow += 2*delta;
            }
            pe->forbidden     &= ~forbidden_edge_mask;
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );/* fix aC edges */
        }
/*exit_case_06:*/
        if ( cur_success ) {
            tot_succes += cur_success;
            /* recalculate InChI from the structure */
            if ( 0 > (ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                            ppt_group_info, ppat_norm, ppat_prep ) ) ) {
                goto exit_function;
            }
            if ( ret = FillOutExtraFixedHDataRestr( pStruct ) ) {
                goto exit_function;
            }
            /*
            if ( ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            */
            cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *v2*/, icr2, &err );
            if ( cmpInChI & IDIF_PROBLEM ) {
                ret = RI_ERR_PROGR; /* severe restore problem */
                goto exit_function;
            }
            if ( err ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            cmpInChI2 = 0;
            memset ( icr2, 0, sizeof(*icr2) );
            if ( iRevrInChI || iOrigInChI ) {
                /* additional mobile-H compare in case of Fixed-H */
                cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *v2*/, icr2, &err );
                if ( cmpInChI & IDIF_PROBLEM ) {
                    ret = RI_ERR_PROGR; /* severe restore problem */
                    goto exit_function;
                }
                if ( err ) {
                    ret = RI_ERR_ALLOC;
                    goto exit_function;
                }
            }
            pStereoRevrs = (pStruct->pOneINChI[0]->StereoIsotopic &&
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[0]->StereoIsotopic->nNumberOfStereoCenters)?
                       pStruct->pOneINChI[0]->StereoIsotopic : pStruct->pOneINChI[0]->Stereo;

            
            pStereo2Revrs = (pStruct->bMobileH == TAUT_YES || !pStruct->pOneINChI[1] ||
                       !pStruct->pOneINChI[1]->nNumberOfAtoms || pStruct->pOneINChI[1]->bDeleted)?
                                NULL:
                       (pStruct->pOneINChI[1]->StereoIsotopic &&
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoBonds +
                       pStruct->pOneINChI[1]->StereoIsotopic->nNumberOfStereoCenters)?
                                pStruct->pOneINChI[1]->StereoIsotopic :
                                pStruct->pOneINChI[1]->Stereo;

        }
    }


exit_function:
    SetForbiddenEdgeMask( pBNS, &FixedStereoEdges, forbidden_stereo_edge_mask  );
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &CurrEdges, EDGE_LIST_FREE );
    AllocEdgeList( &NFlowerEdges, EDGE_LIST_FREE );
    AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_FREE );
    AllocEdgeList( &FixedStereoEdges, EDGE_LIST_FREE );
    AllocEdgeList( &AllRadList, EDGE_LIST_FREE ); /* eliminate memory leak */
    AllocEdgeList( TautMinusEdges+0, EDGE_LIST_FREE );
    AllocEdgeList( TautMinusEdges+1, EDGE_LIST_FREE );

    return ret;
}
#endif
