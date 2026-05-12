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

/*#define CHECK_WIN32_VC_HEAP*/

#include "mode.h"

#if ( READ_INCHI_STRING == 1 )

#include "ichitime.h"
#include "ichicant.h"
#include "ichirvrs.h"

#include "bcf_s.h"

#define INC_ADD_EDGE 64

/****************************************************************************/
int GetPlusMinusVertex( BN_STRUCT *pBNS,
                        ALL_TC_GROUPS *pTCGroups,
                        int bCheckForbiddenPlus,
                        int bCheckForbiddenMinus )
{
    int k, ePlusSuper, eMinusSuper, vPlusSuper, vPlusMinus1 = NO_VERTEX, vPlusMinus2 = NO_VERTEX; /* djb-rwth: removing redundant variables */
    BNS_EDGE *pEdge;
    if (( k = pTCGroups->nGroup[TCG_Plus] ) >= 0 &&
        ( ePlusSuper = pTCGroups->pTCG[k].nForwardEdge ) > 0 &&
        ( vPlusSuper = pTCGroups->pTCG[k].nVertexNumber ) >= pBNS->num_atoms &&
         !( ( pEdge = pBNS->edge + ePlusSuper )->forbidden && bCheckForbiddenPlus ))
    {

        vPlusMinus1 = pEdge->neighbor12 ^ vPlusSuper;
    }
    if (( k = pTCGroups->nGroup[TCG_Minus] ) >= 0 &&
        ( eMinusSuper = pTCGroups->pTCG[k].nForwardEdge ) > 0 &&
        ( pTCGroups->pTCG[k].nVertexNumber ) >= pBNS->num_atoms && /* djb-rwth: removing redundant code */
         !( ( pEdge = pBNS->edge + eMinusSuper )->forbidden && bCheckForbiddenMinus ))
    {

        vPlusMinus2 = pEdge->neighbor12 ^ eMinusSuper;
    }
    if ((bCheckForbiddenPlus && NO_VERTEX == vPlusMinus1) ||
         (bCheckForbiddenMinus && NO_VERTEX == vPlusMinus2)) /* djb-rwth: addressing LLVM warnings */
    {
        return NO_VERTEX;
    }

    return ( NO_VERTEX != vPlusMinus1 ) ? vPlusMinus1 : vPlusMinus2;
}


/****************************************************************************/
int bIsUnsatCarbonInASmallRing( inp_ATOM *at2,
                                VAL_AT *pVA,
                                int iat,
                                BFS_Q *pbfsq,
                                int min_ring_size )
{
    int j, nCurRingSize, nMinRingSize;
    if (min_ring_size < 5)
    {
        /* =C= in a small ring  */
        if (at2[iat].valence == 2 &&
             pVA[iat].cMinRingSize <= 5 &&
             at2[iat].chem_bonds_valence == 4)
        {
            return 1;
        }
    }
    else
    {
        if (at2[iat].valence == 2 &&
             pVA[iat].cMinRingSize &&
             pVA[iat].cMinRingSize <= min_ring_size &&
             at2[iat].chem_bonds_valence == 3)
        {
            return 1;
        }
        nCurRingSize = nMinRingSize = min_ring_size + 1;
        if (( at2[iat].valence == 2 || at2[iat].valence == 3 ) &&
             at2[iat].chem_bonds_valence == at2[iat].valence + 1)
        {
            for (j = 0; j < at2[iat].valence; j++)
            {
                nCurRingSize = is_bond_in_Nmax_memb_ring( at2, iat, j, pbfsq->q,
                                             pbfsq->nAtomLevel,
                                             pbfsq->cSource, (AT_RANK) nMinRingSize /* max ring size */ );
                if (0 < nCurRingSize && nCurRingSize < nMinRingSize)
                {
                    nMinRingSize = nCurRingSize;
                }
            }
            return ( 0 <= nCurRingSize ) ? ( nMinRingSize <= min_ring_size ) : nCurRingSize;
        }
    }

    return 0;
}


/****************************************************************************/
int FixMobileHRestoredStructure( CANON_GLOBALS *pCG,
                                 INCHI_CLOCK *ic,
                                 ICHICONST INPUT_PARMS *ip,
                                 STRUCT_DATA *sd,
                                 BN_STRUCT *pBNS,
                                 BN_DATA *pBD,
                                 StrFromINChI *pStruct,
                                 inp_ATOM *at,
                                 inp_ATOM *at2,
                                 inp_ATOM *at3,
                                 VAL_AT *pVA,
                                 ALL_TC_GROUPS *pTCGroups,
                                 T_GROUP_INFO **ppt_group_info,
                                 inp_ATOM **ppat_norm,
                                 inp_ATOM **ppat_prep,
                                 INChI *pInChI[],
                                 long num_inp,
                                 int bHasSomeFixedH,
                                 int *pnNumRunBNS,
                                 int *pnTotalDelta,
                                 int forbidden_edge_mask,
                                 int forbidden_stereo_edge_mask )
{
    /*--------- process extra or missing Fixed-H on non-tautomeric atoms ------*/
    /* at2 should be the most recently restored atom, Fixed-H */
    int i, j, k, iat, delta, cur_success, ret = 0; /* djb-rwth: removing redundant variables/code */
    CMP2MHINCHI c2i;
    CMP2MHINCHI *pc2i = &c2i;

    EDGE_LIST AllChargeEdges, CurrEdges, CurrEdges2, CurrEdges3, TautEdges, NFlowerEdges, OtherNFlowerEdges, FixedLargeRingStereoEdges;
    EDGE_LIST  *pEdgeList = NULL;

    EdgeIndex e;
    BNS_EDGE  *pe;
    Vertex v1, v2, vPlusMinus;
    BNS_VERTEX *pv1, *pv2;

    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

    int        forbidden_edge_mask_inv = ~forbidden_edge_mask; /* djb-rwth: removing redundant variables */

    INCHI_HEAPCHK

    AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &CurrEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &NFlowerEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &CurrEdges2, EDGE_LIST_CLEAR );
    AllocEdgeList( &CurrEdges3, EDGE_LIST_CLEAR );
    AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &TautEdges, EDGE_LIST_CLEAR );

    /* djb-rwth: removing redundant code */

    if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
    {
        goto exit_function;  /* no fixed-H found */
    }
    /* taut group edges */
    for (i = 0; i < pTCGroups->num_tgroups; i++)
    {
        pv1 = pBNS->vert + ( v1 = pTCGroups->pTCG[i].nVertexNumber ); /* t-group vertex */ /* djb-rwth: ignoring LLVM warning: see comments below */
        for (j = 0; j < pv1->num_adj_edges; j++)
        {
            /* e, pe - tautomeric atom edge; pv2 - endpoint vertex */
            /* Note: pe, pv2, v1 are not used here; they are to show how to traverse t-group */
            pv2 = pBNS->vert + ( pe = pBNS->edge + ( e = pv1->iedge[j] ) )->neighbor1; /* djb-rwth: ignoring LLVM warning: see comments above */
            if ((ret = AddToEdgeList( &TautEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
        }
    }
    /* charge and flower edges */
    for (i = 0; i < pStruct->num_atoms; i++)
    {
        if (( e = pVA[i].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
            ( ret = AddToEdgeList( &AllChargeEdges, e, INC_ADD_EDGE ) ))
        {
            goto exit_function;
        }
        if (( e = pVA[i].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden)
        {
            if ((ret = AddToEdgeList( &AllChargeEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }

            /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
            if (pVA[i].cNumValenceElectrons == 5 && !pVA[i].cMetal && /* N, P, As */
                 NO_VERTEX != ( j = GetChargeFlowerUpperEdge( pBNS, pVA, e ) ))
            {

                if (!pBNS->edge[j].forbidden && pBNS->edge[j].flow)
                {
                    if ((ret = AddToEdgeList( &AllChargeEdges, j, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                    if ((ret = AddToEdgeList( &NFlowerEdges, j, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
                else
                {
                    if ((ret = AddToEdgeList( &OtherNFlowerEdges, j, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_function;
                    }
                }
            }
        }
    }
    if (forbidden_stereo_edge_mask)
    {
        for (i = 0; i < pStruct->num_atoms; i++)
        {
            for (j = 0; j < at2[i].valence; j++)
            {
                if (pBNS->edge[k = pBNS->vert[i].iedge[j]].forbidden == forbidden_stereo_edge_mask)
                {
                    int nMinRingSize = is_bond_in_Nmax_memb_ring( at2, i, j, pStruct->pbfsq->q,
                                                             pStruct->pbfsq->nAtomLevel,
                                                             pStruct->pbfsq->cSource, 99 /* max ring size */ );
                    if (0 < nMinRingSize && ( ret = AddToEdgeList( &FixedLargeRingStereoEdges, k, INC_ADD_EDGE ) ))
                    {
                        goto exit_function;
                    }
                }
            }
        }
    }

    INCHI_HEAPCHK

    if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
    {
        goto exit_function;
    }

    INCHI_HEAPCHK
    if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
    {
        goto exit_function;
    }

    INCHI_HEAPCHK




    if (pc2i->nNumTgInChI == 1 && ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ) &&
            pc2i->nNumTgDBNMinusRevrs + pc2i->nNumTgNHMinusRevrs == 0 && pc2i->nNumTgOMinusInChI &&
            !( pTCGroups->pTCG[0].tg_RestoreFlags & TGRF_MINUS_FIRST ))
    {            /*----------------------------------------------------*/
        /* case 01: restored has -O(-) and does not have N(-) */
        /*          endpoints defined by the original InChI   */
        /*          restored has single taut. group or more   */
        /*          tautomeric endpoints.                     */
        /* Solution: move (-) from endp. -O(-) to endpoints N */
        /*-------
---------------------------------------------*/
        pTCGroups->pTCG[0].tg_RestoreFlags |= TGRF_MINUS_FIRST;
        /* recalculate InChI from the structure */
        if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
            ppt_group_info, ppat_norm, ppat_prep ) ))
        {
            goto exit_function;
        }
        if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function;
        }
        if (!pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed)
        {
            goto exit_function;  /* no fixed-H found */
        }
        if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
        {
            goto exit_function;
        }
        if (!pc2i->bHasDifference)
        {
            goto exit_function; /* nothing to do */
        }
    }
    if (pc2i->nNumTgInChI == 1 && ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ) &&
         pc2i->nNumTgDBNMinusRevrs + pc2i->nNumTgNHMinusRevrs == 0 && pc2i->nNumTgOMinusInChI == 0)
    {
        /*-------------------------------------------------------*/
        /* case 02: restored has no -O(-) and does not have N(-) */
        /*          restored has single taut. group or more      */
        /*          tautomeric endpoints.                        */
        /* Solution: >N-AB=N-  => >N(+)=AB-NH- (add H(+))        */
        /* Solution: >N-AB=NH  => >N(+)=AB-NH2 (add H(+))        */
        /*      SB_N_III  DB_N_III                               */
        /*-------------------------------------------------------*/
        int iat_SB_N_III[MAX_DIFF_MOBH], iat_DB_N_III[MAX_DIFF_MOBH];
        int num_SB_N_III = 0, num_DB_N_III = 0, k1, k2;
        CurrEdges.num_edges = 0;
        cur_success = 0;
        for (i = 0; i < pStruct->num_atoms; i++)
        {
            iat = i;
            if (pVA[iat].cNumValenceElectrons == 5 && pVA[i].cPeriodicRowNumber == 1 &&
                 !at2[iat].endpoint && !at2[iat].charge && !at2[iat].radical)
            {
                if (num_DB_N_III < MAX_DIFF_MOBH && !at2[iat].num_H &&
                     at2[iat].valence == 2 &&
                     at2[iat].chem_bonds_valence == 3 &&
                     !at2[iat].sb_parity[0] &&  /* do not eliminate stereobonds */
                     ( e = pVA[iat].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                     pBNS->edge[e].cap && !pBNS->edge[e].flow)
                {
                    /* -N= */
                    iat_DB_N_III[num_DB_N_III++] = iat;
                }
                else
                {
                    if (num_DB_N_III < MAX_DIFF_MOBH && 1 == at2[iat].num_H &&
                         at2[iat].valence == 1 &&
                         at2[iat].chem_bonds_valence == 2 &&
                         !at2[iat].sb_parity[0] &&  /* do not eliminate stereobonds */
                         ( e = pVA[iat].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                         pBNS->edge[e].cap && !pBNS->edge[e].flow)
                    {
                        /* -N= */
                        iat_DB_N_III[num_DB_N_III++] = iat;
                    }
                    else
                    {
                        if (num_SB_N_III < MAX_DIFF_MOBH && !at2[iat].num_H &&
                                at2[iat].valence == 3 &&
                                at2[iat].chem_bonds_valence == 3 &&
                                ( e = pVA[iat].nCPlusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                                pBNS->edge[e].cap && pBNS->edge[e].flow)
                        {
                            /* -N< */
                            iat_SB_N_III[num_SB_N_III++] = iat;
                            if ((ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_function;
                            }
                        }
                    }
                }
            }
        }
        if (num_DB_N_III && num_SB_N_III)
        {
            EdgeIndex ieMinus;
            BNS_EDGE  *peMinus;
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
            for (i = 0; i < num_DB_N_III && !cur_success; i++)
            {
                iat = iat_DB_N_III[i];
                e = pBNS->edge[k1 = pBNS->vert[iat].iedge[0]].flow ? k1 :
                    pBNS->edge[k2 = pBNS->vert[iat].iedge[1]].flow ? k2 : NO_VERTEX;
                if (e == NO_VERTEX)
                {
                    continue; /* should not happen */
                }
                ieMinus = pVA[iat].nCMinusGroupEdge - 1;
                peMinus = pBNS->edge + ieMinus;
                pe = pBNS->edge + e;
                if (!pe->flow)
                    continue;
                pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                pe->forbidden |= forbidden_edge_mask;     /* fix double bond */
                peMinus->forbidden &= forbidden_edge_mask_inv; /* allow negative charge */
                delta = 1;
                pe->flow -= delta; /* remove (-) from AB-O(-) */
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2 * delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                    (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 2) /* djb-rwth: addressing LLVM warnings */
                {
                    /* Added (-)charge -N= and (+) to -N< => nDeltaCharge == 2 */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if (ret > 0)
                    {
                        /* djb-rwth: removing redundant code */
                        cur_success++; /* 01 */

                        /* eliminate (-) charge and add H */
                        pv1 = pBNS->vert + ( v1 = peMinus->neighbor1 );      /* atom */
                        pv2 = pBNS->vert + ( v2 = peMinus->neighbor12 ^ v1 );/* (=) vertex */ /* djb-rwth: ignoring LLVM warning: consistency of the code */
                        /* effectively eliminate (-) edge by setting its cap=flow= 0 */
                        peMinus->cap--;
                        peMinus->flow--;
                        pv1->st_edge.cap--;
                        pv1->st_edge.flow--;
                        pv2->st_edge.cap--;
                        pv2->st_edge.flow--;
                        pBNS->tot_st_flow -= 2;
                        pBNS->tot_st_cap -= 2;
                        /* add H */
                        pStruct->at[iat].num_H++;
                        /* register total charge increase */
                        pTCGroups->total_charge++;
                        pStruct->nNumRemovedProtonsByRevrs -= 1;
                    }
                }
                else
                {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    peMinus->forbidden |= forbidden_edge_mask;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2 * delta;
                }
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            CurrEdges.num_edges = 0; /* clear current edge list */

            if (cur_success)
            {
                /* djb-rwth: removing redundant code */
                /* recalculate InChI from the structure */
                /* recalculate InChI from the structure */
                if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep ) ))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if (!pc2i->bHasDifference)
                {
                    goto exit_function; /* nothing to do */
                }
            }
        }
    }
    if (pc2i->nNumTgInChI == 1 && ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ) && /* ADP */
        pc2i->nNumTgMInChI == 0 && pc2i->nNumTgNInChI && pc2i->nNumTgOInChI)
    {
        /*-------------------------------------------------------*/
        /* case 03: restored has N and O endpoints, no (-) endp  */
        /* case 04: original has single taut. group or more      */
        /*          tautomeric endpoints.                        */
        /* Solution: 1. Move taut attachment from O to N         */
        /* Solution: 2. Replace the attachment with (-)          */
        /*      SB_N_III  DB_N_III                               */
        /*-------------------------------------------------------*/
        /*
          int iat_SB_N_III[MAX_DIFF_MOBH], iat_DB_N_III[MAX_DIFF_MOBH];
          int num_SB_N_III = 0, num_DB_N_III = 0, k1, k2,
        */
        int itg, j1, j2, bAction = 0;
        BNS_VERTEX *pTg, *pvEndp, *pvEndp2, *pvCent; /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
        Vertex     vEndp, vEndp2, vCent;
        BNS_EDGE   *peTg, *peTg2, *peCent1, *peCent2;
        EdgeIndex  eTg, eTg2;

        CurrEdges.num_edges = 0;
        CurrEdges2.num_edges = 0;
        cur_success = 0;

        /* 1st attempt: -NH-=O => -N(-)-=O  or -N=-OH => -N(-)-=O */
        for (itg = 0; itg < pTCGroups->num_tgroups && !cur_success; itg++)
        {
            pTg = pBNS->vert + pTCGroups->pTCG[itg].nVertexNumber;
            for (i = 0; i < pTg->num_adj_edges && !cur_success; i++)
            {
                pvEndp = pBNS->vert + ( vEndp = ( peTg = pBNS->edge + ( eTg = pTg->iedge[i] ) )->neighbor1 ); /* djb-rwth: ignoring LLVM warning: value used */
                eTg2 = -1;
                if (pVA[vEndp].cNumValenceElectrons == 6 && peTg->cap)
                {
                    /* endpoint -OH or =O found; search for a possible centerpoint */
                    for (j1 = 0; j1 < at2[vEndp].valence && eTg2 < 0; j1++)
                    {
                        peCent1 = pBNS->edge + pvEndp->iedge[j1]; /* edge from O to a centerpoint */
                        pvCent = pBNS->vert + ( vCent = peCent1->neighbor12 ^ vEndp ); /* centerpoint */
                        if (at2[vCent].endpoint || !peCent1->cap ||
                             peCent1->flow + ( peTg->cap == peTg->flow ) != 1)
                        {
                            continue;
                        }
                        /* search for another endpoint, N, around vCent */
                        for (j2 = 0; j2 < at2[vCent].valence; j2++)
                        {
                            peCent2 = pBNS->edge + pvCent->iedge[j2];
                            pvEndp2 = pBNS->vert + ( vEndp2 = peCent2->neighbor12 ^ vCent ); /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
                            if (!peCent2->cap || peCent2->flow + peCent1->flow != 1 ||
                                 at2[vEndp2].endpoint != itg + 1 ||
                                 pVA[vEndp2].cNumValenceElectrons != 5 ||
                                 0 > ( j = pVA[vEndp2].nTautGroupEdge - 1 ) ||
                                 ( peTg2 = pBNS->edge + j )->forbidden ||
                                 peCent2->flow + ( peTg2->cap == peTg2->flow ) != 1)
                            {
                                continue;
                            }
                            eTg2 = j;
                            break; /* found OH-C=N- or O=C-NH- */
                        }
                    }
                }
                if (eTg2 >= 0)
                {
                    /*--------------------------------------------
                                    tg                        tg
                            eTg //\ eTg2              eTg / \\eTg2
                                //  \                     /   \\
                        vEndp HO--C==N vEndp2 -->  vEndp O==C--NH vEndp2
                                ^ ^ ^                     ^ ^ ^
                            eCent1 | eCent2           eCent1 | eCent2
                                    vCent                     vCent

                    additional action: -OH-C=N- => O=C-NH-
                    -------------------------------------------*/
                    if (0 == peTg->cap - peTg->flow && 1 == peTg2->cap - peTg2->flow &&
                         0 == peCent1->flow && 1 == peCent2->flow)
                    {
                        peTg->flow--;          /* 03 prepare */
                        peTg2->flow++;
                        peCent2->flow--;
                        peCent1->flow++;
                        bAction |= 1; /* switched H position */
                    }
                    if (1 == peTg->cap - peTg->flow && 0 == peTg2->cap - peTg2->flow &&
                         1 == peCent1->flow && 0 == peCent2->flow)
                    {
                        /* replace -NH- with -N(-)- */
                        pTCGroups->pTCG[itg].tg_num_H--;
                        pTCGroups->pTCG[itg].tg_num_Minus++;
                        pTCGroups->pTCG[itg].tg_RestoreFlags |= TGRF_MINUS_FIRST;
                        pTCGroups->pTCG[itg].tg_set_Minus = vEndp2 + 1;
                        pStruct->ti.t_group[itg].num[1] ++; /* increment number of (-), keep number of taut attachments */
                        pTCGroups->total_charge--;
                        pTCGroups->tgroup_charge--;
                        pStruct->nNumRemovedProtonsByRevrs += 1;
                        bAction |= 2; /* single NH (at2[vEndp2]) replaced with N(-) */
                        cur_success++; /* 03/04 */
                    }
                }
            }
        }

        if (0 == pc2i->nNumTgNHInChI + pc2i->nNumTgNH2InChI && pc2i->nNumTgOHInChI && !cur_success)
        {
            /* transfer an attachement to N */
            for (itg = 0; itg < pTCGroups->num_tgroups; itg++)
            {
                pTg = pBNS->vert + pTCGroups->pTCG[itg].nVertexNumber;
                for (i = 0; i < pTg->num_adj_edges; i++)
                {
                    pvEndp = pBNS->vert + ( vEndp = ( peTg = pBNS->edge + ( eTg = pTg->iedge[i] ) )->neighbor1 );
                    if (pVA[vEndp].cNumValenceElectrons == 6 &&
                         at2[vEndp].valence == at2[vEndp].chem_bonds_valence &&
                         peTg->flow && peTg->flow == peTg->cap)
                    {
                        /* endpoint -OH found; save the tautomeric group edge */
                        if ((ret = AddToEdgeList( &CurrEdges, eTg, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_function;
                        }
                    }
                    else
                    {
                        if (pVA[vEndp].cNumValenceElectrons == 5 &&
                             pVA[vEndp].cPeriodicRowNumber == 1 &&
                             at2[vEndp].valence + 1 == at2[vEndp].chem_bonds_valence &&
                             peTg->cap && peTg->flow + 1 == peTg->cap)
                        {
                            /* endpoint -N= or =NH found, check for -N=-OH */
                            e = -1;
                            for (j1 = 0; j1 < at2[vEndp].valence && e < 0; j1++)
                            {
                                peCent1 = pBNS->edge + pvEndp->iedge[j1];
                                if (peCent1->flow == 1)
                                {
                                    /* double bond */
                                    pvCent = pBNS->vert + ( vCent = peCent1->neighbor12 ^ vEndp );
                                    if (at2[vCent].endpoint)
                                        continue;
                                    for (j2 = 0; j2 < at2[vCent].valence; j2++)
                                    {
                                        peCent2 = pBNS->edge + pvCent->iedge[j2];
                                        pvEndp2 = pBNS->vert + ( vEndp2 = peCent2->neighbor12 ^ vCent ); /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
                                        if (peCent2->flow || at2[vEndp2].endpoint != itg + 1 ||
                                             pVA[vEndp2].cNumValenceElectrons != 6 ||
                                             0 >= ( e = pVA[vEndp2].nTautGroupEdge - 1 ) ||
                                             pBNS->edge[e].forbidden || !pBNS->edge[e].flow)
                                        {
                                            e = -1;
                                            continue;
                                        }
                                        /*********************/
                                        /* found -N=X-OH     */
                                        /*    vEndp ^ vEndp2 */
                                        /*          vCent    */
                                        /*********************/
                                        /* save this -OH taut edge */
                                        if ((ret = AddToEdgeList( &CurrEdges2, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                                        {
                                            goto exit_function;
                                        }
                                        break;
                                    }
                                }
                            }
                            if (e < 0 && ( ret = AddToEdgeList( &CurrEdges, eTg, INC_ADD_EDGE ) ))
                            {
                                goto exit_function;
                            }
                        }
                    }
                }
            }
            /* rearrange the flows */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            SetForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
            SetForbiddenEdgeMask( pBNS, &CurrEdges2, forbidden_edge_mask );
            pEdgeList = CurrEdges2.num_edges ? &CurrEdges2 : CurrEdges.num_edges ? &CurrEdges : NULL;

            for (i = 0; pEdgeList && i < pEdgeList->num_edges && !cur_success; i++)
            {
                pe = pBNS->edge + pEdgeList->pnEdges[i]; /* pe->flow = 1 <=> -OH */
                if (!pe->flow)
                    continue;
                pv1 = pBNS->vert + ( v1 = pe->neighbor1 );       /* -OH atom */
                pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 ); /* t-group vertex */
                /* locate the t-group */
                for (itg = 0; itg < pTCGroups->num_tgroups; itg++)
                {
                    if (v2 == pTCGroups->pTCG[itg].nVertexNumber)
                    {
                        break;
                    }
                }
                if (itg == pTCGroups->num_tgroups)
                {
                    /* tgroup not found -- should not happen */
                    continue;
                }

                delta = 1;
                pe->flow -= delta; /* add one attachment to  */
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2 * delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                    (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 2) /* djb-rwth: addressing LLVM warning */
                {
                    /* Added (-)charge -N= and (+) to -N< => nDeltaCharge == 2 */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if (ret > 0)
                    {
                        /* djb-rwth: removing redundant code */
                        cur_success++; /* 03 */
                        /* replace -NH- with -N(-)- */
                        pTCGroups->pTCG[itg].tg_num_H--;
                        pTCGroups->pTCG[itg].tg_num_Minus++;
                        pTCGroups->pTCG[itg].tg_RestoreFlags |= TGRF_MINUS_FIRST;
                        pStruct->ti.t_group[itg].num[1] ++;
                        pTCGroups->total_charge--;
                        pTCGroups->tgroup_charge--;
                        pStruct->nNumRemovedProtonsByRevrs += 1;
                        bAction |= 4; /* H in the 1st available NH was replaced with (-) */
                    }
                }
                else
                {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2 * delta;
                }
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
        }
        else
        {
            if (pc2i->nNumTgNHInChI + pc2i->nNumTgNH2InChI && pc2i->nNumTgOInChI && !cur_success)
            {
                /* change an attachement to N from H to (-) */
                for (itg = 0; itg < pTCGroups->num_tgroups && !cur_success; itg++)
                {
                    pTg = pBNS->vert + pTCGroups->pTCG[itg].nVertexNumber;
                    for (i = 0; i < pTg->num_adj_edges && !cur_success; i++)
                    {
                        pvEndp2 = pBNS->vert + ( vEndp2 = ( peTg = pBNS->edge + pTg->iedge[i] )->neighbor1 ); /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
                        if (pVA[vEndp2].cNumValenceElectrons == 5 && pVA[vEndp2].cPeriodicRowNumber == 1 &&
                             at2[vEndp2].valence == at2[vEndp2].chem_bonds_valence &&
                             peTg->flow && peTg->flow == peTg->cap)
                        {
                            /* endpoint -NHn found; change its charge */
                            cur_success++; /* 04 */
                            /* replace -NH- with -N(-)- */
                            pTCGroups->pTCG[itg].tg_num_H--;
                            pTCGroups->pTCG[itg].tg_num_Minus++;
                            pTCGroups->pTCG[itg].tg_RestoreFlags |= TGRF_MINUS_FIRST;
                            pTCGroups->pTCG[itg].tg_set_Minus = vEndp2 + 1;
                            pStruct->ti.t_group[itg].num[1] ++;
                            pTCGroups->total_charge--;
                            pTCGroups->tgroup_charge--;
                            pStruct->nNumRemovedProtonsByRevrs += 1;
                            bAction |= 8; /* manually set (-) charge to NH atom, vEndp2 */
                        }
                    }
                }
            }
        }
        if (cur_success)
        {
            /* djb-rwth: removing redundant code */
            /* recalculate InChI from the structure */
            if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                ppt_group_info, ppat_norm, ppat_prep ) ))
            {
                goto exit_function;
            }
            if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            if (pStruct->One_ti.num_t_groups == 1 && pStruct->One_ti.t_group[0].num[1])
            {
                /* this method did not work: no alt path from N(-) to =O */
                itg = 0;
                if (bAction & ( 8 | 2 ))
                {
                    /* roll back NH -> N(-) replacement; H move from OH to N is not undone */
                    pTCGroups->pTCG[itg].tg_num_H++;
                    pTCGroups->pTCG[itg].tg_num_Minus--;
                    pTCGroups->pTCG[itg].tg_RestoreFlags &= ~TGRF_MINUS_FIRST;
                    pTCGroups->pTCG[itg].tg_set_Minus = 0;
                    pStruct->ti.t_group[itg].num[1] --;
                    pTCGroups->total_charge++;
                    pTCGroups->tgroup_charge++;
                    pStruct->nNumRemovedProtonsByRevrs -= 1;
                    cur_success--;
                }
                else
                {
                    if (bAction & 4)
                    {
                        pTCGroups->pTCG[itg].tg_num_H++;
                        pTCGroups->pTCG[itg].tg_num_Minus--;
                        pTCGroups->pTCG[itg].tg_RestoreFlags &= ~TGRF_MINUS_FIRST;
                        pStruct->ti.t_group[itg].num[1] --;
                        pTCGroups->total_charge++;
                        pTCGroups->tgroup_charge++;
                        pStruct->nNumRemovedProtonsByRevrs -= 1;
                        cur_success--;
                    }
                    else
                    {
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                }
                /* recalculate InChI from the structure */
                if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                    ppt_group_info, ppat_norm, ppat_prep ) ))
                {
                    goto exit_function;
                }
                if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
                if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
            }
            if (!pc2i->bHasDifference)
            {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if (pc2i->nNumTgInChI == 1 && ( pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1 ) && /* ADP */
        pc2i->nNumTgMInChI == 0 && ( pc2i->nNumTgNInChI || pc2i->nNumTgOInChI ) &&
        NO_VERTEX != ( vPlusMinus = GetPlusMinusVertex( pBNS, pTCGroups, 1, 1 ) ))
    {
        /*---------------------------------------------------------------------------*/
        /* case 05: restored has N endpoints, no (-) endpoints                       */
        /*          original has single taut. group or more                          */
        /*          tautomeric endpoints.                                            */
        /* Solution: Find -N< and allow (+) charge change                            */
        /*           Fix all charges and taut attachments exept                      */
        /*           =N- and =O (taut. endpoints)                                    */
        /*           Increment st_edge.cap on (+/-) vertex => add (+) charge to -N<  */
        /*           Increment tot. charge in other places                           */
        /*           Increment t-group st_edge.cap                                   */
        /*           Run BNS                                                         */
        /*                                                                           */
        /*      (+/-)*               (+/-)           Result:                         */
        /*        |                    ||                                            */
        /*        |                    ||            - Added (+) to -N<              */
        /*       (+)super             (+)super       - Added attachment point to O   */
        /*        ||                   |                                             */
        /*        ||          =>       |             To make this attachment H,      */
        /*       (Y)                  (Y)            increment                       */
        /*        |                    ||            pTCGroups->pTCG[itg].tg_num_H   */
        /*        |                    ||                                            */
        /*       (+)hetero            (+)hetero      Technical details:              */
        /*         \\                   \            increase capacities of          */
        /*           N                    N(+)       edges to (+/-) otherwise        */
        /*           |                    ||         flow may not be able to         */
        /*   *(t)--O=R.            (t)==O-R.         increase                        */
        /*                                                                           */
        /*                                                                           */
        /*---------------------------------------------------------------------------*/
        int itg;
        BNS_VERTEX *pTg, *pvEndp; /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
        Vertex     vEndp, vTg;
        BNS_EDGE   *peTg;
        EdgeIndex  eTg;
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];

        CurrEdges.num_edges = 0;
        CurrEdges2.num_edges = 0;
        cur_success = 0;
        /* find -N< and non-taut =N- or =O */
        for (i = 0; i < pStruct->num_atoms; i++)
        {
            iat = nCanon2AtnoRevrs[i];
            /* -N< */
            if (!at2[iat].endpoint && !at2[iat].charge && !at2[iat].radical && !at2[iat].num_H &&
                 pVA[i].cNumValenceElectrons == 5 && pVA[i].cPeriodicRowNumber == 1 &&
                 0 <= ( e = pVA[iat].nCPlusGroupEdge - 1 ) && pBNS->edge[e].flow && !pBNS->edge[e].forbidden)
            {

                if ((ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
            }
        }
        if (!CurrEdges.num_edges)
        {
            goto exit_case_05;
        }
        /* find taut -N= and =O */
        for (itg = 0; itg < pTCGroups->num_tgroups && !cur_success; itg++)
        {
            CurrEdges2.num_edges = 0;
            pTg = pBNS->vert + ( vTg = pTCGroups->pTCG[itg].nVertexNumber );
            for (i = 0; i < pTg->num_adj_edges; i++)
            {
                pvEndp = pBNS->vert + ( vEndp = ( peTg = pBNS->edge + ( eTg = pTg->iedge[i] ) )->neighbor1 ); /* djb-rwth: ignoring LLVM warning: variable used to store initialisation values */
                if (at2[vEndp].charge || at2[vEndp].radical || peTg->cap - peTg->flow != 1)
                {
                    continue;
                }
                /* t-group edges to -N= and =O */
                if ((ret = AddToEdgeList( &CurrEdges2, eTg, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_function;
                }
            }
            if (!CurrEdges2.num_edges)
            {
                goto exit_case_05;
            }
            /* fix all charge edges except -N< and all taut. edges except =O and =N- */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            SetForbiddenEdgeMask( pBNS, &TautEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges2, forbidden_edge_mask );
            delta = 1;
            /* Increment st_edge.cap on (+/-) vertex */
            pBNS->vert[vPlusMinus].st_edge.cap += delta;
            /* Increment st_edge.cap on t-group */
            pTg->st_edge.cap += delta;
            /* total cap count */
            pBNS->tot_st_cap += 2 * delta;

            v1 = vPlusMinus;
            v2 = vTg;

            /* increase capacities of edges to Y  */
            for (i = 0; i < pBNS->vert[vPlusMinus].num_adj_edges; i++)
            {
                j = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i]].neighbor12 ^ vPlusMinus;
                for (k = 0; k < pBNS->vert[j].num_adj_edges; k++)
                {
                    pBNS->edge[pBNS->vert[j].iedge[k]].cap += delta;
                }
            }

            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warnings */
            {
                /* Added (+)charge to -N< => nDeltaCharge == 1 */
                /* Flow change on pe (-)charge edge (atom B-O(-)) is not known to RunBnsTestOnce()) */
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                if (ret > 0)
                {
                    /* djb-rwth: removing redundant code */
                    cur_success++; /* 01 */
                    /* update bookkeeping */
                    pTCGroups->total_charge += delta;
                    pTCGroups->pTCG[itg].edges_cap += delta;
                    pTCGroups->pTCG[itg].tg_num_H += delta;
                    pStruct->nNumRemovedProtonsByRevrs -= delta;
                }
            }
            else
            {
                pBNS->vert[vPlusMinus].st_edge.cap -= delta;
                pTg->st_edge.cap -= delta;
                /*pTCGroups->pTCG[itg].edges_cap     -= delta;*/ /* ???bug??? - commented out 2006-03-22 */
                pBNS->tot_st_cap -= 2 * delta;
                /* decrease capacities of edges to Y  */
                for (i = 0; i < pBNS->vert[vPlusMinus].num_adj_edges; i++)
                {
                    j = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i]].neighbor12 ^ vPlusMinus;
                    for (k = 0; k < pBNS->vert[j].num_adj_edges; k++)
                    {
                        pBNS->edge[pBNS->vert[j].iedge[k]].cap -= delta;
                    }
                }
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &TautEdges, forbidden_edge_mask );
        }
        if (cur_success)
        {
            /* djb-rwth: removing redundant code */
            /* recalculate InChI from the structure */
            if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                ppt_group_info, ppat_norm, ppat_prep ) ))
            {
                goto exit_function;
            }
            if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            if (!pc2i->bHasDifference)
            {
                goto exit_function; /* nothing to do */
            }
        }

    exit_case_05:;
    }

    while (pc2i->nNumDiffMobH && pc2i->nChargeMobHRevrs > pc2i->nChargeMobHInChI)
    {
        /*----------------------------------------------------*/
        /* case 06: restored has extra H attached to -O(-)    */
        /*          while the chrge should be on C, most pro- */
        /*          bably in a small ring.ut. group or more   */
        /*          tautomeric endpoints.                     */
        /* Solution: move (-) from O to C                     */
        /*----------------------------------------------------*/
        int iO, mode;
        EdgeIndex e2;
        BNS_EDGE  *pe2;
        cur_success = 0;
        for (i = 0; !cur_success && i < pc2i->len_c2at; i++)
        {

            if (pc2i->c2at[i].nMobHRevrs == pc2i->c2at[i].nMobHInChI + 1 &&
                 pc2i->c2at[i].nNumHRevrs == pc2i->c2at[i].nMobHInChI &&
                 !pc2i->c2at[i].endptInChI && !pc2i->c2at[i].endptRevrs &&
                 at2[iO = pc2i->c2at[i].atomNumber].charge == -1 &&
                 0 <= ( e = pVA[iO].nCMinusGroupEdge - 1 ) && ( pe = pBNS->edge + e )->flow)
            {
                /* try suitable atoms C */
                /* first look for =C= in a small ring */
                for (mode = 4; !cur_success && mode <= 8; mode++)
                {

                    if (mode == 8)
                        mode = 99;

                    for (iat = 0; !cur_success && iat < pStruct->num_atoms; iat++)
                    {

                        if (!at2[iat].charge && !at2[iat].radical &&
                             pVA[iat].cNumValenceElectrons == 4 &&
                             0 <= ( e2 = pVA[iat].nCMinusGroupEdge - 1 ) && !( pe2 = pBNS->edge + e2 )->flow &&
                             0 < bIsUnsatCarbonInASmallRing( at2, pVA, iat, pStruct->pbfsq, mode ))
                        {

                            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                            /* allow negative charge on the chosen carbon */
                            pe2->forbidden &= forbidden_edge_mask_inv;

                            delta = 1;
                            if (!pe->flow)
                                continue;
                            pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                            pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );
                            pe->flow -= delta;
                            pv1->st_edge.flow -= delta;
                            pv2->st_edge.flow -= delta;
                            pBNS->tot_st_flow -= 2 * delta;

                            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                            if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                                (vPathEnd == v2 && vPathStart == v1) ) && nDeltaCharge == 1) /* djb-rwth: addressing LLVM warning */
                            {
                                /* Added (-)charge to unsaturated C => nDeltaCharge == 2 */
                                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                                if (ret > 0)
                                {
                                    /* djb-rwth: removing redundant code */
                                    cur_success++; /* 01 */
                                    /* djb-rwth: removing redundant code */
                                }
                            }
                            else
                            {
                                pe->forbidden |= forbidden_edge_mask;
                                pe->flow += delta;
                                pv1->st_edge.flow += delta;
                                pv2->st_edge.flow += delta;
                                pBNS->tot_st_flow += 2 * delta;
                            }
                            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                        }
                    }
                }
            }
        }
        if (cur_success)
        {
            /* recalculate InChI from the structure */
            if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                ppt_group_info, ppat_norm, ppat_prep ) ))
            {
                goto exit_function;
            }
            if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            if (!pc2i->bHasDifference)
            {
                goto exit_function; /* nothing to do */
            }
        }
        else
        {
            break;
        }
    }
    if (pc2i->len_c2at && pc2i->nChargeMobHRevrs > pc2i->nChargeMobHInChI)
    {
        /*------------------------------------------------------------------*/
        /* case 07: -NO2 are to be tautomeric but they are not AND          */
        /*          InChI has a SINGLE tautomeric group                     */
        /*                                                                  */
        /*                   (-)O                   (-)O                    */
        /* Solution: convert     \                      \                   */
        /*                        N-X=...-Z(-)   =>      N(+)=X- ...=Z      */
        /*                      //                      /                   */
        /*                     O                    (-)O                    */
        /*                                                                  */
        /*                     O                       O                    */
        /*        or            \\                      \\                  */
        /*                        N-X=...-Z(-)    =>      N=X-  ...=Z       */
        /*                      //                       /                  */
        /*                     O                     (-)O                   */
        /*                                                                  */
        /*                                                                  */
        /*  (a) move (-) from other tautomeric atom to O in O=N-X           */
        /*          or   from other atom that has to be tautomeric          */
        /*               but is not                                         */
        /*  (b) create (+) [ion pair creation] on N as in                   */
        /*                                                                  */
        /*       OH             OH                                          */
        /*      /              /                                            */
        /*  -C=N     =>  =C-N(+)                                            */
        /*     \\             \\                                            */
        /*       O              O                                           */
        /*                                                                  */
        /*------------------------------------------------------------------*/
        int num_DB_O = 0;
        short iat_DB_O[MAX_DIFF_FIXH], iat_NO2[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        /*
        AT_NUMB  *nAtno2CanonRevrs = pStruct->nAtno2Canon[0];
        */
        inp_ATOM *at_Mobile_H_Revrs = ( pStruct->pOne_norm_data[0] &&
                                     pStruct->pOne_norm_data[0]->at ) ? pStruct->pOne_norm_data[0]->at : NULL;

        int iN, one_success;
        BNS_EDGE *peDB_O_Minus;
        int neigh, nNumO, nNumOthers;
#define CHG_SET_WRONG_TAUT_N   0
#define CHG_SET_WRONG_TAUT_O   1
#define CHG_SET_WRONG_TAUT_ALL 2
#define CHG_LAST_SET           2 /* the last index in trying */
#define CHG_SET_O_FIXED        3
#define CHG_SET_NUM            4
        EDGE_LIST ChangeableEdges[CHG_SET_NUM];
        memset( ChangeableEdges, 0, sizeof( ChangeableEdges ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        /* equivalent to AllocEdgeList( &EdgeList, EDGE_LIST_CLEAR ); */
        /*
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        CurrEdges.num_edges = 0; /* clear current edge list */
        cur_success = 0;
        for (i = 0; i < pStruct->num_atoms; i++)
        {
            iat = nCanon2AtnoRevrs[i];
            if ( /* orig. InChI info: taut in orig. InChI =O located in -NO2 that is not taut in Reconstructed InChI */
                 num_DB_O < MAX_DIFF_FIXH &&
                 pVA[iat].cNumValenceElectrons == 6 /* O, S, Se, Te */ &&
                 ( !at2[iat].endpoint /*|| pc2i->c2at[i].nMobHInChI*/ ) &&
                 ( e = pVA[iat].nCMinusGroupEdge - 1 ) >= 0 && !pBNS->edge[e].forbidden &&
                 at2[iat].num_H == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                 /* reversed structure info: */
                 !( at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint ) /*|| pc2i->c2at[i].nMobHRevrs*/ &&
                 !at2[iat].charge &&
                 at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 &&
                 /* find whether it belongs to NO2 */
                 pVA[iN = at2[iat].neighbor[0]].cNumValenceElectrons == 5 &&
                 at2[iN].valence == 3 && ( at2[iN].charge == 0 || at2[iN].charge == 1 ) &&
                 at2[iN].chem_bonds_valence == 5 - at2[iN].charge)
            {
                /* find the second O */
                nNumO = nNumOthers = 0;
                for (k = 0; k < at2[iN].valence; k++)
                {
                    neigh = at2[iN].neighbor[k];
                    if (neigh == iat)
                    {
                        continue;
                    }
                    if (pVA[neigh].cNumValenceElectrons == 6 &&
                         !at2[neigh].endpoint &&
                         !( at_Mobile_H_Revrs && at_Mobile_H_Revrs[neigh].endpoint ) &&
                         at2[neigh].valence == 1 && at2[neigh].num_H == 0 &&
                         at2[neigh].radical == 0 && ( at2[neigh].charge == 0 || at2[neigh].charge == -1 ) &&
                         at2[neigh].chem_bonds_valence - at2[neigh].charge == 2)
                    {
                        nNumO++;
                    }
                    else
                    {
                        if (at2[iN].bond_type[k] == BOND_TYPE_SINGLE &&
                              at2[neigh].valence > 1 &&
                              at2[neigh].valence < at2[neigh].chem_bonds_valence)
                        {
                            nNumOthers++;
                        }
                    }
                }
                if (nNumO != 1 || nNumOthers != 1)
                {
                    continue;
                }
                for (k = 0; k < num_DB_O; k++)
                {
                    if (iat_NO2[k] == iN)
                    {
                        break;
                    }
                }
                if (k == num_DB_O)
                {
                    iat_NO2[num_DB_O] = iN;
                    iat_DB_O[num_DB_O++] = iat;
                }
                /* save the =O (-)-edge to avoid interference */
                if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                {
                    goto exit_case_07;
                }
            }
        }
        if (num_DB_O)
        {
            /* search for falsely tautomeric negatively charged atoms N and O */
            for (i = 0; i < pc2i->len_c2at; i++)
            {
                iat = pc2i->c2at[i].atomNumber;
                if (pc2i->c2at[i].endptRevrs && !pc2i->c2at[i].endptInChI &&
                     pc2i->c2at[i].nAtChargeRevrs == -1 &&
                     0 <= ( e = pVA[iat].nCMinusGroupEdge - 1 ) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                     0 > FindInEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e ))
                {
                    if (pc2i->c2at[i].nValElectr == 6)
                    {
                        if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_WRONG_TAUT_O], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                        {
                            goto exit_case_07;
                        }
                    }
                    else
                        if (pc2i->c2at[i].nValElectr == 5)
                        {
                            if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_WRONG_TAUT_N], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                            {
                                goto exit_case_07;
                            }
                        }
                    if ((ret = AddToEdgeList( &ChangeableEdges[CHG_SET_WRONG_TAUT_ALL], e, INC_ADD_EDGE ))) /* djb-rwth: addressing LLVM warning */
                    {
                        goto exit_case_07;
                    }
                }
            }
            /* ------- finally, try to move charges from O=N --------------*/
            for (i = 0; i < num_DB_O; i++)
            {
                int nDeltaChargeExpected;
                one_success = 0;
                delta = 1;
                iat = iat_DB_O[i];
                peDB_O_Minus = pBNS->edge + ( (long long)pVA[iat].nCMinusGroupEdge - 1 ); /* djb-rwth: cast operator added */
                pe = pBNS->edge + pBNS->vert[iat].iedge[0];

                if (!pe->flow)
                    continue;
                pv1 = pBNS->vert + ( v1 = pe->neighbor1 );
                pv2 = pBNS->vert + ( v2 = pe->neighbor12 ^ v1 );

                pe->forbidden |= forbidden_edge_mask;

                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2 * delta;

                for (k = 0; !one_success && k <= CHG_LAST_SET; k++)
                {
                    if (!ChangeableEdges[k].num_edges)
                    {
                        continue;
                    }
                    nDeltaChargeExpected = 0;

                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    RemoveForbiddenEdgeMask( pBNS, &ChangeableEdges[k], forbidden_edge_mask );
                    /* allow (-) charge to move to N=O */
                    peDB_O_Minus->forbidden &= forbidden_edge_mask_inv;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if (ret == 1 && ( (vPathEnd == v1 && vPathStart == v2) ||
                        (vPathEnd == v2 && vPathStart == v1) ) &&
                                      nDeltaCharge == nDeltaChargeExpected) /* djb-rwth: addressing LLVM warnings */
                    {
                        /* Move (-) charge to =O and remove it an endpoint => nDeltaCharge == 0 */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if (ret > 0)
                        {
                            /* djb-rwth: removing redundant code */
                            one_success++; /* 07 */
                        }
                    }
                    INCHI_HEAPCHK
                }
                cur_success += one_success;

                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                pe->forbidden &= forbidden_edge_mask_inv;

                if (!one_success)
                {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2 * delta;
                }
            }
        }
    exit_case_07:
        for (i = 0; i < CHG_SET_NUM; i++)
        {
            AllocEdgeList( &ChangeableEdges[i], EDGE_LIST_FREE );
        }

        CurrEdges.num_edges = 0; /* clear current edge list */
        if (cur_success)
        {
            /* djb-rwth: removing redundant code */
            /* recalculate InChI from the structure */
            if (0 > ( ret = MakeOneInChIOutOfStrFromINChI2( pCG, ic, ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                ppt_group_info, ppat_norm, ppat_prep ) ))
            {
                goto exit_function;
            }
            if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */
            {
                goto exit_function;
            }
            if (!pc2i->bHasDifference)
            {
                goto exit_function; /* nothing to do */
            }
        }
#undef CHG_SET_NOOH
#undef CHG_SET_WRONG_TAUT
#undef CHG_SET_TAUT
#undef CHG_LAST_SET
#undef CHG_SET_O_FIXED
#undef CHG_SET_NUM
    }

exit_function:
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &CurrEdges, EDGE_LIST_FREE );
    AllocEdgeList( &CurrEdges2, EDGE_LIST_FREE );
    AllocEdgeList( &CurrEdges3, EDGE_LIST_FREE );
    AllocEdgeList( &NFlowerEdges, EDGE_LIST_FREE );
    AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_FREE );
    AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_FREE );
    AllocEdgeList( &TautEdges, EDGE_LIST_FREE );

    return ret;
}

#endif
