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

/*#define CHECK_WIN32_VC_HEAP*/
#include "mode.h"

#if( READ_INCHI_STRING == 1 )

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

/* local types */

/* types for TgDiffHChgFH */
#define fNumRPosChgH 0 /* number of positive charges on endpoints that have H in at2[] */
#define fNumRPosChgU 1 /* number of positive charges on endpoints that have no H in at2[] */
#define fNumRNegChgO 2 /* number of negative charges on O endpoints */
#define fNumRNegChgN 3 /* number of negative charges on N endpoints */
#define fNumRNeutrlH 4 /* number of neutral endp that have H in at2[] */

#define fNumNPosChgH 5 /* number of positive charges on endpoints that have H in atf[] */
#define fNumNPosChgU 6 /* number of positive charges on endpoints that have no H in atf[] */
#define fNumNNegChgO 7 /* number of negative charges on O endpoints */
#define fNumNNegChgN 8 /* number of negative charges on N endpoints */
#define fNumNNeutrlH 9 /* number of neutral endp that have H in atf[] */

#define fNumAllChgT 10 /* total  number of fNum... */

typedef struct tagTgDiffHChgFH {
    short  itg; /* t-group index; endpoint = itg+1 */
    short  nNumHInchi;  /* number of H in t-group from orig. InChI */
    short  nNumHRevrs;  /* number of H in at2[] */
    short  nNumHNorml;  /* number of H in Normalized atfMobile_H_Revrs[] */
    short  nNumMInchi;  /* number of (-) in InChI */
    short  nNumMRevrs;  /* number of (-) in at2[] */
    short  nNumMNorml;  /* number of (-) in atf[] */
    short  nNumPRevrs;  /* number of (+) in at2[] */
    short  nNumPNorml;  /* number of (+) in Normalized atfMobile_H_Revrs[] */
    short n[fNumAllChgT]; /* all numbers */
    short i[fNumAllChgT]; /* all indices */
} TgDiffHChgFH;

/* local prototypes */
static int FillTgDiffHChgFH( TgDiffHChgFH tdhc[], int max_tdhc, inp_ATOM at2[], inp_ATOM atf[],
                            AT_NUMB  *nCanon2AtnoRevrs, VAL_AT *pVA, T_GROUP_INFO *ti, EDGE_LIST *pAtomIndList );


/************************************************************/
int bHas_N_V( inp_ATOM *at2, int num_atoms )
{
    static U_CHAR el_number_N;
    int i, num_found = 0;
    if ( !el_number_N ) {
        el_number_N = get_periodic_table_number( "N" );
    }
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( at2[i].el_number == el_number_N && !at2[i].charge &&
             !at2[i].num_H && !at2[i].radical &&
             at2[i].chem_bonds_valence == 5 &&
             (at2[i].valence==3) ) {
            num_found ++;
        }
    }
    return num_found;
}
/*************************************************************************************/
int FillTgDiffHChgFH( TgDiffHChgFH tdhc[], int max_tdhc, inp_ATOM at2[],
                      inp_ATOM atf[], AT_NUMB  *nCanon2AtnoRevrs, VAL_AT *pVA,
                      T_GROUP_INFO *ti, EDGE_LIST *pAtomIndList )
{

    int i, j, iat, itg, itg_prev, num, itg_out, bOverflow;
    EDGE_LIST IndList;   /* type, itg */
    TgDiffHChgFH cur_tdhc;
    AT_NUMB    *pEndp0;
    inp_ATOM   *at2i, *atfi;
    int         typeR, typeN, type, ret = 0, nCurIndListLen;

    AllocEdgeList( &IndList, EDGE_LIST_CLEAR );
    pAtomIndList->num_edges = 0;
    itg_out = 0;
    bOverflow = 0;
    memset( tdhc, 0, max_tdhc * sizeof(tdhc[0]) );

    for ( itg = 0; itg < ti->num_t_groups; itg ++ ) {
        memset( &cur_tdhc, 0, sizeof(cur_tdhc) );

        cur_tdhc.itg = itg;
        cur_tdhc.nNumHInchi = ti->t_group[itg].num[0] - ti->t_group[itg].num[1];
        cur_tdhc.nNumMInchi = ti->t_group[itg].num[1];
        
        pEndp0 = ti->nEndpointAtomNumber + ti->t_group[itg].nFirstEndpointAtNoPos;
        nCurIndListLen = IndList.num_edges;
        for ( j = 0; j < ti->t_group[itg].nNumEndpoints; j ++ ) {
            i = pEndp0[j];
            iat = nCanon2AtnoRevrs[i];
            
            at2i = at2 + iat;
            atfi = atf + iat;
            
            typeR = typeN = -1;
            if ( at2i->charge == 1 ) {
                if ( at2i->num_H ) {
                    typeR = fNumRPosChgH;
                } else {
                    typeR = fNumRPosChgU;
                }
                cur_tdhc.nNumPRevrs ++;
            } else
            if ( at2i->charge == -1 ) {
                if ( pVA[iat].cNumValenceElectrons == 6) {
                    typeR = fNumRNegChgO;
                } else
                if ( pVA[iat].cNumValenceElectrons == 5) {
                    typeR = fNumRNegChgN;
                }
                cur_tdhc.nNumMRevrs ++;
            } else
            if ( at2i->num_H && at2i->valence == at2i->chem_bonds_valence ) {
                typeR = fNumRNeutrlH;
            }
            cur_tdhc.nNumHRevrs += at2i->num_H;

            if ( atfi->charge == 1 ) {
                if ( atfi->num_H ) {
                    typeN = fNumNPosChgH;
                } else {
                    typeN = fNumNPosChgU;
                }
                cur_tdhc.nNumPNorml ++;
            } else
            if ( atfi->charge == -1 ) {
                if ( pVA[iat].cNumValenceElectrons == 6) {
                    typeN = fNumNNegChgO;
                } else
                if ( pVA[iat].cNumValenceElectrons == 5) {
                    typeN = fNumNNegChgN;
                }
                cur_tdhc.nNumMNorml ++;
            } else
            if ( atfi->num_H && atfi->valence == atfi->chem_bonds_valence ) {
                typeN = fNumNNeutrlH;
            }
            cur_tdhc.nNumHNorml += atfi->num_H;
            if ( at2[iat].charge < 0 || 0 < pVA[iat].nCPlusGroupEdge ) {
                if ( typeR >= 0 && (
                     (ret = AddToEdgeList( &IndList, typeR, INC_ADD_EDGE )) ||
                     (ret = AddToEdgeList( &IndList, itg, INC_ADD_EDGE )) ||
                     (ret = AddToEdgeList( &IndList, iat, INC_ADD_EDGE )) ) ) {
                    goto exit_function;
                }
                if ( typeN >= 0 && (
                     (ret = AddToEdgeList( &IndList, typeN, INC_ADD_EDGE )) ||
                     (ret = AddToEdgeList( &IndList, itg, INC_ADD_EDGE )) ||
                     (ret = AddToEdgeList( &IndList, iat, INC_ADD_EDGE )) ) ) {
                    goto exit_function;
                }
            }

        }
        if ( cur_tdhc.nNumHNorml == cur_tdhc.nNumHInchi &&
             cur_tdhc.nNumMNorml == cur_tdhc.nNumMInchi ) {
            IndList.num_edges = nCurIndListLen; /* t-group seems to be correct */
            continue;
        }
        if ( itg_out < max_tdhc ) {
            tdhc[itg_out ++] = cur_tdhc;
        } else {
            bOverflow |= 1;
            IndList.num_edges = nCurIndListLen;
            break;
        }
    }
    /* fill out atom index list */
    if ( itg_out ) {
        itg_prev = IndList.pnEdges[1]; /* the 1st saved t-group number */
        for ( type = 0; type < fNumAllChgT; type ++ ) {
            j = 0;
            for ( i = 0; i < itg_out; i ++ ) {
                num = 0;
                itg = tdhc[i].itg;
                tdhc[i].i[type] = -999; /* empty */
                while( IndList.pnEdges[j+1] == itg ) {
                    if ( IndList.pnEdges[j] == type ) {
                        if ( !num ++ ) {
                            tdhc[i].i[type] = pAtomIndList->num_edges;
                        }
                        if ( ret = AddToEdgeList( pAtomIndList, IndList.pnEdges[j+2], INC_ADD_EDGE )) {
                            goto exit_function;
                        }
                    }
                    j += 3;
                }
                tdhc[i].n[type] = num;
            }
        }
    }
    ret = itg_out;
exit_function:
    AllocEdgeList( &IndList, EDGE_LIST_FREE );
    return ret;

/*
#undef fNumRPosChgH
#undef fNumRPosChgU
#undef fNumRNegChgO
#undef fNumRNegChgN

#undef fNumNPosChgH
#undef fNumNPosChgU
#undef fNumNNegChgO
#undef fNumNNegChgN

#undef fNumAllChgT    
*/
}

/***********************************************************************************************/
int FixFixedHRestoredStructure(ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
                        StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA,
                        ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ppt_group_info, inp_ATOM **ppat_norm,
                        inp_ATOM **ppat_prep, INChI *pInChI[], long num_inp, int bHasSomeFixedH,
                        int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask, int forbidden_stereo_edge_mask)
{
    /*--------- process extra or missing Fixed-H on non-tautomeric atoms ------*/
    /* at2 should be the most recently restored atom, Fixed-H */
    int i, j, k, delta, num_try, tot_succes, cur_success, ret = 0, bAllowedNFlowerEdges=0, num_zero_ret;
    CMP2FHINCHI c2i;
    CMP2FHINCHI *pc2i = &c2i;

    EDGE_LIST AllChargeEdges, CurrEdges, SFlowerEdges, NFlowerEdges, OtherNFlowerEdges, FixedLargeRingStereoEdges;
    EDGE_LIST AllBondEdges;
    
    EdgeIndex e;
    BNS_EDGE  *pe;
    Vertex v1, v2;
    BNS_VERTEX *pv1, *pv2;

    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

    int        nNumRunBNS = 0, forbidden_edge_mask_inv = ~forbidden_edge_mask;

    INCHI_HEAPCHK

    AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &CurrEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &NFlowerEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &SFlowerEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &AllBondEdges, EDGE_LIST_CLEAR );

    tot_succes = 0;
    
    if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
        goto exit_function;  /* no fixed-H found */
    }

    for ( i = 0; i < pStruct->num_atoms; i ++ ) {
        if ( (e=pVA[i].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
             (ret = AddToEdgeList( &AllChargeEdges, e, INC_ADD_EDGE )) ) {
            goto exit_function;
        }
        if ( (e=pVA[i].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
            if ( ret = AddToEdgeList( &AllChargeEdges, e, INC_ADD_EDGE ) ) {
                goto exit_function;
            }

            /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
            if ( pVA[i].cNumValenceElectrons == 5 && !pVA[i].cMetal && /* N, P, As */
                 NO_VERTEX != (j = GetChargeFlowerUpperEdge( pBNS, pVA, e ))) {

                if ( pBNS->edge[j].forbidden ) {
                    continue;
                }

                if ( pBNS->edge[j].flow ) {
                    if ( ret = AddToEdgeList( &AllChargeEdges, j, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                    if ( ret = AddToEdgeList( &NFlowerEdges, j, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                } else {
                    if ( ret = AddToEdgeList( &OtherNFlowerEdges, j, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                }
            } else
            /* in addition, disallow N(V) creation by forbidding charge flower edge that has flow=1 */
            if ( pVA[i].cNumValenceElectrons == 6 && !pVA[i].cMetal && /* N, P, As */
                 NO_VERTEX != (j = GetChargeFlowerUpperEdge( pBNS, pVA, e ))) {

                if ( pBNS->edge[j].forbidden ) {
                    continue;
                }

                if ( pBNS->edge[j].flow ) {
                    if ( ret = AddToEdgeList( &SFlowerEdges, j, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                }
            }

        }
        for ( j = 0; j < at2[i].valence; j ++ ) {
            k = at2[i].neighbor[j];
            if ( k < i && !pBNS->edge[e=pBNS->vert[i].iedge[j]].forbidden ) {
                if ( ret = AddToEdgeList( &AllBondEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
    }
    if ( forbidden_stereo_edge_mask ) {
        for ( i = 0; i < pStruct->num_atoms; i ++ ) {
            for ( j = 0; j < at2[i].valence; j ++ ) {
                if ( pBNS->edge[k = pBNS->vert[i].iedge[j]].forbidden == forbidden_stereo_edge_mask ) {
                    int nMinRingSize = is_bond_in_Nmax_memb_ring( at2, i, j, pStruct->pbfsq->q,
                                                             pStruct->pbfsq->nAtomLevel,
                                                             pStruct->pbfsq->cSource, 99 /* max ring size */ );
                    if ( 0 < nMinRingSize && (ret = AddToEdgeList( &FixedLargeRingStereoEdges, k, INC_ADD_EDGE ))) {
                        goto exit_function;
                    }
                }
            }
        }
    }

    INCHI_HEAPCHK
    if ( ret = FillOutExtraFixedHDataRestr( pStruct ) ) {
        goto exit_function;
    }
    INCHI_HEAPCHK
    if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
        goto exit_function;
    }

    INCHI_HEAPCHK

    if ( !pc2i->bHasDifference ||
         !pc2i->len_c2at && pc2i->nNumTgRevrs == pc2i->nNumTgInChI &&
         pc2i->nNumEndpRevrs == pc2i->nNumRemHInChI &&
         pc2i->nNumEndpRevrs == pc2i->nNumEndpInChI &&
         !pc2i->nNumTgDiffMinus && !pc2i->nNumTgDiffH ) {
        goto exit_function; /* nothing to do */
    }

    /*goto exit_function;*/ /* debug only*/

    if ( pc2i->len_c2at >= 2 ) {
        /*----------------------------------------------------*/
        /* case 01: restored: O=AB-O(-)  original:  (-)O-AB=O */
        /* FixH:              0    -1                 -1    0 */
        /* MobH:              0     1                  1    0 */
        /*                         non-taut      non-taut     */
        /* O = O, S, Se; charged atoms O are not tautomeric   */
        /* Solution: move (-) from B-O(-) to O=A              */
        /*----------------------------------------------------*/
        int num_DB_O = 0, num_SB_O_Minus = 0, iat;
        short iat_DB_O[MAX_DIFF_FIXH], iat_SB_O_Minus[MAX_DIFF_FIXH];
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( pc2i->c2at[i].nValElectr == 6 /* && !pc2i->c2at[i].endptInChI -- mod#1*/ &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                if ( /* orig. InChI info: */
                     num_SB_O_Minus < MAX_DIFF_FIXH && 
                     pc2i->c2at[i].nFixHInChI ==  0 && pc2i->c2at[i].nMobHInChI == 0 &&
                     /* reversed structure info: */
                     pc2i->c2at[i].nFixHRevrs == -1 && pc2i->c2at[i].nMobHRevrs == 1 &&
                     pc2i->c2at[i].nAtChargeRevrs == -1 && !at2[iat].num_H && /* at2 is Fixed-H */
                     at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 1 ) {
                    iat_SB_O_Minus[num_SB_O_Minus ++] = iat;
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                } else
                if ( /* orig. InChI info: */
                     num_DB_O < MAX_DIFF_FIXH &&
                     pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI ==  1 &&
                     /* reversed structure info: */
                     pc2i->c2at[i].nFixHRevrs ==  0 && pc2i->c2at[i].nMobHRevrs ==  0 &&
                     pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                     at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 ) {
                    iat_DB_O[num_DB_O ++] = iat;
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                }
            }
        }
        if ( num_try = inchi_min( num_SB_O_Minus, num_DB_O ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_O_Minus && cur_success < num_try; i ++ ) {
                iat = iat_SB_O_Minus[i];
                pe   = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta; /* remove (-) from AB-O(-) */
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Added (-)charge to O=AB => nDeltaCharge == -1 */
                    /* Flow change on pe (-)charge edge (atom B-O(-)) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 01 */
                    }
                } else {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
            CurrEdges.num_edges = 0; /* clear current edge list */
        }
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->len_c2at >= 1 ) {
        /*--------------------------------------------------------------*/
        /* case 02: restored: -O(+)=AB-NH2  original:  -O-AB=NH2(+)     */
        /* FixH:               0        0               0      1        */
        /* MobH:               0        2               0      1        */
        /* O = P, As, Sb, O, S, Se, F, Cl, Br, I; not taut. in InChI    */
        /* N = N, O, S, Se, Te; has H; tautomeric or not tautomeric     */
        /* Solution: move (+) from O(+) to NH2                          */
        /*--------------------------------------------------------------*/
        int num_DB_O_Plus = 0, num_SB_NH = 0, iat;
        short iat_DB_O_Plus[MAX_DIFF_FIXH], iat_SB_NH[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : NULL;
        cur_success = 0;
        num_zero_ret = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: =NH2(+), =OH(+) */
                 num_SB_NH < MAX_DIFF_FIXH &&
                 (pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1 ||
                  pc2i->c2at[i].nValElectr == 6 ) /* N, O, S, Se, Te */ &&
                 /*!pc2i->c2at[i].endptInChI &&*/ /* <=== relaxation */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI>0 /*== 1 --modification#2*/ && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                 /* reversed structure info: */
                 pc2i->c2at[i].nFixHRevrs ==  0 && /* pc2i->c2at[i].nMobHRevrs == 0 &&*/
                 pc2i->c2at[i].nAtChargeRevrs == 0 && at2[iat].num_H &&
                 at2[iat].valence == at2[iat].chem_bonds_valence ) {
                iat_SB_NH[num_SB_NH ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom: charge=+1, no H, has double bond, P, As, O, S, Se, Te, F, Cl, Br, I */
                 num_DB_O_Plus < MAX_DIFF_FIXH &&
                 at2[iat].charge == 1 && !at2[iat].num_H &&
                 at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 (pVA[iat].cNumValenceElectrons == 6 || pVA[iat].cNumValenceElectrons == 7 ||
                  pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber > 1) &&
                 /* in orig.InChI: not an endpoint, has no H */
                 !pStruct->endpoint[i] && 
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has (+) edge */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                iat_DB_O_Plus[num_DB_O_Plus ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_DB_O_Plus, num_SB_NH ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
repeat_02_allow_NV:
            for ( i = 0; i < num_SB_NH && cur_success < num_try; i ++ ) {
                iat = iat_SB_NH[i];
                pe   = pBNS->edge + pVA[iat].nCPlusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                
                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == -1 ) {
                    /* Removed charge from O(+) => nDeltaCharge == -1 */
                    /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 02 */
                    }
                } else {
                    num_zero_ret += !ret;
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            if ( num_zero_ret == num_try && !bAllowedNFlowerEdges && NFlowerEdges.num_edges ) {
                RemoveForbiddenEdgeMask( pBNS, &NFlowerEdges, forbidden_edge_mask  );
                bAllowedNFlowerEdges = 1;
                goto repeat_02_allow_NV;
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
            bAllowedNFlowerEdges = 0;
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->len_c2at >= 1 && pc2i->nNumTgRevrs == 1 &&
         (pc2i->nNumEndpRevrs > pc2i->nNumEndpInChI || pc2i->nNumTgInChI > 1) /* ADP in Revrs */ ) {
        /*--------------------------------------------------------------*/
        /* case 03: restored: -N(-)-AB=O    original:  -N=AB-O(-)       */
        /* FixH:               0       0                0     -1        */
        /* MobH:               0       0                0      1        */
        /* O = O, S, Se; N = N;                                         */
        /* restored atoms are tautomeric; original atoms are not taut.  */
        /* restored struct has 1 t-group; original has less endpoints   */
        /*                                and possibly >1 t-groups      */  
        /* Solution: move (-) from N(-) to =O                           */
        /*           these atoms are tautomeric in restored structure   */
        /*--------------------------------------------------------------*/
        int num_SB_N_Minus = 0, num_DB_O = 0, iat;
        short iat_SB_N_Minus[MAX_DIFF_FIXH], iat_DB_O[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        /*
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: -O(-) */
                 num_DB_O < MAX_DIFF_FIXH &&
                 pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                 !pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI ==  1 &&
                 /* reversed structure info: */
                 pc2i->c2at[i].endptRevrs &&
                 pc2i->c2at[i].nFixHRevrs ==  0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                 at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 ) {
                iat_DB_O[num_DB_O ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom N: charge=-1, no H, has no double bond, endpoint */
                 num_SB_N_Minus < MAX_DIFF_FIXH &&
                 at2[iat].charge == -1 && /*!at2[iat].num_H &&*/
                 at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                 at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint &&
                 /* in orig.InChI: not an endpoint, has no H */
                 /* !pStruct->endpoint[i] && */
               /*  
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                */
                 /* has (-) edge */
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                iat_SB_N_Minus[num_SB_N_Minus ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_SB_N_Minus, num_DB_O ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_N_Minus && cur_success < num_try; i ++ ) {
                iat = iat_SB_N_Minus[i];
                pe   = pBNS->edge + pVA[iat].nCMinusGroupEdge-1; /* 2006-03-03: changed from CPlusGroupEdge */
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Added (-) charge to =O => nDeltaCharge == 1 */
                    /* Flow change on pe (-)charge edge (atom -N(-)-) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 03 */
                    }
                } else {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }
    
    if ( pc2i->nNumTgRevrs == 1 && /* pc2i->nNumRemHInChI < 0 &&*/
         (pc2i->nNumEndpRevrs > pc2i->nNumEndpInChI || pc2i->nNumTgInChI > 1) /* ADP in Revrs */ ) {
        /*--------------------------------------------------------------*/
        /* case 03a:restored: -N(-)-AB=O    original:  -N=AB-O(-)       */
        /* FixH:               0       0                0      0        */
        /* MobH:               0       0                0      0        */
        /* O = O, S, Se; N = N;                              taut       */
        /* restored atoms are tautomeric; original atom is; N may be.   */
        /* restored struct has 1 t-group; original has less endpoints   */
        /*                                and possibly >1 t-groups      */  
        /* Solution: move (-) from N(-) to =O                           */
        /*           these atoms are tautomeric in restored structure   */
        /*--------------------------------------------------------------*/
        int num_SB_N_Minus = 0, num_DB_O = 0, iat;
        short iat_SB_N_Minus[MAX_DIFF_FIXH], iat_DB_O[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        S_CHAR   *pnMobHInChI = (pInChI[1] && pInChI[1]->nNum_H)? pInChI[1]->nNum_H :
                                (pInChI[0] && pInChI[0]->nNum_H)? pInChI[0]->nNum_H : NULL;
        S_CHAR   *pnFixHInChI = pStruct->fixed_H;

        cur_success = 0;
        CurrEdges.num_edges = 0;
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom N: charge=-1, no H, has no double bond, endpoint */
                 num_SB_N_Minus < MAX_DIFF_FIXH &&
                 at2[iat].charge == -1 && /*!at2[iat].num_H &&*/
                 at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                 at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint &&
                 /* in orig.InChI: may be an endpoint, has no H */
                 /* has (-) edge */
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                iat_SB_N_Minus[num_SB_N_Minus ++] = iat;
            } else
            if ( num_DB_O < MAX_DIFF_FIXH &&
                 at2[iat].charge == 0 && /*!at2[iat].num_H &&*/
                 at2[iat].valence+1 == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 pVA[iat].cNumValenceElectrons == 6 &&
                 at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint && /* endpoint in Reconstructed */
                 (pStruct->endpoint[i] || /* endpoint or H(+) acceptor in original */
                  pnMobHInChI && pnMobHInChI[i] == 1 && pnFixHInChI && pnFixHInChI[i] == -1 ) &&
                 /* has (-) edge */
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                iat_DB_O[num_DB_O ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }    
        }
        if ( num_try = inchi_min( num_SB_N_Minus, num_DB_O ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            /* allow charge transfer to all found =O */
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_N_Minus && cur_success < num_try; i ++ ) {
                iat = iat_SB_N_Minus[i];
                pe   = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Added (-) charge to =O => nDeltaCharge == 1 */
                    /* Flow change on pe (-)charge edge (atom -N(-)-) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 03a */
                    }
                } else {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->len_c2at >= 1 && pc2i->nNumTgInChI == 1 && /* ADP in InChI */
         (pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1) ) {
        /*--------------------------------------------------------------*/
        /* case 04: restored: OH(+)=AB-O- OH- orig.  HO-AB=O(+)- OH-    */
        /* FixH:               1       0   0          1      0   1      */
        /* MobH:               0       0   1          0      0   0      */
        /*                 non-taut.                taut        taut    */
        /*                                    ADP: one t-group or more endpoints */
        /* O(+) = N, P, As, As, O, S, Se; OH = N, O, S, Se, Te          */
        /* Solution: move (+) from O(+) to NH2                          */
        /*--------------------------------------------------------------*/
        int num_SB_Neutr = 0, num_DB_Charged = 0, iat;
        short iat_SB_Neutr[MAX_DIFF_FIXH], iat_DB_Charged[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        cur_success = 0;
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom: charge=+1, has H, has double bond, N, O, S, Se, Te */
                 num_DB_Charged < MAX_DIFF_FIXH &&
                 at2[iat].charge == 1 && at2[iat].num_H &&
                 at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 (pVA[iat].cNumValenceElectrons == 6 ||
                  pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1) &&
                 /* in orig.InChI: an endpoint, has fixed-H */
                 pStruct->endpoint[i] && 
                 (pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 /*!(nMobHInChI && nMobHInChI[i] ) &&*/
                 /* has (+) edge */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {

                iat_DB_Charged[num_DB_Charged ++] = iat;
                
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            } else
            if ( /* in restored atom: charge=0, has no H, has no double bond, N, P, O, S, Se, Te */
                 num_SB_Neutr < MAX_DIFF_FIXH &&
                 at2[iat].charge == 0 && !at2[iat].num_H &&
                 at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 (pVA[iat].cNumValenceElectrons == 6 ||
                  pVA[iat].cNumValenceElectrons == 5 ) &&
                 /* in orig.InChI: an endpoint, has fixed-H */
                 /* pStruct->endpoint[i] && */
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has (+) edge */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 &&
                 0 == pBNS->edge[e].forbidden ) {

                iat_SB_Neutr[num_SB_Neutr ++] = iat;
                
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_SB_Neutr, num_DB_Charged ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_Neutr && cur_success < num_try; i ++ ) {
                iat = iat_SB_Neutr[i];
                pe   = pBNS->edge + pVA[iat].nCPlusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == -1 ) {
                    /* Removed charge from O(+) => nDeltaCharge == -1 */
                    /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 04 */
                    }
                } else {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->len_c2at > 1 ) {
        /*--------------------------------------------------------------*/
        /* case 05: restored:  O=AB-NH      original:(-)O-AB=NH(+)      */
        /* FixH:               0     0                 -1     1         */
        /* MobH:               0     1                  1     0         */
        /* O = O, S, Se; N = N, O, S, Se, Te; all atoms not tautomeric  */
        /* Solution: Separate charges                                   */
        /*--------------------------------------------------------------*/
        int num_DB_O = 0, num_SB_NH = 0, iat;
        short iat_DB_O[MAX_DIFF_FIXH], iat_SB_NH[MAX_DIFF_FIXH];
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: =NH2(+), =OH(+) */
                 num_SB_NH < MAX_DIFF_FIXH &&
                 (pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1 ||
                  pc2i->c2at[i].nValElectr == 6 ) /* N, O, S, Se, Te */ &&
                 !pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == 1 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                 /* reversed structure info: */
                 pc2i->c2at[i].nFixHRevrs ==  0 && pc2i->c2at[i].nMobHRevrs &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && at2[iat].num_H &&
                 !pc2i->c2at[i].endptRevrs &&
                 at2[iat].valence == at2[iat].chem_bonds_valence ) {
                iat_SB_NH[num_SB_NH ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            } else
            if ( /* orig. InChI info: -O(-) */
                 num_DB_O < MAX_DIFF_FIXH &&
                 (pc2i->c2at[i].nValElectr == 6 ) /* O, S, Se, Te */ &&
                 !pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI ==  1 &&
                 /* reversed structure info: */
                 pc2i->c2at[i].nFixHRevrs ==  0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                 !pc2i->c2at[i].endptRevrs &&
                 at2[iat].valence + 1 == at2[iat].chem_bonds_valence ) {
                iat_DB_O[num_DB_O ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_DB_O, num_SB_NH ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_NH && cur_success < num_try; i ++ ) {
                iat = iat_SB_NH[i];
                pe   = pBNS->edge + pVA[iat].nCPlusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Added charge to =O => nDeltaCharge == 1 */
                    /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 05 */
                    }
                } else {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pStruct->fixed_H && pStruct->endpoint && pc2i->nChargeFixHInChI > 0 && pc2i->nChargeFixHInChI > pc2i->nChargeMobHInChI ) {
        /*----------------------------------------------------------*/
        /* case 06c: restored -NH- or -NH(+)  orig: -NH-            */
        /*  Fixed-H            1       1             0              */
        /*  Mobile-H           0       0             1              */
        /*                     not tautomeric    not tautomeric     */
        /*           has adjacent (+)                               */
        /*           charges                                        */
        /*  Solution: move (+) charges to the -NH- unless it already*/
        /*            N = N, O, S, Se, Te                           */
        /*            has (+) charge blocked by adjacent (+)        */
        /*----------------------------------------------------------*/
        int iat;
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        /*
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
                                      pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds?
                                      pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds : NULL;
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : NULL;
        */
        EDGE_LIST CurChargeEdges;
        EdgeIndex e2;
        cur_success = 0;
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
        CurrEdges.num_edges = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            /* atoms -NH- from which H(+) were removed by the Normalization in orig. InChI */
            iat = pc2i->c2at[i].atomNumber;
            if ( (pc2i->c2at[i].nValElectr == 6 ||
                  pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1) &&
                 !pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                if ( /* orig. InChI info: -NH- */
                     pc2i->c2at[i].nFixHInChI == 1 && pc2i->c2at[i].nMobHInChI == 0 &&
                     /* reversed structure info: */
                     pc2i->c2at[i].nFixHRevrs == 0 && pc2i->c2at[i].nMobHRevrs == 1 && /* was not removed */
                     /*pc2i->c2at[i].nAtChargeRevrs == 0 &&*/ at2[iat].num_H && /* at2 is Fixed-H */
                     at2[iat].valence == at2[iat].chem_bonds_valence ) {
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                }
            }
        }
        for ( i = 0; i < pStruct->num_atoms; i ++ ) {
            /* find adjacent charged atoms */
            iat = nCanon2AtnoRevrs[i];
            if ( pStruct->endpoint[i] || at2[iat].charge != 1 || at2[iat].radical || pVA[iat].cMetal ) {
                continue;
            }
            if ( 0 <= (e=pVA[iat].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden && !pBNS->edge[e].flow && pVA[iat].cNumValenceElectrons >= 5 ) {
                /* positively charged atom */
                for ( j = 0; j < at2[iat].valence; j ++ ) {
                    if ( at2[k=(int)at2[iat].neighbor[j]].charge == 1 && !pVA[k].cMetal &&
                         0 <= (e2=pVA[k].nCPlusGroupEdge-1) && !pBNS->edge[e2].forbidden && !pBNS->edge[e2].flow) {
                        if ( 0 > FindInEdgeList( &CurrEdges, e ) &&
                             0 > FindInEdgeList( &CurChargeEdges, e ) &&
                             ( ret = AddToEdgeList( &CurChargeEdges, e, INC_ADD_EDGE ) ) ) {
                            goto exit_case_06c;
                        }
                        if ( 0 > FindInEdgeList( &CurrEdges, e2 ) &&
                             0 > FindInEdgeList( &CurChargeEdges, e2 ) &&
                             ( ret = AddToEdgeList( &CurChargeEdges, e2, INC_ADD_EDGE ) ) ) {
                            goto exit_case_06c;
                        }
                    }
                }
            }
        }
        if ( num_try = inchi_min( CurrEdges.num_edges, CurChargeEdges.num_edges ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurChargeEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < CurrEdges.num_edges && cur_success < num_try; i ++ ) {
                e = CurrEdges.pnEdges[i];
                pe   = pBNS->edge + e; /* (+)charge edge of -NH- or -OH */
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->flow -= delta; /* add (+) to -NHm */
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == -1 ) {
                    /* Removed (+)charge from -NH- => nDeltaCharge == -1 */
                    /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 06c */
                    }
                } else {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
exit_case_06c:
        CurrEdges.num_edges = 0; /* clear current edge list */
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_FREE );
        if ( ret < 0 ) {
            goto exit_function;
        }
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if (  pc2i->len_c2at >= 2 ) {
        /*------------------------------------------------------------*/
        /* case 06d: restored: XH(+)=-AB-NH    orig.: XH-=AB=NH(+)    */
        /* FixH:                1       1 0          0 1      1       */
        /* MobH:                0    taut 1          1 taut   0       */
        /*                                                            */
        /*                                                            */
        /* N  = N, O, S, Se; atoms N are not tautomeric in orig InChI */
        /* X  = N, O, S, Se, Te, F, Cl, Br, I; atom X is non-taut     */
        /* Solution: move (+) from X  to NH                           */
        /*------------------------------------------------------------*/
        int iat;
        /*
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
                                      pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds?
                                      pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds : 
                                      pStruct->pOne_norm_data[TAUT_NON]->at;
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        EDGE_LIST CurChargeEdges;
        cur_success = 0;
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
        CurrEdges.num_edges = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            /* XH(+) */
            if ( /* reconstructed: non-taut and (+) */
                 (pc2i->c2at[i].nMobHRevrs+1 == pc2i->c2at[i].nFixHRevrs &&
                 pc2i->c2at[i].nFixHRevrs > 0 && !pc2i->c2at[i].endptRevrs &&
                 pc2i->c2at[i].nAtChargeRevrs == 1 &&
                 /* original InChI: non-taut & has H or an endpoint, has Fixed H */
                 (!pc2i->c2at[i].nFixHInChI && pc2i->c2at[i].nMobHInChI == pc2i->c2at[i].nFixHRevrs ||
                 pc2i->c2at[i].nFixHInChI == pc2i->c2at[i].nFixHRevrs && pc2i->c2at[i].endptInChI ))  &&
                 0 <= (e=pVA[iat].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden && !pBNS->edge[e].flow) {

                if (ret = AddToEdgeList( &CurChargeEdges, e, INC_ADD_EDGE )) {
                    goto exit_case_06d;
                }
            } else
            /* -NH- */
            if ( /* original InChI: has H and is not an endpoint */
                (pc2i->c2at[i].nMobHInChI+1 == pc2i->c2at[i].nFixHInChI && 
                 pc2i->c2at[i].nFixHInChI > 0 && !pc2i->c2at[i].endptInChI &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 &&
                 /* reconstructed InChI: non-taut & has H or an endpoint, has Fixed H */
                 (!pc2i->c2at[i].nFixHRevrs && pc2i->c2at[i].nMobHRevrs == pc2i->c2at[i].nFixHInChI ||
                 pc2i->c2at[i].nFixHRevrs == pc2i->c2at[i].nFixHInChI && pc2i->c2at[i].endptRevrs ))  &&
                 0 <= (e=pVA[iat].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden && 
                 pBNS->edge[e].flow) {

                if (ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE )) {
                    goto exit_case_06d;
                }
            }
        }
        if ( num_try = inchi_min( CurrEdges.num_edges, CurChargeEdges.num_edges ) ) {
            /* detected; attempt to fix */
            int bSFlowerEdgesMayBeForbidden = (SFlowerEdges.num_edges > 0);
            int bSFlowerEdgesIsForbidden;
            for ( bSFlowerEdgesIsForbidden = bSFlowerEdgesMayBeForbidden;
                     0 <= bSFlowerEdgesIsForbidden; bSFlowerEdgesIsForbidden -- ) {
                if ( bSFlowerEdgesIsForbidden ) {
                    /* on the 1st pass disallow -S(+)= => =S=, allow only -S(+)= => -S- */
                    SetForbiddenEdgeMask( pBNS, &SFlowerEdges, forbidden_edge_mask );
                }
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurChargeEdges, forbidden_edge_mask  );
                delta = 1;
                for ( i = 0; i < CurrEdges.num_edges && cur_success < num_try; i ++ ) {
                    e = CurrEdges.pnEdges[i];
                    pe   = pBNS->edge + e; /* (+)charge edge of -NH- or -OH */
                    if ( !pe->flow )
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->flow -= delta; /* add (+) to -NHm */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2*delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == -1 ) {
                        /* Removed (+)charge from -NH- => nDeltaCharge == -1 */
                        /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if ( ret > 0 ) {
                            nNumRunBNS ++;
                            cur_success ++; /* 06d */
                        }
                    } else {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2*delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                RemoveForbiddenEdgeMask( pBNS, &SFlowerEdges, forbidden_edge_mask );
            }

        }
exit_case_06d:
        CurrEdges.num_edges = 0; /* clear current edge list */
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_FREE );
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }


    if ( pc2i->len_c2at >= 2 ) {
        /*--------------------------------------------------------*/
        /* case 06: restored: NHn(+)=AB-NHm  orig.: NHn-AB=NHm(+) */
        /* FixH:               1        0            0     1      */
        /* MobH:              n-1       m            n    m-1     */
        /* N = N, O, S, Se; atoms N are not tautomeric            */
        /* Solution: move (+) from NHn(+) to NHn                  */
        /*--------------------------------------------------------*/
        int num_DB_NHn_Plus = 0, num_SB_NHm_Neutr = 0, iat;
        short iat_DB_NHn_Plus[MAX_DIFF_FIXH], iat_SB_NHm_Neutr[MAX_DIFF_FIXH];
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( (pc2i->c2at[i].nValElectr == 6 ||
                  pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1) &&
                 !pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                if ( /* orig. InChI info: NHm */
                     num_SB_NHm_Neutr < MAX_DIFF_FIXH &&
                     pc2i->c2at[i].nFixHInChI == 1 && /*pc2i->c2at[i].nMobHInChI == 0 &&*/
                     /* reversed structure info: */
                     pc2i->c2at[i].nFixHRevrs == 0 && /*pc2i->c2at[i].nMobHRevrs == 1 &&*/
                     pc2i->c2at[i].nAtChargeRevrs == 0 && at2[iat].num_H && /* at2 is Fixed-H */
                     at2[iat].valence == at2[iat].chem_bonds_valence ) {
                    iat_SB_NHm_Neutr[num_SB_NHm_Neutr ++] = iat;
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                } else
                if ( /* orig. InChI info: */
                     num_DB_NHn_Plus < MAX_DIFF_FIXH &&
                     pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI &&*/
                     /* reversed structure info: */
                     pc2i->c2at[i].nFixHRevrs ==  1 && /*pc2i->c2at[i].nMobHRevrs ==  0 &&*/
                     pc2i->c2at[i].nAtChargeRevrs == 1 && at2[iat].num_H &&
                     at2[iat].valence < at2[iat].chem_bonds_valence ) {
                    iat_DB_NHn_Plus[num_DB_NHn_Plus ++] = iat;
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                }
            }
        }
        if ( num_try = inchi_min( num_SB_NHm_Neutr, num_DB_NHn_Plus ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_NHm_Neutr && cur_success < num_try; i ++ ) {
                iat = iat_SB_NHm_Neutr[i];
                pe   = pBNS->edge + pVA[iat].nCPlusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta; /* add (+) to -NHm */
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == -1 ) {
                    /* Removed (+)charge from -NHn => nDeltaCharge == -1 */
                    /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 06 */
                    }
                } else {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( (pc2i->nNumTgInChI > pc2i->nNumTgRevrs && pc2i->nNumTgRevrs == 1 ||
          pc2i->nNumEndpInChI < pc2i->nNumEndpRevrs ) &&
          pStruct->nNumRemovedProtonsMobHInChI == pStruct->One_ti.tni.nNumRemovedProtons &&
          pStruct->fixed_H && pStruct->endpoint && pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds ) {
        /*----------------------------------------------------------*/
        /* case 06a: restored: N'(+)=-AB-NH    orig.: N'-=AB=NH(+)  */
        /* FixH:               0         1            0      1      */
        /* MobH:               0         0            0      0      */
        /*                    single t-group      multiple t-groups */
        /* N  = N, O, S, Se; atoms N are not tautomeric             */
        /* N' = N            atom N' is not tautomeric              */
        /* Solution: move (+) from N' to NH                         */
        /*----------------------------------------------------------*/
        int iat;
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        /*
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        */
        inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
                                      pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds?
                                      pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds : NULL;
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        EDGE_LIST CurChargeEdges;
        cur_success = 0;
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
        CurrEdges.num_edges = 0;
        for ( i = 0; i < pStruct->num_atoms; i ++ ) {
            iat = nCanon2AtnoRevrs[i];
            if ( pStruct->endpoint[i] ) {
                continue;
            }
            /* -NH-, -OH */
            if ( pStruct->fixed_H[i] && !nMobHInChI[i] &&
                 at2[iat].charge == 0 && at2[iat].radical == 0 &&
                 0 <= (e=pVA[iat].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                 (ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) {
                goto exit_case_06a;
            } else
            /* >N(+)= */
            if ( at2[iat].charge == 1 && !at2[iat].num_H &&
                 pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                 atfMobile_H_Revrs && atfMobile_H_Revrs[iat].charge == 0 &&
                 0 <= (e=pVA[iat].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden && !pBNS->edge[e].flow &&
                 (ret = AddToEdgeList( &CurChargeEdges, e, INC_ADD_EDGE ))) {
                goto exit_case_06a;
            }
        }
        if ( num_try = inchi_min( CurrEdges.num_edges, CurChargeEdges.num_edges ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurChargeEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < CurrEdges.num_edges && cur_success < num_try; i ++ ) {
                e = CurrEdges.pnEdges[i];
                pe   = pBNS->edge + e; /* (+)charge edge of -NH- or -OH */
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->flow -= delta; /* add (+) to -NHm */
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == -1 ) {
                    /* Removed (+)charge from -NH- => nDeltaCharge == -1 */
                    /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 06a */
                    }
                } else {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
exit_case_06a:
        CurrEdges.num_edges = 0; /* clear current edge list */
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_FREE );
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }
    if ( (pc2i->nNumTgInChI > pc2i->nNumTgRevrs && pc2i->nNumTgRevrs == 1 ||
          pc2i->nNumEndpInChI < pc2i->nNumEndpRevrs ) &&
          (pStruct->nNumRemovedProtonsMobHInChI == pStruct->One_ti.tni.nNumRemovedProtons ||
           pStruct->nNumRemovedProtonsMobHInChI > pStruct->One_ti.tni.nNumRemovedProtons ) &&
          pStruct->fixed_H && pStruct->endpoint && pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds ) {
        /*----------------------------------------------------------*/
        /* case 06b: restored: X(+)=-AB-NH    orig.: X-=AB=NH(+)    */
        /* FixH:               0        1 1          0      1       */
        /* MobH:               0        0 t          0      0       */
        /*                    single t-group      multiple t-groups */
        /*                                        or no t-groupd    */
        /* N  = N, O, S, Se; atoms N are not tautomeric             */
        /* X  = O, S, Se, Te, F, Cl, Br, I; atom X is not tautomeric*/
        /* Solution: move (+) from X  to NH                         */
        /*----------------------------------------------------------*/
        int iat;
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        /*
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        */
        inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
                                      pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds?
                                      pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds : NULL;
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        EDGE_LIST CurChargeEdges;
        cur_success = 0;
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
        CurrEdges.num_edges = 0;
        for ( i = 0; i < pStruct->num_atoms; i ++ ) {
            iat = nCanon2AtnoRevrs[i];
            if ( pStruct->endpoint[i] ) {
                continue;
            }
            /* -NH-, -OH */
            if ( pStruct->fixed_H[i] && !nMobHInChI[i] &&
                 at2[iat].charge == 0 && at2[iat].radical == 0 &&
                 0 <= (e=pVA[iat].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                 (ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ))) {
                goto exit_case_06b;
            } else
            /* X(+)= */
            if ( at2[iat].charge == 1 && !at2[iat].num_H &&
                 (pVA[iat].cNumValenceElectrons == 6 || pVA[iat].cPeriodicRowNumber == 7) &&
                 atfMobile_H_Revrs && atfMobile_H_Revrs[iat].charge == 1 &&
                 0 <= (e=pVA[iat].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden && !pBNS->edge[e].flow &&
                 (ret = AddToEdgeList( &CurChargeEdges, e, INC_ADD_EDGE ))) {
                goto exit_case_06b;
            }
        }
        if ( num_try = inchi_min( CurrEdges.num_edges, CurChargeEdges.num_edges ) ) {
            /* detected; attempt to fix */
            int bSFlowerEdgesMayBeForbidden = (SFlowerEdges.num_edges > 0);
            int bSFlowerEdgesIsForbidden;
            for ( bSFlowerEdgesIsForbidden = bSFlowerEdgesMayBeForbidden;
                     0 <= bSFlowerEdgesIsForbidden; bSFlowerEdgesIsForbidden -- ) {
                if ( bSFlowerEdgesIsForbidden ) {
                    /* on the 1st pass disallow -S(+)= => =S=, allow only -S(+)= => -S- */
                    SetForbiddenEdgeMask( pBNS, &SFlowerEdges, forbidden_edge_mask );
                }
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurChargeEdges, forbidden_edge_mask  );
                delta = 1;
                for ( i = 0; i < CurrEdges.num_edges && cur_success < num_try; i ++ ) {
                    e = CurrEdges.pnEdges[i];
                    pe   = pBNS->edge + e; /* (+)charge edge of -NH- or -OH */
                    if ( !pe->flow )
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->flow -= delta; /* add (+) to -NHm */
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2*delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == -1 ) {
                        /* Removed (+)charge from -NH- => nDeltaCharge == -1 */
                        /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if ( ret > 0 ) {
                            nNumRunBNS ++;
                            cur_success ++; /* 06b */
                        }
                    } else {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2*delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                RemoveForbiddenEdgeMask( pBNS, &SFlowerEdges, forbidden_edge_mask );
            }

        }
exit_case_06b:
        CurrEdges.num_edges = 0; /* clear current edge list */
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_FREE );
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }



    if ( pc2i->nNumTgInChI > 1 &&
          (pStruct->nNumRemovedProtonsMobHInChI > 0 || pStruct->ti.tni.nNumRemovedProtons > 0 ) &&
          pStruct->fixed_H && pStruct->endpoint && 
          pStruct->pOne_norm_data[TAUT_YES] && pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds ) {
        /*----------------------------------------------------------*/
        /* case 06e:restored: XHn(+)=-AB-YHm  orig.: XHn-=AB=YHm(+) */
        /* FixH:               1          0           1      1      */
        /* MobH:               0          1           t      t      */
        /*                   non-taut atoms       multiple t-groups */
        /*                                                          */
        /* 1. orig. t-group has more H on its endpoints counted     */
        /*          in atf and has no (+) on endpoint that has H    */
        /* 2. orig. t-group has less H on its endpoints counted     */
        /*          in atf and has (+) on endpoint that has H       */
        /*          in reconstructed struct and less H in atf       */
        /* Solution: move (+) from (2) to atom in (1) that has H    */
        /*                                                          */
        /*   tg1  reconstr:   XHn and more H than in orig t-group   */
        /*             atf:   XHn                                   */
        /*   tg2  reconstr:   XHm(+) and less H than in             */
        /*             atf:   XH(m-1)           orig in t-group     */
        /*                                                          */
        /* N  = N, O, S, Se; atoms N are not tautomeric             */
        /* X  = O, S, Se, Te, F, Cl, Br, I; atom X is not tautomeric*/
        /* Solution: move (+) from X  to NH                         */
        /*----------------------------------------------------------*/

        int iat, nNumWrongTg, jjoffs, jj, nNum2RemovePlus, nNum2AddPlus, nNum2MovePlus;
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        /*
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        */
        inp_ATOM *atfMobile_H_Revrs = pStruct->pOne_norm_data[TAUT_YES] &&
                                      pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds?
                                            pStruct->pOne_norm_data[TAUT_YES]->at_fixed_bonds :
                                      pStruct->pOne_norm_data[TAUT_YES] &&
                                      pStruct->pOne_norm_data[TAUT_YES]->at?
                                            pStruct->pOne_norm_data[TAUT_YES]->at : NULL;
        /*
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        EDGE_LIST CurChargeEdges /* source of (+)*/, EndpList;
        TgDiffHChgFH tdhc[MAX_DIFF_FIXH];
        BNS_VERTEX *pv1n, *pv2n;
        BNS_EDGE   *pe1n, *pe2n;
        Vertex      v1n, v2n;

        cur_success = 0;
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_CLEAR );
        AllocEdgeList( &EndpList, EDGE_LIST_CLEAR );
        CurrEdges.num_edges = 0; /* receptors of (+) */
        if ( !atfMobile_H_Revrs ) {
            goto exit_case_06e;
        }
        nNumWrongTg = FillTgDiffHChgFH( tdhc, MAX_DIFF_FIXH, at2, atfMobile_H_Revrs,
                                       nCanon2AtnoRevrs, pVA, &pStruct->ti, &EndpList );
        if ( nNumWrongTg < 1 ) {
            goto exit_case_06e; /* for now only transfer (+) from one Mobile-H group to another */
        }
        nNum2RemovePlus = nNum2AddPlus = nNum2MovePlus = 0;
        for ( i = 0; i < nNumWrongTg; i ++ ) {
            /* detect t-group that has extra (+) on H */
            if ( tdhc[i].nNumHInchi > tdhc[i].nNumHNorml &&
                 tdhc[i].nNumPRevrs > tdhc[i].nNumPNorml && tdhc[i].n[fNumRPosChgH] ) {
                /* count how many (+) to remove */
                /* store XH(+) atom numbers */
                int nNumNeeded = inchi_min( tdhc[i].nNumHInchi-tdhc[i].nNumHNorml, tdhc[i].n[fNumRPosChgH]);
                nNum2RemovePlus += nNumNeeded;
                jjoffs = tdhc[i].i[ fNumRPosChgH ];
                for ( jj = 0; jj < tdhc[i].n[fNumRPosChgH]; jj ++ ) {
                    iat = EndpList.pnEdges[ jjoffs + jj ];
                    e = pVA[iat].nCPlusGroupEdge-1;
                    if ( ret = AddToEdgeList( &CurChargeEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_case_06e;
                    }
                }
            } else
            /* detect t-group that needs (+) on XH to reduce number of H */
            if ( tdhc[i].nNumHInchi < tdhc[i].nNumHNorml && tdhc[i].n[fNumRNeutrlH] ) {
                /* store XH atom numbers */
                int nNumNeeded = inchi_min( tdhc[i].nNumHNorml-tdhc[i].nNumHInchi, tdhc[i].n[fNumRNeutrlH]);
                nNum2AddPlus += nNumNeeded;
                jjoffs = tdhc[i].i[ fNumRNeutrlH ];
                for ( jj = 0; jj < tdhc[i].n[fNumRNeutrlH]; jj ++ ) {
                    iat = EndpList.pnEdges[ jjoffs + jj ];
                    e = pVA[iat].nCPlusGroupEdge-1;
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_case_06e;
                    }
                }
            }
        }
        nNum2MovePlus = inchi_min( nNum2RemovePlus, nNum2AddPlus );
        if ( CurrEdges.num_edges > 0 && CurChargeEdges.num_edges > 0 ) {
            for ( i = 0; 0 < nNum2MovePlus && i < nNumWrongTg; i ++ ) {
                /* detect t-group that has extra (+) on H */
                if ( tdhc[i].nNumHInchi > tdhc[i].nNumHNorml &&
                     tdhc[i].nNumPRevrs > tdhc[i].nNumPNorml && tdhc[i].n[fNumRPosChgH] ) {
                    int nNum2Remove = tdhc[i].nNumHInchi - tdhc[i].nNumHNorml;
                    if ( nNum2Remove < tdhc[i].n[fNumRPosChgH] ) {
                        nNum2Remove = tdhc[i].n[fNumRPosChgH];
                    }
                    /* store XH(+) atom numbers */
                    jjoffs = tdhc[i].i[ fNumRPosChgH ];
                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                    RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
                    for ( jj = 0; 0 < nNum2MovePlus && 0 < nNum2Remove && jj < tdhc[i].n[fNumRPosChgH]; jj ++ ) {
                        iat = EndpList.pnEdges[ jjoffs + jj ];
                        e = pVA[iat].nCPlusGroupEdge-1;
                        pe   = pBNS->edge + pVA[iat].nCPlusGroupEdge-1;
                        if ( pe->flow )
                            continue;
                        pv1 = pBNS->vert + (v1 = pe->neighbor1);
                        pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                        for ( j = pv1->num_adj_edges-1; 0 <= j; j -- ) {
                            pe1n = pBNS->edge + pv1->iedge[j];
                            if ( pe1n->flow && !pe1n->forbidden ) {
                                pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                                break;
                            }
                        }
                        if ( j < 0 )
                            continue; /* not found */
                        
                        for ( j = pv2->num_adj_edges-2; 0 <= j; j -- ) {
                            pe2n = pBNS->edge + pv2->iedge[j];
                            if ( pe2n->flow && !pe2n->forbidden ) {
                                pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                                break;
                            }
                        }
                        if ( j < 0 )
                            continue; /* not found */
                        delta = 1;
                        pe->flow   += delta;
                        pe1n->flow -= delta;
                        pe2n->flow -= delta;
                        pv1n->st_edge.flow -= delta;
                        pv2n->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2*delta;

                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                        if ( ret == 1 && (vPathEnd == v1n && vPathStart == v2n ||
                                          vPathEnd == v2n && vPathStart == v1n) &&
                                          (nDeltaCharge == 0 || nDeltaCharge == 1) ) {
                            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                            if ( ret > 0 ) {
                                nNumRunBNS ++;
                                nNum2Remove --;
                                nNum2MovePlus --;
                                cur_success ++; /* 06e */
                            }
                        } else {
                            pe->flow   -= delta;
                            pe1n->flow += delta;
                            pe2n->flow += delta;
                            pv1n->st_edge.flow += delta;
                            pv2n->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2*delta;
                        }
                        if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                            goto exit_case_06e;
                        }
                    }
                    RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                }
            }
        }
exit_case_06e:
        CurrEdges.num_edges = 0; /* clear current edge list */
        AllocEdgeList( &CurChargeEdges, EDGE_LIST_FREE );
        AllocEdgeList( &EndpList, EDGE_LIST_FREE );
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }



    if ( pc2i->len_c2at >= 1 ) {
        /*--------------------------------------------------------------*/
        /* case 07: restored:  O(-)-AB=O  original:  O=AB-O(-)          */
        /* FixH:               0       0             0     -1           */
        /* MobH:               0       0             0      1           */
        /*                    taut  (non-taut)     (taut) non-taut      */
        /*                    taut  (taut)     (non-taut) non-taut      */   
        /* O = O, S, Se, Te                                             */
        /* Solution: move (-) from O(-)-AB to AB=O                      */
        /*--------------------------------------------------------------*/
        int num_SB_O_Minus = 0, num_DB_O_Neutr = 0, iat;
        short iat_SB_O_Minus[MAX_DIFF_FIXH], iat_DB_O_Neutr[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: -O(-), non-taut */
                 num_DB_O_Neutr < MAX_DIFF_FIXH &&
                 pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                 !pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI ==  1 &&
                 /* reversed structure info: */
                 pc2i->c2at[i].nFixHRevrs ==  0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                 at2[iat].valence < at2[iat].chem_bonds_valence ) {
                iat_DB_O_Neutr[num_DB_O_Neutr ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom: charge=-1, no H, has single bond, O, S, Se, Te */
                 num_SB_O_Minus < MAX_DIFF_FIXH &&
                 at2[iat].charge == -1 && !at2[iat].num_H &&
                 at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 pVA[iat].cNumValenceElectrons == 6 &&
                 at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint &&
                 /* in orig.InChI: not an endpoint, has no H */
                 /*pStruct->endpoint[i] && -- modificatuion#1 */
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has (-) edge */
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                iat_SB_O_Minus[num_SB_O_Minus ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_SB_O_Minus, num_DB_O_Neutr ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_O_Minus && cur_success < num_try; i ++ ) {
                iat = iat_SB_O_Minus[i];
                pe   = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Moved (-) charge to AB=O => nDeltaCharge == 1 */
                    /* Flow change on pe (-)charge edge (O(-)-AB) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 07 */
                    }
                } else {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->len_c2at >= 1 ) {
        /*--------------------------------------------------------------*/
        /* case 07a: restored: O(-)-N(V)B=O  original:  O=N(V)B-O(-)    */
        /* FixH:               0          0             0      -1       */
        /* MobH:               0          0             0       1       */
        /*                non-taut  (non-taut)  non-taut  non-taut      */
        /*                non-taut  (taut)      non-taut  non-taut      */   
        /* O = O, S, Se, Te                                             */
        /* Solution: move (-) from O(-)-AB to AB=O                      */
        /*--------------------------------------------------------------*/
        int num_SB_O_Minus = 0, num_DB_O_Neutr = 0, iat, iN;
        short iat_SB_O_Minus[MAX_DIFF_FIXH], iat_DB_O_Neutr[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: -O(-), non-taut */
                 num_DB_O_Neutr < MAX_DIFF_FIXH &&
                 pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                 !pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == -1 && pc2i->c2at[i].nMobHInChI ==  1 &&
                 /* reversed structure info: */
                 pc2i->c2at[i].nFixHRevrs ==  0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                 at2[iat].valence < at2[iat].chem_bonds_valence ) {
                iat_DB_O_Neutr[num_DB_O_Neutr ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom: charge=-1, no H, has single bond, O, S, Se, Te */
                 num_SB_O_Minus < MAX_DIFF_FIXH &&
                 at2[iat].charge == -1 && !at2[iat].num_H &&
                 at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 pVA[iat].cNumValenceElectrons == 6 &&
                 /*at_Mobile_H_Revrs && !at_Mobile_H_Revrs[iat].endpoint &&*/
                 /* in orig.InChI: not an endpoint, has no H */
                 !pStruct->endpoint[i] &&
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has N(V) neighbor */
                 1 == at2[iat].valence && at2[iN=at2[iat].neighbor[0]].chem_bonds_valence==5 &&
                 !at2[iN].charge && pVA[iN].cNumValenceElectrons == 5 &&
                 /* has (-) edge */
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                iat_SB_O_Minus[num_SB_O_Minus ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_SB_O_Minus, num_DB_O_Neutr ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_O_Minus && cur_success < num_try; i ++ ) {
                iat = iat_SB_O_Minus[i];
                pe   = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Moved (-) charge to AB=O => nDeltaCharge == 1 */
                    /* Flow change on pe (-)charge edge (O(-)-AB) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 07 */
                    }
                } else {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }
    if ( /*(pc2i->len_c2at >= 1 || pc2i->nNumRemHRevrs) &&*/ pc2i->nNumTgInChI == 1 && /* ADP in InChI */
         (pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1) ) {
        /*----------------------------------------------------------------*/
        /* case 08: restored: O(-)-AB=N- OH- orig.   O=AB-N(-)- OH-       */
        /* FixH:               1      0   0          0      0   1         */
        /* MobH:               0      0   1          0      0   0         */
        /*           may be taut or not  non-taut   taut  taut taut       */
        /*                                    ADP: one t-group or more endpoints */
        /* O(-) = S, Se, Te; N = N;                                       */
        /* Solution: move (-) from O(-) to =N-; avoid stereogenic DB on N */
        /*----------------------------------------------------------------*/
        int num_DB_N_Neutr = 0, num_SB_O_Minus = 0, iat;
        short iat_DB_N_Neutr[MAX_DIFF_FIXH], iat_SB_O_Minus[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        cur_success = 0;
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom: charge=-1, has no H, has single bond, O, S, Se, Te */
                 num_SB_O_Minus < MAX_DIFF_FIXH &&
                 at2[iat].charge == -1 && !at2[iat].num_H &&
                 at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 pVA[iat].cNumValenceElectrons == 6 &&
                 /* in orig.InChI: an endpoint, may have fixed-H */
                 pStruct->endpoint[i] && 
                 /*!(pStruct->fixed_H && pStruct->fixed_H[i]) &&*/
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has (-) edge */
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {

                iat_SB_O_Minus[num_SB_O_Minus ++] = iat;
                
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            } else
            if ( /* in restored atom: charge=0, has no H, has double non-stereogenic bond, N */
                 num_DB_N_Neutr < MAX_DIFF_FIXH &&
                 at2[iat].charge == 0 && !at2[iat].num_H && !at2[iat].sb_parity[0] &&
                 at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                 /* in orig.InChI: an endpoint, has no fixed-H */
                 pStruct->endpoint[i] &&
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has (-) edge */
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 &&
                 0 == pBNS->edge[e].forbidden ) {

                iat_DB_N_Neutr[num_DB_N_Neutr ++] = iat;
                
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_DB_N_Neutr, num_SB_O_Minus ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            /* allow stereobonds in rings change */
            if ( forbidden_stereo_edge_mask )
                RemoveForbiddenEdgeMask( pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask );

            delta = 1;
            for ( i = 0; i < num_SB_O_Minus && cur_success < num_try; i ++ ) {
                iat = iat_SB_O_Minus[i];
                pe   = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Moved (-) charge to =N- => nDeltaCharge == 1 */
                    /* Flow change on pe (-)charge edge (atom (-)O-) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 08 */
                    }
                } else {
                    pe->forbidden &= forbidden_edge_mask_inv;
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            if ( forbidden_stereo_edge_mask )
                SetForbiddenEdgeMask( pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->len_c2at >= 2 ) {
        /*--------------------------------------------------------*/
        /* case 09: restored: NH2(+)=C--NH2 orig.: NH2-C(+)-NH2   */
        /* FixH:               2     |  2            0  |   0     */
        /* MobH:               0        0            2      2     */
        /* N = N,            taut      taut     non-taut  non-taut*/
        /* Solution: move (+) from NH2(+) to C                    */
        /*--------------------------------------------------------*/
        int iat;
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( (pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1) &&
                 /* orig. InChI info: */
                 !pc2i->c2at[i].endptInChI &&
                 pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI &&
                 /* reversed structure info: */
                 pc2i->c2at[i].endptRevrs &&
                 pc2i->c2at[i].nFixHRevrs && !pc2i->c2at[i].nMobHRevrs &&
                 pc2i->c2at[i].nAtChargeRevrs == 1 &&
                 at2[iat].valence + 1 == at2[iat].chem_bonds_valence &&
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {
                EdgeIndex eNC = NO_VERTEX, eCPlusC;
                int       iNH2, iatC, iatNH2, icNH2;
                /* found NH2(+)=; locate =C< and find whether it has -NH2 neighbor */
                for ( j = 0; j < at2[iat].valence; j ++ ) {
                    if ( at2[iat].bond_type[j] == BOND_TYPE_DOUBLE )
                        break;
                }
                if ( j == at2[iat].valence )
                    continue;
                eNC = pBNS->vert[iat].iedge[j]; /* edge NH2(+)=C */
                iatC = at2[iat].neighbor[j];
                if ( pVA[iatC].cNumValenceElectrons != 4 || pVA[iatC].cMetal || at2[iatC].charge ||
                     at2[iatC].valence != 3 || at2[iatC].valence+1 != at2[iatC].chem_bonds_valence ||
                     (eCPlusC=pVA[iatC].nCPlusGroupEdge-1) < 0 || pBNS->edge[eCPlusC].forbidden)
                    continue;
                for ( j = 0; j < at2[iatC].valence; j ++ ) {
                    iatNH2 = at2[iatC].neighbor[j];
                    if ( iatNH2 == iat || pVA[iatNH2].cNumValenceElectrons != 5 ||
                         pVA[iatNH2].cPeriodicRowNumber != 1 || !at2[iatNH2].num_H || at2[iatNH2].charge)
                        continue;
                    icNH2 = pStruct->nAtno2Canon[0][iatNH2];
                    for ( iNH2 = 0; iNH2 < pc2i->len_c2at; iNH2 ++ ) {
                        if ( iatNH2 == pc2i->c2at[iNH2].atomNumber )
                            break;
                    }
                    if ( iNH2 == pc2i->len_c2at )
                        continue;

                    if ( (pc2i->c2at[iNH2].nValElectr == 5 && pc2i->c2at[iNH2].nPeriodNum == 1) &&
                         /* orig. InChI info: */
                         !pc2i->c2at[iNH2].endptInChI &&
                         pc2i->c2at[iNH2].nFixHInChI == 0 && pc2i->c2at[iNH2].nMobHInChI &&
                         /* reversed structure info: */
                         pc2i->c2at[iNH2].endptRevrs &&
                         pc2i->c2at[iNH2].nFixHRevrs && !pc2i->c2at[iNH2].nMobHRevrs &&
                         pc2i->c2at[iNH2].nAtChargeRevrs == 0 &&
                         at2[iatNH2].valence == at2[iatNH2].chem_bonds_valence ) {
                        /* we have found NH2(+)=, =C<, and bond between them */
 
                        if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                        if ( ret = AddToEdgeList( &CurrEdges, eCPlusC, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                        SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                        RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
                        delta = 1;

                        pe   = pBNS->edge + eNC;
                        if ( !pe->flow )
                            continue;
                        pv1 = pBNS->vert + (v1 = pe->neighbor1);
                        pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                        pe->forbidden |= forbidden_edge_mask;
                        pe->flow -= delta; /* add (+) to -NHm */
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2*delta;

                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                        if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                          vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 0 ) {
                            /* Removed (+)charge from -NHn => nDeltaCharge == -1 */
                            /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                            if ( ret > 0 ) {
                                nNumRunBNS ++;
                                cur_success ++; /* 09 */
                            }
                        } else {
                            pe->flow += delta;
                            pv1->st_edge.flow += delta;
                            pv2->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2*delta;
                        }
                        INCHI_HEAPCHK
                        pe->forbidden &= forbidden_edge_mask_inv;
                        RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                        CurrEdges.num_edges = 0; /* clear current edge list */
                        break;
                    }
                }
            }
        }
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }


    if ( pc2i->len_c2at >= 2 ) {
        /*--------------------------------------------------------*/
        /* case 10: restored: NH2-X(+)-NH-  orig.: NH2(+)=X-NH-   */
        /* FixH:               0        0            2      1     */
        /* MobH:               2        1            0      0     */
        /* N = N,O,S,Se,Te non-taut  non-taut       taut   taut   */
        /* Solution: move (+) from X(+) to NH2 or NH              */
        /*--------------------------------------------------------*/
        int iat;
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            if ( pc2i->c2at[i].nValue )
                continue;
            iat = pc2i->c2at[i].atomNumber;
            if ( (pc2i->c2at[i].nValElectr == 6 ||
                  pc2i->c2at[i].nValElectr == 5 && pc2i->c2at[i].nPeriodNum == 1) &&
                 /* orig. InChI info: */
                 pc2i->c2at[i].endptInChI &&
                 pc2i->c2at[i].nFixHInChI && !pc2i->c2at[i].nMobHInChI &&
                 /* reversed structure info: */
                 !pc2i->c2at[i].endptRevrs &&
                 !pc2i->c2at[i].nFixHRevrs && pc2i->c2at[i].nMobHRevrs &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 &&
                 at2[iat].valence == at2[iat].chem_bonds_valence &&
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {

                EdgeIndex eCPlusC, eCPlusNH2, bContinue=1;
                int       iNH2, iatC, iatNH2, icNH2, j1, j2;
                BNS_EDGE *pe_iat, *pe_iNH2;
                /* found NH2- locate -X(+) and find whether it has another -NH2 neighbor */
                for ( j1 = 0; j1 < at2[iat].valence && bContinue; j1 ++ ) {
                    if ( at2[iat].bond_type[j1] == BOND_TYPE_SINGLE &&
                         at2[iatC = at2[iat].neighbor[j1]].charge == 1 &&
                         (4 <= pVA[iatC].cNumValenceElectrons && pVA[iatC].cNumValenceElectrons <= 6) &&
                         at2[iatC].valence == at2[iatC].chem_bonds_valence &&
                         (eCPlusC=pVA[iatC].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[eCPlusC].forbidden) {
                        /* found a candidate for X; find another NH2 */
                        for ( j2 = 0; j2 < at2[iatC].valence && bContinue; j2 ++ ) {
                            if ( at2[iatC].bond_type[j2] == BOND_TYPE_SINGLE &&
                                 iat != (iatNH2 = at2[iatC].neighbor[j2]) &&
                                 at2[iatNH2].charge == 0 && at2[iatNH2].num_H &&
                                 (pVA[iatNH2].cNumValenceElectrons==5 || pVA[iatNH2].cNumValenceElectrons==6) &&
                                 at2[iatNH2].valence == at2[iatNH2].chem_bonds_valence &&
                                 (eCPlusNH2=pVA[iatNH2].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[eCPlusNH2].forbidden) {
                                for ( iNH2 = 0; iNH2 < pc2i->len_c2at; iNH2 ++ ) {
                                    if ( iatNH2 != pc2i->c2at[iNH2].atomNumber || pc2i->c2at[iNH2].nValue )
                                        continue;
                                    /* check the second -NH */
                                    icNH2 = pStruct->nAtno2Canon[0][iatNH2]; /* canon number -1 */
                                    if ( /* orig. InChI info: */
                                         pc2i->c2at[iNH2].endptInChI &&
                                         pc2i->c2at[iNH2].nFixHInChI && !pc2i->c2at[iNH2].nMobHInChI &&
                                         /* reversed structure info: */
                                         !pc2i->c2at[iNH2].endptRevrs &&
                                         !pc2i->c2at[iNH2].nFixHRevrs && pc2i->c2at[iNH2].nMobHRevrs &&
                                         pc2i->c2at[iNH2].nAtChargeRevrs == 0 ) {
                                        /* we have found NH-X(+)-NH; remove charge from X(+) */
                                        pe_iat  = pBNS->edge + pBNS->vert[iat].iedge[j1];
                                        pe_iNH2 = pBNS->edge + pBNS->vert[iatC].iedge[j2];
                                        /* pick up one of -NH to move (+) to it */
                                        if ( !pe_iat->forbidden && pBNS->edge[e].flow ) {
                                            pe = pBNS->edge + e;
                                        } else
                                        if ( !pe_iNH2->forbidden && pBNS->edge[eCPlusNH2].flow ) {
                                            pe = pBNS->edge + eCPlusNH2;
                                        } else {
                                            continue; /* none of the two -X(+)- bonds may be changed */
                                        }
                                        if ( ret = AddToEdgeList( &CurrEdges, eCPlusC, INC_ADD_EDGE ) ) {
                                            goto exit_function;
                                        }
                                        SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                                        RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
                                        delta = 1;

                                        pv1 = pBNS->vert + (v1 = pe->neighbor1);
                                        pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                                        /*pe->forbidden |= forbidden_edge_mask;*/
                                        pe->flow -= delta; /* add (+) to -NHm */
                                        pv1->st_edge.flow -= delta;
                                        pv2->st_edge.flow -= delta;
                                        pBNS->tot_st_flow -= 2*delta;

                                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                                        if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                                          vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == -1 ) {
                                            /* Removed (+)charge from -NHn => nDeltaCharge == -1 */
                                            /* Flow change on pe (+)charge edge (atom NHm(+)) is not known to RunBnsTestOnce()) */
                                            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                                            if ( ret > 0 ) {
                                                nNumRunBNS ++;
                                                cur_success ++; /* 10 */
                                                bContinue = 0;
                                                pc2i->c2at[i].nValue    = 1;    /* mark as used */
                                                pc2i->c2at[iNH2].nValue = 1; /* mark as used */
                                            }
                                        } else {
                                            pe->flow += delta;
                                            pv1->st_edge.flow += delta;
                                            pv2->st_edge.flow += delta;
                                            pBNS->tot_st_flow += 2*delta;
                                        }
                                        INCHI_HEAPCHK

                                        /*pe->forbidden &= forbidden_edge_mask_inv;*/
                                        RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                                        CurrEdges.num_edges = 0; /* clear current edge list */
                                        break;
                                    }
                                } /* iNH2: pc2i->c2at[iNH2] cycle */
                            }
                        } /* j2: iatC neighbors cycle */
                    }
                } /* j1: iat neighbors cycle */
            }
        } /* i: pc2i->c2at[i] cycle */
        if ( cur_success ) {
            /*
            for ( i = 0; i < pc2i->len_c2at; i ++ ) {
                pc2i->c2at[i].nValue = 0;
            }
            */
            tot_succes += cur_success;
            /* recalculate InChI from the structure */
            if ( 0 > (ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                            ppt_group_info, ppat_norm, ppat_prep ) ) ) {
                goto exit_function;
            }
            if ( ret = FillOutExtraFixedHDataRestr( pStruct ) ) {
                goto exit_function;
            }
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( /*pc2i->len_c2at >= 1 &&*/ pc2i->nNumTgInChI == 1 && /* ADP in InChI */
         (pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI || pc2i->nNumTgRevrs > 1) ) {
        /*--------------------------------------------------------------*/
        /* case 11: restored: NH(+)=AB-N< OH- orig.  NH-AB=N(+)< OH-    */
        /* FixH:               0       0   0          1      0   1      */
        /* MobH:               1       0   1          0      0   0      */
        /*                 non-taut.                taut        taut    */
        /*                                    ADP: one t-group or more endpoints */
        /* NH(+)= => N, O, S, Se; -N< => N                              */
        /* Solution: move (+) from NH(+) to -N<                         */
        /*--------------------------------------------------------------*/
        int num_SB_Neutr = 0, num_DB_Charged = 0, iat;
        short iat_SB_Neutr[MAX_DIFF_FIXH], iat_DB_Charged[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        cur_success = 0;
        /* search for NH(+)= */
        /* search for -N< */
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom: charge=0, has no H, has no double bond, N only */
                 num_DB_Charged < MAX_DIFF_FIXH &&
                 at2[iat].charge == 1 && at2[iat].num_H &&
                 at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 (pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 ||
                  pVA[iat].cNumValenceElectrons == 6 ) &&
                 /* in orig.InChI: an endpoint, has fixed-H */
                 /*pStruct->endpoint[i] &&*/
                 (pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 /*!(nMobHInChI && nMobHInChI[i] ) &&*/
                 /* has (+) edge */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && 0 == pBNS->edge[e].forbidden ) {

                iat_DB_Charged[num_DB_Charged ++] = iat;
                /*
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
                */
            } else
            if ( /* in restored atom: charge=0, has no H, has no double bond, N only */
                 num_SB_Neutr < MAX_DIFF_FIXH &&
                 at2[iat].charge == 0 && !at2[iat].num_H &&
                 at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 (pVA[iat].cNumValenceElectrons == 5 &&
                  pVA[iat].cPeriodicRowNumber == 1 ) &&
                 /* in orig.InChI: an endpoint, has fixed-H */
                 /*pStruct->endpoint[i] &&*/
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has (+) edge */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && 0 == pBNS->edge[e].forbidden ) {

                iat_SB_Neutr[num_SB_Neutr ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_SB_Neutr, num_DB_Charged ) ) {
            /* detected; attempt to fix */
            BNS_VERTEX *pv1n, *pv2n;
            BNS_EDGE   *pe1n, *pe2n;
            Vertex      v1n, v2n;
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_DB_Charged && cur_success < num_try; i ++ ) {
                iat = iat_DB_Charged[i];
                pe   = pBNS->edge + pVA[iat].nCPlusGroupEdge-1;
                if ( pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                for ( j = pv1->num_adj_edges-1; 0 <= j; j -- ) {
                    pe1n = pBNS->edge + pv1->iedge[j];
                    if ( pe1n->flow && !pe1n->forbidden ) {
                        pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                        break;
                    }
                }
                if ( j < 0 )
                    continue; /* not found */
                
                for ( j = pv2->num_adj_edges-2; 0 <= j; j -- ) {
                    pe2n = pBNS->edge + pv2->iedge[j];
                    if ( pe2n->flow && !pe2n->forbidden ) {
                        pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                        break;
                    }
                }
                if ( j < 0 )
                    continue; /* not found */

                pe->flow   += delta;
                pe1n->flow -= delta;
                pe2n->flow -= delta;
                pv1n->st_edge.flow -= delta;
                pv2n->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1n && vPathStart == v2n ||
                                  vPathEnd == v2n && vPathStart == v1n) &&
                                  (nDeltaCharge == 0 || nDeltaCharge == 1) ) {
                    /* before setting flows the structure could be:
                       [NH+ neigh, v1n]=e1n=[NH+,v1]-pe-[+,v2]=e2n=[another at or its chargeStruct]
                       or
                        
                         [NH+ or ChStr, v1n]=pe1n=[NH+ or ChStr, v1]-pe-[+,v2]=pe2n=[at2 or ChStr, v2n]
                                                                     ^    ^    ^
                                                               NH+(+)edge |  N (+) edge: only
                                                                          |  these are not forbidden
                                                                          |
                                                                   hetero (+) vertex

                        After setting flows (* mark radicals, =pe= is forbidden):
                       
                       *[NH+ or ChStr, v1n]-pe1n-[NH+ or ChStr, v1]=pe=[+,v2]-pe2n-[at2 or ChStr, v2n]*
                                                                     ^    ^    ^
                                                               NH+(+)edge |  N (+) edge: only
                                                                          |  these are not forbidden
                                                                          |
                                                                   hetero (+) vertex

                        Flow in
                        pe1n and pe2n will or will not change, depending on the structure.

                        Consider what happens if pe2n changes. It may only increment.
                        If pe2n flow increments then another (+)edge flow dectrements. If
                        [at2 or ChStr, v2n] is at2 then at2 charge would change from (+) to 0,
                        and another N charge would change from 0 to (+), giving tot. change of
                        number of charges  (-1)+(+1)=0. However, if [at2 or ChStr, v2n] is
                        ChargeStruct then at2 will not be on the alt path and only the creation
                        of another (+) will be detected.
                    */
                    /* Removed charge from O(+) => nDeltaCharge == -1 */
                    /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 11 */
                    }
                } else {
                    pe->flow   -= delta;
                    pe1n->flow += delta;
                    pe2n->flow += delta;
                    pv1n->st_edge.flow += delta;
                    pv2n->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->len_c2at >= 1 && pc2i->nNumTgInChI == 1 &&
         pc2i->nNumRemHInChI >= -1 && /* 2006-03-03 */
         (pc2i->nNumEndpInChI > pc2i->nNumEndpRevrs || pc2i->nNumTgRevrs > 1) /* ADP in InChI */ ) {
        /*--------------------------------------------------------------*/
        /* case 12: restored: O=AB-N<         original: (-)O-AB=N(+)<   */
        /* FixH:              0    0                     0        0     */
        /* MobH:              0    0                     0        0     */
        /*                   non-taut                   taut            */
        /* O = O, S, Se, N; N = N;                                         */
        /* restored atom O is not tautomeric; original atom O is taut.  */
        /* original struct has 1 t-group; restored has less endpoints   */
        /*                             and/or possibly >1 t-groups      */  
        /* Solution: separate charges between O= and -N<                */
        /*           allow moving charge to N(V) to make it N(IV)(+)    */
        /*--------------------------------------------------------------*/
        int bOnly_N_V = 1;
        cur_success = 0;
        while( 1 ) {
            int num_SB_N_Neutr = 0, num_DB_O = 0, iat, num_N_V=0, bN_V;
            short iat_SB_N_Neutr[MAX_DIFF_FIXH], iat_DB_O[MAX_DIFF_FIXH];
            AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
            inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                                 pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
            S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                                   pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
            cur_success = 0;
            for ( i = 0; i < pc2i->len_c2at; i ++ ) {
                iat = pc2i->c2at[i].atomNumber;
                if ( /* orig. InChI info: -O(-) */
                     num_DB_O < MAX_DIFF_FIXH &&
                     (pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ ||
                      pc2i->c2at[i].nValElectr == 5 &&
                      pc2i->c2at[i].nPeriodNum == 1 /* N */ ) &&
                     pc2i->c2at[i].endptInChI &&
                     (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                     pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                     /* reversed structure info: */
                     !pc2i->c2at[i].endptRevrs &&
                     pc2i->c2at[i].nFixHRevrs ==  0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                     pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                     ((pc2i->c2at[i].nValElectr == 6)? 
                           (at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2):
                      (pc2i->c2at[i].nValElectr == 5)?
                           (at2[iat].valence == 2 && at2[iat].chem_bonds_valence == 3):
                           0)) {
                
                    iat_DB_O[num_DB_O ++] = iat;
                    /*
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                    */
                }
            }
            for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                bN_V = 0;
                if ( /* in restored atom N: charge=0, no H, has no double bond, not an endpoint */
                     num_SB_N_Neutr < MAX_DIFF_FIXH &&
                     at2[iat].charge == 0 && !at2[iat].num_H &&
                     (at2[iat].valence == at2[iat].chem_bonds_valence ||
                     (bN_V = at2[iat].valence+2 == at2[iat].chem_bonds_valence)) &&
                     !pVA[iat].cMetal &&
                     pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                     !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint) &&
                     /* in orig.InChI: not an endpoint, has no H */
                     !pStruct->endpoint[i] &&
                     !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                     !(nMobHInChI && nMobHInChI[i]) &&
                     /* has (+) edge */
                     (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {

                    if ( bOnly_N_V && bN_V &&
                         NO_VERTEX != (j = GetChargeFlowerUpperEdge( pBNS, pVA, e )) &&
                         !pBNS->edge[j].forbidden && !pBNS->edge[j].flow ) {
                        if ( !num_N_V ) {
                            /* switch to N(V) only mode */
                            CurrEdges.num_edges = 0;
                            num_SB_N_Neutr = 0;
                        }
                        iat_SB_N_Neutr[num_SB_N_Neutr ++] = iat;
                        num_N_V ++;
                        if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                        if ( ret = AddToEdgeList( &CurrEdges, j, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                    } else
                    if ( !num_N_V ) {
                        iat_SB_N_Neutr[num_SB_N_Neutr ++] = iat;
                        if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                        /* in addition, permit N(V)=>N(IV)(+) change by allowing charge flower edge change flow */
                        if ( bN_V && NO_VERTEX != (j = GetChargeFlowerUpperEdge( pBNS, pVA, e )) &&
                             !pBNS->edge[j].forbidden && !pBNS->edge[j].flow ) {
                            if ( ret = AddToEdgeList( &CurrEdges, j, INC_ADD_EDGE ) ) {
                                goto exit_function;
                            }
                        }
                    }
                }
            }
            if ( num_try = inchi_min( num_SB_N_Neutr, num_DB_O ) ) {
                /* detected; attempt to fix */
                BNS_EDGE *pe_CMinus;
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
                delta = 1;
                for ( i = 0; i < num_DB_O && cur_success < num_try; i ++ ) {
                    iat = iat_DB_O[i];
                    pe_CMinus = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                    pe_CMinus->forbidden &= forbidden_edge_mask_inv;

                    pe   = pBNS->edge + pBNS->vert[iat].iedge[0]; /* double bond O=...*/
                    if ( !pe->flow )
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->forbidden |= forbidden_edge_mask; /* change bond O=X to O(rad)-X(rad) */
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2*delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 2 ) {
                        /* Added (-) charge to =O and (+) charge to N => nDeltaCharge == 2 */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if ( ret > 0 ) {
                            nNumRunBNS ++;
                            cur_success ++; /* 12 */
                        }
                    } else {
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2*delta;
                    }
                    pe->forbidden        &= forbidden_edge_mask_inv; /* allow changes to O=X bond */
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
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
                if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                    goto exit_function;  /* no fixed-H found */
                }
                if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                    goto exit_function;
                }
                if ( !pc2i->bHasDifference ) {
                    goto exit_function; /* nothing to do */
                }
                break;
            } else
            if ( bOnly_N_V ) {
                bOnly_N_V = 0;
            } else {
                break;
            }
        }
    }

    if ( pc2i->nNumTgDiffMinus /*|| pc2i->nNumTgDiffH */ /* no ADP in InChI needed */ ) {
        /*--------------------------------------------------------------*/
        /*                         |                            |       */
        /* case 13: restored: O=AB=N=         original: (-)O-AB-N(+)=   */
        /* FixH:              0    0                     0        0     */
        /* MobH:              0    0                     0        0     */
        /*                        non-taut              taut   non-taut */
        /* O = O, S, Se, N; N = N, P, ...                               */
        /* t-group in original has same num. endpoints                  */
        /*       same num_H and less (-) than in the restored structure */
        /* original atom O is tautomeric, N is not taut in both         */
        /* original struct has 1 t-group; restored has less endpoints   */
        /*                             and/or possibly >1 t-groups      */  
        /* Solution: separate charges between O= and -N<                */
        /*           allow moving charge to N(V) to make it N(IV)(+)    */
        /*--------------------------------------------------------------*/
        int itg;
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;

        S_CHAR   *num_Fixed_H_Revrs = pStruct->pOneINChI[0]->nNum_H_fixed? pStruct->pOneINChI[0]->nNum_H_fixed : NULL;
        S_CHAR   *pnMobHRevrs = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H)?
                                   pStruct->pOneINChI[1]->nNum_H : 
                                (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H)?
                                   pStruct->pOneINChI[0]->nNum_H : NULL;
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        cur_success = 0;
        /* find whether this may help */
        for ( itg = 0; itg < pStruct->ti.num_t_groups && itg < pStruct->One_ti.num_t_groups; itg ++ ) {
            if ( pStruct->ti.t_group[itg].nNumEndpoints == pStruct->One_ti.t_group[itg].nNumEndpoints &&
                 pStruct->ti.t_group[itg].num[0] - pStruct->ti.t_group[itg].num[1] == 
                 pStruct->One_ti.t_group[itg].num[0] - pStruct->One_ti.t_group[itg].num[1] &&
                 pStruct->ti.t_group[itg].num[1] > pStruct->One_ti.t_group[itg].num[1]) {
                /* restored InChI t-group has more (-) and same number of H */

                int num_SB_N_Neutr = 0, num_DB_O = 0, iat;
                short iat_SB_N_Neutr[MAX_DIFF_FIXH], iat_DB_O[MAX_DIFF_FIXH];
                cur_success = 0;
                for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
                    iat = nCanon2AtnoRevrs[i];
                    if ( /* orig. InChI info: -O(-) */
                         num_DB_O < MAX_DIFF_FIXH &&
                         (pVA[i].cNumValenceElectrons == 6 /* O, S, Se, Te */ ) &&
                         pStruct->endpoint[i] == itg+1 &&
                         (e=pVA[i].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                         !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                         !(nMobHInChI && nMobHInChI[i]) &&
                         /* reversed structure info: */
                         /*!pc2i->c2at[i].endptRevrs &&*/
                         !(num_Fixed_H_Revrs && num_Fixed_H_Revrs[iat]) &&
                         !(pnMobHRevrs && pnMobHRevrs[iat]) &&
                         at2[iat].charge == 0 && at2[iat].num_H == 0 &&
                         at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 ) {
            
                        iat_DB_O[num_DB_O ++] = iat;
                        /*
                        if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                        */
                    } else
                    if ( /* in restored atom N: charge=0, no H, has no double bond, not an endpoint */
                         num_SB_N_Neutr < MAX_DIFF_FIXH &&
                         at2[iat].charge == 0 && !at2[iat].num_H &&
                         /*at2[iat].valence == at2[iat].chem_bonds_valence ||*/
                         (at2[iat].valence==4 && at2[iat].chem_bonds_valence==5) &&
                         !pVA[iat].cMetal &&
                         pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber >= 1 &&
                         !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint) &&
                         /* in orig.InChI: not an endpoint, has no H */
                         !pStruct->endpoint[i] &&
                         !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                         !(nMobHInChI && nMobHInChI[i]) &&
                         /* has (+) edge */
                         (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {

                        iat_SB_N_Neutr[num_SB_N_Neutr ++] = iat;
                        if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                            goto exit_function;
                        }
                    }
                }
                if ( num_try = inchi_min( num_SB_N_Neutr, num_DB_O ) ) {
                    /* detected; attempt to fix */
                    BNS_EDGE *pe_CMinus;
                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
                    delta = 1;
                    for ( i = 0; i < num_DB_O && cur_success < num_try; i ++ ) {
                        iat = iat_DB_O[i];
                        pe_CMinus = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                        pe_CMinus->forbidden &= forbidden_edge_mask_inv;

                        pe   = pBNS->edge + pBNS->vert[iat].iedge[0]; /* double bond O=...*/
                        if ( !pe->flow )
                            continue;
                        pv1 = pBNS->vert + (v1 = pe->neighbor1);
                        pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                        pe->forbidden |= forbidden_edge_mask; /* change bond O=X to O(rad)-X(rad) */
                        pe->flow -= delta;
                        pv1->st_edge.flow -= delta;
                        pv2->st_edge.flow -= delta;
                        pBNS->tot_st_flow -= 2*delta;

                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                        if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                          vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 2 ) {
                            /* Added (-) charge to =O and (+) charge to N => nDeltaCharge == 2 */
                            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                            if ( ret > 0 ) {
                                nNumRunBNS ++;
                                cur_success ++; /* 13 */
                            }
                        } else {
                            pe->flow += delta;
                            pv1->st_edge.flow += delta;
                            pv2->st_edge.flow += delta;
                            pBNS->tot_st_flow += 2*delta;
                        }
                        pe->forbidden        &= forbidden_edge_mask_inv; /* allow changes to O=X bond */
                        INCHI_HEAPCHK
                    }
                    RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                }
                CurrEdges.num_edges = 0; /* clear current edge list */
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
                    if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                        goto exit_function;  /* no fixed-H found */
                    }
                    if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                        goto exit_function;
                    }
                    if ( !pc2i->bHasDifference ) {
                        goto exit_function; /* nothing to do */
                    }
                    break;
                }/* else
                if ( bOnly_N_V ) {
                    bOnly_N_V = 0;
                }
                */
                break;
            }
        }
    }

    if ( (pc2i->nNumTgInChI <= 1 &&
        pc2i->nNumRemHInChI > pc2i->nNumRemHRevrs || pc2i->len_c2at) &&
         bHas_N_V( at2, pStruct->num_atoms) ) {
        /*-----------------------------------------------------------------*/
        /*                         |                         |             */
        /* case 14: restored:-N=AB=N=CD-XH original: (-)N-AB-N(+)=CD-XH    */
        /* FixH:              0    0   0/1              0            1     */
        /* MobH:              0    0   1/0              0            0     */
        /*                   non-taut  n/t             any  non     any    */
        /*                                                  taut           */
        /* X = O, S, Se, N; N = N                                          */
        /* t-group in original may have more (-) than in restored          */
        /*       same num_H and less (-) than in the restored structure    */
        /* atom N(V)/N(IV)(+) is not taut in both                          */
        /* The following transformation should be possible:                */
        /*        |                         |                              */
        /*   N=AB=N=CD-XH    ->     (-)N-AB-N-CD=XH(+)                     */
        /* This allows ADP to remove H(+) from -XH                         */
        /* As the result, the original structure has 0 or 1 t-group        */
        /* Solution: separate charges between -N(III)= and  N(V)           */
        /*-----------------------------------------------------------------*/
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;

        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        int num_N_V = 0, iat, i1, i2, i3, e1Flower, e1Plus, e2Plus, e2Minus, e3Plus;
        int max_success = pc2i->nNumRemHInChI - pc2i->nNumRemHRevrs;
        short iat_N_V_Array[MAX_DIFF_FIXH];
        EDGE_LIST iat_X_List, iat_N_III_List;
        AllocEdgeList( &iat_X_List, EDGE_LIST_CLEAR );
        AllocEdgeList( &iat_N_III_List, EDGE_LIST_CLEAR );
        cur_success = 0;
        ret = 0;
        for ( i = 0; i < pStruct->num_atoms; i ++ ) {
            iat = nCanon2AtnoRevrs[i];
            /* search for N(V), 3 bonds */
            if ( /* restored structure */
                 num_N_V < MAX_DIFF_FIXH &&
                 at2[iat].chem_bonds_valence == 5 && at2[iat].valence == 3 &&
                 !at2[iat].charge && !at2[iat].radical &&
                 pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                 !( at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint ) &&
                 !at2[iat].num_H &&
                 (e = pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pBNS->edge[e].flow /* no charge */ &&
                 NO_VERTEX != (j = GetChargeFlowerUpperEdge( pBNS, pVA, e )) && !pBNS->edge[j].forbidden &&
                 !pBNS->edge[j].flow /* neutral, valence=5 */ &&
                 /* orig. InChI */
                 !pStruct->endpoint[i] &&
                 !(nMobHInChI && nMobHInChI[i]) && !pStruct->fixed_H[i] ) {
                iat_N_V_Array[num_N_V ++] = iat;
            } else
            /* search for -N= */
            if ( /* restored structure */
                 at2[iat].chem_bonds_valence == 3 && at2[iat].valence == 2 &&
                 !at2[iat].charge && !at2[iat].radical &&
                 pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                 !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint ) &&
                 !at2[iat].num_H &&
                 (e = pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 !pBNS->edge[e].flow /* no charge */ &&
                 /* orig. InChI */
                 /*!pStruct->endpoint[i] &&*/
                 !(nMobHInChI && nMobHInChI[i]) && !pStruct->fixed_H[i] ) {

                if ( ret = AddToEdgeList( &iat_N_III_List, iat, 32 ) ) {
                    goto exit_case_14;
                }
            } else
            /* search for -OH -NH-, -NH2 */
            if ( /* restored structure */
                 at2[iat].chem_bonds_valence == at2[iat].valence &&
                 !at2[iat].charge && !at2[iat].radical &&
                 (pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 ||
                 pVA[iat].cNumValenceElectrons == 6 ) &&
                 at2[iat].num_H &&
                 (e = pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pBNS->edge[e].flow /* no charge */ &&
                 /* orig. InChI */
                 !(nMobHInChI && nMobHInChI[i]) && pStruct->fixed_H[i] ) {

                if ( ret = AddToEdgeList( &iat_X_List, iat, 32 ) ) {
                    goto exit_case_14;
                }
            }
        }
        if ( !max_success ) {
            max_success = inchi_min( num_N_V, iat_N_III_List.num_edges );
            max_success = inchi_min( max_success, iat_X_List.num_edges );
        }
        if ( num_N_V && iat_N_III_List.num_edges && iat_X_List.num_edges ) {
            for ( i1 = 0; i1 < num_N_V && cur_success < max_success; i1 ++ ) {
                int iat_N_V = iat_N_V_Array[i1];
                if ( NO_VERTEX == iat_N_V ||
                     0 >= (e1Plus = pVA[iat_N_V].nCPlusGroupEdge-1) ||
                     NO_VERTEX == (e1Flower = GetChargeFlowerUpperEdge( pBNS, pVA, e1Plus )) ||
                     1 != pBNS->edge[e1Plus].flow ||
                     0 != pBNS->edge[e1Flower].flow ) {
                    continue;
                }
                for ( i2 = iat_N_III_List.num_edges-1; 0 <= i2 && cur_success < max_success; i2 -- ) {
                    int iat_N_III = iat_N_III_List.pnEdges[i2];
                    if ( NO_VERTEX == iat_N_III ||
                         0 >= (e2Minus = pVA[iat_N_III].nCMinusGroupEdge-1) ||
                         0 >= (e2Plus  = pVA[iat_N_III].nCPlusGroupEdge-1) ||
                         0 != pBNS->edge[e2Minus].flow ||
                         1 != pBNS->edge[e2Plus].flow ) {
                        /* do not consider this atom anymore */
                        iat_N_III_List.pnEdges[i2] = NO_VERTEX;
                        continue;
                    }
                    for ( i3 = iat_X_List.num_edges-1; 0 <= i3 && cur_success < max_success; i3 -- ) {
                        int iat_X = iat_X_List.pnEdges[i3];
                        BNS_VERTEX *pv1n, *pv2n;
                        BNS_EDGE   *pe1n, *pe2n, *pe1Plus, *pe2Minus, *pe3Plus;
                        Vertex      v1n, v2n;
                        ret = 0;
                        if ( NO_VERTEX == iat_X ||
                             0 >= (e3Plus  = pVA[iat_X].nCPlusGroupEdge-1) ||
                             1 != pBNS->edge[e3Plus].flow ) {
                            /* do not consider this atom anymore */
                            iat_X_List.pnEdges[i3] = NO_VERTEX;
                            continue;
                        }
                        /* all is ready to check whether the following applies:
                           forbid changes of all charges and N,P,... flowers
                           allow to change edges: e2Minus, e3Plus
                           Increment flow in e1Flower
                           The result should be: increase in number of charges by 2
                        */
                        pe1Plus  = pBNS->edge + e1Plus;  /* N(V) positive charge edge */
                        pe2Minus = pBNS->edge + e2Minus; /* =N- negative charge edge */
                        pe3Plus  = pBNS->edge + e3Plus;  /* -XH positive charge edge */
                        pe       = pBNS->edge + e1Flower; /* N(V) flower edge */
                        pv1 = pBNS->vert + (v1 = pe->neighbor1);
                        pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);
                        for ( j = pv1->num_adj_edges-1; 0 <= j; j -- ) {
                            pe1n = pBNS->edge + pv1->iedge[j];
                            if ( pe1n->flow && !pe1n->forbidden && pe1n != pe1Plus ) {
                                pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                                break;
                            }
                        }
                        if ( j < 0 )
                            continue; /* not found -- should not happen */
                        for ( j = pv2->num_adj_edges-1; 0 <= j; j -- ) { /* was -2; changed 2006-2-28 12:35pm*/
                            pe2n = pBNS->edge + pv2->iedge[j];
                            if ( pe2n->flow && !pe2n->forbidden && pe2n != pe1Plus ) {
                                pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                                break;
                            }
                        }
                        if ( j < 0 )
                            continue; /* not found -- should not happen */
                        delta = 1;
                        pe->flow           += delta;
                        pe1n->flow         -= delta;
                        pe2n->flow         -= delta;
                        pv1n->st_edge.flow -= delta;
                        pv2n->st_edge.flow -= delta;
                        pBNS->tot_st_flow  -= 2*delta;
                        
                        SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                        SetForbiddenEdgeMask( pBNS, &OtherNFlowerEdges, forbidden_edge_mask );

                        /* allow two charges to change */
                        pe2Minus->forbidden &= forbidden_edge_mask_inv;
                        pe3Plus->forbidden  &= forbidden_edge_mask_inv;
                        /* test #1 */
                        ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                        INCHI_HEAPCHK
                        if ( ret < 0 ) {
                            goto exit_case_14;
                        } else
                        if ( ret == 1 && (vPathEnd == v1n && vPathStart == v2n ||
                                          vPathEnd == v2n && vPathStart == v1n) &&
                                          nDeltaCharge == 2 ) {
                            ; /* success */
                        } else {
                            ret = 0;
                        }
                        /* restore BNS */
                        pe2Minus->forbidden |= forbidden_edge_mask;
                        pe3Plus->forbidden  |= forbidden_edge_mask;
                        pe->flow           -= delta;
                        pe1n->flow         += delta;
                        pe2n->flow         += delta;
                        pv1n->st_edge.flow += delta;
                        pv2n->st_edge.flow += delta;
                        pBNS->tot_st_flow  += 2*delta;
                        if ( ret == 1 ) {
                            /* test #2: check if charge separation is possible */
                            pe->flow           += delta;
                            pe1n->flow         -= delta;
                            pe2n->flow         -= delta;
                            pv1n->st_edge.flow -= delta;
                            pv2n->st_edge.flow -= delta;
                            pBNS->tot_st_flow  -= 2*delta;

                            /* allow two charges (N(V) and N(III)) to change */
                            pe2Minus->forbidden &= forbidden_edge_mask_inv;
                            pe1Plus->forbidden  &= forbidden_edge_mask_inv;
                            /* test #2 */
                            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                              &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                            if ( ret == 1 && (vPathEnd == v1n && vPathStart == v2n ||
                                              vPathEnd == v2n && vPathStart == v1n) &&
                                              nDeltaCharge == 2 ) {
                                /* success; actually change charges */
                                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                                if ( ret > 0 ) {
                                    nNumRunBNS ++;
                                    cur_success ++; /* 14 */
                                }
                            }
                            if ( ret <= 0 ) {
                                /* failed: restore BNS flow */
                                pe->flow           -= delta;
                                pe1n->flow         += delta;
                                pe2n->flow         += delta;
                                pv1n->st_edge.flow += delta;
                                pv2n->st_edge.flow += delta;
                                pBNS->tot_st_flow  += 2*delta;
                            }
                            INCHI_HEAPCHK
                        }
                        RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                        RemoveForbiddenEdgeMask( pBNS, &OtherNFlowerEdges, forbidden_edge_mask  );
                        if ( ret > 0 ) {
                            /* do not repeat for the same atoms */
                            iat_N_V_Array[i1] = NO_VERTEX;
                            iat_N_III_List.pnEdges[i2] = NO_VERTEX;
                            iat_X_List.pnEdges[i3] = NO_VERTEX;
                        }
                        if ( ret < 0 ) {
                            goto exit_case_14;
                        }
                        if ( ret > 0 ) {
                            break;
                        }
                    } /* i3 cycle */
                    if ( ret > 0 ) {
                        break;
                    }
                } /* i2 cycle */
            }
        }
exit_case_14:
        AllocEdgeList( &iat_X_List, EDGE_LIST_FREE );
        AllocEdgeList( &iat_N_III_List, EDGE_LIST_FREE );
        RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        RemoveForbiddenEdgeMask( pBNS, &OtherNFlowerEdges, forbidden_edge_mask  );
        CurrEdges.num_edges = 0; /* clear current edge list */
        if ( ret < 0 ) {
            goto exit_function;
        }
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }


    if ( pc2i->nNumTgMRevrs > pc2i->nNumTgMInChI ||
         pc2i->nNumRemHRevrs < pc2i->nNumRemHInChI ||
         pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI ||
         pc2i->nNumTgInChI <= 1 && pc2i->nNumTgRevrs > pc2i->nNumTgInChI ) {
        /*--------------------------------------------------------------*/
        /* case 15: restored: -(+)O=AB-N<  orig: -O-AB=N(+)<            */
        /* (a) restored t-groups have more (-) than in original InChI   */
        /* (b) Mobile-H    charge: restored > original InChI *and*      */
        /*              removed H: restored < original InChI            */
        /* (c) restored t-groups have less endpnoits than in orig InChI */
        /* O = O, S, Se, Te; N = N                                      */
        /* Solution: move (+) from -O(+)= to -N<                        */
        /*--------------------------------------------------------------*/
        int num_SB_Neutr = 0, num_DB_Charged = 0, iat;
        short iat_SB_Neutr[MAX_DIFF_FIXH], iat_DB_Charged[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        cur_success = 0;
        /* search for -O(+)= */
        /* search for -N< */
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* -O(+)= in restored atom: charge=1, has no H, a double bond */
                 num_DB_Charged < MAX_DIFF_FIXH &&
                 at2[iat].charge == 1 && !at2[iat].num_H &&
                 at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 (pVA[iat].cNumValenceElectrons == 6 ) &&
                 /* in orig.InChI: an endpoint, has fixed-H */
                 /*pStruct->endpoint[i] &&*/
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has (+) edge */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && 0 == pBNS->edge[e].forbidden ) {

                iat_DB_Charged[num_DB_Charged ++] = iat;
                /*
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
                */
            } else
            if ( /* -N< in restored atom: charge=0, has no H, has no double bond, N only */
                 num_SB_Neutr < MAX_DIFF_FIXH &&
                 at2[iat].charge == 0 && !at2[iat].num_H &&
                 at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 (pVA[iat].cNumValenceElectrons == 5 &&
                  pVA[iat].cPeriodicRowNumber == 1 ) &&
                 /* in orig.InChI: an endpoint, has fixed-H */
                 /*pStruct->endpoint[i] &&*/
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has (+) edge */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && 0 == pBNS->edge[e].forbidden ) {

                iat_SB_Neutr[num_SB_Neutr ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_SB_Neutr, num_DB_Charged ) ) {
            /* detected; attempt to fix */
            BNS_VERTEX *pv1n, *pv2n;
            BNS_EDGE   *pe1n, *pe2n;
            Vertex      v1n, v2n;
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_DB_Charged && cur_success < num_try; i ++ ) {
                iat = iat_DB_Charged[i];
                pe   = pBNS->edge + pVA[iat].nCPlusGroupEdge-1;
                if ( pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                for ( j = pv1->num_adj_edges-1; 0 <= j; j -- ) {
                    pe1n = pBNS->edge + pv1->iedge[j];
                    if ( pe1n->flow && !pe1n->forbidden ) {
                        pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                        break;
                    }
                }
                if ( j < 0 )
                    continue; /* not found */
                
                for ( j = pv2->num_adj_edges-1; 0 <= j; j -- ) { /* was -2; changed 2006-2-28 12:35pm*/
                    pe2n = pBNS->edge + pv2->iedge[j];
                    if ( pe2n->flow && !pe2n->forbidden ) {
                        pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                        break;
                    }
                }
                if ( j < 0 )
                    continue; /* not found */

                pe->flow   += delta;
                pe1n->flow -= delta;
                pe2n->flow -= delta;
                pv1n->st_edge.flow -= delta;
                pv2n->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1n && vPathStart == v2n ||
                                  vPathEnd == v2n && vPathStart == v1n) &&
                                  (nDeltaCharge == 0 || nDeltaCharge == 1) ) {
                    /* Moved charge from O(+) to -N< => nDeltaCharge == 1 or 0 if pe2n = -N< charge edge */
                    /* Flow change on pe (+)charge edge (atom NH2) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 15 */
                    }
                } else {
                    pe->flow   -= delta;
                    pe1n->flow += delta;
                    pe2n->flow += delta;
                    pv1n->st_edge.flow += delta;
                    pv2n->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->nNumTgDiffMinus ) {
        /*----------------------------------------------------------------*/
        /* case 16: restored: O=X-NH(-)      orig.:  O(-)-X=NH            */
        /*            t-group: (H,-)                  (2H)                */
        /* O(-) = S, Se, Te; N = N;                                       */
        /* Solution: move (-) from O(-) to -NH(-)                         */
        /*----------------------------------------------------------------*/
        int num_SB_N_Minus = 0, num_DB_O_Neutr = 0, iat, itg;
        short iat_SB_N_Minus[MAX_DIFF_FIXH], iat_DB_O_Neutr[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        cur_success = 0;
        for ( itg = 0; itg < pStruct->ti.num_t_groups && itg < pStruct->One_ti.num_t_groups; itg ++ ) {
            if ( pStruct->ti.t_group[itg].nNumEndpoints != pStruct->One_ti.t_group[itg].nNumEndpoints ||
                 pStruct->ti.t_group[itg].num[1]  >= pStruct->One_ti.t_group[itg].num[1] ) {
                continue;
            }
            CurrEdges.num_edges = num_SB_N_Minus = num_DB_O_Neutr = 0;
            cur_success = 0;
            for ( j = 0, k = pStruct->One_ti.t_group[itg].nFirstEndpointAtNoPos;
                    j < pStruct->One_ti.t_group[itg].nNumEndpoints; j ++ ) {
                i = pStruct->One_ti.nEndpointAtomNumber[k+j]; /* canonical number in restored struct. */
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom: charge=0, has no H, has double bond, O, S, Se, Te */
                     num_DB_O_Neutr < MAX_DIFF_FIXH &&
                     at2[iat].charge == 0 && !at2[iat].num_H &&
                     at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                     pVA[iat].cNumValenceElectrons == 6 &&
                     /* in orig.InChI: an endpoint, may have fixed-H */
                     pStruct->endpoint[i] && 
                     /*!(pStruct->fixed_H && pStruct->fixed_H[i]) &&*/
                     !(nMobHInChI && nMobHInChI[i] ) &&
                     /* has (-) edge */
                     (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {

                    iat_DB_O_Neutr[num_DB_O_Neutr ++] = iat;
                    
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                    
                } else
                if ( /* in restored atom: charge=-1, has H, has double bond, N */
                     num_SB_N_Minus < MAX_DIFF_FIXH &&
                     at2[iat].charge == -1 && at2[iat].num_H &&
                     at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                     pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1 &&
                     /* in orig.InChI: an endpoint, has no fixed-H */
                     pStruct->endpoint[i] &&
                     (pStruct->fixed_H && pStruct->fixed_H[i]) &&
                     !(nMobHInChI && nMobHInChI[i] ) &&
                     /* has (-) edge */
                     (e=pVA[iat].nCMinusGroupEdge-1) >= 0 &&
                     0 == pBNS->edge[e].forbidden ) {

                    iat_SB_N_Minus[num_SB_N_Minus ++] = iat;
                    /*
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                    */
                }
            }
            if ( num_try = inchi_min( num_SB_N_Minus, num_DB_O_Neutr ) ) {
                /* detected; attempt to fix */
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
                /* allow stereobonds in rings change */
                /*
                if ( forbidden_stereo_edge_mask )
                    RemoveForbiddenEdgeMask( pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask );
                */
                delta = 1;
                for ( i = 0; i < num_SB_N_Minus && cur_success < num_try; i ++ ) {
                    iat = iat_SB_N_Minus[i];
                    pe   = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                    if ( !pe->flow )
                        continue;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    /*pe->forbidden |= forbidden_edge_mask;*/
                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2*delta;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                        /* Moved (-) charge to =O => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (atom -NH(-)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if ( ret > 0 ) {
                            nNumRunBNS ++;
                            cur_success ++; /* 16 */
                        }
                    } else {
                        pe->forbidden &= forbidden_edge_mask_inv;
                        pe->flow += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2*delta;
                    }
                    INCHI_HEAPCHK
                }
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                /*
                if ( forbidden_stereo_edge_mask )
                    SetForbiddenEdgeMask( pBNS, &FixedLargeRingStereoEdges, forbidden_stereo_edge_mask );
                */
            }
            CurrEdges.num_edges = 0; /* clear current edge list */
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
                if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                    goto exit_function;  /* no fixed-H found */
                }
                if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                    goto exit_function;
                }
                if ( !pc2i->bHasDifference ) {
                    goto exit_function; /* nothing to do */
                }
            }
        }
    }

    if ( pc2i->nNumRemHInChI < pc2i->nNumRemHRevrs ) {
        /*--------------------------------------------------------------*/
        /* case 17: restored: OH(+)=AB-O-     orig.  HO-AB=O(+)-        */
        /* number of removed H:  n+m                     n              */
        /* OH(+) = N, O, S, Se; -O- = P,As,O,S,Se,Te,F,Cl,Br,I          */
        /* Solution: move (+) from OH(+) to -O-                         */
        /*--------------------------------------------------------------*/
        int num_SB_Neutr = 0, num_DB_Charged = 0, iat;
        short iat_SB_Neutr[MAX_DIFF_FIXH], iat_DB_Charged[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        cur_success = 0;
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom: charge=+1, has H, has double bond, N, O, S, Se, Te */
                 num_DB_Charged < MAX_DIFF_FIXH &&
                 at2[iat].charge == 1 && at2[iat].num_H &&
                 at2[iat].valence < at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 (pVA[iat].cNumValenceElectrons == 6 ||
                  pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1) &&
                 /* has (+) edge */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden ) {

                iat_DB_Charged[num_DB_Charged ++] = iat;
                /*
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
                */
            } else
            if ( /* in restored atom: charge=0, has no H, has no double bond, N, P, O, S, Se, Te */
                 num_SB_Neutr < MAX_DIFF_FIXH &&
                 at2[iat].charge == 0 && !at2[iat].num_H &&
                 at2[iat].valence == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 (pVA[iat].cNumValenceElectrons == 6 || pVA[iat].cNumValenceElectrons == 7 ||
                  pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber > 1 ) &&
                 /* in orig.InChI: not an endpoint */
                 !pStruct->endpoint[i] &&
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 !(nMobHInChI && nMobHInChI[i] ) &&
                 /* has (+) edge */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 &&
                 0 == pBNS->edge[e].forbidden ) {

                iat_SB_Neutr[num_SB_Neutr ++] = iat;
                
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( num_SB_Neutr, num_DB_Charged ) ) {
            BNS_VERTEX *pv1n, *pv2n;
            BNS_EDGE   *pe1n, *pe2n;
            Vertex      v1n, v2n;

            num_try = inchi_min( num_try, pc2i->nNumRemHRevrs-pc2i->nNumRemHInChI);
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_DB_Charged && cur_success < num_try; i ++ ) {
                iat = iat_DB_Charged[i];
                pe   = pBNS->edge + pVA[iat].nCPlusGroupEdge-1;
                if ( pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                for ( j = pv1->num_adj_edges-1; 0 <= j; j -- ) {
                    pe1n = pBNS->edge + pv1->iedge[j];
                    if ( pe1n->flow && !pe1n->forbidden ) {
                        pv1n = pBNS->vert + (v1n = pe1n->neighbor12 ^ v1);
                        break;
                    }
                }
                if ( j < 0 )
                    continue; /* not found */
                
                for ( j = pv2->num_adj_edges-1; 0 <= j; j -- ) { /* was -2; changed 2006-2-28 12:35pm*/
                    pe2n = pBNS->edge + pv2->iedge[j];
                    if ( pe2n->flow && !pe2n->forbidden ) {
                        pv2n = pBNS->vert + (v2n = pe2n->neighbor12 ^ v2);
                        break;
                    }
                }
                if ( j < 0 )
                    continue; /* not found */

                pe->flow   += delta;
                pe1n->flow -= delta;
                pe2n->flow -= delta;
                pv1n->st_edge.flow -= delta;
                pv2n->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1n && vPathStart == v2n ||
                                  vPathEnd == v2n && vPathStart == v1n) &&
                                  (nDeltaCharge == 0 || nDeltaCharge == 1) ) {
                    /* Moved charge from OH(+) to -O- => nDeltaCharge == 1 or 0 if pe2n = -O- charge edge */
                    /* Flow change on pe (+)charge edge (atom OH(+)) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 17 */
                    }
                } else {
                    pe->flow   -= delta;
                    pe1n->flow += delta;
                    pe2n->flow += delta;
                    pv1n->st_edge.flow += delta;
                    pv2n->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( (pc2i->nNumTgInChI && pStruct->endpoint &&
         pc2i->nNumTgMInChI > pc2i->nNumTgMRevrs && pc2i->nNumEndpInChI > pc2i->nNumEndpRevrs ) ) {
        /*-----------------------------------------------------------------*/
        /*                                                                 */
        /* case 18: restored:-N=AB-X                -(-)N-AB-X(+)          */
        /* FixH:              0    0                    0    0             */
        /* MobH:              0    0                    0    0             */
        /*                   non  non                 taut  non            */
        /*                  taut  taut                      taut           */
        /* X = any heteroatom   N=N                                        */
        /* t-group in original has (Hn,-m) in the restored: (Hn,-m+1)      */
        /*       same num_H and more (-) than in the restored structure    */
        /* atom X is not taut in both                                      */
        /* Solution: separate charges between -N(III)= and  X              */
        /*-----------------------------------------------------------------*/
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        int iat, e1, itg, max_success;
        CurrEdges.num_edges = 0;
        cur_success = 0;
        ret = 0;
        /* search for -N= */
        for ( itg = 0; itg < pStruct->ti.num_t_groups && itg < pStruct->One_ti.num_t_groups; itg ++ ) {
            if ( pStruct->ti.t_group[itg].nNumEndpoints <= pStruct->One_ti.t_group[itg].nNumEndpoints ||
                 pStruct->ti.t_group[itg].num[1]  <= pStruct->One_ti.t_group[itg].num[1] ) {
                     continue;
            }
            CurrEdges.num_edges = 0;
            cur_success = 0;
            for ( j = 0, k = pStruct->ti.t_group[itg].nFirstEndpointAtNoPos;
                    j < pStruct->ti.t_group[itg].nNumEndpoints; j ++ ) {
                i = pStruct->ti.nEndpointAtomNumber[k+j]; /* canonical number in restored struct. */
                iat = nCanon2AtnoRevrs[i];
                if ( !pStruct->endpoint[i] || !at_Mobile_H_Revrs || at_Mobile_H_Revrs[iat].endpoint ||
                     pVA[i].cNumValenceElectrons != 5 || pVA[i].cPeriodicRowNumber != 1 ||
                     2 != at2[iat].valence || at2[iat].num_H || at2[iat].radical ||
                     0 <= (e1=pVA[iat].nCPlusGroupEdge-1) && !pBNS->edge[e1].flow ||
                     0 > (e=pVA[iat].nCMinusGroupEdge-1) || pBNS->edge[e].forbidden || pBNS->edge[e].flow ) {
                    continue;
                }
                /* found -N= */
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( !(max_success = CurrEdges.num_edges) ) {
            goto exit_case_18;
        }
        /* search for X */
        for ( i = 0; i < pStruct->num_atoms && cur_success < max_success; i ++ ) {
            iat =  nCanon2AtnoRevrs[i];
            if ( pStruct->endpoint[i] || !pVA[i].cNumValenceElectrons || pVA[i].cNumValenceElectrons == 4 ||
                 at2[iat].num_H || at2[iat].radical ||
                 0 <= (e1=pVA[iat].nCMinusGroupEdge-1) && !pBNS->edge[e1].flow ||
                 0 > (e=pVA[iat].nCPlusGroupEdge-1) || pBNS->edge[e].forbidden || pBNS->edge[e].flow != 1 ) {
                continue;
            }
            /* try to move the charge */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            SetForbiddenEdgeMask( pBNS, &OtherNFlowerEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask );

            pe   = pBNS->edge + e;
            if ( !pe->flow )
                continue;
            pv1 = pBNS->vert + (v1 = pe->neighbor1);
            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

            delta = 1;
            pe->flow -= delta;
            pv1->st_edge.flow -= delta;
            pv2->st_edge.flow -= delta;
            pBNS->tot_st_flow -= 2*delta;

            ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                  &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

            if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                              vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                /* Created (-) charge on -N= => nDeltaCharge == 1 */
                /* Flow change on pe (+)charge edge (atom X) is not known to RunBnsTestOnce()) */
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                if ( ret > 0 ) {
                    nNumRunBNS ++;
                    cur_success ++; /* 18 */
                }
            } else {
                pe->flow += delta;
                pv1->st_edge.flow += delta;
                pv2->st_edge.flow += delta;
                pBNS->tot_st_flow += 2*delta;
            }
            INCHI_HEAPCHK
        }
exit_case_18:
        RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        RemoveForbiddenEdgeMask( pBNS, &OtherNFlowerEdges, forbidden_edge_mask  );
        CurrEdges.num_edges = 0; /* clear current edge list */
        if ( ret < 0 ) {
            goto exit_function;
        }
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }
    if ( pc2i->len_c2at >= 1 ) {
        /*--------------------------------------------------------------*/
        /* case 19 restored:       M--OH   original:  M(-)==OH(+)       */
        /* FixH:               metal  0                      1          */
        /* MobH:                      1                      0          */
        /* O =  O, S, Se, Te; not taut. in InChI                        */
        /* In restored structure has H; tautomeric or not tautomeric    */
        /* Solution: move (+) from -OH to M; charhe on M may vary       */
        /*--------------------------------------------------------------*/
        int iat;
        EdgeIndex eOHPlus, eMPlus, eMMinus, eOMBond;
        BNS_EDGE  *peOHPlus, *peMPlus, *peMMinus, *peOMBond;
        int       iatMetal, ChargeOnMetal, DeltaChargeExpected;
        cur_success = 0;
        num_zero_ret = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: =NH2(+), =OH(+) */
                 (pc2i->c2at[i].nValElectr == 6 ) /* N, O, S, Se, Te */ &&
                 /*!pc2i->c2at[i].endptInChI &&*/ /* <=== relaxation */
                 (e=pVA[iat].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                 pc2i->c2at[i].nFixHInChI == 1 && pc2i->c2at[i].nMobHInChI == 0 &&
                 /* reversed structure info: */
                 pc2i->c2at[i].nFixHRevrs ==  0 && pc2i->c2at[i].nMobHRevrs == 1 &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && at2[iat].num_H &&
                 at2[iat].valence == 1 &&
                 at2[iat].valence == at2[iat].chem_bonds_valence &&
                 /* metal atom */
                 pVA[iatMetal=at2[iat].neighbor[0]].cMetal && 
                 (eMPlus=pVA[iatMetal].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[eMPlus].forbidden &&
                 (eMMinus=pVA[iatMetal].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[eMMinus].forbidden &&
                 !pBNS->edge[eOMBond=pBNS->vert[iat].iedge[0]].forbidden
                 ) {

                /* -OH charge edges */
                if ( ret = AddToEdgeList( &CurrEdges, iat, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( CurrEdges.num_edges ) {
            /* detected; fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            SetForbiddenEdgeMask( pBNS, &NFlowerEdges, forbidden_edge_mask );
            SetForbiddenEdgeMask( pBNS, &AllBondEdges, forbidden_edge_mask );
            for ( i = 0; i < CurrEdges.num_edges; i ++ ) {
                /* v1 is -OH, v2 is adjacent to it Metal */
                iat      = CurrEdges.pnEdges[i];
                iatMetal = at2[iat].neighbor[0];
                peOHPlus = pBNS->edge + (eOHPlus = pVA[iat].nCPlusGroupEdge-1);
                peMPlus  = pBNS->edge + (eMPlus  = pVA[iatMetal].nCPlusGroupEdge-1);
                peMMinus = pBNS->edge + (eMMinus = pVA[iatMetal].nCMinusGroupEdge-1);
                peOMBond = pBNS->edge + (eOMBond =pBNS->vert[iat].iedge[0]);
                /* remove forbidden edge masks */
                peMPlus->forbidden  &= forbidden_edge_mask_inv;
                peMMinus->forbidden &= forbidden_edge_mask_inv;
                peOMBond->forbidden &= forbidden_edge_mask_inv;

                ChargeOnMetal = (peMPlus->cap - peMPlus->flow) - peMMinus->flow;
                if ( 1 == ChargeOnMetal ) {
                    /* We are going to subtract 1 from the charge on Metal */
                    /* Added (+)charge to -OH is not known to RunBnsTestOnce() */
                    DeltaChargeExpected = -1; /* charge will become = 0 */
                } else
                if ( 0 == ChargeOnMetal ) {
                    DeltaChargeExpected = 1; /* charge on Metal will be created */
                } else {
                    DeltaChargeExpected = 0;
                }

                delta = 1;
                pe   = peOHPlus;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->flow -= delta; /* remove (-) from AB-O(-) */
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == DeltaChargeExpected ) {
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 19 */
                    }
                } else {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
                /* set forbidden edge masks back */
                peMPlus->forbidden  |= forbidden_edge_mask;
                peMMinus->forbidden |= forbidden_edge_mask;
                peOMBond->forbidden |= forbidden_edge_mask;
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &NFlowerEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &AllBondEdges, forbidden_edge_mask );

            CurrEdges.num_edges = 0; /* clear current edge list */
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
                if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                    goto exit_function;  /* no fixed-H found */
                }
                if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                    goto exit_function;
                }
                if ( !pc2i->bHasDifference ) {
                    goto exit_function; /* nothing to do */
                }
            }
        }
    }
    if ( pc2i->len_c2at > 1 && pc2i->nNumTgRevrs && pc2i->nNumTgInChI) {
        /*--------------------------------------------------------------*/
        /* case 20: restored:  O(-)-AB=N-   original:   O=AB-N(-)-      */
        /* FixH:               0       0                0     -1        */
        /* MobH:               0       0                0      1        */
        /*                   taut    non-taut       non-taut taut       */
        /*                           or taut                  no H      */
        /*                           no H                               */
        /* O = O, S, Se; N = N, O, S, Se, Te;                           */
        /* restored atoms are taut/non-taut; original are opposite.     */
        /* Solution: move (-) from O(-) to =N-                          */
        /*--------------------------------------------------------------*/
        int num_SB_O_Minus = 0, num_DB_N = 0, iat;
        short iat_SB_O_Minus[MAX_DIFF_FIXH], iat_DB_N[MAX_DIFF_FIXH];
        
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        /*
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        cur_success = 0;
        CurrEdges.num_edges = 0; /* clear current edge list */
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: =O or -N= */
                 num_DB_N < MAX_DIFF_FIXH &&
                 pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pBNS->edge[e].flow == 0 &&
                 pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                 /* if  more than 1 t-group are in orig. InChI then do not move (-) to N */
                 (pc2i->nNumTgInChI == 1 || pc2i->c2at[i].nValElectr == 6) &&
                 /* reversed structure info: */
                 !pc2i->c2at[i].endptRevrs &&
                 pc2i->c2at[i].nFixHRevrs ==  0 && /*pc2i->c2at[i].nMobHRevrs == 0 &&*/
                 pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                 at2[iat].valence + 1 == at2[iat].chem_bonds_valence ) {
                iat_DB_N[num_DB_N ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            } else
            if ( /* orig. InChI info: -O(-) */
                 num_SB_O_Minus < MAX_DIFF_FIXH &&
                 !pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pBNS->edge[e].flow == 1 &&
                 pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                 pc2i->c2at[i].nValElectr == 6 &&
                 /* reversed structure info: */
                 pc2i->c2at[i].endptRevrs &&
                 pc2i->c2at[i].nFixHRevrs ==  0 && pc2i->c2at[i].nMobHRevrs == 0 &&
                 pc2i->c2at[i].nAtChargeRevrs == -1 && !at2[iat].num_H &&
                 at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 1 ) {
                iat_SB_O_Minus[num_SB_O_Minus ++] = iat;
            }
        }
        if ( !num_DB_N ) {
            /* search among N that are tautomeric in both cases */
            for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
                if ( !pStruct->endpoint[i] ) {
                    continue;
                }
                iat = nCanon2AtnoRevrs[i];
                if ( /* in restored atom O: charge=-1, no H, has no double bond, endpoint */
                     num_DB_N < MAX_DIFF_FIXH &&
                     at2[iat].charge == 0 && !at2[iat].num_H &&
                     at2[iat].valence + 1 == at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                     /* in orig.InChI: an endpoint, has no H */
                     !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                     /*!(nMobHInChI && nMobHInChI[i] ) &&*/
                     /* has (-) edge */
                     (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                     !pBNS->edge[e].flow ) {

                    iat_DB_N[num_DB_N ++] = iat;
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                }
            }
        }
        if ( num_try = inchi_min( num_SB_O_Minus, num_DB_N ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_O_Minus && cur_success < num_try; i ++ ) {
                iat = iat_SB_O_Minus[i];
                pe   = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Added (-) charge to =N- => nDeltaCharge == 1 */
                    /* Flow change on pe (-)charge edge (atom -O(-)) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 20 */
                    }
                } else {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }
    if ( pc2i->len_c2at && pc2i->nNumTgRevrs && pc2i->nNumTgHInChI && pStruct->endpoint ) {
        /*--------------------------------------------------------------*/
        /*                      O(-)                      O             */
        /*                      |                         ||            */
        /* case 21: restored: R=S=O         original:   R-S=O           */
        /*                      |                         |             */
        /*                      O(-)                      O(-)          */
        /*                           All O are taut     R is not taut   */
        /*                                                              */
        /* In addition, another atom O that should have been tautomeric */
        /* or has H(+) added in Mobile-H layer is not like that         */
        /* O = O, S, Se;  S=S, Se, Te                                  */
        /* Solution: move (-) from O(-) to =O                           */
        /*           these atoms are tautomeric in restored structure   */
        /*--------------------------------------------------------------*/
        int num_SB_O_Minus = 0, num_DB_O = 0, iat, iS;
        short iat_SB_O_Minus[MAX_DIFF_FIXH], iat_Central[MAX_DIFF_FIXH], iat_DB_O[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        /*
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        CurrEdges.num_edges = 0; /* clear current edge list */
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: =O    */
                 num_DB_O < MAX_DIFF_FIXH &&
                 pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                 (pc2i->c2at[i].endptInChI || pc2i->c2at[i].nMobHInChI) &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                 /* reversed structure info: */
                 !(pc2i->c2at[i].endptRevrs || pc2i->c2at[i].nMobHRevrs) &&
                 pc2i->c2at[i].nFixHRevrs ==  0 &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                 at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 ) {
                iat_DB_O[num_DB_O ++] = iat;
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        for ( i = 0; num_DB_O && i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            if ( !pStruct->endpoint[i] ) {
                continue;
            }
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom O: charge=-1, no H, has no double bond, endpoint */
                 num_SB_O_Minus < MAX_DIFF_FIXH &&
                 at2[iat].charge == -1 && !at2[iat].num_H &&
                 at2[iat].valence == 1 && at2[iat].chem_bonds_valence && !pVA[iat].cMetal &&
                 pVA[iat].cNumValenceElectrons == 6 &&
                 (at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint) &&
                 /* in orig.InChI: an endpoint, has no H */
                 !(pStruct->fixed_H && pStruct->fixed_H[i]) &&
                 /*!(nMobHInChI && nMobHInChI[i] ) &&*/
                 /* has (-) edge */
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pBNS->edge[e].flow ) {
                int nNumTautSB = 0, nNumTautDB = 0, nNumOtherDB = 0, nNumOtherSB = 0, nNumOthers = 0, nNumNegEndp = 0;
                /* traverse neighbors of the centerpoint iS */
                iS = at2[i].neighbor[0];
                for ( j = 0; j < num_SB_O_Minus; j ++ ) {
                    if ( iat_Central[j] == iS )
                        break;
                }
                if ( j < num_SB_O_Minus ) {
                    continue;  /* have already been there */
                }
                for ( j = 0; j < at[iS].valence; j ++ ) {
                    int bond_type = at2[iS].bond_type[j];
                    k = at2[iS].neighbor[j];
                    if ( k == i ) {
                        continue;
                    }
                    if ( pStruct->endpoint[k] == pStruct->endpoint[i] ) {
                        nNumTautSB += ( bond_type == BOND_TYPE_SINGLE );
                        nNumTautDB += ( bond_type == BOND_TYPE_DOUBLE );
                    } else
                    if ( bond_type == BOND_TYPE_DOUBLE ) {
                        nNumOtherDB ++;
                    } else
                    if ( bond_type == BOND_TYPE_SINGLE ) {
                        nNumOtherSB ++;
                    } else {
                        nNumOthers ++;
                    }
                    if ( at2[k].endpoint == at2[i].endpoint && at2[k].valence == 1 &&
                         at2[k].charge    == -1 && pVA[k].cNumValenceElectrons == 6 ) {
                        nNumNegEndp ++;
                    }
                }
                if ( !nNumTautSB ) {
                    continue;
                }
                if ( !( nNumOtherDB && nNumTautDB ) ) {
                    continue; /* ignore */
                }
                
                iat_SB_O_Minus[num_SB_O_Minus] = iat;
                iat_Central[num_SB_O_Minus ++] = iS;
            }
        }
        if ( num_try = inchi_min( num_SB_O_Minus, num_DB_O ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_SB_O_Minus && cur_success < num_try; i ++ ) {
                iat = iat_SB_O_Minus[i];
                pe   = pBNS->edge + pVA[iat].nCMinusGroupEdge-1;
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden |= forbidden_edge_mask;
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Added (-) charge to =O => nDeltaCharge == 1 */
                    /* Flow change on pe (-)charge edge (atom -N(-)-) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 21 */
                    }
                } else {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->len_c2at && pc2i->nNumTgRevrs && pc2i->nNumEndpInChI < pc2i->nNumEndpRevrs ) {
        /*--------------------------------------------------------------*/
        /*                      O                         O             */
        /*                      ||                        ||            */
        /* case 21a:restored: R=S-R' =X     original:   R-S-R' -X(-)    */
        /*                      |                         ||            */
        /*                      O(-)                      O(-)          */
        /*             All O and X are taut      O and X are not taut   */
        /*             it is possible that X is R                       */
        /*                                                              */
        /* O = O, S, Se;  S=S, Se, Te; X = N, O, S, Se, Te              */
        /* Solution: move (-) from O(-) to =X                           */
        /*           these atoms are tautomeric in restored structure   */
        /*--------------------------------------------------------------*/
        int iat, iS;
        /*
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        */
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                    pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        /*
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        EDGE_LIST  OtherSO, CentralS, SOMinus, MinusAcceptord;
        CurrEdges.num_edges = 0; /* clear current edge list */
        AllocEdgeList( &OtherSO, EDGE_LIST_CLEAR );
        AllocEdgeList( &CentralS, EDGE_LIST_CLEAR );
        AllocEdgeList( &SOMinus, EDGE_LIST_CLEAR );
        AllocEdgeList( &MinusAcceptord, EDGE_LIST_CLEAR );
        cur_success = 0;
        if ( !at_Mobile_H_Revrs ) {
            goto exit_case_21a;
        }
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: -X(-)    */
                 /*num_DB_O < MAX_DIFF_FIXH &&*/
                 /*pc2i->c2at[i].nValElectr == 6 */ /* O, S, Se, Te */
                 !pc2i->c2at[i].endptInChI &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                 /* reversed structure info: */
                 (pc2i->c2at[i].endptRevrs || pc2i->c2at[i].nMobHRevrs) &&
                 pc2i->c2at[i].nFixHRevrs ==  0 &&
                 /*pc2i->c2at[i].nAtChargeRevrs == 0 &&*/ !at2[iat].num_H ) {
                if ( pVA[iat].cNumValenceElectrons == 6 && at2[iat].charge == -1 &&
                     pBNS->edge[e].flow &&
                     at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 1 &&
                     pVA[iS=(int)at2[iat].neighbor[0]].cNumValenceElectrons == 6 && pVA[iS].cPeriodicRowNumber > 1 &&
                     at2[iS].valence >= 4 ) {
                    /* a candidate for S in -SO2- */
                    int nNumTautSB = 0, nNumTautDB = 0, nNumOtherDB = 0, nNumOtherSB = 0;
                    int nNumOthers = 0, nNumNegEndp = 0, nNumEndpO = 0;
                    /* check whether we have already found it */
                    if ( 0 <= FindInEdgeList( &CentralS, iS ) ) {
                        continue;
                    }
                    for ( j = 0; j < at[iS].valence; j ++ ) {
                        int bond_type = at2[iS].bond_type[j];
                        k = at2[iS].neighbor[j];
                        if ( k == iat ) {
                            continue;
                        }
                        if ( pc2i->c2at[i].endptRevrs == at_Mobile_H_Revrs[k].endpoint && !at2[k].endpoint ) {
                            nNumTautSB += ( bond_type == BOND_TYPE_SINGLE );
                            nNumTautDB += ( bond_type == BOND_TYPE_DOUBLE );
                            nNumEndpO  += (pVA[k].cNumValenceElectrons == 6 && at2[k].valence == 1);
                        } else
                        if ( bond_type == BOND_TYPE_DOUBLE ) {
                            nNumOtherDB ++;
                        } else
                        if ( bond_type == BOND_TYPE_SINGLE ) {
                            nNumOtherSB ++;
                        } else {
                            nNumOthers ++;
                        }
                        if ( at2[k].endpoint == at2[i].endpoint && at2[k].valence == 1 &&
                             at2[k].charge    == -1 && pVA[k].cNumValenceElectrons == 6 ) {
                            nNumNegEndp ++;
                        }
                    }
                    if ( !nNumEndpO ) {
                        continue;
                    }
                    if ( nNumTautSB + nNumTautDB + nNumOtherDB <= nNumEndpO  ) {
                        continue; /* ignore */
                    }
                    /* collect double bond taut =O */
                    for ( j = 0; j < at[iS].valence; j ++ ) {
                        int bond_type = at2[iS].bond_type[j];
                        k = at2[iS].neighbor[j];
                        if ( pc2i->c2at[i].endptRevrs == at_Mobile_H_Revrs[k].endpoint &&
                             !at2[k].endpoint && pVA[k].cNumValenceElectrons == 6 && at2[k].valence == 1 &&
                             0 <= (e=pVA[k].nCMinusGroupEdge-1) && !pBNS->edge[e].forbidden ) {
                            if ( bond_type == BOND_TYPE_DOUBLE && !at2[k].charge && !pBNS->edge[e].flow) {
                                /* charges to be unchanged */
                                if ( ret = AddToEdgeList( &OtherSO, e, INC_ADD_EDGE ) ) {
                                    goto exit_case_21a;
                                }
                            } else
                            if ( bond_type == BOND_TYPE_SINGLE && at2[k].charge == -1 && pBNS->edge[e].flow ) {
                                /* charges to be removed */
                                if ( ret = AddToEdgeList( &SOMinus, e, INC_ADD_EDGE ) ) {
                                    goto exit_case_21a;
                                }
                            }
                        }
                    }
                    if ( ret = AddToEdgeList( &CentralS, iS, INC_ADD_EDGE ) ) {
                        goto exit_case_21a;
                    }
                } else
                if ( at2[iat].charge == 0 && !pBNS->edge[e].flow &&
                     at2[iat].valence + 1 == at2[iat].chem_bonds_valence ) {
                    /* changeable charges */
                    if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                        goto exit_function;
                    }
                }
            }
        }
        /* remove unchangeable from changeable */
        for ( i = 0; i < OtherSO.num_edges; i ++ ) {
            RemoveFromEdgeListByValue( &CurrEdges, OtherSO.pnEdges[i] );
        }

        if ( num_try = inchi_min( SOMinus.num_edges, CurrEdges.num_edges ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < SOMinus.num_edges && cur_success < num_try; i ++ ) {
                pe   = pBNS->edge + SOMinus.pnEdges[i];
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                /*pe->forbidden |= forbidden_edge_mask;*/
                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 1 ) {
                    /* Added (-) charge to =O => nDeltaCharge == 1 */
                    /* Flow change on pe (-)charge edge (atom -N(-)-) is not known to RunBnsTestOnce()) */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 21a */
                    }
                } else {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
exit_case_21a:
        CurrEdges.num_edges = 0; /* clear current edge list */
        AllocEdgeList( &OtherSO, EDGE_LIST_FREE );
        AllocEdgeList( &CentralS, EDGE_LIST_FREE );
        AllocEdgeList( &SOMinus, EDGE_LIST_FREE );
        AllocEdgeList( &MinusAcceptord, EDGE_LIST_FREE );
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }

    if ( pc2i->len_c2at ) {
        /*------------------------------------------------------------------*/
        /* case 22: restored: N(-)=N(+)=C...=O orig: N#N-N=...-O(-)         */
        /*     im InChI        -O(-) may have H(+) added by Normalization   */
        /*                           or may be tautomeric                   */
        /* Solution: move (-) from N(-) to =O                               */
        /*                                                                  */
        /*------------------------------------------------------------------*/
        int num_DB_O = 0, iat;
        short iat_DB_O[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        int iN2, iC;
        BNS_EDGE *peDB_O_Minus;
        /*
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        CurrEdges.num_edges = 0; /* clear current edge list */
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: =O    */
                 num_DB_O < MAX_DIFF_FIXH &&
                 pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                 (pc2i->c2at[i].endptInChI || pc2i->c2at[i].nMobHInChI) &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                 /* reversed structure info: */
                 !(pc2i->c2at[i].endptRevrs || pc2i->c2at[i].nMobHRevrs) &&
                 pc2i->c2at[i].nFixHRevrs ==  0 &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                 at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 ) {
                iat_DB_O[num_DB_O ++] = iat;
            }
        }
        for ( i = 0; num_DB_O && i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            iat = nCanon2AtnoRevrs[i];
            if ( /* in restored atom O: charge=-1, no H, has no double bond, endpoint */
                 at2[iat].charge == -1 && !at2[iat].num_H &&
                 at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 && !pVA[iat].cMetal &&
                 pVA[iat].cNumValenceElectrons == 5 &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pBNS->edge[e].flow &&
                 !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint) &&
                 pVA[iN2=at2[iat].neighbor[0]].cNumValenceElectrons == 5 &&
                 at2[iat].bond_type[0] == BOND_TYPE_DOUBLE &&
                 at2[iN2].charge == 1 && at2[iN2].valence == 2 && at2[iN2].chem_bonds_valence == 4 &&
                 pVA[iC=at2[iN2].neighbor[at2[iN2].neighbor[0]==iN2]].cNumValenceElectrons == 4 ) {
                
                if ( ret = AddToEdgeList( &CurrEdges, e, INC_ADD_EDGE ) ) {
                    goto exit_function;
                }
            }
        }
        if ( num_try = inchi_min( CurrEdges.num_edges, num_DB_O ) ) {
            /* detected; attempt to fix */
            SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
            RemoveForbiddenEdgeMask( pBNS, &CurrEdges, forbidden_edge_mask  );
            delta = 1;
            for ( i = 0; i < num_DB_O && cur_success < num_try; i ++ ) {
                iat = iat_DB_O[i];

                peDB_O_Minus = pBNS->edge + (pVA[iat].nCMinusGroupEdge-1);
                pe           = pBNS->edge + pBNS->vert[iat].iedge[0];
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden            |= forbidden_edge_mask;
                peDB_O_Minus->forbidden  &= forbidden_edge_mask_inv;

                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) && nDeltaCharge == 0 ) {
                    /* Added (-) charge to =O and removed from =N(-) => nDeltaCharge == 0 */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        cur_success ++; /* 22 */
                    }
                } else {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
                INCHI_HEAPCHK
                pe->forbidden            &= forbidden_edge_mask_inv;
                peDB_O_Minus->forbidden  |= forbidden_edge_mask;
            }
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
        }
        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
    }
    if ( pc2i->len_c2at && pc2i->nNumTgInChI == 1 ) {
        /*------------------------------------------------------------------*/
        /* case 23: -NO2 are to be tautomeric but they are not AND          */
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
        int num_DB_O = 0, iat;
        short iat_DB_O[MAX_DIFF_FIXH], iat_NO2[MAX_DIFF_FIXH];
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        /*
        inp_ATOM *atfMobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at_fixed_bonds)?
                             pStruct->pOne_norm_data[1]->at_fixed_bonds : NULL;
        */
        S_CHAR   *num_Fixed_H_Revrs = pStruct->pOneINChI[0]->nNum_H_fixed? pStruct->pOneINChI[0]->nNum_H_fixed : NULL;
        S_CHAR   *pnMobHRevrs = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H)?
                               pStruct->pOneINChI[1]->nNum_H : 
                            (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H)?
                               pStruct->pOneINChI[0]->nNum_H : NULL;
        int iN, one_success;
        BNS_EDGE *peDB_O_Minus;
        int neigh, nNumO, nNumOthers;
#define CHG_SET_NOOH         0
#define CHG_SET_WRONG_TAUT   1
#define CHG_SET_TAUT         2
#define CHG_LAST_SET         2 /* the last index in trying */
#define CHG_SET_O_FIXED      3
#define CHG_SET_NUM          4
        EDGE_LIST ChangeableEdges[CHG_SET_NUM];
        memset( ChangeableEdges, 0, sizeof(ChangeableEdges) );
        /* equivalent to AllocEdgeList( &EdgeList, EDGE_LIST_CLEAR ); */
        /*
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        CurrEdges.num_edges = 0; /* clear current edge list */
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: taut in orig. InChI =O located in -NO2 that is not taut in Reconstructed InChI */
                 num_DB_O < MAX_DIFF_FIXH &&
                 pc2i->c2at[i].nValElectr == 6 /* O, S, Se, Te */ &&
                 (pc2i->c2at[i].endptInChI /*|| pc2i->c2at[i].nMobHInChI*/) &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == 0 && /*pc2i->c2at[i].nMobHInChI ==  1 &&*/
                 /* reversed structure info: */
                 !(pc2i->c2at[i].endptRevrs /*|| pc2i->c2at[i].nMobHRevrs*/) &&
                 pc2i->c2at[i].nFixHRevrs ==  0 &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                 at2[iat].valence == 1 && at2[iat].chem_bonds_valence == 2 &&
                 /* find whether it belongs to NO2 */
                 pVA[iN=at2[iat].neighbor[0]].cNumValenceElectrons == 5 &&
                 at2[iN].valence == 3 && (at2[iN].charge == 0 || at2[iN].charge == 1) &&
                 at2[iN].chem_bonds_valence == 5 - at2[iN].charge ) {
                /* find the second O */
                nNumO = nNumOthers = 0;
                for ( k = 0; k < at2[iN].valence; k ++ ) {
                    neigh = at2[iN].neighbor[k];
                    if ( neigh == iat ) {
                        continue;
                    }
                    if ( pVA[neigh].cNumValenceElectrons == 6 &&
                         pStruct->endpoint[neigh] &&
                         !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[neigh].endpoint) &&
                         at2[neigh].valence == 1 && at2[neigh].num_H == 0 &&
                         at2[neigh].radical == 0 && (at2[neigh].charge == 0 || at2[neigh].charge == -1) &&
                         at2[neigh].chem_bonds_valence - at2[neigh].charge == 2) {
                        nNumO ++;
                    } else
                    if ( at2[iN].bond_type[k] == BOND_TYPE_SINGLE &&
                         at2[neigh].valence > 1 &&
                         at2[neigh].valence < at2[neigh].chem_bonds_valence ) {
                        nNumOthers ++;
                    }
                }
                if ( nNumO != 1 || nNumOthers != 1 ) {
                    continue;
                }
                for ( k = 0; k < num_DB_O; k ++ ) {
                    if ( iat_NO2[k] == iN ) {
                        break;
                    }
                }
                if ( k == num_DB_O ) {
                    iat_NO2[num_DB_O]     = iN;
                    iat_DB_O[num_DB_O ++] = iat;
                }
                /* save the edge to avoid interference */
                if ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE ) ) {
                    goto exit_case_23;
                }
            }
        }
        if ( num_DB_O ) {
            /* 1. search for =N(=O)-OH; assume =N(+)(-O(-))(-OH) does not happen */
            for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
                /* find O=N(V) */
                iat = nCanon2AtnoRevrs[i];
                if ( !pStruct->endpoint[i] || pVA[i].cNumValenceElectrons != 6 ||
                     at2[iat].valence != 1 || at2[iat].charge ||
                     0 > (e = pVA[iat].nCMinusGroupEdge-1) ||
                     at2[iat].num_H + at2[iat].chem_bonds_valence != 2 ||
                     pVA[iN=at2[iat].neighbor[0]].cNumValenceElectrons != 5 ||
                     0 > (e = pVA[iN].nCPlusGroupEdge-1) ||
                     pBNS->edge[e].forbidden || !pBNS->edge[e].flow ||
                     at2[iN].charge || at2[iN].valence != 3 || at2[iN].chem_bonds_valence != 5) {
                    continue;
                }
                /* find the second O, -OH */
                nNumO = nNumOthers = 0;
                for ( k = 0; k < at2[iN].valence; k ++ ) {
                    neigh = at2[iN].neighbor[k];
                    if ( neigh == iat ) {
                        continue;
                    }
                    if ( pVA[neigh].cNumValenceElectrons == 6 &&
                         pStruct->endpoint[neigh] &&
                         at2[neigh].valence == 1 && at2[neigh].num_H == 1 &&
                         at2[neigh].radical == 0 && (at2[neigh].charge == 0 ) ) {
                        nNumO ++;
                    } else
                    if ( at2[iN].bond_type[k] == BOND_TYPE_DOUBLE &&
                         at2[neigh].valence >= 2 &&
                         at2[neigh].valence < at2[neigh].chem_bonds_valence ) {
                        nNumOthers ++;
                    }
                }
                if ( nNumO != 1 || nNumOthers != 1 ) {
                    continue;
                }
                /* save edges to be changed */
                if ( (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NOOH], e, INC_ADD_EDGE )) ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE ))) {
                    goto exit_case_23;
                }
                if ( NO_VERTEX != (j = GetChargeFlowerUpperEdge( pBNS, pVA, e )) &&
                     (( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NOOH], j, INC_ADD_EDGE ) ) ||
                      ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE ) ))) {
                    goto exit_case_23;
                }
            }
            /* 2. search for (-) atoms that are tautomeric but should not be  */
            /*           or that got H from Normalization but they shouldn't  */
            for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( at2[iat].charge == -1 &&
                     !pStruct->endpoint[i] &&
                     (at_Mobile_H_Revrs && 
                     (at_Mobile_H_Revrs[i].endpoint || at2[iat].num_H < at_Mobile_H_Revrs[i].num_H )) )  {

                    if ( 0 <= (e = pVA[iat].nCMinusGroupEdge-1) &&
                         0 > FindInEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e ) &&
                         !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                         (
                          ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_WRONG_TAUT], e, INC_ADD_EDGE ) ) ||
                          ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE ) )
                         ) ) {
                        goto exit_case_23;
                    }
                } else
                /* negatively charged atom in Reconstructed structure got H(+) from Normalization */
                /* and is not tautomeric; in the original structure it is tautomeric */    
                if ( at2[iat].charge == -1 &&
                     pStruct->endpoint[i] &&
                     !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint) &&
                     (num_Fixed_H_Revrs && num_Fixed_H_Revrs[i] == -1) &&
                     (pnMobHRevrs       && pnMobHRevrs[i]       ==  1) && 
                     pStruct->fixed_H[i] == 0 ) {

                    if ( 0 <= (e = pVA[iat].nCMinusGroupEdge-1) &&
                         0 > FindInEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e ) &&
                         !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                         (
                          ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_WRONG_TAUT], e, INC_ADD_EDGE ) ) ||
                          ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e, INC_ADD_EDGE ) )
                         ) ) {
                        goto exit_case_23;
                    }
                }
            }
            /* 3. Search for (-) atoms that are tautomeric */
            for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
                iat = nCanon2AtnoRevrs[i];
                if ( pStruct->endpoint[i] &&
                     (at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint) &&
                     at2[iat].charge == -1
                    /*&& pVA[i].cNumValenceElectrons == 6*/ ) {
                    if ( 0 <= (e = pVA[iat].nCMinusGroupEdge-1) &&
                         !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                         0 > FindInEdgeList( &ChangeableEdges[CHG_SET_O_FIXED], e )  &&
                         ( ret = AddToEdgeList( &ChangeableEdges[CHG_SET_TAUT], e, INC_ADD_EDGE ) ) ) {
                        goto exit_case_23;
                    }
                }
            }
            /* ------- finally, try to move charges from O=N --------------*/
            for ( i = 0; i < num_DB_O; i ++ ) {
                int nDeltaChargeExpected;
                one_success = 0;
                delta = 1;
                iat = iat_DB_O[i];
                peDB_O_Minus = pBNS->edge + (pVA[iat].nCMinusGroupEdge-1);
                pe =           pBNS->edge + pBNS->vert[iat].iedge[0];
                
                if ( !pe->flow )
                    continue;
                pv1 = pBNS->vert + (v1 = pe->neighbor1);
                pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                pe->forbidden            |= forbidden_edge_mask;

                pe->flow -= delta;
                pv1->st_edge.flow -= delta;
                pv2->st_edge.flow -= delta;
                pBNS->tot_st_flow -= 2*delta;

                for ( k = 0; !one_success && k <= CHG_LAST_SET; k ++ ) {
                    if ( !ChangeableEdges[k].num_edges ) {
                        continue;
                    }
                    nDeltaChargeExpected = (k==CHG_SET_NOOH)? 2 : 0;
                    
                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                    RemoveForbiddenEdgeMask( pBNS, &ChangeableEdges[k], forbidden_edge_mask  );
                    /* allow (-) charge to move to N=O */
                    peDB_O_Minus->forbidden  &= forbidden_edge_mask_inv;
                    
                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) &&
                                      nDeltaCharge == nDeltaChargeExpected ) {
                        /* Move (-) charge to =O and remove it an endpoint => nDeltaCharge == 0 */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if ( ret > 0 ) {
                            nNumRunBNS ++;
                            one_success ++; /* 23 */
                        }
                    }
                    INCHI_HEAPCHK
                }
                cur_success += one_success;
                
                RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                pe->forbidden            &= forbidden_edge_mask_inv;

                if ( !one_success ) {
                    pe->flow += delta;
                    pv1->st_edge.flow += delta;
                    pv2->st_edge.flow += delta;
                    pBNS->tot_st_flow += 2*delta;
                }
            }
        }
exit_case_23:
        for ( i = 0; i < CHG_SET_NUM; i ++ ) {
            AllocEdgeList( &ChangeableEdges[i], EDGE_LIST_FREE );
        }

        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
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

    if ( pc2i->len_c2at && pc2i->nNumTgInChI == 1 ) {
        /*------------------------------------------------------------------*/
        /* case 24: InChI norm. -N(-)-N(+)(IV) => -N=N(V) prevents tauto-   */
        /*          merism on -N(-)- in case of ADP                         */
        /*                                                                  */
        /* Solution: convert       N(V)=N-   ...=X    -> N(IV)(+)-N=...-X(-)*/
        /*                     N(IV)(+)-N(-)-...=X                          */
        /*                                                                  */
        /*      Orig InChI            taut      taut, 1 t-group only(ADP?)  */
        /*   Reconstructed struct   non-taut    possibly not taut           */
        /*                                                                  */
        /*   Details: 1a. store next to N(V) (+)edge its flower edge        */
        /*            1b. store next to N(-) edge NO_VERTEX                 */
        /*            2.  Release (-) edges of other missing endpoints or   */
        /*                all endpoints if no other is missing              */
        /*            3.  Decrement flow on (+) edge                        */
        /*                if flower edge is stored then expect DeltaCharge=2*/
        /*                otherwise DeltaCharge = 0                         */
        /*------------------------------------------------------------------*/
        int iat;
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        inp_ATOM *atf  = (pStruct->pOne_norm_data[1] && pStruct->pOne_norm_data[1]->at_fixed_bonds)?
                            pStruct->pOne_norm_data[1]->at_fixed_bonds : NULL;
        int iN, one_success;
        EdgeIndex  ef, e1;
        BNS_EDGE *pef;
#define CHG_SET_MISSED_TAUT   0
#define CHG_SET_OTHER_TAUT_O  1
#define CHG_SET_OTHER_TAUT_N  2
#define CHG_LAST_SET          2 /* the last index in trying */
#define CHG_SET_NN            3
#define CHG_SET_AVOID         4
#define CHG_SET_NUM           5
        EDGE_LIST ChangeableEdges[CHG_SET_NUM];
        memset( ChangeableEdges, 0, sizeof(ChangeableEdges) );
        /* equivalent to AllocEdgeList( &EdgeList, EDGE_LIST_CLEAR ); */
        /*
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        CurrEdges.num_edges = 0; /* clear current edge list */
        cur_success = 0;
        for ( i = 0; i < pc2i->len_c2at; i ++ ) {
            iat = pc2i->c2at[i].atomNumber;
            if ( /* orig. InChI info: -N=N(V)    */
                 pc2i->c2at[i].nValElectr == 5 /* N or P */ &&
                 (pc2i->c2at[i].endptInChI /* only N */) &&
                 (e1=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e1].forbidden &&
                 pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                 /* reversed structure info: */
                 !pc2i->c2at[i].endptRevrs &&
                 pc2i->c2at[i].nFixHRevrs ==  0 &&
                 pc2i->c2at[i].nAtChargeRevrs == 0 && !at2[iat].num_H &&
                 at2[iat].valence == 2 && at2[iat].chem_bonds_valence == 3 &&
                 /* find whether -N= has =N(V) neighbor; Note: operator comma: (A,B) returns B */
                 (iN = at2[iat].neighbor[at2[iat].bond_type[0] != BOND_TYPE_DOUBLE],
                  pVA[iN].cNumValenceElectrons == 5) &&
                 at2[iN].chem_bonds_valence == 5  && 
                 at2[iN].charge == 0 && !at2[iN].num_H && !at2[iN].radical &&
                 0 <= (e=pVA[iN].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                 0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e )) {

                ef = GetChargeFlowerUpperEdge( pBNS, pVA, e ); /* == NO_VERTEX if N(V) has 4 bonds */
                if ( (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], e, INC_ADD_EDGE ))    ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], ef, INC_ADD_EDGE ))   ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], 1, INC_ADD_EDGE ))    || /* expected nDeltaCharge */
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e1, INC_ADD_EDGE )) ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE )) ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], ef, INC_ADD_EDGE ))) {
                    goto exit_case_24;
                }
                /* mark -N= so that (-) will not be moved to it */
                if ( 0 <= (e = pVA[iat].nCMinusGroupEdge) && !pBNS->edge[e].forbidden &&
                     0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) &&
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ))) {
                    goto exit_case_24;
                }
            } else
            if ( /* orig. InChI info: -N(-)N(IV)(+)    */
                 atf &&
                 pc2i->c2at[i].nValElectr == 5 /* N or P */ &&
                 pc2i->c2at[i].endptInChI /* only N */ &&
                 (e=pVA[iat].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].forbidden &&
                 pc2i->c2at[i].nFixHInChI == 0 && pc2i->c2at[i].nMobHInChI == 0 &&
                 /* reversed structure info: */
                 !pc2i->c2at[i].endptRevrs &&
                 pc2i->c2at[i].nFixHRevrs ==  0 &&
                 pc2i->c2at[i].nAtChargeRevrs == -1 && !at2[iat].num_H &&
                 at2[iat].valence == 2 && at2[iat].chem_bonds_valence == 2 &&
                 atf[iat].valence == 2 && atf[iat].chem_bonds_valence == 3 && 
                 /* find whether -N= has =N(V) neighbor; Note: operator comma: (A,B) returns B */
                 (iN=atf[iat].neighbor[atf[iat].bond_type[0] != BOND_TYPE_DOUBLE],
                  pVA[iN].cNumValenceElectrons == 5) &&
                 at2[iN].charge == 1 && /* double bond neighbor */
                 at2[iN].chem_bonds_valence == 4 &&
                 atf[iN].charge == 0 &&
                 atf[iN].chem_bonds_valence == 5  &&  /* InChI normalization created N(V)=N- out of N(IV)(+)-N(-)- */
                 !at2[iN].num_H && !at2[iN].radical &&
                 0 <= (e=pVA[iat].nCMinusGroupEdge-1) && !pBNS->edge[e].forbidden && pBNS->edge[e].flow &&
                 0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) ) {
                /* save (-) edge */
                if ( (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], e, INC_ADD_EDGE ))    ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], NO_VERTEX, INC_ADD_EDGE ))   ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NN], 1, INC_ADD_EDGE ))   || /* expected nDeltaCharge */
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE ))) {
                    goto exit_case_24;
                }
            }
        }
        if ( !ChangeableEdges[CHG_SET_NN].num_edges ) {
            goto  exit_case_24;
        }
        /* Collect all relevant tautomeric atoms */
        for ( i = 0; i < pStruct->num_atoms; i ++ ) { /* i = canonical number - 1 */
            if ( !pStruct->endpoint[i] ) {
                continue;
            }
            iat = nCanon2AtnoRevrs[i];
            if ( at2[iat].charge || at2[iat].radical || at2[iat].valence == at2[iat].chem_bonds_valence ) {
                continue; /* cannot be an acceptor of (-) */
            }
            if ( 0 > (e=pVA[iat].nCMinusGroupEdge-1) || pBNS->edge[e].forbidden || pBNS->edge[e].flow ) {
                continue;
            }
            if ( 0 <= FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) ) {
                continue; /* has already been used */
            }
            /* missing endpoint */ 
            if ( !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[iat].endpoint) ) {
                if ( 0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) && (
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_MISSED_TAUT], e, INC_ADD_EDGE ))   ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE )))) {
                    goto exit_case_24;
                }
            } else
            /* endpoint O */
            if ( pVA[iat].cNumValenceElectrons == 6 ) {
                 if ( 0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) && (
                      (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_OTHER_TAUT_O], e, INC_ADD_EDGE ))   ||
                      (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE )))){
                    goto exit_case_24;
                 }
            } else
            /* endpoint N */
            if ( pVA[iat].cNumValenceElectrons == 5 ) {
                if ( 0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) && ( 
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_OTHER_TAUT_N], e, INC_ADD_EDGE ))   ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE )))){
                    goto exit_case_24;
                }
            }
        }
        /* ------- finally, try to move charges from -N(-)-N(+) or to N(V) --------------*/
        for ( i = 0; i < ChangeableEdges[CHG_SET_NN].num_edges; i += 3 ) {
            int nDeltaChargeExpected;
            one_success = 0;
            delta = 1;
            pe  = pBNS->edge + ChangeableEdges[CHG_SET_NN].pnEdges[i];
            pef = (NO_VERTEX != ChangeableEdges[CHG_SET_NN].pnEdges[i+1])?
                  pBNS->edge + ChangeableEdges[CHG_SET_NN].pnEdges[i+1] : NULL;
            nDeltaChargeExpected = ChangeableEdges[CHG_SET_NN].pnEdges[i+2];

            if ( !pe->flow )
                continue;
            pv1 = pBNS->vert + (v1 = pe->neighbor1);
            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

            pe->flow -= delta;
            pv1->st_edge.flow -= delta;
            pv2->st_edge.flow -= delta;
            pBNS->tot_st_flow -= 2*delta;

            for ( k = 0; !one_success && k <= CHG_LAST_SET; k ++ ) {
                if ( !ChangeableEdges[k].num_edges ) {
                    continue;
                }
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &ChangeableEdges[k], forbidden_edge_mask  );
                /* allow change of N(V) flower edge */
                if ( pef ) {
                    pef->forbidden  &= forbidden_edge_mask_inv;
                }
                
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) &&
                                  nDeltaCharge == nDeltaChargeExpected ) {
                    /* Move (-) charge to =O and remove it an endpoint => nDeltaCharge == 0 */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        one_success ++; /* 24 */
                    }
                }
                INCHI_HEAPCHK
            }
            cur_success += one_success;
            
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );

            if ( !one_success ) {
                pe->flow += delta;
                pv1->st_edge.flow += delta;
                pv2->st_edge.flow += delta;
                pBNS->tot_st_flow += 2*delta;
            }
        }
exit_case_24:
        for ( i = 0; i < CHG_SET_NUM; i ++ ) {
            AllocEdgeList( &ChangeableEdges[i], EDGE_LIST_FREE );
        }

        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
#undef CHG_SET_NN
#undef CHG_SET_MISSED_TAUT
#undef CHG_SET_OTHER_TAUT_O
#undef CHG_SET_OTHER_TAUT_N
#undef CHG_LAST_SET
#undef CHG_SET_AVOID
#undef CHG_SET_NUM
    }

    /* pStruct->nNumRemovedProtonsMobHInChI == pc2i->nNumRemHInChI */

    if ( pc2i->len_c2at && pc2i->nNumTgInChI == 1 &&
         pc2i->nNumRemHRevrs > pc2i->nNumRemHInChI && 0 > pc2i->nNumRemHInChI &&
         (pc2i->nNumEndpRevrs < pc2i->nNumEndpInChI ||
          pc2i->nNumTgRevrs > pc2i->nNumTgInChI ) ) {
        /*------------------------------------------------------------------*/
        /* case 25: Restored InChI does not have 2 or more added protons    */
        /*                         possibly taut. endpoints are missing     */
        /*                         has -N(-O(-))-O(-) group(s)              */
        /*          Original InChI has only one t-group                     */
        /*                                                                  */
        /* Solution: convert       -N(-O(-))-O(-) -> -N(+)(=O)-O(-)         */
        /*                         and direct 2(-) to the missing taut atoms*/
        /*           at first attempt try to move (-) to N only             */
        /*                                                                  */
        /*------------------------------------------------------------------*/
        int iat;
        AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
        AT_NUMB  *nAtno2CanonRevrs = pStruct->nAtno2Canon[0];
        inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                             pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
        /*
        inp_ATOM *atf  = (pStruct->pOne_norm_data[1] && pStruct->pOne_norm_data[1]->at_fixed_bonds)?
                            pStruct->pOne_norm_data[1]->at_fixed_bonds : NULL;
        */
        int iN, neigh, one_success;
        EdgeIndex  e1, bFirst;
        BNS_EDGE *pef;
#define CHG_SET_MISSED_TAUT_1   0
#define CHG_SET_MISSED_TAUT_ALL 1
#define CHG_SET_OTHER_TAUT_1    2
#define CHG_SET_OTHER_TAUT_ALL  3
#define CHG_LAST_SET            3 /* the last index in trying */
#define CHG_SET_NO_IN_NO2M2     4
#define CHG_SET_AVOID           5
#define CHG_SET_NUM             6
        EDGE_LIST ChangeableEdges[CHG_SET_NUM];
        memset( ChangeableEdges, 0, sizeof(ChangeableEdges) );
        /* equivalent to AllocEdgeList( &EdgeList, EDGE_LIST_CLEAR ); */
        /*
        S_CHAR   *nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H :
                               pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H : 0;
        */
        CurrEdges.num_edges = 0; /* clear current edge list */
        cur_success = 0;
        /* find all -N(-O(-))-O(-) */
        for ( i = 0; i < pStruct->num_atoms; i ++ ) {
            iat = nCanon2AtnoRevrs[i];
            if ( pStruct->endpoint[i] ) {
                if ( 0 > (e=pVA[iat].nCMinusGroupEdge-1) || pBNS->edge[e].forbidden ||
                     0 <= FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) ) {
                    continue;
                }
                bFirst = ( pVA[iat].cNumValenceElectrons == 5 && pc2i->nNumTgInChI == 1 ||
                           pVA[iat].cNumValenceElectrons == 6 && pc2i->nNumTgInChI != 1  );
                    /* many or no t-groups -> try O only first */
                    /* single t-group -> try only N first */
                if ( !(at_Mobile_H_Revrs && at_Mobile_H_Revrs[i].endpoint) ) {
                    /* missed tautomeric endpoint */
                    if ( bFirst && 
                         (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_MISSED_TAUT_1], e, INC_ADD_EDGE ))) {
                        goto exit_case_25;
                    }
                    if (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_MISSED_TAUT_ALL], e, INC_ADD_EDGE )) {
                        goto exit_case_25;
                    }
                }
                if ( bFirst && 
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_OTHER_TAUT_1], e, INC_ADD_EDGE ))) {
                    goto exit_case_25;
                }
                if (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_OTHER_TAUT_ALL], e, INC_ADD_EDGE )) {
                    goto exit_case_25;
                }
                if (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE )) {
                    goto exit_case_25;
                }
            } else
            if ( at2[iat].valence == 1 && at2[iat].charge == -1 &&
                 pVA[iat].cNumValenceElectrons == 6 &&
                 pVA[iN=at2[iat].neighbor[0]].cNumValenceElectrons == 5 && /* -O(-) */
                 !pStruct->endpoint[nAtno2CanonRevrs[iN]] &&
                 at2[iN].valence == 3 && at2[iN].chem_bonds_valence == 3 &&
                 !at2[iN].charge && !at2[iN].radical && 
                 0 <= (e=pVA[iN].nCPlusGroupEdge-1) && !pBNS->edge[e].forbidden &&
                 pBNS->edge[e].flow && /* NPlus edge */
                 0 <= (e1 = pVA[iat].nCMinusGroupEdge-1) && !pBNS->edge[e1].forbidden &&
                 pBNS->edge[e1].flow &&  /* OMinus edge */
                 0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e ) &&
                 0 > FindInEdgeList( &ChangeableEdges[CHG_SET_AVOID], e1 )) {
                /* found >N-O(-) */
                int nNumO = 0, nNumOthers = 0;
                for ( k = 0; k < at2[iN].valence; k ++ ) {
                    neigh = at2[iN].neighbor[k];
                    if ( neigh == iat ) {
                        continue;
                    }
                    if ( pVA[neigh].cNumValenceElectrons == 6 &&
                         !pStruct->endpoint[neigh] &&
                         at2[neigh].valence == 1 && at2[neigh].num_H == 0 &&
                         at2[neigh].radical == 0 && at2[neigh].charge == -1 &&
                         at2[neigh].chem_bonds_valence == 1 ) {
                        nNumO ++;
                    } else
                    if ( at2[iN].bond_type[k] == BOND_TYPE_SINGLE &&
                         at2[neigh].valence > 1 &&
                         at2[neigh].valence < at2[neigh].chem_bonds_valence ) {
                        nNumOthers ++;
                    }
                }
                if ( nNumO != 1 && nNumOthers != 1 ) {
                    continue;
                }
                /* save charge edges: NPlus first, OMinus second */
                if ( (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NO_IN_NO2M2], e, INC_ADD_EDGE )) ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_NO_IN_NO2M2], e1, INC_ADD_EDGE ))   ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e, INC_ADD_EDGE )) ||
                     (ret = AddToEdgeList( &ChangeableEdges[CHG_SET_AVOID], e1, INC_ADD_EDGE ))) {
                    goto exit_case_25;
                }
            }
        }
        if ( !ChangeableEdges[CHG_SET_NO_IN_NO2M2].num_edges ||
             !ChangeableEdges[CHG_SET_OTHER_TAUT_ALL].num_edges ) {
            goto exit_case_25;
        }
        /* ------- finally, try to move charges from -NO2(2-) or to tautomeric endpoints ----*/
        for ( i = 0; i < ChangeableEdges[CHG_SET_NO_IN_NO2M2].num_edges; i += 2 ) {
            int nDeltaChargeExpected = 3; 
            /* change flow on O(-) to make it neutral; 3 new charges will be created:
               N(+), and two (-) on InChI endpoints
               alternatively, if we change flow on N to make N(+) then O(-) will
               be nutralized (-1 charge) and two (-) charges on taut. endpoints will be
               created (+2); the total change in this case would be (-1)+(+2) = +1
            */
            one_success = 0;
            delta = 1;
            pe  = pBNS->edge + ChangeableEdges[CHG_SET_NO_IN_NO2M2].pnEdges[i+1]; /* O(-) edge */
            pef = pBNS->edge + ChangeableEdges[CHG_SET_NO_IN_NO2M2].pnEdges[i]; /* >N- (+) edge */

            if ( !pe->flow )
                continue;
            pv1 = pBNS->vert + (v1 = pe->neighbor1);
            pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

            pe->flow -= delta;
            pv1->st_edge.flow -= delta;
            pv2->st_edge.flow -= delta;
            pBNS->tot_st_flow -= 2*delta;

            for ( k = 0; !one_success && k <= CHG_LAST_SET; k ++ ) {
                if ( !ChangeableEdges[k].num_edges ) {
                    continue;
                }
                SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );
                RemoveForbiddenEdgeMask( pBNS, &ChangeableEdges[k], forbidden_edge_mask  );
                /* allow change of N(V) flower edge */
                pef->forbidden  &= forbidden_edge_mask_inv;
                
                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) &&
                                  nDeltaCharge == nDeltaChargeExpected ) {
                    /* Move (-) charge to =O and remove it an endpoint => nDeltaCharge == 0 */
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    if ( ret > 0 ) {
                        nNumRunBNS ++;
                        one_success ++; /* 24 */
                    }
                }
                INCHI_HEAPCHK
            }
            cur_success += one_success;
            
            RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask );

            if ( !one_success ) {
                pe->flow += delta;
                pv1->st_edge.flow += delta;
                pv2->st_edge.flow += delta;
                pBNS->tot_st_flow += 2*delta;
            }
        }
exit_case_25:
        for ( i = 0; i < CHG_SET_NUM; i ++ ) {
            AllocEdgeList( &ChangeableEdges[i], EDGE_LIST_FREE );
        }

        CurrEdges.num_edges = 0; /* clear current edge list */
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
            if ( !pInChI[0]->nNum_H_fixed && !pStruct->pOneINChI[0]->nNum_H_fixed ) {
                goto exit_function;  /* no fixed-H found */
            }
            if ( ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ) ) {
                goto exit_function;
            }
            if ( !pc2i->bHasDifference ) {
                goto exit_function; /* nothing to do */
            }
        }
#undef CHG_SET_NN
#undef CHG_SET_MISSED_TAUT
#undef CHG_SET_OTHER_TAUT_O
#undef CHG_SET_OTHER_TAUT_N
#undef CHG_LAST_SET
#undef CHG_SET_AVOID
#undef CHG_SET_NUM
    }


exit_function:
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &CurrEdges, EDGE_LIST_FREE );
    AllocEdgeList( &NFlowerEdges, EDGE_LIST_FREE );
    AllocEdgeList( &SFlowerEdges, EDGE_LIST_FREE );
    AllocEdgeList( &OtherNFlowerEdges, EDGE_LIST_FREE );
    AllocEdgeList( &FixedLargeRingStereoEdges, EDGE_LIST_FREE );
    AllocEdgeList( &AllBondEdges, EDGE_LIST_FREE );
    return ret < 0? ret : (pc2i->bHasDifference && tot_succes);
}
#endif
