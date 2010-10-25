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

/*^^^ */
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
/*^^^ */

/********************** Forbid carbon charge edges ***********************************/
int ForbidCarbonChargeEdges( BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups, EDGE_LIST *pCarbonChargeEdges, int forbidden_edge_mask  )
{
#define MAX_NUM_CARBON_CHARGE_EDGES 2
    int nType, i, k, ret;
    BNS_EDGE   *pEdge;
    if ( ret = AllocEdgeList( pCarbonChargeEdges, MAX_NUM_CARBON_CHARGE_EDGES ) ) {
        goto exit_function;
    }
    pCarbonChargeEdges->num_edges = 0;
    for ( i = 0; i < MAX_NUM_CARBON_CHARGE_EDGES; i ++ ) {
        switch( i ) {
        case 0:
            nType = TCG_Plus_C0;
            break;
        case 1:
            nType = TCG_Minus_C0;
            break;
        default:
            ret = RI_ERR_PROGR;
            goto exit_function;
        }
        if ( (k = pTCGroups->nGroup[nType]) >= 0 ) {
            k = pTCGroups->pTCG[k].nForwardEdge;
            if ( k > 0 ) {
                pEdge = pBNS->edge + k;
                if ( !(pEdge->forbidden & forbidden_edge_mask) ) {
                    pEdge->forbidden |= forbidden_edge_mask;
                    if ( ret = AddToEdgeList( pCarbonChargeEdges, k, 0 ) ) {
                        goto exit_function;
                    }
                }
            } else {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
        }
    }
    ret = pCarbonChargeEdges->num_edges;
exit_function:
    return ret;
#undef MAX_NUM_CARBON_CHARGE_EDGES
}
/******************************************************************************************************/
int ForbidNintrogenPlus2BondsInSmallRings( BN_STRUCT *pBNS, inp_ATOM *at, int num_at,
                                           VAL_AT *pVA, int min_ring_size, ALL_TC_GROUPS *pTCGroups,
                                           EDGE_LIST *pNplus2BondsEdges, int forbidden_edge_mask  )
{
    int i, j, ret;
    BNS_EDGE   *e;

    ret = 0;
        /* --- forbid edges that allow to make =N(+)= or #N(+)- in small ring */
    for ( i = 0; i < num_at; i ++ ) {
        if ( at[i].valence == 2 &&
             !at[i].num_H && !at[i].endpoint &&
             pVA[i].cNumValenceElectrons == 5 &&
             pVA[i].cPeriodicRowNumber == 1 &&
             !pVA[i].cMaxFlowToMetal && pVA[i].nCPlusGroupEdge > 0 &&
             pVA[i].cnListIndex > 0 && cnList[pVA[i].cnListIndex-1].bits == cn_bits_MNP &&
             pVA[i].cMinRingSize && pVA[i].cMinRingSize <= min_ring_size ) {

            e = pBNS->edge + (j = pVA[i].nCPlusGroupEdge - 1);
            if ( !(e->forbidden & forbidden_edge_mask) ) {
                e->forbidden |= forbidden_edge_mask;
                if ( ret = AddToEdgeList( pNplus2BondsEdges, j, 128 ) ) {
                    goto exit_function;
                }
            }
        }
    }
    ret = 0;
exit_function:
    return ret;
}

/*************************************************************************************************
Problem: Formula in InChI from the reversed structure has less H than in the input InChI
Solutions:

(a)   |                        |                 
     -B(-)-NH-=..-=N(+)<   => -B(-)-NH(+)=-..=-N<  (H is not removed from the ion pair)
      |                        |                 

                  |                      |   
(b)  >N(+)=-=...-=N-NH     =>  >N-=-...=-N(+)-NH  (charge from onium cannot be moved to remove H+)
                  |                      |   
*************************************************************************************************/
int FixLessHydrogenInFormula( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
                              inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask )
{
    int iBPlus=NO_VERTEX, iNV=NO_VERTEX, iNH = NO_VERTEX, neigh;
    EDGE_LIST NewlyFixedEdges;
    int ret, i, j;
    int num_at = pStruct->num_atoms;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    /* for RunBnsTestOnce */
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;

    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_CLEAR );
    if ( ret = AllocEdgeList( &NewlyFixedEdges, 2*num_at ) ) {
        goto exit_function;
    }
    for ( i = 0; i < num_at; i ++ ) {
        if ( (j = pVA[i].nCMinusGroupEdge-1) >= 0 ) {
            if ( ret = AddToEdgeList( &NewlyFixedEdges, j, 0 )) {
                goto exit_function;
            }
            pBNS->edge[j].forbidden |= forbidden_edge_mask;
        }
        if ( (j = pVA[i].nCPlusGroupEdge-1) >= 0 ) {
            if ( ret = AddToEdgeList( &NewlyFixedEdges, j, 0 )) {
                goto exit_function;
            }
            pBNS->edge[j].forbidden |= forbidden_edge_mask;
        }
    }
    /* extra H has been removed; check non-tautomeric atoms */
    for ( i = 0; i < num_at; i ++ ) {
        if ( !at2[i].endpoint && !pVA[i].cMetal &&
              pVA[i].cNumValenceElectrons == 5 && pVA[i].cPeriodicRowNumber == 1 &&
             at2[i].num_H == atf[i].num_H + 1) {
            /* H was removed from N */
            iNH = i;
            break;
        }
    }
    if ( 0 <= iNH && iNH < num_at ) {
        /* check neighbors for  |                 |
                          (a)  -B(+)-  or  (b)   =N-
                                |                 |
        */
        for ( j = 0; j < at2[i].valence; j ++ ) {
            neigh = at2[iNH].neighbor[j];
            if ( at2[neigh].valence == 4 ) {
                if ( at2[neigh].charge == -1 && at2[neigh].chem_bonds_valence == 4 &&
                     !at2[neigh].radical && !at[neigh].num_H ) {
                    iBPlus = neigh;
                }
            }
        }
    }
    if ( 0 <= iNH && iNH < num_at ) {
        int bond_type_at2;
        int bond_type_atf;
        int num_bonds_in_path = 0;
        int delta = -1, nxt = iNH, prv = NO_VERTEX, nxt_is_NPlus;
        /* the changed bond to the dehydrogenated atom H should have greater order */
        /* delta = (new bond order in atf[]) - (restored bond order in at2[]) */
        nxt_is_NPlus = 0;
        do {
            i = nxt;
            nxt = NO_VERTEX;
            delta = -delta;
            for ( j = 0; j < at2[i].valence; j ++ ) {
                bond_type_at2 = at2[i].bond_type[j] & BOND_TYPE_MASK; /* restored bond */
                bond_type_atf = atf[i].bond_type[j] & BOND_TYPE_MASK; /* normalized bond */
                nxt_is_NPlus  = 0;
                if ( (bond_type_atf - bond_type_at2 == delta || bond_type_atf == BOND_ALT12NS) &&
                     BOND_TYPE_SINGLE <= bond_type_at2 + delta && bond_type_at2 + delta <= BOND_TYPE_TRIPLE &&
                     !at2[(int)at2[i].neighbor[j]].cFlags ) {
                    prv = i;
                    nxt = at2[i].neighbor[j];
                    nxt_is_NPlus = at2[nxt].charge == 1 && atf[nxt].charge == 0 &&
                                   pVA[nxt].cNumValenceElectrons == 5 && pVA[nxt].cPeriodicRowNumber == 1;
                    at2[i].cFlags |= 1;  /* avoid cycling */
                    num_bonds_in_path ++;
                    if ( delta == -1 && at2[prv].valence == 4 && at2[prv].chem_bonds_valence == 5 &&
                         !at2[prv].charge && !at2[prv].radical && pVA[prv].cNumValenceElectrons == 5 &&
                         pVA[prv].nCPlusGroupEdge > 0 ) {
                        iNV = prv;
                    }
                    if ( at2[nxt].charge != atf[nxt].charge ) {
                        if ( (at2[nxt].charge == 1 || atf[nxt].charge == 1) &&
                              pVA[nxt].nCPlusGroupEdge > 0 ) {
                            pBNS->edge[pVA[nxt].nCPlusGroupEdge-1].forbidden &= inv_forbidden_edge_mask; 
                        }
                        if ( (at2[nxt].charge == -1 || atf[nxt].charge == -1) &&
                              pVA[nxt].nCMinusGroupEdge > 0 ) {
                            pBNS->edge[pVA[nxt].nCMinusGroupEdge-1].forbidden &= inv_forbidden_edge_mask; 
                        }
                    }
                    break; /* found */
                }
            }
        } while ( nxt >= 0 && !( nxt_is_NPlus && delta == -1 ) );
        for ( i = 0; i < num_at; i ++ ) {
            at2[i].cFlags = 0;
        }
        if ( nxt >= 0 && nxt_is_NPlus && delta == -1 ) {
            /* a simple alt path from NH-= to =N(+) has been found */
            if ( iBPlus || iNV ) {
                /* move (+) charge from N(+) to iNV or, if iBPlus, then to iNH */
                if ( iNV >= 0 && (j = pVA[iNV].nCPlusGroupEdge-1) > 0 && pBNS->edge[j].flow > 0 ||
                     iNH >= 0 && (j = pVA[iNH].nCPlusGroupEdge-1) > 0 && pBNS->edge[j].flow > 0 ) {
                    int          ieFlower;
                    BNS_EDGE    *pe  = pBNS->edge + j, *peFlower = NULL;
                    Vertex      v1   = pe->neighbor1;
                    Vertex      v2   = v1 ^ pe->neighbor12;
                    BNS_VERTEX  *pv1 = pBNS->vert + v1;
                    BNS_VERTEX  *pv2 = pBNS->vert + v2;

                    delta = 1;
                    /* prevent conversion of >N(+)= into N(V) neutral */
                    ieFlower = GetChargeFlowerUpperEdge( pBNS, pVA, pVA[nxt].nCPlusGroupEdge-1 );
                    if ( ieFlower >= 0 ) {
                        peFlower = pBNS->edge + ieFlower;
                        if ( peFlower->flow == delta ) {
                            peFlower->forbidden |= forbidden_edge_mask;
                            if ( ret = AddToEdgeList( &NewlyFixedEdges, ieFlower, 0 )) {
                                goto exit_function;
                            }
                        }
                    }
                    pe->forbidden |= forbidden_edge_mask;
                    pe->flow          -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2*delta;
                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                    if ( ret < 0 ) {
                        goto exit_function;
                    }
                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) &&
                                      nDeltaCharge <= 0  /* charge moving to this atom disappers*/ ) {
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        (*pnNumRunBNS) ++;
                        if ( ret < 0 ) {
                            goto exit_function;
                        } else
                        if ( ret == 1 ) {
                            *pnTotalDelta += ret;
                        } else {
                            ret = RI_ERR_PROGR;
                            goto exit_function;
                        }
                    } else {
                        ret = 0;
                        pe->flow          += delta;
                        pv1->st_edge.flow += delta;
                        pv2->st_edge.flow += delta;
                        pBNS->tot_st_flow += 2*delta;
                    }

                }

            }
        }
    }
exit_function:
    /* remove bond fixation */
    RemoveForbiddenEdgeMask( pBNS, &NewlyFixedEdges, forbidden_edge_mask );
    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_FREE );
    return ret;
}
/***********************************************************************************************


    X=Y-O(-)  => X(-)-Y=O


************************************************************************************************/
int FixMoreHydrogenInFormula( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
                              inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask )
{
    int iNH = NO_VERTEX, neigh, neigh2;
    EDGE_LIST NewlyFixedEdges;
    int ret, i, j, k, k2, delta;
    int num_at = pStruct->num_atoms;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    Vertex v1, v2;
    /* for RunBnsTestOnce */
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
    BNS_EDGE *pe, *pe2;

    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_CLEAR );
    if ( ret = AllocEdgeList( &NewlyFixedEdges, 2*num_at ) ) {
        goto exit_function;
    }
    /* fix all charges */
    for ( i = 0; i < num_at; i ++ ) {
        if ( (j = pVA[i].nCMinusGroupEdge-1) >= 0 ) {
            if ( ret = AddToEdgeList( &NewlyFixedEdges, j, 0 )) {
                goto exit_function;
            }
            pBNS->edge[j].forbidden |= forbidden_edge_mask;
        }
        if ( (j = pVA[i].nCPlusGroupEdge-1) >= 0 ) {
            if ( ret = AddToEdgeList( &NewlyFixedEdges, j, 0 )) {
                goto exit_function;
            }
            pBNS->edge[j].forbidden |= forbidden_edge_mask;
        }
    }

    /* H(+) has been added to -O(-); check non-tautomeric atoms */
    for ( i = 0; i < num_at; i ++ ) {
        if ( !(pStruct->bMobileH? at2[i].endpoint : pStruct->endpoint[i]) && !pVA[i].cMetal &&
             at2[i].num_H + 1 == atf[i].num_H &&      /* normalization added H ??? What would happen in Fixed-H case?*/
             (k = pVA[i].nCMinusGroupEdge-1) >= 0 &&
             pBNS->edge[k].flow == 1 &&               /* atom had (-) charge before preprocessing */
             at2[i].charge == -1 && atf[i].charge == 0 && /* and has no charge after preprocessing */
             at2[i].valence == 1 && at2[i].chem_bonds_valence == 1 && /* connected by a single bond */
             pVA[i].cNumValenceElectrons == 6 &&     /* atom is O, S, Se, Te */
             at2[neigh=at2[i].neighbor[0]].chem_bonds_valence > at2[neigh].valence
             /* atom's single neighbor has multiple bond(s)*/
            ) {
            /* H(+) was added to O in Y=X-O(-), where X is the only neighbor of O, X=neigh, Y=neigh2 */
            iNH = i;
            for ( j = 0; j < at2[neigh].valence; j ++ ) {
                neigh2 = at2[neigh].neighbor[j];
                if ( neigh2 != iNH && !at2[neigh2].endpoint &&
                     !pBNS->edge[(int)pBNS->vert[neigh].iedge[j]].forbidden &&
                     4 <= pVA[neigh2].cNumValenceElectrons &&
                          pVA[neigh2].cNumValenceElectrons <= 5 && /* neig2 is C or N */
                     (k2 = pVA[neigh2].nCMinusGroupEdge-1) >= 0 &&
                     !pBNS->edge[k2].flow /* negative charge may be moved to neigh2 */ ) {
                    break;
                }
            }
            if ( j < at2[neigh].valence ) {
                delta = 1;
                pe  = pBNS->edge + k;  /* -O(-) negative charge edge; flow = 1 */ 
                pe2 = pBNS->edge + k2; /* X charge edge; flow = 0 */
                v1  = pe->neighbor1;
                v2  = pe->neighbor12 ^ v1;
                pe->flow                    -= delta;
                pBNS->vert[v1].st_edge.flow -= delta;
                pBNS->vert[v2].st_edge.flow -= delta;
                pBNS->tot_st_flow           -= 2*delta;
                pe2->forbidden              &= inv_forbidden_edge_mask; /* allow the charge to move */

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                if ( ret < 0 ) {
                    goto exit_function;
                }
                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) &&
                                  nDeltaCharge <= 1 ) {
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    (*pnNumRunBNS) ++;
                    if ( ret < 0 ) {
                        goto exit_function;
                    } else
                    if ( ret ) {
                        *pnTotalDelta += ret;
                    } else {
                        ret = RI_ERR_PROGR;
                    }
                    break;
                } else {
                    /* the attempt has failed; restore the flow */
                    ret = 0;
                    pe->flow                    += delta;
                    pBNS->vert[v1].st_edge.flow += delta;
                    pBNS->vert[v2].st_edge.flow += delta;
                    pBNS->tot_st_flow           += 2*delta;
                }
            }
        }
    }
exit_function:
    /* remove bond fixation */
    RemoveForbiddenEdgeMask( pBNS, &NewlyFixedEdges, forbidden_edge_mask );
    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_FREE );
    return ret;
}
#if( FIX_ADD_PROTON_FOR_ADP == 1 )
/******************************************************************************************************/
int FixAddProtonForADP( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
                              inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, ICR *picr,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask )
{
    int iBPlus=NO_VERTEX, iNV=NO_VERTEX, iNH = NO_VERTEX, neigh, neigh2;
    EDGE_LIST NewlyFixedEdges;
    int ret, i, j, k, k2, delta;
    int num_at = pStruct->num_atoms;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    Vertex v1, v2;
    /* for RunBnsTestOnce */
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
    BNS_EDGE *pe, *pe2;

    ret = 0;
    /*
    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_CLEAR );
    
    for ( i = 0; i < num_at; i ++ ) {
        if ( at2[i].radical == RADICAL_DOUBLET && at2[i].endpoint ) {
            pStruct->bExtract |= EXTRACT_STRUCT_NUMBER;
            ret = 1;
            break;
        }
    }
    */
    return ret;
}
#endif
/******************************************************************************************************
      OH              OH
     /               /
  -NH      =>   -NH(+)        to eliminate false tautomerism. S(IV) or N(V) or P(V) may be a centerpoint
     \\              \
      O               O(-)
*******************************************************************************************************/
int FixRemoveExtraTautEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
                              inp_ATOM *at2, inp_ATOM *atf, inp_ATOM *atn, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, ICR *picr,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask )
{
    EDGE_LIST NewlyFixedEdges;
    int ret, i, j, k, delta, centerpoint, endpoint1, endpoint2;
    int num_at = pStruct->num_atoms;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    Vertex v1, v2;
    /* for RunBnsTestOnce */
    Vertex     vPathStart, vPathEnd;
    int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
    BNS_EDGE *pe, *pe2;

    ret = 0;

    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_CLEAR );
    if ( ret = AllocEdgeList( &NewlyFixedEdges, 2*num_at ) ) {
        goto exit_function;
    }
    /* fix all charges */
    for ( i = 0; i < num_at; i ++ ) {
        if ( (j = pVA[i].nCMinusGroupEdge-1) >= 0 ) {
            if ( ret = AddToEdgeList( &NewlyFixedEdges, j, 0 )) {
                goto exit_function;
            }
            pBNS->edge[j].forbidden |= forbidden_edge_mask;
        }
        if ( (j = pVA[i].nCPlusGroupEdge-1) >= 0 ) {
            if ( ret = AddToEdgeList( &NewlyFixedEdges, j, 0 )) {
                goto exit_function;
            }
            pBNS->edge[j].forbidden |= forbidden_edge_mask;
        }
    }

    for ( i = 0; i < picr->num_endp_in1_only; i ++ ) {
        endpoint1 = picr->endp_in1_only[i]-1;
        if ( at2[endpoint1].valence == at2[endpoint1].chem_bonds_valence ||
             pVA[endpoint1].nCMinusGroupEdge <= 0 ) {
            continue;
        }
        /* find centerpoint */
        for ( j = 0; j < at2[endpoint1].valence; j ++ ) {
            if ( BOND_TYPE_DOUBLE == ( BOND_TYPE_MASK & at2[endpoint1].bond_type[j] ) ) {
                centerpoint = at2[endpoint1].neighbor[j];
                if ( at2[centerpoint].charge || pVA[centerpoint].nCPlusGroupEdge <= 0 ||
                     !is_centerpoint_elem( at2[centerpoint].el_number ) ) {
                    continue;
                }
                /* -- the centerpoint as depicted has no ChargeStruct flower ---
                m = GetChargeFlowerUpperEdge( pBNS, pVA, pVA[centerpoint].nCPlusGroupEdge-1 );
                if ( m < 0 || pBNS->edge[m].flow ) {
                    continue;
                }
                */
                /* find 2nd endpoint */
                for ( k = 0; k < at2[centerpoint].valence; k ++ ) {
                    if ( BOND_TYPE_SINGLE != ( BOND_TYPE_MASK & at2[centerpoint].bond_type[k] ) ) {
                        continue;
                    }
                    endpoint2 = at2[centerpoint].neighbor[k];
                    if ( !at2[endpoint2].endpoint && atn[endpoint2].endpoint ) {
                        break;
                    }
                }
                if ( k == at2[centerpoint].valence ) {
                    continue;
                }
                /* the centerpoint and two extra endpoints have been found */
                pe = pBNS->edge + pVA[centerpoint].nCPlusGroupEdge - 1;
                if ( !pe->flow  ) {
                    continue;
                }
                pe2 = pBNS->edge + pVA[endpoint1].nCMinusGroupEdge - 1;
                if ( pe2->flow ) {
                    continue;
                }
                delta = 1;
                v1 = pe->neighbor1;
                v2 = pe->neighbor12 ^ v1;
                pe->flow                    -= delta;
                pBNS->vert[v1].st_edge.flow -= delta;
                pBNS->vert[v2].st_edge.flow -= delta;
                pBNS->tot_st_flow           -= 2*delta;
                pe2->forbidden              &= inv_forbidden_edge_mask; /* allow the charge to move */

                ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                      &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                if ( ret < 0 ) {
                    goto exit_function;
                }
                if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                  vPathEnd == v2 && vPathStart == v1) &&
                                  nDeltaCharge <= 1 ) {
                    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                    (*pnNumRunBNS) ++;
                    if ( ret < 0 ) {
                        goto exit_function;
                    } else
                    if ( ret ) {
                        *pnTotalDelta += ret;
                    } else {
                        ret = RI_ERR_PROGR;
                    }
                    goto exit_function;
                } else {
                    ret = 0;
                    pe->flow                    += delta;
                    pBNS->vert[v1].st_edge.flow += delta;
                    pBNS->vert[v2].st_edge.flow += delta;
                    pBNS->tot_st_flow           += 2*delta;
                    pe2->forbidden              |= forbidden_edge_mask;
                }
            }
        }
    }

exit_function:
    /* remove bond fixation */
    RemoveForbiddenEdgeMask( pBNS, &NewlyFixedEdges, forbidden_edge_mask );
    AllocEdgeList( &NewlyFixedEdges, EDGE_LIST_FREE );
    return ret;
}
/*************************************************************************/
int  FillOutExtraFixedHDataRestr( StrFromINChI *pStruct )
{
    int i, j, k, len, ret = 0;
    AT_NUMB *pNum;
    for ( i = 0; i < TAUT_NUM; i ++ ) {
        if ( pStruct->pOneINChI_Aux[i] ) {
            pNum = (pStruct->pOneINChI_Aux[i]->nIsotopicOrigAtNosInCanonOrd &&
                    pStruct->pOneINChI_Aux[i]->nIsotopicOrigAtNosInCanonOrd[0])?
                               pStruct->pOneINChI_Aux[i]->nIsotopicOrigAtNosInCanonOrd:
                   (pStruct->pOneINChI_Aux[i]->nOrigAtNosInCanonOrd &&
                    pStruct->pOneINChI_Aux[i]->nOrigAtNosInCanonOrd[0])?
                               pStruct->pOneINChI_Aux[i]->nOrigAtNosInCanonOrd : NULL;
        } else {
            pNum = NULL;
        }
        if ( pNum ) {
            len = pStruct->num_atoms * sizeof(pStruct->nCanon2Atno[0][0]);
            if ( !pStruct->nCanon2Atno[i] &&
                 !(pStruct->nCanon2Atno[i] = (AT_NUMB *)inchi_malloc( len )) ||
                 !pStruct->nAtno2Canon[i] &&
                 !(pStruct->nAtno2Canon[i] = (AT_NUMB *)inchi_malloc( len ))) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
            
            INCHI_HEAPCHK

            memcpy( pStruct->nCanon2Atno[i], pNum, len ); /* ??? the next for(...) fills it out */

            INCHI_HEAPCHK

            for ( j = 0; j < pStruct->num_atoms; j ++ ) {
                k = pNum[j]-1; /* atom number */
                pStruct->nCanon2Atno[i][j] = (AT_NUMB)k;
                pStruct->nAtno2Canon[i][k] = (AT_NUMB)j;
                INCHI_HEAPCHK
            }
        } else
        if ( !i ) {
            ret = RI_ERR_PROGR;
            goto exit_function;
        } else {
            if ( pStruct->nCanon2Atno[i] ) {
                inchi_free( pStruct->nCanon2Atno[i] );
                pStruct->nCanon2Atno[i] = NULL;
            }
            INCHI_HEAPCHK
            if ( pStruct->nAtno2Canon[i] ) {
                inchi_free( pStruct->nAtno2Canon[i] );
                pStruct->nAtno2Canon[i] = NULL;
            }
            INCHI_HEAPCHK
        }
    }

exit_function:
    return ret;
}
/*************************************************************************/
int  FillOutExtraFixedHDataInChI( StrFromINChI *pStruct, INChI *pInChI[] )
{
    int ret = 0;
    /*--- allocate memory for Mobile/Fixed-H data from the input InChI ---*/
    if ( NULL == pStruct->endpoint ) {
        pStruct->endpoint = (AT_NUMB *)inchi_calloc(pStruct->num_atoms, sizeof(pStruct->endpoint[0]));
    } else {
        memset( pStruct->endpoint, 0, pStruct->num_atoms * sizeof(pStruct->endpoint[0] ) );
    }
    if ( NULL == pStruct->fixed_H ) {
        pStruct->fixed_H = (S_CHAR *)inchi_malloc(pStruct->num_atoms * sizeof(pStruct->fixed_H[0]));
    }
    if ( !pStruct->endpoint || !pStruct->fixed_H ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    /*--- fill out Mobile/Fixed-H data from the input InChI ---*/
    GetTgroupInfoFromInChI( &pStruct->ti, NULL, pStruct->endpoint, pInChI[1] );
    if ( pInChI[0]->nNum_H_fixed ) {
        memcpy( pStruct->fixed_H, pInChI[0]->nNum_H_fixed, pStruct->num_atoms * sizeof(pStruct->fixed_H[0]) );
    } else {
        memset( pStruct->fixed_H, 0, pStruct->num_atoms * sizeof(pStruct->fixed_H[0]) );
    }

exit_function:
    return ret;
}
/***********************************************************************************************/
int FillOutCMP2FHINCHI( StrFromINChI *pStruct, inp_ATOM *at2, VAL_AT *pVA, INChI *pInChI[], CMP2FHINCHI *pc2i )
{
    int       ret = 0, i, j;
    int       bFixHRevrsExists  = pInChI[1] && pInChI[1]->nNumberOfAtoms > 0 && !pInChI[1]->bDeleted;
    inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[1] &&
                                 pStruct->pOne_norm_data[1]->at)? pStruct->pOne_norm_data[1]->at : NULL;
    S_CHAR   *num_Fixed_H_Revrs = pStruct->pOneINChI[0]->nNum_H_fixed? pStruct->pOneINChI[0]->nNum_H_fixed : NULL;
    /* atom number in structure that produced original InChI is atom number in all inp_ATOM *atoms */
    /* atom number in structure that produced restored InChI is in nAtomRevrs[]: */
    AT_NUMB  *nAtno2CanonRevrs = pStruct->nAtno2Canon[0];
    S_CHAR   *pnMobHInChI = (pInChI[1] && pInChI[1]->nNum_H)? pInChI[1]->nNum_H :
                            (pInChI[0] && pInChI[0]->nNum_H)? pInChI[0]->nNum_H : NULL;
    S_CHAR   *pnMobHRevrs = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H)?
                               pStruct->pOneINChI[1]->nNum_H : 
                            (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H)?
                               pStruct->pOneINChI[0]->nNum_H : NULL;
    int     nNumTgHInChI, nNumTgMInChI, nNumTgHRevrs, nNumTgMRevrs;
    memset( pc2i, 0, sizeof(*pc2i) );
    pc2i->nNumTgInChI = pStruct->ti.num_t_groups;
    pc2i->nNumTgRevrs = pStruct->One_ti.num_t_groups;
    pc2i->bHasDifference |= pc2i->nNumTgInChI != pc2i->nNumTgRevrs;

    pc2i->nNumRemHInChI = pStruct->nNumRemovedProtonsMobHInChI;
    pc2i->nNumRemHRevrs = pStruct->One_ti.tni.nNumRemovedProtons;
    pc2i->bHasDifference |= pc2i->nNumRemHInChI != pc2i->nNumRemHRevrs;

    pc2i->bFixedHLayerExistsRevrs = bFixHRevrsExists;
    pc2i->bHasDifference |= !bFixHRevrsExists;

    for ( i = 0; i < pStruct->ti.num_t_groups && i < pStruct->One_ti.num_t_groups; i ++ ) {
        nNumTgHInChI = pStruct->ti.t_group[i].num[0] - pStruct->ti.t_group[i].num[1];
        nNumTgMInChI = pStruct->ti.t_group[i].num[1];
        nNumTgHRevrs = pStruct->One_ti.t_group[i].num[0] - pStruct->One_ti.t_group[i].num[1];
        nNumTgMRevrs = pStruct->One_ti.t_group[i].num[1];

        pc2i->bHasDifference |= nNumTgHInChI != nNumTgHRevrs;
        pc2i->bHasDifference |= nNumTgMInChI != nNumTgMRevrs;

        if ( pStruct->ti.t_group[i].nNumEndpoints ==
             pStruct->One_ti.t_group[i].nNumEndpoints ) {
            
            if ( nNumTgHInChI != nNumTgHRevrs ) {
                pc2i->nNumTgDiffH ++;
            }
            if ( nNumTgMInChI != nNumTgMRevrs ) {
                pc2i->nNumTgDiffMinus ++;
            }
        }
        pc2i->bHasDifference |= pStruct->ti.t_group[i].nNumEndpoints !=
                                pStruct->One_ti.t_group[i].nNumEndpoints;

        pc2i->nNumTgHInChI += nNumTgHInChI;
        pc2i->nNumTgMInChI += nNumTgMInChI;
        pc2i->nNumTgHRevrs += nNumTgHRevrs;
        pc2i->nNumTgMRevrs += nNumTgMRevrs;

    }
    for ( ; i < pStruct->ti.num_t_groups; i ++ ) {
        nNumTgHInChI = pStruct->ti.t_group[i].num[0] - pStruct->ti.t_group[i].num[1];
        nNumTgMInChI = pStruct->ti.t_group[i].num[1];
        pc2i->nNumTgHInChI += nNumTgHInChI;
        pc2i->nNumTgMInChI += nNumTgMInChI;
        pc2i->bHasDifference |= 1;
    }
    for ( ; i < pStruct->One_ti.num_t_groups; i ++ ) {
        nNumTgHRevrs = pStruct->One_ti.t_group[i].num[0] - pStruct->One_ti.t_group[i].num[1];
        nNumTgMRevrs = pStruct->One_ti.t_group[i].num[1];
        pc2i->nNumTgHRevrs += nNumTgHRevrs;
        pc2i->nNumTgMRevrs += nNumTgMRevrs;
        pc2i->bHasDifference |= 1;

    }
    for ( i = j  = 0; i < pStruct->num_atoms; i ++ ) {
        /* i = original InChI canonical number - 1 */
        /* k = atom number from InChI created out of restored Fixed-H structure */
        int iCanonRevrs = nAtno2CanonRevrs[i];
        int endptInChI  = pStruct->endpoint[i]; /* endpoint in InChI */
        int endptRevrs  = at_Mobile_H_Revrs? at_Mobile_H_Revrs[i].endpoint : 0;
        int nFixHInChI  = pStruct->fixed_H[i];
        int nFixHRevrs  = num_Fixed_H_Revrs? num_Fixed_H_Revrs[iCanonRevrs]:0;
        int nMobHInChI  = pnMobHInChI? pnMobHInChI[i]:0;
        int nMobHRevrs =  pnMobHRevrs? pnMobHRevrs[iCanonRevrs]:0;
        if ( /*(!endptInChI || !endptRevrs) &&*/ (nFixHInChI != nFixHRevrs ) ||
             (!endptInChI != !endptRevrs) || nMobHInChI != nMobHRevrs ) {
             /* in InChI or reversed InChI atom[i] is not tautomeric */
             /* and number of fixed-H on the atom[i] differs */
            if ( j >= MAX_DIFF_FIXH ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            pc2i->c2at[j].endptInChI = endptInChI;
            pc2i->c2at[j].endptRevrs = endptRevrs;
            pc2i->bHasDifference |= !endptInChI != !endptRevrs;
            pc2i->c2at[j].atomNumber = i;
            pc2i->c2at[j].nValElectr = pVA[i].cNumValenceElectrons;
            pc2i->c2at[j].nPeriodNum = pVA[i].cPeriodicRowNumber;
            pc2i->c2at[j].nFixHInChI = nFixHInChI;
            pc2i->c2at[j].nFixHRevrs = nFixHRevrs;
            pc2i->bHasDifference |= nFixHInChI != nFixHRevrs;
            pc2i->c2at[j].nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H[i] :
                                       pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H[i] : 0;
            pc2i->c2at[j].nMobHRevrs = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H)?
                                           pStruct->pOneINChI[1]->nNum_H[iCanonRevrs] : 
                                       (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H)?
                                           pStruct->pOneINChI[0]->nNum_H[iCanonRevrs] : 0;
            pc2i->nNumDiffMobH   += (nMobHInChI != nMobHRevrs && !endptRevrs && !endptInChI);
            pc2i->bHasDifference |= nMobHInChI != nMobHRevrs;
            pc2i->c2at[j].nNumHRevrs = at2[i].num_H;
            pc2i->c2at[j].nAtChargeRevrs = at2[i].charge;
            j ++;
        }
        pc2i->nNumEndpInChI += (endptInChI != 0); 
        pc2i->nNumEndpRevrs += (endptRevrs != 0);

        if ( !pVA[i].cMetal ) {
            pc2i->nChargeFixHRevrsNonMetal += at2[i].charge;
            pc2i->nChargeMobHRevrsNonMetal += at_Mobile_H_Revrs? at_Mobile_H_Revrs[i].charge : 0;
        }

        /*pStruct->bExtract |= EXTRACT_STRUCT_NUMBER;*/
    }
    pc2i->nChargeFixHInChI = pInChI[0]? pInChI[0]->nTotalCharge : 0;
    pc2i->nChargeMobHInChI = pInChI[1]? pInChI[1]->nTotalCharge : 0;

    pc2i->nChargeMobHRevrs = pStruct->pOneINChI[1]? pStruct->pOneINChI[1]->nTotalCharge : 
                             pStruct->pOneINChI[0]? pStruct->pOneINChI[0]->nTotalCharge : 0;
    pc2i->nChargeFixHRevrs = pStruct->pOneINChI[0]? pStruct->pOneINChI[0]->nTotalCharge : 0;
    
    pc2i->bHasDifference |= pc2i->nChargeFixHInChI != pc2i->nChargeFixHRevrs;
    pc2i->bHasDifference |= pc2i->nChargeMobHInChI != pc2i->nChargeMobHRevrs;

exit_function:
    pc2i->len_c2at = j;

    return ret;
}
/***********************************************************************************************/
int FillOutCMP2MHINCHI( StrFromINChI *pStruct, ALL_TC_GROUPS *pTCGroups, inp_ATOM *at2,
                        VAL_AT *pVA, INChI *pInChI[], CMP2MHINCHI *pc2i )
{
    int       ret = 0, i, j, iat;
    int       bFixHRevrsExists  = pInChI[1] && pInChI[1]->nNumberOfAtoms > 0 && !pInChI[1]->bDeleted;
    inp_ATOM *at_Mobile_H_Revrs = (pStruct->pOne_norm_data[0] &&
                                 pStruct->pOne_norm_data[0]->at)? pStruct->pOne_norm_data[0]->at : NULL;
    /* atom number in structure that produced original InChI is atom number in all inp_ATOM *atoms */
    /* atom number in structure that produced restored InChI is in nAtomRevrs[]: */
    AT_NUMB  *nCanon2AtnoRevrs = pStruct->nCanon2Atno[0];
    AT_NUMB  *nAtno2CanonRevrs = pStruct->nAtno2Canon[0];
    S_CHAR   *pnMobHInChI = (pInChI[0] && pInChI[0]->nNum_H)? pInChI[0]->nNum_H : NULL;
    S_CHAR   *pnMobHRevrs = (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H)?
                               pStruct->pOneINChI[0]->nNum_H : NULL;
    int     nNumTgHInChI, nNumTgMInChI, nNumTgHRevrs, nNumTgMRevrs;
    memset( pc2i, 0, sizeof(*pc2i) );
    pc2i->nNumTgInChI = pStruct->ti.num_t_groups;
    pc2i->nNumTgRevrs = pStruct->One_ti.num_t_groups;
    pc2i->bHasDifference |= pc2i->nNumTgInChI != pc2i->nNumTgRevrs;

    pc2i->nNumRemHInChI = pStruct->nNumRemovedProtonsMobHInChI;
    pc2i->nNumRemHRevrs = pStruct->One_ti.tni.nNumRemovedProtons;
    /*pc2i->bHasDifference |= pc2i->nNumRemHInChI != pc2i->nNumRemHRevrs;*/

    pc2i->bFixedHLayerExistsRevrs = bFixHRevrsExists;
    /*pc2i->bHasDifference |= !bFixHRevrsExists;*/

    for ( i = 0; i < pStruct->ti.num_t_groups; i ++ ) {
        int jFst = pStruct->ti.t_group[i].nFirstEndpointAtNoPos;
        int jNum = pStruct->ti.t_group[i].nNumEndpoints;
        int is_N, is_O; 
        for ( j = 0; j < jNum; j ++ ) {
            iat = pStruct->ti.nEndpointAtomNumber[jFst + j];
            is_N = pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1;
            is_O = pVA[iat].cNumValenceElectrons == 6;
            if ( is_N + is_O != 1 ) {
                return RI_ERR_SYNTAX;
            }
            pc2i->nNumTgNInChI += is_N;
            pc2i->nNumTgOInChI += is_O;
            if ( at2[iat].chem_bonds_valence == at2[iat].valence ) {
                /* donor */
                if ( is_N ) {
                    /* N */
                    pc2i->nNumTgNHInChI      += at2[iat].charge ==  0 && at2[iat].num_H == 1;
                    pc2i->nNumTgNH2InChI     += at2[iat].charge ==  0 && at2[iat].num_H == 2;
                    pc2i->nNumTgNMinusInChI  += at2[iat].charge == -1 && at2[iat].num_H == 0;
                    pc2i->nNumTgNHMinusInChI += at2[iat].charge == -1 && at2[iat].num_H == 1;
                } else {
                    /* O, S, Se, Te */
                    pc2i->nNumTgOHInChI      += at2[iat].charge ==  0 && at2[iat].num_H == 1;
                    pc2i->nNumTgOMinusInChI  += at2[iat].charge == -1 && at2[iat].num_H == 0;
                }
            } else
            if ( at2[iat].chem_bonds_valence == at2[iat].valence+1 ) {
                /* donor */
                if ( is_N ) {
                    /* N */
                    pc2i->nNumTgDBNHInChI      += at2[iat].charge ==  0 && at2[iat].num_H == 1;
                    pc2i->nNumTgDBNMinusInChI  += at2[iat].charge == -1 && at2[iat].num_H == 0;
                    pc2i->nNumTgDBNInChI       += at2[iat].charge ==  0 && at2[iat].num_H == 0;
                } else {
                    /* O, S, Se, Te */
                    pc2i->nNumTgDBOInChI       += at2[iat].charge ==  0 && at2[iat].num_H == 0;
                }
            }
        }
    }
    for ( i = 0; i < pStruct->One_ti.num_t_groups; i ++ ) {
        int jFst = pStruct->One_ti.t_group[i].nFirstEndpointAtNoPos;
        int jNum = pStruct->One_ti.t_group[i].nNumEndpoints;
        int is_N, is_O; 
        for ( j = 0; j < jNum; j ++ ) {
            iat = nCanon2AtnoRevrs[(int)pStruct->One_ti.nEndpointAtomNumber[jFst + j]];
            is_N = pVA[iat].cNumValenceElectrons == 5 && pVA[iat].cPeriodicRowNumber == 1;
            is_O = pVA[iat].cNumValenceElectrons == 6;
            if ( is_N + is_O != 1 ) {
                return RI_ERR_PROGR;
            }
            pc2i->nNumTgNRevrs += is_N;
            pc2i->nNumTgORevrs += is_O;
            if ( at2[iat].chem_bonds_valence == at2[iat].valence ) {
                /* donor */
                if ( is_N ) {
                    /* N */
                    pc2i->nNumTgNHRevrs      += at2[iat].charge ==  0 && at2[iat].num_H == 1;
                    pc2i->nNumTgNH2Revrs     += at2[iat].charge ==  0 && at2[iat].num_H == 2;
                    pc2i->nNumTgNMinusRevrs  += at2[iat].charge == -1 && at2[iat].num_H == 0;
                    pc2i->nNumTgNHMinusRevrs += at2[iat].charge == -1 && at2[iat].num_H == 1;
                } else {
                    /* O, S, Se, Te */
                    pc2i->nNumTgOHRevrs      += at2[iat].charge ==  0 && at2[iat].num_H == 1;
                    pc2i->nNumTgOMinusRevrs  += at2[iat].charge == -1 && at2[iat].num_H == 0;
                }
            } else
            if ( at2[iat].chem_bonds_valence == at2[iat].valence+1 ) {
                /* donor */
                if ( is_N ) {
                    /* N */
                    pc2i->nNumTgDBNHRevrs      += at2[iat].charge ==  0 && at2[iat].num_H == 1;
                    pc2i->nNumTgDBNMinusRevrs  += at2[iat].charge == -1 && at2[iat].num_H == 0;
                    pc2i->nNumTgDBNRevrs       += at2[iat].charge ==  0 && at2[iat].num_H == 0;
                } else {
                    /* O, S, Se, Te */
                    pc2i->nNumTgDBORevrs       += at2[iat].charge ==  0 && at2[iat].num_H == 0;
                }
            }
        }
    }

    for ( i = 0; i < pStruct->ti.num_t_groups && i < pStruct->One_ti.num_t_groups; i ++ ) {
        nNumTgHInChI = pStruct->ti.t_group[i].num[0] - pStruct->ti.t_group[i].num[1];
        nNumTgMInChI = pStruct->ti.t_group[i].num[1];
        nNumTgHRevrs = pStruct->One_ti.t_group[i].num[0] - pStruct->One_ti.t_group[i].num[1];
        nNumTgMRevrs = pStruct->One_ti.t_group[i].num[1];

        pc2i->bHasDifference |= nNumTgHInChI != nNumTgHRevrs;
        pc2i->bHasDifference |= nNumTgMInChI != nNumTgMRevrs;

        if ( pStruct->ti.t_group[i].nNumEndpoints ==
             pStruct->One_ti.t_group[i].nNumEndpoints ) {
            
            if ( nNumTgHInChI != nNumTgHRevrs ) {
                pc2i->nNumTgDiffH ++;
            }
            if ( nNumTgMInChI != nNumTgMRevrs ) {
                pc2i->nNumTgDiffMinus ++;
            }
        }
        pc2i->bHasDifference |= pStruct->ti.t_group[i].nNumEndpoints !=
                                pStruct->One_ti.t_group[i].nNumEndpoints;

        pc2i->nNumTgHInChI += nNumTgHInChI;
        pc2i->nNumTgMInChI += nNumTgMInChI;
        pc2i->nNumTgHRevrs += nNumTgHRevrs;
        pc2i->nNumTgMRevrs += nNumTgMRevrs;

    }
    for ( ; i < pStruct->ti.num_t_groups; i ++ ) {
        nNumTgHInChI = pStruct->ti.t_group[i].num[0] - pStruct->ti.t_group[i].num[1];
        nNumTgMInChI = pStruct->ti.t_group[i].num[1];
        pc2i->nNumTgHInChI += nNumTgHInChI;
        pc2i->nNumTgMInChI += nNumTgMInChI;
        pc2i->bHasDifference |= 1;
    }
    for ( ; i < pStruct->One_ti.num_t_groups; i ++ ) {
        nNumTgHRevrs = pStruct->One_ti.t_group[i].num[0] - pStruct->One_ti.t_group[i].num[1];
        nNumTgMRevrs = pStruct->One_ti.t_group[i].num[1];
        pc2i->nNumTgHRevrs += nNumTgHRevrs;
        pc2i->nNumTgMRevrs += nNumTgMRevrs;
        pc2i->bHasDifference |= 1;

    }
    for ( i = j  = 0; i < pStruct->num_atoms; i ++ ) {
        /* i = original InChI canonical number - 1 */
        /* k = atom number from InChI created out of restored Fixed-H structure */
        int iCanonRevrs = nAtno2CanonRevrs[i];
        int endptInChI  = at2[i].endpoint; /* endpoint in InChI */
        int endptRevrs  = at_Mobile_H_Revrs? at_Mobile_H_Revrs[i].endpoint : 0;
        int nMobHInChI  = pnMobHInChI? pnMobHInChI[i]:0;
        int nMobHRevrs =  pnMobHRevrs? pnMobHRevrs[iCanonRevrs]:0;
        if ( (!endptInChI != !endptRevrs) || nMobHInChI != nMobHRevrs ) {
             /* in InChI or reversed InChI atom[i] is not tautomeric */
             /* and number of fixed-H on the atom[i] differs */
            if ( j >= MAX_DIFF_FIXH ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            pc2i->c2at[j].endptInChI = endptInChI;
            pc2i->c2at[j].endptRevrs = endptRevrs;
            pc2i->bHasDifference |= !endptInChI != !endptRevrs;
            pc2i->c2at[j].atomNumber = i;
            pc2i->c2at[j].nValElectr = pVA[i].cNumValenceElectrons;
            pc2i->c2at[j].nPeriodNum = pVA[i].cPeriodicRowNumber;
            pc2i->c2at[j].nMobHInChI = pInChI[1] && pInChI[1]->nNum_H? pInChI[1]->nNum_H[i] :
                                       pInChI[0] && pInChI[0]->nNum_H? pInChI[0]->nNum_H[i] : 0;
            pc2i->c2at[j].nMobHRevrs = (pStruct->pOneINChI[1] && pStruct->pOneINChI[1]->nNum_H)?
                                           pStruct->pOneINChI[1]->nNum_H[iCanonRevrs] : 
                                       (pStruct->pOneINChI[0] && pStruct->pOneINChI[0]->nNum_H)?
                                           pStruct->pOneINChI[0]->nNum_H[iCanonRevrs] : 0;
            
            pc2i->nNumDiffMobH   += (nMobHInChI != nMobHRevrs && !endptRevrs && !endptInChI);
            pc2i->bHasDifference |= (nMobHInChI != nMobHRevrs);
            pc2i->c2at[j].nNumHRevrs = at2[i].num_H;
            pc2i->c2at[j].nAtChargeRevrs = at2[i].charge;
            j ++;
        }
        pc2i->nNumEndpInChI += (endptInChI != 0); 
        pc2i->nNumEndpRevrs += (endptRevrs != 0);

        if ( !pVA[i].cMetal ) {
            pc2i->nChargeMobHRevrsNonMetal += (at_Mobile_H_Revrs && !at_Mobile_H_Revrs[i].endpoint)? at_Mobile_H_Revrs[i].charge : 0;
        }


        /*pStruct->bExtract |= EXTRACT_STRUCT_NUMBER;*/
    }
    pc2i->nChargeMobHRevrsNonMetal += pTCGroups->tgroup_charge;

    pc2i->nChargeMobHInChI = pInChI[0]? pInChI[0]->nTotalCharge : 0;

    pc2i->nChargeMobHRevrs = pStruct->pOneINChI[0]? pStruct->pOneINChI[0]->nTotalCharge : 0;
    
    pc2i->bHasDifference |= pc2i->nChargeMobHInChI != pc2i->nChargeMobHRevrs;

exit_function:
    pc2i->len_c2at = j;

    return ret;
}
/******************************************************************************************************/
int NormalizeAndCompare(ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
                        StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA,
                        ALL_TC_GROUPS *pTCGroups, INChI *pInChI[], long num_inp, int bHasSomeFixedH,
                        int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask, int forbidden_stereo_edge_mask)
{
    int i;
    int err;
    ICR icr, icr2;
    int num_norm_endpoints, num_endpoints, num_norm_t_groups, num_mobile, num_norm_mobile, ret = 0;
#if( bRELEASE_VERSION == 0 )
#ifndef INCHI_LIBRARY
    const char *szCurHdr = (ip->pSdfValue && ip->pSdfValue[0])? ip->pSdfValue : "???";
    int         iComponent = pTCGroups->iComponent;
#endif
#endif
    T_GROUP_INFO *t_group_info = NULL;
    inp_ATOM     *at_norm = NULL; /* normalized */
    inp_ATOM     *at_prep = NULL; /* preprocessed */
    INCHI_MODE  cmpInChI, cmpInChI2;
    int         nDeltaPrev, nDeltaCur;
    int         iOrigInChI, iRevrInChI;


    /***********************************************************/
    /* normalize and create one component InChI                */
    /***********************************************************/
    ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                          &t_group_info, &at_norm, &at_prep );
    if ( ret < 0 ) {
#if( bRELEASE_VERSION == 0 )
#ifndef INCHI_LIBRARY
        fprintf( stdout, "\nERROR in MakeOneInchi-1: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                 szCurHdr? szCurHdr: "???", iComponent, pStruct->iInchiRec? 'R':'D', pStruct->iMobileH?'M':'F', ret);
#endif
#endif
        goto exit_function;
    }
    if ( pStruct->bMobileH == TAUT_NON ) {
        /* these indexes are used to compare Mobile-H InChI */
        iOrigInChI = (pInChI[1] && pInChI[1]->nNumberOfAtoms && !pInChI[1]->bDeleted)? 1 : 0;
        iRevrInChI = (pStruct->pOneINChI[1] &&pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted)? 1 : 0;
    } else {
        iOrigInChI = 0;
        iRevrInChI = 0;
    }

    /************************************************************/
    /* compare                                                  */
    /************************************************************/
    if ( pStruct->iMobileH == TAUT_NON && (ret = FillOutExtraFixedHDataRestr( pStruct )) ) {
        goto exit_function;
    }
    cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *a2*/, &icr, &err );
    if ( cmpInChI & IDIF_PROBLEM ) {
        ret = RI_ERR_PROGR; /* severe restore problem */
        goto exit_function;
    }
    if ( err ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    /********** InChI from restored structure has LESS hydrogen atoms ******************************/
    if ( (cmpInChI & IDIF_LESS_H) && at_prep && 0 < (nDeltaCur = icr.tot_num_H2 - icr.tot_num_H1) ) {
        do {
            ret =  FixLessHydrogenInFormula( pBNS, pBD, pStruct, at, at2, at_prep, pVA, pTCGroups,
                                             pnNumRunBNS, pnTotalDelta, forbidden_edge_mask );
            if ( ret < 0 ) {
                goto exit_function;
            }
            if ( ret ) {
                /* Probably success. The changes are in pBNS. Create new InChI out of the new restored structure */
                ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                      &t_group_info, &at_norm, &at_prep );
                if ( ret < 0 ) {
#if( bRELEASE_VERSION == 0 )
#ifndef INCHI_LIBRARY
                    fprintf( stdout, "\nERROR in MakeOneInchi-2: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                             szCurHdr? szCurHdr: "???", iComponent, pStruct->iInchiRec? 'R':'D', pStruct->iMobileH?'M':'F', ret);
#endif
#endif
                    goto exit_function;
                }
                /* compare new InChI to the original InChI */
                if ( pStruct->bMobileH == TAUT_NON ) {
                    iRevrInChI = (pStruct->pOneINChI[1] &&pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted)? 1 : 0;
                } else {
                    iRevrInChI = 0;
                }
                if ( pStruct->iMobileH == TAUT_NON && (ret = FillOutExtraFixedHDataRestr( pStruct )) ) {
                    goto exit_function;
                }
                cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL, &icr, &err );
                nDeltaPrev = nDeltaCur;
                nDeltaCur  = icr.tot_num_H2 - icr.tot_num_H1;
            } else {
                break;
            }
        } while( (cmpInChI & IDIF_LESS_H) && at_prep && nDeltaCur && nDeltaCur < nDeltaPrev );
    }
    /********** InChI from restored structure has MORE hydrogen atoms ******************************/
    if ( (cmpInChI & IDIF_MORE_H) && at_prep && 0 < (nDeltaCur = icr.tot_num_H1 - icr.tot_num_H2) ) {
        do {
            ret =  FixMoreHydrogenInFormula( pBNS, pBD, pStruct, at, at2, at_prep, pVA, pTCGroups,
                                             pnNumRunBNS, pnTotalDelta, forbidden_edge_mask );
            if ( ret < 0 ) {
                goto exit_function;
            }
            if ( ret ) {
                /* Probably success. The changes are in pBNS. Create new InChI out of the new restored structure */
                ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                      &t_group_info, &at_norm, &at_prep );
                if ( ret < 0 ) {
#if( bRELEASE_VERSION == 0 )
#ifndef INCHI_LIBRARY
                    fprintf( stdout, "\nERROR in MakeOneInchi-3: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                             szCurHdr? szCurHdr: "???", iComponent, pStruct->iInchiRec? 'R':'D', pStruct->iMobileH?'M':'F', ret);
#endif
#endif
                    goto exit_function;
                }
                /* compare new InChI to the original InChI */
                if ( pStruct->bMobileH == TAUT_NON ) {
                    iRevrInChI = (pStruct->pOneINChI[1] &&pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted)? 1 : 0;
                } else {
                    iRevrInChI = 0;
                }
                if ( pStruct->iMobileH == TAUT_NON && (ret = FillOutExtraFixedHDataRestr( pStruct )) ) {
                    goto exit_function;
                }
                cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL, &icr, &err );
                nDeltaPrev = nDeltaCur;
                nDeltaCur  = icr.tot_num_H1 - icr.tot_num_H2;
            } else {
                break;
            }
        } while( (cmpInChI & IDIF_MORE_H) && at_prep && nDeltaCur && nDeltaCur < nDeltaPrev );
    }
    /***************** Fix non-taut atoms normalized to tautomeric endpoints ***********************/
    if ( (cmpInChI & IDIF_EXTRA_TG_ENDP) && at_norm && 0 < (nDeltaCur = icr.num_endp_in1_only) ) {
        do {
            ret = FixRemoveExtraTautEndpoints( pBNS, pBD, pStruct, at, at2, at_prep, at_norm, pVA, pTCGroups, &icr,
                                               pnNumRunBNS, pnTotalDelta, forbidden_edge_mask );
            if ( ret < 0 ) {
                goto exit_function;
            }
            if ( ret ) {
                /* Probably success. The changes are in pBNS. Create new InChI out of the new restored structure */
                ret = MakeOneInChIOutOfStrFromINChI2( ip, sd, pBNS, pStruct, at, at2, at3, pVA, pTCGroups,
                                                      &t_group_info, &at_norm, &at_prep );
                if ( ret < 0 ) {
#if( bRELEASE_VERSION == 0 )
#ifndef INCHI_LIBRARY
                    fprintf( stdout, "\nERROR in MakeOneInchi-4: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                             szCurHdr? szCurHdr: "???", iComponent, pStruct->iInchiRec? 'R':'D', pStruct->iMobileH?'M':'F', ret);
#endif
#endif
                    goto exit_function;
                }
                /* compare new InChI to the original InChI */
                if ( pStruct->bMobileH == TAUT_NON ) {
                    iRevrInChI = (pStruct->pOneINChI[1] &&pStruct->pOneINChI[1]->nNumberOfAtoms && !pStruct->pOneINChI[1]->bDeleted)? 1 : 0;
                } else {
                    iRevrInChI = 0;
                }
                if ( pStruct->iMobileH == TAUT_NON && (ret = FillOutExtraFixedHDataRestr( pStruct )) ) {
                    goto exit_function;
                }
                cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL, &icr, &err );
                nDeltaPrev = nDeltaCur;
                nDeltaCur  = icr.num_endp_in1_only;
            } else {
                break;
            }
        } while( (cmpInChI & IDIF_EXTRA_TG_ENDP) && at_norm && nDeltaCur && nDeltaCur < nDeltaPrev );
    }
    /************************ case of Fixed-H ******************************************************/

    if ( pStruct->bMobileH == TAUT_NON ) {
        int num_tries = 0;
        do {
            if ( 0 > (ret = FixFixedHRestoredStructure(ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA, pTCGroups,
                                                  &t_group_info, &at_norm, &at_prep, pInChI,
                                                  num_inp, bHasSomeFixedH, pnNumRunBNS, pnTotalDelta, forbidden_edge_mask,
                                                  forbidden_stereo_edge_mask) ) ) {
                goto exit_function;
            }
        } while( num_tries ++ < 2 && ret > 0 );
    }
    /************************ case of Fixed-H ******************************************************/
    if ( pStruct->bMobileH == TAUT_YES ) {
        if ( 0 > (ret = FixMobileHRestoredStructure(ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA, pTCGroups,
                                              &t_group_info, &at_norm, &at_prep, pInChI,
                                              num_inp, bHasSomeFixedH, pnNumRunBNS, pnTotalDelta, forbidden_edge_mask,
                                              forbidden_stereo_edge_mask) ) ) {
            goto exit_function;
        }
    }
    /**********************************************************************************************/
    /* stereo */
    cmpInChI = CompareReversedINChI2( pStruct->pOneINChI[0], pInChI[0], pStruct->pOneINChI_Aux[0], NULL /*INChI_Aux *a2*/, &icr, &err );
    if ( cmpInChI & IDIF_PROBLEM ) {
        ret = RI_ERR_PROGR; /* severe restore problem */
        goto exit_function;
    }
    if ( err ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    cmpInChI2 = 0;
    memset ( &icr2, 0, sizeof(icr2) );
    if ( iRevrInChI || iOrigInChI ) {
        /* additional mobile-H compare in case of Fixed-H */
        cmpInChI2 = CompareReversedINChI2( pStruct->pOneINChI[iRevrInChI], pInChI[iOrigInChI], pStruct->pOneINChI_Aux[iRevrInChI], NULL /*INChI_Aux *a2*/, &icr2, &err );
        if ( cmpInChI & IDIF_PROBLEM ) {
            ret = RI_ERR_PROGR; /* severe restore problem */
            goto exit_function;
        }
        if ( err ) {
            ret = RI_ERR_ALLOC;
            goto exit_function;
        }
    }
    ret = FixRestoredStructureStereo( cmpInChI, &icr, cmpInChI2, &icr2,
                                          ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA, pTCGroups,
                                          &t_group_info, &at_norm, &at_prep, pInChI,
                                          num_inp, pnNumRunBNS, pnTotalDelta, forbidden_edge_mask,
                                          forbidden_stereo_edge_mask);

    if ( ret < 0 ) {
        goto exit_function;
    }
#if( FIX_ADD_PROTON_FOR_ADP == 1 )
    /************************ check and fix ADP by adding a proton (dummy) *************************/
    if ( cmpInChI && pTCGroups->num_tgroups && pBNS->tot_st_cap > pBNS->tot_st_flow ) {
        ret = FixAddProtonForADP( pBNS, pBD, pStruct, at, at2, at_prep, pVA, pTCGroups, &icr,
                                  pnNumRunBNS, pnTotalDelta, forbidden_edge_mask );
        if ( ret < 0 ) {
            goto exit_function;
        }
    }
#endif
     /* moved to MakeOneInChIOutOfStrFromINChI():
       pStruct->nNumRemovedProtons = (pStruct->iMobileH == TAUT_YES)? pStruct->One_ti.tni.nNumRemovedProtons : 0;
     */

    /* count endpoints */
    num_endpoints      = 0;
    num_norm_endpoints = 0;
    num_norm_t_groups  = 0;
    num_mobile         = 0;
    num_norm_mobile    = 0;
    at_norm            = pStruct->pOne_norm_data[0]->at;
    for ( i = 0; i < pTCGroups->num_tgroups; i ++ ) {
        num_endpoints += pTCGroups->pTCG[i].num_edges;
        num_mobile    += pTCGroups->pTCG[i].tg_num_H + pTCGroups->pTCG[i].tg_num_Minus;
    }

    if ( t_group_info ) {
        /* after canonicalization, t_group_info->t_group[i].num[0] = number of H   */
        /*                         t_group_info->t_group[i].num[1] = number of (-) */
        for ( i = 0; i < t_group_info->num_t_groups; i ++ ) {
            if ( t_group_info->t_group[i].num[0] ) {
                num_norm_t_groups ++;
                num_norm_endpoints += t_group_info->t_group[i].nNumEndpoints;
                num_norm_mobile    += t_group_info->t_group[i].num[0]+t_group_info->t_group[i].num[1];
            }
        }
    }
#if( bRELEASE_VERSION == 0 )
#ifndef INCHI_LIBRARY
    if ( num_norm_t_groups != pTCGroups->num_tgroups || num_norm_endpoints != num_endpoints ) {
        /* need aggressive (de)protonation */
        /* pStruct->bExtract |= EXTRACT_STRUCT_NUMBER; */
        fprintf( stdout, "NORMCOMP: %s comp=%d %c%c: InChI/NormRvrs NumTg=%d/%d NumEndp=%d/%d\n",
                 (*ip).pSdfValue, (*pTCGroups).iComponent,
                 pStruct->iInchiRec? 'R':'D', pStruct->iMobileH?'M':'F',
                 pTCGroups->num_tgroups, num_norm_t_groups,
                 num_endpoints, num_norm_endpoints ); 
        
    }
#endif
#endif

exit_function:
    
    for( i = 0; i < TAUT_NUM; i ++ ) {
        Free_INChI( &pStruct->pOneINChI[i] );
        Free_INChI_Aux( &pStruct->pOneINChI_Aux[i] );
        FreeInpAtomData( pStruct->pOne_norm_data[i] );
        if ( pStruct->pOne_norm_data[i] ) {
            inchi_free( pStruct->pOne_norm_data[i] );
            pStruct->pOne_norm_data[i] = NULL;
        }
    }
    free_t_group_info( &pStruct->One_ti );
    return ret;
}
/******************************************************************************************************/
/* Find A=X< where all bonds to X except A=X are marked as stereogenic; temporary allow stereobonds   */
/* change and make A=X bonds single                                                                   */
int CheckAndRefixStereobonds(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                              inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int forbidden_edge_stereo     = BNS_EDGE_FORBIDDEN_MASK;
    int inv_forbidden_edge_stereo = ~forbidden_edge_stereo;

    int i, k, ne, j1, j2, num_wrong, num_fixed;
    int ret2, retBNS, ret;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    EDGE_LIST FixedEdges, WrongEdges, CarbonChargeEdges;

    BNS_EDGE   *pEdge;
    Vertex      v1, v2;
    BNS_VERTEX *pv1, *pv2;

    ret = 0;

    /* to simplify, prepare new at[] from pBNS */
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        return ret;
    }

    num_wrong = 0;
    /* find wrong double bonds */
    for ( i = 0; i < num_at; i ++ ) {
        if ( at2[i].valence == 3 &&
             at2[i].chem_bonds_valence - at2[i].valence == 1 &&
             at2[i].sb_parity[0] && at2[i].sb_parity[1] && !at2[i].sb_parity[2] &&
             (at2[i].bond_type[j1=(int)at2[i].sb_ord[0]] & BOND_TYPE_MASK) == BOND_TYPE_SINGLE &&
             (at2[i].bond_type[j2=(int)at2[i].sb_ord[1]] & BOND_TYPE_MASK) == BOND_TYPE_SINGLE  &&
             j1 != j2 ) {

            num_wrong ++;
        }
    }
    if ( !num_wrong ) {
        return 0;
    }
    num_fixed = 0;
    for ( i = 0; i < pBNS->num_bonds; i ++ ) {
        pEdge = pBNS->edge + i;
        if ( pEdge->forbidden & forbidden_edge_stereo ) {
            num_fixed ++;
        }
    }

    /* there may be no fixed stereo bonds at all, see #87607 */
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &FixedEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &WrongEdges, EDGE_LIST_CLEAR );

    /* do not goto exit_function before reaching this point: EdgeLists have not been initiated */

    if ( 0 > (ret = ForbidCarbonChargeEdges( pBNS, pTCGroups, &CarbonChargeEdges, forbidden_edge_mask ))) {
        goto exit_function;
    }
    if ( (ret = AllocEdgeList( &FixedEdges, num_fixed )) ||
         (ret = AllocEdgeList( &WrongEdges, num_wrong )) ) {
        goto exit_function;
    }
    /* collect wrong double bonds and set flow=0 */
    for ( i = 0; i < num_at && WrongEdges.num_edges < num_wrong; i ++ ) {
        if ( at2[i].valence == 3 &&
             at2[i].chem_bonds_valence - at2[i].valence == 1 &&
             at2[i].sb_parity[0] && at2[i].sb_parity[1] && !at2[i].sb_parity[2] &&
             (at2[i].bond_type[j1=(int)at2[i].sb_ord[0]] & BOND_TYPE_MASK) == BOND_TYPE_SINGLE &&
             (at2[i].bond_type[j2=(int)at2[i].sb_ord[1]] & BOND_TYPE_MASK) == BOND_TYPE_SINGLE &&
             j1 != j2 ) {
            switch ( j1 + j2 ) {
            case 1: /* 0, 1 */
                k = 2;
                break;
            case 2: /* 0, 2 */
                k = 1;
                break;
            case 3: /* 1, 2 */
                k = 0;
                break;
            default:
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            ne = pBNS->vert[i].iedge[k];
            pEdge = pBNS->edge + ne;
            v1 = pEdge->neighbor1;
            v2 = pEdge->neighbor12 ^ v1;
            pv1 = pBNS->vert + v1;
            pv2 = pBNS->vert + v2;

            if ( !pEdge->flow ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            pEdge->flow --;
            pEdge->forbidden |= forbidden_edge_mask;
            pv1->st_edge.flow --;
            pv2->st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            if ( ret = AddToEdgeList( &WrongEdges, ne, 0 )) {
                goto exit_function;
            }
        }
    }
    /* remove forbidden mark from stereo bonds (unfix stereo bonds) */
    for ( i = 0; i < pBNS->num_bonds && FixedEdges.num_edges < num_fixed; i ++ ) {
        pEdge = pBNS->edge + i;
        if ( pEdge->forbidden & forbidden_edge_stereo ) {
            pEdge->forbidden &= inv_forbidden_edge_stereo;
            FixedEdges.pnEdges[FixedEdges.num_edges ++] = i;
        }
    }
    /* Run BNS to move charges and rearrange bond orders */
    retBNS = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
    (*pnNumRunBNS) ++;
    if ( retBNS < 0 ) {
        goto exit_function;
    } else
    if ( retBNS > 0 ) {
        *pnTotalDelta += retBNS;
    }
    /* remove forbidden_edge_mask and set forbidden_edge_stereo */
    RemoveForbiddenEdgeMask( pBNS, &WrongEdges, forbidden_edge_mask );
    /* allow carbon charges to change */
    RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, forbidden_edge_mask );
    /* fix previously unfixed stereo bonds */
    SetForbiddenEdgeMask( pBNS, &FixedEdges, forbidden_edge_stereo );
    /* Run BNS again in case not all edge flows are maximal */
    ret2 = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
    (*pnNumRunBNS) ++;
    if ( ret2 < 0 ) {
        goto exit_function;
    } else
    if ( ret2 > 0 ) {
        *pnTotalDelta += retBNS;
    }
    ret = retBNS;

exit_function:

    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &FixedEdges, EDGE_LIST_FREE );
    AllocEdgeList( &WrongEdges, EDGE_LIST_FREE );

    return ret;
}
/******************************************************************************************************/
/* Find and eliminate false Mobile-H groups: Cl(=O)3(-O(-)) => Cl(-)(=O)4                             */
int MoveChargeToRemoveCenerpoints(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                              inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int i, j, neigh, num_endpoints, num_success;
    int num_donors, num_acceptors, bond_type, num_donors_O, num_acceptors_O, is_centerpoint_N, num_known_endpoints, num_wrong_neigh;
    int ret2, ret_forbid_edges, ret, delta;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int forbidden_edge_test = BNS_EDGE_FORBIDDEN_TEST;
    int bPossiblyIgnore = pStruct->charge >= 0 && (!pTCGroups->num_tgroups || pStruct->iMobileH == TAUT_NON && pStruct->ti.num_t_groups);
    S_CHAR MobileChargeNeigh[MAXVAL], DoubleBondAcceptors[MAXVAL], DoubleBondNotONeigh[MAXVAL];
    int    numMobileChargeNeigh, numDoubleBondAcceptors, numDoubleBondNotONeigh, numOtherDoubleBondOAcceptors=0;
    EDGE_LIST ChargeListAllExcept_DB_O;


    BNS_EDGE   *pEdgeMinus, *pe;
    Vertex      v1m, v2m;
    BNS_VERTEX *pv1m, *pv2m;
    ret = 0;
    num_success = 0;

    /* count O(+)H, N(+)H */

    /*
    if ( pStruct->charge >= 0 && (!pTCGroups->num_tgroups || pStruct->iMobileH == TAUT_NON && pStruct->ti.num_t_groups) ) {
        goto exit_function;
    }
    */
    if ( ret = AllocEdgeList( &ChargeListAllExcept_DB_O, EDGE_LIST_CLEAR ) ) {
        goto exit_function;
    }


    /* to simplify, prepare new at[] from pBNS */
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
#if( FIND_RING_SYSTEMS == 1 )
    ret2 = MarkRingSystemsInp( at2, num_at, 0 );
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
#endif
    /* mark bonds that cannot be tautomeric; do not forget to remove the marks later */
    ret_forbid_edges = SetForbiddenEdges( pBNS, at2, num_at, forbidden_edge_test );
    if ( ret_forbid_edges < 0 ) {
        ret = ret_forbid_edges;
        goto exit_function;
    }

    for ( i = 0; i < num_at; i ++ ) {
        if ( pVA[i].cNumValenceElectrons != 4 && /* not C, Si, Ge */
             !(pVA[i].nTautGroupEdge || pStruct->iMobileH == TAUT_NON && pStruct->endpoint && pStruct->endpoint[i] ) &&
             !at2[i].num_H && !at2[i].charge && at2[i].valence >= 2 &&
             at2[i].valence < at2[i].chem_bonds_valence &&
             is_centerpoint_elem( at2[i].el_number ) ) {
            
            is_centerpoint_N = (pVA[i].cNumValenceElectrons == 5 && (pVA[i].cPeriodicRowNumber == 1 || pVA[i].cMetal));
            /* look at the neighbors */
            numMobileChargeNeigh = numDoubleBondAcceptors = numDoubleBondNotONeigh = num_donors = num_acceptors = 0;
            num_donors_O = num_acceptors_O = 0;
            num_known_endpoints = num_wrong_neigh = 0;
            for ( j = 0, num_endpoints = 0; j < at2[i].valence; j ++ ) {
                neigh = at2[i].neighbor[j];
                if ( (at2[neigh].endpoint || pStruct->iMobileH == TAUT_NON && pStruct->endpoint && pStruct->endpoint[neigh]) || at2[neigh].charge > 0 ) {
                    num_known_endpoints ++;
                    continue;
                }
                if ( pBNS->edge[pBNS->vert[i].iedge[j]].forbidden & forbidden_edge_test ) {
                    continue;
                }
                bond_type = at2[i].bond_type[j] & BOND_TYPE_MASK;
                if ( bond_type > BOND_TYPE_DOUBLE ) {
                    num_wrong_neigh ++;
                    continue;
                }
                if ( at2[neigh].num_H && bond_type == BOND_TYPE_SINGLE ) {
                    break;  /* not this case */
                }
                if ( at2[neigh].chem_bonds_valence - at2[neigh].charge
                     != get_endpoint_valence( at2[neigh].el_number ) ) {
                    if ( bond_type == BOND_TYPE_DOUBLE && pVA[neigh].cNumValenceElectrons != 6 ) {
                        DoubleBondNotONeigh[numDoubleBondNotONeigh ++] = j;
                    }
                    continue;
                }
                if ( at2[neigh].charge == -1 && bond_type == BOND_TYPE_SINGLE &&
                     (pVA[neigh].nCMinusGroupEdge < 1 || pBNS->edge[pVA[neigh].nCMinusGroupEdge-1].flow != 1)  ) {
                    break;
                }
                switch( bond_type ) {
                case BOND_TYPE_SINGLE:
                    if ( at2[neigh].charge != -1 || pVA[neigh].nCMinusGroupEdge <= 0 ) {
                        num_wrong_neigh ++;
                        continue;
                    }
                    num_donors ++;
                    num_donors_O += (pVA[neigh].cNumValenceElectrons == 6 && pVA[neigh].cPeriodicRowNumber <= 4);
                    MobileChargeNeigh[numMobileChargeNeigh ++] = j;
                    break;
                case BOND_TYPE_DOUBLE:
                    if ( at2[neigh].charge ) {
                        num_wrong_neigh ++;
                        continue;
                    }
                    DoubleBondAcceptors[numDoubleBondAcceptors ++]  = j;
                    num_acceptors ++;
                    num_acceptors_O += (pVA[neigh].cNumValenceElectrons == 6 && pVA[neigh].cPeriodicRowNumber <= 4);
                }
            }
            if ( j != at2[i].valence || !num_donors || !num_acceptors ) {
                continue;
            }
            /* special case NOn(-) */
            if ( is_centerpoint_N && (num_donors == num_donors_O) && (num_acceptors == num_acceptors_O) ) {
                continue;
            }
            if ( pStruct->iMobileH == TAUT_NON && num_donors == numDoubleBondNotONeigh ) {
                /* fix all charges except on =O */
                Vertex     vPathStart, vPathEnd;
                int        nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
                int k, e, num_MovedCharges = 0;

                if ( !ChargeListAllExcept_DB_O.num_edges ) {
                    numOtherDoubleBondOAcceptors = 0;
                    for ( k = 0; k < num_at; k ++ ) {
                        if ( 1 == at2[k].valence && pBNS->edge[pBNS->vert[k].iedge[0]].flow &&
                             !pBNS->edge[pBNS->vert[k].iedge[0]].forbidden &&
                             !((e=pVA[k].nCMinusGroupEdge-1) >= 0 && pBNS->edge[e].flow) &&
                             !((e=pVA[k].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[e].flow) &&
                             /* 0 == at2[k].charge && */
                             pVA[k].cNumValenceElectrons == 6 && !pVA[k].cMetal &&
                             pStruct->endpoint && pStruct->endpoint[k] ||
                             pStruct->fixed_H && pStruct->fixed_H[k] ) {
                            numOtherDoubleBondOAcceptors ++;  /* do not fix this minus edge */
                        } else
                        if ( (e=pVA[k].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[e].flow &&
                             !pBNS->edge[e].forbidden &&
                             ( ret = AddToEdgeList( &ChargeListAllExcept_DB_O, e, 64 )) ) {
                            goto exit_function;
                        }
                        if ( (e=pVA[k].nCPlusGroupEdge-1) >= 0 &&
                             !pBNS->edge[e].forbidden &&
                             ( ret = AddToEdgeList( &ChargeListAllExcept_DB_O, e, 64 )) ) {
                            goto exit_function;
                        }
                    }
                }
                /* fix double bonds to non-O neighbors connected by double bonds;
                   we will try to make these bons single */
                for ( k = 0; k < numDoubleBondNotONeigh; k ++ ) {
                    e = pBNS->vert[i].iedge[(int)DoubleBondNotONeigh[k]];
                    if ( !pBNS->edge[e].forbidden &&
                         (ret = AddToEdgeList( &ChargeListAllExcept_DB_O, e, 64 ))) {
                        goto exit_function;
                    }
                }
                /* attempt to make DoubleBondNotONeigh[] single */
                SetForbiddenEdgeMask( pBNS, &ChargeListAllExcept_DB_O, forbidden_edge_mask);
                for ( k = 0; k < numDoubleBondNotONeigh && num_MovedCharges < numMobileChargeNeigh; k ++ ) {
                    pe = pBNS->edge + pBNS->vert[i].iedge[(int)DoubleBondNotONeigh[k]];
                    delta = 1;
                    if ( pe->flow != delta )
                        continue;
                    pv1m = pBNS->vert + (v1m = pe->neighbor1);
                    pv2m = pBNS->vert + (v2m = pe->neighbor12 ^ v1m);
                    pv1m->st_edge.flow -= delta;
                    pv2m->st_edge.flow -= delta;
                    pe->flow           -= delta;
                    pBNS->tot_st_flow  -= 2*delta;
                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                    if ( ret < 0 ) {
                        goto exit_function;
                    }
                    if ( ret == 1 && (vPathEnd == v1m && vPathStart == v2m ||
                                      vPathEnd == v2m && vPathStart == v1m) &&
                                      nDeltaCharge == 0  /* (-) moving from one to another atom*/ ) {
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        (*pnNumRunBNS) ++;
                        if ( ret < 0 ) {
                            goto exit_function;
                        } else
                        if ( ret == 1 ) {
                            *pnTotalDelta += ret;
                            num_MovedCharges ++;
                        } else {
                            ret = RI_ERR_PROGR;
                            goto exit_function;
                        }
                    } else {
                        ret = 0;
                        pv1m->st_edge.flow += delta;
                        pv2m->st_edge.flow += delta;
                        pe->flow           += delta;
                        pBNS->tot_st_flow  += 2*delta;
                    }
                }
                RemoveForbiddenEdgeMask( pBNS, &ChargeListAllExcept_DB_O, forbidden_edge_mask);
            } else
            if ( !bPossiblyIgnore || !num_known_endpoints && !num_wrong_neigh && (num_acceptors_O + num_donors_O) >=3  ) {
                /* remove negative charges from the neighbors */
                pBNS->vert[i].st_edge.cap += num_donors; /* enough to make all bonds to donors double */
                pBNS->tot_st_cap          += num_donors;
                pVA[i].cInitCharge        -= num_donors; /* work no matter what are known charge/valence */
                for ( j = 0; j < numMobileChargeNeigh; j ++ ) {
                    neigh = at2[i].neighbor[ (int)MobileChargeNeigh[j] ];
                    pEdgeMinus = pBNS->edge + (pVA[neigh].nCMinusGroupEdge-1);
                    v1m   = pEdgeMinus->neighbor1;
                    v2m   = pEdgeMinus->neighbor12 ^ v1m;
                    pv1m  = pBNS->vert + v1m;
                    pv2m  = pBNS->vert + v2m;
                    delta = pEdgeMinus->flow;
                    pv1m->st_edge.flow -= delta;
                    pv2m->st_edge.flow -= delta;
                    if ( IS_BNS_VT_C_GR( pv1m->type ) ) {
                        /* irreversible change to ChargeStruct */
                        pv1m->st_edge.cap -= delta;
                    } else
                    if ( IS_BNS_VT_C_GR( pv2m->type ) ) {
                        /* irreversible change to ChargeStruct */
                        pv2m->st_edge.cap -= delta;
                    } else {
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                    pBNS->tot_st_cap  -=   delta;
                    pBNS->tot_st_flow -= 2*delta;
                    pEdgeMinus->flow  -= delta;
                }
                ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                (*pnNumRunBNS) ++;
                if ( ret < 0 ) {
                    goto exit_function;
                } else
                if ( ret == num_donors ) {
                    *pnTotalDelta += ret;
                    num_success ++;
                    /*pStruct->bExtract |= EXTRACT_STRUCT_NUMBER;*/
                } else {
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                }
            }
        }
    }
    if ( ret_forbid_edges ) {
        /* remove the marks */
        RemoveForbiddenBondFlowBits( pBNS, forbidden_edge_test );
    }
    ret = num_success;
exit_function:
    AllocEdgeList( &ChargeListAllExcept_DB_O, EDGE_LIST_FREE );
    return ret;
}
/******************************************************************************************************/
/* Find and eliminate cases when Mobile H endpoint has radical on it (typical for wrong P(VI)(=O)3OH  */
int MakeSingleBondsMetal2ChargedHeteroat(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                              inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int i;
    
    int ret2, ret, pass;
    int num_at = pStruct->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;

    int         j, k;
    int        cur_num_edges;
    BNS_EDGE  *e;
    Vertex     v1, v2;

    EdgeIndex *pFixedEdges;
    int        nNumEdgesToFix;
    
    ret = 0;

    /* to simplify, prepare new at[] from pBNS */
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }

    pFixedEdges = NULL;

    nNumEdgesToFix = 0; /* cpunt nNumEdgesToFix only when pass==0 */
    cur_num_edges  = 0; /* count cur_num_edges  only when pass==1; at the end they must be equal */
    for ( pass = 0; pass < 2; pass ++ ) {
        if ( pass ) {
            /* 2nd pass: allocate edge storage */
            if ( !nNumEdgesToFix ) {
                break; /* nothing to do */
            }
            pFixedEdges = (EdgeIndex *)inchi_malloc(nNumEdgesToFix * sizeof( pFixedEdges[0] ) );
            if ( !pFixedEdges ) {
                ret = RI_ERR_ALLOC;
                goto exit_function;
            }
        }
        for ( i = 0; i < num_at; i ++ ) {
            int neigh;
            if ( pVA[i].cMetal ) {
                for ( j = 0; j < at2[i].valence; j ++ ) {
                    neigh = at2[i].neighbor[j];
                    if ( pVA[neigh].cNumValenceElectrons == 4 &&
                         pVA[neigh].cPeriodicRowNumber   == 1  ) {
                        continue; /* ignore carbon */
                    }
                    if ( at2[i].bond_type[j] > BOND_TYPE_SINGLE && at2[neigh].charge &&
                         !pVA[neigh].cMetal && pVA[neigh].cnListIndex > 0 ) {
                        int cnBits = at2[neigh].charge > 0? MAKE_CN_BITS(cn_bits_N, cn_bits_P, 0, 0) :
                                                            MAKE_CN_BITS(cn_bits_N, cn_bits_M, 0, 0);
                        int atBits = cnList[pVA[neigh].cnListIndex-1].bits;
                        for ( k = 0; k < MAX_NUM_CN_BITS-1; k ++, atBits >>= cn_bits_shift ) { /* ??? */
                            if ( (atBits & cnBits) == cnBits ) {
                                break;
                            }
                        }
                        if ( k == MAX_NUM_CN_BITS-1 ) {
                            continue;
                        }
                        if ( pass == 0 ) {
                            nNumEdgesToFix ++;
                        } else {
                            pFixedEdges[ cur_num_edges ++ ] = pBNS->vert[i].iedge[j];
                        }
                    }
                }
            }
        }
    }
    
    /* restore the initial structures */
    memcpy( at2, at, (num_at + num_deleted_H)*sizeof(at2[0]));

    if ( nNumEdgesToFix && pFixedEdges ) {
        if ( nNumEdgesToFix != cur_num_edges ) {
            ret = RI_ERR_PROGR;
            goto exit_function;
        }
        /* change edge flow, fix the edges, and run BNS */
        for ( i = 0; i < nNumEdgesToFix; i ++ ) {
            e = pBNS->edge + pFixedEdges[i];
            v1 = e->neighbor1;
            v2 = e->neighbor12 ^ v1;
            e->flow --;
            e->forbidden |= forbidden_edge_mask;
            pBNS->vert[v1].st_edge.flow --;
            pBNS->vert[v2].st_edge.flow --;
            pBNS->tot_st_flow -= 2;
            (*pnTotalDelta) -= 2;
        }
        /* Run BNS allowing to change any charges */
        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
        (*pnNumRunBNS) ++;
        if ( ret < 0 ) {
            goto exit_function;
        } else {
            (*pnTotalDelta) += ret;
        }
        /* unfix the edges */
        for ( i = 0; i < nNumEdgesToFix; i ++ ) {
            e = pBNS->edge + pFixedEdges[i];
            e->forbidden &= inv_forbidden_edge_mask;
        }
        if ( ret < 2 * nNumEdgesToFix ) {
            /* not all fixes succeeded */
            ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
            (*pnNumRunBNS) ++;
            if ( ret < 0 ) {
                goto exit_function;
            } else {
                (*pnTotalDelta) += ret;
            }
        }
    }
    if ( pFixedEdges ) {
        inchi_free( pFixedEdges );
        pFixedEdges = NULL;
    }


exit_function:
    return ret;
}
/**************************************************************************/
/* In Reconnected structure change 'salt bonds' to 'coordination bonds    */
/* for example, M-O-C=  ->  M(+)-O(-)-C=                                  */
/* Defect: instead of NH2-C=O(+)-M it will restore NH2(+)=C-O(-)-M(+)     */
/* However, in this release metal-organic compounds do not get much care  */
int SaltBondsToCoordBonds(BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                          inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                          int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)
{
    int i;
    
    int ret2, ret, cur_success;
    int num_at = pStruct->num_atoms;
    int num_edges = pBNS->num_bonds + 2 * pBNS->num_atoms;
    int num_deleted_H = pStruct->num_deleted_H;
    int len_at = num_at + num_deleted_H;
    int inv_forbidden_edge_mask = ~forbidden_edge_mask;
    EDGE_LIST AllChargeEdges;

    int         j, k, n;
    BNS_EDGE  *pe, *pePlusMetal, *peMinusO;
    BNS_VERTEX *pv1, *pv2, *pvO, *pvM;
    Vertex     v1, v2, vPlusMinus;

    EdgeIndex  ie, iePlusMetal, ieMinusO;
    
    Vertex     vPathStart, vPathEnd;
    int        delta, nPathLen, nDeltaH, nDeltaCharge, nNumVisitedAtoms;
    
    ret = 0;
    cur_success = 0;
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_CLEAR );

    if ( pStruct->iInchiRec == INCHI_BAS || !pStruct->pSrm->bMetalAddFlower || pStruct->pSrm->nMetalMinBondOrder ) {
        goto exit_function;
    }

    /* to simplify, prepare new at[] from pBNS */
    memcpy( at2, at, len_at*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
    if ( ret2 < 0 ) {
        ret = ret2;
        goto exit_function;
    }
    for ( i = 0; i < num_at; i ++ ) {
        if ( bIsMetalSalt( at2, i ) ) {
            if ( !AllChargeEdges.num_edges ) {
                /*--------- one-time action: fix all bonds, charges, taut. group edges ------------*/
                for ( j = 0; j < num_at; j ++ ) {
                    /* all bonds */
                    for ( k = 0; k < at2[j].valence; k ++ ) {
                        n = at2[j].neighbor[k];
                        if ( n < j && !pBNS->edge[ie = pBNS->vert[j].iedge[k]].forbidden &&
                             ( ret = AddToEdgeList( &AllChargeEdges, ie, num_edges ) ) ) {
                            goto exit_function;
                        }
                    }
                    /* charge edges */
                    if ( (ie=pVA[j].nCMinusGroupEdge-1) >= 0 && !pBNS->edge[ie].forbidden &&
                         (ret = AddToEdgeList( &AllChargeEdges, ie, num_edges ) ) ) {
                        goto exit_function;
                    }
                    if ( (ie=pVA[j].nCPlusGroupEdge-1) >= 0 && !pBNS->edge[ie].forbidden &&
                         ( ret = AddToEdgeList( &AllChargeEdges, ie, num_edges ) ) ) {
                            goto exit_function;
                    }
                }
                /* taut group edges */
                for ( j = 0; j < pTCGroups->num_tgroups; j ++ ) {
                    pv1  = pBNS->vert + (v1=pTCGroups->pTCG[j].nVertexNumber); /* t-group vertex */
                    for ( k = 0; k < pv1->num_adj_edges; k ++ ) {
                        /* ie, pe - tautomeric atom edge; pv2 - endpoint vertex */
                        /* Note: pe, pv2, v1 are not used here; they are to show how to traverse t-group */
                        pv2 = pBNS->vert + (pe = pBNS->edge + (ie=pv1->iedge[k]))->neighbor1;
                        if ( ret = AddToEdgeList( &AllChargeEdges, ie, num_edges ) ) {
                            goto exit_function;
                        }
                    }
                }
                /*---------------------------------------------------------------*/
            }
            /* replace all single bonds to neutral neighbors with zero-order bonds
               allow neighbor charge change to (-1) and metal atom charge increment +1 */
            for ( k = 0; k < at2[i].valence; k ++ ) {
                n  = at2[i].neighbor[k];
                pe = pBNS->edge + pBNS->vert[i].iedge[k];
                if ( at2[n].charge || at2[i].bond_type[k] != BOND_TYPE_SINGLE ) {
                    continue;
                }
                iePlusMetal = pVA[i].nCPlusGroupEdge-1;
                ieMinusO    = pVA[n].nCMinusGroupEdge-1;

                if ( pe->flow != 1 || pe->forbidden || iePlusMetal < 0 ) {
                    continue;
                }
                pePlusMetal = pBNS->edge + iePlusMetal;
                if ( pePlusMetal->flow <= 0 ) {
                    continue; /* to add (+) to metal this flow must be decremented */
                }
                if ( ieMinusO >= 0 ) {
                    /* usually does not happen */
                    peMinusO    = pBNS->edge + ieMinusO;

                    if ( peMinusO->flow || pePlusMetal->forbidden || peMinusO->forbidden ) {
                        continue;
                    }

                    /* decrement bond order to 0 */
                    delta = 1;
                    pv1 = pBNS->vert + (v1 = pe->neighbor1);
                    pv2 = pBNS->vert + (v2 = pe->neighbor12 ^ v1);

                    pe->flow -= delta;
                    pv1->st_edge.flow -= delta;
                    pv2->st_edge.flow -= delta;
                    pBNS->tot_st_flow -= 2*delta;

                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                    pePlusMetal->forbidden &= inv_forbidden_edge_mask;
                    peMinusO->forbidden &= inv_forbidden_edge_mask;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );

                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) /*&& nDeltaCharge > 0*/ ) {
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
                    RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                } else
                if ( NO_VERTEX != (vPlusMinus = GetPlusMinusVertex( pBNS, pTCGroups, 1, 1 ) ) ) {
                    /* manually add (-) charge to O and (+) charge to metal */
                    /* decrement bond order to 0 */
                    /*---------------------------------------------------------------------------*/
                    /*                                                                           */
                    /*      (+/-)*               (+/-)           Result:                         */
                    /*        |                    ||                                            */
                    /*        |                    ||            - Added (+) to M                */
                    /*       (+)super             (+)super       - Incremented bond M-O          */
                    /*        ||                   |                                             */
                    /*        ||          =>       |             To make this attachment H,      */
                    /*       (Y)                  (Y)            increment                       */
                    /*        |                    ||            pTCGroups->pTCG[itg].tg_num_H   */
                    /*        |                    ||                                            */
                    /*       (+)metal             (+)hetero      Technical details:              */
                    /*         \\                   \            increase capacities of          */
                    /*           M                    M(+)       edges to (+/-) otherwise        */
                    /*           |                    ||         flow may not be able to         */
                    /*          -O*                -O-O          increase                        */
                    /*                                                                           */
                    /*   After that change M=O bond order from 2 to 0                            */
                    /*---------------------------------------------------------------------------*/
                    int i1, j1, k1;
                    delta = 1;
                    pvO = pBNS->vert + n;
                    pvM = pBNS->vert + i;
                    /* Increment st_edge.cap on (+/-) vertex */
                    pBNS->vert[vPlusMinus].st_edge.cap += delta;
                    /* Increment st_edge.cap on O */
                    pvO->st_edge.cap                   += delta;
                    /* increment cap on M-O edge */
                    pe->cap += delta;
                    /* total cap count */
                    pBNS->tot_st_cap                   += 2*delta;

                    v1 = vPlusMinus;
                    v2 = n; /* atom O */

                    /* increase capacities of edges to Y  */
                    for ( i1 = 0; i1 < pBNS->vert[vPlusMinus].num_adj_edges; i1 ++ ) {
                        j1 = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i1]].neighbor12 ^ vPlusMinus;
                        for ( k1 = 0; k1 < pBNS->vert[j1].num_adj_edges; k1 ++ ) {
                            pBNS->edge[pBNS->vert[j1].iedge[k1]].cap += delta;
                        }
                    }
                    SetForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                    pePlusMetal->forbidden &= inv_forbidden_edge_mask;
                    pe->forbidden          &= inv_forbidden_edge_mask;

                    ret = RunBnsTestOnce( pBNS, pBD, pVA, &vPathStart, &vPathEnd, &nPathLen,
                                          &nDeltaH, &nDeltaCharge, &nNumVisitedAtoms );
                    cur_success = 0;
                    if ( ret == 1 && (vPathEnd == v1 && vPathStart == v2 ||
                                      vPathEnd == v2 && vPathStart == v1) /*&& nDeltaCharge == 1*/ ) {
                        /* Added (+)charge to -N< => nDeltaCharge == 1 */
                        /* Flow change on pe (-)charge edge (atom B-O(-)) is not known to RunBnsTestOnce()) */
                        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
                        if ( ret > 0 ) {
                            (*pnNumRunBNS) ++;
                            cur_success ++; /* 01 */
                        }
                    }
                    if ( cur_success ) {
                        /* set bond M=O order = 0 */
                        if ( pe->flow != 2*delta ) {
                            ret = RI_ERR_PROGR;
                            goto exit_function;
                        }
                        /* reduce pe bond order by 2*delta */
                        pe->flow            -= 2*delta;
                        pvO->st_edge.cap    -= 2*delta;
                        pvO->st_edge.flow   -= 2*delta;
                        pvM->st_edge.flow   -= 2*delta;
                        pvM->st_edge.cap    -= 2*delta;
                        pBNS->tot_st_cap    -= 3*delta;
                        pBNS->tot_st_flow   -= 4*delta;
                        /* fix M-O bond order to zero */
                        pe->cap             -= 2*delta;
                        /* add fixed (-) charge to O */
                        pVA[n].cInitCharge -= delta;
                    } else {
                        /* failed */
                        pBNS->vert[vPlusMinus].st_edge.cap -= delta;
                        pvO->st_edge.cap                   -= delta;
                        /*pTCGroups->pTCG[itg].edges_cap     -= delta;*/ /* ???bug??? - commented out 2006-03-22 */
                        pBNS->tot_st_cap                   -= 2*delta;
                        /* decrease capacities of edges to Y  */
                        for ( i1 = 0; i1 < pBNS->vert[vPlusMinus].num_adj_edges; i1 ++ ) {
                            j1 = pBNS->edge[pBNS->vert[vPlusMinus].iedge[i1]].neighbor12 ^ vPlusMinus;
                            for ( k1 = 0; k1 < pBNS->vert[j1].num_adj_edges; k1 ++ ) {
                                pBNS->edge[pBNS->vert[j1].iedge[k1]].cap -= delta;
                            }
                        }
                    }
                    RemoveForbiddenEdgeMask( pBNS, &AllChargeEdges, forbidden_edge_mask  );
                }
            }
        }
    }

exit_function:
    AllocEdgeList( &AllChargeEdges, EDGE_LIST_FREE );
    return ret;
}
#if ( KEEP_METAL_EDGE_FLOW == 1 )
/******************************************************************************************************/
int ForbidMetalCarbonEdges( BN_STRUCT *pBNS, inp_ATOM *at, int num_at, VAL_AT *pVA,
                           ALL_TC_GROUPS *pTCGroups, EDGE_LIST *pMetalCarbonEdges, int forbidden_edge_mask  )
{

    int i, j, neigh, nNumEdgeMetalCarbon = 0, pass = 0, ret = 0;
    BNS_VERTEX *pVert, *pNeigh;
    BNS_EDGE   *pEdge;

    /* count carbon-metal edges */
    
    if ( pTCGroups->num_metal_atoms ) {
fill_ForbiddenEdgesMetalCarbon:
        for ( i = 0; i < num_at; i ++ ) {
            if ( pVA[i].cMetal && pVA[i].cNumBondsToMetal ) {
                pVert = pBNS->vert + i;
                for ( j = 0; j < pVert->num_adj_edges; j ++ ) {
                    pEdge = pBNS->edge + pVert->iedge[j];
                    neigh = pEdge->neighbor12 ^ i;
                    pNeigh = pBNS->vert + neigh;
                    if ( !IS_BNS_VT_ATOM(pNeigh->type) )
                        continue;
                    if ( at[neigh].endpoint )
                        continue;
                    if ( pVA[neigh].cNumValenceElectrons == 4 && pVA[neigh].cPeriodicRowNumber == 1 &&
                         pNeigh->st_edge.cap >= at[neigh].valence+1 ) {
                        if ( pass ) {
                            if ( ret = AddToEdgeList( pMetalCarbonEdges, pVert->iedge[j], 0 ) ) {
                                goto exit_function;
                            }
                            pEdge->forbidden |= forbidden_edge_mask;
                        } else {
                            nNumEdgeMetalCarbon ++;
                        }
                    }

                }
            }
        }
        if ( !pass && nNumEdgeMetalCarbon ) {
            if ( ret = AllocEdgeList( pMetalCarbonEdges, nNumEdgeMetalCarbon ) ) {
                goto exit_function;
            }
            pass ++;
            goto fill_ForbiddenEdgesMetalCarbon;
        }
    }
exit_function:
    return ret;
}
#endif
/******************************************************************************************************/
int RunBnsRestore1( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
                    StrFromINChI *pStruct, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, INChI *pInChI[],
                    long num_inp, int bHasSomeFixedH )
{
    int        nNumRunBNS = 0;
    
    EDGE_LIST CarbonChargeEdges, MetalCarbonEdges, Nplus2BondsEdges;

    int nTotalDelta = 0, ret = 0, tot_num_fixes;
    inp_ATOM *at          = pStruct->at;
    inp_ATOM *at2         = NULL; /* restored structure */
    inp_ATOM *at3         = NULL; /* structure for calculating one InChI */
    int     num_at        = pStruct->num_atoms;
    int     num_deleted_H = pStruct->num_deleted_H;
#ifdef _DEBUG
    int ret2;
#endif

#if ( KEEP_METAL_EDGE_FLOW == 1 )    
    BNS_VERTEX *pVert, *pNeigh;
    int         j, neigh;
#endif

    /* Edge lista initialization */
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &MetalCarbonEdges, EDGE_LIST_CLEAR );
    AllocEdgeList( &Nplus2BondsEdges, EDGE_LIST_CLEAR );

    if ( pStruct->iMobileH == TAUT_NON &&
         ( ret = FillOutExtraFixedHDataInChI( pStruct, pInChI ) ) ) {
        goto exit_function;
    }

    if ( !at2 && !(at2 = (inp_ATOM *)inchi_malloc((num_at + num_deleted_H)*sizeof(at2[0]))) ||
         !at3 && !(at3 = (inp_ATOM *)inchi_malloc((num_at + num_deleted_H)*sizeof(at3[0])))) {
        return RI_ERR_ALLOC;
    }

    if ( 0 > (ret = ForbidCarbonChargeEdges( pBNS, pTCGroups, &CarbonChargeEdges, BNS_EDGE_FORBIDDEN_TEMP ))) {
        goto exit_function;
    }

#if ( KEEP_METAL_EDGE_FLOW == 1 )    
    /* count edges of -C(IV)< carbons connected to metals */
    if ( 0 > (ret = ForbidMetalCarbonEdges( pBNS, at, num_at, pVA, pTCGroups, &MetalCarbonEdges, BNS_EDGE_FORBIDDEN_TEMP ))) {
        goto exit_function;
    }
#endif
    if ( 0 > (ret = ForbidNintrogenPlus2BondsInSmallRings( pBNS, at, num_at, pVA, 6,
                                           pTCGroups, &Nplus2BondsEdges, BNS_EDGE_FORBIDDEN_TEMP ) ) ) {
        goto exit_function;
    }

    /*********** Run BNS #1: no charge on carbons and =N= ***************/
    if ( Nplus2BondsEdges.num_edges ) {
        /* Run BNS leaving carbon charges unchanged */
        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
        nNumRunBNS ++;
        if ( ret < 0 ) {
            goto exit_function;
        } else {
            nTotalDelta += ret;
        }
        RemoveForbiddenEdgeMask( pBNS, &Nplus2BondsEdges, BNS_EDGE_FORBIDDEN_TEMP );
        AllocEdgeList( &Nplus2BondsEdges, EDGE_LIST_FREE );
    }
#ifdef _DEBUG
    /* debug only */
    memcpy( at2, at, (pStruct->num_atoms + pStruct->num_deleted_H)*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
#endif    
    /*************************** extend min ring size to 8 ****************************/
    if ( 0 > (ret = ForbidNintrogenPlus2BondsInSmallRings( pBNS, at, num_at, pVA, 8,
                                           pTCGroups, &Nplus2BondsEdges, BNS_EDGE_FORBIDDEN_TEMP ) ) ) {
        goto exit_function;
    }
    if ( Nplus2BondsEdges.num_edges ) {
        /* Run BNS leaving carbon charges unchanged */
        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
        nNumRunBNS ++;
        if ( ret < 0 ) {
            goto exit_function;
        } else {
            nTotalDelta += ret;
        }
        RemoveForbiddenEdgeMask( pBNS, &Nplus2BondsEdges, BNS_EDGE_FORBIDDEN_TEMP );
        AllocEdgeList( &Nplus2BondsEdges, EDGE_LIST_FREE );
    }
#ifdef _DEBUG
    /* debug only */
    memcpy( at2, at, (pStruct->num_atoms + pStruct->num_deleted_H)*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
#endif    
    /*******************************************************************/
    if ( CarbonChargeEdges.num_edges > 0 ) {
        /* Run BNS leaving carbon charges unchanged */
        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
        nNumRunBNS ++;
        if ( ret < 0 ) {
            goto exit_function;
        } else {
            nTotalDelta += ret;
        }
        RemoveForbiddenEdgeMask( pBNS, &CarbonChargeEdges, BNS_EDGE_FORBIDDEN_TEMP );
        AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    }
#ifdef _DEBUG
    /* debug only */
    memcpy( at2, at, (pStruct->num_atoms + pStruct->num_deleted_H)*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
#endif    
    /*******************************************************************/
    if ( MetalCarbonEdges.num_edges > 0 ) {
        /* Run BNS leaving carbon charges unchanged */
        ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
        nNumRunBNS ++;
        if ( ret < 0 ) {
            goto exit_function;
        } else {
            nTotalDelta += ret;
        }
        RemoveForbiddenEdgeMask( pBNS, &MetalCarbonEdges, BNS_EDGE_FORBIDDEN_TEMP );
        AllocEdgeList( &MetalCarbonEdges, EDGE_LIST_FREE );
    }
    /*******************************************************************/
    /* Run BNS allowing to change any charges */
    ret = RunBnsRestoreOnce( pBNS, pBD, pVA, pTCGroups );
    nNumRunBNS ++;
    if ( ret < 0 ) {
        goto exit_function;
    } else {
        nTotalDelta += ret;
    }
#ifdef _DEBUG
    /* debug only */
    memcpy( at2, at, (pStruct->num_atoms + pStruct->num_deleted_H)*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
#endif

#if ( BNS_RAD_SEARCH == 1 )
    /******************************************************************/
    /* move unfulfilled 'radicals' from ChargeStruct to atoms         */
    /* and set change charges of affected atoms to fit total charge   */
     ret = MoveRadToAtomsAddCharges( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups, BNS_EDGE_FORBIDDEN_TEMP );
    if ( ret < 0 ) {
        goto exit_function;
    }
#endif
    /**************************************************************/
    /**************************************************************/
    /*****           fix restore inconsistencies              *****/
    /**************************************************************/
    /**************************************************************/
#ifdef _DEBUG
    /* debug only */
    memcpy( at2, at, (pStruct->num_atoms + pStruct->num_deleted_H)*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
#endif    

    /* rearrange (+) and (-) edges flow so that there is no (+)flow=0 and (-)flow=1 */
    ret = RearrangePlusMinusEdgesFlow( pBNS, pBD, pVA, pTCGroups, BNS_EDGE_FORBIDDEN_TEMP );
    if ( ret < 0 ) {
        goto exit_function;
    }

    /*****************************************************************/
    /*       Increment zero order metal bonds to heteroatoms         */
    /*****************************************************************/
    ret = IncrementZeroOrderBondsToHeteroat( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                             &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    
#ifdef _DEBUG
    /* debug only */
    memcpy( at2, at, (pStruct->num_atoms + pStruct->num_deleted_H)*sizeof(at2[0]));
    pStruct->at = at2;
    ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
    pStruct->at = at;
#endif

#if (MOVE_CHARGES_FROM_HETEREO_TO_METAL == 1 )
    /*****************************************************************/
    /* move charges from heteroatoms to metal atoms                  */
    /*****************************************************************/
    ret = MoveChargeFromHeteroatomsToMetals( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                             &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
#endif
    /***********************************************************************
            NH2                NH2
               \                  \
                C==S(+)-   =>      C(+)-S-   where NH2 are not tautomeric
               /                  /
            NH2                NH2
    ************************************************************************/
    ret = MovePlusFromS2DiaminoCarbon( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                       &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /*****************************************************************/
    /*       Avoid charge separation on heteroatoms                  */
    /*****************************************************************/
    ret = EliminateChargeSeparationOnHeteroatoms( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                                  &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP, 0);
    if ( ret < 0 ) {
        goto exit_function;
    }
    if ( ret ) {
        /*charge separation remains; allow changes of stereobonds in a ring and try again */
        ret = EliminateChargeSeparationOnHeteroatoms( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                                      &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP,
                                                      BNS_EDGE_FORBIDDEN_MASK);
        if ( ret < 0 ) {
            goto exit_function;
        }
    }
    /*****************************************************************/
    /*         convert N#N(+)-N= into N(-)=N(+)=N-                   */
    /*****************************************************************/
    ret = RestoreNNNgroup( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                           &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /*****************************************************************/
    /*     convert Metal(q)-N(-)-O(-) Metal(q-2)-N=O (local change)  */
    /*****************************************************************/
    ret = FixMetal_Nminus_Ominus( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                           &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /*****************************************************************/
    /*         convert N(-)=C= into N#C-         -                   */
    /*****************************************************************/
    ret = RestoreCyanoGroup( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                           &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /*****************************************************************/
    /*         convert C(+)#N(+)- into C(-)#N(+)-                    */
    /*****************************************************************/
    ret = RestoreIsoCyanoGroup( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                           &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /*****************************************************************/
    /*         eliminate =N(V)= if possible                          */
    /*                    |                                          */
    /*****************************************************************/
    ret = EliminateNitrogen5Val3Bonds(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                      &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }

    /*****************************************************************/
    /*                    |      |                                   */
    /*         convert   -S- to =S= if possible                      */
    /*                    |      |                                   */
    /*****************************************************************/
    ret = Convert_SIV_to_SVI(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                             &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }

    /*****************************************************************/
    /*                  =N(+)=O     =N-O(-)                          */
    /*         convert           => if possible                      */
    /*                  Metal(q)    Metal(q+2)                       */
    /*****************************************************************/
    ret = PlusFromDB_N_DB_O_to_Metal(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                             &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }

    /*****************************************************************/
    /*  forbidden edges prevents required in InChI tautomerism       */
    /*  incorrectly restored mobile H mix separate tautomeric groups */
    /*  because an edge may not become forbidden                     */
    /* note: removes this 'forbidden_edge' bit from ALL edges        */
    /*****************************************************************/
    ret = MoveMobileHToAvoidFixedBonds( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                      &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);

    if ( ret < 0 ) {
        goto exit_function;
    }
    /**************************************************************************/
    /* 2. Mobile H endpoint has radical on it (typical for wrong P(VI)(=O)3OH */
    tot_num_fixes = 0;
    if ( pStruct->iMobileH==TAUT_NON ) {
        ret = RemoveRadFromMobileHEndpointFixH( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                          &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    } else {
        ret = RemoveRadFromMobileHEndpoint( pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                      &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    }
    if ( ret < 0 ) {
        goto exit_function;
    }
    tot_num_fixes += ret;
    /**************************************************************/
    /* make bonds between a charged heteroatom and a metal single */
    ret = MakeSingleBondsMetal2ChargedHeteroat(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                      &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /**************************************************************/
    /* move (+) charges to >N- and other centerpoints             */
    ret = MoveChargeToMakeCenerpoints(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                      &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }

    /**************************************************************************/
    /* Find and eliminate false Mobile-H groups: Cl(=O)3(-O(-)) => Cl(-)(=O)4 */
    ret = MoveChargeToRemoveCenerpoints(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                      &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /**************************************************************************/
    /* Find A=X< where all bonds to X except A=X are marked as stereogenic    */
    /* make bonds A=X single                                                  */
    ret = CheckAndRefixStereobonds(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                                 &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /**************************************************************************/
    /* In Reconnected structure change 'salt bonds' to 'coordination bonds    */
    /* for example, M-O-C=  ->  M(+)-O(-)-C=                                  */
    /* Defect: instead of NH2-C=O(+)-M it will restore NH2(+)=C-O(-)-M(+)     */
    /* However, in this release metal-organic compounds do not get much care  */
    ret = SaltBondsToCoordBonds(pBNS, pBD, pStruct, at, at2, pVA, pTCGroups,
                               &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /**************************************************************************/
    /* Normalize the structure and compare t-groups and stereobonds           */
    ret = NormalizeAndCompare(ip, sd, pBNS, pBD, pStruct, at, at2, at3, pVA, pTCGroups, pInChI, num_inp, bHasSomeFixedH,
                              &nNumRunBNS, &nTotalDelta, BNS_EDGE_FORBIDDEN_TEMP, BNS_EDGE_FORBIDDEN_MASK);
    if ( ret < 0 ) {
        goto exit_function;
    }
    /**************************************************************************/
    /* Create InChI out of the restored structure                             */


    /*ret = nTotalDelta;*/

exit_function:
    pStruct->at  = at;
    pStruct->at2 = at2;
    at2 = NULL;
    AllocEdgeList( &CarbonChargeEdges, EDGE_LIST_FREE );
    AllocEdgeList( &MetalCarbonEdges, EDGE_LIST_FREE );
    AllocEdgeList( &Nplus2BondsEdges, EDGE_LIST_FREE );
    if ( at2 ) {
        inchi_free( at2 );
    }
    if ( at3 ) {
        inchi_free( at3 );
    }

    return ret;
}

/******************************************************************************************************/
int RestoreAtomMakeBNS( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, StrFromINChI *pStruct, int iComponent,
                        int iAtNoOffset, INChI *pInChI[], const char *szCurHdr, long num_inp, int bHasSomeFixedH )
{
    int i, j, ret = 0, ret2;
    /*int nDelta, nTotalDelta;*/
    VAL_AT   *pVA = NULL;
    VAL_AT    va1;
    int    num_at = pStruct->num_atoms;
    inp_ATOM *at  = pStruct->at;
    ALL_TC_GROUPS   TCGroups;
    ALL_TC_GROUPS *pTCGroups = &TCGroups;
    int            nAddEdges2eachAtom = 2, nAddVertices = 0;

    BFS_Q bfsq;

    /* BNS creation */
    BN_STRUCT     *pBNS = NULL;
    BN_DATA       *pBD  = NULL;
    int            nNum_changed_bonds = 0;
    int            bTreatMoreAtomsAsMetals = 0, bSecondPassNewMetals=0;
    int            nMaxAddAtoms = 2, nMaxAddEdges = 2, max_altp = BN_MAX_ALTP;

    memset( pTCGroups, 0, sizeof(pTCGroups[0]) );
    for ( i = 0; i < NUM_TCGROUP_TYPES; i ++ ) {
        pTCGroups->nGroup[i] = TCG_None; /* unassigned */
    }
    pTCGroups->iComponent = iComponent;
    pTCGroups->iAtNoOffset = iAtNoOffset;

    if ( num_at == 1 ) {
        /* single atom -- no bonds to restore */
        inp_ATOM *at2 = (inp_ATOM *)inchi_malloc(sizeof(at2[0])*(pStruct->num_atoms+pStruct->num_deleted_H));
        inp_ATOM *at3 = (inp_ATOM *)inchi_malloc(sizeof(at3[0])*(pStruct->num_atoms+pStruct->num_deleted_H));
        pStruct->at2 = at2;
        at[0].charge = pInChI[0]->nTotalCharge;
        if ( at2 ) {
            memcpy( at2, at, sizeof(at2[0])*(pStruct->num_atoms+pStruct->num_deleted_H));
        }
        if ( !at2 || !at3 ) {
            if ( at3 ) inchi_free( at3 );
            return RI_ERR_ALLOC;
        }
        ret = MakeOneInChIOutOfStrFromINChI( ip, sd, pStruct, pStruct->at2, at3, pTCGroups );
        /* clean up */
        for( i = 0; i < TAUT_NUM; i ++ ) {
            Free_INChI( &pStruct->pOneINChI[i] );
            Free_INChI_Aux( &pStruct->pOneINChI_Aux[i] );
            FreeInpAtomData( pStruct->pOne_norm_data[i] );
            if ( pStruct->pOne_norm_data[i] ) {
                inchi_free( pStruct->pOne_norm_data[i] );
                pStruct->pOne_norm_data[i] = NULL;
            }
        }
        free_t_group_info( &pStruct->One_ti );
        inchi_free( at3 );

        return ret;
    }
    
    AllocBfsQueue( &bfsq, BFS_Q_CLEAR, 0 );
    if ( !(pVA = (VAL_AT *) inchi_calloc( num_at, sizeof( pVA[0] ) ) ) ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    pStruct->pVA = pVA;
    memset( &va1, 0, sizeof(va1) );
    pTCGroups->total_charge = pInChI[0]->nTotalCharge;
    if ( 0 > ( ret = AllocBfsQueue( &bfsq, num_at, 0 /* min ring size undefined */ ) ) ) {
        goto exit_function;
    }
    pStruct->pbfsq = &bfsq;
    
    if ( pStruct->iMobileH == TAUT_NON && pInChI[1] && pInChI[1]->nNumberOfAtoms > 1 &&
        ( ret = FillOutpStructEndpointFromInChI( pInChI[1], &pStruct->endpoint )) ) {
        goto exit_function;
    }
    
    /* mark metal atoms; find min ring sizes for atoms that have 2 bonds */
    for ( i = 0; i < num_at; i ++ ) {
        pVA[i].cNumValenceElectrons = get_sp_element_type( at[i].el_number, &j );
        pVA[i].cPeriodicRowNumber = j;
        pVA[i].cPeriodicNumber    = at[i].el_number;
        pVA[i].cNumValenceElectrons --; /* = -1 d- and f- metals, 0 for H, 1 for Na, 2 for Mg,.. = (ATYPE_Xx-1)  */

        if ( is_el_a_metal( at[i].el_number ) ) {
            if ( pStruct->pSrm->bStereoRemovesMetalFlag ) {
                /* treat metal as non-metal if it is stereogenic or has a stereobond */
                pVA[i].cMetal = !( at[i].p_parity || at[i].sb_parity[0] );
            } else {
                pVA[i].cMetal = 1;
            }
        }
        if ( at[i].valence == 2 && !at[i].num_H ) {
            pVA[i].cMinRingSize = is_bond_in_Nmax_memb_ring( at, i, 0, bfsq.q, bfsq.nAtomLevel,
                                                             bfsq.cSource, 99 /* max ring size */ );
        } else {
            pVA[i].cMinRingSize = 0;
        }
    }
    /* AllocBfsQueue( &bfsq, BFS_Q_FREE, 0 ); */

repeat_for_new_metals:
    /* set valences for the first time; find ChargeValence structures for each atom */
    for ( i = 0; i < num_at; i ++ ) {
        /* get additional fictitious atoms information */
        pVA[i].cInitFreeValences = 0;
        ret = GetAtomRestoreInfo( at, i, pVA, pStruct->pSrm, pStruct->bMobileH, pStruct->endpoint );
        if ( ret < 0 ) {
            goto exit_function;
        }
        if ( ret == TREAT_ATOM_AS_METAL && !bSecondPassNewMetals && !pVA[i].cMetal ) {
            if ( pStruct->pSrm->bStereoRemovesMetalFlag ) {
                /* treat metal as non-metal if it is stereogenic or has a stereobond */
                pVA[i].cMetal = !( at[i].p_parity || at[i].sb_parity[0] );
            } else {
                pVA[i].cMetal = 1;
            }
            if ( pVA[i].cMetal ) {
                bTreatMoreAtomsAsMetals ++;
            }
        }
        pTCGroups->charge_on_atoms += pVA[i].cInitCharge;
    }
    if ( bTreatMoreAtomsAsMetals && !bSecondPassNewMetals ) {
        for ( i = 0; i < num_at; i ++ ) {
            /* clear all members of pVA[i] except two */
            pTCGroups->charge_on_atoms -= pVA[i].cInitCharge;
            va1.cMetal               = pVA[i].cMetal;
            va1.cMinRingSize         = pVA[i].cMinRingSize;
            va1.cNumValenceElectrons = pVA[i].cNumValenceElectrons;
            va1.cPeriodicRowNumber   = pVA[i].cPeriodicRowNumber;
            va1.cPeriodicNumber      = pVA[i].cPeriodicNumber;
            pVA[i]           = va1;
        }
        bSecondPassNewMetals = 1;
        goto repeat_for_new_metals;
    }

    /* count atoms, bonds, additional edges and vertices in ChargeValence structures and t-groups */
    ret = nCountBnsSizes( at, num_at,  nAddEdges2eachAtom, nAddVertices, &pStruct->ti,
                          pVA, pStruct->pSrm, pTCGroups );
    if ( ret < 0 ) {
        goto exit_function;
    }
    
    /* find and count groups; add counts of all other vertices to be created */
    ret = nAddSuperCGroups( pTCGroups );
    if ( ret < 0 ) {
        goto exit_function;
    }

    /* create the BNS and fill it with all real atoms */
    pBNS = AllocateAndInitTCGBnStruct( pStruct, pVA, pTCGroups,
                                       nMaxAddAtoms, nMaxAddEdges, max_altp, &nNum_changed_bonds );
    if ( !pBNS ) {
        ret = BNS_OUT_OF_RAM;
        goto exit_function;
    }
    /* add t-groups to the BNS */
    ret = AddTGroups2TCGBnStruct( pBNS, pStruct, pVA, pTCGroups, nMaxAddEdges );
    if ( ret < 0 ) {
        goto exit_function;
    }

    /* add c-groups to the BNS; adjust charges */
    ret = AddCGroups2TCGBnStruct( pBNS, pStruct, pVA, pTCGroups, nMaxAddEdges );
    if ( ret < 0 ) {
        goto exit_function;
    }

    /* allocate BNData */
    pBD = AllocateAndInitBnData( pBNS->max_vertices + pBNS->max_vertices/2 );
    if ( !pBD ) {
        ret = BNS_OUT_OF_RAM;
        goto exit_function;
    }
    CheckBnsConsistency( pStruct, pBNS, pVA, pTCGroups, 0 );
    
    /* restore bonds & charges */
    ret = RunBnsRestore1( ip, sd, pBNS, pBD, pStruct, pVA, pTCGroups, pInChI, num_inp, bHasSomeFixedH );
    if ( ret < 0 ) {
        goto exit_function;
    }

    ret = CheckBnsConsistency( pStruct, pBNS, pVA, pTCGroups, 1 );
#if( bRELEASE_VERSION == 0 )
#ifndef INCHI_LIBRARY
    if ( ret ) {
        fprintf( stdout, "Msg for: %ld %s comp=%d %c%c\n", num_inp, (szCurHdr && szCurHdr[0])? szCurHdr : "", iComponent, pStruct->iInchiRec? 'R':'D', pStruct->iMobileH?'M':'F' ); 
    }
    if ( pStruct->iMobileH == TAUT_YES && pStruct->nNumRemovedProtons ) {
        fprintf( stdout, "REMOVED_PROTONS%+d %ld %s\n", pStruct->nNumRemovedProtons, num_inp, (szCurHdr && szCurHdr[0])? szCurHdr : "" ); 
        /*pStruct->bExtract |= EXTRACT_STRUCT_NUMBER;*/
    }
    if ( pStruct->bExtract & EXTRACT_STRUCT_NUMBER ) {
        fprintf( stdout, "EXTRACT: %ld: %s\n", num_inp, (szCurHdr && szCurHdr[0])? szCurHdr : "" );
    }
#endif
#endif
    {  /* create the final structure in pStruct->at2 */
        inp_ATOM *at_tmp = pStruct->at;
        pStruct->at = pStruct->at2;
        memcpy( pStruct->at, at_tmp, sizeof(pStruct->at[0])*(pStruct->num_atoms + pStruct->num_deleted_H) );
        ret2 = CopyBnsToAtom( pStruct, pBNS, pVA, pTCGroups, 1 );
        pStruct->at2 = pStruct->at;
        pStruct->at = at_tmp;
        if ( ret2 < 0 ) {
            ret = ret2;
        }
    }

exit_function:

    pStruct->pbfsq = NULL;
    AllocBfsQueue( &bfsq, BFS_Q_FREE, 0 );

    pBD = DeAllocateBnData( pBD );
    pBNS = DeAllocateBnStruct( pBNS );
    /*
    if ( pVA ) inchi_free( pVA );
    */
    if ( pTCGroups->pTCG ) inchi_free( pTCGroups->pTCG );

    return ret;
}
/******************************************************************************************************/
int OneInChI2Atom( ICHICONST INPUT_PARMS *ip_inp, STRUCT_DATA *sd, const char *szCurHdr, long num_inp,
                   StrFromINChI *pStruct, int iComponent, int iAtNoOffset, int bHasSomeFixedH, INChI *pInChI[])
{
    int ret;
    INPUT_PARMS *ip, ip_loc;

    ip_loc = *ip_inp;
    ip     = &ip_loc;

    sd->pStrErrStruct[0] = '\0';
    ret = RestoreAtomConnectionsSetStereo( pStruct, iComponent, iAtNoOffset, pInChI[0], pInChI[1]);
    if ( ret < 0 ) {
        goto exit_function;
    }
    ret = SetStereoBondTypesFrom0DStereo( pStruct, pInChI[0]);
    if ( ret < 0 ) {
        goto exit_function;
    }
    ret = ReconcileAllCmlBondParities( pStruct->at, pStruct->num_atoms, 0 );
    if ( ret < 0 ) {
        goto exit_function;
    }
    /* main InChI restore function */
    ret = RestoreAtomMakeBNS( ip, sd, pStruct, iComponent, iAtNoOffset, pInChI, szCurHdr, num_inp, bHasSomeFixedH );

#ifndef INCHI_ANSI_ONLY
    if ( (pStruct->num_inp_actual>0? pStruct->num_inp_actual : num_inp) >= ip->first_struct_number &&
        ( (/*ret > 0 &&*/ ip->bDisplayIfRestoreWarnings ) && pStruct->pXYZ ) ) {
        inchiTime     ulTStart;
        InchiTimeGet( &ulTStart );
        DisplayRestoredComponent( pStruct, iComponent, iAtNoOffset, pInChI[0], szCurHdr );
        sd->ulStructTime -= InchiTimeElapsed( &ulTStart ); /* subtract display time */
    }
#endif
    if ( ret < 0 ) {
        goto exit_function;
    }
    if ( (pStruct->num_inp_actual? pStruct->num_inp_actual: num_inp) >= ip->first_struct_number && ret >= 0 ) {
        /* remove t-group markings and increment zero-order bonds,
           otherwise MakeInChIOutOfStrFromINChI2() woild fail */
        /* --- moved to MakeInChIOutOfStrFromINChI2 ---
        IncrZeroBondsAndClearEndpts(pStruct->at2, pStruct->num_atoms, iComponent+1);
        CopySt2At( pStruct->at2, pStruct->st, pStruct->num_atoms );
        */
        /* include all restored structure features in pStruct->at2 */
        /* make full InChI out of pStruct->at2, pStruct->num_atoms */
        /***************************************************************************************/
        /* !!! pStruct->One_InChI etc. were removed at the exit from NormalizeAndCompare() !!! */
        /***************************************************************************************/
        if ( bHasSomeFixedH && pStruct->iInchiRec == INCHI_REC && pStruct->iMobileH == TAUT_YES &&
             !pStruct->bFixedHExists && !(ip->nMode & REQ_MODE_BASIC) ) {
            /* reconnected components without Fixed-H layer may produce 'tautomeric' fragments like Cl(-) */
            ip->nMode |= REQ_MODE_BASIC;
        }
        ret = MakeInChIOutOfStrFromINChI2( ip, sd, pStruct, iComponent, iAtNoOffset, num_inp );
        if ( ret >= 0 ) {
            ;
        } 
#if( bRELEASE_VERSION == 0 )
#ifndef INCHI_LIBRARY
        else {
            fprintf( stdout, "\nERROR in MakeInChI-1: %ld %s Comp:%d %c%c Err:%d\n", num_inp,
                     szCurHdr? szCurHdr: "???", iComponent, pStruct->iInchiRec? 'R':'D', pStruct->iMobileH?'M':'F', ret);
        }
#endif
#endif
    }


exit_function:
    return ret;
}
/********************************************************************************************/
int MakeProtonComponent( StrFromINChI *pStruct, int iComponent, int num_prot )
{
    inp_ATOM *at = NULL;
    int        i;

    if ( num_prot <= 0 ) {
        return 0;
    }
    /* allocate */
    pStruct->at  = (inp_ATOM *) inchi_calloc( num_prot, sizeof(pStruct->at[0]) );
    pStruct->at2 = (inp_ATOM *) inchi_calloc( num_prot, sizeof(pStruct->at2[0]) );
    if ( !pStruct->at || !pStruct->at2 ) {
        return 0;
    }
    /* create protons */
    at = pStruct->at;
    /* fill out proton atom info */
    for ( i = 0; i < num_prot; i ++ ) {
        strcpy( at[i].elname, "H" );
        at[i].el_number = EL_NUMBER_H;
        at[i].orig_at_number = i+1;
        /*
        at[i].orig_compt_at_numb = i + 1;
        at[i].component = i + 1;
        */
        at[i].charge    = 1;
    }
    memcpy( pStruct->at2, at, num_prot * sizeof(pStruct->at2[0]) );
    pStruct->bDeleted = 0;
    pStruct->num_atoms = num_prot;
    pStruct->bMobileH  = TAUT_YES;
    pStruct->iMobileH  = TAUT_YES;
    return num_prot; 
}
/********************************************************************************************/
int AddRemProtonsInRestrStruct( ICHICONST INPUT_PARMS *ip_inp,  STRUCT_DATA *sd, long num_inp,
                                int bHasSomeFixedH,
                                StrFromINChI *pStruct, int num_components,
                                StrFromINChI *pStructR, int num_componentsR,
                                NUM_H *nProtonsToBeRemovedByNormFromRevrs, int *recmet_change_balance )
{   /* on entry and exit, all at[i].num_H do not include isotopic H  and explicit terminal H are connected */
    int  iComp, q, ret = 0;
    int      num_atoms, tot_num_at, num_deleted_H, num_tg, num_changed, num_deleted_components;
    inp_ATOM *at;
    INPUT_PARMS *ip, ip_loc;
    int      num_prot = *nProtonsToBeRemovedByNormFromRevrs;
    int      delta_recmet_prot, num_prot_prev, bAccumulateChanges=0, nNumProtAddedByRevrs;
    INChI_Aux *pINChI_Aux;
    INCHI_MODE bNormalizationFlags;
    int        nChargeRevrs, nChargeInChI;

    if ( !num_prot ) {
        return 0;
    }
    delta_recmet_prot = 0;
    num_changed       = 0;
    num_deleted_components = 0;
    ip_loc = *ip_inp;
    ip     = &ip_loc;
    /*----------------------------------------------------------------------------------
    nLink < 0 && num_componentsR > 0 => This is a Disconnected structure component; it is
                                        same as already processed reconnected one
                                        Do no preicess it

    nLink > 0 && num_componentsR > 0 => This is a Disconnected structure component;
    (should not happen)                 It it is a result of (nLink-1)th Reconeected
                                        component disconnection (NOT IMPLEMENTED YET)

    nLink = 0                        => Process this component. It is either a reconnected
                                        component, or a result of a disconnection (for now)

    nLink > 0 && num_componentsR = 0 => This is a Reconnected component that is same as
                                        a disconnected one that will not be processed.
                                        Process and save charge delta.
    -----------------------------------------------------------------------------------*/

    for ( iComp = 0; iComp < num_components && num_prot; iComp ++ ) {
        bAccumulateChanges = 0;
        if ( pStruct[iComp].nLink < 0 && num_componentsR > 0 ) {
            /* check */
            q = -(pStruct[iComp].nLink+1);
            if ( !pStructR || !num_componentsR || q >= num_componentsR || pStructR[q].nLink != (iComp+1) ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            continue; /* Disconnected structure component has already been processed as a Reconnected one */
        }

        at           = pStruct[iComp].at2;
        num_atoms    = pStruct[iComp].num_atoms;
        tot_num_at   = pStruct[iComp].num_atoms+(num_deleted_H=pStruct[iComp].num_deleted_H);
        bAccumulateChanges = ( pStruct[iComp].nLink > 0 && !num_componentsR );
        nChargeRevrs = pStruct[iComp].nChargeRevrs;
        nChargeInChI = pStruct[iComp].nChargeInChI;
        num_deleted_components += (0 != pStruct[iComp].bDeleted);
        if ( !at || !num_atoms ) {
            continue;
        }
        /* find whether it is a reconnected structure */
        q = bRevInchiComponentExists( pStruct+iComp, INCHI_REC, TAUT_YES, 0 )? INCHI_REC : INCHI_BAS;
        /*
        q = pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC] &&
            pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC][0][TAUT_YES] &&
            pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC][0][TAUT_YES]->nNumberOfAtoms? INCHI_REC : INCHI_BAS;
        */
        pINChI_Aux = pStruct[iComp].RevInChI.pINChI_Aux[q][0][TAUT_YES]; /* 0 = 1st component in RevInChI */
        /*nNumProtAddedByRevrs = pINChI_Aux->nNumRemovedProtons;*/
        nNumProtAddedByRevrs = -pStruct[iComp].nNumRemovedProtonsByRevrs;
        bNormalizationFlags  = pINChI_Aux->bNormalizationFlags;
        num_tg               = pINChI_Aux->nNumberOfTGroups;


        /* disconnect all explicit H and add the number of implicit iso H and all explicit terminal H to the number of implicit H */
        if ( 0 > ( ret = DisconnectedConnectedH( at, num_atoms, num_deleted_H ) ) ) {
            goto exit_function;
        }
        num_prot_prev = num_prot;
        ret = AddRemoveProtonsRestr( at, num_atoms, &num_prot, nNumProtAddedByRevrs,
                                     bNormalizationFlags, num_tg, nChargeRevrs, nChargeInChI );
        
        pStruct[iComp].bPostProcessed = ret;
        num_changed += (ret > 0);
        if ( ret < 0 ) {
            goto exit_function;
        }
        if ( ret > 0 ) {
            /* recalculate InChI; it will reconnect at */
            StrFromINChI *pStruct1 = pStruct + iComp;
            INCHI_MODE    nMode = ip->nMode;
            FreeAllINChIArrays( pStruct1->RevInChI.pINChI,
                                pStruct1->RevInChI.pINChI_Aux,
                                pStruct1->RevInChI.num_components );

            if ( bHasSomeFixedH && pStruct1->iInchiRec == INCHI_REC && pStruct1->iMobileH == TAUT_YES &&
                 !pStruct1->bFixedHExists && !(ip->nMode & REQ_MODE_BASIC) ) {
                /* reconnected components without Fixed-H layer may produce 'tautomeric' fragments like Cl(-) */
                ip->nMode |= REQ_MODE_BASIC;
            }
            /* calls ConnectDisconnectedH(...): subtracts number of implicit iso H from implicit H */
            ret = MakeInChIOutOfStrFromINChI2( ip, sd, pStruct1, 0, 0, num_inp );
            ip->nMode = nMode;
            if ( ret < 0 ) {
                goto exit_function;
            }
        } else {
            /* reconnect disconnected terminal H and subtracts number of implicit iso H from implicit H */
            if ( 0 > ( ret = ConnectDisconnectedH( at, num_atoms, num_deleted_H ) ) ) {
                goto exit_function;
            }
        }
        if ( bAccumulateChanges && recmet_change_balance ) {
            /* processed Reconnected layer component that is also present in Disconnected layer */
            delta_recmet_prot += num_prot - num_prot_prev;
        }
    }

    iComp = num_components-1;
    if ( !bHasSomeFixedH && num_prot > 0 && 1 == num_deleted_components && iComp >= 0 && pStruct[iComp].bDeleted ) {
        /* add bare protons to the deleted Mobile-H component; undelete the component */
        num_prot_prev = num_prot;
        if ( !MakeProtonComponent( pStruct+iComp, iComp, num_prot ) ) {
            goto exit_function;
        } else {
            /* recalculate InChI; it will reconnect at */
            StrFromINChI *pStruct1 = pStruct + iComp;
            INCHI_MODE    nMode = ip->nMode;
            num_changed ++;
            num_prot = 0;
            FreeAllINChIArrays( pStruct1->RevInChI.pINChI,
                                pStruct1->RevInChI.pINChI_Aux,
                                pStruct1->RevInChI.num_components );

            if ( bHasSomeFixedH && pStruct1->iInchiRec == INCHI_REC && pStruct1->iMobileH == TAUT_YES &&
                 !pStruct1->bFixedHExists && !(ip->nMode & REQ_MODE_BASIC) ) {
                /* reconnected components without Fixed-H layer may produce 'tautomeric' fragments like Cl(-) */
                ip->nMode |= REQ_MODE_BASIC;
            }
            /* Although MakeInChIOutOfStrFromINChI2() calls ConnectDisconnectedH(...) */
            /* to subtracts number of implicit iso H from implicit H */
            /* this CANNOT have any effect on the deleted H component */
            ret = MakeInChIOutOfStrFromINChI2( ip, sd, pStruct1, 0, 0, num_inp );
            ip->nMode = nMode;
            if ( ret < 0 ) {
                goto exit_function;
            }
            if ( bAccumulateChanges && recmet_change_balance ) {
                /* processed Reconnected layer component that is also present in Disconnected layer */
                delta_recmet_prot += num_prot - num_prot_prev;
            }
        }
    }
    *nProtonsToBeRemovedByNormFromRevrs = num_prot;
    if ( recmet_change_balance ) {
        *recmet_change_balance = delta_recmet_prot;
    }

exit_function:
    return ret < 0? ret : num_changed;
}
/**********************************************************************************/
int AddRemIsoProtonsInRestrStruct( ICHICONST INPUT_PARMS *ip_inp,  STRUCT_DATA *sd, long num_inp, int bHasSomeFixedH,
                                StrFromINChI *pStruct, int num_components,
                                StrFromINChI *pStructR, int num_componentsR,
                                NUM_H pProtonBalance[], NUM_H recmet_change_balance[] )
{   /* on entry and exit, all at[i].num_H do not include isotopic H and explicit terminal H are connected */
    int  iComp, q, k, ret = 0, bNotEmpty;
    int      num_atoms, tot_num_at, num_deleted_H, num_tg, num_changed;
    inp_ATOM *at;
    NUM_H    num_prot[NUM_H_ISOTOPES], delta_recmet_prot[NUM_H_ISOTOPES], num_prot_prev[NUM_H_ISOTOPES];
    int      bAccumulateChanges;
    INChI_Aux *pINChI_Aux;
    INChI     *pINChI;
    INCHI_MODE bNormalizationFlags;
    INPUT_PARMS *ip, ip_loc;

    ip_loc = *ip_inp;
    ip     = &ip_loc;

    memcpy( num_prot, pProtonBalance, sizeof(num_prot) );
    for ( bNotEmpty=0, k = 0; k < NUM_H_ISOTOPES; k ++ ) {
        bNotEmpty |= num_prot[k];
    }
    if ( !bNotEmpty ) {
        return 0;
    }
    memset ( delta_recmet_prot, 0, sizeof(delta_recmet_prot));
    num_changed       = 0;
    /*----------------------------------------------------------------------------------
    nLink < 0 && num_componentsR > 0 => This is a Disconnected structure component; it is
                                        same as already processed reconnected one
                                        Do no preicess it

    nLink > 0 && num_componentsR > 0 => This is a Disconnected structure component;
    (should not happen)                 It it is a result of (nLink-1)th Reconeected
                                        component disconnection (NOT IMPLEMENTED YET)

    nLink = 0                        => Process this component. It is either a reconnected
                                        component, or a result of a disconnection (for now)

    nLink > 0 && num_componentsR = 0 => This is a Reconnected component that is same as
                                        a disconnected one that will not be processed.
                                        Process and save charge delta.
    -----------------------------------------------------------------------------------*/

    for ( iComp = 0; iComp < num_components && num_prot; iComp ++ ) {
        bAccumulateChanges = 0;
        if ( pStruct[iComp].nLink < 0 && num_componentsR > 0 ) {
            /* check */
            q = -(pStruct[iComp].nLink+1);
            if ( !pStructR || !num_componentsR || q >= num_componentsR || pStructR[q].nLink != (iComp+1) ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            continue; /* Disconnected structure component has already been processed as a Reconnected one */
        }

        at         = pStruct[iComp].at2;
        num_atoms  = pStruct[iComp].num_atoms;
        tot_num_at = pStruct[iComp].num_atoms+(num_deleted_H=pStruct[iComp].num_deleted_H);
        bAccumulateChanges = ( pStruct[iComp].nLink > 0 && !num_componentsR );
        
        if ( !at || !num_atoms ) {
            continue;
        }
        /* find whether it is a reconnected structure */
        q = pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC] &&
            pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC][0][TAUT_YES] &&
            pStruct[iComp].RevInChI.pINChI_Aux[INCHI_REC][0][TAUT_YES]->nNumberOfAtoms? INCHI_REC : INCHI_BAS;

        pINChI_Aux = pStruct[iComp].RevInChI.pINChI_Aux[q][0][TAUT_YES]; /* 0 = 1st component in RevInChI */
        pINChI     = pStruct[iComp].RevInChI.pINChI[q][0][TAUT_YES]; /* 0 = 1st component in RevInChI */
        bNormalizationFlags  = pINChI_Aux->bNormalizationFlags;
        num_tg               = pINChI_Aux->nNumberOfTGroups;
        memcpy( num_prot_prev, num_prot, sizeof(num_prot_prev) );

        /* pass CONNECTED explicit H to AddRemoveIsoProtonsRestr() for isotopic H addition */
        ret = AddRemoveIsoProtonsRestr( at, num_atoms, num_prot, num_tg );
        
        pStruct[iComp].bPostProcessed |= ret;
        num_changed += (ret > 0);
        if ( ret < 0 ) {
            goto exit_function;
        }
        if ( ret > 0 ) {
            StrFromINChI *pStruct1 = pStruct+iComp;
            INCHI_MODE    nMode = ip->nMode;
            /* recalculate InChI; MakeInChIOutOfStrFromINChI2() will reconnect explicit H */
            /* disconnect all explicit H and add the number of implicit iso H and all explicit terminal H to the number of implicit H */
            if ( 0 > ( ret = DisconnectedConnectedH( at, num_atoms, num_deleted_H ) ) ) {
                goto exit_function;
            }
            FreeAllINChIArrays( pStruct1->RevInChI.pINChI,
                                pStruct1->RevInChI.pINChI_Aux,
                                pStruct1->RevInChI.num_components );
            if ( bHasSomeFixedH && pStruct1->iInchiRec == INCHI_REC && pStruct1->iMobileH == TAUT_YES &&
                 !pStruct1->bFixedHExists && !(ip->nMode & REQ_MODE_BASIC) ) {
                /* reconnected components without Fixed-H layer may produce 'tautomeric' fragments like Cl(-) */
                ip->nMode |= REQ_MODE_BASIC;
            }
            /* input: disconnected explicit H, output: connected explicit H */
            ret = MakeInChIOutOfStrFromINChI2( ip, sd, pStruct1, 0, 0, num_inp );
            ip->nMode = nMode;
            if ( ret < 0 ) {
                goto exit_function;
            }
        }
        /* the following was commented out 2007-08-28 by DT. Reason: it's a bug since H must be already connected */
        /* else {
            if ( 0 > ( ret = ConnectDisconnectedH( at, num_atoms, num_deleted_H ) ) ) {
                goto exit_function;
            }
        } */
        if ( bAccumulateChanges ) {
            /* processed Reconnected layer component that is also present in Disconnected layer */
            for ( k = 0; k < NUM_H_ISOTOPES; k ++ ) {
                delta_recmet_prot[k] += num_prot[k] - num_prot_prev[k];
            }
        }
    }

    memcpy ( pProtonBalance, num_prot, sizeof(num_prot) );
    if ( recmet_change_balance ) {
        memcpy ( recmet_change_balance, delta_recmet_prot, sizeof(delta_recmet_prot) );
    }
exit_function:
    return ret < 0? ret : num_changed;
}

#endif
