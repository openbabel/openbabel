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

#include "inpdef.h"
#include "extr_ct.h"
#include "inpdef.h"
#include "ichitaut.h"
#include "ichinorm.h"
#include "ichicant.h"
#include "ichicomn.h"

#include "ichicomp.h"

#include "util.h"

#include "ichi_bns.h"


/* local prototypes */
int SetTautomericBonds( inp_ATOM *at, int nNumBondPos, T_BONDPOS *BondPos );
int CompRankTautomer(const void* a1, const void* a2 );
int RegisterEndPoints( T_GROUP_INFO *t_group_info, /* T_GROUP *t_group, int *pnum_t, int max_num_t,*/
                       T_ENDPOINT *EndPoint, int nNumEndPoints, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi
                       , struct BalancedNetworkStructure *pBNS );
int cmpTGroupNumber( const void *a1, const void *a2 );
int comp_candidates( const void *a1, const void *a2 );
int MoveEndpoint( inp_ATOM *at, S_CANDIDATE *s_candidate, AT_NUMB endpoint, AT_NUMB *nTGroupNewNumbers,
                  AT_NUMB *nTGroupPosition, int nNewTGroupOrd, T_GROUP_INFO *t_group_info);

int FindAccessibleEndPoints( T_ENDPOINT *EndPoint, int *nNumEndPoints, T_BONDPOS *BondPos, int *nNumBondPos,
                         struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD,
                         inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi, int taut_mode );

/* bits for GetChargeType */

#define C_SUBTYPE_CHARGED     0
#define C_SUBTYPE_p_DONOR     1  /* new */
#define C_SUBTYPE_p_ACCEPT    2  /* new */
#define C_SUBTYPE_H_ACCEPT    4
#define C_SUBTYPE_H_DONOR     8
#define C_SUBTYPE_NEUTRAL    16


/* internal stack array size */
#define MAX_STACK_ARRAY_LEN 127
#define MAX_TGROUP_ARRAY_LEN 127

/* local prototypes */
int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype );
int GetNeutralRepsIfNeeded( AT_NUMB *pri, AT_NUMB *prj, inp_ATOM *at, int num_atoms, T_ENDPOINT *EndPoint, int nNumEndPoints, C_GROUP_INFO *cgi );
int bCanBeACPoint( inp_ATOM *at, S_CHAR cCharge, S_CHAR cChangeValence, S_CHAR neutral_bonds_valence,
                   S_CHAR neutral_valence, S_CHAR nEndpointValence, S_CHAR *cChargeSubtype );
int CmpCCandidates( const void *a1, const void *a2 );
int RegisterCPoints( C_GROUP *c_group, int *pnum_c, int max_num_c, T_GROUP_INFO *t_group_info,
                     int point1, int point2, int ctype, inp_ATOM *at, int num_atoms );
int GetSaltChargeType( inp_ATOM *at, int at_no, T_GROUP_INFO *t_group_info, int *s_subtype );
int GetOtherSaltChargeType( inp_ATOM *at, int at_no, T_GROUP_INFO *t_group_info, int *s_subtype, int bAccept_O );
int MergeSaltTautGroupsBlind( inp_ATOM *at, int s_type, int num_atoms, S_GROUP_INFO *s_group_info, int nNumCandidates,
                          T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS );
int ConnectSaltTGroups2SuperTGroup( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info, int nNumCandidates,
                          T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS, int *nNewTGroupNumber, int *vertSuperTGroup );
int bDoNotMergeNonTautAtom(inp_ATOM *at, int at_no);
int GetOtherSaltType( inp_ATOM *at, int at_no, int *s_subtype );


/****************************************************************/
/*  tautomers: Sorting globals */
AT_RANK        *pn_tRankForSort;

/*************************************************************************************/
int is_centerpoint_elem( U_CHAR el_number )
{
    static U_CHAR el_numb[12];
    static int len;
    int i;
    if ( !el_numb[0] && !len ) {
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "C" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "N" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "P" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "S" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "I" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "As" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "Sb" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "Se" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "Te" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "Cl" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "Br" );
    }
    for ( i = 0; i < len; i ++ ) {
        if ( el_numb[i] == el_number ) {
            return 1;
        }
    }
    return 0;
}
#if ( KETO_ENOL_TAUT == 1 )  /* post v.1 feature */
/*************************************************************************************/
int is_centerpoint_elem_KET( U_CHAR el_number )
{
    static U_CHAR el_numb[1];
    static int len;
    int i;
    if ( !el_numb[0] && !len ) {
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "C" );
    }
    for ( i = 0; i < len; i ++ ) {
        if ( el_numb[i] == el_number ) {
            return 1;
        }
    }
    return 0;
}
#endif
/*************************************************************************************/
int is_centerpoint_elem_strict( U_CHAR el_number )
{
    static U_CHAR el_numb[6];
    static int len;
    int i;
    if ( !el_numb[0] && !len ) {
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "C" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "N" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "P" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "As" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "Sb" );
    }
    for ( i = 0; i < len; i ++ ) {
        if ( el_numb[i] == el_number ) {
            return 1;
        }
    }
    return 0;
}
/*************************************************************************************/
int get_endpoint_valence( U_CHAR el_number )
{
    static U_CHAR el_numb[6];
    static int len, len2;
    int i;
    if ( !el_numb[0] && !len ) {
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "O" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "S" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "Se" );
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "Te" );
        len2 = len;
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "N" );
    }
    for ( i = 0; i < len; i ++ ) {
        if ( el_numb[i] == el_number ) {
            return i < len2? 2 : 3;
        }
    }
    return 0;
}
#if( KETO_ENOL_TAUT == 1 )  /* post v.1 feature */
/*************************************************************************************/
int get_endpoint_valence_KET( U_CHAR el_number )
{
    static U_CHAR el_numb[2];
    static int len, len2;
    int i;
    if ( !el_numb[0] && !len ) {
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "O" );
        len2 = len;
        el_numb[len++] = (U_CHAR)get_periodic_table_number( "C" );
    }
    for ( i = 0; i < len; i ++ ) {
        if ( el_numb[i] == el_number ) {
            return i < len2? 2 : 4;
        }
    }
    return 0;
}
#endif
/********************************************************************************************************/
int AddAtom2num(  AT_RANK num[], inp_ATOM *atom, int at_no, int bSubtract )
{  /*  bSubtract: 0=> add, 1=>subtract, 2=> fill */
    inp_ATOM *at = atom + at_no;
    int       k;
    int       nMobile  =  (at->charge == -1);
    if ( bSubtract == 1 ) {
        /* 1: subtract */
        num[1] -= nMobile;
        nMobile += at->num_H;
        num[0] -= nMobile;
        for ( k = 0; k < T_NUM_ISOTOPIC; k ++ ) {
            /*  T (3H isotope) first because it has higher weight */
            num[T_NUM_NO_ISOTOPIC+k] -= at->num_iso_H[NUM_H_ISOTOPES-k-1];
        }
    } else {
        if ( bSubtract == 2 ) {
            /* fill */
            memset( num, 0, (T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC)*sizeof(num[0]) );
        }
        /* else (0): add */
        num[1] += nMobile;
        nMobile += at->num_H;
        num[0] += nMobile;
        for ( k = 0; k < T_NUM_ISOTOPIC; k ++ ) {
            /*  T (3H isotope) first because it has higher weight */
            num[T_NUM_NO_ISOTOPIC+k] += at->num_iso_H[NUM_H_ISOTOPES-k-1];
        }
    }
    return nMobile;
}
/********************************************************************************************************/
void AddAtom2DA( AT_RANK num_DA[], inp_ATOM *atom, int at_no, int bSubtract )
{   /*  bSubtract: 0=> add, 1=>subtract, 2=> fill */
    inp_ATOM *at = atom + at_no;
    int       nDelta, nAcidic_O;

    if (at->charge < -1 || at->charge == 1 && !at->c_point || at->charge > 1 )
        return;

    nDelta = ( bSubtract == 1 )? -1 : 1;

    /* "Acidic" O, S, Se, Te recognition */
    if ( at->at_type & ATT_ACIDIC_CO ) {
        nAcidic_O = nDelta;
    } else {
        nAcidic_O = 0;
    }

    if ( bSubtract == 2 ) { /* 2: fill, otherwise add */
        memset( num_DA, 0, TG_NUM_DA * sizeof(num_DA[0]) );
    }
    if ( at->charge <= 0 && at->valence == at->chem_bonds_valence ||
         /* neutral or negative donor */
         at->charge > 0  && at->valence + 1 == at->chem_bonds_valence
         /* positively charged donor */ 
       ) {
        if ( at->charge < 0 ) {
            num_DA[TG_Num_dM] += nDelta;
            num_DA[TG_Num_dO] += nAcidic_O;
        } else
        if ( at->num_H ) {
            num_DA[TG_Num_dH] += nDelta;
            num_DA[TG_Num_dO] += nAcidic_O;
        }
    } else
    if ( at->charge <= 0 && at->valence + 1 == at->chem_bonds_valence ||
         at->charge > 0  && at->valence + 2 == at->chem_bonds_valence ) {
        /* acceptor */
        if ( at->charge < 0 ) {
            num_DA[TG_Num_aM] += nDelta;
        } else
        if ( at->num_H ) {
            num_DA[TG_Num_aH] += nDelta;
        } else {
            num_DA[TG_Num_aO] += nAcidic_O; /* acidic O-acceptor has no H or charge */
        }
    }
    return;
}
/********************************************************************************************************/
int AddEndPoint( T_ENDPOINT *pEndPoint, inp_ATOM *at, int iat )
{        
    pEndPoint->nAtomNumber  = iat;
    pEndPoint->nEquNumber   = 0;
    pEndPoint->nGroupNumber = at[iat].endpoint;
    if ( at[iat].endpoint ) {
        /* already an endpoint */
        memset( pEndPoint->num, 0, sizeof(pEndPoint->num) );
    } else {
        /* not an endpoint yet, make it an endpoint */
        AddAtom2num( pEndPoint->num, at, iat, 2 );  /* fill */
        AddAtom2DA( pEndPoint->num_DA, at, iat, 2 );
        /*
        nMobile  = pEndPoint->num[1] = (at[iat].charge == -1);
        nMobile  = pEndPoint->num[0] = at[iat].num_H + nMobile;
        for ( k = 0; k < T_NUM_ISOTOPIC; k ++ ) {
            pEndPoint->num[T_NUM_NO_ISOTOPIC+k] = at[iat].num_iso_H[NUM_H_ISOTOPES-k-1];
        }
        */
    }
    return 0;
}
/********************************************************************************************************/
int nGetEndpointInfo( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif )
{
    int  nEndpointValence;
    int  nMobile;
    S_CHAR cChargeSubtype;

    if ( atom[iat].radical && atom[iat].radical != RADICAL_SINGLET )
        return 0; /* a radical */
    if ( !(nEndpointValence = get_endpoint_valence( atom[iat].el_number )) )
        return 0; /* not an endpoint */
    if ( nEndpointValence <= atom[iat].valence )
        return 0; /* not an endpoint, for example >N(+)< or >N<  or >O(+)- or >O- or >N- or -O- */

    if ( atom[iat].charge == -1 || atom[iat].charge == 0 ) {
        /* not a positive charge-point */
        if ( nEndpointValence < atom[iat].chem_bonds_valence )
            return 0; /* abnormal valence > standard endpoint valence */
        nMobile = atom[iat].num_H + (atom[iat].charge == -1);
        if ( nMobile + atom[iat].chem_bonds_valence != nEndpointValence )
            return 0; /* non-standard endpoint valence */
        switch ( atom[iat].chem_bonds_valence - atom[iat].valence ) {
        case 0:
            eif->cDonor    = 1;
            eif->cAcceptor = 0;
            break;
        case 1:
            eif->cDonor    = 0;
            eif->cAcceptor = 1;
            break;
        default:
            return 0;
        }
        eif->cMobile = nMobile;
        eif->cNeutralBondsValence = nEndpointValence-nMobile;
        eif->cMoveableCharge = 0;
#if ( KETO_ENOL_TAUT == 1 )
        eif->cKetoEnolCode = 0;
#endif
        return nEndpointValence;
    } else
    if ( atom[iat].c_point &&
         0 <= GetChargeType( atom, iat, &cChargeSubtype ) &&
         ((int)cChargeSubtype & (C_SUBTYPE_H_ACCEPT|C_SUBTYPE_H_DONOR))
        ) {
        /* charge-point */
        if ( cChargeSubtype & C_SUBTYPE_H_ACCEPT ) {
            eif->cDonor    = 0;
            eif->cAcceptor = 1;
        } else
        if ( cChargeSubtype & C_SUBTYPE_H_DONOR ) {
            eif->cDonor    = 1;
            eif->cAcceptor = 0;
        } else {
            return 0;
        }
        eif->cMobile = atom[iat].num_H;
        eif->cNeutralBondsValence = nEndpointValence-atom[iat].num_H;
        eif->cMoveableCharge = atom[iat].charge;
#if ( KETO_ENOL_TAUT == 1 )
        eif->cKetoEnolCode = 0;
#endif
        return nEndpointValence;
    }
    return 0;
}
#if ( KETO_ENOL_TAUT == 1 )  /* post v.1 feature */
/********************************************************************************************************/
int nGetEndpointInfo_KET( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif )
{
    int  nEndpointValence;
    int  nMobile;
    S_CHAR cChargeSubtype;

    /*
    static U_CHAR el_number_O, el_number_C;

    if ( !el_number_O ) {
        el_number_O = (U_CHAR)get_periodic_table_number( "O" );
        el_number_C = (U_CHAR)get_periodic_table_number( "C" );
    }
    */
    if ( atom[iat].radical && atom[iat].radical != RADICAL_SINGLET )
        return 0; /* a radical */
    if ( !(nEndpointValence = get_endpoint_valence_KET( atom[iat].el_number )) )
        return 0; /* not an endpoint; only O and C can be an endpoint for keto-enol tautomerism */
    if ( nEndpointValence <= atom[iat].valence )
        return 0; /* not an endpoint, for example >N(+)< or >N<  or >O(+)- or >O- or >N- or -O- */
    if ( nEndpointValence == 4 && atom[iat].valence < 2 )
        return 0; /* exclude O==C--CH3  <=> HO--C==CH2 */
    if ( nEndpointValence == 2 && atom[iat].valence > 1 )
        return 0; /* exclude --O--C==CH-- */

    if ( atom[iat].charge == -1 || atom[iat].charge == 0 ) {
        /* not a positive charge-point */
        if ( nEndpointValence < atom[iat].chem_bonds_valence )
            return 0; /* abnormal valence > standard endpoint valence */
        nMobile = atom[iat].num_H + (atom[iat].charge == -1);
        if ( nMobile + atom[iat].chem_bonds_valence != nEndpointValence )
            return 0; /* non-standard endpoint valence */
        switch ( atom[iat].chem_bonds_valence - atom[iat].valence ) {
        case 0:
            eif->cDonor    = 1;
            eif->cAcceptor = 0;
            break;
        case 1:
            eif->cDonor    = 0;
            eif->cAcceptor = 1;
            break;
        default:
            return 0;
        }
        eif->cMobile = nMobile;
        eif->cNeutralBondsValence = nEndpointValence-nMobile;
        eif->cMoveableCharge = 0;
        eif->cKetoEnolCode   = (nEndpointValence == 2)? 1 : (nEndpointValence == 4)? 2 : 0;
        return nEndpointValence;
    } else
    if ( atom[iat].c_point &&
         0 <= GetChargeType( atom, iat, &cChargeSubtype ) &&
         ((int)cChargeSubtype & (C_SUBTYPE_H_ACCEPT|C_SUBTYPE_H_DONOR))
        ) {
        /* charge-point; currently only O for keto-enol tautomerism */
        if ( cChargeSubtype & C_SUBTYPE_H_ACCEPT ) {
            eif->cDonor    = 0;
            eif->cAcceptor = 1;
        } else
        if ( cChargeSubtype & C_SUBTYPE_H_DONOR ) {
            eif->cDonor    = 1;
            eif->cAcceptor = 0;
        } else {
            return 0;
        }
        eif->cMobile = atom[iat].num_H;
        eif->cNeutralBondsValence = nEndpointValence-atom[iat].num_H;
        eif->cMoveableCharge = atom[iat].charge;
        eif->cKetoEnolCode   = (nEndpointValence == 2)? 1 : (nEndpointValence == 4)? 2 : 0;
        return nEndpointValence;
    }
    return 0;
}
#endif
/********************************************************************************************************/
/*  RegisterEndPoints ret>0 => new registration happened, 0 => no changes, -1 => program error (debug) */
int RegisterEndPoints( T_GROUP_INFO *t_group_info, /* T_GROUP *t_group, int *pnum_t, int max_num_t,*/
                       T_ENDPOINT *EndPoint, int nNumEndPoints, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi
                       , struct BalancedNetworkStructure *pBNS )
{
    T_GROUP  *t_group   =  t_group_info->t_group;
    int      *pnum_t    = &t_group_info->num_t_groups;
    int       max_num_t =  t_group_info->max_num_t_groups;
    int       nNumZeroEqu, nNumNewTGroups;
    AT_NUMB   group, prev_group, prev_eqnum, nNextGroupNumber, nLeastGroupNumber; 
    int       nNumGroups, num_t, difference;
    int       i, j, k, ret;
    AT_NUMB   nNewTgNumberStackArray[MAX_STACK_ARRAY_LEN+1];
    AT_NUMB   nGroupNumberStackArray[MAX_STACK_ARRAY_LEN+1];
    AT_NUMB   nGroupNewNumberStackArray[MAX_STACK_ARRAY_LEN+1];
    AT_NUMB  *nNewTgNumber    = nNewTgNumberStackArray;
    AT_NUMB  *nGroupNumber    = nGroupNumberStackArray;
    AT_NUMB  *nGroupNewNumber = nGroupNewNumberStackArray;

    if ( nNumEndPoints <= 0 )
        return 0; /* nothing to do */
    num_t            = *pnum_t;
    difference       = 0;
    nNextGroupNumber = 0;
    nNumZeroEqu      = 0;
    ret              = 0;
    /* find max group number; increment it to obtain next available group number */
    for ( i = 0; i < num_t; i ++ ) {
        if ( nNextGroupNumber < t_group[i].nGroupNumber )
            nNextGroupNumber = t_group[i].nGroupNumber;
    }
    nNextGroupNumber ++;

    /* find min non-zero group number nLeastGroupNumber;
       count zero EndPoint[i].nEquNumber
       if all EndPoint[i].nGroupNumber are equal and non-zero then exit: nothing to do.
    */
    nLeastGroupNumber = nNextGroupNumber;
    prev_group        = EndPoint[0].nGroupNumber;
    prev_eqnum        = EndPoint[0].nEquNumber;
    for ( i = j = k = 0; i < nNumEndPoints; i ++ ) {
        if ( group = EndPoint[i].nGroupNumber ) {
            if ( group < nLeastGroupNumber ) {
                nLeastGroupNumber = group;
            }
        }
        j += (prev_group == EndPoint[i].nGroupNumber); /* count endpoints that belong to the 1st group */
        k += (prev_eqnum == EndPoint[i].nEquNumber);   /* count endpoints that belongo to a group equivalent to the 1st group */
        nNumZeroEqu += !EndPoint[i].nEquNumber;        /* count endpoints that have been processed by FindAccessibleEndPoints() */
    }
    if ( j == nNumEndPoints && prev_group && k == nNumEndPoints ) {
        /* all endpoints already belong to one t-group;
           the last comparison is not needed for now 
           because EndPoint[i].nEquNumber cannot make
           endpont partitioning finer
        */
        return 0;
    }
         
    nNumNewTGroups   = 0;

    if ( !nNumZeroEqu ) {
        /* EndPoint[] has been processed by FindAccessibleEndPoints;
         * equal EndPoint[i].nEquNumber mark endpoints belonging to
         * the same t-group
         * Since now the next available t-group number, nNextGroupNumber,
         * is known,replace fict. IDs assigned by FindAccessibleEndPoints
         * with correct new t-group numbers.
         */
        for ( i = 0; i < nNumEndPoints; i ++ ) {
            if ( (group = EndPoint[i].nEquNumber) >= nNextGroupNumber ) {
                /* replace fict. IDs assigned by FindAccessibleEndPoints() with new t-group numbers */
                /* these fict. IDs have values = (num_atoms+1), (num_atoms+2),...; they may be non-contiguous */
                for ( j = 0; j < nNumNewTGroups; j ++ ) {
                    if ( group == nGroupNewNumber[j] )
                        break;
                }
                if ( j == nNumNewTGroups ) {
                    /* found new fict. ID = group */
                    if ( j == MAX_STACK_ARRAY_LEN && nGroupNewNumber == nGroupNewNumberStackArray ) {
                        /* stack array overflow; allocate more memory than may be needed */
                        nGroupNewNumber = (AT_NUMB *)inchi_malloc(nNumEndPoints*sizeof(nGroupNewNumber[0]));
                        if ( !nGroupNewNumber ) {
                            ret = -1;
                            goto exit_function;
                        }
                        memcpy( nGroupNewNumber, nGroupNewNumberStackArray, nNumNewTGroups*sizeof(nGroupNewNumber[0]));
                    }
                    /* save newly found fict. t-group ID to compare to the next values of EndPoint[].nEquNumber  */
                    nGroupNewNumber[j] = group;
                    nNumNewTGroups ++;
                }
                EndPoint[i].nEquNumber = nNextGroupNumber + j;
            }
        } /* after this point the values just stored in nGroupNewNumber[] will not
             be used. However, the obtained nNumNewTGroups value will be used */
    } else
    if ( nNumZeroEqu == nNumEndPoints ) {
        /* EndPoint[] has NOT been processed by FindAccessibleEndPoints;
           all atoms and t-groups to which endpoints belong should be merged into a single t-group
         */
        if ( nLeastGroupNumber == nNextGroupNumber ) {
            /* flag to create a new t-group: none of the found
             *  endpoints belong to an already known t-group
             */
            nNumNewTGroups = 1; /* otherwise 0 */
        }
        /* All EndPoint[*].nEquNumber are zeroes. All endpoints will
         * belong to one new or old t-group; its ID is nLeastGroupNumber.
         * Set  EndPoint[i].nEquNumber = nLeastGroupNumber;
         */
        for ( i = 0; i < nNumEndPoints; i ++ ) {
            EndPoint[i].nEquNumber = nLeastGroupNumber;
        }
    } else {
        ret = -1; /* program error: only some of EndPoint[i].nEquNumber are zero */ /*   <BRKPT> */
        goto exit_function;
    }

    if ( nNumNewTGroups ) {
        /* create new nNumNewTGroups t-group(s) */
        if ( num_t + nNumNewTGroups > max_num_t ) {
            ret = -1; /* found too many t-groups */ /*   <BRKPT> */
            goto exit_function;
        }
        /* initialize new t-group(s) */
        memset( t_group + num_t, 0, nNumNewTGroups * sizeof(t_group[0]) );
        for ( i = 0; i < nNumNewTGroups; i ++ ) {
            t_group[num_t+i].nGroupNumber = nNextGroupNumber + i;
        }
    }

    /* At this point:
     *   EndPoint[i].nGroupNumber == 0   => the endpoint atom does not belong to a t-group yet
     *   EndPoint[i].nGroupNumber  > 0   => current t-group ID of the endpoint atom
     *   EndPoint[i].nEquNumber      -->    new ID of a tautomeric group of this endpoint atom
     *   EndPoint[i].nAtomNumber     -->    number of the endpoint atom
     */

     nNumGroups = 0; /* counts the groups to be renumbered */
     for ( i = j = 0; i < nNumEndPoints; i ++ ) {
         if ( group = EndPoint[i].nGroupNumber ) {
             if ( group == EndPoint[i].nEquNumber ) {
                continue; /* ignore: the endpoint belongs to the same t-group as before */
             }
             /* save information for renumbering of the existing t-groups */
             for ( j = 0; j < nNumGroups; j ++ ) {
                 if ( group == nGroupNumber[j] ) {
                     if ( EndPoint[i].nEquNumber != nGroupNewNumber[j] ) {
                         ret = -1; /* program error */ /*   <BRKPT> */
                         goto exit_function;
                     }
                     break;
                 }
             }
             if ( j == nNumGroups ) {
                 /* discovered a new t-group number; store it together with its nEquNumber */
                 if ( j == MAX_STACK_ARRAY_LEN ) {
                    if ( nGroupNewNumber == nGroupNewNumberStackArray ) {
                        nGroupNewNumber = (AT_NUMB *)inchi_malloc(nNumEndPoints*sizeof(nGroupNewNumber[0]));
                        if ( !nGroupNewNumber ) {
                            ret = -1;
                            goto exit_function;
                        }
                        memcpy( nGroupNewNumber, nGroupNewNumberStackArray, nNumGroups*sizeof(nGroupNewNumber[0]));
                    }
                    if ( nGroupNumber == nGroupNumberStackArray ) {
                        nGroupNumber = (AT_NUMB *)inchi_malloc(nNumEndPoints*sizeof(nGroupNumber[0]));
                        if ( !nGroupNumber ) {
                            ret = -1;
                            goto exit_function;
                        }
                        memcpy( nGroupNumber, nGroupNumberStackArray, nNumGroups*sizeof(nGroupNumber[0]));
                    }
                 }

                 nGroupNumber[j] = group;  /* old t-group ID */
                 nGroupNewNumber[j] = EndPoint[i].nEquNumber; /* new t-group ID */
                 nNumGroups ++;
             }
         } else {
            /* add a new endpoint to the newly created or previously existing t-groups */
             group = EndPoint[i].nEquNumber;
             if ( group >= nNextGroupNumber ) {
                 /* get index of a new t-group from equ number */
                 j = num_t + group - nNextGroupNumber; /* newly assigned IDs are contiguous */
             } else {
                 /* old t-group */
                 if ( j >= num_t || group != t_group[j].nGroupNumber ) {
                     /* search only if j is not a needed group index */
                     for ( j = 0; j < num_t; j ++ ) {
                         if ( group == t_group[j].nGroupNumber )
                             break;
                     }
                     if ( j == num_t ) {
                         ret = -1; /* program error: t-group not found */ /*   <BRKPT> */
                         goto exit_function;
                     }
                 }
             }
             /* add aton to existing or new t-group */
             t_group[j].nNumEndpoints ++;
             for ( k = 0; k < (int)(sizeof(t_group->num)/sizeof(t_group->num[0])); k ++ )
                 t_group[j].num[k]    += EndPoint[i].num[k];
             for ( k = 0; k < (int)(sizeof(t_group->num_DA)/sizeof(t_group->num_DA[0])); k ++ )
                 t_group[j].num_DA[k] += EndPoint[i].num_DA[k];
             /* mark endpoint */
             at[EndPoint[i].nAtomNumber].endpoint = group;
             difference ++;
         }
     }

     difference += nNumGroups;
     num_t      += nNumNewTGroups;
     if ( !difference ) {
         ret = 0; /* nothing to do. Not necessarily a program error: happens if all EndPoint[i].nGroupNumber==EndPoint[i].nEquNumber  */
         goto exit_function;
     }
     
     if ( nNumGroups ) {
         /* prepare for renumbering: find max t-group number */
         for ( i = 0, nNextGroupNumber = 0; i <  num_t; i ++ ) {
             if ( nNextGroupNumber < t_group[i].nGroupNumber ) {
                 nNextGroupNumber = t_group[i].nGroupNumber;
             }
         }
     }
     /* renumber and merge t-groups */
     for ( i = 0; i < nNumGroups; i ++ ) {
         int i1, i2;
         AT_NUMB group1 = nGroupNumber[i];
         AT_NUMB group2 = nGroupNewNumber[i];
         /* add group1 to group2, then delete group1. */
         for ( j = 0, i1 = i2 = -1; j < num_t && (i1 < 0 || i2 < 0); j ++ ) {
             if ( i1 < 0 && group1 == t_group[j].nGroupNumber )
                 i1 = j;
             if ( i2 < 0 && group2 == t_group[j].nGroupNumber )
                 i2 = j;
         }
         if ( i1 < 0 || i2 < 0 ) {
             ret = -1; /* program error */ /*   <BRKPT> */
             goto exit_function;
         }
         /*  add t_group[i1] to t_group[i2] and remove t_group[i1] */
         for ( k = 0; k < (int)(sizeof(t_group->num)/sizeof(t_group->num[0])); k ++ )
             t_group[i2].num[k]    += t_group[i1].num[k];
         for ( k = 0; k < (int)(sizeof(t_group->num_DA)/sizeof(t_group->num_DA[0])); k ++ )
             t_group[i2].num_DA[k] += t_group[i1].num_DA[k];
         t_group[i2].nNumEndpoints += t_group[i1].nNumEndpoints;
         num_t --;
         if ( num_t > i1 ) {
             memmove( t_group+i1, t_group+i1+1, ( num_t - i1)*sizeof(t_group[0]) );
         }
     }

     if ( nNumGroups ) {
         /* there are groups to merge */
         if ( nNextGroupNumber >= MAX_STACK_ARRAY_LEN ) {
             nNewTgNumber = (AT_NUMB *)inchi_malloc((nNextGroupNumber+1)*sizeof(*nNewTgNumber));
             if ( !nNewTgNumber ) {
                 ret = -1;
                 goto exit_function; /* error: out of RAM */
             }
         }
         memset( nNewTgNumber, 0, (nNextGroupNumber+1)*sizeof(*nNewTgNumber) );
         for ( i = 0; i < num_t; i ++ ) {
             nNewTgNumber[t_group[i].nGroupNumber] = i+1; /* new t-group numbers */
         }
         for ( j = 0; j < nNumGroups; j ++ ) {
             if ( !nNewTgNumber[nGroupNumber[j]] && nNewTgNumber[nGroupNewNumber[j]] ) {
                     nNewTgNumber[nGroupNumber[j]] = nNewTgNumber[nGroupNewNumber[j]];
             } else {
                 ret = -1; /* program error: all new numbers must have been marked */
                 goto exit_function;
             }
         }
         /* renumber t-groups */
         for ( i = 0; i < num_t; i ++ ) {
             t_group[i].nGroupNumber = nNewTgNumber[t_group[i].nGroupNumber];
         }
#if( bRELEASE_VERSION != 1 )
         /* Check: debug only */
         for ( i = 1; i < num_t; i ++ ) {
             if ( 1 != t_group[i].nGroupNumber - t_group[i-1].nGroupNumber ) {
                 ret = -1; /* debug */
                 goto exit_function;
             }
         }
#endif
         /* renumber endpoints */
         for ( i = 0; i < num_atoms; i ++ ) {
             if ( group = at[i].endpoint ) {
                 if ( !(at[i].endpoint = nNewTgNumber[group]) || nNextGroupNumber <= nNewTgNumber[group] ) {
                     ret = -1; /* program error */
                     goto exit_function;
                 }
             }
         }
     }
     if ( nNewTgNumber != nNewTgNumberStackArray ) {
         inchi_free( nNewTgNumber );
         nNewTgNumber = nNewTgNumberStackArray;
     }
     if ( nGroupNumber != nGroupNumberStackArray ) {
         inchi_free(nGroupNumber);
         nGroupNumber = nGroupNumberStackArray;
     }
     if ( nGroupNewNumber != nGroupNewNumberStackArray ) {
         inchi_free( nGroupNewNumber );
         nGroupNewNumber = nGroupNewNumberStackArray;
     }
     if ( !t_group_info->tGroupNumber ) {
        t_group_info->tGroupNumber = (AT_NUMB *)inchi_malloc(2*max_num_t*sizeof(t_group_info->tGroupNumber[0])); 
        if ( !t_group_info->tGroupNumber ) {
            ret = -1;
            goto exit_function;
        }
     }
     /* fill out t-group index 2004-02-27 */
     memset( t_group_info->tGroupNumber, 0, 2*max_num_t*sizeof(t_group_info->tGroupNumber[0]) );
     for ( i = 0; i < num_t; i ++ ) {
         if ( t_group[i].nNumEndpoints && t_group[i].nGroupNumber )
            t_group_info->tGroupNumber[t_group[i].nGroupNumber] = i+1;
     }

     if ( pBNS && (pBNS->tot_st_cap == pBNS->tot_st_flow || ALWAYS_ADD_TG_ON_THE_FLY) ) {
        T_GROUP_INFO tgi;
        int ret_bns;
        memset( &tgi, 0, sizeof(tgi) );
        tgi.num_t_groups = num_t;
        tgi.t_group = t_group;
#if( KETO_ENOL_TAUT == 1 )
        tgi.bTautFlags |= (t_group_info->bTautFlags & TG_FLAG_KETO_ENOL_TAUT); /* needed in AddTGroups2BnStruct() */
#endif
        /* reinitialize BN Structure */
        ret_bns = ReInitBnStruct( pBNS, at, num_atoms, 0 );
        if ( IS_BNS_ERROR( ret_bns ) ) {
            return ret_bns;
        }
        if ( *pBNS->pbTautFlags & TG_FLAG_MOVE_POS_CHARGES ) {
            /* set new charge groups */
            ret_bns = AddCGroups2BnStruct( pBNS, at, num_atoms, cgi );
            if ( IS_BNS_ERROR( ret_bns ) ) {
                return ret_bns;
            }
        }
        /* set new tautomeric groups */
        ret_bns = AddTGroups2BnStruct( pBNS, at, num_atoms, &tgi );
        if ( IS_BNS_ERROR( ret_bns ) ) {
            return ret_bns;
        }
     }

     *pnum_t = num_t;
     return difference;

exit_function:
     if ( nNewTgNumber != nNewTgNumberStackArray ) {
         inchi_free( nNewTgNumber );
     }
     if ( nGroupNumber != nGroupNumberStackArray ) {
         inchi_free(nGroupNumber);
     }
     if ( nGroupNewNumber != nGroupNewNumberStackArray ) {
         inchi_free( nGroupNewNumber );
     }
     return ret;
}
/*******************************************************************************************************
 * change non-alternating and non-tautomeric bonds
 * (that is, single and double bonds) to tautomeric
 */
int SetTautomericBonds( inp_ATOM *at, int nNumBondPos, T_BONDPOS *BondPos )
{
    int k, n;
    for ( k = n = 0; k < nNumBondPos; k ++ ) {
        int neighbor_index = BondPos[k].neighbor_index;
        int center         = BondPos[k].nAtomNumber;
        int bond_mark      = at[center].bond_type[neighbor_index];
        int bond_type      = bond_mark & ~BOND_MARK_ALL;
        int neighbor;
#if( REPLACE_ALT_WITH_TAUT == 1 )
        if ( bond_type != BOND_TAUTOM )
#else
        if ( bond_type != BOND_ALTERN && bond_type != BOND_TAUTOM )
#endif
        {
            int ii;
            /*  change bond type to BOND_TAUTOM presering higher bits marks */
            bond_type = (bond_mark & BOND_MARK_ALL) | BOND_TAUTOM;
            /*  change center-neighbor bond */
            at[center].bond_type[neighbor_index] = bond_type;
            neighbor = at[center].neighbor[neighbor_index];
            for ( ii = 0; ii < at[neighbor].valence; ii ++ ) {
                if ( at[neighbor].neighbor[ii] == center ) {
                    /*  neighbor-center bond found */
                    at[neighbor].bond_type[ii] = bond_type;
                    break;
                }
            }
            n ++;
        }
    }
    return n;
}

/********************************************************************************************************/
int GetNeutralRepsIfNeeded( AT_NUMB *pri, AT_NUMB *prj, inp_ATOM *at, int num_atoms, T_ENDPOINT *EndPoint, int nNumEndPoints, C_GROUP_INFO *cgi )
{
    AT_NUMB ri = *pri;
    AT_NUMB rj = *prj;
    int     i, k;
    AT_NUMB c_point, endpoint, r;

    if ( (c_point = at[ri].c_point) && (c_point == at[rj].c_point) && 
         (at[ri].charge == 1 || at[rj].charge == 1) && cgi && cgi->num_c_groups > 0 ) {
        /* at[ri] and at[rj] belong to the same charge group, at least one is charged */
        /* MS VC++ 2005 reports unreachable code here ??? */
        for ( k = 0; k < cgi->num_c_groups; k ++ ) {
            if ( cgi->c_group[k].nGroupNumber == c_point ) {
                /* cgi->c_group[k] is found to be this charge group */
                if ( cgi->c_group[k].num_CPoints - cgi->c_group[k].num[0] < 2 ) {
                    /* Only one neutral in the c-group: we will not be able to neutralize both
                       when looking for the alt path to discover the tautomerism.
                       Therefore we need to find a neutral t-group representative */
                    /* at[rj] */
                    if ( endpoint = at[rj].endpoint ) {
                        for ( i = 0; i < nNumEndPoints; i ++ ) {
                            if ( (r=EndPoint[i].nAtomNumber) == *prj )
                                continue; /* ignore at[*prj] */
                            if ( at[r].endpoint != endpoint )
                                continue; /* at[r] does not belong to the same t-group as at[*prj]; ignore the atom */
                            if ( !at[r].c_point ) {
                                rj = r; /* found a neutral t-group representative */
                                break;
                            }
                            if ( at[r].c_point != c_point && c_point == at[rj].c_point ) {
                                /* replace only once because of (c_point == at[rj].c_point) condition  */
                                rj = r;
                            }
                        }
                        if ( rj == *prj /*&& at[ri].endpoint*/ ) {
                            /* !!! "&& at[ri].endpoint": only between 2 t-groups 2004-02-27;
                                the change disabled due to undiscovered yet possibility of ambiguity*/
                            /* no replacement has been found in EndPoint[]; try all atoms in the t-group */
                            for ( i = 0; i < num_atoms; i ++ ) {
                                if ( at[i].endpoint != endpoint )
                                    continue;
                                if ( i == (int)*prj )
                                    continue;
                                if ( !at[i].c_point ) {
                                    rj = (AT_NUMB)i; /* found neutral t-group representative */
                                    break;
                                }
                                if ( at[i].c_point != c_point && c_point == at[rj].c_point ) {
                                    /* replace only once */
                                    rj = (AT_NUMB)i;
                                }
                            }
                        }
                    }
                    /* at[ri] */
                    if ( endpoint = at[ri].endpoint ) {
                        for ( i = 0; i < nNumEndPoints; i ++ ) {
                            if ( (r=EndPoint[i].nAtomNumber) == *pri )
                                continue;
                            if ( at[r].endpoint != endpoint )
                                continue;
                            if ( !at[r].c_point ) {
                                ri = r; /* found neutral t-group representative */
                                break;
                            }
                            if ( at[r].c_point != c_point && c_point == at[ri].c_point &&
                                 at[r].c_point != at[rj].c_point ) {
                                /* replace only once */
                                ri = r;
                            }
                        }
                        if ( ri == *pri && at[rj].endpoint ) {
                            /* !!! "&& at[rj].endpoint": only between 2 t-groups 2004-02-27;
                               the change disabled due to undiscovered yet possibility of ambiguity */
                            for ( i = 0; i < num_atoms; i ++ ) {
                                if ( at[i].endpoint != endpoint )
                                    continue;
                                if ( i == (int)*pri )
                                    continue;
                                if ( !at[i].c_point ) {
                                    ri = (AT_NUMB)i; /* found neutral t-group representative */
                                    break;
                                }
                                if ( at[i].c_point != c_point && c_point == at[ri].c_point &&
                                     at[i].c_point != at[rj].c_point) {
                                    /* replace only once */
                                    ri = (AT_NUMB)i;
                                }
                            }
                        }
                    }

                }
            }
            break;
        }
        *prj = rj;
        *pri = ri;
    }
    return 0;
}

/********************************************************************************************************/
int FindAccessibleEndPoints( T_ENDPOINT *EndPoint, int *nNumEndPoints, T_BONDPOS *BondPos, int *nNumBondPos,
                         struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD,
                         inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi, int taut_mode )
{
    AT_NUMB nTGroupRepresenative[MAXVAL], nTGroupEqu[MAXVAL], nTGEndPointNo[MAXVAL], ri, rj;
    AT_NUMB nCurTGroupNumber, nMaxTGroupNumber, nNumTgroupNumbers, nMaxEquNumber;
    int   i, j, k, nNumDiffTGroupNumbers = 0, nNumFoundEqu, nErr;
    
    if ( *nNumEndPoints != *nNumBondPos )
        return 0;
    /* collect all group numbers. Fill EndPoint[i].nEquNumber */
    for ( i = 0; i < *nNumEndPoints; i ++ ) {
        nCurTGroupNumber = EndPoint[i].nEquNumber = EndPoint[i].nGroupNumber; /* initial equivalence */
        if ( nCurTGroupNumber ) {
            /* found endpoint that already belongs to a t-group */
            for ( j = 0; j < nNumDiffTGroupNumbers; j ++ ) {
                if ( nTGroupEqu[j] == nCurTGroupNumber )
                    break;
            }
            if ( j == nNumDiffTGroupNumbers ) {
                nTGroupRepresenative[nNumDiffTGroupNumbers] = EndPoint[i].nAtomNumber;
                nTGroupEqu[nNumDiffTGroupNumbers]           = EndPoint[i].nGroupNumber;
                nTGEndPointNo[nNumDiffTGroupNumbers]        = i;
                nNumDiffTGroupNumbers ++;
            }
        }
    }


    /* check whether each pair belongs to the same t-group and establish the equivalence(s) */
    for ( i = 0, nNumFoundEqu=0; i < nNumDiffTGroupNumbers; i ++ ) {
        for ( j = i+1; j < nNumDiffTGroupNumbers; j ++ ) {
            ri = nTGroupRepresenative[i];
            rj = nTGroupRepresenative[j];
            /* both at[ri] and at[rj] are known to belong to tautomeric groups */
            GetNeutralRepsIfNeeded( &ri, &rj, at, num_atoms, EndPoint, *nNumEndPoints, cgi );
            nErr = bExistsAnyAltPath( pBNS, pBD, at, num_atoms, ri, rj, taut_mode );
            if ( IS_BNS_ERROR(nErr) )
                return nErr;
            if ( 0 == nErr )
                continue; /* alt path between at[ri] and at[rj] not found */
            nCurTGroupNumber = inchi_min( nTGroupEqu[i], nTGroupEqu[j] );
            nMaxTGroupNumber = inchi_max( nTGroupEqu[i], nTGroupEqu[j] );
            for ( k = 0; k < nNumDiffTGroupNumbers; k ++ ) {
                if ( nTGroupEqu[k]==nMaxTGroupNumber ) {
                    nTGroupEqu[k] = nCurTGroupNumber;
                    nNumFoundEqu ++;
                }
            }
            for ( k = 0; k < *nNumEndPoints; k ++ ) {
                if ( EndPoint[k].nEquNumber == nMaxTGroupNumber ) {
                    EndPoint[k].nEquNumber = nCurTGroupNumber;
                }
            }
        }
    }
    if ( nNumFoundEqu ) {
        /* leave in only non-equivalent representatives */
        for ( i = 1, k = 0; i < nNumDiffTGroupNumbers; i ++ ) {
            for ( j = 0; j < i; j ++ ) {
                if ( nTGroupEqu[j] == nTGroupEqu[i] ) {
                    nTGroupEqu[i] = 0;  /* i > j; mark equivalent for removal*/
                    break;
                }
            }
        }
        for ( i = j = 0; i < nNumDiffTGroupNumbers; i ++ ) {
            if ( nTGroupEqu[i] ) {
                if ( i != j ) { /* remove the marked */
                    nTGroupEqu[j]           = nTGroupEqu[i];
                    nTGroupRepresenative[j] = nTGroupRepresenative[i];
                    nTGEndPointNo[j]        = nTGEndPointNo[i];
                }
                j ++;
            }
        }
        nNumDiffTGroupNumbers = j; /* number of known t-group representatives */
    }
    /* collect endpoints that have not been assigned to t-groups */
    for ( i = 0, j = nNumDiffTGroupNumbers; i < *nNumEndPoints; i ++ ) {
        if ( EndPoint[i].nEquNumber )
            continue;
        nTGroupEqu[j] = 0;
        nTGroupRepresenative[j] = EndPoint[i].nAtomNumber;
        nTGEndPointNo[j]        = i;
        j ++;

    }
    nNumTgroupNumbers = j;
    nMaxEquNumber = num_atoms + 1; /* impossible atom or t-group number */

    /* check whether each pair belongs to the same group and establish the equivalence(s) */
    for ( i = 0, nNumFoundEqu=0; i < nNumTgroupNumbers; i ++ ) {
        for ( j = i+1; j < nNumTgroupNumbers; j ++ ) {
            if ( nTGroupEqu[i] != nTGroupEqu[j] && (i>=nNumDiffTGroupNumbers || j>=nNumDiffTGroupNumbers) ||
                 /* equivalence of a t-group and a non-t-group atom */
                 !nTGroupEqu[i] && !nTGroupEqu[j]
                 /* equivalence of two non-t-group atoms */
               ) {
                ri = nTGroupRepresenative[i];
                rj = nTGroupRepresenative[j];
                /*------------------------------!!!---------------------------------------------  
                    Explanation why GetNeutralRepsIfNeeded() may need to be changed 2004-02-27
                    The change has been disabled due to undiscovered yet possibility of ambiguity
                    to search for neutral only among EndPoint[] in case taut-not_taut pairs
                
                    Counterexample:   O=C-NH(+)=C-NH2
                                      1   2       3

                    Has already been found: 2-3 (+)-charge exchange
                                            1-2 tautomerism (charge removed to 3)
                    Now testing: 2-3 tautomerism. If not commented out,
                       GetNeutralRepsIfNeeded() would replace 2-3 test with 1-3 test because:
                        o Charge group has only one neutral and both 2 and 3 belong to it,
                          therefore we cannot neutralize both; search for neutral representative;
                        o Since 1 and 2 belong to the same t-group and 1 is neutral,
                          test 1-3 instead of 2-3.
                    This breaks our condition:
                        Test tautomeric H movement only between neutral atoms.
                 -----------------------------------------------------------------------------*/
                GetNeutralRepsIfNeeded( &ri, &rj, at, num_atoms, EndPoint, *nNumEndPoints, cgi );
                
                nErr = bExistsAnyAltPath( pBNS, pBD, at, num_atoms, ri, rj, taut_mode );
                if ( IS_BNS_ERROR(nErr) )
                    return nErr;
                if ( nErr <= 0 )
                    continue;
                if ( nTGroupEqu[i] && nTGroupEqu[j] ) {
                    /* found equivalence of two t-groups; at least one of them must be a new one */
                    nCurTGroupNumber = inchi_min( nTGroupEqu[i], nTGroupEqu[j] );
                    nMaxTGroupNumber = inchi_max( nTGroupEqu[i], nTGroupEqu[j] );
                    for ( k = 0; k < nNumTgroupNumbers; k ++ ) {
                        if ( nTGroupEqu[k]==nMaxTGroupNumber ) {
                            nTGroupEqu[k] = nCurTGroupNumber;
                            nNumFoundEqu ++;
                        }
                    }
                    for ( k = 0; k < *nNumEndPoints; k ++ ) {
                        if ( EndPoint[k].nEquNumber == nMaxTGroupNumber ) {
                            EndPoint[k].nEquNumber = nCurTGroupNumber;
                        }
                    }
                } else
                if ( nTGroupEqu[i] ) { /* extend existing t-group */
                    nTGroupEqu[j] = nTGroupEqu[i];
                    EndPoint[nTGEndPointNo[j]].nEquNumber = nTGroupEqu[i];

                } else
                if ( nTGroupEqu[j] ) { /* extend existing t-group */
                    nTGroupEqu[i] = nTGroupEqu[j];
                    EndPoint[nTGEndPointNo[i]].nEquNumber = nTGroupEqu[j];

                } else { /* establis a new t-group */
                    nTGroupEqu[i] =
                    nTGroupEqu[j] = nMaxEquNumber; /* assign a fict. ID to establish equivalence */
                    EndPoint[nTGEndPointNo[i]].nEquNumber =
                    EndPoint[nTGEndPointNo[j]].nEquNumber = nMaxEquNumber;
                    nMaxEquNumber ++;
                }
            }
        }
    }
    /* eliminate endpoints and bonds that do not belong to t-group(s) 
       (they have not been found connected by an alt path to any other endpoint)
    */
    for ( i = 0, j = 0; i < *nNumEndPoints; i ++ ) {
        if ( EndPoint[i].nEquNumber ) {
#if( IGNORE_SINGLE_ENDPOINTS == 1 )  /* 1-28-2003 */
            for ( k = 0, nNumFoundEqu = 0; k < *nNumEndPoints; k ++ ) {
                nNumFoundEqu += (EndPoint[i].nEquNumber == EndPoint[k].nEquNumber);
            }
            if ( nNumFoundEqu <= 1 ) { /* one time it is equal to itself when i == k above */
                /* if EndPoint[i] is not "equivalent" to any other EndPoint then ignore it */
                continue;
            }
#endif
            if ( i != j ) { /* save endpoints that are found to be connected to other endpoints by alt paths */
                EndPoint[j] = EndPoint[i];
                BondPos[j]  = BondPos[i];
            }
            j ++;
        }
    }

#if( IGNORE_SINGLE_ENDPOINTS != 1 )  /* 1-28-2003 */
    /* Do not allow a centerpoint to have only one tautomeric bond */
    /* Hack: we may have only one centerpoint */
    /* BondPos[*].nAtomNumber are centerpoints */
    if ( j == 1 ) {
        /* check if there exist other centerpoint neighbors
         * connected to it by another tautomeric-bond
         */
        for ( i = 0, k = 0; i < at[BondPos[0].nAtomNumber].valence; i ++ ) {
            k += ( i != BondPos[0].neighbor_index &&
                   BOND_TAUTOM == (at[BondPos[0].nAtomNumber].bond_type[i] & ~BOND_MARK_ALL));
        }
        if ( !k ) {
            j = 0;
        }
    }
#endif

    *nNumEndPoints = *nNumBondPos = j;
    return j;
    
}

/*#if( MOVE_CHARGES == 1 ) */  /* { */
/********************************************************************************************************/

/**********************************************/
/*                                            */
/* definitions for positive ion recognition   */
/*                                            */
/**********************************************/

typedef struct tagChargeType {  /* meaning see in bCanBeACPoint() */
    char    elname[3];
    S_CHAR  charge;
    S_CHAR  neutral_valence;
    S_CHAR  neutral_bonds_valence; /* valence of a neutral atom */
    S_CHAR  cChangeValence; /* charge increases valence by this value */
    S_CHAR  cChargeType;    /* different types are treated separately */
    S_CHAR  num_bonds;      /* added 02-06-2005 */
} CHARGE_TYPE;

CHARGE_TYPE CType[] = {
    { "N\0",  1, 3, 3, 1, 0, 0 },
    { "P\0",  1, 3, 3, 1, 1, 0 },
#if( ADD_MOVEABLE_O_PLUS == 1 )
    { "O\0",  1, 2, 2, 1, 2, 2 }, /* added 02-06-2005 */
    { "S\0",  1, 2, 2, 1, 3, 2 }, /* added 03-18-2005 */
    { "Se",   1, 2, 2, 1, 4, 2 }, /* added 03-18-2005 */
    { "Te",   1, 2, 2, 1, 5, 2 }, /* added 03-18-2005 */
#endif
};

/* bits */

#define C_SUBTYPE_CHARGED     0
#define C_SUBTYPE_p_DONOR     1  /* new */
#define C_SUBTYPE_p_ACCEPT    2  /* new */
#define C_SUBTYPE_H_ACCEPT    4
#define C_SUBTYPE_H_DONOR     8
#define C_SUBTYPE_NEUTRAL    16

/* make sure any C_SUBTYPE_CHARGED_... < any C_SUBTYPE_NEUTRAL_... */
/* charged */
#define C_SUBTYPE_CHARGED_NON_TAUT          (C_SUBTYPE_CHARGED)
#define C_SUBTYPE_CHARGED_p_DONOR           (C_SUBTYPE_CHARGED|C_SUBTYPE_p_DONOR)
#define C_SUBTYPE_CHARGED_H_ACCEPT          (C_SUBTYPE_CHARGED|C_SUBTYPE_H_ACCEPT)
#define C_SUBTYPE_CHARGED_H_ACCEPT_p_DONOR  (C_SUBTYPE_CHARGED|C_SUBTYPE_H_ACCEPT|C_SUBTYPE_p_DONOR)
#define C_SUBTYPE_CHARGED_H_DONOR           (C_SUBTYPE_CHARGED|C_SUBTYPE_H_DONOR |C_SUBTYPE_p_DONOR)
/* neutral */
#define C_SUBTYPE_NEUTRAL_NON_TAUT          (C_SUBTYPE_NEUTRAL)
#define C_SUBTYPE_NEUTRAL_H_ACCEPT          (C_SUBTYPE_NEUTRAL|C_SUBTYPE_H_ACCEPT)
#define C_SUBTYPE_NEUTRAL_H_ACCEPT_p_ACCEPT (C_SUBTYPE_NEUTRAL|C_SUBTYPE_H_ACCEPT|C_SUBTYPE_p_ACCEPT)
#define C_SUBTYPE_NEUTRAL_H_DONOR           (C_SUBTYPE_NEUTRAL|C_SUBTYPE_H_DONOR)

#define NUM_C_TYPES  (int)(sizeof( CType )/sizeof(CType[0]))


/********************************************************************************************************/
int bCanBeACPoint( inp_ATOM *at, S_CHAR cCharge, S_CHAR cChangeValence, S_CHAR neutral_bonds_valence,
                   S_CHAR neutral_valence, S_CHAR nEndpointValence, S_CHAR *cChargeSubtype )
{
    int nChangeValence;
    int nNumBonds;
    int nBondsValence;
    int bNegCharge = (at->charge == -1);  /* add fict. bonds to (-) 2004-02-24*/
    if ( at->charge == cCharge && at->valence == at->chem_bonds_valence && at->num_H ) {
        /* proton donors candidates >NH(+)-, >NH2(+), -NH3(+), >OH(+), -OH2(+) */
        /* charged, added p-transfer -- 01-28-2004 */
        nChangeValence = at->charge * cChangeValence; /* +1 or -1; currently only +1 */
        nBondsValence  = at->chem_bonds_valence + at->num_H;
        if ( nBondsValence == neutral_bonds_valence + nChangeValence && nEndpointValence ) {
            *cChargeSubtype = C_SUBTYPE_CHARGED_p_DONOR; /* ignore Phosphorus p-donors for now */
        }
        return 0;
    } else
    if ( at->charge == cCharge && at->valence < at->chem_bonds_valence ) {
        /* the requirement at->valence < at->chem_bonds_valence rejects
           candidates >NH(+)-, >NH2(+), -NH3(+), >N(+)<, >OH(+), -OH2(+), >O(+)-
           Moveable charge requires double bonds; these ions have no double bonds
         */

        /* charged */
        nChangeValence = at->charge * cChangeValence; /* +1 or -1; currently only +1 */
        nBondsValence  = at->chem_bonds_valence + at->num_H;
        nNumBonds      = at->valence + at->num_H;
        if ( nBondsValence == neutral_bonds_valence + nChangeValence ) { /* known valence */
            if ( nNumBonds == neutral_valence ) {
                /* non-tautomeric: >N(+)=, =O(+)-
                   possibly tautomeric donor: =NH(+)-, =NH2(+), =OH(+) */
                if ( at->valence == neutral_valence || !nEndpointValence ) {
                    /* non-tautomeric: >N(+)=, =O(+)-; any suitable P+: >P(+)=, =PH(+)-, =PH2(+) */
                    *cChargeSubtype = C_SUBTYPE_CHARGED_NON_TAUT;
                } else {
                    /* possibly tautomeric donor: =NH(+)-, =NH2(+), =OH(+) */
                    *cChargeSubtype = C_SUBTYPE_CHARGED_H_DONOR;
                }
                return 1;
            }
            if ( nNumBonds == neutral_valence - 1 ) {
                /* possibly tutomeric acceptor: =N(+)=, #N(+)-, #NH(+), #O(+) */
                if ( nEndpointValence ) {
                    *cChargeSubtype = at->num_H? C_SUBTYPE_CHARGED_H_ACCEPT_p_DONOR : C_SUBTYPE_CHARGED_H_ACCEPT;
                } else {
                    /* =P(+)=, #P(+)-, #PH(+) */
                    *cChargeSubtype = C_SUBTYPE_CHARGED_NON_TAUT;
                }
                return 1; /* charge type, charged */
            }
        }

    } else
    if ( at->charge == 0 || bNegCharge ) {
        /* neutral atom or anion, all bonds are single */
        nBondsValence  = at->chem_bonds_valence + at->num_H + bNegCharge; /* add fict. bonds to (-) 2004-02-24*/
        nNumBonds      = at->valence + at->num_H + bNegCharge; /* add fict. bonds to (-) 2004-02-24*/
        if ( nBondsValence == neutral_bonds_valence ) {
            if ( nNumBonds == neutral_valence ) {
                /* only single bonds: >N-, >NH, -NH2, -O-, -OH, >P- >PH -PH2 */
                /*                    >N(-), -NH(-), -O(-). >P(-) -PH(-) */
                if ( at->valence == neutral_valence || !nEndpointValence ) {
                    /* >N-, -O-, any P(3 single bonds): >P- >PH -PH2  */
                    *cChargeSubtype = C_SUBTYPE_NEUTRAL_NON_TAUT;
                } else
                if ( at->valence < neutral_valence /*&& nEndpointValence */ ) {
                    /* num_H > 0: >NH -NH2 -OH */
                    /* num_H = 0: none C_SUBTYPE_NEUTRAL_H_ACCEPT for now */
                    *cChargeSubtype = at->num_H? C_SUBTYPE_NEUTRAL_H_DONOR: C_SUBTYPE_NEUTRAL_H_ACCEPT;
                } else {
                    return 0;
                }
                return 1; /* charge type, neutral */
            }
            if ( nNumBonds == neutral_valence - 1 ) {
                /* possibly tautomeric acceptor =N-, =NH, =O or non-taut =P-, =PH */
                if ( nEndpointValence ) {
                    /* =N-,  =NH, =O  */
                    *cChargeSubtype = C_SUBTYPE_NEUTRAL_H_ACCEPT_p_ACCEPT;
                } else {
                    /* =P-, =PH */
                    *cChargeSubtype = C_SUBTYPE_NEUTRAL_NON_TAUT;
                }
                return 1; /* charge type, (+) => neutral */
            }
        }
    }
    return 0;
}
/********************************************************************************************************/
int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype )
{
    int i, n;
    S_CHAR    nEndpointValence;
    inp_ATOM *at = atom + iat;

    *cChargeSubtype = 0;
    /* ignore ion pairs and charges != 1 */
    if ( abs(at->charge) == 1 ) {
        for ( i = 0; i < at->valence; i ++ ) {
            n = at->neighbor[i];
            /* allow negatively charged tautomeric neighbors 2004-02-26 */
            if ( abs(atom[n].charge + at->charge) < abs(atom[n].charge - at->charge) && !atom[n].endpoint ) {
                return -1; /* charges have different signs */
            }
        }
    } else
    if ( at->charge ) {
        return -1; /* abs(charge) != 1 */
    }
    /* find candidates */
    for ( i = 0; i < NUM_C_TYPES; i ++ ) {
        if ( !strcmp( at->elname, CType[i].elname ) &&
             (!CType[i].num_bonds || CType[i].num_bonds==at->valence && at->nNumAtInRingSystem >= 5) ) {
            nEndpointValence = (S_CHAR)get_endpoint_valence(at->el_number );
            if ( bCanBeACPoint( at, CType[i].charge, CType[i].cChangeValence, CType[i].neutral_bonds_valence,
                                CType[i].neutral_valence, nEndpointValence, cChargeSubtype ) ) {
                return CType[i].cChargeType;
            }
        }
    }
    return -1;
}
/********************************************************************************************************/
int CmpCCandidates( const void *a1, const void *a2 )
{
    const C_CANDIDATE *c1 = (const C_CANDIDATE *)a1;
    const C_CANDIDATE *c2 = (const C_CANDIDATE *)a2;
    int ret;
    if ( ret = (int)c1->type - (int)c2->type )
        return ret;
    if ( ret = (int)c1->subtype - (int)c2->subtype )
        return ret;
    ret = (int)c1->atnumber - (int)c2->atnumber;
    return ret;
}
/********************************************************************************************************/
int RegisterCPoints( C_GROUP *c_group, int *pnum_c, int max_num_c, T_GROUP_INFO *t_group_info,
                     int point1, int point2, int ctype, inp_ATOM *at, int num_atoms )
{
    int num_c = *pnum_c, i, i1, i2;
    AT_NUMB nGroupNumber = 0, nNewGroupNumber;


    if ( at[point1].c_point == at[point2].c_point ) {
        if ( at[point1].c_point )
            return 0;
        memset( c_group+num_c, 0, sizeof(c_group[0]) );
        if ( num_c < max_num_c ) {
            c_group[num_c].num[0] = CHARGED_CPOINT(at,point1) + CHARGED_CPOINT(at, point2);
            c_group[num_c].num_CPoints += 2;
            c_group[num_c].cGroupType   = ctype;
            /* get next available c-group number */
            for ( i = 0; i < num_c; i ++ ) {
                if ( nGroupNumber < c_group[i].nGroupNumber )
                    nGroupNumber = c_group[i].nGroupNumber;
            }
            nGroupNumber ++;
            c_group[num_c].nGroupNumber =
            at[point1].c_point          =
            at[point2].c_point          = nGroupNumber;
            *pnum_c = num_c+1;
            /* count protons */
            if ( at[point1].num_H ) {
                c_group[num_c].num[1] ++;
            } else
            if ( at[point2].num_H ) {
                c_group[num_c].num[1] ++;
            } else
            if ( (at[point1].endpoint || at[point2].endpoint) && t_group_info && t_group_info->t_group && t_group_info->num_t_groups ) {
            /* !!! add later !!! */
            }


            return 1;
        }
        return BNS_CPOINT_ERR; /* overflow */
    }
    if ( at[point1].c_point > at[point2].c_point ) {
        /* make sure at[point1].c_point < at[point2].c_point */
        i = point1;
        point1 = point2;
        point2 = i;
    }
    if ( !at[point1].c_point ) {
        /* add a new c-endpoint to an existing c-group */
        nGroupNumber = at[point2].c_point;
        for ( i = 0; i < num_c; i ++ ) {
            if ( nGroupNumber == c_group[i].nGroupNumber ) {
                at[point1].c_point = at[point2].c_point;
                c_group[i].num_CPoints ++;
                c_group[i].num[0] += CHARGED_CPOINT(at,point1);
                return 1;
            }
        }
        return BNS_CPOINT_ERR; /* program error: c-group not found */
    } else {
        /* merge two c-groups */
        nNewGroupNumber = at[point1].c_point;
        nGroupNumber    = at[point2].c_point;
        for ( i = 0, i1=i2=-1; i < num_c && (i1 < 0 || i2 < 0); i ++ ) {
            if ( nNewGroupNumber == c_group[i].nGroupNumber ) {
                i1 = i;
                continue;
            }
            if ( nGroupNumber    == c_group[i].nGroupNumber ) {
                i2 = i;
                continue;
            }
        }
        if ( i1 < 0 || i2 < 0 ) {
            return BNS_CPOINT_ERR; /* at least one not found */
        }

        c_group[i1].num[0]      += c_group[i2].num[0];
        c_group[i1].num_CPoints += c_group[i2].num_CPoints;
        num_c --;
        if ( num_c > i2 ) {
            memmove( c_group+i2, c_group+i2+1, ( num_c - i2)*sizeof(c_group[0]) );
        }
        *pnum_c = num_c;
        /* renumber c-groups */
        for ( i = 0; i < num_c; i ++ ) {
            if ( c_group[i].nGroupNumber > nGroupNumber ) {
                c_group[i].nGroupNumber --;
            }
        }
        /* renumber c-points */
        for ( i = 0; i < num_atoms; i ++ ) {
            if ( at[i].c_point > nGroupNumber ) {
                at[i].c_point --;
            } else
            if ( at[i].c_point == nGroupNumber ) {
                at[i].c_point = nNewGroupNumber;
            }
        }
        return 1;
    }
}

/********************************************************************************************************/
int MarkChargeGroups(inp_ATOM *at, int num_atoms, 
                     C_GROUP_INFO *c_group_info, T_GROUP_INFO *t_group_info,
                     struct BalancedNetworkStructure *pBNS, 
                     struct BalancedNetworkData *pBD)
{   
    int nNumChanges = 0;

    if ( c_group_info && c_group_info->c_candidate && c_group_info->max_num_candidates > 0 ) {
        int i, i1, i2, i3, j, num_tested;
        C_CANDIDATE *c_candidate       = c_group_info->c_candidate;
        int          nMaxNumCandidates = c_group_info->max_num_candidates;
        int          nNumCandidates    = c_group_info->num_candidates;
        S_CHAR       c_type, c_subtype;
        int          iat1, iat2, ret, nDelta;

        if ( nNumCandidates == -1 ) 
        {
            nNumCandidates = 0; /* 2004-02-26 they could appear after t-group discovery */
            /*return 0;*/
        }
        if ( nNumCandidates == 0 ) 
        {
            for ( i = 0, nNumCandidates = 0; i < num_atoms; i ++ ) 
            {
                if ( 0 <= (c_type = GetChargeType( at, i, &c_subtype )) ) 
                {
                    if ( nNumCandidates >= nMaxNumCandidates ) 
                    {
                        return BNS_VERT_EDGE_OVFL;
                    }
                    c_candidate[nNumCandidates].atnumber = i;
                    c_candidate[nNumCandidates].type     = c_type;
                    c_candidate[nNumCandidates].subtype  = c_subtype;
                    nNumCandidates ++;
                }
            }
            if ( nNumCandidates <= 1 ) 
            {
                c_group_info->num_candidates = -1; /* no candidate exists */
                return 0;
            }
        }
        /* sorting keys: (1) atom type (N,P); (2) uncharged=16/charged=0; (3) other;
           atom-charged-N .... i1
              ...
           atom-charged-N
           atom-neutral-N .... i2
              ...
           atom-neutral-N
           atom-charged-P .... i3 ... i1
              ...
           atom-charged-P
           atom-neutral-P ........... i2
              ...
           atom-neutral-P
           end.           ........... i3
        */
        qsort(c_candidate, nNumCandidates, sizeof(c_candidate[0]), CmpCCandidates);
    
        i1 = 0;
        num_tested = 0;
        nDelta     = 0;

        while ( i1 < nNumCandidates ) 
        {
    
            /* the the first charged candidate of a new atom type */
            for (; i1 < nNumCandidates && (c_candidate[i1].subtype & C_SUBTYPE_NEUTRAL); i1 ++)
                ;
            if ( i1 == nNumCandidates )
                break; /* not found */
        
            /* bypass other charged candidates of the same atom type */
            for ( i2 = i1+1; i2 < nNumCandidates &&
                            c_candidate[i2].type    == c_candidate[i1].type &&
                            !(c_candidate[i2].subtype & C_SUBTYPE_NEUTRAL); i2++ )
                ;
            if ( i2 == nNumCandidates )
                break; /* no neutral candidates */
        
            /* find next to the last neutral candidate of the same atom type */
            for ( i3 = i2;  i3 < nNumCandidates &&
                            c_candidate[i3].type    == c_candidate[i1].type; i3 ++ )
                ; 
        
            if ( i3 == i2 ) 
            {
                /* no neutral candidates found */
                if ( i2 < nNumCandidates ) 
                {
                    i1 = i3;
                    continue;  /* move to the next atom type */
                }
                break; /* nothing more to do */
            }

            /* found charged candidates: i1...i2-1;  neutral candidates: i2...i3-1 */
            for ( i = i1; i < i2; i ++ ) 
            {
                iat1 = c_candidate[i].atnumber;
                for ( j = i2; j < i3; j ++ ) 
                {
                    /* check alt path at[iat1]=-=-...-at[iat2]; at[iat1] is charged, at[iat2] is neutral */
                    num_tested ++;
                    iat2 = c_candidate[j].atnumber;
                    if ( at[iat1].c_point && at[iat1].c_point == at[iat2].c_point )
                        continue;
                    ret = bExistsAltPath( pBNS, pBD, NULL, at, num_atoms, iat1, iat2, ALT_PATH_MODE_CHARGE );
                    if ( IS_BNS_ERROR( ret ) ) 
                    {
                        return ret;
                    }
                    if ( ret & 1 ) 
                    {
                        nDelta       = (ret & ~3) >> 2;
                        nNumChanges += (ret & 2);
                        ret = RegisterCPoints( c_group_info->c_group, &c_group_info->num_c_groups,
                                               c_group_info->max_num_c_groups, t_group_info,
                                               iat1, iat2, c_candidate[i1].type, at, num_atoms );
                        if ( IS_BNS_ERROR( ret ) ) 
                        {
                            return ret;
                        }
                        if ( nDelta ) 
                        {
                            goto quick_exit;
                        }
                    }
                }
            }
            i1 = i3;
        }
quick_exit:
        if ( c_group_info->num_candidates == 0 ) 
        {
            /* first time: initialize */
            c_group_info->num_candidates = num_tested? nNumCandidates : -1; /* no candidate exists */
        }
        
    }
    return nNumChanges;
}

/********************************************************************************************************/
int GetSaltChargeType(inp_ATOM *at, int at_no, T_GROUP_INFO *t_group_info, int *s_subtype )
{
    static int el_number_C  = 0;
    static int el_number_O  = 0;
    static int el_number_S  = 0;
    static int el_number_Se = 0;
    static int el_number_Te = 0;

/* 
   type (returned value):
    -1 => ignore
     0 => oxygen
   subtype:
     1 = SALT_DONOR_H   => has H 
     2 = SALT_DONOR_Neg => has (-) charge
     4 = SALT_ACCEPTOR  => may be an acceptor of H or (-), but not necessarily

   O-atom should be:
     - a terminal atom
     - connected to unsaturated, uncharged, non-radical atom C that has chemical valence 4:
     H-donors:             =CH-OH, =C(-X)-OH
     possible H-acceptors: -CH=O, >C=O
     H-acceptors are true if O is tautomeric
*/
    int iC, tg, i, type;
    /* one-time initialization */
    if ( !el_number_O ) {
        el_number_C  = get_periodic_table_number( "C" );
        el_number_O  = get_periodic_table_number( "O" );
        el_number_S  = get_periodic_table_number( "S" );
        el_number_Se = get_periodic_table_number( "Se" );
        el_number_Te = get_periodic_table_number( "Te" );
    }
    *s_subtype = 0; /* initialize the output */
    /* check whether it is a candidate */
    if ( at[at_no].valence != 1 ||
         at[at_no].radical && at[at_no].radical != RADICAL_SINGLET ||
         at[at_no].charge < -1 ||
         at[at_no].charge > 0 && !at[at_no].c_point ) {
        return -1;
    }
    
    if ( at[at_no].el_number == el_number_O  ||
         at[at_no].el_number == el_number_S  ||
         at[at_no].el_number == el_number_Se ||
         at[at_no].el_number == el_number_Te ) {
        type = 0;  /* terminal oxygen atom, needs more to be checked... */
    } else {
        type = -1; /* ignore this atom */
    }

    if ( type < 0 ||
         at[at_no].chem_bonds_valence + at[at_no].num_H !=
         get_el_valence(at[at_no].el_number, at[at_no].charge, 0) ) {
        return -1; /* non-standard valence or not an oxygen */
    }
    
    iC = at[at_no].neighbor[0];

#if ( SALT_WITH_PROTONS == 1 )
    if ( at[iC].el_number != el_number_C ||
         at[iC].chem_bonds_valence + at[iC].num_H != 4 || /* allow =C(H)-OH or -C(H)=O */
         at[iC].charge         ||
         at[iC].radical && at[iC].radical != RADICAL_SINGLET ||
         at[iC].valence == at[iC].chem_bonds_valence ) {
        return -1; /* oxigen is connected to a wrong atom */
    }
#else
    if ( at[iC].el_number != el_number_C ||
         at[iC].num_H ||
         at[iC].chem_bonds_valence != 4 ||  /* allow only no H on C */
         at[iC].charge         ||
         at[iC].radical && at[iC].radical != RADICAL_SINGLET ||
         at[iC].valence == at[iC].chem_bonds_valence ) {
        return -1; /* oxigen is connected to a wrong atom */
    }
#endif
    if ( (tg = at[at_no].endpoint) && t_group_info && t_group_info->t_group ) {
        /* O-atom is in a tautomeric group */
        for ( i = 0; i < t_group_info->num_t_groups; i ++ ) {
            if ( tg == t_group_info->t_group[i].nGroupNumber ) {
                /*
                t_group_info->t_group[i].num[0] = number of attached H-atoms and negative charges
                t_group_info->t_group[i].num[1] = number of attached negative charges
                */
                if ( t_group_info->t_group[i].num[0] > t_group_info->t_group[i].num[1] ) {
                    *s_subtype |= SALT_DONOR_H; /* has H */
                }
                if ( t_group_info->t_group[i].num[1] ) {
                    *s_subtype |= SALT_DONOR_Neg; /* has (-) */
                }
                *s_subtype |= SALT_ACCEPTOR; /* there is always an acceptor in a t-group */
                return type;
            }
        }
        return -1; /* error: t-group not found */
    }
    /* O is not not in a tautomeric group */
    /* assume valence(O-) < valence(O) < valence(O+) */
    if ( at[at_no].charge == -1 ) {
        *s_subtype |= SALT_DONOR_Neg; /* has (-) */
    }
    if ( at[at_no].charge <= 0 && at[at_no].num_H ) {
        *s_subtype |= SALT_DONOR_H; /* has H */
    }
    if ( at[at_no].charge == 0 && at[at_no].chem_bonds_valence == 2 ) {
        *s_subtype |= SALT_ACCEPTOR;
    }
    /* since O cannot be a charge point, the following cannot happen: */
    if ( at[at_no].charge == 1 && at[at_no].c_point && at[at_no].chem_bonds_valence == 2 && at[at_no].num_H ) {
        *s_subtype |= SALT_DONOR_H; /* has H */
    }
    return type;
}
/********************************************************************************************************/
int bDoNotMergeNonTautAtom(inp_ATOM *at, int at_no)
{
    static int el_number_N  = 0;

    if ( !el_number_N ) {
        el_number_N  = get_periodic_table_number( "N" );
    }
    if ( at[at_no].el_number == el_number_N )
    {
        return 1;
    }
    return 0;
}
/********************************************************************************************************/
int GetOtherSaltChargeType( inp_ATOM *at, int at_no, T_GROUP_INFO *t_group_info, int *s_subtype, int bAccept_O )
{
   /* static int el_number_C  = 0; */
   /* static int el_number_N  = 0; */
    static int el_number_O  = 0;
    static int el_number_S  = 0;
    static int el_number_Se = 0;
    static int el_number_Te = 0;

/* 
   type (returned value):
    -1 => ignore
     1 => not an oxygen
   subtype:
     1 = SALT_DONOR_H   => has H 
     2 = SALT_DONOR_Neg => has (-) charge
     4 = SALT_ACCEPTOR  => may be an acceptor of H or (-), but not necessarily

   the atom should be:
     - a tautomeric endpoint atom
     - connected to possible centerpoint atom

   another description of the atom searched here:

   any possibly tautomeric atom adjacent to a possibly centerpoint
      that has at least one double bond (possibly if positively charged);
   if eif.cAcceptor then the bond between the atom and the centerpoint must be possibly double
   if eif.cAcceptor then the bond must be possibly single
   Donors that belong to a t-group are also acceptors


*/
    int tg, i, j, type, endpoint_valence, num_centerpoints, bond_type, centerpoint;
    ENDPOINT_INFO eif;
    /* one-time initialization */
    if ( !el_number_O && !bAccept_O ) {
       /* el_number_C  = get_periodic_table_number( "C" ); */
       /* el_number_N  = get_periodic_table_number( "N" ); */
        el_number_O  = get_periodic_table_number( "O" );
        el_number_S  = get_periodic_table_number( "S" );
        el_number_Se = get_periodic_table_number( "Se" );
        el_number_Te = get_periodic_table_number( "Te" );
    }
    *s_subtype = 0; /* initialize the output */
    if ( !bAccept_O /* only N */ && 
         (at[at_no].el_number == el_number_O  ||
          at[at_no].el_number == el_number_S  ||
          at[at_no].el_number == el_number_Se ||
          at[at_no].el_number == el_number_Te ) ) {
        return -1; /* we are not looking for oxygen here */
    }

    type = 1;
    if ( !(endpoint_valence = nGetEndpointInfo( at, at_no, &eif )) ) {
        return -1; /* not a possible endpoint */
    } else {
        /* at[at_no] is not not in a tautomeric group; use eif previously filled out by nGetEndpointInfo */
        /* check whether there is adjacent atom-candidate for a centerpoint */
        num_centerpoints = 0;
        for ( j = 0; j < at[at_no].valence; j ++ ) {
            bond_type   = (int)at[at_no].bond_type[j] & BOND_TYPE_MASK;
            centerpoint = (int)at[at_no].neighbor[j];  /*  a centerpoint candidate */
            if ( ( eif.cAcceptor && (bond_type == BOND_DOUBLE  || 
                                     bond_type == BOND_ALTERN  || /* possibly double */
                                     bond_type == BOND_ALT12NS ||
                                     bond_type == BOND_TAUTOM   )  ||
                   eif.cDonor    && (bond_type == BOND_SINGLE  ||
                                     bond_type == BOND_ALTERN  || /* possibly single */
                                     bond_type == BOND_ALT12NS ||
                                     bond_type == BOND_TAUTOM   )  ) &&
                   (at[centerpoint].chem_bonds_valence >  at[centerpoint].valence ||
                   /* check for possible endpoint added 2004-02-24 */
                    at[centerpoint].chem_bonds_valence == at[centerpoint].valence &&
                    (at[centerpoint].endpoint || at[centerpoint].c_point) /* tautomerism or charge may increment at[centerpoint].chem_bonds_valence*/ ) && 
                   is_centerpoint_elem( at[centerpoint].el_number ) ) {
                num_centerpoints ++;
                break; /* at least one possibly centerpoint neighbor has been found */
            }
        }
        if ( !num_centerpoints ) {
            return -1;
        }
        /* moved here from just after "type = 1;" line 2004-02-26 */
        if ( (tg = at[at_no].endpoint) && t_group_info && t_group_info->t_group ) {
            /* atom is in a tautomeric group */
            for ( i = 0; i < t_group_info->num_t_groups; i ++ ) {
                if ( tg == t_group_info->t_group[i].nGroupNumber ) {
                    /*
                    t_group_info->t_group[i].num[0] = number of attached H-atoms and negative charges
                    t_group_info->t_group[i].num[1] = number of attached negative charges
                    */
                    if ( t_group_info->t_group[i].num[0] > t_group_info->t_group[i].num[1] ) {
                        *s_subtype |= SALT_DONOR_H; /* has H */
                    }
                    if ( t_group_info->t_group[i].num[1] ) {
                        *s_subtype |= SALT_DONOR_Neg; /* has (-) */
                    }
                    *s_subtype |= SALT_ACCEPTOR; /* there is always an acceptor in a t-group */
                    return type;
                }
            }
            return -1; /* error: t-group not found */
        }

        if ( eif.cAcceptor ) {
            *s_subtype |= SALT_ACCEPTOR;
        }
        if ( eif.cDonor ) {
            if ( at[at_no].charge == -1 ) {
                *s_subtype |= SALT_DONOR_Neg; /* has (-) */
            }
            if ( at[at_no].num_H ) {
                *s_subtype |= SALT_DONOR_H; /* has H */
            }
        }
    }
    return type;
}
/********************************************************************************************************/
int GetOtherSaltType( inp_ATOM *at, int at_no, int *s_subtype )
{
    static int el_number_C  = 0;
   /* static int el_number_N  = 0; */
   /* static int el_number_O  = 0; */
    static int el_number_S  = 0;
    static int el_number_Se = 0;
    static int el_number_Te = 0;

/* 
   type (returned value):
    -1 => ignore
     2 => found:                           SH
          proton donor     -CH2-SH, >CH-SH, >C<    S(-)
          proton acceptor  -CH2-S(-), >CH-S(-), >C<
   subtype:
     1 = SALT_DONOR_H   => has H 
     2 = SALT_DONOR_Neg => has (-) charge
     4 = SALT_ACCEPTOR  => may be an acceptor of H or (-), but not necessarily

   non-O-atom should be:
     - a tautomeric endpoint atom
     - connected to possible middle point atom
*/
    int type, endpoint_valence, bond_type, centerpoint;
    ENDPOINT_INFO eif;

    if ( at[at_no].valence != 1 || at[at_no].chem_bonds_valence != 1 ||
         1 != (at[at_no].num_H==1) + (at[at_no].charge==-1) ) {
        return -1;
    }
    /* one-time initialization */
    if ( !el_number_S ) {
        el_number_C  = get_periodic_table_number( "C" );
       /* el_number_N  = get_periodic_table_number( "N" ); */
       /* el_number_O  = get_periodic_table_number( "O" ); */
        el_number_S  = get_periodic_table_number( "S" );
        el_number_Se = get_periodic_table_number( "Se" );
        el_number_Te = get_periodic_table_number( "Te" );
    }
    *s_subtype = 0; /* initialize the output */
    if ( !(at[at_no].el_number == el_number_S  ||
           at[at_no].el_number == el_number_Se ||
           at[at_no].el_number == el_number_Te ) ) {
        return -1; /* we are not looking for oxygen here */
    }

    type = 2; /* non-tautomeric p-donor or acceptor: C-SH, C-S(-) */

    if ( !(endpoint_valence = nGetEndpointInfo( at, at_no, &eif )) ||
         eif.cMoveableCharge && !at[at_no].c_point || !eif.cDonor || eif.cAcceptor ) {
        return -1; /* not a possible -SH or -S(-) */
    } else {
        /* at[at_no] is not not in a tautomeric group; use eif previously filled out by nGetEndpointInfo */
        /* check whether there is adjacent atom-candidate for a centerpoint */
        centerpoint = (int)at[at_no].neighbor[0];
        bond_type   = (int)at[at_no].bond_type[0] & BOND_TYPE_MASK;
        if ( at[centerpoint].el_number != el_number_C ||
             at[centerpoint].charge ||
             at[centerpoint].radical && at[centerpoint].radical != RADICAL_SINGLET ||
             at[centerpoint].valence != at[centerpoint].chem_bonds_valence ) {
            return -1; /* not a carbon with all single bonds */
        }
        if ( at[at_no].num_H == 1 ) {
            *s_subtype |= SALT_p_DONOR;
        } else
        if ( at[at_no].charge == -1 ) {
            *s_subtype |= SALT_p_ACCEPTOR;
        } else {
            return -1;
        }
    }
    return type;
}

/********************************************************************************************************/
/* new version: merge all, check alt paths, then unmerge unreachable O-atoms if any */
/* Check for oxygen negative charge-H tautomerism (Salts)
   allowed long-range tautomerism; more than one H or (-) can be moved, for example:
   HO-C=C-O(-)         O=C-C=O   
     /   \              /   \    
    R     R            R     R   
    |     |       =>   |     |   
    R'    R'           R'    R'  
     \   /              \   /    
    O=C-C=O           HO-C=C-O(-)

    To check:

                          |             |
     -add all possible HO-C=, O=C, (-)O-C= (including all containing O t-groups) into one t-group;
     -temporarily disconnect one of previously not belonging to any t-group O-atoms from the one t-group;
     -find whether there is an alt path allowing H or (-) to migrate
      from the temp. disconnected O to any one left in the group.
      If the alt path does not exist then the temp. disconnected atom does not
      participate in the H/(-) migrartion and it will be unmarked/unmerged.

*/
/********************************************************************************************************/
int comp_candidates( const void *a1, const void *a2 )
{
    const S_CANDIDATE *s1 = (const S_CANDIDATE *)a1;
    const S_CANDIDATE *s2 = (const S_CANDIDATE *)a2;
    int ret;
    if ( s1->type >= 0 /* enabled < */  && s2->type < 0 /* disabled */ )
        return -1; /* enabled goes first */
    if ( s1->type <  0 /* disabled > */ && s2->type >= 0 /* enabled */ )
        return 1;
    if ( s1->endpoint && !s2->endpoint )
        return -1; /* tautomeric goes first; only tautomeric may be disabled */
    if ( !s1->endpoint && s2->endpoint )
        return 1; /* tautomeric goes first; only tautomeric may be disabled */
    if ( s1->endpoint && s2->endpoint && (ret = (int)s1->endpoint - (int)s2->endpoint) ) {
        return ret;
    }
    return (int)s1->atnumber - (int)s2->atnumber;
}
/********************************************************************************************************/
int MarkSaltChargeGroups2 ( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info,
                          T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD )
{
/* BNS_EDGE_FORBIDDEN_TEMP */
#define ALT_PATH_FOUND    (MAX_ATOMS+1)
#define NO_ENDPOINT       (MAX_ATOMS+2)  /* the two defines must be different */
#define DISABLE_CANDIDATE 10
#define cPAIR(a,b) cPair[a+b*nNumLeftCandidates]
#define ACCEPTOR_PAIR 1
#define DONOR_PAIR    2

    int nNumChanges = 0, nNumOtherChanges = 0, nNumAcidicChanges = 0, nTotNumChanges = 0;
    S_CHAR      *cPair    = NULL;
    T_ENDPOINT  *EndPoint = NULL;
    if ( s_group_info && s_group_info->s_candidate && s_group_info->max_num_candidates > 0 ) {
        int i, j, i1, j1;
        S_CANDIDATE *s_candidate         = s_group_info->s_candidate;
        int          nMaxNumCandidates   = s_group_info->max_num_candidates;
        int          nNumCandidates      = s_group_info->num_candidates;
        int          nNumOtherCandidates = s_group_info->num_other_candidates;
        int          nNumPOnlyCandidates = s_group_info->num_p_only_candidates;
        int          nNumLeftCandidates  = 0;
        int          nNumMarkedCandidates = 0;
        int          s_type, s_subtype;
        int          ret, nDelta;
        int          bHardAddedRemovedProtons = t_group_info && (t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT);

        int          s_subtype_all = 0;
        int          nDonorPairs, nAcceptorPairs, nCurDonorPairs, nCurAcceptorPairs, bAlreadyTested;
/*
        ENDPOINT_INFO    eif;
*/

#if( IGNORE_TGROUP_WITHOUT_H == 1 )
        int          bTGroupHasNegativeChargesOnly = 1;
#endif
        /*return 0;*/ /* debug only */
        
        i1 = -1;

        if ( nNumCandidates <= -2 || !t_group_info || !t_group_info->t_group ) {
            return 0;
        }
        /*************************************************************************/
        /* find all candidates including those with differen s_type (other type) */
        /*************************************************************************/
        for ( i = 0, nNumCandidates = nNumOtherCandidates = nNumPOnlyCandidates = 0; i < num_atoms; i ++ ) {
            if ( 0 == (s_type = GetSaltChargeType( at, i, t_group_info, &s_subtype )) ||
                 /* -C=O or =C-OH, O = S, Se, Te */
                 1 == (s_type = GetOtherSaltChargeType( at, i, t_group_info, &s_subtype, 1/* bAccept_O*/ ))  ||
                 /* =Z-MH or -Z=M, Z = centerpoint, M = endpoint, other than above */
                 2 == (s_type = GetOtherSaltType( at, i, &s_subtype ) ) ||
                 ( bHardAddedRemovedProtons && 4 == (s_type = bIsHardRemHCandidate( at, i, &s_subtype ) ) )
                 /* >C-SH, >C-S(-); S=S,Se,Te */
               ) {

                if ( nNumCandidates >= nMaxNumCandidates ) {
                    return BNS_VERT_EDGE_OVFL;
                }
                s_candidate[nNumCandidates].atnumber = i;
                s_candidate[nNumCandidates].type     = s_type;
                s_candidate[nNumCandidates].subtype  = s_subtype;
                s_candidate[nNumCandidates].endpoint = at[i].endpoint;
                nNumCandidates ++;
                nNumOtherCandidates += (1 == s_type);
                s_subtype_all                        |= s_subtype;
                i1 = i; /* save a representative of a tautomeric group */
            }
        }

        if ( nNumCandidates <= 1 ||  /* TG_FLAG_ALLOW_NO_NEGTV_O <=> CHARGED_SALTS_ONLY=0 */
             !(s_subtype_all & SALT_ACCEPTOR) ||
             (((t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O) ||
               (t_group_info->bTautFlagsDone & TG_FLAG_FOUND_SALT_CHARGES_DONE) ||
               (t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT)) ?
                !(s_subtype_all & (SALT_DONOR)):
                (!(s_subtype_all & SALT_DONOR_Neg) || nNumOtherCandidates == nNumCandidates ))
            ) {
            s_group_info->num_candidates = 0; /* no candidate exists */
            return 0;
        }
        if ( !(s_subtype_all & (SALT_DONOR_Neg) ) ) {
            t_group_info->bTautFlagsDone |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
        }

        /************************************************************************************/
        /* Mark redundant candidates so that only one candidate from one t-group is left in */
        /************************************************************************************/
        for ( i = 0; i < nNumCandidates; i ++ ) {
            if ( 2 == s_candidate[nNumCandidates].type ) {
                s_candidate[i].type -= DISABLE_CANDIDATE; /* disable >C-SH candidates */
                nNumLeftCandidates ++; /* count rejected */
                continue;
            }
            if ( s_candidate[i].endpoint ) {
                for ( j = i-1; 0 <= j; j -- ) {
                    if ( s_candidate[i].endpoint == s_candidate[j].endpoint ) {
                        s_candidate[i].type -= DISABLE_CANDIDATE; /* disable subsequent redundant */
                        nNumLeftCandidates ++; /* count rejected */
                        break;
                    }
                }
            }
        }
        nNumLeftCandidates = nNumCandidates - nNumLeftCandidates; /* subtract num. rejected from the total */
        s_group_info->num_candidates = 0; /* reinit next time */
        /*********************************************************************/
        /* reorder so that all disabled are at the end, tautomeric are first */
        /*********************************************************************/
        qsort ( s_candidate, nNumCandidates, sizeof(s_candidate[0]), comp_candidates );
        cPair = (S_CHAR *)inchi_calloc( nNumLeftCandidates*nNumLeftCandidates, sizeof(cPair[0]) );
        if ( !cPair ) {
            /*printf("BNS_OUT_OF_RAM-6\n");*/
            nTotNumChanges = BNS_OUT_OF_RAM;
            goto quick_exit;
        }
        nDonorPairs = nAcceptorPairs = 0;
        /**********************************************************************/
        /* Find whether we have at least one donor pair and one acceptor pair */
        /**********************************************************************/
        for ( i = 0; i < nNumLeftCandidates; i ++ ) {
            nCurDonorPairs = nCurAcceptorPairs = 0;
            for ( j = 0; j <= i; j ++ ) {
                if ( i == j && !s_candidate[i].endpoint ) {
                    continue;  /* same non-taut atom. However, success for i==j means     *
                                * that the whole tautomeric group may donate or accept 2H */
                }
                /* check for acceptor pair */
                if ( (s_candidate[i].subtype & SALT_ACCEPTOR) && (s_candidate[j].subtype & SALT_ACCEPTOR) &&
                     (ret = bExistsAltPath( pBNS, pBD, NULL, at, num_atoms, s_candidate[i].atnumber,
                                            s_candidate[j].atnumber, ALT_PATH_MODE_ADD2H_TST ))) {
                    if ( IS_BNS_ERROR( ret ) ) {
                        nTotNumChanges = ret;
                        goto quick_exit;
                    }
                    if ( ret & 1 ) {
                        nDelta       = (ret & ~3) >> 2;
                        /*nNumChanges += (ret & 2);*/
                        if ( nDelta ) {
                            /* alt path unleashed previously localized radicals and they annihilated */
                            nNumChanges = 0;
                            nTotNumChanges = BNS_RADICAL_ERR;
                            goto quick_exit;
                        }
                        cPAIR(i,j) |= ACCEPTOR_PAIR; /* the result: mark the pair */
                        /*cPAIR(j,i) |= ACCEPTOR_PAIR;*/
                    }
                }
                /* check for donor pair */
                if ( (s_candidate[i].subtype & SALT_DONOR) && (s_candidate[j].subtype & SALT_DONOR) &&
                     (ret = bExistsAltPath( pBNS, pBD, NULL, at, num_atoms, s_candidate[i].atnumber,
                                            s_candidate[j].atnumber, ALT_PATH_MODE_REM2H_TST ))) {
                    if ( IS_BNS_ERROR( ret ) ) {
                        nTotNumChanges = ret;
                        goto quick_exit;
                    }
                    if ( ret & 1 ) {
                        nDelta       = (ret & ~3) >> 2;
                        /*nNumChanges += (ret & 2);*/
                        if ( nDelta ) {
                            /* alt path unleashed previously localized radicals and they annihilated */
                            nNumChanges = 0;
                            nTotNumChanges = BNS_RADICAL_ERR;
                            goto quick_exit;
                        }
                        cPAIR(i,j) |= DONOR_PAIR; /* the result: mark the pair */
                        /*cPAIR(j,i) |= ACCEPTOR_PAIR;*/
                    }
                }
                /* since the results will be used later to change bonds, check only now */
                /* when both results for (i,j) have been obtained. */
                if ( cPAIR(i,j) & ACCEPTOR_PAIR ) {
                    nCurAcceptorPairs ++;
                    if ( nDonorPairs ) {
                        /* find donor pair (i1,j1) such that i!=i1, i!=j1, j!=i1, j!=j1 */
                        for ( i1 = 0; i1 < i; i1 ++ ) {
                            for ( j1 = 0; j1 <= i1; j1 ++ ) {
                                /* here always j1 < i && i1 < i therefore we do not compare i to i1 or j1 */
                                if ( j1 != j && i1 != j && (cPAIR(i1,j1) & DONOR_PAIR) ) { 
                                    /* both the donor and the acceptor pairs have been found */
                                    goto bFound2Pairs;
                                }
                            }
                        }
                    }
                }
                if ( cPAIR(i,j) & DONOR_PAIR ) {
                    nCurDonorPairs ++;
                    if ( nAcceptorPairs ) {
                        /* find acceptor pair (i1,j1) such that i!=i1, i!=j1, j!=i1, j!=j1 */
                        for ( i1 = 0; i1 < i; i1 ++ ) {
                            for ( j1 = 0; j1 <= i1; j1 ++ ) {
                                /* here always j1 < i && i1 < i therefore we do not compare i to i1 or j1 */
                                if ( j1 != j && i1 != j && (cPAIR(i1,j1) & ACCEPTOR_PAIR) ) { 
                                    /* both the donor and the acceptor pairs have been found */
                                    goto bFound2Pairs;
                                }
                            }
                        }
                    }
                }
            }
            nDonorPairs    += nCurDonorPairs;
            nAcceptorPairs += nCurAcceptorPairs;
        }
        /* nothing has been found */
        nNumChanges = 0;
        inchi_free( cPair );
        cPair = NULL;
        goto quick_exit;


        /* both the donor and the acceptor pairs have been found */
bFound2Pairs:
        /* first, try already found pairs */
        i1 = i;
        j1 = j;

        /* Find all possible donor and acceptor pairs */
        nNumMarkedCandidates = 0;
        for ( i = 0; i < nNumLeftCandidates; i ++ ) {
            nCurDonorPairs = nCurAcceptorPairs = 0;
            for ( j = 0; j <= i; j ++ ) {
                bAlreadyTested = (i < i1 || i == i1 && j <= j1);
                if ( bAlreadyTested && (cPAIR(i,j) & ACCEPTOR_PAIR) || !bAlreadyTested ) {
                    /* checking for acceptor pair */
                    if ( (s_candidate[i].subtype & SALT_ACCEPTOR) && (s_candidate[j].subtype & SALT_ACCEPTOR) &&
                         (ret = bExistsAltPath( pBNS, pBD, NULL, at, num_atoms, s_candidate[i].atnumber,
                                                s_candidate[j].atnumber, ALT_PATH_MODE_ADD2H_CHG ))) {
                        if ( IS_BNS_ERROR( ret ) ) {
                            nTotNumChanges = ret;
                            goto quick_exit;
                        }
                        if ( ret & 1 ) {
                            nDelta       = (ret & ~3) >> 2;
                            nNumChanges += (ret & 2);
                            if ( nDelta ) {
                                /* alt path unleashed previously localized radicals and they annihilated */
                                nNumChanges = 0;
                                nTotNumChanges = BNS_RADICAL_ERR;
                                goto quick_exit;
                            }
                            cPAIR(i,j) |= ACCEPTOR_PAIR;
                            /*cPAIR(j,i) |= ACCEPTOR_PAIR;*/
                            nCurAcceptorPairs += !bAlreadyTested;
                            if ( !(s_candidate[i].subtype & SALT_SELECTED) ) {
                                s_candidate[i].subtype |= SALT_SELECTED;
                                nNumMarkedCandidates ++;
                                if ( !s_candidate[i].endpoint && s_candidate[i].type ) {
                                    nNumOtherChanges ++;
                                } else {
                                    nNumAcidicChanges ++;
                                }
                            }
                            if ( !(s_candidate[j].subtype & SALT_SELECTED) ) {
                                s_candidate[j].subtype |= SALT_SELECTED;
                                nNumMarkedCandidates ++;
                                if ( !s_candidate[j].endpoint && s_candidate[j].type ) {
                                    nNumOtherChanges ++;
                                } else {
                                    nNumAcidicChanges ++;
                                }
                            }
                        }
                    }
                }
                if ( bAlreadyTested && (cPAIR(i,j) & DONOR_PAIR) || !bAlreadyTested ) {
                    /* checking for donor pair */
                    if ( (s_candidate[i].subtype & SALT_DONOR) && (s_candidate[j].subtype & SALT_DONOR) &&
                         (ret = bExistsAltPath( pBNS, pBD, NULL, at, num_atoms, s_candidate[i].atnumber,
                                                s_candidate[j].atnumber, ALT_PATH_MODE_REM2H_CHG ))) {
                        if ( IS_BNS_ERROR( ret ) ) {
                            nTotNumChanges = ret;
                            goto quick_exit;
                        }
                        if ( ret & 1 ) {
                            nDelta       = (ret & ~3) >> 2;
                            nNumChanges += (ret & 2);
                            if ( nDelta ) {
                                /* alt path unleashed previously localized radicals and they annihilated */
                                nNumChanges = 0;
                                nTotNumChanges = BNS_RADICAL_ERR;
                                goto quick_exit;
                            }
                            cPAIR(i,j) |= DONOR_PAIR;
                            /*cPAIR(j,i) |= ACCEPTOR_PAIR;*/
                            nCurDonorPairs += !bAlreadyTested;
                            if ( !(s_candidate[i].subtype & SALT_SELECTED) ) {
                                s_candidate[i].subtype |= SALT_SELECTED;
                                nNumMarkedCandidates ++;
                                if ( !s_candidate[i].endpoint && s_candidate[i].type ) {
                                    nNumOtherChanges ++;
                                } else {
                                    nNumAcidicChanges ++;
                                }
                            }
                            if ( !(s_candidate[j].subtype & SALT_SELECTED) ) {
                                s_candidate[j].subtype |= SALT_SELECTED;
                                nNumMarkedCandidates ++;
                                if ( !s_candidate[j].endpoint && s_candidate[j].type ) {
                                    nNumOtherChanges ++;
                                } else {
                                    nNumAcidicChanges ++;
                                }
                            }
                        }
                    }
                }
            }
            nDonorPairs    += nCurDonorPairs;
            nAcceptorPairs += nCurAcceptorPairs;
        }
        inchi_free( cPair );
        cPair = NULL;

        if ( nNumMarkedCandidates  ) {
            EndPoint = (T_ENDPOINT *)inchi_calloc( nNumMarkedCandidates, sizeof(EndPoint[0]));
            if ( !EndPoint ) {
                /*printf("BNS_OUT_OF_RAM-7\n");*/
                nTotNumChanges = BNS_OUT_OF_RAM;
                goto quick_exit;
            }
            for ( i = 0, j = 0; i < nNumLeftCandidates; i ++ ) {
                if ( s_candidate[i].subtype & SALT_SELECTED ) {
                    s_candidate[i].subtype ^= SALT_SELECTED; /* remove the flag */
                    if ( j < nNumMarkedCandidates ) {
                        i1 = s_candidate[i].atnumber; /* save a representative of the t-group to be created */
                        AddEndPoint( EndPoint+j, at, i1 );
                    }
                    j ++;
                }
            }
            if ( j != nNumMarkedCandidates ) {
                nTotNumChanges = BNS_PROGRAM_ERR;
                goto quick_exit;
            }
            /* merge all marked atoms and their t-groups into one t-group */
            ret = RegisterEndPoints( t_group_info, EndPoint, nNumMarkedCandidates, at, num_atoms, c_group_info, pBNS );
            if (  ret == -1 ) {
                ret = BNS_PROGRAM_ERR;
            }
            if ( ret < 0 ) {
                nTotNumChanges = ret;
                goto quick_exit;
            }
            nTotNumChanges += (ret > 0);
            inchi_free( EndPoint );
            EndPoint = NULL;

            if ( nNumMarkedCandidates ) {
                for ( i = nNumLeftCandidates; i < nNumCandidates; i ++ ) {
                    s_candidate[i].type += DISABLE_CANDIDATE;
                    j1 = s_candidate[i].atnumber;
                    if ( at[j1].endpoint == at[i1].endpoint ) {
                        if ( !s_candidate[i].endpoint && s_candidate[i].type ) {
                            nNumOtherChanges ++;
                        } else {
                            nNumAcidicChanges ++;
                        }
                    }
                }
            } else {
                for ( i = nNumLeftCandidates; i < nNumCandidates; i ++ ) {
                    s_candidate[i].type += DISABLE_CANDIDATE;
                }
            }

            /* find whether the new t-group have any movable H */
            for ( i = 0, bTGroupHasNegativeChargesOnly = 0; i < t_group_info->num_t_groups; i ++ ) {
                if ( t_group_info->t_group[i].nGroupNumber == at[i1].endpoint &&
                     t_group_info->t_group[i].num[0] == t_group_info->t_group[i].num[1] ) {
                    bTGroupHasNegativeChargesOnly = 1;
                    break;
                }
            }
        }
        nTotNumChanges = ( nTotNumChanges > 0);

#if( IGNORE_TGROUP_WITHOUT_H == 1 )
        if ( nTotNumChanges && bTGroupHasNegativeChargesOnly ) {
            nTotNumChanges = 2;  /* means no moveable H has been affected */
        }
#endif
    }
    
quick_exit:
    if ( nNumOtherChanges && nTotNumChanges == 1 ) {
        nTotNumChanges = 5; /* not only acidic atoms merged */
    }
    if ( cPair ) {
        inchi_free( cPair );
        /*cPair = NULL;*/
    }
    if ( EndPoint ) {
        inchi_free ( EndPoint );
        /*EndPoint = NULL;*/
    }
    return nTotNumChanges; /* 0=>no changes, 1=>new salt tautomerism found, 2=>only new charge tautomerism found */
#undef ALT_PATH_FOUND
#undef NO_ENDPOINT
}
/********************************************************************************************************/
/* regular one-path version: find alt paths then merge */
/* Check for oxygen negative charge-H tautomerism (Salts)
   allowed long-range tautomerism; only one H or (-) can be moved, for example:
   HO-C=X-Y=Z-...-C=O  => O=C-X=Y-Z=...=C-OH
*/

#if ( SALT_WITH_PROTONS == 1 )

#define MAX_LOCAL_TGNUM 0 /* was 32; disable since it has not been used */

#if ( MAX_LOCAL_TGNUM > 0 )
typedef struct tagTGroupData {
    S_SHORT nGroupNumber; /* t-group number from t_group_info->t_group->nGroupNumber */
    S_SHORT nGroupIndex;  /* TGroupData[nGroupNumber]nGroupIndex = index of t_group in t_group_info */
    S_SHORT nDonorM;      /* number of endpoint-donors that have negative charge (Minus) */
    S_SHORT nDonorH;      /* number of endpoint-donors that have only H */
    S_SHORT nAccepM;      /* number of endpoint-acceptors that have negative charge (Minus) */
    S_SHORT nAccepH;      /* number of endpoint-acceptors that have H and no negative charge */
    S_SHORT nAccep0;      /* number of endpoint-acceptors that have no H and no negative charge */ 
    S_SHORT nDonorA;      /* number of acidic endpoint-donors */
    S_SHORT nAccepS;      /* number of acidic endpoint-acceptors */
} TGroupData;
#endif
/********************************************************************************************************/
int MarkSaltChargeGroups ( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info,
                          T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD )
{
    
    int nNumChanges = 0, nTotNumChanges = 0;
    if ( s_group_info && s_group_info->s_candidate && s_group_info->max_num_candidates > 0 ) {
        int i, i1, i2, j, j1, j2, jj, ii1, ii2, jj1, jj2, /*k,*/ num_tested;
        S_CANDIDATE *s_candidate         = s_group_info->s_candidate;
        int          nMaxNumCandidates   = s_group_info->max_num_candidates;
        int          nNumCandidates      = s_group_info->num_candidates;
        int          nNumOtherCandidates = s_group_info->num_other_candidates;
        int          nNumPOnlyCandidates = s_group_info->num_p_only_candidates;
        int          s_type, s_subtype;
        int          ret, nDelta, /*nMobile,*/ err = 0;
        int          s_subtype_all = 0;
        int          nGroupNumber;
        T_ENDPOINT   EndPoint[2];
#if ( MAX_LOCAL_TGNUM > 0 )
        TGroupData   tgData[MAX_LOCAL_TGNUM];
        TGroupData   *ptgData = tgData;
#endif
        if ( nNumCandidates <= -1 || !t_group_info || !t_group_info->t_group ) {
            return 0;
        }

        /* count t-groups */
        for ( i = 0, nGroupNumber = 0; i < t_group_info->num_t_groups; i ++ ) {
            if ( nGroupNumber < t_group_info->t_group[i].nGroupNumber ) {
                nGroupNumber = t_group_info->t_group[i].nGroupNumber; /* max. t-group number */
            }
        }
#if ( MAX_LOCAL_TGNUM > 0 )
        /* prepare memory */
        if ( nGroupNumber >= MAX_LOCAL_TGNUM ) {
            if ( !( ptgData = (TGroupData*)inchi_calloc( nGroupNumber+1, sizeof(TGroupData) ) ) ) {
                err = BNS_OUT_OF_RAM;
                goto quick_exit;
            }
        } else {
            memset( ptgData, 0, sizeof(tgData) );
        }
        ptgData[0].nGroupIndex = -1; /* data for non-tautomeric atoms */
        for ( i = 0, nGroupNumber = 0; i < t_group_info->num_t_groups; i ++ ) {
            if ( nGroupNumber = t_group_info->t_group[i].nGroupNumber ) {
                ptgData[nGroupNumber].nGroupIndex = i;
                ptgData[i].nGroupNumber = nGroupNumber;
            }
        }
#endif
        nNumCandidates = 0; /* always recalculate 2004-03-22 */
        num_tested = 0;

        if ( nNumCandidates == 0 ) {
            for ( i = 0, nNumCandidates = nNumOtherCandidates = nNumPOnlyCandidates = 0; i < num_atoms; i ++ ) {
                if ( 0 == (s_type = GetSaltChargeType( at, i, t_group_info, &s_subtype )) ||
                     /* -C=O or =C-OH, O = S, Se, Te */
#if( INCL_NON_SALT_CANDIDATATES == 1 )            
                     1 == (s_type = GetOtherSaltChargeType( at, i, t_group_info, &s_subtype, 1 )) ||
                     /* =Z-MH or -Z=M, Z = centerpoint, M = endpoint, other than above */
#endif
                     2 == (s_type = GetOtherSaltType( at, i, &s_subtype ) ) 
                     /* >C-SH, >C-S(-); S=S,Se,Te */
                   ) {

                    if ( nNumCandidates >= nMaxNumCandidates ) {
                        err = BNS_VERT_EDGE_OVFL;
                        goto quick_exit;
                    }
                    s_candidate[nNumCandidates].atnumber = i;
                    s_candidate[nNumCandidates].type     = s_type;
                    s_candidate[nNumCandidates].subtype  = s_subtype;
                    s_candidate[nNumCandidates].endpoint = at[i].endpoint;
                    nNumCandidates ++;
                    nNumOtherCandidates += (1 == s_type);
                    nNumPOnlyCandidates += (2 == s_type);
                    s_subtype_all                        |= s_subtype;
                    /*i1 = i;*/ /* save a representative of a tautomeric group */
                }
            }

            /* changes: TG_FLAG_ALLOW_NO_NEGTV_O replaced CHARGED_SALTS_ONLY==0 */
            if ( nNumCandidates <= 1 ||
                !(s_subtype_all & SALT_ACCEPTOR) ||
                 (((t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O)||
                   (t_group_info->bTautFlagsDone & TG_FLAG_FOUND_SALT_CHARGES_DONE) ||
                   (t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT)) ?
                    !(s_subtype_all & (SALT_DONOR_Neg | SALT_DONOR_H)):
                    (!(s_subtype_all & SALT_DONOR_Neg) || nNumOtherCandidates==nNumCandidates))
               ) {
                s_group_info->num_candidates = -1; /* no candidate exists */
                goto quick_exit;
            }
            if ( !(s_subtype_all & (SALT_DONOR_Neg) ) ) {
                t_group_info->bTautFlagsDone |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
            }
        } else {
            for ( i = 0; i < nNumCandidates; i ++ ) {
                i1 = s_candidate[i].atnumber;
                if ( 0 <= (s_type = GetSaltChargeType( at, i1, t_group_info, &s_subtype ))
#if( INCL_NON_SALT_CANDIDATATES == 1 )            
                     || 0 < (s_type = GetOtherSaltChargeType( at, i1, t_group_info, &s_subtype, 1 /* bAccept_O*/ ))
#endif
                    ) {
                    s_candidate[nNumCandidates].type     = s_type;
                    s_candidate[nNumCandidates].subtype  = s_subtype;
                    s_candidate[nNumCandidates].endpoint = at[i1].endpoint;
                }
            }
        }
        /* Look for alt paths connecting:
           SALT_DONOR_Neg to SALT_ACCEPTOR  : long distance migration of negative charges
           SALT_DONOR_H   to SALT_ACCEPTOR  : long distance migration of H-atoms
        */
        do {
            nNumChanges = 0;
            for ( i1 = 0; i1 < nNumCandidates; i1 ++ ) {
                j1 = s_candidate[i1].atnumber; 
                for ( i2 = i1+1; i2 < nNumCandidates; i2 ++ ) {
                    /* prev. approach: do not test if both candidates are not "salt-type". Disabled 2004-03-18
                    if ( s_candidate[i1].type && s_candidate[i2].type )
                        continue;
                    */
                    j2 = s_candidate[i2].atnumber;
                    if ( at[j1].endpoint && at[j1].endpoint == at[j2].endpoint ) {
                        continue;
                    }
                    for ( j = 0; j < 2; j ++ ) {
                        if ( j ) {
                            ii1 = i2; /* candidate 1 (donor)    ordering number */
                            ii2 = i1; /* candidate 2 (acceptor) ordering number */
                            jj1 = j2; /* candidate 1 (donor)    atom number */
                            jj2 = j1; /* candidate 2 (acceptor) atom number */
                        } else {      /* transposition */
                            ii1 = i1; /* candidate 1 (donor)    ordering number */
                            ii2 = i2; /* candidate 2 (acceptor) ordering number */
                            jj1 = j1; /* candidate 1 (donor)    atom number     */
                            jj2 = j2; /* candidate 2 (acceptor) atom number     */
                        }

                        if ( ( s_candidate[ii1].subtype & (SALT_DONOR_Neg | SALT_DONOR_H) ) &&
                             ( s_candidate[ii2].subtype & SALT_ACCEPTOR ) ) {
                            ret = bExistsAltPath( pBNS, pBD, NULL, at, num_atoms, jj2, jj1, ALT_PATH_MODE_4_SALT );
                            num_tested ++;
                            if ( IS_BNS_ERROR( ret ) ) {
                                err = ret;
                                goto quick_exit;
                            }
                            if ( ret & 1 ) {
                                nDelta       = (ret & ~3) >> 2;
                                nNumChanges += (ret & 2);
                                for ( i = 0; i < 2; i ++ ) {
                                    jj = i? jj2 : jj1;
                                    AddEndPoint( EndPoint+i, at, jj );
                                }
                                /* add/merge taut groups and reinit pBNS in the fly */
                                ret = RegisterEndPoints(  t_group_info,
                                                          EndPoint, 2, at, num_atoms, c_group_info, pBNS );
                                if (  ret == -1 ) {
                                    ret = BNS_PROGRAM_ERR;
                                }
                                if ( ret < 0 ) {
                                    err = ret;
                                    goto quick_exit;
                                }
                                if ( nDelta ) {
                                    err = BNS_RADICAL_ERR;
                                    goto quick_exit;
                                }
                                nNumChanges += (ret > 0);
                                break; /* avoid redundant repetition */
                            }
                        }
                    }
                }
            }
            nTotNumChanges += nNumChanges;
        } while ( num_tested && nNumChanges );
    
quick_exit:
        if ( !err ) {
            nTotNumChanges += nNumChanges; /* nNumChanges != 0 only in case of 'goto quick_exit' */
            if ( s_group_info->num_candidates == 0 ) {
                /* first time: initialize */
                s_group_info->num_candidates = num_tested? nNumCandidates : -1; /* no candidate exists */
            }
        } else {
            nTotNumChanges = err;
        }
#if ( MAX_LOCAL_TGNUM > 0 )
        if ( ptgData != tgData ) {
            inchi_free( ptgData );
        }
#endif        
    }
    return nTotNumChanges;
}
#else
/********************************************************************************************************/
int MarkSaltChargeGroups ( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info,
                          T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD )
{
    
    int nNumChanges = 0, nTotNumChanges = 0;
    if ( s_group_info && s_group_info->s_candidate && s_group_info->max_num_candidates > 0 ) {
        int i, i1, i2, j, j1, j2, jj, ii1, ii2, jj1, jj2, k, num_tested;
        S_CANDIDATE *s_candidate         = s_group_info->s_candidate;
        int          nMaxNumCandidates   = s_group_info->max_num_candidates;
        int          nNumCandidates      = s_group_info->num_candidates;
        int          nNumOtherCandidates = s_group_info->num_other_candidates;
        int          s_type, s_subtype;
        int          ret, nDelta, nMobile;
        int          s_subtype_all = 0;
        T_ENDPOINT   EndPoint[2];

        if ( nNumCandidates <= -1 || !t_group_info || !t_group_info->t_group ) {
            return 0;
        } else
        if ( nNumCandidates == 0 ) {
            for ( i = 0, nNumCandidates = nNumOtherCandidates = 0; i < num_atoms; i ++ ) {
                if ( 0 <= (s_type = GetSaltChargeType( at, i, t_group_info, &s_subtype )) ) {
                    if ( nNumCandidates >= nMaxNumCandidates ) {
                        return BNS_VERT_EDGE_OVFL;
                    }
                    s_candidate[nNumCandidates].atnumber = i;
                    s_candidate[nNumCandidates].type     = s_type;
                    s_candidate[nNumCandidates].subtype  = s_subtype;
                    s_candidate[nNumCandidates].endpoint = at[i].endpoint;
                    nNumCandidates ++;
                    s_subtype_all                        |= s_subtype;
                    /*i1 = i;*/ /* save a representative of a tautomeric group */
                }
#if( INCL_NON_SALT_CANDIDATATES == 1 )            
                else  /* new */
                if ( 0 < (s_type = GetOtherSaltChargeType( at, i, t_group_info, &s_subtype, 1 /* bAccept_O*/ )) ) {
                    if ( nNumCandidates >= nMaxNumCandidates ) {
                        return BNS_VERT_EDGE_OVFL;
                    }
                    s_candidate[nNumCandidates].atnumber = i;
                    s_candidate[nNumCandidates].type     = s_type;
                    s_candidate[nNumCandidates].subtype  = s_subtype;
                    s_candidate[nNumCandidates].endpoint = at[i].endpoint;
                    nNumCandidates ++;
                    nNumOtherCandidates ++;
                    s_subtype_all                        |= s_subtype;
                }
#endif
            }

            /* changes: TG_FLAG_ALLOW_NO_NEGTV_O replaced CHARGED_SALTS_ONLY==0 */
            if ( nNumCandidates <= 1 || nNumOtherCandidates == nNumCandidates || 
                 ((t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O) ?
                    !(s_subtype_all & (SALT_DONOR_Neg | SALT_DONOR_H)):
                    !(s_subtype_all & SALT_DONOR_Neg)) ||
                 !(s_subtype_all & SALT_ACCEPTOR)) {
                s_group_info->num_candidates = -1; /* no candidate exists */
                return 0;
            }
            if ( !(s_subtype_all & (SALT_DONOR_Neg) ) ) {
                t_group_info->bTautFlagsDone |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
            }
        } else {
            for ( i = 0; i < nNumCandidates; i ++ ) {
                i1 = s_candidate[i].atnumber;
                if ( 0 <= (s_type = GetSaltChargeType( at, i1, t_group_info, &s_subtype ))
#if( INCL_NON_SALT_CANDIDATATES == 1 )            
                     || 0 < (s_type = GetOtherSaltChargeType( at, i1, t_group_info, &s_subtype, 1 /* bAccept_O*/ ))
#endif
                    ) {
                    s_candidate[nNumCandidates].type     = s_type;
                    s_candidate[nNumCandidates].subtype  = s_subtype;
                    s_candidate[nNumCandidates].endpoint = at[i1].endpoint;
                }
            }
        }
        /* Look for alt paths connecting:
           SALT_DONOR_Neg to SALT_ACCEPTOR  : long distance migration of negative charges
           SALT_DONOR_H   to SALT_ACCEPTOR  : long distance migration of H-atoms
        */
        num_tested = 0;
        do {
            nNumChanges = 0;
            for ( i1 = 0; i1 < nNumCandidates; i1 ++ ) {
                j1 = s_candidate[i1].atnumber; 
                for ( i2 = i1+1; i2 < nNumCandidates; i2 ++ ) {
                    if ( s_candidate[i1].type && s_candidate[i2].type )
                        continue; /* both candidates are not "salt-type" */
                    j2 = s_candidate[i2].atnumber;
                    if ( at[j1].endpoint && at[j1].endpoint == at[j2].endpoint ) {
                        continue;
                    }
                    for ( j = 0; j < 2; j ++ ) {
                        if ( j ) {
                            ii1 = i2; /* candidate 1 (donor)    ordering number */
                            ii2 = i1; /* candidate 2 (acceptor) ordering number */
                            jj1 = j2; /* candidate 1 (donor)    atom number */
                            jj2 = j1; /* candidate 2 (acceptor) atom number */
                        } else {      /* transposition */
                            ii1 = i1; /* candidate 1 (donor)    ordering number */
                            ii2 = i2; /* candidate 2 (acceptor) ordering number */
                            jj1 = j1; /* candidate 1 (donor)    atom number     */
                            jj2 = j2; /* candidate 2 (acceptor) atom number     */
                        }

                        if ( ( s_candidate[ii1].subtype & (SALT_DONOR_Neg | SALT_DONOR_H) ) &&
                             ( s_candidate[ii2].subtype & SALT_ACCEPTOR ) ) {
                            ret = bExistsAltPath( pBNS, pBD, NULL, at, num_atoms, jj2, jj1, ALT_PATH_MODE_4_SALT );
                            num_tested ++;
                            if ( IS_BNS_ERROR( ret ) ) {
                                return ret;
                            }
                            if ( ret & 1 ) {
                                nDelta       = (ret & ~3) >> 2;
                                nNumChanges += (ret & 2);
                                for ( i = 0; i < 2; i ++ ) {
                                    jj = i? jj2 : jj1;
                                    EndPoint[i].nAtomNumber  = jj;
                                    EndPoint[i].nEquNumber   = 0;
                                    EndPoint[i].nGroupNumber = at[jj].endpoint;
                                    if ( at[jj].endpoint ) {
                                        memset( EndPoint[i].num, 0, sizeof(EndPoint[i].num) );
                                    } else {
                                        AddAtom2num( EndPoint[i].num, at, jj, 2 ); /* fill out */
                                        AddAtom2DA( EndPoint[i].num_DA, at, jj, 2 );
                                        /*
                                        nMobile  = EndPoint[i].num[1] = (at[jj].charge == -1);
                                        nMobile  = EndPoint[i].num[0] = at[jj].num_H + nMobile;
                                        for ( k = 0; k < T_NUM_ISOTOPIC; k ++ ) {
                                            EndPoint[i].num[T_NUM_NO_ISOTOPIC+k] = at[jj].num_iso_H[NUM_H_ISOTOPES-k-1];
                                        }
                                        */
                                    }
                                }
                                /* add/merge taut groups and reinit pBNS */
                                ret = RegisterEndPoints(  t_group_info,
                                                          EndPoint, 2, at, num_atoms, c_group_info, pBNS );
                                if (  ret < 0 ) {
                                    return ret;
                                }
                                nNumChanges += (ret > 0);
                                if ( nDelta ) {
                                    goto quick_exit;
                                }
                                break; /* avoid redundant repetition */
                            }
                        }
                    }
                }
            }
            nTotNumChanges += nNumChanges;
        } while ( num_tested && nNumChanges );
    
quick_exit:
        nTotNumChanges += nNumChanges; /* nNumChanges != 0 only in case of 'goto quick_exit' */
        if ( s_group_info->num_candidates == 0 ) {
            /* first time: initialize */
            s_group_info->num_candidates = num_tested? nNumCandidates : -1; /* no candidate exists */
        }
        
    }
    return nTotNumChanges;
}
#endif

/*****************************************************************************/
int MergeSaltTautGroups( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info,
                          T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS )
{
    /* count candidates to be connected: exclude pure donors that do not belong to any t-group */
    AT_NUMB    nCurTGroupNumber;
    int        i, j, /*k,*/ ret, iat, /*nMobile,*/ nMinNumEndpoints;
    int        s_subtype_all, s_subtype_taut;
    int        nMaxNumCandidates, nNumCandidates, nNumCandidates2;
    T_ENDPOINT EndPointStackArray[MAX_STACK_ARRAY_LEN]; /* will be reallocated if too short */
    T_ENDPOINT  *EndPoint = EndPointStackArray;


    if ( !s_group_info || !s_group_info->s_candidate || /*s_group_info->num_candidates <= 0 ||*/
         !t_group_info || !t_group_info->t_group || !c_group_info ) {
        return 0;
    }
    nMinNumEndpoints = 0;
    nMaxNumCandidates = s_group_info->max_num_candidates;
    nCurTGroupNumber = MAX_ATOMS;  /* impossible t-group number */
    s_subtype_all = s_subtype_taut = 0;
    /* collect tautomeric acidic O and previously non-tautomeric C-OH, C-SH, C-O(-), C-S(-)  */
    /* find whether previously found tautomeric atoms have both mobile H and (-) */
    if ( 1 || (s_group_info->num_candidates < 0) ) {
        /* can be only -O(-)  and -OH */
        int          s_type, s_subtype;
        S_CANDIDATE *s_candidate       = s_group_info->s_candidate;
        for ( i = 0, nNumCandidates = nNumCandidates2 = 0; i < num_atoms; i ++ ) {
            s_subtype = 0;
            if ( 0 == (s_type = GetSaltChargeType( at, i, t_group_info, &s_subtype )) ||
                 /* -C=O or =C-OH, O = S, Se, Te */
                 
                 /*(t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT) &&*/
                 1 == (s_type = GetOtherSaltChargeType( at, i, t_group_info, &s_subtype, 1/* bAccept_O*/ ))  ||
                 /* =Z-MH or -Z=M, Z = centerpoint, M = endpoint, other than above. M may be N */

                 2 == (s_type = GetOtherSaltType( at, i, &s_subtype )) ||
                 /* >C-SH, >C-S(-); S=S,Se,Te */
                 
                 /* other proton donor or acceptor */
                 bHasAcidicHydrogen( at, i)    && ((s_type=3), (s_subtype = SALT_p_DONOR)) ||
                 bHasAcidicMinus( at, i)       && ((s_type=3), (s_subtype = SALT_p_ACCEPTOR))
               ) {

                if ( nNumCandidates >= nMaxNumCandidates ) {
                    return BNS_VERT_EDGE_OVFL;
                }
                if ( at[i].endpoint ) {
                    s_subtype_taut |= s_subtype;
                } else
                if ( bDoNotMergeNonTautAtom(at, i) ) {
                    continue; /* ignore non-tautomeric N */
                }
                if ( !( s_subtype & SALT_DONOR_ALL ) ||
                     (s_subtype & SALT_ACCEPTOR) && !at[i].endpoint ) {
                    continue;  /* do not include non-taut acceptors like -C=O */
                }
                s_candidate[nNumCandidates].atnumber = i;
                s_candidate[nNumCandidates].type     = s_type;
                s_candidate[nNumCandidates].subtype  = s_subtype;
                s_candidate[nNumCandidates].endpoint = at[i].endpoint;
                nNumCandidates ++;
                s_subtype_all  |= s_subtype;
            }
        }
        /*
         Forced merging occurs upon:
         ===========================
                (t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O) or
                (t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT)


         Allow forced merging in cases:
           {t-groups}  (H, (-)}  {H, (-), t-groups}


         Normal salt merging in cases:
            (H, (-)} {H, (-), t-groups},
           
         Cannot merge H into t-groups if no (-) is present
         */


        if ( (t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O) ||
             (t_group_info->bTautFlagsDone & TG_FLAG_FOUND_SALT_CHARGES_DONE) ||
             (t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT) ) {
            /* force merge even though no negative charges are present */
            if ( nNumCandidates <= 1 ||
                 (!(s_subtype_all & SALT_DONOR_Neg2) || !(s_subtype_all & SALT_DONOR_H2)) &&
                 !t_group_info->num_t_groups ) {
                s_group_info->num_candidates = -1; /* no candidate exists */
                return 0;
            }
        } else {
            /* normal salt mode: merge if both -XH and -X(-) are present */
            if ( nNumCandidates <= 1 ||
                 (!(s_subtype_all & SALT_DONOR_Neg2) || !(s_subtype_all & SALT_DONOR_H2)) ) {
                s_group_info->num_candidates = -1; /* no candidate exists */
                return 0;
            }
        }
        /* -- old code --
        if ( nNumCandidates <= 1 ||
             (((t_group_info->bTautFlags & TG_FLAG_ALLOW_NO_NEGTV_O) ||
               (t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT)) ?
                 !(s_subtype_all & SALT_DONOR_ALL):
                 !(s_subtype_all & SALT_DONOR_Neg2)
             ) 
           ) {
            s_group_info->num_candidates = -1;
            return 0;
        }
        */
        if ( !(s_subtype_all & (SALT_DONOR_Neg2) ) ) {
            t_group_info->bTautFlagsDone |= TG_FLAG_ALLOW_NO_NEGTV_O_DONE;
        }
        s_group_info->num_candidates = nNumCandidates;
    }

    for ( i = 0; i < s_group_info->num_candidates; i ++ ) {
        iat = s_group_info->s_candidate[i].atnumber;
        if ( (s_group_info->s_candidate[i].subtype & SALT_ACCEPTOR) && !at[iat].endpoint ) {
            continue; /* should not happen */
        }
        s_subtype_all |= s_group_info->s_candidate[i].subtype;
        if ( at[iat].endpoint != nCurTGroupNumber || !at[iat].endpoint ) {
            nMinNumEndpoints ++;
        }
        nCurTGroupNumber = (int)at[iat].endpoint;
    }
    if ( nMinNumEndpoints <= 1 ) {
        return 0; /* too few endpoints */
    }
    
    /* make sure we have enough memory */
    if ( nMinNumEndpoints > MAX_STACK_ARRAY_LEN ) {
        if ( !(EndPoint = (T_ENDPOINT *)inchi_calloc( nMinNumEndpoints, sizeof(EndPoint[0]) ) ) ) {
            /*printf("BNS_OUT_OF_RAM-8\n");*/
            return BNS_OUT_OF_RAM;
        }
    }
    
    nCurTGroupNumber = MAX_ATOMS;  /* impossible t-group number */
    for ( i = j = 0; i < s_group_info->num_candidates; i ++ ) {
        iat = s_group_info->s_candidate[i].atnumber;
        if ( s_group_info->s_candidate[i].subtype == SALT_ACCEPTOR && !at[iat].endpoint ) {
            continue;
        }
        if ( at[iat].endpoint != nCurTGroupNumber || !at[iat].endpoint ) {
            AddEndPoint( EndPoint+j, at, iat );
            j ++;
        }
        nCurTGroupNumber = (int)at[iat].endpoint;
    }

    ret = RegisterEndPoints(  t_group_info,
                              EndPoint, j, at, num_atoms, c_group_info, pBNS );
    if (  ret == -1 ) {
        ret = BNS_PROGRAM_ERR;
    }

    if ( EndPoint != EndPointStackArray ) {
        inchi_free( EndPoint );
    }

    return ret;
}

/*****************************************************************************/
int MakeIsotopicHGroup(inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info, 
                       T_GROUP_INFO *t_group_info)
{
    /* all tautomeric atoms and all possible H+ donors and acceptors that have H */
    int        i, j, k, n, bHasH, tg, nError=0;
    int        s_subtype_all, s_subtype_taut;
    int        nMaxNumCandidates, nNumCandidates, nNumNonTautCandidates;


    if ( !s_group_info || !s_group_info->s_candidate || /*s_group_info->num_candidates <= 0 ||*/
         !t_group_info || !t_group_info->t_group ) {
        return 0;
    }
    nMaxNumCandidates = s_group_info->max_num_candidates;
    s_subtype_all = s_subtype_taut = 0;
    memset( t_group_info->num_iso_H, 0, sizeof(t_group_info->num_iso_H) );
    if ( 1 || (s_group_info->num_candidates < 0) ) {
        int          s_type, s_subtype;
        S_CANDIDATE *s_candidate       = s_group_info->s_candidate;
        for ( i = 0, nNumCandidates = nNumNonTautCandidates = 0; i < num_atoms; i ++ ) {
            s_subtype = 0;
            s_type    = 0;
            if ( at[i].endpoint ) {
                if ( (tg = t_group_info->tGroupNumber[at[i].endpoint]) &&
                     at[i].endpoint == t_group_info->t_group[tg-=1].nGroupNumber ) {
                    bHasH = (int)t_group_info->t_group[tg].num[0] - (int)t_group_info->t_group[tg].num[1];
                } else {
                    nError = BNS_PROGRAM_ERR;
                    break;
                }
            } else {
                bHasH = (int)at[i].num_H;
            }
            if ( bHasH && at[i].endpoint || /* tautomeric atoms */

                 /* non-tautomeric heteroatoms that 
                    (a) have H and 
                    (b) may be donors of H
                    therefore may exchange isotopic-non-isotopic H */
                 bHasH &&
                 (0 == (s_type = GetSaltChargeType( at, i, t_group_info, &s_subtype )) ||
                 /* -C=O or =C-OH, O = S, Se, Te */
                 
                 /*(t_group_info->tni.bNormalizationFlags & FLAG_FORCE_SALT_TAUT) &&*/
                 1 == (s_type = GetOtherSaltChargeType( at, i, t_group_info, &s_subtype, 1/* bAccept_O*/ ))  ||
                 /* =Z-MH or -Z=M, Z = centerpoint, M = endpoint, other than above. M may be N */

                 2 == (s_type = GetOtherSaltType( at, i, &s_subtype )) ||
                 /* >C-SH, >C-S(-); S=S,Se,Te */
                 
                 /* other proton donor or acceptor */
                 bHasAcidicHydrogen( at, i) && ((s_type=3), (s_subtype = SALT_p_DONOR)) ||
                 bHasAcidicMinus( at, i)    && ((s_type=3), (s_subtype = SALT_p_ACCEPTOR)) ||
                 bHasOtherExchangableH (at, i) && ((s_type=3), (s_subtype = SALT_DONOR_H)) )
                 
               ) {

                if ( nNumCandidates >= nMaxNumCandidates ) {
                    return BNS_VERT_EDGE_OVFL;
                }
                s_candidate[nNumCandidates].atnumber = i;
                s_candidate[nNumCandidates].type     = s_type;
                s_candidate[nNumCandidates].subtype  = s_subtype;
                s_candidate[nNumCandidates].endpoint = at[i].endpoint;
                nNumCandidates ++;
                nNumNonTautCandidates += !at[i].endpoint;
                s_subtype_all  |= s_subtype;
            }
        }
        if ( nError ) {
            return nError;
        }
        if ( nNumCandidates > 0 ) {
            t_group_info->nIsotopicEndpointAtomNumber = (AT_NUMB *)inchi_calloc( nNumNonTautCandidates+1, sizeof(t_group_info->nIsotopicEndpointAtomNumber[0]));
            t_group_info->nIsotopicEndpointAtomNumber[0] = nNumNonTautCandidates;
            for ( i = 0, n = 1; i < nNumCandidates; i ++ ) {
                k = s_candidate[i].atnumber;
                if ( !at[k].endpoint ) {
                    t_group_info->nIsotopicEndpointAtomNumber[n++] = k;
                }
                for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
                    t_group_info->num_iso_H[j] += at[k].num_iso_H[j];
                }
                at[k].cFlags |= AT_FLAG_ISO_H_POINT;
            }
            t_group_info->nNumIsotopicEndpoints = nNumNonTautCandidates+1;
        }
    }
    return nNumCandidates;
}

/*#else*/ /* } DISCONNECT_SALTS == 0 */

/**********************************************************************************
                       Charges and tautomeric endpoints (N only)
 **********************************************************************************

 H = number of possibly moveable hydrogen atoms
 C = possibly moveable positive charge

 - = single bond
 = = double bond
 # = triple bond

+-----------------------------------------------------------------------------+
|ca-| H    | edges to t-   | 1 bond         |  2 bonds       |  3 bonds    *) |
|se | C    | and c-groups  | (valence)      |  (valence)     |  (valence)     |
| # |      | (edges flow)  |                |                |                |
+---|------+---------------+----------------+----------------+----------------|
| 1 | H=0  | --  (1)       | =NH     (3)    |  =N-     (3)   |   >N-     (3)  |
|   | C=0  | ==            |                |                |                |
+---|------+---------------+----------------+----------------+----------------|
| 2 | H=1  | ==  (2)       | -NH2    (3)    |  -NH-    (3)   |   none         |
|   | C=0  | ==            |                |                |                |
+---|------+---------------+----------------+----------------+----------------|
| 3 | H=0  | --  (0)       | #NH(+)  (4)    |  =N(+)=  (4) +)|   >N(+)=  (4)  |
|   | C=1  | --            | (prohibited    |                |                |
|   |      |               | by edge cap)   |                |                |
+---|------+---------------+----------------+----------------+----------------|
| 4 | H=1  | ==  (1)       | =NH2(+) (4)  +)|  =NH(+)- (4) +)|   >NH(+)- (4)  |     
|   | C=1  | --            |                |                |                |
+---+-------------------------------------------------------------------------+

  *) Cannot be a tautomeric endpoint

  +) The three charged types of atoms [=N(+)=, =NH(+)-, =NH2(+)] should be
     checked for possible H-tautomerism. Other types in the marked by *)
     column should not be checked as long as H(+) exchange is not considered
     tautomeric.

  Other possibilities:  -NH3(+)  >NH2(+)  >N(+)<  cannot be H-tautomeric endpoints.

  Case #1 (H=0, C=0) and #4 (H=1,C=0) is indistinguishable from the
  viewpoint of edges flow and capacities except for flow from N to (+) vertex.

  Without taking precautions H(+) can be transferred

  from =NH2(+)  to =NH,
  from =NH(+)-  to =N-,
  from >NH(+)-  to >N-

  or to any other appropriate atom that has a lone electron pair and bonds
  will not change. In this case no bond must be marked as tautomeric.

  For this reason before attempting to transfer H from one endpoint to
  another the charges on the two atoms should be set to zero by
  forcing zero flow from each of atoms to the (+)-vertices if the
  atoms belong to a c-group.

 **********************************************************************************/


/********************************************************************************************************/
/*   MarkTautomerGroups: do not identify positively charged N as endpoints for now */
int MarkTautomerGroups( inp_ATOM *at, int num_atoms, T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info
                        , struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD )
{
    int i, j, k, m, endpoint_valence, centerpoint, endpoint, bond_type, nMobile, num_changes=0, tot_changes=0;
    T_ENDPOINT EndPoint[MAXVAL];
    T_BONDPOS  BondPos[MAXVAL];
    AT_NUMB    nGroupNumber;
    int        bDiffGroups;
    int  nNumEndPoints, nNumBondPos, nNumPossibleMobile;
    int  bTautBond, bNonTautBond, bAltBond;
    int  nNumDonor, nNumAcceptor, bPossiblyEndpoint;
    T_GROUP *t_group;
    int *pnum_t, max_num_t, bIgnoreIsotopic;
    ENDPOINT_INFO eif1, eif2;
    int nErr = 0;
#define ALLOWED_EDGE(PBNS, IAT,IBOND)  ( !(PBNS) || !(PBNS)->edge || !(PBNS)->vert || !(PBNS)->edge[(PBNS)->vert[IAT].iedge[IBOND]].forbidden)
#define ACTUAL_ORDER(PBNS, IAT,IBOND, BTYPE)  ( ((PBNS) && (PBNS)->edge && (PBNS)->vert &&\
    ((BTYPE)==BOND_ALT_123 || (BTYPE)==BOND_ALT_13 || (BTYPE)==BOND_ALT_23))? (PBNS)->edge[(PBNS)->vert[IAT].iedge[IBOND]].flow+BOND_TYPE_SINGLE:(BTYPE))


    if ( !t_group_info || !(t_group_info->bTautFlags & TG_FLAG_TEST_TAUT__ATOMS) )
        return 0;
    /*  initial t_group allocation */
    if ( !t_group_info->t_group && !t_group_info->max_num_t_groups ) {
        INCHI_MODE bTautFlags     = t_group_info->bTautFlags;       /*  save initial setting */
        INCHI_MODE bTautFlagsDone = t_group_info->bTautFlagsDone;   /*  save previous findings, if any */
        TNI       tni           = t_group_info->tni;
        AT_NUMB   *tGroupNumber = t_group_info->tGroupNumber;
        bIgnoreIsotopic = t_group_info->bIgnoreIsotopic;
        memset( t_group_info, 0, sizeof(*t_group_info) );
        t_group_info->bIgnoreIsotopic = bIgnoreIsotopic; /*  restore initial setting */
        t_group_info->bTautFlags      = bTautFlags;
        t_group_info->bTautFlagsDone  = bTautFlagsDone;
        t_group_info->tni             = tni;
        t_group_info->tGroupNumber    = tGroupNumber;
        t_group_info->max_num_t_groups = num_atoms/2+1; /*  upper limit */
        if (!(t_group_info->t_group = (T_GROUP*)inchi_calloc(t_group_info->max_num_t_groups, sizeof(t_group[0])))) {
            return (t_group_info->max_num_t_groups = -1); /*  failed, out of RAM */
        }
    }
    /*  check if t_group_info exists */
    if ( !t_group_info->t_group || !t_group_info->max_num_t_groups )
        return 0;

    if ( 0 > t_group_info->max_num_t_groups )
        return t_group_info->max_num_t_groups;

    pnum_t          = &t_group_info->num_t_groups; /*  number of found tautomer endpoint groups */
    t_group         =  t_group_info->t_group;
    max_num_t       =  t_group_info->max_num_t_groups;
    bIgnoreIsotopic =  t_group_info->bIgnoreIsotopic;
    /*  1-3 tautomers */
    for ( i = 0; i < num_atoms; i ++ ) {
        /*  find possible endpoint Z = at[i] */
        if ( endpoint_valence = nGetEndpointInfo( at, i, &eif1 ) ) {
            /*  1st endpoint candidate found. Find centerpoint candidate */
            for ( j = 0; j < at[i].valence; j ++ ) {
                bond_type   = (int)at[i].bond_type[j] & ~BOND_MARK_ALL;
#if( FIX_BOND23_IN_TAUT == 1 )
                bond_type = ACTUAL_ORDER(pBNS,i,j,bond_type);
#endif
                centerpoint = (int)at[i].neighbor[j];  /*  a centerpoint candidate */
                if ( (bond_type == BOND_DOUBLE ||
                      bond_type == BOND_ALTERN ||
                      bond_type == BOND_ALT12NS ||
                      bond_type == BOND_TAUTOM) && is_centerpoint_elem( at[centerpoint].el_number ) 
                      && ALLOWED_EDGE(pBNS, i, j)
                      ) {
                    /*  test a centerpoint candidate. */
                    /*  find all endpoints including at[i] and store them into EndPoint[] */
                    nNumPossibleMobile =  0;
                    nGroupNumber       =  (AT_NUMB)num_atoms; /*  greater than any tautomeric group number */
                    bDiffGroups        = -1;         /*  ignore the first difference */
                    nNumDonor = nNumAcceptor = 0;
                    for ( k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k ++ ) {
                        endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
                        bond_type    = (int)at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
#if( FIX_BOND23_IN_TAUT == 1 )
                        bond_type = ACTUAL_ORDER(pBNS,centerpoint,k,bond_type);
#endif
                        bTautBond    =
                        bNonTautBond =
                        bAltBond     =
                        bPossiblyEndpoint = 0;
                        if ( !ALLOWED_EDGE(pBNS, centerpoint, k) ) {
                            continue;
                        } else
                        if ( bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_TAUTOM ) {
                            bTautBond = 1;
#if( REPLACE_ALT_WITH_TAUT == 1 )
                            bAltBond  = (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS);
#endif
                        } else
                        if ( bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE )
                            bNonTautBond = 1;
                        else
                            continue;

                        if ( !(endpoint_valence = nGetEndpointInfo( at, endpoint, &eif1 )) )
                            continue; /*  not an endpoint element or can't have mobile groups */
                        /*  save information about the found possible tautomeric endpoint */
                        /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
                        nMobile  =
                        AddAtom2num( EndPoint[nNumEndPoints].num, at, endpoint, 2 ); /* fill out */
                        AddAtom2DA( EndPoint[nNumEndPoints].num_DA, at, endpoint, 2 );
                        /* --- why is isitopic info missing ? -- see below
                        nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
                        nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
                        */
                        if ( bNonTautBond ) {
                            m = (bond_type == BOND_SINGLE && (nMobile || at[endpoint].endpoint));
                            nNumDonor         += m;
                            bPossiblyEndpoint += m;
                            m = (bond_type == BOND_DOUBLE );
                            nNumAcceptor      += m;
                            bPossiblyEndpoint += m;
                        } else {
                            /*  tautomeric or alternating bond */
                            m = (0 != at[endpoint].endpoint || eif1.cDonor );
                            nNumDonor         += m;
                            bPossiblyEndpoint += m;
                            m = ( at[endpoint].endpoint ||
                                  eif1.cNeutralBondsValence > at[endpoint].valence );
                            nNumAcceptor      += m;
                            bPossiblyEndpoint += m;
                        }
                        if ( !bPossiblyEndpoint )
                            continue;
                        EndPoint[nNumEndPoints].nGroupNumber  = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
                        EndPoint[nNumEndPoints].nEquNumber    = 0;
                        EndPoint[nNumEndPoints].nAtomNumber   = (AT_NUMB)endpoint;
                        if ( nGroupNumber != at[endpoint].endpoint ) {
                            bDiffGroups ++;
                            nGroupNumber  = at[endpoint].endpoint;
                        }
                        
                        /*  save positions of all, not only possibly tautomeric bonds */
#if( REPLACE_ALT_WITH_TAUT != 1 )
                        if ( bNonTautBond || bAltBond ) {
#endif
                            BondPos[nNumBondPos].nAtomNumber    = (AT_NUMB)centerpoint;
                            BondPos[nNumBondPos].neighbor_index = (AT_NUMB)k; /* bond ordering number; used to change bonds to tautomeric only  */
                            nNumBondPos ++;
#if( REPLACE_ALT_WITH_TAUT != 1 )
                        }
#endif
                        /*  mobile group is possible if (a) the endpoint has a mobile group or */
                        /*                              (b) the centerpoint is adjacent to another endpoint */
                        nNumPossibleMobile += (nMobile>0 || at[endpoint].endpoint);
                        nNumEndPoints ++;
                    }
                    if ( nNumEndPoints > 1 && nNumPossibleMobile && nNumDonor && nNumAcceptor ) {
                        /*
                         * a tautomeric group has been found
                         *
                         * at this point:
                         * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
                         * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
                         * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
                         * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
                         */

                        nErr=FindAccessibleEndPoints( EndPoint, &nNumEndPoints, BondPos, &nNumBondPos,
                                                 pBNS, pBD, at, num_atoms, c_group_info, ALT_PATH_MODE_TAUTOM );
                        if ( IS_BNS_ERROR(nErr) ) {
                            return nErr;
                        }
                        nErr = 0;

                        if ( nNumEndPoints > 0 ) {
                            if ( !nGroupNumber || bDiffGroups > 0 ) {
                                num_changes = RegisterEndPoints( t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS );
                                if ( num_changes == -1 ) {
                                    nErr = CT_TAUCOUNT_ERR;
                                }
                                if ( num_changes < 0 ) {
                                    nErr = num_changes;
                                }
                                if ( nErr )
                                    goto exit_function;
                                tot_changes += (num_changes>0);
                            }
                            if ( nNumBondPos > 0 ) {
                                /*  some of the bonds have not been marked as tautomeric yet */
                                num_changes = SetTautomericBonds( at, nNumBondPos, BondPos );
                                tot_changes += (num_changes>0);
                            }
                        }
                    }
                }
            }
        }
    }
#if ( KETO_ENOL_TAUT == 1 )  /***** post v.1 feature *****/
    if ( t_group_info->bTautFlags & TG_FLAG_KETO_ENOL_TAUT ) {
        /* 1,3 keto-enol tautomerism */
        for ( i = 0; i < num_atoms; i ++ ) {
            /*  find possible endpoint Z = at[i] */
            if ( endpoint_valence = nGetEndpointInfo_KET( at, i, &eif1 ) ) {
                /*  1st endpoint candidate found. Find centerpoint candidate */
                for ( j = 0; j < at[i].valence; j ++ ) {
                    bond_type   = (int)at[i].bond_type[j] & ~BOND_MARK_ALL;
#if( FIX_BOND23_IN_TAUT == 1 )
                    bond_type = ACTUAL_ORDER(pBNS,i,j,bond_type);
#endif
                    centerpoint = (int)at[i].neighbor[j];  /*  a centerpoint candidate */
                    if ( (bond_type == BOND_DOUBLE ||
                          bond_type == BOND_ALTERN ||
                          bond_type == BOND_ALT12NS ||
                          bond_type == BOND_TAUTOM) &&
                          is_centerpoint_elem_KET( at[centerpoint].el_number ) &&
                          !at[centerpoint].charge && !at[centerpoint].radical &&
                          /* only normal carbon is allowed */
                          4 == at[centerpoint].chem_bonds_valence + at[centerpoint].num_H
                          && ALLOWED_EDGE(pBNS, i, j)
                          ) {
                        int num_O = 0;
                        int num_C = 0;
                        /*  test a centerpoint candidate. */
                        /*  find all endpoints including at[i] and store them into EndPoint[] */
                        nNumPossibleMobile =  0;
                        nGroupNumber       =  (AT_NUMB)num_atoms; /*  greater than any tautomeric group number */
                        bDiffGroups        = -1;         /*  ignore the first difference */
                        nNumDonor = nNumAcceptor = 0;
                        for ( k = 0, nNumEndPoints = 0, nNumBondPos = 0; k < at[centerpoint].valence; k ++ ) {
                            endpoint = at[centerpoint].neighbor[k]; /*  endpoint candidate */
                            bond_type    = (int)at[centerpoint].bond_type[k] & ~BOND_MARK_ALL;
#if( FIX_BOND23_IN_TAUT == 1 )
                            bond_type = ACTUAL_ORDER(pBNS,centerpoint,k,bond_type);
#endif
                            bTautBond    =
                            bNonTautBond =
                            bAltBond     =
                            bPossiblyEndpoint = 0;
                            if ( !ALLOWED_EDGE(pBNS, centerpoint, k) ) {
                                continue;
                            } else
                            if ( bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS || bond_type == BOND_TAUTOM ) {
                                bTautBond = 1;
#if( REPLACE_ALT_WITH_TAUT == 1 )
                                bAltBond  = (bond_type == BOND_ALTERN || bond_type == BOND_ALT12NS);
#endif
                            } else
                            if ( bond_type == BOND_SINGLE || bond_type == BOND_DOUBLE )
                                bNonTautBond = 1;
                            else
                                continue;

                            if ( !(endpoint_valence = nGetEndpointInfo_KET( at, endpoint, &eif2 )) ) {
                                continue;
                            }
                            /*
                            if ( 3 != eif1.cKetoEnolCode + eif2.cKetoEnolCode && endpoint != i )
                                continue;
                            */
                            /*  save information about the found possible tautomeric endpoint */
                            /*  2 = T_NUM_NO_ISOTOPIC non-isotopic values */
                            nMobile  =
                            AddAtom2num( EndPoint[nNumEndPoints].num, at, endpoint, 2 ); /* fill out */
                            AddAtom2DA( EndPoint[nNumEndPoints].num_DA, at, endpoint, 2 );
                            /* --- why is isitopic info missing ? -- see below
                            nMobile  = EndPoint[nNumEndPoints].num[1] = (at[endpoint].charge == -1);
                            nMobile  = EndPoint[nNumEndPoints].num[0] = at[endpoint].num_H + nMobile;
                            */
                            if ( bNonTautBond ) {
                                m = (bond_type == BOND_SINGLE && (nMobile || at[endpoint].endpoint));
                                nNumDonor         += m;
                                bPossiblyEndpoint += m;
                                m = (bond_type == BOND_DOUBLE );
                                nNumAcceptor      += m;
                                bPossiblyEndpoint += m;
                            } else {
                                /*  tautomeric or alternating bond */
                                m = (0 != at[endpoint].endpoint || eif1.cDonor );
                                nNumDonor         += m;
                                bPossiblyEndpoint += m;
                                m = ( at[endpoint].endpoint ||
                                      eif1.cNeutralBondsValence > at[endpoint].valence );
                                nNumAcceptor      += m;
                                bPossiblyEndpoint += m;
                            }
                            if ( !bPossiblyEndpoint )
                                continue;

                            num_O += (endpoint_valence == 2);
                            num_C += (endpoint_valence == 4);

                            EndPoint[nNumEndPoints].nGroupNumber  = at[endpoint].endpoint; /* =0 if it is an endpoint for the 1st time */
                            EndPoint[nNumEndPoints].nEquNumber    = 0;
                            EndPoint[nNumEndPoints].nAtomNumber   = (AT_NUMB)endpoint;
                            if ( nGroupNumber != at[endpoint].endpoint ) {
                                bDiffGroups ++;
                                nGroupNumber  = at[endpoint].endpoint;
                            }
                            
                            /*  save positions of all, not only possibly tautomeric bonds */
#if( REPLACE_ALT_WITH_TAUT != 1 )
                            if ( bNonTautBond || bAltBond ) {
#endif
                                BondPos[nNumBondPos].nAtomNumber    = (AT_NUMB)centerpoint;
                                BondPos[nNumBondPos].neighbor_index = (AT_NUMB)k; /* bond ordering number; used to change bonds to tautomeric only  */
                                nNumBondPos ++;
#if( REPLACE_ALT_WITH_TAUT != 1 )
                            }
#endif
                            /*  mobile group is possible if (a) the endpoint has a mobile group or */
                            /*                              (b) the centerpoint is adjacent to another endpoint */
                            nNumPossibleMobile += (nMobile>0 || at[endpoint].endpoint);
                            nNumEndPoints ++;
                        }
                        if ( nNumEndPoints > 1 && nNumPossibleMobile && nNumDonor && nNumAcceptor && num_O==1 && num_C ) {
                            /*
                             * a tautomeric group has been found
                             *
                             * at this point:
                             * nGroupNumber = 0 if all endpoints belong to a newly discovered tautomeric group
                             * bDiffGroups  > 0 if at least 2 tautomeric groups are to be merged (one of them can be new)
                             * case (nGroupNumber != 0 && bDiffGroups = 0 ) ignored because all endpoints belong to the same known t-group
                             * case (nGroupNumber != 0 && bDiffGroups < 0 ) cannot happen
                             */

                            nErr=FindAccessibleEndPoints( EndPoint, &nNumEndPoints, BondPos, &nNumBondPos,
                                                     pBNS, pBD, at, num_atoms, c_group_info, ALT_PATH_MODE_TAUTOM_KET );
                            if ( IS_BNS_ERROR(nErr) ) {
                                return nErr;
                            }
                            nErr = 0;

                            if ( nNumEndPoints > 0 ) {
                                if ( !nGroupNumber || bDiffGroups > 0 ) {
                                    num_changes = RegisterEndPoints( t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS );
                                    if ( num_changes == -1 ) {
                                        nErr = CT_TAUCOUNT_ERR;
                                    }
                                    if ( num_changes < 0 ) {
                                        nErr = num_changes;
                                    }
                                    if ( nErr )
                                        goto exit_function;
                                    tot_changes += (num_changes>0);
                                }
                                if ( nNumBondPos > 0 ) {
                                    /*  some of the bonds have not been marked as tautomeric yet */
                                    num_changes = SetTautomericBonds( at, nNumBondPos, BondPos );
                                    tot_changes += (num_changes>0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif  /* KETO_ENOL_TAUT */

#if( TAUT_OTHER == 1 ) /* { */
    if ( !tot_changes ) {
#define MAX_ALT_PATH_LEN 8
        int nMaxLenDfsPath = MAX_ALT_PATH_LEN;
        int i1, i2;
        AT_RANK *nDfsPathPos = (AT_RANK  *)inchi_calloc( num_atoms, sizeof(nDfsPathPos[0]) );
        DFS_PATH DfsPath[MAX_ALT_PATH_LEN];
        int      ret;
        if ( !nDfsPathPos || !DfsPath ) {
            tot_changes = CT_OUT_OF_RAM;  /*   <BRKPT> */
            goto free_memory;
        }
#if( TAUT_15_NON_RING      == 1 ) /***** post v.1 feature *****/
        if ( t_group_info->bTautFlags & TG_FLAG_1_5_TAUT ) {
            /*  1,5 tautomerism; one of the endpoints should no be on a ring  */
            /*
                  O                OH                 O      
                  ||               |                  ||     
                  A--pos-          A--pos-            A--pos-
                 /   sib-        //   sib-     ?     /   sib-
                C    ly          C    ly            CH   ly
                \\   a     <-->   \   a       <-->   \   a   
                  B--ring          B--ring            B--ring
                  |                ||                 ||     
                  NH               N                  N

               Note: few recent modifications now allow the terminal N be in a ring, too
             */
            for ( i1 = 0; i1 < num_atoms; i1 ++ ) {
                /*  find possible endpoint Z = at[i1] */
                if ( !(endpoint_valence = nGetEndpointInfo( at, i1, &eif1 ) ) /*||
                     at[i1].nNumAtInRingSystem > 1*/ ) {
                    continue; /*  not a possibly endpoint */
                }

                if ( 1 ) {
                    nNumEndPoints = 0;
                    nNumBondPos   = 0;

                    ret = nGet15TautInAltPath( at, i1, nDfsPathPos,
                                                  DfsPath,  nMaxLenDfsPath,
                                                  EndPoint, sizeof(EndPoint)/sizeof(EndPoint[0]),
                                                  BondPos, sizeof(BondPos)/sizeof(BondPos[0]),
                                                  &nNumEndPoints, &nNumBondPos, 
                                                  pBNS, pBD, num_atoms);
                    if ( ret > 0 ) {
                        if ( nNumEndPoints ) {
                            num_changes = RegisterEndPoints( t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
                            if ( num_changes == -1 ) {
                                nErr = CT_TAUCOUNT_ERR;
                            }
                            if ( num_changes < 0 ) {
                                nErr = num_changes;
                            }
                            if ( nErr )
                                goto free_memory;
                            tot_changes += (num_changes > 0);
                        }
                        if ( nNumBondPos ) {
                            tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
                        }
                    } else
                    if ( IS_BNS_ERROR( ret ) ) {
                        nErr = ret;
                        goto free_memory;
                    }
                }
            }
        }
#endif
#if( TAUT_4PYRIDINOL_RINGS == 1 )
        /*  6-member rings */
        /*
              O              OH             OH   
              ||             |              |    
             /  \          //  \           /  \\ 
            ||   ||  <-->  |    ||  <-->  ||   | 
             \  /          \\  /           \  // 
              NH             N              N    
         */
        for ( i1 = 0; i1 < num_atoms; i1 ++ ) {
            /*  find possible endpoint Z = at[i1] */
            if ( 3 != (endpoint_valence = nGetEndpointInfo( at, i1, &eif1 ) ) ||
                 2 != at[i1].valence ) {
                continue; /*  not a nitrogen atom or a wrong valence */
            }

            if ( at[i1].nNumAtInRingSystem >= 6 ) {
                nNumEndPoints = 0;
                nNumBondPos   = 0;

                ret = nGet15TautIn6MembAltRing( at, i1, nDfsPathPos,
                                              DfsPath,  nMaxLenDfsPath,
                                              EndPoint, sizeof(EndPoint)/sizeof(EndPoint[0]),
                                              BondPos, sizeof(BondPos)/sizeof(BondPos[0]),
                                              &nNumEndPoints, &nNumBondPos, 
                                              pBNS, pBD, num_atoms);
                if ( ret > 0 ) {
                    if ( nNumEndPoints ) {
                        num_changes = RegisterEndPoints( t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
                        if ( num_changes == -1 ) {
                            nErr = CT_TAUCOUNT_ERR;
                        }
                        if ( num_changes < 0 ) {
                            nErr = num_changes;
                        }
                        if ( nErr )
                            goto free_memory;
                        tot_changes += (num_changes > 0);
                    }
                    if ( nNumBondPos ) {
                        tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
                    }
                } else
                if ( IS_BNS_ERROR( ret ) ) {
                    nErr = ret;
                    goto free_memory;
                }
            }
        }
#endif /* TAUT_4PYRIDINOL_RINGS */
#if( TAUT_PYRAZOLE_RINGS == 1 )
        /* 5-member rings:

            Z               Z   
          /  \\           //  \ 
         X     Y  <-->   X     Y
         \\   /          \    //
          N--NH           HN--N  
                       
             ^             ^
             search for these NH
        */
        /*  5-member rings (pyrazole derivatives): look for the neighboring N */
        for ( i1 = 0; i1 < num_atoms; i1 ++ ) {
            if ( 2 == at[i1].valence &&
                 at[i1].nNumAtInRingSystem >= 5 &&
                 3 == (endpoint_valence = nGetEndpointInfo( at, i1, &eif1 ))
               ) {
                nMobile = at[i1].num_H + (at[i1].charge == -1);
                for ( j = 0; j < at[i1].valence; j ++ ) {
                    int nMobile2, endpoint_valence2;
                    i2 = at[i1].neighbor[j];

                    /*  may be important */
                    if ( i2 >= i1 )
                        continue; /*  do not try same pair 2 times */

                    if ( at[i2].nRingSystem != at[i1].nRingSystem )
                        continue;
                    
                    bond_type = (at[i1].bond_type[j] & ~BOND_MARK_ALL);
                    if ( bond_type != BOND_SINGLE &&
                         bond_type != BOND_TAUTOM &&
                         bond_type != BOND_ALT12NS &&
                         bond_type != BOND_ALTERN ||  /* added 1-15-2002 */
                         2 != at[i2].valence      ||
                         3 != (endpoint_valence2 = nGetEndpointInfo( at, i2, &eif2 ) ) ) {
                        continue; /*  not a nitrogen atom or a wrong valence or not a single bond */
                    }
                    nMobile2 = at[i2].num_H + (at[i2].charge == -1);  /*  number of mobile groups */
#if( TAUT_IGNORE_EQL_ENDPOINTS == 1 )
                      if ( at[i1].endpoint && at[i1].endpoint == at[i2].endpoint )
                          continue; /* atoms already belong to the same t-group */
#endif
                    if ( !at[i1].endpoint && !at[i2].endpoint && 1!=nMobile + nMobile2 )
                        continue;

                    ret = nGet12TautIn5MembAltRing( at, i1, j, nDfsPathPos,
                                                  DfsPath,  nMaxLenDfsPath,
                                                  EndPoint, sizeof(EndPoint)/sizeof(EndPoint[0]),
                                                  BondPos, sizeof(BondPos)/sizeof(BondPos[0]),
                                                  &nNumEndPoints, &nNumBondPos 
                                                 , pBNS, pBD, num_atoms);
                    if ( ret > 0 ) {
                        if ( nNumEndPoints ) {
                            num_changes = RegisterEndPoints( t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
                            if ( num_changes == -1 ) {
                                nErr = CT_TAUCOUNT_ERR;
                            }
                            if ( num_changes < 0 ) {
                                nErr = num_changes;
                            }
                            if ( nErr )
                                goto free_memory;
                            tot_changes += (num_changes > 0);
                        }
                        if ( nNumBondPos ) {
                            tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
                        }
                    } else
                    if ( IS_BNS_ERROR( ret ) ) {
                        nErr = ret;
                        goto free_memory;
                    }
                }
            }
        }
#endif /* TAUT_PYRAZOLE_RINGS */
#if ( TAUT_TROPOLONE_7 == 1 || TAUT_TROPOLONE_5 == 1 ) /* { */
        /********************************************************
         *                                         A  B
         *                                         | ||
         * 7-member rings (tropolones): look for M=Q--R--ZH,
         *                                       ^ ^  ^  ^
         *                               endpoint1 i1 i2 endpoint2
         * where A-Q-R=B belong to a 7-member alt. (except Q-R bond) ring: ..=A-(Q-R)=B-..
         * Bond Q-R should be single or tautomeric or alternating
         * M=Q and R-ZH should be chain (non-ring) bonds
         * Same for 5-member rings
         */
        for ( i1 = 0; i1 < num_atoms; i1 ++ ) {
            if ( at[i1].nNumAtInRingSystem >=
#if( TAUT_TROPOLONE_5 == 1 )
                  5
#else
                  7 
#endif
                 &&
                 bIsCenterPointStrict( at, i1 ) &&
#if( TAUT_RINGS_ATTACH_CHAIN == 1 )        
                 at[i1].bCutVertex &&
#endif
                 at[i1].valence == 3 && !at[i1].endpoint ) {
                int nMobile1, endpoint1, endpoint1_valence, bond_type1;
                int nMobile2, endpoint2, endpoint2_valence, bond_type2;
                for ( j = 0; j < at[i1].valence; j ++ ) {
                    i2 = at[i1].neighbor[j];
                    /*
                         // may be important
                    if ( i2 > i1 )
                        continue; // do not try same pair 2 times
                    */
                    if ( at[i2].nRingSystem != at[i1].nRingSystem ||
                         !bIsCenterPointStrict( at, i2 ) ||
#if( TAUT_RINGS_ATTACH_CHAIN == 1 )        
                         !at[i2].bCutVertex ||
#endif                         
                         at[i2].valence != 3 || at[i2].endpoint )
                        continue;
                    bond_type = (at[i1].bond_type[j] & ~BOND_MARK_ALL);
                    if ( bond_type != BOND_SINGLE &&
                         bond_type != BOND_TAUTOM &&
                         bond_type != BOND_ALT12NS &&
                         bond_type != BOND_ALTERN ) {
                        continue; /*  not a single bond between Q-R */
                    }
                    /*  find endpoints */
                    for ( k = 0; k < at[i1].valence; k ++ ) {
                        endpoint1 = at[i1].neighbor[k];
                        if ( endpoint1 == i2 )
                            continue; /*  j == k */
                        if ( !(endpoint1_valence = nGetEndpointInfo( at, endpoint1, &eif1 ) ) )
                            continue; /*  not an endpoint1 element or can't have mobile groups */
#if( TAUT_RINGS_ATTACH_CHAIN == 1 )        
                        if ( at[endpoint1].nRingSystem == at[i1].nRingSystem )
                            continue;
#endif
                        nMobile1  = at[endpoint1].num_H + (at[endpoint1].charge == -1);  /*  number of mobile groups */
                        if ( nMobile1 + at[endpoint1].chem_bonds_valence != endpoint1_valence )
                            continue; /*  abnormal endpoint1 valence; ignore. */
                        bond_type1 = (at[i1].bond_type[k] & ~BOND_MARK_ALL);
                        
                        if ( bond_type1 != BOND_SINGLE &&
                             bond_type1 != BOND_DOUBLE &&
                             bond_type1 != BOND_TAUTOM &&
                             bond_type1 != BOND_ALT12NS &&
                             bond_type1 != BOND_ALTERN )
                            continue;
                        
                        for ( m = 0; m < at[i2].valence; m ++ ) {
                            endpoint2 = at[i2].neighbor[m];
                            if ( endpoint2 == i1 )
                                continue;
                            if ( !(endpoint2_valence = nGetEndpointInfo( at, endpoint2, &eif2 )) )
                                continue; /*  not an endpoint2 element or can't have mobile groups */
#if( TAUT_RINGS_ATTACH_CHAIN == 1 )        
                            if ( at[endpoint2].nRingSystem == at[i2].nRingSystem )
                                continue;
#endif
                            nMobile2  = at[endpoint2].num_H + (at[endpoint2].charge == -1);  /*  number of mobile groups */
                            bond_type2 = (at[i2].bond_type[m] & ~BOND_MARK_ALL);
                            
                            if ( bond_type2 != BOND_SINGLE &&
                                 bond_type2 != BOND_DOUBLE &&
                                 bond_type2 != BOND_TAUTOM &&
                                 bond_type2 != BOND_ALT12NS &&
                                 bond_type2 != BOND_ALTERN )
                                continue;
                            
                            /*  final test for possible tautomerism */
                            nMobile = 0;
                            
                            if ( ALLOWED_EDGE(pBNS, i1, k) && ALLOWED_EDGE(pBNS, i2, m) ) {
                            
                                /*  can mobile group move from 1 to 2? */
                                nMobile += (at[endpoint1].endpoint || nMobile1) &&  /*  from endpoint1 */
                                           (bond_type1 != BOND_DOUBLE)   &&
                                       
                                            (at[endpoint2].endpoint ||          /*  to endpoint2 */
                                            eif2.cNeutralBondsValence > at[endpoint2].valence ) &&
                                           (bond_type2 != BOND_SINGLE); 


                                /*  can mobile group move from 2 to 1? */
                                nMobile += (at[endpoint2].endpoint || nMobile2) &&  /*  from endpoint2 */
                                           (bond_type2 != BOND_DOUBLE)   && /*changed from BOND_SINGLE 2004-02-26 */
                                       
                                            (at[endpoint1].endpoint ||          /*  to endpoint1 */
                                            eif1.cNeutralBondsValence > at[endpoint1].valence ) &&
                                           (bond_type1 != BOND_SINGLE);
                            }
                            if ( !nMobile )
                                continue;
                            
                            if ( bond_type1 == bond_type2 &&
                                 (bond_type1 == BOND_SINGLE || bond_type1 == BOND_DOUBLE) )
                                continue;
                            /* -- old --
                            if ( !at[endpoint1].endpoint && !at[endpoint2].endpoint && 1 != nMobile1 + nMobile2 )
                                continue;
                            */
                            /* -- new --

                            if ( !at[endpoint1].endpoint && !at[endpoint2].endpoint ) {
                                if ( !(bond_type1 == BOND_SINGLE || bond_type1 == BOND_DOUBLE) ||
                                     !(bond_type2 == BOND_SINGLE || bond_type2 == BOND_DOUBLE) ) {
                                    // at this point bond_type1 != bond_type2
                                    continue;
                                }
                                if ( bond_type1 == BOND_SINGLE && !nMobile1 ||
                                     bond_type2 == BOND_SINGLE && !nMobile2 ||
                                     0 == nMobile1 + nMobile2 ) {
                                    continue;
                                }
                            }
                            */
#if ( TAUT_TROPOLONE_7 == 1 )
                            if ( at[i1].nNumAtInRingSystem >= 7 ) {
                                ret = nGet14TautIn7MembAltRing( at, i1, j, k, m, nDfsPathPos,
                                                              DfsPath,  nMaxLenDfsPath,
                                                              EndPoint, sizeof(EndPoint)/sizeof(EndPoint[0]),
                                                              BondPos, sizeof(BondPos)/sizeof(BondPos[0]),
                                                              &nNumEndPoints, &nNumBondPos, 
                                                              pBNS, pBD, num_atoms);
                                if ( ret > 0 ) {
                                    if ( nNumEndPoints ) {
                                        num_changes = RegisterEndPoints( t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
                                        if ( num_changes == -1 ) {
                                            nErr = CT_TAUCOUNT_ERR;
                                        }
                                        if ( num_changes < 0 ) {
                                            nErr = num_changes;
                                        }
                                        if ( nErr )
                                            goto free_memory;
                                        tot_changes += (num_changes > 0);
                                    }
                                    if ( nNumBondPos ) {
                                        tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
                                    }
                                } else
                                if ( IS_BNS_ERROR( ret ) ) {
                                    nErr = ret;
                                    goto free_memory;
                                }
                            }
#endif

#if ( TAUT_TROPOLONE_5 == 1 )
                            if ( at[i1].nNumAtInRingSystem >= 5 ) {
                                ret = nGet14TautIn5MembAltRing( at, i1, j, k, m, nDfsPathPos,
                                                              DfsPath,  nMaxLenDfsPath,
                                                              EndPoint, sizeof(EndPoint)/sizeof(EndPoint[0]),
                                                              BondPos, sizeof(BondPos)/sizeof(BondPos[0]),
                                                              &nNumEndPoints, &nNumBondPos, 
                                                              pBNS, pBD, num_atoms);
                                if ( ret > 0 ) {
                                    if ( nNumEndPoints ) {
                                        num_changes = RegisterEndPoints( t_group_info, EndPoint, nNumEndPoints, at, num_atoms, c_group_info, pBNS);
                                        if ( num_changes == -1 ) {
                                            nErr = CT_TAUCOUNT_ERR;
                                        }
                                        if ( num_changes < 0 ) {
                                            nErr = num_changes;
                                        }
                                        if ( nErr )
                                            goto free_memory;
                                        tot_changes += (num_changes > 0);
                                    }
                                    if ( nNumBondPos ) {
                                        tot_changes += ( 0 < SetTautomericBonds( at, nNumBondPos, BondPos ) );
                                    }
                                } else
                                if ( IS_BNS_ERROR( ret ) ) {
                                    nErr = ret;
                                    goto free_memory;
                                }
                            }
#endif
                        }
                    }
                }
            }
        }
#endif /* } TAUT_TROPOLONE */
free_memory:
        if ( nDfsPathPos ) {
            inchi_free( nDfsPathPos );
        }
#undef MAX_ALT_PATH_LEN
    }
#endif  /* } FIND_RING_SYSTEMS */
exit_function:
    return nErr < 0? nErr : tot_changes;
}

/******************************************************************************/
int free_t_group_info( T_GROUP_INFO *t_group_info )
{
    if ( t_group_info ) {
        if ( t_group_info->t_group ) {
            inchi_free( t_group_info->t_group );
        }
        if ( t_group_info->nEndpointAtomNumber ) {
            inchi_free( t_group_info->nEndpointAtomNumber );
        }
        if ( t_group_info->tGroupNumber ) {
            inchi_free( t_group_info->tGroupNumber );
        }
        if ( t_group_info->nIsotopicEndpointAtomNumber ) {
            inchi_free( t_group_info->nIsotopicEndpointAtomNumber );
        }
        memset( t_group_info, 0, sizeof(*t_group_info));
    }
    return 0;
}

/*******************************************************************************/
/**/
int make_a_copy_of_t_group_info( T_GROUP_INFO *t_group_info, T_GROUP_INFO *t_group_info_orig )
{
    int err = 0, len;
    free_t_group_info( t_group_info );
    if ( t_group_info_orig && t_group_info ) {
        if ( (len=t_group_info_orig->max_num_t_groups) > 0 ) {
            if (t_group_info->t_group =
                (T_GROUP*)inchi_malloc( len * sizeof(t_group_info->t_group[0]))) {
                memcpy(t_group_info->t_group,
                       t_group_info_orig->t_group,
                       len * sizeof(t_group_info->t_group[0]));
            } else {
                err ++;
            }
        }
        if ( (len = t_group_info_orig->nNumEndpoints) > 0 ) {
            if (t_group_info->nEndpointAtomNumber =
                (AT_NUMB*)inchi_malloc( len * sizeof(t_group_info->nEndpointAtomNumber[0]))) {
                memcpy(t_group_info->nEndpointAtomNumber,
                       t_group_info_orig->nEndpointAtomNumber,
                       len * sizeof(t_group_info->nEndpointAtomNumber[0]));
            } else {
                err ++;
            }
        }
        if ( (len = t_group_info_orig->num_t_groups) > 0 ) {
            if (t_group_info->tGroupNumber =
                (AT_NUMB*)inchi_malloc( len * TGSO_TOTAL_LEN * sizeof(t_group_info->tGroupNumber[0]))) {
                memcpy(t_group_info->tGroupNumber,
                       t_group_info_orig->tGroupNumber,
                       len * TGSO_TOTAL_LEN * sizeof(t_group_info->tGroupNumber[0]));
            } else {
                err ++;
            }
        }
        if ( (len = t_group_info_orig->nNumIsotopicEndpoints) > 0 ) {
            if (t_group_info->nIsotopicEndpointAtomNumber =
                (AT_NUMB*)inchi_malloc( len * sizeof(t_group_info->nIsotopicEndpointAtomNumber[0]))) {
                memcpy(t_group_info->nIsotopicEndpointAtomNumber,
                       t_group_info_orig->nIsotopicEndpointAtomNumber,
                       len * sizeof(t_group_info->nIsotopicEndpointAtomNumber[0]));
            } else {
                err ++;
            }
        }
        if ( !err ) {
            t_group_info->nNumEndpoints              = t_group_info_orig->nNumEndpoints;   
            t_group_info->num_t_groups               = t_group_info_orig->num_t_groups;    
            t_group_info->max_num_t_groups           = t_group_info_orig->max_num_t_groups;
            t_group_info->bIgnoreIsotopic            = t_group_info_orig->bIgnoreIsotopic;
            t_group_info->nNumIsotopicEndpoints      = t_group_info_orig->nNumIsotopicEndpoints;
            t_group_info->tni                        = t_group_info_orig->tni;
            /*
            t_group_info->nNumRemovedExplicitH       = t_group_info_orig->nNumRemovedExplicitH;
            t_group_info->nNumRemovedProtons         = t_group_info_orig->nNumRemovedProtons;
            t_group_info->bNormalizationFlags        = t_group_info_orig->bNormalizationFlags;
            */
            /*
            t_group_info->bHardAddedRemovedProtons   = t_group_info_orig->bHardAddedRemovedProtons;
            t_group_info->bSimpleAddedRemovedProtons = t_group_info_orig->bSimpleAddedRemovedProtons;
            t_group_info->nNumCanceledCharges        = t_group_info_orig->nNumCanceledCharges;
            */
        }
        t_group_info->bTautFlags         = t_group_info_orig->bTautFlags;
        t_group_info->bTautFlagsDone     = t_group_info_orig->bTautFlagsDone;
    }
    return err;
}
/*******************************************************************************/
/*  set tautomer group isotopic sort keys */
int set_tautomer_iso_sort_keys( T_GROUP_INFO *t_group_info )
{
    T_GROUP       *t_group;
    T_GROUP_ISOWT Mult = 1;
    int     i, j, num_t_groups, num_iso_t_groups = 0;
    if ( !t_group_info || !(t_group = t_group_info->t_group) ||
         0 >= (num_t_groups = t_group_info->num_t_groups) || t_group_info->nNumIsotopicEndpoints )
        return 0;
    for ( i = 0; i < num_t_groups; i ++ ) {
        t_group[i].iWeight = 0;
        j = T_NUM_ISOTOPIC - 1;
        Mult = 1;
        do {
            t_group[i].iWeight += Mult * (T_GROUP_ISOWT)t_group[i].num[T_NUM_NO_ISOTOPIC+j];
        } while ( --j >= 0 && (Mult *= T_GROUP_ISOWT_MULT) );
        num_iso_t_groups += (t_group[i].iWeight != 0);
    }
    return num_iso_t_groups;
}

/******************************************************************************
 *
 *  Fill t_group_info with information necessary to fill out tautomer part
 *  of the linear connection table record.
 *  Note: on input, t_group_info should contain information created by MarkTautomerGroups()
 *        No previous t_group_info adjustment due to throwing out disconnected parts of
 *        the chemical structure is needed.
 *
 *  Note2: throws out t_groups containing negative charges only (IGNORE_TGROUP_WITHOUT_H==1)
 *         (leave their tautomeric bonds unchanged)
 *  Note3: removes negative charges from other tautomeric groups
 *         and adjust counts of mobile atoms if permitted         (REMOVE_TGROUP_CHARGE==1)
 */
int CountTautomerGroups( sp_ATOM *at, int num_atoms, T_GROUP_INFO *t_group_info )
{
    int i, j, ret = 0, nNumEndpoints, max_t_group, num_groups_noH;

    AT_NUMB    nGroupNumber, nNewGroupNumber, *nCurrEndpointAtNoPos = NULL;
    
    T_GROUP   *t_group;
    int        num_t;
    /* int bIgnoreIsotopic, max_num_t; */
    AT_NUMB   *nTautomerGroupNumber        = NULL;
    AT_NUMB   *nEndpointAtomNumber         = NULL;
    AT_NUMB   *tGroupNumber  = NULL;
    
    if ( !t_group_info || !t_group_info->t_group || 0 >= t_group_info->max_num_t_groups ) {
        return 0; /* empty t-groups */
    }
    num_t           =  t_group_info->num_t_groups;
    t_group         =  t_group_info->t_group;
    /*
      max_num_t       =  t_group_info->max_num_t_groups;
      bIgnoreIsotopic =  t_group_info->bIgnoreIsotopic;
     */
    num_groups_noH  = 0;

    /* the following 2 arrays are to be rebuilt here */
    if ( t_group_info->nEndpointAtomNumber ) {
        inchi_free ( t_group_info->nEndpointAtomNumber );
        t_group_info->nEndpointAtomNumber = NULL;
    }
    if ( t_group_info->tGroupNumber ) {
        inchi_free ( t_group_info->tGroupNumber );
        t_group_info->tGroupNumber = NULL;
    }
    /*  find max_t_group */
    for ( i = 0, max_t_group = 0; i < t_group_info->num_t_groups; i ++ ) {
        if ( max_t_group < t_group[i].nGroupNumber )
            max_t_group = t_group[i].nGroupNumber;
    }
    /*  allocate memory for temp storage of numbers of endpoints  */
    if ( max_t_group &&
         !(nTautomerGroupNumber = (AT_NUMB*) inchi_calloc( max_t_group+1, sizeof(nTautomerGroupNumber[0]) ) /*temp*/ ) ) {
        goto err_exit_function; /*  program error: out of RAM */ /*   <BRKPT> */
    }
    
    /*  count endpoints for each tautomer group */
    for ( i = 0, nNumEndpoints = 0; i < num_atoms; i ++ ) {
        if ( (j = at[i].endpoint) == 0 )
            continue;
        if ( j > max_t_group ) /*  debug only */
            goto err_exit_function; /*  program error */ /*   <BRKPT> */
        nTautomerGroupNumber[j] ++;
        nNumEndpoints ++;
    }
    
    if ( !nNumEndpoints ) {
        goto exit_function; /*  not a tautomer */
    }

    /*  allocate temporary array */
    if ( !(nEndpointAtomNumber  = (AT_NUMB*) inchi_calloc( nNumEndpoints, sizeof(nEndpointAtomNumber[0]) ) ) ||
         !(nCurrEndpointAtNoPos = (AT_NUMB*) inchi_calloc( num_t, sizeof(nCurrEndpointAtNoPos[0]) ) /*temp*/ ) ) {
        goto err_exit_function; /*   program error: out of RAM */ /*   <BRKPT> */
    }
    /*
     * Remove missing endpoints from t_group. Since only one
     * disconnected part is processed, some endpoints groups may have disappeared.
     * Mark t_groups containing charges only for subsequent removal
     */
    for ( i = 0, nNewGroupNumber = 0; i < num_t; /*i ++*/ ) {
        int bNoH = 0, nNumH;
        nGroupNumber  = t_group[i].nGroupNumber;
        for ( j = 1, nNumH = t_group[i].num[0]; j < T_NUM_NO_ISOTOPIC; j ++ ) {
            nNumH -= (int)t_group[i].num[j];
        }
        if ( t_group[i].nNumEndpoints != nTautomerGroupNumber[(int)nGroupNumber]
#if( IGNORE_TGROUP_WITHOUT_H == 1 )
             || (bNoH = (t_group[i].num[0]==t_group[i].num[1]))  /* only for (H,-) t-groups; (+) t-groups are not removed */
#endif
           ) {
            if ( !nTautomerGroupNumber[(int)nGroupNumber] || bNoH ) {
                /*  the group belongs to another disconnected part of the structure or has only charges */
                /*  Remove the group */
                num_t --;
                if ( i < num_t )
                    memmove( t_group+i, t_group+i+1, (num_t-i)*sizeof(t_group[0]) );
                if ( bNoH ) {
                    /*  group contains no mobile hydrogen atoms, only charges. Prepare to remove it. */
                    nTautomerGroupNumber[(int)nGroupNumber] = 0;
                    num_groups_noH ++;
                }
                /*i --;*/
            } else {
                /*  different number of endpoints */
                goto err_exit_function; /*  program error */ /*   <BRKPT> */
            }
        } else {
            /*  renumber t_group and prepare to renumber at[i].endpoint */
            nTautomerGroupNumber[(int)nGroupNumber] =
            t_group[i].nGroupNumber                 = ++nNewGroupNumber; /*  = i+1 */
            /*  get first group atom orig. number position in the nEndpointAtomNumber[] */
            /*  and in the tautomer endpoint canon numbers part of the connection table */
            t_group[i].nFirstEndpointAtNoPos = nCurrEndpointAtNoPos[i]  =
                i? (t_group[i-1].nFirstEndpointAtNoPos+t_group[i-1].nNumEndpoints) : 0;
            t_group[i].num[0] = nNumH;
#if( REMOVE_TGROUP_CHARGE == 1 )
            t_group[i].num[1]  = 0;  /* remove only (-) charges */
#endif
            /* -- wrong condition. Disabled.
            if ( t_group[i].nGroupNumber != i + 1 ) { // for debug only
                goto err_exit_function; // program error
            }
            */
            i ++;
        }
    }
    if ( num_t != nNewGroupNumber ) { /*  for debug only */
        goto err_exit_function; /*  program error */ /*   <BRKPT> */
    }
    
    /*  check if any tautomer group was left */
    if ( !nNewGroupNumber ) {
        if ( !num_groups_noH )
            goto err_exit_function; /*  program error: not a tautomer */ /*   <BRKPT> */
        else
            goto exit_function;
    }
    /*
     * an array for tautomer group sorting later, at the time of storing Connection Table
     * Later the sorting consists out of 2 steps:
     * 1) Sort t_group[i].nNumEndpoints endpoint atom ranks within each endpoint group
     *    starting from t_group[i].nFirstEndpointAtNoPos; i = 0..t_group_info->num_t_groups-1
     * 2) Sort the groups indexes t_group_info->tGroupNumber[]
     */
    if ( !(tGroupNumber=
           (AT_NUMB*)inchi_calloc(nNewGroupNumber*TGSO_TOTAL_LEN, sizeof(tGroupNumber[0])))) {
        goto err_exit_function; /*  out of RAM */
    }
    for ( i = 0; i < nNewGroupNumber; i ++ ) {
        tGroupNumber[i] = (AT_NUMB)i; /*  initialization: original t_group number = (at[i]->endpoint-1) */
    }
    /*
     * renumber endpoint atoms and save their orig. atom 
     * numbers for filling out the tautomer part of the LinearCT.
     * nCurrEndpointAtNoPos[j] is an index of the atom number in the nEndpointAtomNumber[]
     */
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( j = (int)at[i].endpoint ) {
            j = (int)(at[i].endpoint = nTautomerGroupNumber[j])-1; /*  new t_group number */
            if ( j >= 0 ) { /*  j=-1 in case of no mobile hydrogen atoms (charges only), group being removed */
                if ( nCurrEndpointAtNoPos[j] >=   /*  debug only */
                     t_group[j].nFirstEndpointAtNoPos+t_group[j].nNumEndpoints ) {
                    goto err_exit_function; /*  program error */ /*   <BRKPT> */
                }
                nEndpointAtomNumber[(int)nCurrEndpointAtNoPos[j] ++] = (AT_NUMB)i;
            } else {
                nNumEndpoints --; /*  endpoint has been removed */
            }
        }
    }
    t_group_info->num_t_groups               = nNewGroupNumber;
    t_group_info->nNumEndpoints              = nNumEndpoints;
    t_group_info->nEndpointAtomNumber        = nEndpointAtomNumber;
    t_group_info->tGroupNumber               = tGroupNumber; /* only the 1st segment filled */
    inchi_free ( nTautomerGroupNumber );
    inchi_free ( nCurrEndpointAtNoPos );
    return nNumEndpoints + T_GROUP_HDR_LEN * nNewGroupNumber + 1; /*  nLenLinearCTTautomer */

err_exit_function:
    ret = CT_TAUCOUNT_ERR;
exit_function:
    /*  release allocated memory; set "no tautomeric group" */
    if ( nEndpointAtomNumber )
        inchi_free ( nEndpointAtomNumber );
    if ( nTautomerGroupNumber )
        inchi_free ( nTautomerGroupNumber );
    if ( tGroupNumber )
        inchi_free ( tGroupNumber );
    if ( nCurrEndpointAtNoPos )
        inchi_free ( nCurrEndpointAtNoPos );
    t_group_info->nNumEndpoints = 0;
    t_group_info->num_t_groups  = 0;
    if ( !ret && ((t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT) ||
                   t_group_info->nNumIsotopicEndpoints>1 && (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))) ) {
        ret = 1; /* only protons have been (re)moved or neitralization happened */
    }
    return ret;
}
#if( READ_INCHI_STRING == 1 )
#if( INCLUDE_NORMALIZATION_ENTRY_POINT == 1 )
/********************************************************************************/
int CountTautomerGroupsInpAt( inp_ATOM *at, int num_atoms, T_GROUP_INFO *t_group_info )
{
    int i, j, ret = 0, nNumEndpoints, max_t_group, num_groups_noH;

    AT_NUMB    nGroupNumber, nNewGroupNumber, *nCurrEndpointAtNoPos = NULL;
    
    T_GROUP   *t_group;
    int        num_t;
    /* int bIgnoreIsotopic, max_num_t; */
    AT_NUMB   *nTautomerGroupNumber        = NULL;
    AT_NUMB   *nEndpointAtomNumber         = NULL;
    AT_NUMB   *tGroupNumber  = NULL;
    
    if ( !t_group_info || !t_group_info->t_group || 0 >= t_group_info->max_num_t_groups ) {
        return 0; /* empty t-groups */
    }
    num_t           =  t_group_info->num_t_groups;
    t_group         =  t_group_info->t_group;
    /*
      max_num_t       =  t_group_info->max_num_t_groups;
      bIgnoreIsotopic =  t_group_info->bIgnoreIsotopic;
     */
    num_groups_noH  = 0;

    /* the following 2 arrays are to be rebuilt here */
    if ( t_group_info->nEndpointAtomNumber ) {
        inchi_free ( t_group_info->nEndpointAtomNumber );
        t_group_info->nEndpointAtomNumber = NULL;
    }
    if ( t_group_info->tGroupNumber ) {
        inchi_free ( t_group_info->tGroupNumber );
        t_group_info->tGroupNumber = NULL;
    }
    /*  find max_t_group */
    for ( i = 0, max_t_group = 0; i < t_group_info->num_t_groups; i ++ ) {
        if ( max_t_group < t_group[i].nGroupNumber )
            max_t_group = t_group[i].nGroupNumber;
    }
    /*  allocate memory for temp storage of numbers of endpoints  */
    if ( max_t_group &&
         !(nTautomerGroupNumber = (AT_NUMB*) inchi_calloc( max_t_group+1, sizeof(nTautomerGroupNumber[0]) ) /*temp*/ ) ) {
        goto err_exit_function; /*  program error: out of RAM */ /*   <BRKPT> */
    }
    
    /*  count endpoints for each tautomer group */
    for ( i = 0, nNumEndpoints = 0; i < num_atoms; i ++ ) {
        if ( (j = at[i].endpoint) == 0 )
            continue;
        if ( j > max_t_group ) /*  debug only */
            goto err_exit_function; /*  program error */ /*   <BRKPT> */
        nTautomerGroupNumber[j] ++;
        nNumEndpoints ++;
    }
    
    if ( !nNumEndpoints ) {
        goto exit_function; /*  not a tautomer */
    }

    /*  allocate temporary array */
    if ( !(nEndpointAtomNumber  = (AT_NUMB*) inchi_calloc( nNumEndpoints, sizeof(nEndpointAtomNumber[0]) ) ) ||
         !(nCurrEndpointAtNoPos = (AT_NUMB*) inchi_calloc( num_t, sizeof(nCurrEndpointAtNoPos[0]) ) /*temp*/ ) ) {
        goto err_exit_function; /*   program error: out of RAM */ /*   <BRKPT> */
    }
    /*
     * Remove missing endpoints from t_group. Since only one
     * disconnected part is processed, some endpoints groups may have disappeared.
     * Mark t_groups containing charges only for subsequent removal
     */
    for ( i = 0, nNewGroupNumber = 0; i < num_t; /*i ++*/ ) {
        int bNoH = 0, nNumH;
        nGroupNumber  = t_group[i].nGroupNumber;
        for ( j = 1, nNumH = t_group[i].num[0]; j < T_NUM_NO_ISOTOPIC; j ++ ) {
            nNumH -= (int)t_group[i].num[j];
        }
        if ( t_group[i].nNumEndpoints != nTautomerGroupNumber[(int)nGroupNumber]
#if( IGNORE_TGROUP_WITHOUT_H == 1 )
             || (bNoH = (t_group[i].num[0]==t_group[i].num[1]))  /* only for (H,-) t-groups; (+) t-groups are not removed */
#endif
           ) {
            if ( !nTautomerGroupNumber[(int)nGroupNumber] || bNoH ) {
                /*  the group belongs to another disconnected part of the structure or has only charges */
                /*  Remove the group */
                num_t --;
                if ( i < num_t )
                    memmove( t_group+i, t_group+i+1, (num_t-i)*sizeof(t_group[0]) );
                if ( bNoH ) {
                    /*  group contains no mobile hydrogen atoms, only charges. Prepare to remove it. */
                    nTautomerGroupNumber[(int)nGroupNumber] = 0;
                    num_groups_noH ++;
                }
                /*i --;*/
            } else {
                /*  different number of endpoints */
                goto err_exit_function; /*  program error */ /*   <BRKPT> */
            }
        } else {
            /*  renumber t_group and prepare to renumber at[i].endpoint */
            nTautomerGroupNumber[(int)nGroupNumber] =
            t_group[i].nGroupNumber                 = ++nNewGroupNumber; /*  = i+1 */
            /*  get first group atom orig. number position in the nEndpointAtomNumber[] */
            /*  and in the tautomer endpoint canon numbers part of the connection table */
            t_group[i].nFirstEndpointAtNoPos = nCurrEndpointAtNoPos[i]  =
                i? (t_group[i-1].nFirstEndpointAtNoPos+t_group[i-1].nNumEndpoints) : 0;
            t_group[i].num[0] = nNumH;
#if( REMOVE_TGROUP_CHARGE == 1 )
            t_group[i].num[1]  = 0;  /* remove only (-) charges */
#endif
            /* -- wrong condition. Disabled.
            if ( t_group[i].nGroupNumber != i + 1 ) { // for debug only
                goto err_exit_function; // program error
            }
            */
            i ++;
        }
    }
    if ( num_t != nNewGroupNumber ) { /*  for debug only */
        goto err_exit_function; /*  program error */ /*   <BRKPT> */
    }
    
    /*  check if any tautomer group was left */
    if ( !nNewGroupNumber ) {
        if ( !num_groups_noH )
            goto err_exit_function; /*  program error: not a tautomer */ /*   <BRKPT> */
        else
            goto exit_function;
    }
    /*
     * an array for tautomer group sorting later, at the time of storing Connection Table
     * Later the sorting consists out of 2 steps:
     * 1) Sort t_group[i].nNumEndpoints endpoint atom ranks within each endpoint group
     *    starting from t_group[i].nFirstEndpointAtNoPos; i = 0..t_group_info->num_t_groups-1
     * 2) Sort the groups indexes t_group_info->tGroupNumber[]
     */
    if ( !(tGroupNumber=
           (AT_NUMB*)inchi_calloc(nNewGroupNumber*TGSO_TOTAL_LEN, sizeof(tGroupNumber[0])))) {
        goto err_exit_function; /*  out of RAM */
    }
    for ( i = 0; i < nNewGroupNumber; i ++ ) {
        tGroupNumber[i] = (AT_NUMB)i; /*  initialization: original t_group number = (at[i]->endpoint-1) */
    }
    /*
     * renumber endpoint atoms and save their orig. atom 
     * numbers for filling out the tautomer part of the LinearCT.
     * nCurrEndpointAtNoPos[j] is an index of the atom number in the nEndpointAtomNumber[]
     */
    for ( i = 0; i < num_atoms; i ++ ) {
        if ( j = (int)at[i].endpoint ) {
            j = (int)(at[i].endpoint = nTautomerGroupNumber[j])-1; /*  new t_group number */
            if ( j >= 0 ) { /*  j=-1 in case of no mobile hydrogen atoms (charges only), group being removed */
                if ( nCurrEndpointAtNoPos[j] >=   /*  debug only */
                     t_group[j].nFirstEndpointAtNoPos+t_group[j].nNumEndpoints ) {
                    goto err_exit_function; /*  program error */ /*   <BRKPT> */
                }
                nEndpointAtomNumber[(int)nCurrEndpointAtNoPos[j] ++] = (AT_NUMB)i;
            } else {
                nNumEndpoints --; /*  endpoint has been removed */
            }
        }
    }
    t_group_info->num_t_groups               = nNewGroupNumber;
    t_group_info->nNumEndpoints              = nNumEndpoints;
    t_group_info->nEndpointAtomNumber        = nEndpointAtomNumber;
    t_group_info->tGroupNumber               = tGroupNumber; /* only the 1st segment filled */
    inchi_free ( nTautomerGroupNumber );
    inchi_free ( nCurrEndpointAtNoPos );
    return nNumEndpoints + T_GROUP_HDR_LEN * nNewGroupNumber + 1; /*  nLenLinearCTTautomer */

err_exit_function:
    ret = CT_TAUCOUNT_ERR;
exit_function:
    /*  release allocated memory; set "no tautomeric group" */
    if ( nEndpointAtomNumber )
        inchi_free ( nEndpointAtomNumber );
    if ( nTautomerGroupNumber )
        inchi_free ( nTautomerGroupNumber );
    if ( tGroupNumber )
        inchi_free ( tGroupNumber );
    if ( nCurrEndpointAtNoPos )
        inchi_free ( nCurrEndpointAtNoPos );
    t_group_info->nNumEndpoints = 0;
    t_group_info->num_t_groups  = 0;
    if ( !ret && ((t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT) ||
                   t_group_info->nNumIsotopicEndpoints>1 && (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))) ) {
        ret = 1; /* only protons have been (re)moved or neitralization happened */
    }
    return ret;
}
#endif
#endif
/**************************************************************
 * tautomers: Compare for sorting
 ******************************************************************/
/*  Compare for sorting Ranks only */
/*  Globals: pn_tRankForSort */
int CompRankTautomer(const void* a1, const void* a2 )
{
    int ret = (int)pn_tRankForSort[(int)(*(const AT_RANK*)a1)] -
              (int)pn_tRankForSort[(int)(*(const AT_RANK*)a2)];
    return ret;
}
/*********************************************************************/
int SortTautomerGroupsAndEndpoints( T_GROUP_INFO *t_group_info, int num_atoms, int num_at_tg, AT_RANK *nRank )
{
    int i, nFirstEndpointAtNoPos, nNumEndpoints;
    AT_NUMB  *nEndpointAtomNumber;
    int       num_t_groups = num_at_tg - num_atoms;
    T_GROUP   *t_group = NULL;
    /*  check if sorting is required */
     
    if ( num_t_groups <= 0 || t_group_info->nNumEndpoints < 2 ) {
         return 0; /*  no tautomer data */
    }
    t_group = t_group_info->t_group;
    /*  sort endpoints within the groups */
    for ( i = 0; i < num_t_groups; i ++ ) {
        if ( t_group[i].nNumEndpoints < 2 )
            continue;  /*  program error; should not happen */ /*   <BRKPT> */
        /*  set globals for sorting */
        nFirstEndpointAtNoPos = t_group[i].nFirstEndpointAtNoPos;
        nNumEndpoints         = t_group[i].nNumEndpoints;
        if ( nNumEndpoints + nFirstEndpointAtNoPos > t_group_info->nNumEndpoints ) { /*  for debug only */
            return CT_TAUCOUNT_ERR; /*  program error */ /*   <BRKPT> */
        }
        nEndpointAtomNumber = t_group_info->nEndpointAtomNumber+(int)nFirstEndpointAtNoPos;
        pn_tRankForSort = nRank;
        insertions_sort( nEndpointAtomNumber, nNumEndpoints, sizeof(nEndpointAtomNumber[0]), CompRankTautomer);
    }
    /*  sort the tautomeric groups according to their ranks only
        (that is, ignoring the isotopic composition of the mobile groups and ranks of the endpoints) */
    if ( t_group_info->num_t_groups > 1 ) {
        /*  set globals for sorting */
        /*  a hack: the ranks of all tautomeric groups are */
        /*  located at nRank[num_atoms..num_at_tg-1] */
        pn_tRankForSort    = nRank+num_atoms;
        /*  sort */
        /*  ordering numbers to sort : t_group_info->tGroupNumber; */
        insertions_sort( t_group_info->tGroupNumber, num_t_groups,
                         sizeof(t_group_info->tGroupNumber[0]), CompRankTautomer);
    }
    return t_group_info->num_t_groups;
}
