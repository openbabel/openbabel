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
#include <ctype.h>
#include <math.h>

/*^^^ */
/* #define CHECK_WIN32_VC_HEAP */
#include "mode.h"

#if ( READ_INCHI_STRING == 1 )

#include "ichicomp.h"
#include "ichi.h"
#include "ichitime.h"
#include "ichierr.h"
#include "util.h"
#include "strutil.h"

/* reverse InChI */
#include "ichimain.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichister.h"
#include "strutil.h"
#include "ichisize.h"
#include "ichiring.h"
#include "ichinorm.h"

#include "ichirvrs.h"
#include "inchicmp.h"

/******************************************************************************************************/
int InChI2Atom( ICHICONST INPUT_PARMS *ip,  STRUCT_DATA *sd, const char *szCurHdr, long num_inp,
                StrFromINChI *pStruct, int iComponent, int iAtNoOffset, int  bI2A_Flag, int bHasSomeFixedH, InpInChI *OneInput)
{
    int           iINChI   = (bI2A_Flag & I2A_FLAG_RECMET)? INCHI_REC : INCHI_BAS;
    int           bMobileH = (bI2A_Flag & I2A_FLAG_FIXEDH)? TAUT_NON  : TAUT_YES;
    INChI        *pInChI[TAUT_NUM];
    int           ret = 0;

    memset( pInChI, 0, sizeof(pInChI) );
    /* disconnected or reconnected */
    if ( iINChI == INCHI_REC ) {
        if ( !OneInput->nNumComponents[iINChI][TAUT_YES] ) {
            iINChI = INCHI_BAS;
        }
    }
    if ( iComponent >= OneInput->nNumComponents[iINChI][TAUT_YES] ) {
        return 0; /* component does not exist */
    }
    /* mobile or fixed H */
    pStruct->bFixedHExists = 0;
    if ( bMobileH == TAUT_NON ) {
        if ( !OneInput->nNumComponents[iINChI][bMobileH] ) {
            /* only one InChI exists (no mobile H) */
            bMobileH = TAUT_YES;
        }
    } 
    
    if ( iComponent >= OneInput->nNumComponents[iINChI][bMobileH] ) {
        return 0; /* component does not exist */
    }
    /* pointer to the InChI that is going to be reversed */
    pInChI[0] = &OneInput->pInpInChI[iINChI][bMobileH][iComponent];
    pStruct->bMobileH = bMobileH;
    pStruct->iINCHI   = iINChI;
    /* deleted component only in case Mobile-H and compound contains only protons */
    if ( pInChI[0]->bDeleted ) {
        return 0; /* deleted component, presumably H(+) */
    }

    if ( bMobileH == TAUT_NON && OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons ) {
        pStruct->nNumRemovedProtonsMobHInChI = 
            OneInput->nNumProtons[iINChI][TAUT_YES].pNumProtons[iComponent].nNumRemovedProtons;
    }

    if ( bMobileH == TAUT_NON || bMobileH == TAUT_YES && OneInput->pInpInChI[iINChI][TAUT_NON] &&
         OneInput->pInpInChI[iINChI][TAUT_NON][iComponent].nNumberOfAtoms > 0 &&
        !OneInput->pInpInChI[iINChI][TAUT_NON][iComponent].bDeleted ) {
        pStruct->bFixedHExists = 1;
    }
    if ( bMobileH == TAUT_NON && iComponent < OneInput->nNumComponents[iINChI][TAUT_YES] &&
         OneInput->pInpInChI[iINChI][TAUT_YES] &&
         OneInput->pInpInChI[iINChI][TAUT_YES][iComponent].nNumberOfAtoms > 0 &&
         !OneInput->pInpInChI[iINChI][TAUT_YES][iComponent].bDeleted ) {
        /* pointer to the Mobile-H InChI if we are reversing Fixed-H InChI */
        pInChI[1] = &OneInput->pInpInChI[iINChI][TAUT_YES][iComponent];
    }
    pStruct->num_inp_actual = OneInput->num_inp;
    ret = OneInChI2Atom( ip, sd, szCurHdr, num_inp, pStruct, iComponent, iAtNoOffset, bHasSomeFixedH, pInChI);
    return ret; /* same interpretation as in ProcessOneStructure ??? */
}

/*******************************************************************/
void RemoveFixHInChIIdentical2MobH( InpInChI *pOneInput )
{
    int iInchiRec, cur_num_comp, k;
    /* eliminate Fixed-H InChI that are exactly came as the corresponding Mobile-H structures */
    for ( iInchiRec = 0; iInchiRec < INCHI_NUM; iInchiRec ++ ) {
        cur_num_comp = inchi_min(pOneInput->nNumComponents[iInchiRec][TAUT_YES],
                                 pOneInput->nNumComponents[iInchiRec][TAUT_NON]);
        for ( k = 0; k < cur_num_comp; k ++ ) {
            if ( !CompareReversedINChI( pOneInput->pInpInChI[iInchiRec][TAUT_YES]+k,
                                        pOneInput->pInpInChI[iInchiRec][TAUT_NON]+k, NULL, NULL ) ) {
                Free_INChI_Members( pOneInput->pInpInChI[iInchiRec][TAUT_NON]+k );
                memset( pOneInput->pInpInChI[iInchiRec][TAUT_NON]+k, 0, sizeof(pOneInput->pInpInChI[0][0][0]) );
            }
        }
    }
}
/*******************************************************************/
int MarkDisconectedIdenticalToReconnected ( InpInChI *pOneInput )
{
    /* mark Disconnected InChI components that are exactly came as Reconnected ones */
    /* Disconnected will have a negative number of the reconnected component */
    /* Reconnected will have a positive number of the disconnected component */
    int k1, k2, num_marked = 0;
    for ( k1 = 0; k1 < inchi_max(pOneInput->nNumComponents[INCHI_BAS][TAUT_YES],
                                 pOneInput->nNumComponents[INCHI_BAS][TAUT_NON]); k1 ++ ) {
    for ( k2 = 0; k2 < inchi_max(pOneInput->nNumComponents[INCHI_REC][TAUT_YES],
                                 pOneInput->nNumComponents[INCHI_REC][TAUT_NON]); k2 ++ ) {
        int eqM = ( k1 < pOneInput->nNumComponents[INCHI_BAS][TAUT_YES] &&
                    k2 < pOneInput->nNumComponents[INCHI_REC][TAUT_YES] &&
                   !pOneInput->pInpInChI[INCHI_REC][TAUT_YES][k2].nLink && /* already linked */
                   !pOneInput->pInpInChI[INCHI_BAS][TAUT_YES][k1].bDeleted &&
                    pOneInput->pInpInChI[INCHI_BAS][TAUT_YES][k1].nNumberOfAtoms &&
                    pOneInput->pInpInChI[INCHI_BAS][TAUT_YES][k1].nNumberOfAtoms ==
                    pOneInput->pInpInChI[INCHI_REC][TAUT_YES][k2].nNumberOfAtoms &&
                   !pOneInput->pInpInChI[INCHI_REC][TAUT_YES][k2].bDeleted &&
                    !CompareReversedINChI( pOneInput->pInpInChI[INCHI_REC][TAUT_YES]+k2,
                                           pOneInput->pInpInChI[INCHI_BAS][TAUT_YES]+k1,
                                           NULL, NULL ));
        int isF1 = (k1 < pOneInput->nNumComponents[INCHI_BAS][TAUT_NON] &&
                    0 == pOneInput->pInpInChI[INCHI_BAS][TAUT_NON][k1].bDeleted &&
                    0  < pOneInput->pInpInChI[INCHI_BAS][TAUT_NON][k1].nNumberOfAtoms );
        int isF2 = (k2 < pOneInput->nNumComponents[INCHI_REC][TAUT_NON] &&
                    0 == pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].bDeleted &&
                    0  < pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].nNumberOfAtoms );
        int eqF =  isF1 && isF2 &&
                   !pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].nLink &&
                    pOneInput->pInpInChI[INCHI_BAS][TAUT_NON][k1].nNumberOfAtoms ==
                    pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].nNumberOfAtoms &&
                    !CompareReversedINChI( pOneInput->pInpInChI[INCHI_REC][TAUT_NON]+k2,
                                           pOneInput->pInpInChI[INCHI_BAS][TAUT_NON]+k1,
                                           NULL, NULL );
        if ( eqM && (!isF1 && !isF2 || eqF ) ) {
            pOneInput->pInpInChI[INCHI_BAS][TAUT_YES][k1].nLink = -(k2+1);
            pOneInput->pInpInChI[INCHI_REC][TAUT_YES][k2].nLink =  (k1+1);
            if ( eqF ) {
                pOneInput->pInpInChI[INCHI_BAS][TAUT_NON][k1].nLink = -(k2+1);
                pOneInput->pInpInChI[INCHI_REC][TAUT_NON][k2].nLink =  (k1+1);
            }
            num_marked ++;
            break; /* equal InChI has been deleted from the disconnected layer, get next k1 */
        }
    }
    }
    return num_marked;

}
/**************************************************************/
void SetUpSrm( SRM *pSrm )
{
    /* structure restore parms !!!!! */
    memset( pSrm, 0, sizeof(pSrm[0]) );
    pSrm->bFixStereoBonds      = FIX_STEREO_BOND_ORDER;
    pSrm->nMetal2EndpointMinBondOrder  = 1;
    pSrm->nMetal2EndpointInitEdgeFlow  = 0;
    if ( METAL_FREE_CHARGE_VAL == 1 ) { 
        pSrm->bMetalAddFlower      = 1;
        /* the next 3 parameters: */
        /* 0, 0, 0 => all bonds 0, no init radical on metal */
        /* 0, 0, 1 => all bonds 0,    init radical on metal */
        /* 0, 1, 0 => wrong */
        /* 0, 1, 1 => all bonds 1, no init radical on metal */
        /* 1, 0, 1 => min bond order 1, all bonds to metal have order 1 */
        /* 1, 1, 0 => wrong */
        /* 1, 1, 1 => wrong */
        pSrm->nMetalMinBondOrder   = 0; 
        pSrm->nMetalInitEdgeFlow   = 1;  
        pSrm->nMetalInitBondOrder  = 1; 
        pSrm->bStereoRemovesMetalFlag = pSrm->bFixStereoBonds;
        pSrm->nMetalFlowerParam_D     = 16;
        pSrm->nMetalMaxCharge_D       = 16;
    } else {
        pSrm->bMetalAddFlower      = 0;
        pSrm->nMetalMinBondOrder   = 1;
        pSrm->nMetalInitEdgeFlow   = 0;
        pSrm->nMetalInitBondOrder  = 1;
        pSrm->bStereoRemovesMetalFlag = pSrm->bFixStereoBonds;
        pSrm->nMetalFlowerParam_D     = 16;
        pSrm->nMetalMaxCharge_D       = 0;
    }
    /*
    pSrm->nMetalInitBondOrder  = pSrm->nMetalMinBondOrder 
                             + pSrm->nMetalInitEdgeFlow;
    */
    pSrm->nMetal2EndpointInitBondOrder = pSrm->nMetal2EndpointMinBondOrder 
                                     + pSrm->nMetal2EndpointInitEdgeFlow;

}
/**************************************************************************************/
int MergeStructureComponents( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, long num_inp, char *szCurHdr,
                         ICHICONST SRM *pSrm, int bReqNonTaut, StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM],
                         InpInChI *pOneInput )
{
    int iInchiRec, iMobileH, iAlternH, num_components, tot_just_atoms, tot_removed_H, tot_atoms, cur_nA, cur_nH;
    int k, i, j, ret, iCurAtomOffs, iNxtAtomOffs, iCurDelHOffs, iNxtDelHOffs, len, len2, iShiftH, icomp;
    int *nAtomOffs=NULL, *nDelHOffs=NULL;
    StrFromINChI *pStruct1;
    inp_ATOM *at=NULL, *a;

    ret = 0;
    pOneInput->num_atoms = 0;
    /* select highest detail level */
    if ( num_components = pOneInput->nNumComponents[INCHI_REC][TAUT_NON] ) {
        iInchiRec = INCHI_REC;
        iMobileH  = TAUT_NON;
    } else
    if ( num_components = pOneInput->nNumComponents[INCHI_REC][TAUT_YES] ) {
        iInchiRec = INCHI_REC;
        iMobileH  = TAUT_YES;
    } else
    if ( num_components = pOneInput->nNumComponents[INCHI_BAS][TAUT_NON] ) {
        iInchiRec = INCHI_BAS;
        iMobileH  = TAUT_NON;
    } else
    if ( num_components = pOneInput->nNumComponents[INCHI_BAS][TAUT_YES] ) {
        iInchiRec = INCHI_BAS;
        iMobileH  = TAUT_YES;
    } else {
        return 0; /* no components available */
    }

    nAtomOffs = (int*) inchi_malloc((num_components+1) * sizeof(nAtomOffs[0]));
    nDelHOffs = (int*) inchi_malloc((num_components+1) * sizeof(nDelHOffs[0]));
    if ( !nAtomOffs || !nDelHOffs ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    /* count number of atoms and removed H */
    tot_just_atoms = tot_removed_H = tot_atoms = 0;
    iAlternH = (iMobileH==TAUT_NON && pOneInput->nNumComponents[iInchiRec][TAUT_YES])? TAUT_YES : -1;
    nAtomOffs[0] = nDelHOffs[0] = 0;
    for ( k = 0; k < num_components; k ++ ) {
        pStruct1 = pStruct[iInchiRec][iMobileH][k].num_atoms? pStruct[iInchiRec][iMobileH]+k :
                   iAlternH>=0 &&
                   pStruct[iInchiRec][iAlternH][k].num_atoms? pStruct[iInchiRec][iAlternH]+k : NULL;
        if ( !pStruct1 || !pStruct1->at2 || !pStruct1->num_atoms || pStruct1->bDeleted ) {
            cur_nA = cur_nH = 0;
        } else {
            cur_nA = pStruct1->num_atoms;
            cur_nH = pStruct1->num_deleted_H;
        }
        nAtomOffs[k+1] = nAtomOffs[k] + cur_nA;
        nDelHOffs[k+1] = nDelHOffs[k] + cur_nH;
    }
    tot_just_atoms = nAtomOffs[num_components];
    /* shift all H to the end */
    for ( k = 0; k <= num_components; k ++ ) {
        nDelHOffs[k] += tot_just_atoms;
    }
    tot_atoms = nDelHOffs[num_components];

    /* merge atoms together: 1. Allocate */
    if ( NULL == (at = (inp_ATOM *) inchi_malloc( (tot_atoms+1) * sizeof(at[0]) ) ) ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    if ( !tot_atoms ) {
        ret = 0;
        goto exit_function; /* empty structure */
    }
    /* merge atoms together: 2. Copy */
    for ( k = 0; k < num_components; k ++ ) {
        pStruct1 = pStruct[iInchiRec][iMobileH][k].num_atoms? pStruct[iInchiRec][iMobileH]+k :
                   iAlternH>=0 &&
                   pStruct[iInchiRec][iAlternH][k].num_atoms? pStruct[iInchiRec][iAlternH]+k : NULL;
        if ( len = nAtomOffs[k+1] - nAtomOffs[k] ) {
            memcpy( at + nAtomOffs[k], pStruct1->at2, len * sizeof(at[0]) );
            if ( len2 = nDelHOffs[k+1] - nDelHOffs[k] ) {
                memcpy( at + nDelHOffs[k], pStruct1->at2+len, len2 * sizeof(at[0]) );
            }
        }
    }
    /* merge atoms together: 3. Update atom numbers */
    icomp = 0;
    for ( k = 0; k < num_components; k ++ ) {
        iCurAtomOffs = nAtomOffs[k];
        iNxtAtomOffs = nAtomOffs[k+1];
        iCurDelHOffs = nDelHOffs[k];
        iNxtDelHOffs = nDelHOffs[k+1];
        len = nAtomOffs[k+1] - nAtomOffs[k]; /* number of atoms in a component excluding explicit H */
        iShiftH      = iCurDelHOffs - len;
        if ( !len ) {
            continue;
        }
        icomp ++; /* current component number */
        /* update atoms */
        for ( i = iCurAtomOffs; i < iNxtAtomOffs; i ++ ) {
            
            a = at+i;

            a->endpoint = 0;
            a->bAmbiguousStereo = 0;
            a->at_type = 0;
            a->bCutVertex = 0;
            a->bUsed0DParity = 0;
            a->cFlags = 0;
            a->nBlockSystem = 0;
            a->nNumAtInRingSystem = 0;
            a->nRingSystem = 0;
           
            for ( j = 0; j < a->valence; j ++ ) {
                if ( a->neighbor[j] < len ) {
                    a->neighbor[j] += iCurAtomOffs; /* atom */
                } else {
                    a->neighbor[j] += iShiftH;      /* explicit H */
                }
            }
            a->orig_at_number += iCurAtomOffs;
            a->component = icomp;
            if ( a->p_parity ) {
                for ( j = 0; j < MAX_NUM_STEREO_ATOM_NEIGH; j ++ ) {
                    if ( a->p_orig_at_num[j] <= len ) {
                        /* originally, orig_at_num = atom_index+1, therefore <= instead of < */
                        a->p_orig_at_num[j] += iCurAtomOffs;
                    } else {
                        a->p_orig_at_num[j] += iShiftH;
                    }
                }
            }
            for ( j = 0; j < MAX_NUM_STEREO_BONDS && a->sb_parity[j]; j ++ ) {
                if ( a->sn_orig_at_num[j] <= len ) {
                    /* originally, orig_at_num = atom_index+1, therefore <= instead of < */
                    a->sn_orig_at_num[j] += iCurAtomOffs;
                } else {
                    a->sn_orig_at_num[j] += iShiftH;
                }
            }
        }
        /* update fixed-H */
        for ( i = iCurDelHOffs; i < iNxtDelHOffs; i ++ ) {
            a = at+i;
            a->neighbor[0]    += iCurAtomOffs;
            a->orig_at_number += iShiftH;

        }
    }
    /* save the results */
    pOneInput->atom      = at;
    pOneInput->num_atoms = tot_atoms;
    at = NULL;

exit_function:
    if ( at )        inchi_free( at );  /* in case of failure */
    if ( nAtomOffs ) inchi_free( nAtomOffs );
    if ( nDelHOffs ) inchi_free( nDelHOffs );
    return ret;
}
#ifndef COMPILE_ANSI_ONLY
static PER_DRAW_PARMS pdp;
/******************************************************************************************************/
int DisplayAllRestoredComponents( inp_ATOM *at, int num_at, const char *szCurHdr )
{
    int    ret;
    char     szTitle[512];
    DRAW_PARMS dp;
    TBL_DRAW_PARMS tdp;
    if ( num_at <= 0 ) {
        return 0;
    }
    memset( &dp, 0, sizeof(dp));
    memset( &tdp, 0, sizeof(tdp) );
    //memset( &pdp, 0, sizeof(pdp) );
    dp.sdp.tdp       = &tdp;
    dp.pdp           = &pdp;
    dp.sdp.nFontSize = -9;
    sprintf( szTitle, "All Components of Restored %s Structure", szCurHdr? szCurHdr : "(No structure name)");
    ret = DisplayStructure( at, num_at, 0 /* nNumDeletedH*/, 0 /*bAdd_DT_to_num_H*/,
                      0 /*nNumRemovedProtons*/, NULL /*NUM_H *nNumRemovedProtonsIsotopic*/,
                      1 /*int bIsotopic*/, 0 /*bTautomeric*/,
                      NULL /* pINChI */, NULL /* INChI_Aux **cur_INChI_Aux*/,
                      0 /*bAbcNumbers*/, &dp, 0 /*INCHI_MODE nMode*/, szTitle );
    return 0;
}
/******************************************************************************************************/
int DisplayOneRestoredComponent( StrFromINChI *pStruct, inp_ATOM *at,
                                 int iComponent, int nNumComponents, int bMobileH,
                                 const char *szCurHdr )
{
    int    ret, k;
    int    num_at        = pStruct->num_atoms;
    XYZ_COORD *pxyz      = pStruct->pXYZ;
    char     szTitle[512];
    DRAW_PARMS dp;
    TBL_DRAW_PARMS tdp;
    int         iInchiRec = pStruct->iInchiRec;
    int         iMobileH  = pStruct->iMobileH;
    INChI     **pInChI = NULL;
    INChI_Aux **pAux   = NULL;
    int         nNumRemovedProtons         = pAux? pAux[iMobileH]->nNumRemovedProtons : 0;
    NUM_H      *nNumRemovedProtonsIsotopic = pAux? pAux[iMobileH]->nNumRemovedIsotopicH : NULL;


    if ( num_at <= 0 || !pxyz ) {
        return 0;
    }
    if ( iInchiRec && !pStruct->RevInChI.pINChI_Aux[iInchiRec][0] ) {
        iInchiRec = 0;
    }
    k = iMobileH;
    if ( !bRevInchiComponentExists( pStruct, iInchiRec, k, 0 ) ) {
        k = ALT_TAUT(k);
    }
    pInChI = pStruct->RevInChI.pINChI[iInchiRec][0];
    pAux   = pStruct->RevInChI.pINChI_Aux[iInchiRec][0];
    

    memset( &dp, 0, sizeof(dp));
    memset( &tdp, 0, sizeof(tdp) );
    //memset( &pdp, 0, sizeof(pdp) );
    dp.sdp.tdp       = &tdp;
    dp.pdp           = &pdp;
    dp.sdp.nFontSize = -9;
    sprintf( szTitle, "Restored %s Component %d of %d %c%c",
                      szCurHdr? szCurHdr : "(No structure name)", iComponent+1, nNumComponents,
                      pStruct->iInchiRec? 'R':'D', pStruct->iMobileH?'M':'F' );
    ret = DisplayStructure( at, num_at, 0 /* nNumDeletedH*/, 0 /*bAdd_DT_to_num_H*/,
                      nNumRemovedProtons, /*NULL*/ nNumRemovedProtonsIsotopic,
                      1 /*int bIsotopic*/, k,
                      pInChI, pAux,
                      0 /*bAbcNumbers*/, &dp, 0 /*INCHI_MODE nMode*/, szTitle );
    return 0;
}
/******************************************************************************************************/
int DisplayRestoredComponent( StrFromINChI *pStruct, int iComponent, int iAtNoOffset, INChI *pInChI, const char *szCurHdr )
{
    int    i, ret;
    int    num_at        = pStruct->num_atoms;
    int    num_deleted_H = pStruct->num_deleted_H;
    inp_ATOM *atom       = pStruct->at2;
    XYZ_COORD *pxyz      = pStruct->pXYZ;
    inp_ATOM *at         = NULL;
    char     szTitle[512];
    DRAW_PARMS dp;
    TBL_DRAW_PARMS tdp;
    if ( !atom || num_at <= 0 || !pxyz ) {
        return 0;
    }
    at = (inp_ATOM *)inchi_calloc( num_at + num_deleted_H, sizeof(at[0]) );
    if ( !at ) {
        return RI_ERR_ALLOC;
    }
    memcpy( at, atom, (num_at + num_deleted_H) * sizeof(at[0]) );
    for ( i = 0; i < num_at; i ++ ) {
        at[i].x = pxyz[i].xyz[0]; 
        at[i].y = pxyz[i].xyz[1]; 
        at[i].z = pxyz[i].xyz[2]; 
    }
    memset( &dp, 0, sizeof(dp));
    memset( &tdp, 0, sizeof(tdp) );
    //memset( &pdp, 0, sizeof(pdp) );
    dp.sdp.tdp       = &tdp;
    dp.pdp           = &pdp;
    dp.sdp.nFontSize = -9;
    sprintf( szTitle, "DBG Restored %s Component %d %c%c", szCurHdr? szCurHdr : "(No structure name)", iComponent+1, pStruct->iInchiRec? 'R':'D', pStruct->iMobileH?'M':'F' );
    ret = DisplayStructure( at, num_at, 0 /* nNumDeletedH*/, 0 /*bAdd_DT_to_num_H*/,
                      0 /*nNumRemovedProtons*/, NULL /*NUM_H *nNumRemovedProtonsIsotopic*/,
                      1 /*int bIsotopic*/, 0 /*bTautomeric*/,
                      &pInChI, NULL /* INChI_Aux **cur_INChI_Aux*/,
                      0 /*bAbcNumbers*/, &dp, 0 /*INCHI_MODE nMode*/, szTitle );
    inchi_free( at );
    return 0;
}
/**************************************************************************************/
int DisplayStructureComponents( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, long num_inp, char *szCurHdr,
                         ICHICONST SRM *pSrm, int bReqNonTaut, StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM],
                         InpInChI *pOneInput )
{
    int iInchiRec, iMobileH, iCurMobH, iAlternH, num_components, tot_just_atoms, tot_removed_H, tot_atoms, cur_nA, cur_nH;
    int k, i, j, ret, iCurAtomOffs, iNxtAtomOffs, iCurDelHOffs, iNxtDelHOffs, len, len2, iShiftH, icomp;
    int *nAtomOffs=NULL, *nDelHOffs=NULL, bNoCoord=0, iNewCoord=0, nNewCoord=0;
    double x_max=-1.0e16, x_min = 1.0e16, y_max=-1.0e16, y_min=1.0e16, delta = 0.0;
    StrFromINChI *pStruct1;
    inp_ATOM *at=NULL, *a;

    if (!ip->bDisplayCompositeResults && !ip->bDisplay ) {
        return 0;
    }

    ret = 0;
    pOneInput->num_atoms = 0;
    /* select highest detail level */
    if ( num_components = pOneInput->nNumComponents[INCHI_REC][TAUT_NON] ) {
        iInchiRec = INCHI_REC;
        iMobileH  = TAUT_NON;
    } else
    if ( num_components = pOneInput->nNumComponents[INCHI_REC][TAUT_YES] ) {
        iInchiRec = INCHI_REC;
        iMobileH  = TAUT_YES;
    } else
    if ( num_components = pOneInput->nNumComponents[INCHI_BAS][TAUT_NON] ) {
        iInchiRec = INCHI_BAS;
        iMobileH  = TAUT_NON;
    } else
    if ( num_components = pOneInput->nNumComponents[INCHI_BAS][TAUT_YES] ) {
        iInchiRec = INCHI_BAS;
        iMobileH  = TAUT_YES;
    } else {
        return 0; /* no components available */
    }
    for ( k = 0; k < num_components; k ++ ) {
        if ( pStruct[iInchiRec][iMobileH][k].bDeleted )
            break;
    }
    num_components = k;

    nAtomOffs = (int*) inchi_malloc((num_components+1) * sizeof(nAtomOffs[0]));
    nDelHOffs = (int*) inchi_malloc((num_components+1) * sizeof(nDelHOffs[0]));
    if ( !nAtomOffs || !nDelHOffs ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    /* count number of atoms and removed H */
    tot_just_atoms = tot_removed_H = tot_atoms = 0;
    iAlternH = (iMobileH==TAUT_NON && pOneInput->nNumComponents[iInchiRec][TAUT_YES])? TAUT_YES : -1;
    nAtomOffs[0] = nDelHOffs[0] = 0;
    for ( k = 0; k < num_components; k ++ ) {
        pStruct1 = pStruct[iInchiRec][iMobileH][k].num_atoms? pStruct[iInchiRec][iMobileH]+k :
                   iAlternH>=0 &&
                   pStruct[iInchiRec][iAlternH][k].num_atoms? pStruct[iInchiRec][iAlternH]+k : NULL;
        if ( !pStruct1 || !pStruct1->at2 || !pStruct1->num_atoms ) {
            cur_nA = cur_nH = 0;
        } else {
            cur_nA = pStruct1->num_atoms;
            cur_nH = pStruct1->num_deleted_H;
            if ( cur_nA && !pStruct1->pXYZ ) {
                if ( !k ) {
                    ret = 0; /* no coordinates available */
                    goto exit_function;
                } else {
                    bNoCoord ++;
                }
            }
        }
        nAtomOffs[k+1] = nAtomOffs[k] + cur_nA;
        nDelHOffs[k+1] = nDelHOffs[k] + cur_nH;
    }
    tot_just_atoms = nAtomOffs[num_components];
    /* shift all H to the end */
    for ( k = 0; k <= num_components; k ++ ) {
        nDelHOffs[k] += tot_just_atoms;
    }
    tot_atoms = nDelHOffs[num_components];

    /* merge atoms together: 1. Allocate */
    if ( NULL == (at = (inp_ATOM *) inchi_malloc( (tot_atoms+1) * sizeof(at[0]) ) ) ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }
    if ( !tot_atoms ) {
        ret = 0;
        goto exit_function; /* empty structure */
    }
    /* merge atoms together: 2. Copy */
    for ( k = 0; k < num_components; k ++ ) {
        pStruct1 = pStruct[iInchiRec][iMobileH][k].num_atoms? pStruct[iInchiRec][iCurMobH=iMobileH]+k :
                   iAlternH>=0 &&
                   pStruct[iInchiRec][iAlternH][k].num_atoms? pStruct[iInchiRec][iCurMobH=iAlternH]+k : NULL;
        if ( len = nAtomOffs[k+1] - nAtomOffs[k] ) {
            XYZ_COORD *pxyz = pStruct1->pXYZ;
            len2 = nDelHOffs[k+1] - nDelHOffs[k]; /* do not separate H from the atom: we will not need them */
            iCurAtomOffs = nAtomOffs[k];
            a = at + iCurAtomOffs;
            memcpy( a, pStruct1->at2, (len+len2) * sizeof(at[0]) );
            DisconnectedConnectedH( a, len, len2 );
            if ( pxyz ) {
                for ( i = 0; i < len; i ++ ) {
                    a[i].x = pxyz[i].xyz[0];
                    x_max = inchi_max( x_max, pxyz[i].xyz[0] );
                    x_min = inchi_min( x_min, pxyz[i].xyz[0] );
                    a[i].y = pxyz[i].xyz[1];
                    y_max = inchi_max( y_max, pxyz[i].xyz[1] );
                    y_min = inchi_min( y_min, pxyz[i].xyz[1] );
                    a[i].z = pxyz[i].xyz[2];
                    nNewCoord ++;
                }
            } else {
                if ( !iNewCoord ) {
                    if ( !nNewCoord ) {
                        ret = 0;
                        goto exit_function; /* empty structure */
                    }
                    delta = inchi_max(x_max - x_min, y_max - y_min);
                    if ( delta == 0.0 ) {
                        delta = 0.5 * (x_max+x_min);
                        if ( delta == 0.0 )
                            delta = 1.0;
                    } else {
                        delta /= sqrt( (double)(nNewCoord+1) );
                    }
                }
                for ( i = 0; i < len; i ++ ) {
                    a[i].x = x_max + delta;
                    a[i].y = y_max - iNewCoord * delta;
                    a[i].z = 0.0;
                    iNewCoord ++;
                }
                if ( pStruct1->pXYZ = (XYZ_COORD *)inchi_calloc(len, sizeof(pStruct1->pXYZ[0]) ) ) {

                    for ( i = 0; i < len; i ++ ) {
                        pStruct1->pXYZ[i].xyz[0] = a[i].x;
                        pStruct1->pXYZ[i].xyz[1] = a[i].y;
                        pStruct1->pXYZ[i].xyz[2] = 0.0;
                    }
                }
            }
            if ( ip->bDisplay || ip->bDisplayCompositeResults && 1 == num_components ) {
               DisplayOneRestoredComponent( pStruct1, a, k, num_components, iCurMobH, szCurHdr );
            }
            if ( !pxyz && pStruct1->pXYZ ) {
                inchi_free( pStruct1->pXYZ );
                pStruct1->pXYZ = NULL;
            }
        }
    }
    /* merge atoms together: 3. Update atom numbers */
    icomp = 0;
    if ( ip->bDisplayCompositeResults && num_components > 1 ) {
        for ( k = 0; k < num_components; k ++ ) {
            /* display each restored component if requested */
            iCurAtomOffs = nAtomOffs[k];
            iNxtAtomOffs = nAtomOffs[k+1];
            iCurDelHOffs = nDelHOffs[k];
            iNxtDelHOffs = nDelHOffs[k+1];
            len = nAtomOffs[k+1] - nAtomOffs[k]; /* number of atoms in a component excluding explicit H */
            iShiftH      = iCurDelHOffs - len;
            if ( !len ) {
                continue;
            }
            icomp ++; /* current component number */
            /* update atoms */
            for ( i = iCurAtomOffs; i < iNxtAtomOffs; i ++ ) {
                a = at+i;
                for ( j = 0; j < a->valence; j ++ ) {
                    if ( a->neighbor[j] < len ) {
                        a->neighbor[j] += iCurAtomOffs; /* atom */
                    } else {
                        ret = RI_ERR_PROGR;  /* explicit H */
                        goto exit_function;
                    }
                }
                a->orig_at_number += iCurAtomOffs;
            }
        }
        tot_atoms = nAtomOffs[num_components];
        DisplayAllRestoredComponents( at, tot_atoms, szCurHdr );

    }

exit_function:
    if ( at )        inchi_free( at );  /* in case of failure */
    if ( nAtomOffs ) inchi_free( nAtomOffs );
    if ( nDelHOffs ) inchi_free( nDelHOffs );
    return ret;
}
#endif
/**************************************************************************************/
int AllInchiToStructure( ICHICONST INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, long num_inp, char *szCurHdr,
                         ICHICONST SRM *pSrm, int bHasSomeFixedH, StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM],
                         InpInChI *pOneInput )
{
    int iInchiRec, iMobileH, cur_num_comp, bCurI2A_Flag, k, ret, num_err;
    INPUT_PARMS *ip, ip_loc;
    STRUCT_DATA *sd, sd_loc;
    long          ulProcessingTime = 0;
    inchiTime     ulTStart;

    InchiTimeGet( &ulTStart );
    ip  = &ip_loc;
    *ip = *ip_inp;
    sd  = &sd_loc;
    memset( sd, 0, sizeof(*sd));
    sd->ulStructTime = sd_inp->ulStructTime;
    ret = 0;
    num_err = 0;
    for ( iInchiRec = 0; iInchiRec < INCHI_NUM; iInchiRec ++ ) { /* Disconnected/Connected */
        for ( iMobileH = 0; iMobileH < TAUT_NUM; iMobileH ++ ) { /* Mobile/Fixed H */
            cur_num_comp = pOneInput->nNumComponents[iInchiRec][iMobileH];
            if ( !cur_num_comp ) {
                continue;
            }
            /* allocate memory for all existing components */
            pStruct[iInchiRec][iMobileH] = (StrFromINChI *)inchi_calloc( cur_num_comp, sizeof(pStruct[0][0][0]));
            if ( !pStruct[iInchiRec][iMobileH] ) {
                ret = RI_ERR_ALLOC;
                goto exit_error;
            }
            /* set conversion mode */
            bCurI2A_Flag = (iMobileH? 0: I2A_FLAG_FIXEDH) | (iInchiRec? I2A_FLAG_RECMET : 0);
            if ( iMobileH ) {
                ip->nMode &= ~REQ_MODE_BASIC;
            } else {
                ip->nMode |= REQ_MODE_BASIC;
            }
            /* InChI --> structure conversion for all components except duplicated */
            for ( k = 0; k < cur_num_comp; k ++ ) { /* components */
                if ( !iMobileH && !pOneInput->pInpInChI[iInchiRec][iMobileH][k].nNumberOfAtoms ||
                     pOneInput->pInpInChI[iInchiRec][iMobileH][k].bDeleted ||
                     pOneInput->pInpInChI[iInchiRec][iMobileH][k].nLink < 0 ) {

                    pStruct[iInchiRec][iMobileH][k].nLink = pOneInput->pInpInChI[iInchiRec][iMobileH][k].nLink;
                    pStruct[iInchiRec][iMobileH][k].bDeleted = pOneInput->pInpInChI[iInchiRec][iMobileH][k].bDeleted;
                    continue; /* do not create a structure out of an unavailable
                                 Fixed-H InChI or out of the one present in Reconnected layer */
#ifdef NEVER  /* a wrong attempt to process deleted components here */
                    if ( pStruct[iInchiRec][iMobileH][k].nLink = pOneInput->pInpInChI[iInchiRec][iMobileH][k].nLink ) {
                        continue; /* do not create a structure out of an unavailable
                                     Fixed-H InChI or out of the one present in Reconnected layer */
                    } else
                    if ( iMobileH && pOneInput->pInpInChI[iInchiRec][iMobileH][k].nNumberOfAtoms &&
                         pOneInput->pInpInChI[iInchiRec][iMobileH][k].bDeleted &&
                         pOneInput->pInpInChI[iInchiRec][iMobileH][0].bDeleted ) {
                        /* all components are protons */
                        ;
                    } else {
                        continue;
                    }
#endif
                }
                if ( bHasSomeFixedH && iMobileH && k < pOneInput->nNumComponents[iInchiRec][TAUT_NON] &&
                     pOneInput->pInpInChI[iInchiRec][TAUT_NON][k].nNumberOfAtoms ) {
                    continue; /* do not process Mobile-H if Fixed-H is requested and exists */
                }
                pStruct[iInchiRec][iMobileH][k].pSrm      = pSrm;
                pStruct[iInchiRec][iMobileH][k].iInchiRec = iInchiRec;
                pStruct[iInchiRec][iMobileH][k].iMobileH  = iMobileH;

                /****************************************************/
                /*                                                  */
                /* Convert InChI of one component into a Structure  */
                /*                                                  */
                /****************************************************/

                ret = InChI2Atom( ip, sd, szCurHdr, num_inp, pStruct[iInchiRec][iMobileH]+k, k,
                                   0 /* AtNoOffset*/, bCurI2A_Flag, bHasSomeFixedH, pOneInput );
                pStruct[iInchiRec][iMobileH][k].nLink = pOneInput->pInpInChI[iInchiRec][iMobileH][k].nLink;
                if ( ret < 0 ) {
#if ( bRELEASE_VERSION != 1 )
#ifndef TARGET_API_LIB
                    /* !!! Conversion Error -- Ignore for now !!! */
                    fprintf( stdout, "%ld %s Conversion failed: %d, %c%c comp %d\n",
                        num_inp, szCurHdr? szCurHdr : "Struct", ret, iInchiRec? 'R':'D', iMobileH? 'M':'F', k+1); 
#endif
#endif
                    if ( ret == CT_USER_QUIT_ERR ) {
                        goto exit_error;
                    }
                    pStruct[iInchiRec][iMobileH][k].nError = ret;
                    ret = 0; /* force to ignore the errors for now !!!! */
                    num_err ++;
                }
            }
        }
    }
exit_error:
    ulProcessingTime += InchiTimeElapsed( &ulTStart );
    sd->ulStructTime += ulProcessingTime;
    return ret<0? ret : num_err;
}
/**************************************************************************************/
int AddProtonAndIsoHBalanceToMobHStruct( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd,
                                         long num_inp, int bHasSomeFixedH, char *szCurHdr,
                             StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM], InpInChI *pOneInput)
{
    COMPONENT_REM_PROTONS nToBeRemovedByNormFromRevrs[INCHI_NUM];
    int                   nRemovedByNormFromRevrs[INCHI_NUM];
    int                   nRemovedByRevrs[INCHI_NUM];

    int   nDeltaFromDisconnected = 0, nRemovedProtonsByNormFromRevrs, nRemovedProtonsByRevrs, num_changes = 0;
    NUM_H nIsoDeltaFromDisconnected[NUM_H_ISOTOPES];
    int iInchiRec, i, k, k1, ret = 0;
    int  nChargeInChI, nChargeRevrs;

    if ( bHasSomeFixedH ) {
        return 0; /* 2005-03-01 */
    }

    /* num protons removed by InChI Normalization from the original structure */
    for ( i = 0; i < INCHI_NUM; i ++ ) {
        nToBeRemovedByNormFromRevrs[i].nNumRemovedProtons = pOneInput->nNumProtons[i][TAUT_YES].nNumRemovedProtons;
        for ( k = 0; k < NUM_H_ISOTOPES; k ++ ) {
            nToBeRemovedByNormFromRevrs[i].nNumRemovedIsotopicH[k] = pOneInput->nNumProtons[i][TAUT_YES].nNumRemovedIsotopicH[k];
        }
    }
    /* accumulate here num. protons removed by the normalization from the reversed structure */
    nRemovedByNormFromRevrs[INCHI_BAS] =
    nRemovedByNormFromRevrs[INCHI_REC] = 0;
    nRemovedByRevrs[INCHI_REC] =
    nRemovedByRevrs[INCHI_BAS] = 0;
    /* protons added/removed by InChI Normalization to/from Restored Structure might have been added by StructureRestore */
    for ( iInchiRec = 0; iInchiRec < INCHI_NUM; iInchiRec ++ ) {
        for ( k = 0; k < pOneInput->nNumComponents[iInchiRec][TAUT_YES]; k ++ ) {
            if ( !bInpInchiComponentExists( pOneInput, iInchiRec, TAUT_YES, k ) ) {
                continue;
            }
            nRemovedProtonsByNormFromRevrs = 0; /* Num protons removed from the Restored Structure by InChI Normalization */
            nRemovedProtonsByRevrs         = 0; /* Num protons removed by the Reconstruction from the Restored Structure */
            if ( iInchiRec == INCHI_REC || iInchiRec == INCHI_BAS && (k1=pStruct[iInchiRec][TAUT_YES][k].nLink) >= 0 ) {
                
                REV_INCHI  *pRevInChI   = &pStruct[iInchiRec][TAUT_YES][k].RevInChI;
                INChI_Aux  **pINChI_Aux2 = pRevInChI->pINChI_Aux[iInchiRec][0]; /* component 0*/
                INChI      **pINChI_Revr = pRevInChI->pINChI[iInchiRec][0];
                INChI       *pINChI_Orig = pOneInput->pInpInChI[iInchiRec][TAUT_YES]+k;
                nChargeRevrs = pINChI_Revr? pINChI_Revr[TAUT_YES]->nTotalCharge : NO_VALUE_INT;
                nChargeInChI = pINChI_Orig->nTotalCharge;
                if ( pINChI_Aux2 ) {
                    nRemovedProtonsByNormFromRevrs = pINChI_Aux2[TAUT_YES]->nNumRemovedProtons;
                }
                nRemovedProtonsByRevrs = pStruct[iInchiRec][TAUT_YES][k].nNumRemovedProtonsByRevrs;
                pStruct[iInchiRec][TAUT_YES][k].nChargeRevrs = nChargeRevrs;
                pStruct[iInchiRec][TAUT_YES][k].nChargeInChI = nChargeInChI;
            } else
            if ( 0 <= ( k1 = -(1+pStruct[iInchiRec][TAUT_YES][k].nLink) ) ) {
                REV_INCHI  *pRevInChI   = &pStruct[INCHI_REC][TAUT_YES][k1].RevInChI;
                INChI_Aux  **pINChI_Aux2 = pRevInChI->pINChI_Aux[INCHI_BAS][0]; /* component 0 */
                INChI      **pINChI_Revr = pRevInChI->pINChI[INCHI_BAS][0];
                INChI       *pINChI_Orig = pOneInput->pInpInChI[INCHI_REC][TAUT_YES]+k1;
                nChargeRevrs = pINChI_Revr? pINChI_Revr[TAUT_YES]->nTotalCharge : NO_VALUE_INT;
                nChargeInChI = pINChI_Orig->nTotalCharge;
                if ( pINChI_Aux2 ) {
                    nRemovedProtonsByNormFromRevrs = pINChI_Aux2[TAUT_YES]->nNumRemovedProtons;
                }
                /* this component cannot be disconnected because it is same as in reconnected layer */
                nRemovedProtonsByRevrs = pStruct[INCHI_REC][TAUT_YES][k1].nNumRemovedProtonsByRevrs;
                pStruct[iInchiRec][TAUT_YES][k1].nChargeRevrs = nChargeRevrs;
                pStruct[iInchiRec][TAUT_YES][k1].nChargeInChI = nChargeInChI;
            }
            /* how many protons (to be removed by InChI Normalization) to add = 
               (proton balance in InChI} - 
               {number of protons known to be removed by InChI Normalization from Reconstructed structure} */
            nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedProtons -= nRemovedProtonsByNormFromRevrs;
            nRemovedByNormFromRevrs[iInchiRec] += nRemovedProtonsByNormFromRevrs;
            nRemovedByRevrs[iInchiRec]         += nRemovedProtonsByRevrs;
            pStruct[iInchiRec][TAUT_YES][k].nRemovedProtonsByNormFromRevrs = nRemovedProtonsByNormFromRevrs;
        }
    }

    /* Since fixed-H layer is missing we need to add proton balance to the components */
    memset( nIsoDeltaFromDisconnected, 0, sizeof(nIsoDeltaFromDisconnected) );
    for ( iInchiRec = INCHI_REC; INCHI_BAS <= iInchiRec; iInchiRec -- ) {
        /*
        if ( !pOneInput->nNumComponents[iInchiRec][TAUT_NON] &&
              pOneInput->nNumComponents[iInchiRec][TAUT_YES] ) {
        */
             int bHasRecMobH = (iInchiRec==INCHI_BAS && pOneInput->nNumComponents[INCHI_REC][TAUT_YES]);
             /* bHasRecMobH means all components that could not be disconnected are in reconnected part */
             if ( iInchiRec==INCHI_BAS ) {
                 /* second pass: common structures have been changed */
                 nToBeRemovedByNormFromRevrs[INCHI_BAS].nNumRemovedProtons += nDeltaFromDisconnected;
             }
             /* after proton removal InChI is recalculated */

             ret = AddRemProtonsInRestrStruct( ip, sd, num_inp, bHasSomeFixedH, pStruct[iInchiRec][TAUT_YES],
                                     pOneInput->nNumComponents[iInchiRec][TAUT_YES],
                                     bHasRecMobH? pStruct[INCHI_REC][TAUT_YES] : NULL,
                                     bHasRecMobH? pOneInput->nNumComponents[INCHI_REC][TAUT_YES]:0,
                                     &nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedProtons,
                                     (iInchiRec==INCHI_REC)?&nDeltaFromDisconnected : NULL);
             if ( ret < 0 ) {
                 goto exit_function;
             }
             num_changes += ret;
        /*
        }
        */
    }
    /* if fixed-H layer is missing then we need to add isotopic exchangeable proton balance to the components */
    for ( iInchiRec = INCHI_REC; INCHI_BAS <= iInchiRec; iInchiRec -- ) {
        /*
        if ( !pOneInput->nNumComponents[iInchiRec][TAUT_NON] &&
              pOneInput->nNumComponents[iInchiRec][TAUT_YES] ) {
        */
             int bHasRecMobH = (iInchiRec==INCHI_BAS && pOneInput->nNumComponents[INCHI_REC][TAUT_YES]);
             /* bHasRecMobH means all components that could not be disconnected are in reconnected part */
             if ( iInchiRec==INCHI_BAS ) {
                 /* second pass: common structures have been changed */
                 for ( k = 0; k < NUM_H_ISOTOPES; k ++ ) {
                     nToBeRemovedByNormFromRevrs[INCHI_BAS].nNumRemovedIsotopicH[k] += nIsoDeltaFromDisconnected[k];
                 }
             }
             /* after proton removal InChI is recalculated */
             ret = AddRemIsoProtonsInRestrStruct( ip, sd, num_inp, bHasSomeFixedH, pStruct[iInchiRec][TAUT_YES],
                                     pOneInput->nNumComponents[iInchiRec][TAUT_YES],
                                     bHasRecMobH? pStruct[INCHI_REC][TAUT_YES] : NULL,
                                     bHasRecMobH? pOneInput->nNumComponents[INCHI_REC][TAUT_YES]:0,
                                     nToBeRemovedByNormFromRevrs[iInchiRec].nNumRemovedIsotopicH,
                                     (iInchiRec==INCHI_REC)?nIsoDeltaFromDisconnected : NULL);
             if ( ret < 0 ) {
                 goto exit_function;
             }
             num_changes += ret;
        /*
        }
        */
    }

exit_function:
    return ret;
}

/*************************************************************/
void FreeStrFromINChI( StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM], int nNumComponents[INCHI_NUM][TAUT_NUM] )
{
    int iInchiRec, iMobileH, cur_num_comp, k, j;
    StrFromINChI *pStruct1;
    for ( iInchiRec = 0; iInchiRec < INCHI_NUM; iInchiRec ++ ) {
        for ( iMobileH = 0; iMobileH < TAUT_NUM; iMobileH ++ ) {
            cur_num_comp = nNumComponents[iInchiRec][iMobileH];
            if ( !cur_num_comp || !(pStruct1=pStruct[iInchiRec][iMobileH]) ) {
                continue;
            }
            for ( k = 0; k < cur_num_comp; k ++ ) {
                if ( pStruct1[k].at ) {
                    inchi_free(pStruct1[k].at);
                }
                if ( pStruct1[k].at2 ) {
                    inchi_free(pStruct1[k].at2);
                }
                if ( pStruct1[k].st ) {
                    inchi_free(pStruct1[k].st);
                }
                if ( pStruct1[k].pVA ) {
                    inchi_free(pStruct1[k].pVA);
                }
                /*
                if ( pStruct1[k].ti.t_group ) {
                    inchi_free( pStruct1[k].ti.t_group );
                }
                */
                if ( pStruct1[k].pXYZ ) {
                    inchi_free(pStruct1[k].pXYZ);
                }
                /*==== begin ====*/
                free_t_group_info( &pStruct1[k].ti );
                if ( pStruct1[k].endpoint ) {
                    inchi_free(pStruct1[k].endpoint);
                }
                if ( pStruct1[k].fixed_H ) {
                    inchi_free(pStruct1[k].fixed_H);
                }
                for ( j = 0; j < TAUT_NUM; j ++ ) {
                    if ( pStruct1[k].nAtno2Canon[j] )
                        inchi_free( pStruct1[k].nAtno2Canon[j] );
                    if ( pStruct1[k].nCanon2Atno[j] )
                        inchi_free( pStruct1[k].nCanon2Atno[j] );
                }
                /*===== end ======*/
                /*  free INChI memory */
                FreeAllINChIArrays( pStruct1[k].RevInChI.pINChI,
                                    pStruct1[k].RevInChI.pINChI_Aux,
                                    pStruct1[k].RevInChI.num_components );
#ifdef NEVER
                /* don't do that: these are just pointers to OneInput structure members */
                Free_INChI( &pStruct1[k].pINChI );
                Free_INChI_Aux( &pStruct1[k].pINChI_Aux );
                if ( pStruct1[k].inp_norm_data ) {
                    FreeInpAtomData( pStruct1[k].inp_norm_data );
                    inchi_free( pStruct1[k].inp_norm_data );
                }
#endif
            }
            inchi_free(pStruct[iInchiRec][iMobileH]);
            pStruct[iInchiRec][iMobileH] = NULL;
        }
    }
}
/********************************************************************/
void FreeInpInChI( InpInChI *pOneInput )
{
    int iINChI, k, j;
    for ( iINChI = 0; iINChI < INCHI_NUM; iINChI ++ ) {
        for ( j = 0; j < TAUT_NUM; j ++ ) {
            if ( pOneInput->pInpInChI[iINChI][j] ) {
                for ( k = 0; k < pOneInput->nNumComponents[iINChI][j]; k ++ ) {
                    Free_INChI_Members( &pOneInput->pInpInChI[iINChI][j][k] );
                }
                inchi_free(pOneInput->pInpInChI[iINChI][j]);
                pOneInput->pInpInChI[iINChI][j] = NULL;
            }
            if ( pOneInput->nNumProtons[iINChI][j].pNumProtons ) {
                inchi_free( pOneInput->nNumProtons[iINChI][j].pNumProtons );
                pOneInput->nNumProtons[iINChI][j].pNumProtons = NULL;
            }
        }
    }
    if ( pOneInput->atom ) inchi_free(pOneInput->atom);
    memset( pOneInput, 0, sizeof(*pOneInput) );
}

/***********************************************************************************************/
int CompareAllOrigInchiToRevInChI(StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM], InpInChI *pOneInput, int bReqNonTaut,
                                  long num_inp, char *szCurHdr)
{
    int i, iInchiRec, iMobileH, iMobileHpStruct, num_components, iComponent, ret=0;
    COMPONENT_REM_PROTONS nCurRemovedProtons, nNumRemovedProtons;
    INChI     *pInChI[TAUT_NUM];
    INCHI_MODE  CompareInchiFlags[TAUT_NUM];
    memset( pOneInput->CompareInchiFlags[0], 0, sizeof(pOneInput->CompareInchiFlags[0]) );
    memset( &nNumRemovedProtons, 0, sizeof(nNumRemovedProtons) );
    
    /* do we have reconnected InChI ?*/
    iInchiRec = INCHI_REC;
    iMobileH  = TAUT_NON;
    if ( !pOneInput->nNumComponents[iInchiRec][TAUT_YES] && !pOneInput->nNumComponents[iInchiRec][TAUT_NON] ) {
        iInchiRec = INCHI_BAS;
    }
    /* do we have Mobile or Fixed-H ? */
    if ( !pOneInput->nNumComponents[iInchiRec][TAUT_NON] || !bReqNonTaut ) {
        iMobileH = TAUT_YES;  /* index for pOneInput */
    }
    /* if a restored structure has Fixed-H InChI then its mobile-H restored InChI is in Fixed-H pStruct */
    num_components = pOneInput->nNumComponents[iInchiRec][iMobileH];
    for ( iComponent = 0; iComponent < num_components; iComponent ++ ) {
        int bMobileH = iMobileH;
        pInChI[0]     = pInChI[1]     = NULL;
        if ( pOneInput->pInpInChI[iInchiRec][bMobileH][iComponent].nNumberOfAtoms &&
             !pOneInput->pInpInChI[iInchiRec][bMobileH][iComponent].bDeleted ) {
            /* the requested InChI layer exists */
            pInChI[0]     = &pOneInput->pInpInChI[iInchiRec][bMobileH][iComponent];
            if ( bMobileH == TAUT_NON ) {
                pInChI[1] = &pOneInput->pInpInChI[iInchiRec][TAUT_YES][iComponent];
            }
        } else
        if ( bMobileH == TAUT_NON &&
             pOneInput->pInpInChI[iInchiRec][TAUT_YES][iComponent].nNumberOfAtoms &&
             !pOneInput->pInpInChI[iInchiRec][TAUT_YES][iComponent].bDeleted ) {
            /* the requested Fixed-H InChI layer does not exist; however, the Mobile-H does exist */
            bMobileH = TAUT_YES; /* only Mobile-H is available */
            pInChI[0] = &pOneInput->pInpInChI[iInchiRec][bMobileH][iComponent];
        }
        memset( CompareInchiFlags, 0, sizeof(CompareInchiFlags) );
        memset( &nCurRemovedProtons, 0, sizeof(nCurRemovedProtons) );
        iMobileHpStruct = 
#if ( bRELEASE_VERSION == 0 )
#ifndef TARGET_API_LIB
        /* legacy: reproduce old output */
        OldPrintCompareOneOrigInchiToRevInChI(pStruct[iInchiRec][bMobileH]+iComponent, pInChI, bMobileH,
                                              iComponent, num_inp, szCurHdr);
#endif
#endif
        /* one component comparison result bits */
        ret = CompareOneOrigInchiToRevInChI( pStruct[iInchiRec][bMobileH]+iComponent, pInChI, bMobileH, iComponent,
                                             num_inp, szCurHdr, &nCurRemovedProtons, CompareInchiFlags);
        if ( ret >= 0 ) {
            /* no errors encountered -> accumulate removed protons from individual Mobile-H layers of components */
            nNumRemovedProtons.nNumRemovedProtons += nCurRemovedProtons.nNumRemovedProtons;
            for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) {
                nNumRemovedProtons.nNumRemovedIsotopicH[i] += nCurRemovedProtons.nNumRemovedIsotopicH[i];
            }
            /* accumulate compare bits */
            for ( i = 0; i < TAUT_NUM; i ++ ) {
                pOneInput->CompareInchiFlags[0][i] |= CompareInchiFlags[i];
            }
        } else {
            goto exit_function;
        }
    }
    if ( iMobileH == TAUT_YES ) {
        if ( pOneInput->nNumProtons[iInchiRec][iMobileH].pNumProtons ) {
            ret = RI_ERR_PROGR; /* in Mobile-H case proton balances are split between compoments */
        } else {
            /*   num removed protons in orig. InChI      num removed protons in restored InChi */
            if ( nNumRemovedProtons.nNumRemovedProtons != pOneInput->nNumProtons[iInchiRec][iMobileH].nNumRemovedProtons ) {
                /* restored structure InChI has less or more removed protons */
                pOneInput->CompareInchiFlags[0][TAUT_YES] |= INCHIDIFF_MOBH_PROTONS;
#if ( bRELEASE_VERSION == 0 )
                /* debug output only */
                {
                    int num_H_AddedByRevrs = pOneInput->nNumProtons[iInchiRec][iMobileH].nNumRemovedProtons
                                              - nNumRemovedProtons.nNumRemovedProtons;
                    fprintf( stdout, "COMPARE_INCHI: %ld: %s %cM: Proton balance (Diff: %d, RevrsRem=%d)\n",
                        num_inp, szCurHdr? szCurHdr : "Struct", iInchiRec? 'R':'D',
                        pOneInput->nNumProtons[iInchiRec][iMobileH].nNumRemovedProtons,num_H_AddedByRevrs);
                }
#endif
            }
            for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) {
                if ( nNumRemovedProtons.nNumRemovedIsotopicH[i] != pOneInput->nNumProtons[iInchiRec][TAUT_YES].nNumRemovedIsotopicH[i] ) {
                    pOneInput->CompareInchiFlags[0][TAUT_YES] |= INCHIDIFF_MOB_ISO_H;
#if ( bRELEASE_VERSION == 0 )
                    /* debug output only */
                    {
                    int num_H_AddedByRevrs = pOneInput->nNumProtons[iInchiRec][TAUT_YES].nNumRemovedIsotopicH[i]
                      - nNumRemovedProtons.nNumRemovedIsotopicH[i];
                    fprintf( stdout, "COMPARE_INCHI: %ld: %s %cM: Iso Xchg %dH balance (Diff: %d, RevrsRem=%d)\n",
                        num_inp, szCurHdr? szCurHdr : "Struct", iInchiRec? 'R':'D', i+1,
                        pOneInput->nNumProtons[iInchiRec][TAUT_YES].nNumRemovedIsotopicH[i],num_H_AddedByRevrs);
                    }
#endif
                }
            }
        }
    }

exit_function:
    return ret;
}
/***********************************************************************************************/
int CompareAllDisconnectedOrigInchiToRevInChI(StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM],
                                              InpInChI *pOneInput, int bHasSomeFixedH,
                                              long num_inp, char *szCurHdr)
{
    int i, k, m, n, iInChI, iMobileH, bMobileH, ifk;
    int num_components_D, num_components_R;
    int nNumCompHaveSeparateProtons_D, nNumCompHaveSeparateProtons_R;
    int num_fragments_D, num_fragments_R, num_fragments_DR, num_fragments, iComponent, ret;
    int ifInChI, ifMobileH, bfMobileH, nLink;
    COMPONENT_REM_PROTONS nNumRemovedProtons_D;     /* removed from the disconnected layer of the Input InChI */
    COMPONENT_REM_PROTONS nNumRemovedProtons_D_all; /* if only totals are avalable */
    COMPONENT_REM_PROTONS nNumRemovedProtons_R; /* removed from disconnected layer of the reconstructed struct */
    COMPONENT_REM_PROTONS nNumRemovedProtons_R_all;
    INCHI_MODE  CompareInchiFlags[TAUT_NUM];
    StrFromINChI *pStruct1;
    INChI_Aux *pINChI_Aux;
    INCHI_SORT *pINChISort1 = NULL; /* from reversed structure */
    INCHI_SORT *pINChISort2 = NULL; /* original input InChI */
    int        nNumNonTaut1=0, nNumNonTaut2=0;

    ret = 0;
    memset( pOneInput->CompareInchiFlags[1], 0, sizeof(pOneInput->CompareInchiFlags[1]) );

    /* count components that are not subject to disconnection */
    if ( !pOneInput->nNumComponents[INCHI_REC][TAUT_YES] &&
         !pOneInput->nNumComponents[INCHI_REC][TAUT_NON] ) {
        return 0; /* nothing to do */
    }

    memset( &nNumRemovedProtons_D, 0, sizeof(nNumRemovedProtons_D) );
    memset( &nNumRemovedProtons_R, 0, sizeof(nNumRemovedProtons_R) );
    memset( &nNumRemovedProtons_D_all, 0, sizeof(nNumRemovedProtons_D_all) );
    memset( &nNumRemovedProtons_R_all, 0, sizeof(nNumRemovedProtons_R_all) );
    memset( CompareInchiFlags, 0, sizeof(CompareInchiFlags) );

    num_components_D = inchi_max( pOneInput->nNumComponents[INCHI_BAS][TAUT_YES],
                                  pOneInput->nNumComponents[INCHI_BAS][TAUT_NON] );
    num_components_R = inchi_max( pOneInput->nNumComponents[INCHI_REC][TAUT_YES],
                                  pOneInput->nNumComponents[INCHI_REC][TAUT_NON] );
    /***********************************************************************************************/
    /* InpInChI: count fragments -- disconnected components that do not match reconnected          */
    /* Accumulate removed H and isotopic H from ALL Fixed-H disconnected components except deleted */
    /* This segment collects info from the original InChI                                          */
    /***********************************************************************************************/
    /*---- Original InChI ----*/
    num_fragments_D = 0;
    iInChI   = INCHI_BAS;
    iMobileH = bHasSomeFixedH? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
    nNumCompHaveSeparateProtons_D = 0;

    /* in case of Mobile-H components here are the proton totals from the original InChI disconn. layer */
    nNumRemovedProtons_D.nNumRemovedProtons = pOneInput->nNumProtons[iInChI][TAUT_YES].nNumRemovedProtons;
    memcpy( nNumRemovedProtons_D.nNumRemovedIsotopicH, 
            pOneInput->nNumProtons[iInChI][TAUT_YES].nNumRemovedIsotopicH,
            sizeof(nNumRemovedProtons_D.nNumRemovedIsotopicH) ); /* total for the disconnected layer */

    for ( k = 0; k < num_components_D; k ++ ) {
        bMobileH = iMobileH;
        if ( !bInpInchiComponentExists( pOneInput, iInChI, bMobileH, k ) ) {
            if ( bInpInchiComponentExists( pOneInput, iInChI, TAUT_YES, k ) ) {
                bMobileH = TAUT_YES;
            } else {
                continue; /* component is missing ??? */
            }
        }
        if ( 0 > (nLink = pOneInput->pInpInChI[iInChI][bMobileH][k].nLink) ) {
            /* component in Disconnected layer is linked to the identical one in the Reconnected layer */
            if ( pOneInput->nNumProtons[INCHI_REC][TAUT_YES].pNumProtons ) {
                nNumCompHaveSeparateProtons_D ++;
                nLink = -(1+nLink);
                nNumRemovedProtons_D.nNumRemovedProtons += pOneInput->nNumProtons[INCHI_REC][TAUT_YES].pNumProtons[nLink].nNumRemovedProtons;
                for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                    nNumRemovedProtons_D.nNumRemovedIsotopicH[m] += pOneInput->nNumProtons[INCHI_REC][TAUT_YES].pNumProtons[nLink].nNumRemovedIsotopicH[m];
                }
            }
            continue; /* same as reconnected */
        }
        /* component in the reconnected layer that was disconnected */
        nNumNonTaut2 += (bMobileH == TAUT_NON);
        if ( pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons ) {
            nNumCompHaveSeparateProtons_D ++;
            nNumRemovedProtons_D.nNumRemovedProtons += pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons[k].nNumRemovedProtons;
            for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                nNumRemovedProtons_D.nNumRemovedIsotopicH[m] += pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons[k].nNumRemovedIsotopicH[m];
            }
        }
        num_fragments_D ++; /* number of disconnected fragments from original reconnected structure */
    }
    /* in case of Mobile-H components here are the proton totals from the original InChI */
    /*
    nNumRemovedProtons_D_all.nNumRemovedProtons = pOneInput->nNumProtons[iInChI][TAUT_YES].nNumRemovedProtons;
    memcpy( nNumRemovedProtons_D_all.nNumRemovedIsotopicH, 
            pOneInput->nNumProtons[iInChI][TAUT_YES].nNumRemovedIsotopicH,
            sizeof(nNumRemovedProtons_D_all.nNumRemovedIsotopicH) );

    */
    /****************************************************************************************************/
    /* count fragments in reconstructed reconnected structure                                           */
    /* accumulate removed H and isotopic H from ALL reconstructed reconnected components except deleted */
    /* This segment collects info from the reconstructed structure InChI                                */
    /****************************************************************************************************/
    /*---- InChI from the reconstructed reconnected structure ----*/
    num_fragments_R = 0;
    iInChI   = INCHI_REC;
    iMobileH = bHasSomeFixedH? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
    nNumCompHaveSeparateProtons_R = 0;
    for ( k = 0; k < num_components_R; k ++ ) {
        bMobileH = iMobileH;
        if ( !bInpInchiComponentExists( pOneInput, iInChI, bMobileH, k ) ) {
            if ( bInpInchiComponentExists( pOneInput, iInChI, TAUT_YES, k ) ) {
                bMobileH = TAUT_YES;
            } else {
                continue; /* component is missing ??? (Deleted proton in Mobile-H layer) */
            }
        }
        if ( 0 < pOneInput->pInpInChI[iInChI][bMobileH][k].nLink ) {
            /* this reconstructed reconnected component was NOT DISCONNECTED */
            /* same component is in the disconnected layer, it has no metal atoms or is an isolated metal atom */
            pStruct1 = pStruct[iInChI][bMobileH]+k;
            ifMobileH = TAUT_YES;  /* Mobile-H Aux_Info contains number removed protons */
            ifInChI   = INCHI_BAS; /* this component cannot be reconnected */
            ifk       = 0;         /* 0th component since it is InChI of a single component */
            /* The statement in the following line is *WRONG*, component number mixed with bMobileH:  */
            /* in RevInchi, when only Mobile-H is present then its only non-NULL InChI has index 0==TAUT_NON */
            if ( bRevInchiComponentExists( pStruct1, ifInChI, ifMobileH, ifk ) ) {
                /* count protons */
                pINChI_Aux = pStruct1->RevInChI.pINChI_Aux[ifInChI][ifk][ifMobileH];
                if ( pINChI_Aux ) {
                    nNumRemovedProtons_R.nNumRemovedProtons += pINChI_Aux->nNumRemovedProtons;
                    for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                        nNumRemovedProtons_R.nNumRemovedIsotopicH[m] += pINChI_Aux->nNumRemovedIsotopicH[m];
                    }
                }
            }
            nNumCompHaveSeparateProtons_R += bRevInchiComponentExists( pStruct1, ifInChI, ALT_TAUT(ifMobileH), ifk );
            continue; /* same as disconnected, has no metal atoms */
        }
        /* this reconstructed reconnected component WAS DISCONNECTED; check its fragments */
        /* it does not have same component in the disconnected layer */
        pStruct1 = pStruct[iInChI][bMobileH]+k;
        num_fragments = pStruct1->RevInChI.num_components[INCHI_BAS];
        ifInChI = INCHI_BAS; /* disconnected layer */
        ifMobileH = bHasSomeFixedH? TAUT_NON : TAUT_YES;
        for ( ifk = 0; ifk < num_fragments; ifk ++ ) {
            bfMobileH = ifMobileH;
            if ( !bRevInchiComponentExists( pStruct1, ifInChI, bfMobileH, ifk ) ) {
                if ( bRevInchiComponentExists( pStruct1, ifInChI, TAUT_YES, ifk ) ) {
                    bfMobileH = TAUT_YES;
                } else {
                    continue; /* fragment does not exist ??? */
                }
            }
            nNumNonTaut1           += (bfMobileH == TAUT_NON);
            nNumCompHaveSeparateProtons_R += (bfMobileH == TAUT_NON);
            /* count protons from fragments made by metal disconnection */
            pINChI_Aux = pStruct1->RevInChI.pINChI_Aux[ifInChI][ifk][TAUT_YES];
            if ( pINChI_Aux ) {
                nNumRemovedProtons_R.nNumRemovedProtons += pINChI_Aux->nNumRemovedProtons;
                for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                    nNumRemovedProtons_R.nNumRemovedIsotopicH[m] += pINChI_Aux->nNumRemovedIsotopicH[m];
                }
            }
            num_fragments_R ++; /* number of disconnected fragments from reconstructed reconnected structure */
        }
    }
    /*---------------- special treatment of the last reconstructed component -----------------*/
    /*---------------- this may contain separate protons added by the reconstruction ---------*/
    k = num_components_R - 1;
    pStruct1 = pStruct[iInChI][iMobileH]+k;
    if ( iMobileH == TAUT_YES && !bHasSomeFixedH &&
         bInpInchiComponentDeleted( pOneInput, iInChI, iMobileH, k ) &&
         (num_fragments = pStruct1->RevInChI.num_components[INCHI_BAS]) ) {

        ifInChI = INCHI_BAS; /* disconnected layer */
        ifMobileH = TAUT_YES;
        for ( ifk = 0; ifk < num_fragments; ifk ++ ) {
            bfMobileH = ifMobileH;
            if ( !bRevInchiComponentDeleted( pStruct1, ifInChI, bfMobileH, ifk ) ) {
                continue; /* fragment does exist ??? Should not happen */
            }
            /*
            nNumNonTaut1           += (bfMobileH == TAUT_NON);
            nNumCompHaveSeparateProtons_R += (bfMobileH == TAUT_NON);
            */
            /* count protons from fragments made by metal disconnection */
            pINChI_Aux = pStruct1->RevInChI.pINChI_Aux[ifInChI][ifk][TAUT_YES];
            if ( pINChI_Aux ) {
                nNumRemovedProtons_R.nNumRemovedProtons += pINChI_Aux->nNumRemovedProtons;
                for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
                    nNumRemovedProtons_R.nNumRemovedIsotopicH[m] += pINChI_Aux->nNumRemovedIsotopicH[m];
                }
            }
            /*num_fragments_R ++;*/ /* number of disconnected fragments from reconstructed reconnected structure */
        }
    }



    num_fragments_DR = inchi_max( num_fragments_D, num_fragments_R );
    /* in case of correct reconstruction, num_fragments_D, num_fragments_R */

    if ( !num_fragments_DR ) {
        return 0; /* no component was disconnected */
    }
    if ( num_fragments_D != num_fragments_R ) {
        for ( i = 0; i < TAUT_NUM; i ++ ) {
            if ( pOneInput->nNumComponents[INCHI_BAS][i] ) {
                pOneInput->CompareInchiFlags[1][i] |= INCHIDIFF_PROBLEM;
            }
        }
        return 1; /* severe error */
    }


    pINChISort1 = (INCHI_SORT *)inchi_calloc(num_fragments_DR, sizeof(pINChISort1[0]));
    pINChISort2 = (INCHI_SORT *)inchi_calloc(num_fragments_DR, sizeof(pINChISort2[0]));
    if ( !pINChISort1 || !pINChISort2 ) {
        ret = RI_ERR_ALLOC;
        goto exit_function;
    }

    /* accumulate original InChI of fragments -- disconnected components that do not match reconnected */
    iInChI   = INCHI_BAS;
    iMobileH = bHasSomeFixedH? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
    for ( k = n = 0; k < num_components_D; k ++ ) {
        bMobileH = iMobileH;
        if ( !bInpInchiComponentExists( pOneInput, iInChI, bMobileH, k ) ) {
            if ( bInpInchiComponentExists( pOneInput, iInChI, TAUT_YES, k ) ) {
                bMobileH = TAUT_YES;
            } else {
                continue; /* component is missing ??? (Deleted proton in Mobile-H layer) */
            }
        }
        if ( 0 > pOneInput->pInpInChI[iInChI][bMobileH][k].nLink ) {
            continue; /* same as reconnected */
        }
        /* the component exists in disconnected layer of the orig. InChI only: it is a fragment */
        pINChISort2[n].pINChI[bMobileH] = pOneInput->pInpInChI[iInChI][bMobileH] + k;
        if ( bMobileH == TAUT_NON && 
             (bInpInchiComponentExists( pOneInput, iInChI, TAUT_YES, k ) ||
              bInpInchiComponentDeleted( pOneInput, iInChI, TAUT_YES, k ) ) ) {
            pINChISort2[n].pINChI[TAUT_YES] = pOneInput->pInpInChI[iInChI][TAUT_YES] + k;
        }
        /* the last sort key is a number of removed protons */
        pINChISort2[n].ord_number = pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons?
                     pOneInput->nNumProtons[iInChI][TAUT_YES].pNumProtons[k].nNumRemovedProtons : 0;
        pINChISort2[n].n1 = k;  /* orig. InChI disconnected layer component number */
        pINChISort2[n].n2 = -1; /* no fragment index */
        n ++;     
    }

    /* accumulate fragments from the reconstructed structure */
    iInChI   = INCHI_REC;
    iMobileH = bHasSomeFixedH? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
    for ( k = n = 0; k < num_components_R; k ++ ) {
        bMobileH = iMobileH;
        if ( !bInpInchiComponentExists( pOneInput, iInChI, bMobileH, k ) ) {
            if ( bInpInchiComponentExists( pOneInput, iInChI, TAUT_YES, k ) ) {
                bMobileH = TAUT_YES;
            } else {
                continue; /* component is missing ??? (Deleted proton in Mobile-H layer) */
            }
        }
        /* the reconstructed structure */
        if ( 0 < pOneInput->pInpInChI[iInChI][bMobileH][k].nLink ) {
            continue; /* same as disconnected, has no metal atoms */
        }
        /* this reconstructed structure was disconnected */
        pStruct1 = pStruct[iInChI][bMobileH]+k;
        num_fragments = pStruct1->RevInChI.num_components[INCHI_BAS];
        ifInChI = INCHI_BAS;
        ifMobileH = bHasSomeFixedH? TAUT_NON : TAUT_YES;
        for ( i = 0; i < num_fragments; i ++ ) {
            bfMobileH = ifMobileH;
            if ( !bRevInchiComponentExists( pStruct1, ifInChI, bfMobileH, i ) ) {
                if ( bRevInchiComponentExists( pStruct1, ifInChI, TAUT_YES, i ) ) {
                    bfMobileH = TAUT_YES;
                } else {
                    continue; /* component is missing ??? */
                }
            }
            pINChISort1[n].pINChI[bfMobileH] =  pStruct1->RevInChI.pINChI[ifInChI][i][bfMobileH];
            if ( bfMobileH == TAUT_NON /*&& bRevInchiComponentExists( pStruct1, ifInChI, TAUT_YES, i )*/ ) {
                pINChISort1[n].pINChI[TAUT_YES] = pStruct1->RevInChI.pINChI[ifInChI][i][TAUT_YES];
                /* remove Fixed-H InChI if is is identical to Mobile-H */
                /* do it exactly same way the identical components were removed from InpInChI */
                if ( !CompareReversedINChI( pINChISort1[n].pINChI[bfMobileH],
                                            pINChISort1[n].pINChI[TAUT_YES], NULL, NULL ) ) {
                    pINChISort1[n].pINChI[bfMobileH] = NULL; /* remove Fixed-H layer */
                } else {
                    pINChISort1[n].ord_number = pStruct1->RevInChI.pINChI_Aux[ifInChI][i][TAUT_YES]->nNumRemovedProtons;
                }
            }

            pINChISort1[n].n1 = k;  /* reconstructed reconnected structure component index */
            pINChISort1[n].n2 = i;  /* index of a fragment made out of this component */
            n ++;
        }
    }
    
    /* sort fragment InChI before comparing them */
    qsort( pINChISort1, num_fragments_D, sizeof(pINChISort1[0]), CompINChITaut2 );
    qsort( pINChISort2, num_fragments_R, sizeof(pINChISort2[0]), CompINChITaut2 );
    
    /* compare fragments -- components present in disconnected layer only */
    for ( iComponent = 0; iComponent < num_fragments_DR; iComponent ++ ) {
        INChI *pInChI1[TAUT_NUM]; /* from reversed structure */
        INChI *pInChI2[TAUT_NUM]; /* original input InChI */
        for ( i = 0; i < TAUT_NUM; i ++ ) {
            pInChI1[i] = pINChISort1[iComponent].pINChI[i];
            pInChI2[i] = pINChISort2[iComponent].pINChI[i];
        }
        CompareTwoPairsOfInChI( pInChI1, pInChI2, !bHasSomeFixedH, CompareInchiFlags );
    }
    
    if ( /*nNumNonTaut1 && nNumNonTaut2 &&*/ bHasSomeFixedH ) {
        if ( nNumCompHaveSeparateProtons_D || nNumCompHaveSeparateProtons_R ) {
            /* for each component, compare number removed protons */
            /* comparison does not make sense if Disconnected Fixed-H layer is not present */
            for ( iComponent = 0; iComponent < num_fragments_DR; iComponent ++ ) {
                NUM_H   nNumRemovedIsotopicH1[NUM_H_ISOTOPES];
                NUM_H   nNumRemovedIsotopicH2[NUM_H_ISOTOPES];

                memset( nNumRemovedIsotopicH1, 0, sizeof(nNumRemovedIsotopicH1) );
                memset( nNumRemovedIsotopicH2, 0, sizeof(nNumRemovedIsotopicH2) );
                /* compare removed protons */
                if ( pINChISort1[iComponent].ord_number != pINChISort2[iComponent].ord_number ) {
                    CompareInchiFlags[TAUT_YES] |= INCHIDIFF_MOBH_PROTONS; /* diff number of removed protons */
                }
                /* also compare removed isotopic atoms H */
                k = pINChISort2[iComponent].n1; /* input InChI, OneInput */
                if ( pOneInput->nNumProtons[INCHI_BAS][TAUT_YES].pNumProtons ) {
                    memcpy( nNumRemovedIsotopicH2,
                        pOneInput->nNumProtons[INCHI_BAS][TAUT_YES].pNumProtons[k].nNumRemovedIsotopicH,
                        sizeof( nNumRemovedIsotopicH2 ) );
                }
               /* get fragments of reconstructed structure removed protons info */
                k = pINChISort1[iComponent].n1; /* restored component number */
                i = pINChISort1[iComponent].n2; /* subcomponent number */
                iInChI   = INCHI_REC;
                iMobileH = bHasSomeFixedH? !pOneInput->nNumComponents[iInChI][TAUT_NON] : TAUT_YES;
                bMobileH = iMobileH;
                if ( !bInpInchiComponentExists( pOneInput, iInChI, bMobileH, k ) ) {
                    if ( bInpInchiComponentExists( pOneInput, iInChI, TAUT_YES, k ) ) {
                        bMobileH = TAUT_YES;
                    } else {
                       goto compare_iso_H;
                    }
                }
                if ( pOneInput->pInpInChI[iInChI][bMobileH][k].nLink ) {
                    continue;
                    /*
                    ret = RI_ERR_PROGR;
                    goto exit_function;
                    */
                }
                pStruct1 = pStruct[iInChI][bMobileH]+k;
                num_fragments = pStruct1->RevInChI.num_components[INCHI_BAS];
                ifInChI = INCHI_BAS;
                ifMobileH = bHasSomeFixedH? TAUT_NON : TAUT_YES;
                if ( i < num_fragments ) {
                    bfMobileH = ifMobileH;
                    if ( !bRevInchiComponentExists( pStruct1, ifInChI, bfMobileH, i ) ) {
                        if ( bRevInchiComponentExists( pStruct1, ifInChI, TAUT_YES, i ) ) {
                            bfMobileH = TAUT_YES;
                        } else {
                            goto compare_iso_H;
                        }
                    }
                    memcpy( nNumRemovedIsotopicH1,
                        pStruct1->RevInChI.pINChI_Aux[ifInChI][i][TAUT_YES]->nNumRemovedIsotopicH,
                        sizeof( nNumRemovedIsotopicH1 ) );
                }
compare_iso_H:
                if ( memcmp( nNumRemovedIsotopicH1, nNumRemovedIsotopicH2, sizeof( nNumRemovedIsotopicH1 ) ) ) {
                    CompareInchiFlags[TAUT_YES] |= INCHIDIFF_REM_ISO_H;
                }
            }
        }
    } else
    /*if ( !nNumNonTaut1 && !nNumNonTaut2 || !bHasSomeFixedH )*/ {
        /* compare totals for removed protons and isotopic H */
        if ( pOneInput->nNumProtons[INCHI_BAS][TAUT_YES].nNumRemovedProtons !=
             nNumRemovedProtons_R.nNumRemovedProtons ) {
            CompareInchiFlags[TAUT_YES] |= INCHIDIFF_MOBH_PROTONS;
        }
        if ( memcmp( pOneInput->nNumProtons[INCHI_BAS][TAUT_YES].nNumRemovedIsotopicH,
                     nNumRemovedProtons_R.nNumRemovedIsotopicH,
                     sizeof( nNumRemovedProtons_R.nNumRemovedIsotopicH ) ) ) {
                CompareInchiFlags[TAUT_YES] |= INCHIDIFF_REM_ISO_H;
        }
    }

    if ( !nNumNonTaut1 == !nNumNonTaut2 ) {
        ; /* difference if(nNumNonTaut1 != nNumNonTaut2) will be caught in InChI comparison */
    } else
    if ( nNumNonTaut1 ) {
        /* reconstructed has Fixed-H while the original has not: extra Fixed-H layer */
        CompareInchiFlags[TAUT_YES] |= INCHIDIFF_WRONG_TAUT;
    } else {
        /* the original InChI has Fixed-H while the reconstructed one has not: missing Fixed-H layer */
        CompareInchiFlags[TAUT_YES] |= INCHIDIFF_NO_TAUT;
    }
    for ( i = 0; i < TAUT_NUM; i ++ ) {
        pOneInput->CompareInchiFlags[1][i] |= CompareInchiFlags[i];
    }

    /* compare totals */
    if ( nNumRemovedProtons_R.nNumRemovedProtons != nNumRemovedProtons_D.nNumRemovedProtons ) {
        CompareInchiFlags[TAUT_YES] |= INCHIDIFF_MOBH_PROTONS; /* diff number of removed protons */
    }
    if ( memcmp( nNumRemovedProtons_R.nNumRemovedIsotopicH, 
                 nNumRemovedProtons_D.nNumRemovedIsotopicH,
                 sizeof( nNumRemovedProtons_D.nNumRemovedIsotopicH ) ) ) {
        CompareInchiFlags[TAUT_YES] |= INCHIDIFF_REM_ISO_H;
    }

exit_function:

    if ( pINChISort1 ) inchi_free( pINChISort1 );
    if ( pINChISort2 ) inchi_free( pINChISort2 );

    return ret;
}
/******************************************************************************************************/
int CompareTwoPairsOfInChI( INChI *pInChI1[TAUT_NUM], INChI *pInChI2[TAUT_NUM],
                            int bMobileH, INCHI_MODE CompareInchiFlags[] )
{
    int iMobileH, err=0;
    INCHI_MODE cmp;
    for ( iMobileH = 0; iMobileH < TAUT_NUM; iMobileH ++ ) {
        if ( !pInChI1[iMobileH] != !pInChI2[iMobileH] ) {
            if ( iMobileH == TAUT_NON &&
                 pInChI1[TAUT_YES] && pInChI1[TAUT_YES] ) {
                CompareInchiFlags[iMobileH] |= INCHIDIFF_COMP_HLAYER;
            } else {
                CompareInchiFlags[iMobileH] |= INCHIDIFF_COMP_NUMBER;
            }
            continue;
        }
        if ( pInChI1[iMobileH] && pInChI2[iMobileH] ) {
            cmp = CompareReversedINChI3( pInChI1[iMobileH], pInChI2[iMobileH], NULL, NULL, &err );
            if ( cmp ) {
                CompareInchiFlags[iMobileH] |= cmp;
            }
        }
    }
    return err;
}
/******************************************************************************************************/
int CompareOneOrigInchiToRevInChI(StrFromINChI *pStruct, INChI *pInChI[TAUT_NUM], int bMobileH, int iComponent,
                                  long num_inp, char *szCurHdr,
                                  COMPONENT_REM_PROTONS *nCurRemovedProtons, INCHI_MODE CompareInchiFlags[])
{
    int ret = pStruct->RevInChI.nRetVal, err=0;
    INCHI_MODE cmp;
    if ( ret == _IS_OKAY || ret == _IS_WARNING ) {
        /* ignore bMobileH for now */
        int i, i0, b /* created type */, b0 /* requested type*/, j, k;
        /* pINChI[iINCHI][iComponent][bTaut] */
        /* i0 = requested Rec/Disconnected: 1/0 */
        /* i  = what InChI creaded out of the restored structure */
        /* b0 = requested Mobile/Fixed-H: 1/0 */
        /* b  = what InChI creaded out of the restored structure */
        i = i0 = pStruct->iINCHI;
        b = b0 = pStruct->iMobileH;
        if ( i == INCHI_REC && !pStruct->RevInChI.num_components[i] ) {
            i = INCHI_BAS;
        }
        if ( b == TAUT_NON && (!pStruct->RevInChI.pINChI[i] ||
                               !pStruct->RevInChI.pINChI[i][0][b] ||
                               !pStruct->RevInChI.pINChI[i][0][b]->nNumberOfAtoms ) ) {
            b = TAUT_YES;
        }
        if ( pStruct->bDeleted && (!pInChI[0] || pInChI[0]->bDeleted ) ) {
            return 0;
        }

        if ( pStruct->RevInChI.num_components[i] > 1 && 
             !pStruct->RevInChI.pINChI[i][1][b]->bDeleted ||
             pStruct->RevInChI.num_components[i] < 1 ) {
            CompareInchiFlags[bMobileH] |= INCHIDIFF_COMP_NUMBER;
        }
        if ( b != b0 || b != bMobileH || b0 != bMobileH || i > i0 ) {
            /* do not print messages about TAUT_YES instead of TAUT_NON */
            CompareInchiFlags[bMobileH] |= INCHIDIFF_COMP_HLAYER;
        }

        if ( pStruct->RevInChI.num_components[i] ) {
            /* compare InChI from restored structure; '0' in [i][0][b] is the first component */
            if ( b == TAUT_YES && pStruct->RevInChI.pINChI[i][0][b]->bDeleted && (!pInChI[0] || pInChI[0]->bDeleted ) ) {
                /* the 1st component is made out of proton(s) and the input component is missing or also a proton */
                cmp = 0;
            } else {
                cmp = CompareReversedINChI3( pStruct->RevInChI.pINChI[i][0][b], pInChI[0], NULL, NULL, &err );
                if ( cmp ) {
                    CompareInchiFlags[bMobileH] |= cmp;
                }
            }
            if ( b == b0 && b == TAUT_NON ) {
                if ( pStruct->RevInChI.pINChI[i][0][TAUT_YES] &&
                     !pStruct->RevInChI.pINChI[i][0][TAUT_YES]->bDeleted ||
                     pInChI[1] && !pInChI[1]->bDeleted ) {

                    /* in addition to fixed-H also compare mobile-H InChI */
                    cmp = CompareReversedINChI3( pStruct->RevInChI.pINChI[i][0][TAUT_YES], pInChI[1], NULL, NULL, &err );
                    if ( cmp ) {
                        CompareInchiFlags[TAUT_YES] |= cmp;
                    }
                }
                /* compare removed H */
                if ( pStruct->nNumRemovedProtonsMobHInChI != pStruct->RevInChI.pINChI_Aux[i][0][TAUT_YES]->nNumRemovedProtons ) {
                    CompareInchiFlags[TAUT_YES] |= INCHIDIFF_MOBH_PROTONS;
                }
            }
            memset( nCurRemovedProtons, 0, sizeof(*nCurRemovedProtons) );
            for ( k = 0; k < pStruct->RevInChI.num_components[i]; k ++ ) {
                if ( !k || pStruct->RevInChI.pINChI[i][k][TAUT_YES]->bDeleted ) {
                    /* get removed protons from the 1st component; add othere only if they are deleted protons */
                    nCurRemovedProtons->nNumRemovedProtons += pStruct->RevInChI.pINChI_Aux[i][k][TAUT_YES]->nNumRemovedProtons;
                    for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
                        nCurRemovedProtons->nNumRemovedIsotopicH[j] += pStruct->RevInChI.pINChI_Aux[i][k][TAUT_YES]->nNumRemovedIsotopicH[j];
                    }
                }
            }
        }
    } else {
        CompareInchiFlags[bMobileH] |= INCHIDIFF_STR2INCHI_ERR;
    }
    return err;
}
/*************************************************************************************/
INCHI_MODE CompareReversedStereoINChI3( INChI_Stereo *s1/* InChI from reversed struct */, INChI_Stereo *s2 /* input InChI */, ICR *picr)
{
    int ret = 0;
    int j1, j2, num_eq, num_dif, num_extra_undf, num_miss_undf, num_in1_only, num_in2_only;
    int bAddSb = !(picr->num_sb_undef_in1_only + picr->num_sb_in1_only + picr->num_sb_in2_only);
    int bAddSc = !(picr->num_sc_undef_in1_only + picr->num_sc_in1_only + picr->num_sc_in2_only);
    
    int nNumSc1 = s1? s1->nNumberOfStereoCenters : 0;
    int nNumSc2 = s2? s2->nNumberOfStereoCenters : 0;
    int nNumSb1 = s1? s1->nNumberOfStereoBonds   : 0;
    int nNumSb2 = s2? s2->nNumberOfStereoBonds   : 0;
    
    if ( (nNumSc1 || nNumSc1) &&
         ( nNumSc1 != nNumSc2 ||
           memcmp( s1->nNumber,  s2->nNumber,  nNumSc1*sizeof(s1->nNumber[0] ) ) ||
           memcmp( s1->t_parity, s2->t_parity, nNumSc1*sizeof(s1->t_parity[0]) ) ) ) {

        num_eq = num_dif = num_extra_undf = num_miss_undf = num_in1_only = num_in2_only = 0;
        for ( j1 = j2 = 0; j1 < nNumSc1 && j2 < nNumSc2; ) {
            if ( s1->nNumber[j1] ==  s2->nNumber[j2] ) {
                if ( s1->t_parity[j1] == s2->t_parity[j2] ) {
                    num_eq ++;
                } else {
                    num_dif ++;
                }
                j1 ++;
                j2 ++;
            } else
            if ( s1->nNumber[j1] < s2->nNumber[j2] ) {
                num_in1_only ++;
                if ( s1->t_parity[j1] == AB_PARITY_UNDF ) {
                    num_extra_undf ++;
                }
                if ( bAddSc ) {
                    if ( picr->num_sc_in1_only < ICR_MAX_SC_IN1_ONLY )
                        picr->sc_in1_only[picr->num_sc_in1_only ++] = j1;
                    if ( s1->t_parity[j1] == AB_PARITY_UNDF ) {
                        if ( picr->num_sc_undef_in1_only < ICR_MAX_SC_UNDF )
                            picr->sc_undef_in1_only[picr->num_sc_undef_in1_only ++] = j1;
                    }
                }
                j1 ++;
            } else {
                num_in2_only ++;
                if ( s2->t_parity[j2] == AB_PARITY_UNDF ) {
                    num_miss_undf ++;
                }
                if ( bAddSc ) {
                    if ( picr->num_sc_in2_only < ICR_MAX_SC_IN2_ONLY )
                        picr->sc_in2_only[picr->num_sc_in2_only ++] = j2;
                    if ( s2->t_parity[j2] == AB_PARITY_UNDF ) {
                        if ( picr->num_sc_undef_in2_only < ICR_MAX_SC_UNDF )
                            picr->sc_undef_in2_only[picr->num_sc_undef_in2_only ++] = j1;
                    }
                }
                j2 ++;
            }
        }
        while ( j1 < nNumSc1 ) {
            if ( s1->t_parity[j1] == AB_PARITY_UNDF ) {
                num_extra_undf ++;
            }
            num_in1_only ++;
            if ( bAddSc ) {
                if ( picr->num_sc_in1_only < ICR_MAX_SC_IN1_ONLY )
                    picr->sc_in1_only[picr->num_sc_in1_only ++] = j1;
                if ( s1->t_parity[j1] == AB_PARITY_UNDF ) {
                    if ( picr->num_sc_undef_in1_only < ICR_MAX_SC_UNDF )
                        picr->sc_undef_in1_only[picr->num_sc_undef_in1_only ++] = j1;
                }
            }
            j1 ++;
        }
        while ( j2 < nNumSc2 ) {
            if ( s2->t_parity[j2] == AB_PARITY_UNDF ) {
                num_miss_undf ++;
            }
            num_in2_only ++;
            if ( bAddSc ) {
                if ( picr->num_sc_in2_only < ICR_MAX_SC_IN2_ONLY )
                    picr->sc_in2_only[picr->num_sc_in2_only ++] = j2;
            }
            j2 ++;
        }
        if ( num_dif ) {
            ret |= INCHIDIFF_SC_PARITY; 
        }
        if ( num_in1_only ) {
            if ( num_extra_undf ) {
                ret |= INCHIDIFF_SC_EXTRA_UNDF;
            }
            if ( num_in1_only != num_extra_undf ) {
                ret |= INCHIDIFF_SC_EXTRA;
            }
        }
        if ( num_in2_only ) {
            if ( num_miss_undf ) {
                ret |= INCHIDIFF_SC_MISS_UNDF;
            }
            if ( num_in2_only != num_miss_undf ) {
                ret |= INCHIDIFF_SC_MISS;
            }
        }
    }
    if ( s1 && s2 && (s2->nCompInv2Abs != 2) && s1->nCompInv2Abs != s2->nCompInv2Abs && s1->nCompInv2Abs && s2->nCompInv2Abs ) {
        ret |= INCHIDIFF_SC_INV; /* 2007-07-13 DT: added (s2->nCompInv2Abs != 2) to fix bug reoprted by Yerin on 2007/02/28 */
                                 /* Bug description: falsely reported "Stereo centers/allenes: Falsely inverted" for /S2 or /S3 */
    }

    if ( (nNumSb1 || nNumSb2 ) &&
         (nNumSb1 != nNumSb2 ||
          memcmp( s1->nBondAtom1, s2->nBondAtom1, nNumSb1*sizeof(s1->nBondAtom1[0]) ) ||
          memcmp( s1->nBondAtom2, s2->nBondAtom2, nNumSb1*sizeof(s1->nBondAtom2[0]) ) ||
          memcmp( s1->b_parity,   s2->b_parity,   nNumSb1*sizeof(s1->b_parity[0]) ) ) ) {

        num_eq = num_dif = num_extra_undf = num_miss_undf = num_in1_only = num_in2_only = 0;
        for ( j1 = j2 = 0; j1 < nNumSb1 && j2 < nNumSb2; ) {
            if ( s1->nBondAtom1[j1] ==  s2->nBondAtom1[j2] &&
                 s1->nBondAtom2[j1] ==  s2->nBondAtom2[j2] ) {
                if ( s1->b_parity[j1] == s2->b_parity[j2] ) {
                    num_eq ++;
                } else {
                    num_dif ++;
                }
                j1 ++;
                j2 ++;
            } else
            if ( s1->nBondAtom1[j1] <  s2->nBondAtom1[j2] ||
                 s1->nBondAtom1[j1] == s2->nBondAtom1[j2] && s1->nBondAtom2[j1] <  s2->nBondAtom2[j2]) {
                num_in1_only ++;
                if ( s1->b_parity[j1] == AB_PARITY_UNDF ) {
                    num_extra_undf ++;
                }
                if ( bAddSb ) {
                    if ( picr->num_sb_in1_only < ICR_MAX_SB_IN1_ONLY )
                        picr->sb_in1_only[picr->num_sb_in1_only ++] = j1;
                    if ( s1->b_parity[j1] == AB_PARITY_UNDF ) {
                        if ( picr->num_sb_undef_in1_only < ICR_MAX_SB_UNDF )
                            picr->sb_undef_in1_only[picr->num_sb_undef_in1_only ++] = j1;
                    }
                }
                j1 ++;
            } else {
                num_in2_only ++;
                if ( s2->b_parity[j2] == AB_PARITY_UNDF ) {
                    num_miss_undf ++;
                }
                if ( bAddSb ) {
                    if ( picr->num_sb_in2_only < ICR_MAX_SB_IN2_ONLY )
                        picr->sb_in2_only[picr->num_sb_in2_only ++] = j2;
                    if ( s2->b_parity[j2] == AB_PARITY_UNDF ) {
                        if ( picr->num_sb_undef_in2_only < ICR_MAX_SB_UNDF )
                            picr->sb_undef_in2_only[picr->num_sb_undef_in2_only ++] = j1;
                    }
                }
                j2 ++;
            }
        }
        while ( j1 < nNumSb1 ) {
            num_in1_only ++;
            if ( s1->b_parity[j1] == AB_PARITY_UNDF ) {
                num_extra_undf ++;
            }
            if ( bAddSb ) {
                if ( picr->num_sb_in1_only < ICR_MAX_SB_IN1_ONLY )
                    picr->sb_in1_only[picr->num_sb_in1_only ++] = j1;
                if ( s1->b_parity[j1] == AB_PARITY_UNDF ) {
                    if ( picr->num_sb_undef_in1_only < ICR_MAX_SB_UNDF )
                        picr->sb_undef_in1_only[picr->num_sb_undef_in1_only ++] = j1;
                }
            }
            j1 ++;
        }
        while ( j2 < nNumSb2 ) {
            num_in2_only ++;
            if ( s2->b_parity[j2] == AB_PARITY_UNDF ) {
                num_miss_undf ++;
            }
            if ( bAddSb ) {
                if ( picr->num_sb_in2_only < ICR_MAX_SB_IN2_ONLY )
                    picr->sb_in2_only[picr->num_sb_in2_only ++] = j2;
                if ( s2->b_parity[j2] == AB_PARITY_UNDF ) {
                    if ( picr->num_sb_undef_in2_only < ICR_MAX_SB_UNDF )
                        picr->sb_undef_in2_only[picr->num_sb_undef_in2_only ++] = j1;
                }
            }
            j2 ++;
        }
        if ( num_dif ) {
            ret |= INCHIDIFF_SB_PARITY; 
        }
        if ( num_in1_only ) {
            if ( num_extra_undf ) {
                ret |= INCHIDIFF_SB_EXTRA_UNDF;
            }
            if ( num_in1_only != num_extra_undf ) {
                ret |= INCHIDIFF_SB_EXTRA;
            }
        }
        if ( num_in2_only ) {
            if ( num_miss_undf ) {
                ret |= INCHIDIFF_SB_MISS_UNDF;
            }
            if ( num_in2_only != num_miss_undf ) {
                ret |= INCHIDIFF_SB_MISS;
            }
        }
    }

    return ret;
}
/*********************************************************************************************************/
INCHI_MODE CompareReversedINChI3( INChI *i1 /* InChI from reversed struct */, INChI *i2 /* input InChI */,
                                  INChI_Aux *a1, INChI_Aux *a2, int *err )
{
    INCHI_MODE ret = 0;
    INChI_Stereo *Stereo1=NULL, *Stereo2=NULL;
    int  n1, n2, m, j, j1, j2, ret2, num_H1, num_H2;
    ICR icr;
    ICR *picr = &icr;
    
    *err = 0;

    memset( picr, 0, sizeof(*picr) );

    if ( i1 == NULL && i2 == NULL )
        return 0;
    if ( (i1 == NULL) ^ (i2 == NULL) ) {
        ret |= INCHIDIFF_PROBLEM; /* one InChI exists while another doesn't */
        goto exit_function;
    }
    
    if ( i1->nErrorCode == i2->nErrorCode ) {
        if ( i1->nErrorCode ) {
            ret |= INCHIDIFF_PROBLEM; /* both InChI have same error codes */
            goto exit_function;
        }
    } else {
        ret |= INCHIDIFF_PROBLEM; /* at least one InChI has an error code */
        goto exit_function;
    }
    
    if ( i1->nNumberOfAtoms != i2->nNumberOfAtoms ) {
        ret |= INCHIDIFF_NUM_AT;
        goto exit_function;
    }
    if ( i1->nNumberOfAtoms > 0 ) {
        if ( memcmp( i1->nAtom, i2->nAtom, i1->nNumberOfAtoms*sizeof(i1->nAtom[0]) ) ) {
            ret |= INCHIDIFF_ATOMS;
            goto exit_function;
        }
        /* INCHIDIFF_NON_TAUT_H,  INCHIDIFF_MORE_FH, INCHIDIFF_LESS_FH */
        if ( memcmp( i1->nNum_H, i2->nNum_H, i1->nNumberOfAtoms*sizeof(i1->nNum_H[0]) ) ) {
            ret |= INCHIDIFF_POSITION_H;
            for ( j1 = 0; j1 < i1->nNumberOfAtoms; j1 ++ ) {
                if ( i1->nNum_H[j1] != i2->nNum_H[j1] && picr->num_diff_pos_H < ICR_MAX_DIFF_FIXED_H ) {
                    picr->diff_pos_H_at[picr->num_diff_pos_H] = j1;
                    picr->diff_pos_H_nH[picr->num_diff_pos_H] = i1->nNum_H[j1] - i2->nNum_H[j1];
                    picr->num_diff_pos_H ++;
                }
            }
        }
        /* fixed H */
        if ( i1->nNum_H_fixed || i2->nNum_H_fixed ) {
            int bHasFixedH1 = 0, bHasFixedH2 = 0, i;
            if ( i1->nNum_H_fixed ) {
                for ( i = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( i1->nNum_H_fixed[i] ) {
                        bHasFixedH1 ++;
                    }
                }
            }
            if ( i2->nNum_H_fixed ) {
                for ( i = 0; i < i2->nNumberOfAtoms; i ++ ) {
                    if ( i2->nNum_H_fixed[i] ) {
                        bHasFixedH2 ++;
                    }
                }
            }
            if ( bHasFixedH1 && !bHasFixedH2 ) {
                for ( i = j = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( i1->nNum_H_fixed[i] ) {
                        if ( j < ICR_MAX_DIFF_FIXED_H ) {
                            picr->fixed_H_at1_more[j] = i;
                            picr->fixed_H_nH1_more[j] = i1->nNum_H_fixed[i];
                            j ++;
                        }
                    }
                }
                picr->num_fixed_H1_more = j;
                ret |= INCHIDIFF_MORE_FH; /* Extra Fixed-H */
            } else
            if ( !bHasFixedH1 && bHasFixedH2 ) {
                for ( i = j = 0; i < i2->nNumberOfAtoms; i ++ ) {
                    if ( i2->nNum_H_fixed[i] ) {
                        if ( j < ICR_MAX_DIFF_FIXED_H ) {
                            picr->fixed_H_at2_more[j] = i;
                            picr->fixed_H_nH2_more[j] = i2->nNum_H_fixed[i];
                            j ++;
                        }
                    }
                }
                picr->num_fixed_H2_more = j;
                ret |= INCHIDIFF_LESS_FH; /* Missed Fixed-H */
            } else
            if ( bHasFixedH1 && bHasFixedH2 &&
                 memcmp( i1->nNum_H_fixed, i2->nNum_H_fixed, i1->nNumberOfAtoms*sizeof(i1->nNum_H_fixed[0]) ) ) {
                for ( i = j1 = j2 = 0; i < i1->nNumberOfAtoms; i ++ ) {
                    if ( i1->nNum_H_fixed[i] > i2->nNum_H_fixed[i] ) {
                        if ( j1 < ICR_MAX_DIFF_FIXED_H ) {
                            picr->fixed_H_at1_more[j1] = i;
                            picr->fixed_H_nH1_more[j1] = i1->nNum_H_fixed[i] - i2->nNum_H_fixed[i];
                            j1 ++;
                        }
                    } else
                    if ( i1->nNum_H_fixed[i] < i2->nNum_H_fixed[i] ) {
                        if ( j2 < ICR_MAX_DIFF_FIXED_H ) {
                            picr->fixed_H_at2_more[j2] = i;
                            picr->fixed_H_nH2_more[j2] = i2->nNum_H_fixed[i] - i1->nNum_H_fixed[i];
                            j2 ++;
                        }
                    }
                }
                ret |= (j1? INCHIDIFF_MORE_FH:0) | (j2? INCHIDIFF_LESS_FH:0);
                picr->num_fixed_H1_more = j1;
                picr->num_fixed_H2_more = j2;
            }
        }
    }
    /* compare formulas and H */
    num_H1 = 0;
    num_H2 = 0;
    ret2 = CompareHillFormulasNoH( i1->szHillFormula, i2->szHillFormula, &num_H1, &num_H2 );
    picr->tot_num_H1 = num_H1;
    picr->tot_num_H2 = num_H2;
    if ( ret2 ) {
        ret |= INCHIDIFF_NUM_EL;
        goto exit_function;
    }
    if ( num_H1 > num_H2 ) {
        ret |= INCHIDIFF_MORE_H;
    }
    if ( num_H1 < num_H2 ) {
        ret |= INCHIDIFF_LESS_H;
    }

    if ( i1->lenConnTable != i2->lenConnTable ) {
        ret |= INCHIDIFF_CON_LEN;
        goto exit_function;
    } else
    if ( i1->lenConnTable > 0 && memcmp( i1->nConnTable, i2->nConnTable, i1->lenConnTable*sizeof(i1->nConnTable[0]) ) ) {
        ret |= INCHIDIFF_CON_TBL;
        goto exit_function;
    }
    /* output special cases: different number of t-groups, different sizes of t-groups, different endpoints */
    /* in isotopic or deprotonated cases i1->lenTautomer == 1 && i1->nTautomer[0] = 0 */
/*
    if ( i1->lenTautomer != i2->lenTautomer && (i1->lenTautomer > 1 || i2->lenTautomer > 1) ) {
        ret |=  INCHIDIFF_TAUT_LEN; 
    }
*/
    /* compare number of t-groups */
    n1 = i1->lenTautomer? i1->nTautomer[0] : 0;
    n2 = i2->lenTautomer? i2->nTautomer[0] : 0;
    if ( !n1 && n2 ) {
        ret |= INCHIDIFF_NO_TAUT;
    } else
    if ( n1 && !n2 ) {
        ret |= INCHIDIFF_WRONG_TAUT;
    } else
    if ( n1 == 1 && n2 > 1 ) {
        ret |= INCHIDIFF_SINGLE_TG;
    } else
    if ( n1 > 1 && n2 == 1 ) {
        ret |= INCHIDIFF_MULTIPLE_TG;
    } else
    if ( n1 != n2 ) {
        ret |= INCHIDIFF_NUM_TG;
    }
    if ( n1 || n2 ) {
        /* number of endpoints */
        int num1 = 0, num2 = 0, num_M1=0, num_M2=0;
        int len, num_eq, num_in1_only, num_in2_only;
        AT_NUMB *pe1 = (AT_NUMB *) inchi_malloc( (i1->lenTautomer+1) * sizeof(pe1[0]) );
        AT_NUMB *pe2 = (AT_NUMB *) inchi_malloc( (i2->lenTautomer+1) * sizeof(pe2[0]) );
        num_H1 = num_H2=0;
        /* collect endpoints, H, (-) */
        if ( !pe1 || !pe2 ) {
            if ( pe1 ) inchi_free( pe1 );
            if ( pe2 ) inchi_free( pe2 );
            *err = RI_ERR_ALLOC; /* allocation error */
            goto exit_function;
        }
        for ( m = 1; m < i1->lenTautomer; m += len ) {
            len = i1->nTautomer[m ++];
            num_H1 += i1->nTautomer[m];
            num_M1 += i1->nTautomer[m+1];
            for ( j = 2; j < len; j ++ ) {
                pe1[num1 ++] = i1->nTautomer[m + j];
            }
        }
        for ( m = 1; m < i2->lenTautomer; m += len ) {
            len = i2->nTautomer[m ++];
            num_H2 += i2->nTautomer[m];
            num_M2 += i2->nTautomer[m+1];
            for ( j = 2; j < len; j ++ ) {
                pe2[num2 ++] = i2->nTautomer[m + j];
            }
        }
        picr->num_taut_H1 = num_H1;
        picr->num_taut_H2 = num_H2;
        picr->num_taut_M1 = num_M1;
        picr->num_taut_M2 = num_M2;
        /* sort endpoints */
        insertions_sort_AT_NUMB( pe1, num1 );
        insertions_sort_AT_NUMB( pe2, num2 );
        /* compare */
        /*
        if ( num1 < num2 ) {
            ret |= INCHIDIFF_LESS_TG_ENDP;
        } else
        if ( num1 > num2 ) {
            ret |= INCHIDIFF_MORE_TG_ENDP;
        }
        */
        /* compare all */
        num_eq = num_in1_only = num_in2_only = 0;
        for ( j1 = j2 = 0; j1 < num1 && j2 < num2; ) {
            if( pe1[j1] == pe2[j2] ) {
                j1 ++;
                j2 ++;
                num_eq ++;
            } else
            if ( pe1[j1] < pe2[j1] ) {
                if ( picr->num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY ) {
                    picr->endp_in1_only[picr->num_endp_in1_only ++] = pe1[j1];
                }
                j1 ++;
                num_in1_only ++;
            } else {
                if ( picr->num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY ) {
                    picr->endp_in2_only[picr->num_endp_in2_only ++] = pe2[j2];
                }
                j2 ++;
                num_in2_only ++;
            }
        }
        while ( j1 < num1 ) {
            if ( picr->num_endp_in1_only < ICR_MAX_ENDP_IN1_ONLY ) {
                picr->endp_in1_only[picr->num_endp_in1_only ++] = pe1[j1];
            }
            j1 ++;
            num_in1_only ++;
        }
        while ( j2 < num2 ) {
            if ( picr->num_endp_in2_only < ICR_MAX_ENDP_IN2_ONLY ) {
                picr->endp_in2_only[picr->num_endp_in2_only ++] = pe2[j2];
            }
            j2 ++;
            num_in2_only ++;
        }
        if ( num_in1_only ) {
            ret |= INCHIDIFF_EXTRA_TG_ENDP;
        }
        if ( num_in2_only ) {
            ret |= INCHIDIFF_MISS_TG_ENDP;
        }
        if ( !num_in1_only && !num_in2_only && num_eq ) {
           ; /* same t-groups endpoints */
        } else {
           ret |= INCHIDIFF_DIFF_TG_ENDP;
        }
        inchi_free( pe1 );
        inchi_free( pe2 );

    }

    if ( (i1->lenTautomer > 1 && i2->lenTautomer > 1) &&
         ( i1->lenTautomer != i2->lenTautomer ||
         memcmp( i1->nTautomer, i2->nTautomer, i1->lenTautomer*sizeof(i1->nTautomer[0]) ) ) )
        ret |= INCHIDIFF_TG;

    if ( i1->nNumberOfIsotopicAtoms != i2->nNumberOfIsotopicAtoms ) {
        ret |= INCHIDIFF_NUM_ISO_AT;
    } else
    if ( i1->nNumberOfIsotopicAtoms > 0 && memcmp( i1->IsotopicAtom, i2->IsotopicAtom, i1->nNumberOfIsotopicAtoms*sizeof(i1->IsotopicAtom[0]) ) )
        ret |= INCHIDIFF_ISO_AT;
    if ( i1->nTotalCharge != i2->nTotalCharge )
        ret |= INCHIDIFF_CHARGE;
    if ( a1 && a1->nNumRemovedProtons && (!a2 || a2->nNumRemovedProtons != a1->nNumRemovedProtons) ) {
        ret |= INCHIDIFF_REM_PROT;
    }
    if ( a1 && (!a2 || 
         a2->nNumRemovedIsotopicH[0] != a1->nNumRemovedIsotopicH[0] ||
         a2->nNumRemovedIsotopicH[1] != a1->nNumRemovedIsotopicH[1] ||
         a2->nNumRemovedIsotopicH[2] != a1->nNumRemovedIsotopicH[2]) ) {
        ret |= INCHIDIFF_REM_ISO_H;
    }

/*
    if ( i1->nPossibleLocationsOfIsotopicH && i2->nPossibleLocationsOfIsotopicH ) {
        if ( i1->nPossibleLocationsOfIsotopicH[0] != i2->nPossibleLocationsOfIsotopicH[0] ||
             memcmp(i1->nPossibleLocationsOfIsotopicH, i2->nPossibleLocationsOfIsotopicH,
                    sizeof(i1->nPossibleLocationsOfIsotopicH[0])*i1->nPossibleLocationsOfIsotopicH[0]) )
            return 18;
    } else
    if ( !i1->nPossibleLocationsOfIsotopicH != !i2->nPossibleLocationsOfIsotopicH ) {
        return 19;
    }
*/
    if ( i1->StereoIsotopic &&
         i1->StereoIsotopic->nNumberOfStereoBonds + i1->StereoIsotopic->nNumberOfStereoCenters ) {
        Stereo1 = i1->StereoIsotopic;
    } else {
        Stereo1 = i1->Stereo;
    }
    if ( i2->StereoIsotopic &&
         i2->StereoIsotopic->nNumberOfStereoBonds + i2->StereoIsotopic->nNumberOfStereoCenters ) {
        Stereo2 = i2->StereoIsotopic;
    } else {
        Stereo2 = i2->Stereo;
    }
    ret |= CompareReversedStereoINChI3( Stereo1, Stereo2, picr );

exit_function:

    picr->flags = ret;

    return ret;
}
/* message group names */
CMP_INCHI_MSG_GROUP CompareInchiMsgsGroup[] = {
{IDGRP_ERR,     " Error:"},
{IDGRP_H,       " Hydrogens:"},
{IDGRP_MOB_GRP, " Mobile-H groups:"},
{IDGRP_ISO_AT,  " Isotopic:"},
{IDGRP_CHARGE,  " Charge(s):"},
{IDGRP_PROTONS, " Proton balance:"},
{IDGRP_ISO_H,   " Exchangeable isotopic H:"},
{IDGRP_SC,      " Stereo centers/allenes:"},
{IDGRP_SB,      " Stereobonds/cumulenes:"},
{IDGRP_HLAYER,  " Fixed-H layer:"},
{IDGRP_COMP,    " Number of components:"},
{IDGRP_CONV_ERR," Conversion encountered:"},
{IDGRP_ZERO,    ""}
};
/* messages */
CMP_INCHI_MSG  CompareInchiMsgs[] = {
{INCHIDIFF_PROBLEM      ,IDGRP_ERR,     " Wrong result"                   }, /*0x00000001,  severe: at least one InChI does not exist */
{INCHIDIFF_POSITION_H   ,IDGRP_H,       " Locations or number"            }, /*0x00000002,  difference in non-taut {Mobile-H} or all H {Fixed-H} location/number */
{INCHIDIFF_MORE_FH      ,IDGRP_H,       " Fixed-H"                        }, /*0x00000004,  extra fixed H */
{INCHIDIFF_LESS_FH      ,IDGRP_H,       " Fixed-H"                        }, /*0x00000004,  missing fixed H */
{INCHIDIFF_MORE_H       ,IDGRP_H,       " Number"                         }, /*0x00000008,  formulas differ in number of H */
{INCHIDIFF_LESS_H       ,IDGRP_H,       " Number"                         }, /*0x00000008,  formulas differ in number of H */
{INCHIDIFF_NO_TAUT      ,IDGRP_MOB_GRP, " Missing"                        }, /*0x00000010,  restored structure has no taut groups while the original InChI has some */
{INCHIDIFF_WRONG_TAUT   ,IDGRP_MOB_GRP, " Falsely present"                }, /*0x00000020,  restored has tautomerism while the original does not have it */
{INCHIDIFF_SINGLE_TG    ,IDGRP_MOB_GRP, " One instead of multiple"        }, /*0x00000040,  restored has 1 taut. group while the original InChI has multiple tg */
{INCHIDIFF_MULTIPLE_TG  ,IDGRP_MOB_GRP, " Multiple instead of one"        }, /*0x00000080,  restored has multiple tg while the original InChI has only one tg */
{INCHIDIFF_EXTRA_TG_ENDP,IDGRP_MOB_GRP, " Attachment points"              }, /*0x00000100,  extra tautomeric endpoint{s} in restored structure */
{INCHIDIFF_MISS_TG_ENDP ,IDGRP_MOB_GRP, " Attachment points"              }, /*0x00000100,  one or more tg endpoint is not in the restored structure */
{INCHIDIFF_DIFF_TG_ENDP ,IDGRP_MOB_GRP, " Attachment points"              }, /*0x00000100,  lists of tg endpoints are different */
{INCHIDIFF_NUM_TG       ,IDGRP_MOB_GRP, " Number"                         }, /*0x00000200,  different number of tautomeric groups */
{INCHIDIFF_TG           ,IDGRP_MOB_GRP, " Do not match"                   }, /*0x00000200,  different tautomeric groups */
{INCHIDIFF_NUM_ISO_AT   ,IDGRP_ISO_AT,  " Atoms do not match"             }, /*0x00000400,  ?severe: restored struct. has different number of isotopic atoms */
{INCHIDIFF_ISO_AT       ,IDGRP_ISO_AT,  " Atoms do not match"             }, /*0x00000400,  ?severe: restored struct. has different locations/isotopes of isotopic atoms */
{INCHIDIFF_REM_ISO_H    ,IDGRP_ISO_H,   " Does not match for a component" }, /*0x00000800,  isotopic H removed */
{INCHIDIFF_MOB_ISO_H    ,IDGRP_ISO_H,   " Do not match"                   }, /*0x00001000,  different number of mobile exchangeable isotopic H */
{INCHIDIFF_CHARGE       ,IDGRP_CHARGE,  " Do not match"                   }, /*0x00002000,  restored structure has different charge */
{INCHIDIFF_REM_PROT     ,IDGRP_PROTONS, " Does not match for a component" }, /*0x00004000,  proton{s} removed/added from the restored structure */
{INCHIDIFF_MOBH_PROTONS ,IDGRP_PROTONS, " Does not match"                 }, /*0x00008000,  different proton balance */
{INCHIDIFF_SC_INV       ,IDGRP_SC,      " Falsely inverted"               }, /*0x00010000,  restores structure has different inversion stereocenter mark */
{INCHIDIFF_SC_PARITY    ,IDGRP_SC,      " Wrong parity"                   }, /*0x00020000,  restored structure has stereoatoms or allenes with different parity */
{INCHIDIFF_SC_EXTRA_UNDF,IDGRP_SC,      " Extra undefined"                }, /*0x00040000,  restored structure has extra undefined stereocenter{s} */
{INCHIDIFF_SC_EXTRA     ,IDGRP_SC,      " Extra known"                    }, /*0x00080000,  restored structure has extra stereocenter{s} */
{INCHIDIFF_SC_MISS_UNDF ,IDGRP_SC,      " Missing undefined"              }, /*0x00100000,  restored structure has not some undefined stereocenter{s} */
{INCHIDIFF_SC_MISS      ,IDGRP_SC,      " Missing known"                  }, /*0x00200000,  restored structure has not some stereocenters that are not undefined */
{INCHIDIFF_SB_PARITY    ,IDGRP_SB,      " Wrong parity"                   }, /*0x00400000,  restored structure has stereobonds or cumulenes with different parity */
{INCHIDIFF_SB_EXTRA_UNDF,IDGRP_SB,      " Extra undefined"                }, /*0x00800000,  restored structure has extra undefined stereobond{s} */
{INCHIDIFF_SB_EXTRA     ,IDGRP_SB,      " Missing known"                  }, /*0x01000000,  restored structure has extra stereobond{s} */
{INCHIDIFF_SB_MISS_UNDF ,IDGRP_SB,      " Missing undefined"              }, /*0x02000000,  restored structure has not some undefined stereocenters */
{INCHIDIFF_SB_MISS      ,IDGRP_SB,      " Missing known"                  }, /*0x04000000,  restored structure has not some stereobonds that are not undefined */
{INCHIDIFF_COMP_HLAYER  ,IDGRP_HLAYER,  " Missing or extra"               }, /*0x08000000,  Restored component has Mobile-H layer instead of both Mobile-H & Fixed-H or both instead of one */
{INCHIDIFF_COMP_NUMBER  ,IDGRP_COMP,    " Does not match"                 }, /*0x10000000,  wrong number of components */
{INCHIDIFF_STR2INCHI_ERR,IDGRP_CONV_ERR," Error"                          },  /*0x20000000   Restored structure to InChI conversion error */
{INCHIDIFF_ZERO         ,IDGRP_ZERO,    ""                                }
};

/*************************************************************************/
int AddOneMsg( char *szMsg, int used_len, int tot_len, const char *szAddMsg, const char *szDelim )
{
    const char ellip[] = "...";
    int len = strlen( szAddMsg );
    int len_delim = (used_len && szDelim)? strlen(szDelim) : 0;
    int len_to_copy;
    if ( len + len_delim + used_len < tot_len ) {
        if ( len_delim ) {
            strcpy( szMsg+used_len, szDelim );
            used_len += len_delim;
        }
        strcpy( szMsg+used_len, szAddMsg );
        used_len += len;
    } else
    if ( (len_to_copy = (tot_len - used_len - len_delim - (int)sizeof(ellip))) > 10 ) {
        if ( len_delim ) {
            strcpy( szMsg+used_len, szDelim );
            used_len += len_delim;
        }
        strncpy( szMsg+used_len, szAddMsg, len_to_copy );
        used_len += len_to_copy;
        strcpy( szMsg+used_len, ellip );
        used_len += sizeof( ellip ) - 1;
    }
    return used_len;
}
/*************************************************************************/
int FillOutCompareMessage( char *szMsg, int nLenMsg, INCHI_MODE bits[] )
{
    int bMobileH, k, n, len = strlen( szMsg );
    int iPrevGrpIdx, iCurGrpIdx, bFound;
    INCHI_MODE bit;
    static const char *hdr = " Problems/mismatches:";
    char szOneMsg[256];
    if ( bits[TAUT_YES] || bits[TAUT_NON] ) {
        if ( !strstr( szMsg, hdr ) ) {
            len = AddOneMsg( szMsg, len, nLenMsg, hdr, NULL );
        }
        for ( bMobileH = TAUT_YES; 0 <= bMobileH; bMobileH -- ) {
            if ( bits[bMobileH] ) {
                strcpy( szOneMsg, bMobileH==TAUT_YES? " Mobile-H(" : " Fixed-H(" );
                len = AddOneMsg( szMsg, len, nLenMsg, szOneMsg, NULL );
            }
            bit = 1;
            iPrevGrpIdx = -1;
            do {
                if ( bit & bits[bMobileH] ) {
                    /* search for the message */
                    bFound = 0;
                    for ( k = 0; CompareInchiMsgs[k].nBit != INCHIDIFF_ZERO && !bFound; k ++ ) {
                        if ( bit & (INCHI_MODE)CompareInchiMsgs[k].nBit ) {
                            /* message found */
                            for ( n = 0; CompareInchiMsgsGroup[n].nGroupID != IDGRP_ZERO; n ++ ) {
                                if ( CompareInchiMsgsGroup[n].nGroupID == CompareInchiMsgs[k].nGroupID ) {
                                    iCurGrpIdx = n;
                                    if ( iCurGrpIdx != iPrevGrpIdx ) {
                                        if ( iPrevGrpIdx >= 0 ) {
                                            len = AddOneMsg( szMsg, len, nLenMsg, ";", NULL );
                                        }
                                        len = AddOneMsg( szMsg, len, nLenMsg, CompareInchiMsgsGroup[iCurGrpIdx].szGroupName, NULL );
                                    }
                                    len = AddOneMsg( szMsg, len, nLenMsg, CompareInchiMsgs[k].szMsg, iCurGrpIdx == iPrevGrpIdx? ",":NULL );
                                    iPrevGrpIdx = iCurGrpIdx;
                                    bFound = 1;
                                    break;
                                }
                            }
                        }
                    }
                }
                bit <<= 1;
            } while ( bit );
            if ( bits[bMobileH] ) {
                len = AddOneMsg( szMsg, len, nLenMsg, ")", NULL );
            }
        }
    }
    return len;
}

#endif
