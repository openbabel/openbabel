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
#include <stdarg.h>
/* #include <varargs.h> */
#include <errno.h>
#include <limits.h>



#include "mode.h"       /* moved from below, suggestion by David Mosenkis */

#include "ichitime.h"

#ifndef COMPILE_ANSI_ONLY
#include <conio.h>
#endif

#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichidrp.h"
#include "ichierr.h"
#include "ichimain.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichi_io.h"

#ifdef TARGET_LIB_FOR_WINCHI
#include "ichi_lib.h"
#endif
#include "inchi_api.h"

#include "ichicomp.h"

#if ( ADD_CMLPP == 1 )
#include "readcml.hpp"
#include "debug.h"
#endif


/* for DisplayTheWholeStructure() */

#define COMP_ORIG_0_MAIN  0x0001
#define COMP_ORIG_0_RECN  0x0002
#define COMP_PREP_0_MAIN  0x0004
#define COMP_PREP_0_RECN  0x0008
#define COMP_ORIG_1_MAIN  0x0010
#define COMP_ORIG_1_RECN  0x0020


/* local prototypes */
int GetProcessingWarningsOneINChI(INChI *pINChI, INP_ATOM_DATA *inp_norm_data, char *pStrErrStruct);
int  GetProcessingWarnings(INChI *cur_INChI[], INP_ATOM_DATA **inp_norm_data, STRUCT_DATA *sd);
int DisplayTheWholeStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, 
                              INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file,
                              ORIG_ATOM_DATA *orig_inp_data, long num_inp, int iINChI, int bShowStruct, int bINCHI_LIB_Flag );
int DuplicateOrigAtom( ORIG_ATOM_DATA *new_orig_atom, ORIG_ATOM_DATA *orig_atom );
int bCheckUnusualValences( ORIG_ATOM_DATA *orig_at_data, int bAddIsoH,  char *pStrErrStruct );
int CreateCompositeNormAtom( COMP_ATOM_DATA *composite_norm_data, INP_ATOM_DATA2 *all_inp_norm_data,
                             PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, int num_components, INCHI_MODE nMode );
int DetectInputINChIFileType( FILE **inp_file, INPUT_PARMS *ip, const char *fmode );


/*  callback */
int (*ConsoleQuit)(void) = NULL; /*  Console user issued CTRL+C etc. */
int (*UserAction)(void)  = NULL; /*  callback */

#ifdef TARGET_LIB_FOR_WINCHI
void (*FWPRINT) (const char * format, va_list argptr )=NULL;
void (*DRAWDATA) ( struct DrawData * pDrawData) = NULL;
int (*DRAWDATA_EXISTS) ( int nComponent, int nType, int bReconnected ) = NULL;
struct DrawData * (*GET_DRAWDATA) ( int nComponent, int nType, int bReconnected ) = NULL;
#endif

#if ( TEST_RENUMB_ATOMS == 1 ) /* { */
/************************************************/
/*    atoms renumbering -- for testing only     */
/************************************************/
typedef struct tagRenumbData {
    PINChI2         ren_INChI2[1];
    PINChI_Aux2     ren_INChI_Aux[1];
    INP_ATOM_DATA  orig_inp_cur_data;
    INP_ATOM_DATA  saved_inp_cur_data;
#if ( TEST_RENUMB_ATOMS_SAVE_LONGEST == 1 || TEST_RENUMB_SWITCH == 1 )
    INP_ATOM_DATA  longest_inp_cur_data;
#endif
    INP_ATOM_DATA  ren_inp_norm_data1, ren_inp_norm_data2;
    INP_ATOM_DATA  *ren_inp_norm_data[2];
    int            ren_counter;
    int            num_taut, num_non_taut, num_taut0, num_non_taut0;
    AT_RANK       *new_ord;
    int            nRet2, c1, c2, nComp, bRenumbErr;
    unsigned long  ulCurTimeNorm0, ulCurTimeCanon0, ulCurTimeNorm1, ulCurTimeCanon1;
    unsigned long  ulCurTimeNorm, ulCurTimeCanon, ulMaxTimeNorm, ulMaxTimeCanon;
    unsigned long  ulMaxTime, ulCurTime, ulCurTime0, ulCurTime1;
#if ( bRELEASE_VERSION == 0 )
    int bExtract;
#endif
} RENUMB_DATA;

int RenumberingTestInit( RENUMB_DATA *pRenumbData, INP_ATOM_DATA *inp_cur_data );
int RenumberingTestUninit( RENUMB_DATA *pRenumbData );
int RenumberingTest( PINChI2 *pICh, PINChI_Aux2 *pINChI_Aux, ORIG_ATOM_DATA *orig_inp_data, int iINChI,
                     RENUMB_DATA *pRenumbData, INP_ATOM_DATA *inp_cur_data, INP_ATOM_DATA **inp_norm_data,
                     STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *prb_file,
                     int i, long num_inp, NORM_CANON_FLAGS *pncFlags);
/*
int RenumberingTest( INChI *pINChI[][TAUT_NUM], INChI_Aux *pINChI_Aux[][TAUT_NUM], int iINChI,
                     RENUMB_DATA *pRenumbData, INP_ATOM_DATA *inp_cur_data, INP_ATOM_DATA **inp_norm_data,
                     STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, INCHI_IOSTREAM *log_file, int i, long num_inp);
*/
#endif /* } TEST_RENUMB_ATOMS */



#ifndef COMPILE_ANSI_ONLY
/********************************************************************/
void FillTableParms( SET_DRAW_PARMS *sdp, INChI **cur_INChI, INChI_Aux **cur_INChI_Aux,
                     INCHI_MODE nMode, int bShowIsotopic, int indx )
{
    TBL_DRAW_PARMS *tdp = sdp->tdp;
    char   (*ReqShownFound)[TDP_NUM_PAR] = tdp->ReqShownFound;
    int  i, j;
    INChI_Stereo *Stereo;
    int          bShowTaut = (cur_INChI && cur_INChI[indx]->lenTautomer > 0)? 1 : 0;
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    int bRelRac = 0 != (nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO ));
#endif
    if ( !cur_INChI || !cur_INChI_Aux ) {
        sdp->tdp->bDrawTbl = 0;
        sdp->bOrigAtom     = 1;
        return;
    }

    /*  Displayed */
    ReqShownFound[ilSHOWN][itBASIC]    =  bShowTaut?     'T':'\0';
    ReqShownFound[ilSHOWN][itISOTOPIC] =  bShowIsotopic? 'I':'\0';
/*
    ReqShownFound[ilSHOWN][itBASIC]    =  bShowTaut?     'T':'B';
    ReqShownFound[ilSHOWN][itISOTOPIC] =  bShowIsotopic? 'I':'N';
 */
    i = indx;
    if ( cur_INChI[i] ) {
        Stereo = bShowIsotopic? cur_INChI[i]->StereoIsotopic : cur_INChI[i]->Stereo;
    } else {
        Stereo = NULL;
    }
#if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    if ( Stereo && ( 0 < Stereo->nNumberOfStereoBonds ||
                     0 < Stereo->nNumberOfStereoCenters-bRelRac ) ) {
        ReqShownFound[ilSHOWN][itSTEREO] = 'S';
        if ( Stereo->nNumberOfStereoCenters && Stereo->nCompInv2Abs == -1 &&
             ( nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO ) ) ) {
            if ( Stereo->nNumberOfStereoCenters < 2 && !Stereo->nNumberOfStereoBonds ) {
                ReqShownFound[ilSHOWN][itSTEREO] = '\0';
            } else
            if ( Stereo->nNumberOfStereoCenters >= 2 ) {
                ReqShownFound[ilSHOWN][itSTEREO] =  's'; /* shown Inverted stereo */
            }
        }
#else  /* REL_RAC_STEREO_IGN_1_SC == 0 */
    if ( Stereo && ( Stereo->nNumberOfStereoBonds || Stereo->nNumberOfStereoCenters ) ) {
        ReqShownFound[ilSHOWN][itSTEREO] = 'S';
        if ( Stereo->nNumberOfStereoCenters && Stereo->nCompInv2Abs == -1 &&
             ( nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO ) ) ) {
            /*
            if ( Stereo->nNumberOfStereoCenters < 2 && !Stereo->nNumberOfStereoBonds ) {
                ReqShownFound[ilSHOWN][itSTEREO] = '\0';
            } else
            if ( Stereo->nNumberOfStereoCenters >= 2 ) {
            */
                ReqShownFound[ilSHOWN][itSTEREO] =  's'; /* shown Inverted stereo */
            /*
            }
            */
        }
#endif /* REL_RAC_STEREO_IGN_1_SC */
    } else {
        ReqShownFound[ilSHOWN][itSTEREO] = '\0';
    }
    /*
    ReqShownFound[ilSHOWN][itSTEREO]   =
        (bShowIsotopic? (cur_INChI[i] && cur_INChI[i]->StereoIsotopic &&
                         (cur_INChI[i]->StereoIsotopic->nNumberOfStereoBonds ||
                          cur_INChI[i]->StereoIsotopic->nNumberOfStereoCenters) )
                        :
                        (cur_INChI[i] && cur_INChI[i]->Stereo &&
                         (cur_INChI[i]->Stereo->nNumberOfStereoBonds ||
                          cur_INChI[i]->Stereo->nNumberOfStereoCenters) )
        ) ? 'S':'\0';
    */

    /* remove zeroes between chars */
    for ( i = j = 0; i < TDP_NUM_PAR; i ++ ) {
        if ( ReqShownFound[ilSHOWN][i] >= ' ' ) {
            ReqShownFound[ilSHOWN][j++] = ReqShownFound[ilSHOWN][i];
        }
    }
    i = j;
    for ( ; i < TDP_NUM_PAR; i ++ ) {
        ReqShownFound[ilSHOWN][i] = '\0';
    }

    sdp->tdp->bDrawTbl = j? 1 : 0;
    sdp->bOrigAtom     = 0;
}
/********************************************************************/
void FillCompositeTableParms( SET_DRAW_PARMS *sdp, AT_NUMB StereoFlags,
                     INCHI_MODE nMode, int bShowIsotopic, int bShowTaut )
{
    TBL_DRAW_PARMS *tdp = sdp->tdp;
    char    (*ReqShownFound)[TDP_NUM_PAR] = tdp->ReqShownFound;
    int  i, j;

    /*  Displayed */
    ReqShownFound[ilSHOWN][itBASIC]    =  bShowTaut?     'T':'\0';
    ReqShownFound[ilSHOWN][itISOTOPIC] =  bShowIsotopic? 'I':'\0';
/*
    ReqShownFound[ilSHOWN][itBASIC]    =  bShowTaut?     'T':'B';
    ReqShownFound[ilSHOWN][itISOTOPIC] =  bShowIsotopic? 'I':'N';
 */
    if ( StereoFlags & INF_STEREO ) {
        ReqShownFound[ilSHOWN][itSTEREO] = 'S';
        if ( (StereoFlags & INF_STEREO_INV) &&
            ( nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO ) ) ) {
            if (StereoFlags & (INF_STEREO_REL | INF_STEREO_RAC) ) {
                ReqShownFound[ilSHOWN][itSTEREO] = 's';
            } else {
                ReqShownFound[ilSHOWN][itSTEREO] =  '\0'; /* shown Inverted stereo */
            }
        }
    } else {
        ReqShownFound[ilSHOWN][itSTEREO] = '\0';
    }
    /*
    ReqShownFound[ilSHOWN][itSTEREO]   =
        (bShowIsotopic? (cur_INChI[i] && cur_INChI[i]->StereoIsotopic &&
                         (cur_INChI[i]->StereoIsotopic->nNumberOfStereoBonds ||
                          cur_INChI[i]->StereoIsotopic->nNumberOfStereoCenters) )
                        :
                        (cur_INChI[i] && cur_INChI[i]->Stereo &&
                         (cur_INChI[i]->Stereo->nNumberOfStereoBonds ||
                          cur_INChI[i]->Stereo->nNumberOfStereoCenters) )
        ) ? 'S':'\0';
    */

    /* remove zeroes between chars */
    for ( i = j = 0; i < TDP_NUM_PAR; i ++ ) {
        if ( ReqShownFound[ilSHOWN][i] >= ' ' ) {
            ReqShownFound[ilSHOWN][j++] = ReqShownFound[ilSHOWN][i];
        }
    }
    i = j;
    for ( ; i < TDP_NUM_PAR; i ++ ) {
        ReqShownFound[ilSHOWN][i] = '\0';
    }

    sdp->tdp->bDrawTbl = j? 1 : 0;
    sdp->bOrigAtom     = 0;
}
#endif
/* IchiParm.c was here */
/*******************************************************************/
#ifndef COMPILE_ANSI_ONLY
#ifndef TARGET_LIB_FOR_WINCHI
/*******************************************************************/
int DisplayStructure( inp_ATOM *at, int num_at, int num_removed_H, int bAdd_DT_to_num_H,
                      int nNumRemovedProtons, NUM_H *nNumRemovedProtonsIsotopic,
                      int bIsotopic, int j /*bTautomeric*/,
                      INChI **cur_INChI, INChI_Aux **cur_INChI_Aux,
                      int bAbcNumbers, DRAW_PARMS *dp, INCHI_MODE nMode, char *szTitle )
{
    INF_ATOM_DATA inf_data = {NULL,};
    int err = -1;
    if ( CreateInfoAtomData( &inf_data, num_at, 1 ) ) {
        err = 0;
        FillOutInfAtom( at, &inf_data, num_at, num_removed_H, bAdd_DT_to_num_H,
                        nNumRemovedProtons, nNumRemovedProtonsIsotopic, bIsotopic,
                        cur_INChI?cur_INChI[j]:NULL,
                        cur_INChI_Aux?cur_INChI_Aux[j]:NULL, bAbcNumbers, nMode);
        FillTableParms( &dp->sdp, cur_INChI, cur_INChI_Aux, nMode, bIsotopic, j );
        err = DisplayInputStructure( szTitle, at, &inf_data, num_at, dp );
        FreeInfoAtomData( &inf_data );
    }
    return err;
}

/*******************************************************************/
int DisplayCompositeStructure( COMP_ATOM_DATA *composite_norm_data, int bIsotopic, int bTautomeric,
                      PINChI2 *pINChI2, PINChI_Aux2 *pINChI_Aux2,
                      int bAbcNumbers, DRAW_PARMS *dp, INCHI_MODE nMode, char *szTitle )
{
    INF_ATOM_DATA inf_data;
    int err = -1, ret;
    memset( &inf_data, 0, sizeof(inf_data) );
    if ( CreateInfoAtomData( &inf_data, (composite_norm_data+bTautomeric)->num_at,
                              (composite_norm_data+bTautomeric)->num_components ) ) {
        ret = FillOutCompositeCanonInfAtom(composite_norm_data, &inf_data,
                                 bIsotopic, bTautomeric,
                                 pINChI2, pINChI_Aux2, bAbcNumbers, nMode);
        if ( !ret ) {
            goto exit_function; /* error */
        }
        if ( bTautomeric == TAUT_INI ) {
            /*
            FillOutInfAtom( (composite_norm_data+bTautomeric)->at, &inf_data, (composite_norm_data+bTautomeric)->num_at,
                            (composite_norm_data+bTautomeric)->num_removed_H, bAdd_DT_to_num_H,
                            (composite_norm_data+bTautomeric)->nNumRemovedProtons,
                            (composite_norm_data+bTautomeric)->nNumRemovedProtonsIsotopic, bIsotopic,
                            NULL, NULL, bAbcNumbers, nMode);
            */
            ;
        } else {
            /* real check for tautomeric components 02-04-2005 */
            int m, nNumTautComponents = 0;
            if ( 1 == bTautomeric ) {
                for ( m = 0; m < composite_norm_data[TAUT_YES].num_components; m ++ ) {
                    if ( !pINChI2[m][TAUT_YES] )
                        continue;
                    if ( pINChI2[m][TAUT_YES]->bDeleted || pINChI2[m][TAUT_YES]->lenTautomer > 0 )
                        nNumTautComponents ++;
                }
            }
            FillCompositeTableParms( &dp->sdp, inf_data.StereoFlags, nMode, bIsotopic, nNumTautComponents );
        }
        err = DisplayInputStructure( szTitle, (composite_norm_data+bTautomeric)->at, &inf_data, (composite_norm_data+bTautomeric)->num_at, dp );
        FreeInfoAtomData( &inf_data );
    }
exit_function:
    return err;
}
#endif
#endif
/************************************************/
const char *ErrMsg( int nErrorCode )
{
    const char *p;
    static char szErrMsg[64];
    switch( nErrorCode ) {
        case 0:                      p = "";                      break;
        case CT_OVERFLOW:            p = "ARRAY OVERFLOW";        break;
        case CT_LEN_MISMATCH:        p = "LENGTH_MISMATCH";       break;
        case CT_OUT_OF_RAM:          p = "Out of RAM";            break;
        case CT_RANKING_ERR:         p = "RANKING_ERR";           break;
        case CT_ISOCOUNT_ERR:        p = "ISOCOUNT_ERR";          break;
        case CT_TAUCOUNT_ERR:        p = "TAUCOUNT_ERR";          break;
        case CT_ISOTAUCOUNT_ERR:     p = "ISOTAUCOUNT_ERR";       break;
        case CT_MAPCOUNT_ERR:        p = "MAPCOUNT_ERR";          break;
        case CT_TIMEOUT_ERR:         p = "Time limit exceeded";   break;
        case CT_ISO_H_ERR:           p = "ISO_H_ERR";             break;
        case CT_STEREOCOUNT_ERR:     p = "STEREOCOUNT_ERR";       break;
        case CT_ATOMCOUNT_ERR:       p = "ATOMCOUNT_ERR";         break;
        case CT_STEREOBOND_ERROR:    p = "STEREOBOND_ERR";        break;
        case CT_USER_QUIT_ERR:       p = "User requested termination"; break;
        case CT_REMOVE_STEREO_ERR:   p = "REMOVE_STEREO_ERR";     break;
        case CT_CALC_STEREO_ERR:     p = "CALC_STEREO_ERR";       break;
        case CT_STEREO_CANON_ERR:    p = "STEREO_CANON_ERR";      break;
        case CT_CANON_ERR:           p = "CANON_ERR";             break;
        case CT_WRONG_FORMULA:       p = "Wrong or missing chemical formula";  break;
        /*case CT_CANON_ERR2:          p = "CT_CANON_ERR2";         break;*/
        case CT_UNKNOWN_ERR:         p = "UNKNOWN_ERR";           break;
        case BNS_RADICAL_ERR:        p = "Cannot process free radical center"; break;
        case BNS_ALTBOND_ERR:        p = "Cannot process aromatic bonds";      break;

        default:
            if ( nErrorCode > CT_UNKNOWN_ERR ) {
                sprintf( szErrMsg, "No description(%d)", nErrorCode );
                p = szErrMsg;
            } else {
                sprintf( szErrMsg, "UNKNOWN_ERR(%d)", CT_UNKNOWN_ERR - nErrorCode );
                p = szErrMsg;
            }
            break;
    }
    return p;
}
/***********************************************************************************/
#ifndef COMPILE_ANSI_ONLY /* { */
/***********************************************************************************/
int SaveEquComponentsInfoAndSortOrder ( int iINChI, INCHI_SORT *pINChISort[TAUT_NUM], int *num_components,
                                        ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                                        COMP_ATOM_DATA composite_norm_data[TAUT_NUM+1],
                                        int bCompareComponents )
{
    int nRet = 0, i, k, nNumDeleted;
    /* equivalent components and sorting order */
    /* bCompareComponents: bit = 1 => compare */
    /*                     bit = 2 => compare non-isotopic */
    /*                     bit = 4 => compare non-tautomeric  */
    int bCompareIsotopic, bCompareTaut, bCompareAlt;
    ORIG_ATOM_DATA *inp_data = NULL;

    if ( num_components[iINChI] <= 1 )
        return 0;
#ifdef TARGET_LIB_FOR_WINCHI
    if ( !DRAWDATA )
        return 0;
#endif
    if ( !(bCompareComponents & CMP_COMPONENTS) )
        return 0;
    bCompareIsotopic = !(bCompareComponents & CMP_COMPONENTS_NONISO);
    bCompareTaut     =  (bCompareComponents & CMP_COMPONENTS_NONTAUT) ? TAUT_NON : TAUT_YES;
    bCompareAlt      =  ALT_TAUT(bCompareTaut);
    if ( num_components[iINChI] > 1 ) {
        if ( prep_inp_data[iINChI].bSavedInINCHI_LIB[iINChI] && prep_inp_data[iINChI].bPreprocessed[iINChI] ) {
            inp_data       = prep_inp_data+iINChI;
        } else
        if ( orig_inp_data->bSavedInINCHI_LIB[iINChI] && !orig_inp_data->bPreprocessed[iINChI] ) {
            inp_data       = orig_inp_data;
        } else {
            inp_data       = NULL;
        }
        if ( inp_data && !inp_data->nEquLabels && !prep_inp_data[iINChI].nSortedOrder ) {
            int i1, i2, nSet;
            AT_NUMB nAtNo;
            AT_NUMB nNumAtoms = (AT_NUMB)inp_data->num_inp_atoms;
            if ( (prep_inp_data[iINChI].nSortedOrder =
                                 (AT_NUMB *)inchi_calloc(num_components[iINChI]+1,
                                                    sizeof(prep_inp_data[0].nSortedOrder[0])))) {
                inp_data->nNumEquSets = 0;
                for ( i1 = 0, nSet = 0; i1 < num_components[iINChI]; i1 = i2 ) {
                    nNumDeleted = (pINChISort[bCompareTaut][i1].pINChI[bCompareTaut] && pINChISort[bCompareTaut][i1].pINChI[bCompareTaut]->bDeleted);
                    for ( i2 = i1+1; i2 < num_components[iINChI]; i2 ++ ) {
                        /* isotopic/non-isotopic comparison does not separate equivalent components */
                        if ( CompINChI2( pINChISort[bCompareTaut]+i1, pINChISort[bCompareTaut]+i2, bCompareTaut, bCompareIsotopic ) ) {
                            break;
                        } else {
                            nNumDeleted += (pINChISort[bCompareTaut][i2].pINChI[bCompareTaut] && pINChISort[bCompareTaut][i2].pINChI[bCompareTaut]->bDeleted);
                        }
                    }
                    if ( i2 - i1 - nNumDeleted > 1 ) {
                        if ( inp_data->nEquLabels ||
                             (inp_data->nEquLabels = (AT_NUMB *)inchi_calloc(inp_data->num_inp_atoms+1,
                                                 sizeof(inp_data->nEquLabels[0]))) ) {
                            nSet ++; /* found i2-i1 equivalent components && memory has been allocated */
                            for ( i = i1; i < i2; i ++ ) {
                                INChI_Aux *pINChI_Aux;
                                if (pINChISort[bCompareTaut][i].pINChI[bCompareTaut] && pINChISort[bCompareTaut][i].pINChI[bCompareTaut]->bDeleted)
                                    continue;
                                pINChI_Aux = (pINChISort[bCompareTaut][i].pINChI_Aux[bCompareTaut] &&
                                             pINChISort[bCompareTaut][i].pINChI_Aux[bCompareTaut]->nNumberOfAtoms)?
                                                pINChISort[bCompareTaut][i].pINChI_Aux[bCompareTaut]:
                                             (pINChISort[bCompareTaut][i].pINChI_Aux[bCompareAlt] &&
                                             pINChISort[bCompareTaut][i].pINChI_Aux[bCompareAlt]->nNumberOfAtoms)?
                                                pINChISort[bCompareTaut][i].pINChI_Aux[bCompareAlt]:
                                             (INChI_Aux *)NULL;
                                if ( pINChI_Aux && pINChI_Aux->nOrigAtNosInCanonOrd ) {
                                    for ( k = 0; k < pINChI_Aux->nNumberOfAtoms; k ++ ) {
                                        if ( (nAtNo = pINChI_Aux->nOrigAtNosInCanonOrd[k]) &&
                                              nAtNo <= nNumAtoms ) {
                                            inp_data->nEquLabels[nAtNo-1] = nSet;
                                        }
                                    }
                                }
                            }
                        } else {
                            return CT_OUT_OF_RAM;
                        }
                    }
                }
                nRet |= nSet? 1:0;
            } else {
                return CT_OUT_OF_RAM;
            }
            inp_data->nNumEquSets = nSet;
            /* output order */
            prep_inp_data[iINChI].nSortedOrder[0] = 0;
            for ( i1 = 0; i1 < num_components[iINChI]; i1 ++ ) {
                prep_inp_data[iINChI].nSortedOrder[i1+1] = pINChISort[TAUT_YES][i1].ord_number+1;
            }
#ifdef TARGET_LIB_FOR_WINCHI /* { */
            if ( DRAWDATA && GET_DRAWDATA && inp_data->nNumEquSets > 0 && inp_data->nEquLabels ) {
                int    nType = inp_data->bPreprocessed[iINChI]?
                                        COMPONENT_ORIGINAL_PREPROCESSED :
                                        COMPONENT_ORIGINAL;
                struct DrawData *pDrawData = GET_DRAWDATA( 0, nType, iINChI);
                if ( pDrawData && pDrawData->pWindowData && !pDrawData->pWindowData->nEquLabels ) {
                    /* copy equivalence data from inp_data to pDrawData->pWindowData */
                    if ( inp_data->nEquLabels &&
                         (pDrawData->pWindowData->nEquLabels = (AT_NUMB *)inchi_calloc(inp_data->num_inp_atoms,
                                                                          sizeof(inp_data->nEquLabels[0])))) {
                        memcpy( pDrawData->pWindowData->nEquLabels, inp_data->nEquLabels,
                                 inp_data->num_inp_atoms * sizeof(inp_data->nEquLabels[0]));
                        pDrawData->pWindowData->nNumEquSets  = inp_data->nNumEquSets;
                        pDrawData->pWindowData->nCurEquLabel = 0;
                    }
                }
            }
#endif  /* } TARGET_LIB_FOR_WINCHI */
        }
    }
    return nRet;
}

/************************************************************************************************/
int DisplayTheWholeCompositeStructure( INPUT_PARMS *ip, STRUCT_DATA *sd, long num_inp, int iINChI,
                                       PINChI2 *pINChI2, PINChI_Aux2 *pINChI_Aux2,
                                       ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                                       COMP_ATOM_DATA composite_norm_data[TAUT_NUM+1] )
{
    ORIG_ATOM_DATA *inp_data = NULL;
    int jj, j, k, err = 0, nNumIntermediateTaut = 0, bDisplayTaut;
    char szTitle[256];
    int nNumTautComponents, m;

    int bCompareIsotopic   = !(ip->bCompareComponents & CMP_COMPONENTS_NONISO);
    int bCompareTaut       =  (ip->bCompareComponents & CMP_COMPONENTS_NONTAUT) ? TAUT_NON : TAUT_YES;

    if ( ip->bCompareComponents & CMP_COMPONENTS ) {
        if ( prep_inp_data[iINChI].bSavedInINCHI_LIB[iINChI] && prep_inp_data[iINChI].bPreprocessed[iINChI] ) {
            inp_data       = prep_inp_data+iINChI;
        } else
        if ( orig_inp_data->bSavedInINCHI_LIB[iINChI] && !orig_inp_data->bPreprocessed[iINChI] ) {
            inp_data       = orig_inp_data;
        }
    }
    /**************************************************************************
     * display from one up to 4 structure pictures-results for all components *
     * Enable buttons:                                                        *
     * BN (non-tautomeric non-isotopic): inp_norm_data[0]->bExists            *
     * TN (tautomeric non-isotopic):     inp_norm_data[1]->bExists            *
     * BI (non-tautomeric isotopic):     inp_norm_data[0]->bExists &&         *
     *                                   inp_norm_data[0]->bHasIsotopicLayer  *
     * TI (tautomeric isotopic):         inp_norm_data[1]->bExists &&         *
     *                                   inp_norm_data[1]->bHasIsotopicLayer  *
     **************************************************************************/
    for ( jj = 0; ip->bDisplayCompositeResults && !sd->bUserQuitComponentDisplay && jj <= TAUT_INI; jj ++ ) {
    /*for ( j = 0; ip->bDisplayCompositeResults && !sd->bUserQuitComponentDisplay && j <= TAUT_INI; j ++ )*/
        j = (jj==0)? TAUT_NON : (jj==1)? TAUT_INI : (jj==2)? TAUT_YES : -1;
        if ( j < 0 )
            continue;
        if ( composite_norm_data[j].bExists && composite_norm_data[j].num_components > 1 ) {
            bDisplayTaut = (!(ip->nMode & REQ_MODE_BASIC) && !j)? -1 : j;
            nNumTautComponents = 0;
            if ( bDisplayTaut ) {
                /* find whether the structure is actually tautomeric */
                for ( m = 0; m < composite_norm_data[TAUT_YES].num_components; m ++ ) {
                    if ( !pINChI2[m][TAUT_YES] )
                        continue;
                    if ( pINChI2[m][TAUT_YES]->bDeleted || pINChI2[m][TAUT_YES]->lenTautomer > 0 )
                        nNumTautComponents ++;
                }
            }
            for ( k = 0; k <= composite_norm_data[j].bHasIsotopicLayer && !sd->bUserQuitComponentDisplay; k ++ ) {
                /*  added number of components, added another format for a single component case - DCh */
                int bMobileH = (bDisplayTaut>0 && nNumTautComponents);
                sprintf( szTitle, "%s Structure #%ld%s%s.%s%s%s%s%s",
                              j == TAUT_INI? "Preprocessed":"Result for", num_inp,
                              bMobileH? ", mobile H":
                              bDisplayTaut==0?", fixed H":"",
                              /*j? ", mobile H":", fixed H",*/
                              k? ", isotopic":"",
                              SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), iINChI? " (Reconnected)":"");
#ifndef TARGET_LIB_FOR_WINCHI
                /****** Display composite Result structure **************/
                nNumIntermediateTaut += (j == TAUT_INI );  /* display TAUT_INI (preprocessed) only once */
                if ( j != TAUT_INI || nNumIntermediateTaut == 1 ) {
                    err = DisplayCompositeStructure( composite_norm_data, j==TAUT_INI? 1:k /* bIsotopic*/,
                                                   j/*tautomeric*/,
                                                   j==TAUT_INI? NULL:pINChI2, j==TAUT_INI? NULL:pINChI_Aux2,
                                                   ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
                }
                if ( sd->bUserQuitComponentDisplay = (err==ESC_KEY) ) {
                    break;
                }

                if ( inp_data && inp_data->nEquLabels && inp_data->nNumEquSets && !sd->bUserQuitComponentDisplay &&
                     ((j == bCompareTaut || bCompareTaut && j == TAUT_INI) || 
                      bCompareTaut && !composite_norm_data[bCompareTaut].bExists) &&
                     (k == bCompareIsotopic || 
                      bCompareIsotopic && !composite_norm_data[j].bHasIsotopicLayer) ) {
                    AT_NUMB         nEquSet;
                    int             bDisplaySaved = ip->bDisplay;
                    /****** Display Equ Sets of composite Result structure **************/
                    for ( nEquSet = 1; nEquSet <= inp_data->nNumEquSets; nEquSet ++ ) {
                        sprintf( szTitle, "Equ set %d of %d, %s Structure #%ld%s%s.%s%s%s%s%s",
                                      nEquSet, inp_data->nNumEquSets,
                                      j == TAUT_INI? "Preprocessed":"Result for",
                                      num_inp,
                                      (bDisplayTaut>0 && nNumTautComponents)? ", mobile H": bDisplayTaut==0?", fixed H":"",
                                      /*j? ", mobile H":", fixed H",*/
                                      k? ", isotopic":"",
                                      SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), iINChI? " (Reconnected)":"");
                        ip->dp.nEquLabels   = inp_data->nEquLabels;
                        ip->dp.nCurEquLabel = nEquSet;
                        ip->dp.nNumEquSets  = inp_data->nNumEquSets;
                        ip->bDisplay = 1; /* force display if it was not requested */
                        err = DisplayCompositeStructure( composite_norm_data, k, j,
                                               pINChI2, pINChI_Aux2,
                                               ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
                        ip->dp.nEquLabels   = NULL;
                        ip->dp.nCurEquLabel = 0;
                        ip->dp.nNumEquSets  = 0;
                        ip->bDisplay = bDisplaySaved; /* restore display option */

                        if ( sd->bUserQuitComponentDisplay = (err==ESC_KEY) ) {
                            break;
                        }
                    }
                }
#else
                if(DRAWDATA && j <= TAUT_YES)
                {
                    struct DrawData vDrawData;
                    vDrawData.pWindowData = CreateWinDataComposite_( composite_norm_data, k, j,
                                                                     pINChI2, pINChI_Aux2,
                                                                     ip->bAbcNumbers, &ip->dp, ip->nMode);
                    /* vDrawData.pWindowData = CreateWinData_( composite_norm_data[j].at, composite_norm_data[j].num_at,
                                         k, j, pINChI[i], pINChI_Aux[i],ip->bAbcNumbers, &ip->dp, ip->nMode ); */
                    if( vDrawData.pWindowData != NULL )
                    {
                        int nType;
                        vDrawData.nComponent = 0;
                        if( j == 0 )
                            nType = (k == 0) ? COMPONENT_BN: COMPONENT_BI;
                        else
                            nType = (k == 0) ? COMPONENT_TN: COMPONENT_TI;
                        vDrawData.nType        = nType;
                           vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                        vDrawData.szTitle              = _strdup(szTitle);
                        vDrawData.pWindowData->szTitle = _strdup(szTitle);
                        if ( inp_data && inp_data->nEquLabels && inp_data->nNumEquSets &&
                             (j == bCompareTaut     || bCompareTaut     && !composite_norm_data[bCompareTaut].bExists) &&
                             (k == bCompareIsotopic || bCompareIsotopic && !composite_norm_data[j].bHasIsotopicLayer)  &&
                             (vDrawData.pWindowData->nEquLabels = (AT_NUMB *)inchi_calloc(inp_data->num_inp_atoms,
                                                                                  sizeof(inp_data->nEquLabels[0])))) {
                            memcpy( vDrawData.pWindowData->nEquLabels, inp_data->nEquLabels,
                                     inp_data->num_inp_atoms * sizeof(inp_data->nEquLabels[0]));
                            vDrawData.pWindowData->nNumEquSets  = inp_data->nNumEquSets;
                            vDrawData.pWindowData->nCurEquLabel = 0;
                        }
                        DRAWDATA(&vDrawData);
                    }
                } else
                if(DRAWDATA && GET_DRAWDATA && j == TAUT_INI)
                {
                    struct DrawData vDrawData;
                    struct DrawData *pDrawData;

                    if ( !(ip->bCompareComponents & CMP_COMPONENTS) ||
                           (ip->bCompareComponents & CMP_COMPONENTS_NONTAUT) ||
                          !k != !composite_norm_data[j].bHasIsotopicLayer ) {

                        continue;
                    }
                    /*
                    vDrawData.pWindowData = CreateWinDataComposite_( composite_norm_data, k, j,
                                                                     pINChI2, pINChI_Aux2,
                                                                     ip->bAbcNumbers, &ip->dp, ip->nMode);
                    */
                    vDrawData.pWindowData = CreateWinDataComposite_( composite_norm_data, 1 /*k*/, j,
                                                                     NULL, NULL,
                                                                     ip->bAbcNumbers, &ip->dp, ip->nMode);
                    if( vDrawData.pWindowData != NULL )
                    {
                        int nType = COMPONENT_ORIGINAL_PREPROCESSED;
                        pDrawData = GET_DRAWDATA( 0, nType, iINChI);
                        if ( pDrawData  ) {
                            FreeDrawData( pDrawData );
                            pDrawData->pWindowData = vDrawData.pWindowData;
                            vDrawData.pWindowData  = NULL;
                        } else {
                            pDrawData = &vDrawData;
                        }

                    /* vDrawData.pWindowData = CreateWinData_( composite_norm_data[j].at, composite_norm_data[j].num_at,
                                        k, j, pINChI[i], pINChI_Aux[i],ip->bAbcNumbers, &ip->dp, ip->nMode ); */
                        pDrawData->nComponent   = 0;
                        pDrawData->nType        = nType;
                           pDrawData->bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                        pDrawData->szTitle              = _strdup(szTitle);
                        pDrawData->pWindowData->szTitle = _strdup(szTitle);
                        if ( inp_data && inp_data->nEquLabels && inp_data->nNumEquSets &&
                             /*(j == bCompareTaut     || bCompareTaut     && !composite_norm_data[bCompareTaut].bExists) &&*/
                             /*(k == bCompareIsotopic || bCompareIsotopic && !composite_norm_data[j].bHasIsotopicLayer)  &&*/
                             (pDrawData->pWindowData->nEquLabels = (AT_NUMB *)inchi_calloc(inp_data->num_inp_atoms,
                                                                                  sizeof(inp_data->nEquLabels[0])))) {
                            memcpy( pDrawData->pWindowData->nEquLabels, inp_data->nEquLabels,
                                     inp_data->num_inp_atoms * sizeof(inp_data->nEquLabels[0]));
                            pDrawData->pWindowData->nNumEquSets  = inp_data->nNumEquSets;
                            pDrawData->pWindowData->nCurEquLabel = 0;
                        }
                        if ( pDrawData == &vDrawData ) {
                            DRAWDATA(pDrawData);  /* there was no prepocessed structure */
                        }
                    }
                }
#endif
            }
        }
    }
    return err;
}

#endif /*  }COMPILE_ANSI_ONLY */



/***********************************************************************************/
/* pINChI[INCHI_BAS] refers to either disconnected or original structure;            */
/*                  num_components[INCHI_BAS] > 0 if there was input structure      */
/***********************************************************************************/
/* pINChI[INCHI_REC] refers to the reconnected structure,                            */
/*                  and only if the input structure has been disconnected, that is,*/
/*                  num_components[INCHI_REC] > 0                                   */
/***********************************************************************************/
int SortAndPrintINChI(INCHI_IOSTREAM *output_file, 
                      char *pStr, int nStrLen, 
                      INCHI_IOSTREAM *log_file,
                      INPUT_PARMS *ip, 
                      ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                      COMP_ATOM_DATA composite_norm_data[INCHI_NUM][TAUT_NUM+1],
                      ORIG_STRUCT *pOrigStruct, int num_components[INCHI_NUM],
                      int num_non_taut[INCHI_NUM], int num_taut[INCHI_NUM],
                      INCHI_MODE bTautFlags[INCHI_NUM], INCHI_MODE bTautFlagsDone[INCHI_NUM],
                      NORM_CANON_FLAGS *pncFlags, long num_inp,
                      PINChI2 *pINChI[INCHI_NUM], 
                      PINChI_Aux2 *pINChI_Aux[INCHI_NUM], 
                      int *pSortPrintINChIFlags, unsigned char save_opt_bits)
{
    INCHI_SORT *pINChISort[INCHI_NUM][TAUT_NUM];
    int       j, i, k, k1, ret, iINChI, max_num_components;
    INCHI_MODE nMode;
    int       bDisconnectedCoord = (0 != (bTautFlagsDone[0] & TG_FLAG_DISCONNECT_COORD_DONE));
    int bINChIOutputOptions0, bCurOption, bINChIOutputOptionsCur, bEmbedReconnected, bAnnInXmlBrackets;
    static const char szAnnHdr[] = "InChI ANNOTATED CONTENTS";
    int ikflag = 0;

    ret = 1;
    for ( i = 0; i < INCHI_NUM; i ++ ) {
        for ( k = 0; k < TAUT_NUM; k ++ ) {
            bTautFlags[i]     |= pncFlags->bTautFlags[i][k];
            bTautFlagsDone[i] |= pncFlags->bTautFlagsDone[i][k];
        }
    }
    nMode = ip->nMode;
    if ( !(nMode & (REQ_MODE_BASIC|REQ_MODE_TAUT)) ) {
        nMode |= (REQ_MODE_BASIC|REQ_MODE_TAUT);
    }

    max_num_components = 0;
    for ( j = 0; j < INCHI_NUM; j ++ ) {
        if ( max_num_components < num_components[j] )
            max_num_components = num_components[j];
    }
    if ( max_num_components <= 0 )
        max_num_components = 1;

    for ( j = 0, i = 0; j < INCHI_NUM; j ++ ) {
        if ( num_components[j] ) {
            for ( k1 = 0; k1 < TAUT_NUM; k1 ++ ) {
                pINChISort[j][k1] = (INCHI_SORT *)inchi_calloc(max_num_components, sizeof(pINChISort[0][0][0]) );
                i += !pINChISort[j][k1]; /* number of failed allocatons */
            }
        } else {
            for ( k1 = 0; k1 < TAUT_NUM; k1 ++ ) {
                pINChISort[j][k1] = NULL; /* keep BC happy */
            }
        }
    }
    if ( i ) {
        ret = CT_OUT_OF_RAM;
        goto exit_function;
    }


    for ( j = 0; j < INCHI_NUM; j ++ ) {

        if ( !num_components[j] ) {
            continue;
        }

        iINChI = j;

#if ( OUTPUT_CONNECTED_METAL_ONLY == 1 ) /* test: output connected as the only one INChI */
        if ( INCHI_BAS == j && num_components[INCHI_REC] ) {
            j = INCHI_REC;
        }
#endif

        /*j = INCHI_BAS; <- for debug only */
        /* for only normal or disconnected coord compounds */
        /* (j=0=INCHI_BAS => normal or disconnected, j=1=INCHI_REC => reconnected */
        for ( k1 = 0; k1 < TAUT_NUM; k1 ++ ) {
            for ( i = 0; i < num_components[j]; i ++ ) {
                for ( k = 0; k < TAUT_NUM; k ++ ) {
                    pINChISort[j][k1][i].pINChI[k]     = pINChI[j][i][k];
                    pINChISort[j][k1][i].pINChI_Aux[k] = pINChI_Aux[j][i][k];
                }
                pINChISort[j][k1][i].ord_number = i;
            }
        }
        /* sort component INChIs */
        for ( k1 = 0; k1 < TAUT_NUM; k1 ++ ) {
            switch ( k1 ) {
            case TAUT_NON:
                qsort( pINChISort[j][k1], num_components[j], sizeof(pINChISort[0][0][0]), CompINChINonTaut2 );
                break;
            case TAUT_YES:
                qsort( pINChISort[j][k1], num_components[j], sizeof(pINChISort[0][0][0]), CompINChITaut2 );
                break;
            }
        }
#ifndef COMPILE_ANSI_ONLY
/* find equivalent and wINChI display order; use requested in ip->bCompareComponents comparison */
        ret = SaveEquComponentsInfoAndSortOrder ( iINChI, pINChISort[j], num_components, orig_inp_data, prep_inp_data,
#if ( FIX_DALKE_BUGS == 1 )
                                                  composite_norm_data? composite_norm_data[j]:NULL,
#else
                                                  composite_norm_data[j],
#endif
                                                  ip->bCompareComponents );
        if ( RETURNED_ERROR( ret ) ) {
            ret = 0;
            goto exit_function;
        } else {
            ret = 1;
        }
#endif
    }
    
    if ( !( ip->bINChIOutputOptions & INCHI_OUT_PRINT_OPTIONS ) ) {
        /* prepare InChI from the structures obtained by reversing InChI for returning to the caller */
        for ( j = 0; j < INCHI_NUM; j ++ ) {
            if ( !num_components[j] ) {
                continue;
            }
            /* pINChI[iINCHI][iComponent][bTaut] */
            /* j  = disconnected/connected */
            /* k1 = sort order for Mobile or Fixed H */
            k1 = TAUT_YES; /* in Mobile H order */
            /* store components in Mobile H order */
            
            for ( i = 0; i < num_components[j]; i ++ ) {

                if ( pINChISort[j][k1][i].pINChI[TAUT_NON] &&
                    !pINChISort[j][k1][i].pINChI[TAUT_YES] ) {
                    /* make sure Mobile-H is always present */
                    for ( k = 0; k < TAUT_NUM; k ++ ) {
                        pINChI[j][i][k]     = pINChISort[j][k1][i].pINChI[ALT_TAUT(k)];
                        pINChI_Aux[j][i][k] = pINChISort[j][k1][i].pINChI_Aux[ALT_TAUT(k)];
                    }
                } else {

                    for ( k = 0; k < TAUT_NUM; k ++ ) {
                        pINChI[j][i][k]     = pINChISort[j][k1][i].pINChI[k];
                        pINChI_Aux[j][i][k] = pINChISort[j][k1][i].pINChI_Aux[k];
                    }
                }
            }
        }

    } else {
        
        /* print inchi string(s) */


        bINChIOutputOptions0 = ip->bINChIOutputOptions & ~INCHI_OUT_PRINT_OPTIONS;

        bEmbedReconnected    = ip->bINChIOutputOptions & INCHI_OUT_EMBED_REC;

        for ( i = 0; i < 4; i ++ ) {
            switch( i ) {
            case 0:
                bCurOption = INCHI_OUT_XML;
                break;
            case 1:
                bCurOption = INCHI_OUT_PLAIN_TEXT;
                break;
            case 2:
                bCurOption = INCHI_OUT_PLAIN_TEXT_COMMENTS;
                break;
            case 3:
                bCurOption = INCHI_OUT_XML_TEXT_COMMENTS;
                break;
            default:
                continue;
            }
            if ( ip->bINChIOutputOptions & bCurOption ) {
                bAnnInXmlBrackets = 0;
                if ( i == 1 ) {
                    ;/*bEmbedReconnected = 0;*/
                }
                if ( i == 3 ) {
                    bCurOption = INCHI_OUT_XML; /* xml output as annotation */
                }
                bINChIOutputOptionsCur = bINChIOutputOptions0 | bCurOption;
                switch ( i ) {
                case 0:
                case 1:
                    /* output INChI */
                    bINChIOutputOptionsCur |= bEmbedReconnected;
                    break;
                case 2:
                case 3:
                    /* output annotation */
                    bAnnInXmlBrackets = (i == 2 && (ip->bINChIOutputOptions & INCHI_OUT_XML ));
                    if ( bAnnInXmlBrackets ) {
                        inchi_ios_print( output_file, "\n<%s>\n", szAnnHdr );
                    } else {
                        inchi_ios_print( output_file, "\n==== %s ====\n", szAnnHdr );
                    }
                    bINChIOutputOptionsCur |= bEmbedReconnected;
                    bINChIOutputOptionsCur &= ~INCHI_OUT_TABBED_OUTPUT;
                    break;
                default:
                    continue;
                }

#ifdef TARGET_LIB_FOR_WINCHI
                if ( ikflag==0 )
                    output_file->type = INCHI_IOSTREAM_STRING;
#endif


                ret &= OutputINChI2(pStr, nStrLen, 
                                    pINChISort, INCHI_BAS /*iINChI*/, 
                                    pOrigStruct,
                                    bDisconnectedCoord, OUT_TN, 
                                    bINChIOutputOptionsCur, 
                                    0 != (bINChIOutputOptionsCur & INCHI_OUT_XML),
                                    ip->bAbcNumbers, ip->bCtPredecessors, 
                                    ip->bNoStructLabels, 
                                    num_components, num_non_taut, num_taut,
                                    output_file, log_file, num_inp,
                                    ip->pSdfLabel,ip->pSdfValue, ip->lSdfId, 
                                    pSortPrintINChIFlags,
                                    save_opt_bits);

                if ( ret &&  !(bINChIOutputOptionsCur & INCHI_OUT_EMBED_REC) ) 
                {
                    ret &= OutputINChI2(pStr, nStrLen, pINChISort, INCHI_REC /*iINChI*/, 
                                        pOrigStruct,
                                        bDisconnectedCoord, OUT_TN, 
                                        bINChIOutputOptionsCur, 
                                        0 != (bINChIOutputOptionsCur & INCHI_OUT_XML),
                                        ip->bAbcNumbers, ip->bCtPredecessors, 
                                        ip->bNoStructLabels, 
                                        num_components, num_non_taut, num_taut,
                                        output_file, log_file, num_inp,
                                        ip->pSdfLabel,ip->pSdfValue, ip->lSdfId, 
                                        pSortPrintINChIFlags,
                                        save_opt_bits);
                }

#ifdef TARGET_LIB_FOR_WINCHI
                /* always calculate InChIKey */
                ikflag++;
                if (ikflag==1) 
                {
                    if (ret)
                    {                    
                        char ik_string[256];    /*^^^ Resulting InChIKey string */
                        int ik_ret=0;           /*^^^ InChIKey-calc result code */
                        int xhash1, xhash2;
                        char szXtra1[256], szXtra2[256];
                        size_t slen = output_file->s.nUsedLength;
                        char *buf = NULL;
                        extract_inchi_substring(&buf, output_file->s.pStr, slen);            
                        inchi_ios_flush(output_file);
                        output_file->type = INCHI_IOSTREAM_FILE;
                        /* calculate and print InChIKey */
                        if (NULL!=buf)
                        {
                            xhash1 = xhash2 = 0;
                            if ( ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1 ) ||
                                 ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2 ) )
                                 xhash1 = 1;
                            if ( ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA2 ) ||
                                 ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2 ) )
                                 xhash2 = 1;                
                            ik_ret = GetINCHIKeyFromINCHI(buf, xhash1, xhash2, ik_string, szXtra1, szXtra2);
                            inchi_free(buf);
                        }
                        else
                            ik_ret = 3;                     
    
                        if (ik_ret==INCHIKEY_OK)   
                        {
                            /* NB: correctly treat tabbed output with InChIKey & hash extensions */                
                            char csep = '\n';
                                if ( ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT ) 
                            csep = '\t';
                            inchi_ios_print(output_file, "InChIKey=%-s",ik_string);
                            if ( xhash1 )
                                inchi_ios_print(output_file, "%cXHash1=%-s",csep,szXtra1);
                            if ( xhash2 )
                                inchi_ios_print(output_file, "%cXHash2=%-s",csep,szXtra2);
                            inchi_ios_print(output_file, "\n");
                        }
                        else            
                        {
                            inchi_ios_print(log_file, "Warning (Could not compute InChIKey: ", num_inp);    
                        }
                            
                        /*inchi_ios_flush(output_file);
                        inchi_ios_flush2(log_file, stderr);*/
                    }
                    else
                    {
                        inchi_ios_flush(output_file);
                        output_file->type = INCHI_IOSTREAM_FILE;
                    }

                }

#endif

                if ( bAnnInXmlBrackets ) {
                    inchi_ios_print( output_file, "</%s>\n\n", szAnnHdr );
                }
                if ( !ret ) {
                    break;
                }
            }
        }
    }

exit_function:
    for ( j = 0; j < INCHI_NUM; j ++ ) {
        for ( k1 = 0, i = 0; k1 < TAUT_NUM; k1 ++ ) {
            if ( pINChISort[j][k1] ) {
                inchi_free( pINChISort[j][k1] );
            }
        }
    }
    ret = ret? 0 : _IS_FATAL;



    return ret;
}
/**********************************************************************************/
void FreeAllINChIArrays( PINChI2 *pINChI[INCHI_NUM], PINChI_Aux2 *pINChI_Aux[INCHI_NUM], int num_components[INCHI_NUM] )
{
    int k;
    for ( k = 0; k < INCHI_NUM; k ++ ) {
        FreeINChIArrays( pINChI[k], pINChI_Aux[k], num_components[k] );
        num_components[k] = 0;
        if ( pINChI[k] ) {
            inchi_free( pINChI[k] );
            pINChI[k] = NULL;
        }
        if ( pINChI_Aux[k] ) {
            inchi_free( pINChI_Aux[k] );
            pINChI_Aux[k] = NULL;
        }

    }
}
/**********************************************************************************/
void FreeINChIArrays( PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, int num_components )
{
    int i, k;
    /* release allocated memory */
    if ( pINChI ) {
        for ( i = 0; i < num_components; i ++ ) {
            for ( k = 0; k < TAUT_NUM; k ++ ) {
                Free_INChI( &pINChI[i][k] );
                /*
                inchi_free( pINChI[i][k] );
                pINChI[i][k] = NULL;
                */
            }
        }
    }
    if ( pINChI_Aux ) {
        for ( i = 0; i < num_components; i ++ ) {
            for ( k = 0; k < TAUT_NUM; k ++ ) {
                Free_INChI_Aux( &pINChI_Aux[i][k] );
                /*
                inchi_free( pINChI_Aux[i][k] );
                pINChI_Aux[i][k] = NULL;
                */
            }
        }
    }
}


/**********************************************
 * output " L=V" or " L missing" or ""
 * The fprintf format string must contain %s%s%s%s
 */

const char gsMissing[] = "is missing";
const char gsEmpty[]   = "";
const char gsSpace[]   = " ";
const char gsEqual[]   = "=";

#ifndef TARGET_API_LIB
/*********************************************************************************************************/
void SplitTime( unsigned long ulTotalTime, int *hours, int *minutes, int *seconds, int *mseconds )
{
        *mseconds = (int)(ulTotalTime % 1000);
        ulTotalTime /= 1000;
        *seconds = (int)(ulTotalTime % 60);
        ulTotalTime /= 60;
        *minutes = (int)(ulTotalTime % 60);
        ulTotalTime /= 60;
        *hours = (int)(ulTotalTime);
}
/*********************************************************************************************************/
int ReadTheStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_IOSTREAM  *inp_file, ORIG_ATOM_DATA *orig_inp_data,
                      /* for CML:*/ int inp_index, int *out_index )
{
    inchiTime     ulTStart;
    int           nRet = 0, nRet2 = 0;
    int           bGetOrigCoord = !(ip->bINChIOutputOptions & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO));
    INCHI_MODE InpAtomFlags = 0;  /* reading Molfile may set FLAG_INP_AT_CHIRAL bit */

    /* vABParityUnknown holds actual value of an internal constant signifying       */
    /* unknown parity: either the same as for undefined parity (default==standard)  */
    /*  or a specific one (non-std; requested by SLUUD switch).                     */
    int vABParityUnknown = AB_PARITY_UNDF;
    if ( 0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO) ) 
    {
        /* Make labels for unknown and undefined stereo different */
        vABParityUnknown = AB_PARITY_UNKN;
    }

    memset( sd, 0, sizeof(*sd) );
    switch ( ip->nInputType ) {
    case INPUT_MOLFILE:
    case INPUT_SDFILE:
        if ( orig_inp_data ) {
            if ( ip->pSdfValue && ip->pSdfValue[0] ) {
                /* Added 07-29-2003 to avoid inheriting exact value from prev. structure
                   and to make reference to a (bad) structure with unknown ID Value */
                char *p, *q;  /* q shadows prev declaration of const char *q */
                int  n;
                if ( (p = strrchr( ip->pSdfValue, '+' )) &&
                     '[' == *(p-1) && 0 < (n=strtol(p+1,&q,10)) && q[0] && ']'==q[0] && !q[1] ) {
                    sprintf( p+1, "%d]", n+1 );
                } else {
                    strcat( ip->pSdfValue, " [+1]" );
                }
            }
            InchiTimeGet( &ulTStart );
            sd->fPtrStart = (inp_file->f == stdin)? -1 : ftell( inp_file->f );
            /*  read the original structure */
            nRet2 = MolfileToOrigAtom( inp_file->f, orig_inp_data, ip->bMergeAllInputStructures, bGetOrigCoord, ip->bDoNotAddH,
                               ip->pSdfLabel, ip->pSdfValue, &ip->lSdfId, &ip->lMolfileNumber,
                               &InpAtomFlags, &sd->nStructReadError, sd->pStrErrStruct );


            if ( !ip->bGetSdfileId || ip->lSdfId == 999999) ip->lSdfId = 0;
            if ( !ip->bGetMolfileNumber || ip->lMolfileNumber < 0 ) ip->lMolfileNumber = 0;
            sd->fPtrEnd = (inp_file->f == stdin)? -1 : ftell( inp_file->f );
            sd->ulStructTime += InchiTimeElapsed( &ulTStart );
#if ( bRELEASE_VERSION == 0 )
            sd->bExtract |= orig_inp_data->bExtract;
#endif
            /* 2004-11-16: added Molfile Chiral Flag Mode */
            /* *****************************************************************************
             * Chiral flags are set in: 
             * - RunICHI.c #1610 -- ReadTheStructure()     -- cInChI, wInChI (here)
             * - e_IchiMain.c #273 -- main()               -- C example of calling InChI dll
             * - inchi_dll.c  #1662 -- ExtractOneStructure -- InChI dll code 
             *******************************************************************************/   
            /* 1. Highest precedence: Chiral Flag set by the user */
            if ( ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL ) {
                InpAtomFlags = FLAG_INP_AT_CHIRAL; /* forced by the user */
            } else
            if ( ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL ) {
                InpAtomFlags = FLAG_INP_AT_NONCHIRAL; /* forced by the user */
            } else
            if ( (InpAtomFlags & FLAG_INP_AT_CHIRAL) && (InpAtomFlags && FLAG_INP_AT_NONCHIRAL) ) {
                InpAtomFlags &= ~FLAG_INP_AT_NONCHIRAL;
            }
            /* save requested flags in the AuxInfo */
            sd->bChiralFlag &= ~( FLAG_INP_AT_CHIRAL | FLAG_INP_AT_NONCHIRAL );
            sd->bChiralFlag |= InpAtomFlags & ( FLAG_INP_AT_CHIRAL | FLAG_INP_AT_NONCHIRAL );
            /* quick fix: modify ip->nMode on the fly */
            /* 2. The user requested both Stereo AND Chiral flag */
            if ( (ip->nMode & REQ_MODE_CHIR_FLG_STEREO) && (ip->nMode & REQ_MODE_STEREO) ) {
                if ( InpAtomFlags & FLAG_INP_AT_CHIRAL ) {
                    /* structure has chiral flag or the user said it is chiral */
                    ip->nMode &= ~(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO);
                    sd->bChiralFlag |= FLAG_INP_AT_CHIRAL; /* write AuxInfo as chiral */
                } else {
                    ip->nMode &= ~REQ_MODE_RACEMIC_STEREO;
                    ip->nMode |=  REQ_MODE_RELATIVE_STEREO;
                    sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL; /* write AuxInfo as explicitly not chiral */
                }
            }
        } else {
            /*  read the next original structure */
            int           nStructReadError=0;
            if ( !ip->bMergeAllInputStructures ) {
                nRet2 = MolfileToOrigAtom( inp_file->f, NULL, 0, 0, 0,
                                       NULL, NULL, NULL, NULL,
                                       NULL, &nStructReadError, NULL );
                if ( nRet2 <= 0 && 10 < nStructReadError && nStructReadError < 20 ) {
                    return _IS_EOF;
                }
            } else {
                return _IS_EOF;
            }
        }
        break;
    case INPUT_INCHI_XML:
    case INPUT_INCHI_PLAIN:
        if ( orig_inp_data ) {
            if ( ip->pSdfValue && ip->pSdfValue[0] ) {
                /* Added 07-29-2003 to avoid inheriting exact value from prev. structure
                   and to make reference to a (bad) structure with unknown ID Value */
                char *p, *q;
                int  n;
                if ( (p = strrchr( ip->pSdfValue, '+' )) &&
                     '[' == *(p-1) && 0 < (n=strtol(p+1,&q,10)) && q[0] && ']'==q[0] && !q[1] ) {
                    sprintf( p+1, "%d]", n+1 );
                } else {
                    strcat( ip->pSdfValue, " [+1]" );
                }
            }
            InchiTimeGet( &ulTStart );
            sd->fPtrStart = (inp_file->f == stdin)? -1 : ftell( inp_file->f );
            /*  read the original structure */
            nRet2 = INChIToOrigAtom( inp_file, orig_inp_data, ip->bMergeAllInputStructures,
                               bGetOrigCoord, ip->bDoNotAddH, vABParityUnknown,
                               ip->nInputType, ip->pSdfLabel, ip->pSdfValue, &ip->lMolfileNumber,
                               &InpAtomFlags, &sd->nStructReadError, sd->pStrErrStruct );
            /*if ( !ip->bGetSdfileId || ip->lSdfId == 999999) ip->lSdfId = 0;*/
            sd->fPtrEnd = (inp_file->f == stdin)? -1 : ftell( inp_file->f );

            sd->ulStructTime += InchiTimeElapsed( &ulTStart );
#if ( bRELEASE_VERSION == 0 )
            sd->bExtract |= orig_inp_data->bExtract;
#endif
            /* 2004-11-16: added Molfile Chiral Flag Mode */
            if ( ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL ) {
                InpAtomFlags = FLAG_INP_AT_CHIRAL; /* forced by the user */
            } else
            if ( ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL ) {
                InpAtomFlags = FLAG_INP_AT_NONCHIRAL; /* forced by the user */
            } else
            if ( (InpAtomFlags & FLAG_INP_AT_CHIRAL) && (InpAtomFlags && FLAG_INP_AT_NONCHIRAL) ) {
                InpAtomFlags &= ~FLAG_INP_AT_NONCHIRAL;
            }
            sd->bChiralFlag |= InpAtomFlags; /* copy chiral flag to AuxInfo */
            /* quick fix: modify ip->nMode on the fly */
            if ( (ip->nMode & REQ_MODE_CHIR_FLG_STEREO) && (ip->nMode & REQ_MODE_STEREO) ) {
                if ( InpAtomFlags & FLAG_INP_AT_CHIRAL ) {
                    ip->nMode &= ~(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO);
                } else {
                    ip->nMode &= ~REQ_MODE_RACEMIC_STEREO;
                    ip->nMode |=  REQ_MODE_RELATIVE_STEREO;
                }
            }
        } else {
            /*  read the next original structure */
            int           nStructReadError=0;
            if ( !ip->bMergeAllInputStructures ) {
                nRet2 = INChIToOrigAtom( inp_file, NULL, 0, 0, 0, 0, 
                                    ip->nInputType, NULL, NULL, NULL, NULL, &nStructReadError, NULL );
                if ( nRet2 <= 0 && 10 < nStructReadError && nStructReadError < 20 ) {
                    return _IS_EOF;
                }
            } else {
                return _IS_EOF;
            }
        }
        break;

#if ( ADD_CMLPP == 1 )
    /* BILLY 8/6/04 */
    case INPUT_CMLFILE:
        if ( orig_inp_data ) {

            InchiTimeGet( &ulTStart );
            /*
            if ( inp_index >= 0 ) {
                sd->fPtrStart = inp_index;
            } else {
                sd->fPtrStart = GetCmlStructIndex();
            }
            */
            sd->fPtrStart = -1; /* disable "CopyMOLfile() for CML input files */
            sd->fPtrEnd = -1;
            /*  read the original structure */
            nRet = CmlfileToOrigAtom( inp_file->f, orig_inp_data, ip->bMergeAllInputStructures,
                               bGetOrigCoord, ip->bDoNotAddH, inp_index, out_index,
                               ip->pSdfLabel, ip->pSdfValue, &ip->lSdfId,
                               &sd->nStructReadError, sd->pStrErrStruct );


            sd->ulStructTime += InchiTimeElapsed( &ulTStart );
#if ( bRELEASE_VERSION == 0 )
            sd->bExtract |= orig_inp_data->bExtract;
#endif
        } else {
            /*  read the next original structure */
            int nStructReadError=0;
            if ( !ip->bMergeAllInputStructures ) {
                nRet2 = CmlfileToOrigAtom( inp_file->f, NULL, 0, 0, 0, inp_index, out_index,
                                       NULL, NULL, NULL, &nStructReadError, NULL );

                if ( nRet2 <= 0 && 10 < nStructReadError && nStructReadError < 20 ) {
                    return _IS_EOF;
                }
            } else {
                return _IS_EOF;
            }
        }
        break;
#endif

    default:
        nRet = _IS_FATAL; /*  wrong file type */
    }
    return nRet;
}
#endif
/*****************************************************************************************************/
int TreatReadTheStructureErrors(  STRUCT_DATA *sd, INPUT_PARMS *ip, int nLogMask,
                                  INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file, 
                                  INCHI_IOSTREAM *prb_file, /*^^^ was: INCHI_IOSTREAM */
                                  ORIG_ATOM_DATA *orig_inp_data, long *num_inp, char *pStr, int nStrLen )
{
    int nRet = _IS_OKAY;
    /*  End of file */
    if ( 10 < sd->nStructReadError && sd->nStructReadError < 20 ) {
        if ( sd->pStrErrStruct[0] ) {
            inchi_ios_eprint( log_file, "%s inp structure #%ld: End of file.%s%s%s%s    \n", sd->pStrErrStruct, *num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
        }
        inchi_ios_eprint( log_file, "End of file detected after structure #%ld.   \n", *num_inp-1 );
        nRet = _IS_EOF;
        goto exit_function; /*  end of file */
    }

    /*(*num_inp) ++;*/

    /*  Skipping the structures */
    if ( *num_inp < ip->first_struct_number ) {

#if ( !defined(TARGET_API_LIB) && !defined(TARGET_EXE_STANDALONE) )
/*^^^ #ifndef TARGET_API_LIB */
        if ( log_file->f != stderr ) {
            inchi_fprintf( stderr, "\rSkipping structure #%ld.%s%s%s%s...", *num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue));
        }
#endif
        nRet = sd->nErrorType = _IS_SKIP;
        goto exit_function;
    }

    sd->nErrorType = GetInpStructErrorType( ip, sd->nStructReadError, sd->pStrErrStruct, orig_inp_data->num_inp_atoms );

    /*  init xml output */
    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && !ip->bXmlStarted ) {
        OutputINChIXmlRootStartTag( output_file );
        ip->bXmlStarted ++;
    }
    /*  init xml structure block */
    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && !sd->bXmlStructStarted ) {
        if ( !OutputINChIXmlStructStartTag( output_file, pStr, 1, nStrLen, ip->bNoStructLabels,
                                           *num_inp, ip->pSdfLabel, ip->pSdfValue ) ) {
            inchi_ios_eprint( log_file, "Cannot create start xml tag for structure #%ld.%s%s%s%s Terminating.\n", *num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            sd->bXmlStructStarted = -1;
            nRet = _IS_FATAL;
            goto exit_function;
        }
        sd->bXmlStructStarted ++;
    }

    /*  Fatal error */
    if ( sd->nErrorType == _IS_FATAL ) {
        if ( nLogMask & LOG_MASK_FATAL )
            inchi_ios_eprint( log_file, "Fatal Error %d (aborted; %s) inp structure #%ld.%s%s%s%s\n",
                    sd->nStructReadError, sd->pStrErrStruct, *num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
#if ( bRELEASE_VERSION == 1 || EXTR_FLAGS == 0 )
        if ( prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && !ip->bSaveAllGoodStructsAsProblem ) {
            CopyMOLfile(inp_file->f, sd->fPtrStart, sd->fPtrEnd, prb_file->f, *num_inp);
        }
#endif
        /* goto exit_function; */
    }
    /*  Non-fatal errors: do not produce INChI */
    if ( sd->nErrorType == _IS_ERROR ) {  /*  70 => too many atoms */
        if ( nLogMask & LOG_MASK_ERR )
            inchi_ios_eprint( log_file, "Error %d (no %s; %s) inp structure #%ld.%s%s%s%s\n",
                    sd->nStructReadError, (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY)?"Molfile":INCHI_NAME,
                    sd->pStrErrStruct, *num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
#if ( bRELEASE_VERSION == 1 || EXTR_FLAGS == 0 )
        if ( prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && !ip->bSaveAllGoodStructsAsProblem) {
            CopyMOLfile(inp_file->f, sd->fPtrStart, sd->fPtrEnd, prb_file->f, *num_inp);
        }
#endif
    }

    /*  Warnings: try to produce INChI */
    if ( sd->nErrorType == _IS_WARNING ) {
        if ( nLogMask & LOG_MASK_WARN )
            inchi_ios_eprint( log_file, "Warning: (%s) inp structure #%ld.%s%s%s%s\n",
                    sd->pStrErrStruct, *num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
    }

    /*  xml error/warning processing; close xml struct block if error */
    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML)
#ifdef TARGET_LIB_FOR_WINCHI
         || (ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW) && (ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT)
#endif
        ) {
        if ( sd->nErrorType != _IS_OKAY && sd->nErrorType != _IS_WARNING ) {
            sd->nErrorType =
                ProcessStructError( output_file, log_file, /*sd->nStructReadError,*/
                         sd->pStrErrStruct, sd->nErrorType, &sd->bXmlStructStarted, *num_inp, ip, pStr, nStrLen );
        }
    }
exit_function:
    if ( nRet <= _IS_OKAY && sd->nErrorType > 0 ) {
        nRet = sd->nErrorType;
    }
    return nRet;
}
/******************************************************************************************************/
int GetOneComponent( STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file,
                     INP_ATOM_DATA *inp_cur_data,
                     ORIG_ATOM_DATA *orig_inp_data, int i, long num_inp, char *pStr, int nStrLen )
{
    inchiTime ulTStart;
    InchiTimeGet( &ulTStart );
    CreateInpAtomData( inp_cur_data, orig_inp_data->nCurAtLen[i], 0 );
    inp_cur_data->num_at = ExtractConnectedComponent( orig_inp_data->at, orig_inp_data->num_inp_atoms, i+1, inp_cur_data->at );
    sd->ulStructTime += InchiTimeElapsed( &ulTStart );

    /*  error processing */
    if ( inp_cur_data->num_at <= 0 || orig_inp_data->nCurAtLen[i] != inp_cur_data->num_at ) {
        /*  log error message */
        AddMOLfileError(sd->pStrErrStruct, "Cannot extract Component");
        inchi_ios_eprint( log_file, "%s #%d structure #%ld.%s%s%s%s\n", sd->pStrErrStruct, i+1, num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue));
        sd->nErrorCode = inp_cur_data->num_at < 0? inp_cur_data->num_at : (orig_inp_data->nCurAtLen[i] != inp_cur_data->num_at)? CT_ATOMCOUNT_ERR : CT_UNKNOWN_ERR;
        /* num_err ++; */
        sd->nErrorType = _IS_ERROR;
        if ( (ip->bINChIOutputOptions & INCHI_OUT_XML)
#ifdef TARGET_LIB_FOR_WINCHI
             || (ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW) && (ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT)
#endif
           ) {
            /*  xml error message */
            sd->nErrorType = ProcessStructError( output_file, log_file, /*sd->nErrorCode,*/ sd->pStrErrStruct,
                                            sd->nErrorType, &sd->bXmlStructStarted, num_inp, ip, pStr, nStrLen );
        }
    }
    return sd->nErrorType;
}
/*******************************************************************************************/
int GetProcessingWarningsOneINChI(INChI *pINChI, INP_ATOM_DATA *inp_norm_data, char *pStrErrStruct)
{
    int j;
    int nAmbiguousStereoAtoms, nAmbiguousStereoBonds;
    nAmbiguousStereoAtoms = 0;
    nAmbiguousStereoBonds = 0;

    if ( inp_norm_data->at ) {
        for ( j = 0; j < pINChI->nNumberOfAtoms; j ++ ) {
            if ( inp_norm_data->at[j].bAmbiguousStereo & (AMBIGUOUS_STEREO_ATOM | AMBIGUOUS_STEREO_ATOM_ISO) ) {
                nAmbiguousStereoAtoms ++;
            }
            if ( inp_norm_data->at[j].bAmbiguousStereo & (AMBIGUOUS_STEREO_BOND | AMBIGUOUS_STEREO_BOND_ISO) ) {
                nAmbiguousStereoBonds ++;
            }
        }
        if ( nAmbiguousStereoAtoms ) {
            AddMOLfileError(pStrErrStruct, "Ambiguous stereo:");
            AddMOLfileError(pStrErrStruct, "center(s)");
        }
        if ( nAmbiguousStereoBonds ) {
            AddMOLfileError(pStrErrStruct, "Ambiguous stereo:");
            AddMOLfileError(pStrErrStruct, "bond(s)");
        }
    }
    return (nAmbiguousStereoAtoms || nAmbiguousStereoBonds);
}
/*******************************************************************************************/
int  GetProcessingWarnings(INChI *cur_INChI[], INP_ATOM_DATA **inp_norm_data, STRUCT_DATA *sd)
{
    int i, ret = 0;
    for (i = 0; i < TAUT_NUM; i ++ ) {
        if ( cur_INChI[i] && cur_INChI[i]->nNumberOfAtoms>0 ) {
            ret |= GetProcessingWarningsOneINChI(cur_INChI[i], inp_norm_data[i], sd->pStrErrStruct);
        }
    }
    return ret;
}

/*******************************************************************************************/
int CreateOneComponentINChI( STRUCT_DATA *sd, INPUT_PARMS *ip, INP_ATOM_DATA *inp_cur_data, ORIG_ATOM_DATA *orig_inp_data,
                            PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, int iINChI,
                            int i, long num_inp, INP_ATOM_DATA **inp_norm_data, NORM_CANON_FLAGS *pncFlags, 
                            INCHI_IOSTREAM *log_file )
{
    inchiTime     ulTStart, ulTEnd, *pulTEnd = NULL;
    int           k, num_at, ret = 0;
    int           bOrigCoord;
    INCHI_MODE     bTautFlags     = ip->bTautFlags;
    INCHI_MODE     bTautFlagsDone = (ip->bTautFlagsDone | sd->bTautFlagsDone[INCHI_BAS]);
    INChI       *cur_INChI[TAUT_NUM];
    INChI_Aux   *cur_INChI_Aux[TAUT_NUM];
    long          lElapsedTime;
    /*
    PINChI2     *pINChI     = pINChI2[iINChI];
    PINChI_Aux2 *pINChI_Aux = pINChI_Aux2[iINChI];
    */
    InchiTimeGet( &ulTStart );
    bOrigCoord = !(ip->bINChIOutputOptions & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO));

    for ( k = 0; k < TAUT_NUM; k ++ ) {
        cur_INChI[k]      = NULL;
        cur_INChI_Aux[k]  = NULL;
    }
    /*  allocate memory for non-tautimeric (k=0) and tautomeric (k=1) results */
    for ( k = 0; k < TAUT_NUM; k ++ ) {
        int nAllocMode = (k==TAUT_YES? REQ_MODE_TAUT:0) |
                         (bTautFlagsDone & ( TG_FLAG_FOUND_ISOTOPIC_H_DONE |
                                             TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE ))?
                         (ip->nMode & REQ_MODE_ISO):0;

        if ( (k==TAUT_NON && (ip->nMode & REQ_MODE_BASIC )) ||
             (k==TAUT_YES && (ip->nMode & REQ_MODE_TAUT ))     ) {
            /*  alloc INChI and INChI_Aux */
            cur_INChI[k]     = Alloc_INChI( inp_cur_data->at, inp_cur_data->num_at, &inp_cur_data->num_bonds,
                                          &inp_cur_data->num_isotopic, nAllocMode );
            cur_INChI_Aux[k] = Alloc_INChI_Aux( inp_cur_data->num_at,
                                          inp_cur_data->num_isotopic, nAllocMode, bOrigCoord );
            if ( cur_INChI_Aux[k] ) {
                cur_INChI_Aux[k]->bIsIsotopic = inp_cur_data->num_isotopic;
            }
            /*  alloc memory for the output structure: non-tautomeric and tautomeric (for displaying) */
            CreateInpAtomData( inp_norm_data[k], inp_cur_data->num_at, k );
        } else {
            FreeInpAtomData( inp_norm_data[k] );
        }
    }
    lElapsedTime = InchiTimeElapsed( &ulTStart );
    if ( ip->msec_MaxTime ) {
        ip->msec_LeftTime -= lElapsedTime;
    }
    sd->ulStructTime += lElapsedTime;


/*^^^#if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined( TARGET_API_LIB ) ) */
#if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined( TARGET_API_LIB ) && !defined(TARGET_EXE_STANDALONE) )
#if ( TEST_RENUMB_ATOMS != 1 )
    /*  log file / console output */
    if ( log_file->f && log_file->f != stderr ) { /* NULL log_file now ignored. 11-23-2005 */
        if ( ip->bDisplay )
            inchi_ios_eprint( log_file, "Component #%d structure #%ld.%s%s%s%s...\n", i+1, num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
        else
            inchi_fprintf( stderr, "Component #%d structure #%ld.%s%s%s%s...\r", i+1, num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
    }
#endif
#endif
    /******************************************************
     *
     *  Get one component canonical numberings, etc.
     *
     ******************************************************/

    /*
     * Create_INChI() return value:
     * num_at <= 0: error code
     * num_at >  0: number of atoms (excluding terminal hydrogen atoms)
     * inp_norm_data[0] => non-tautomeric, inp_norm_data[1] => tautomeric
     */
    InchiTimeGet( &ulTStart );
    if ( ip->msec_MaxTime ) {
        ulTEnd = ulTStart;
        pulTEnd = &ulTEnd;
        if ( ip->msec_LeftTime > 0 ) {
            InchiTimeAddMsec( pulTEnd, ip->msec_LeftTime );
        }
    }
    num_at = Create_INChI( cur_INChI, cur_INChI_Aux, orig_inp_data/* not used */, inp_cur_data->at,
                          inp_norm_data,
                          inp_cur_data->num_at,
                          ip->nMode, &bTautFlags, &bTautFlagsDone, pulTEnd, NULL, sd->pStrErrStruct);
    SetConnectedComponentNumber( inp_cur_data->at, inp_cur_data->num_at, i+1 ); /*  normalization alters structure component number */
    for ( k = 0; k < TAUT_NUM; k ++ ) {
        if ( cur_INChI_Aux[k] && cur_INChI_Aux[k]->nNumberOfAtoms > 0 ) {
            pncFlags->bNormalizationFlags[iINChI][k] |= cur_INChI_Aux[k]->bNormalizationFlags;
            pncFlags->bTautFlags[iINChI][k]          |= cur_INChI_Aux[k]->bTautFlags;
            pncFlags->bTautFlagsDone[iINChI][k]      |= cur_INChI_Aux[k]->bTautFlagsDone;
            pncFlags->nCanonFlags[iINChI][k]         |= cur_INChI_Aux[k]->nCanonFlags;
        }
    }

    /*  Detect errors */
    if ( num_at < 0 ) {
        sd->nErrorCode = num_at;
    } else
    if ( num_at == 0 ) {
        sd->nErrorCode = -1;
    } else
    if ( cur_INChI[TAUT_NON] && cur_INChI[TAUT_NON]->nErrorCode ) {
        /*  non-tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_NON]->nErrorCode;
    } else
    if ( cur_INChI[TAUT_YES] && cur_INChI[TAUT_YES]->nErrorCode ) {
        /*  tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_YES]->nErrorCode;
    }
#if ( bRELEASE_VERSION == 0 )
    if ( cur_INChI[TAUT_NON] ) sd->bExtract |= cur_INChI[TAUT_NON]->bExtract;
    if ( cur_INChI[TAUT_YES] ) sd->bExtract |= cur_INChI[TAUT_YES]->bExtract;
    if ( (TG_FLAG_TEST_TAUT3_SALTS_DONE & bTautFlagsDone) ) {
        sd->bExtract |= EXTR_TEST_TAUT3_SALTS_DONE;
    }
#endif
    /*  detect and store stereo warnings */
    if ( !sd->nErrorCode ) {
        GetProcessingWarnings(cur_INChI, inp_norm_data, sd);
    }

    lElapsedTime = InchiTimeElapsed( &ulTStart );
    if ( ip->msec_MaxTime ) {
        ip->msec_LeftTime -= lElapsedTime;
    }
    sd->ulStructTime += lElapsedTime;
#ifndef TARGET_API_LIB
    /*  Display the results */
    if ( ip->bDisplay )
        eat_keyboard_input();
#endif
    /*  a) No matter what happened save the allocated INChI pointers */
    /*  save the INChI of the current component */

    InchiTimeGet( &ulTStart );
    for ( k = 0; k < TAUT_NUM; k ++ ) {
        pINChI[i][k]     = cur_INChI[k];
        pINChI_Aux[i][k] = cur_INChI_Aux[k];

        cur_INChI[k]     = NULL;
        cur_INChI_Aux[k] = NULL;
    }

    /*  b) Count one component structure and/or INChI results only if there was no error */
    /*     Set inp_norm_data[j]->num_removed_H = number of removed explicit H           */

    if ( !sd->nErrorCode ) {

        /*  find where the current processed structure is located */
        int cur_is_in_non_taut = (pINChI[i][TAUT_NON] && pINChI[i][TAUT_NON]->nNumberOfAtoms>0);
        int cur_is_in_taut     = (pINChI[i][TAUT_YES] && pINChI[i][TAUT_YES]->nNumberOfAtoms>0);
        int cur_is_non_taut = (cur_is_in_non_taut && 0 == pINChI[i][TAUT_NON]->lenTautomer) ||
                              (cur_is_in_taut     && 0 == pINChI[i][TAUT_YES]->lenTautomer);
        int cur_is_taut     = cur_is_in_taut     && 0 <  pINChI[i][TAUT_YES]->lenTautomer;
        /*
        sd->bTautFlags[iINChI]     |= bTautFlags;
        sd->bTautFlagsDone[iINChI] |= bTautFlagsDone;
        */
        if ( cur_is_non_taut + cur_is_taut ) {
            /*  count tautomeric and non-tautomeric components of the structures */
            int j1 = cur_is_in_non_taut? TAUT_NON:TAUT_YES;
            int j2 = cur_is_in_taut?     TAUT_YES:TAUT_NON;
            int j;
            sd->num_non_taut[iINChI] += cur_is_non_taut;
            sd->num_taut[iINChI]     += cur_is_taut;
            for ( j = j1; j <= j2; j ++ ) {
                int bIsotopic = (pINChI[i][j]->nNumberOfIsotopicAtoms ||
                                 pINChI[i][j]->nNumberOfIsotopicTGroups ||
                                 (pINChI[i][j]->nPossibleLocationsOfIsotopicH && pINChI[i][j]->nPossibleLocationsOfIsotopicH[0]>1));
                if ( j == TAUT_YES ) {
                    bIsotopic |= (0 < pINChI_Aux[i][j]->nNumRemovedIsotopicH[0] + 
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[1] +
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[2]);
                }
                inp_norm_data[j]->bExists = 1; /*  j=0: non-taut exists, j=1: taut exists */
                inp_norm_data[j]->bHasIsotopicLayer = bIsotopic;
                /*inp_norm_data[j]->num_removed_H = inp_norm_data[j]->num_at - num_at;*/
            }
        }
    }
/*
    return (sd->nErrorCode==CT_OUT_OF_RAM || sd->nErrorCode==CT_USER_QUIT_ERR)? _IS_FATAL :
            sd->nErrorCode? _IS_ERROR : 0;
*/
    if ( sd->nErrorCode==CT_OUT_OF_RAM || sd->nErrorCode==CT_USER_QUIT_ERR ) {
        ret = _IS_FATAL;
    } else
    if ( sd->nErrorCode ) {
        ret = _IS_ERROR;
    }
    lElapsedTime = InchiTimeElapsed( &ulTStart );
    if ( ip->msec_MaxTime ) {
        ip->msec_LeftTime -= lElapsedTime;
    }
    sd->ulStructTime += lElapsedTime;
    return ret;
}
/****************************************************************************************************/
int TreatCreateOneComponentINChIError(STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data,
                                     int i, long num_inp,
                                     INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file, 
                                     INCHI_IOSTREAM *prb_file, /*^^^ was: INCHI_IOSTREAM */
                                     char *pStr, int nStrLen )
{
    if ( sd->nErrorCode ) {
        AddMOLfileError(sd->pStrErrStruct, ErrMsg(sd->nErrorCode) );
        inchi_ios_eprint( log_file, "Error %d (%s) structure #%ld component %d.%s%s%s%s\n",
            sd->nErrorCode, sd->pStrErrStruct, num_inp, i+1, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
        sd->nErrorType = (sd->nErrorCode==CT_OUT_OF_RAM || sd->nErrorCode==CT_USER_QUIT_ERR)? _IS_FATAL : _IS_ERROR;
        if ( (ip->bINChIOutputOptions & INCHI_OUT_XML)
#ifdef TARGET_LIB_FOR_WINCHI
             || (ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW) && (ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT)
#endif
            ) {
            sd->nErrorType = ProcessStructError( output_file, log_file, /*sd->nErrorCode,*/ sd->pStrErrStruct,
                                            sd->nErrorType, &sd->bXmlStructStarted, num_inp, ip, pStr, nStrLen );
            /*  save the problem structure */
            if ( prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && !ip->bSaveAllGoodStructsAsProblem ) {
                CopyMOLfile(inp_file->f, sd->fPtrStart, sd->fPtrEnd, prb_file->f, num_inp);
            }
        } else {
            /*  save the problem structure */
            if ( sd->nErrorCode && prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && !ip->bSaveAllGoodStructsAsProblem ) {
                CopyMOLfile(inp_file->f, sd->fPtrStart, sd->fPtrEnd, prb_file->f, num_inp);
            }
        }
    }
/*^^^ #ifndef TARGET_API_LIB */
#if ( !defined( TARGET_API_LIB ) && !defined(TARGET_EXE_STANDALONE) )
    /*  print the logfile record */
    if ( log_file->f && log_file->f != stderr && (sd->ulStructTime >= 1000 || sd->nErrorCode) ) {
        fprintf( log_file->f, "%10lu msec structure #%ld.%s%s%s%s (%d component%s, %d atom%s, error=%d).\n",
                sd->ulStructTime, num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue),
                orig_inp_data->num_components, orig_inp_data->num_components==1?"":"s",
                orig_inp_data->num_inp_atoms, orig_inp_data->num_inp_atoms==1?"":"s", sd->nErrorCode );
    }
#endif
    return sd->nErrorType;
}
/****************************************************************************************************/
int TreatCreateINChIWarning(STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, long num_inp,
                                     INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file, 
                                     INCHI_IOSTREAM *prb_file, /*^^^ was: INCHI_IOSTREAM */
                                     char *pStr, int nStrLen )
{
#if ( bRELEASE_VERSION == 0 && (EXTR_FLAGS || EXTR_MASK) )
    if ( EXTR_MASK? ((sd->bExtract & EXTR_MASK) == EXTR_FLAGS) : (sd->bExtract & EXTR_FLAGS) ) {
        char szMsg[64];
        sprintf( szMsg, "ExtractStruct.code=0x%X", sd->bExtract);
        AddMOLfileError(sd->pStrErrStruct, szMsg);
    }
#endif
    if ( !sd->nErrorCode && sd->pStrErrStruct[0] ) {
        inchi_ios_eprint( log_file, "Warning (%s) structure #%ld.%s%s%s%s\n",
            sd->pStrErrStruct, num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
        sd->nErrorType = _IS_WARNING;
        if ( (ip->bINChIOutputOptions & INCHI_OUT_XML)
#ifdef TARGET_LIB_FOR_WINCHI
             || (ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW) && (ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT)
#endif
            ) {
            sd->nErrorType = ProcessStructError( output_file, log_file, /*sd->nErrorCode,*/ sd->pStrErrStruct,
                                            sd->nErrorType, &sd->bXmlStructStarted, num_inp, ip, pStr, nStrLen );
        }
        /*  save the structure as a problem structure if requested */
        if ( ip->bSaveWarningStructsAsProblem && !ip->bSaveAllGoodStructsAsProblem &&
             prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd ) {
            CopyMOLfile(inp_file->f, sd->fPtrStart, sd->fPtrEnd, prb_file->f, num_inp);
        }
#if ( bRELEASE_VERSION == 0 )
        /*  otherwise extract the structure as a problem structure if requested */
        else
        if ( (EXTR_MASK? ((sd->bExtract & EXTR_MASK) == EXTR_FLAGS) : (sd->bExtract & EXTR_FLAGS)) && !ip->bSaveAllGoodStructsAsProblem &&
             prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd ) {
            CopyMOLfile(inp_file->f, sd->fPtrStart, sd->fPtrEnd, prb_file->f, num_inp);
        }
#endif
    }
#if ( bRELEASE_VERSION != 1 && bOUTPUT_ONE_STRUCT_TIME == 1 )
#ifndef TARGET_API_LIB
    if ( log_file && log_file != stderr ) {
        fprintf( log_file, "%10lu msec structure %1dD #%ld.%s%s%s%s (%d component%s, %d atom%s, error=%d).\n",
                sd->ulStructTime, orig_inp_data->num_dimensions, num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue),
                orig_inp_data->num_components, orig_inp_data->num_components==1?"":"s",
                orig_inp_data->num_inp_atoms, orig_inp_data->num_inp_atoms==1?"":"s", sd->nErrorCode );
    }
#else
    if ( log_file ) {
        inchi_ios_eprint( log_file, "%10lu msec structure %1dD #%ld.%s%s%s%s (%d component%s, %d atom%s, error=%d).\n",
                sd->ulStructTime, orig_inp_data->num_dimensions, num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue),
                orig_inp_data->num_components, orig_inp_data->num_components==1?"":"s",
                orig_inp_data->num_inp_atoms, orig_inp_data->num_inp_atoms==1?"":"s", sd->nErrorCode );
    }
#endif
#endif
    return sd->nErrorType;
}
/*******************************************************************************************/
int DuplicateOrigAtom( ORIG_ATOM_DATA *new_orig_atom, ORIG_ATOM_DATA *orig_atom )
{
    inp_ATOM  *at                 = NULL;
    AT_NUMB   *nCurAtLen          = NULL;
    AT_NUMB   *nOldCompNumber     = NULL;

    if ( new_orig_atom->at && new_orig_atom->num_inp_atoms >= orig_atom->num_inp_atoms ) {
        at             = new_orig_atom->at;
    } else {
        at             = (inp_ATOM *)inchi_calloc(orig_atom->num_inp_atoms+1, sizeof(at[0]));
    }
    if ( new_orig_atom->nOldCompNumber && new_orig_atom->num_components >= orig_atom->num_components ) {
        nCurAtLen      = new_orig_atom->nCurAtLen;
    } else {
        nCurAtLen      = (AT_NUMB *)inchi_calloc(orig_atom->num_components+1, sizeof(nCurAtLen[0]));
    }
    if ( new_orig_atom->nCurAtLen && new_orig_atom->num_components >= orig_atom->num_components ) {
        nOldCompNumber = new_orig_atom->nOldCompNumber;
    } else {
        nOldCompNumber = (AT_NUMB *)inchi_calloc(orig_atom->num_components+1, sizeof(nOldCompNumber[0]));
    }

    if ( at && nCurAtLen && nOldCompNumber ) {
        /* copy */
        if ( orig_atom->at )
            memcpy( at, orig_atom->at, orig_atom->num_inp_atoms * sizeof(new_orig_atom->at[0]) );
        if ( orig_atom->nCurAtLen )
            memcpy( nCurAtLen, orig_atom->nCurAtLen, orig_atom->num_components*sizeof(nCurAtLen[0]) );
        if ( orig_atom->nOldCompNumber )
            memcpy( nOldCompNumber, orig_atom->nOldCompNumber, orig_atom->num_components*sizeof(nOldCompNumber[0]) );
        /* deallocate */
        if ( new_orig_atom->at && new_orig_atom->at != at )
            inchi_free( new_orig_atom->at );
        if ( new_orig_atom->nCurAtLen && new_orig_atom->nCurAtLen != nCurAtLen )
            inchi_free( new_orig_atom->nCurAtLen );
        if ( new_orig_atom->nOldCompNumber && new_orig_atom->nOldCompNumber != nOldCompNumber )
            inchi_free( new_orig_atom->nOldCompNumber );

        *new_orig_atom                = *orig_atom;
        new_orig_atom->at             = at;
        new_orig_atom->nCurAtLen      = nCurAtLen;
        new_orig_atom->nOldCompNumber = nOldCompNumber;
        /* data that are not to be copied */
        new_orig_atom->nNumEquSets    = 0;
        memset(new_orig_atom->bSavedInINCHI_LIB, 0, sizeof(new_orig_atom->bSavedInINCHI_LIB));
        memset(new_orig_atom->bPreprocessed,    0, sizeof(new_orig_atom->bPreprocessed));
        /* arrays that are not to be copied */
        new_orig_atom->szCoord        = NULL;
        new_orig_atom->nEquLabels     = NULL;
        new_orig_atom->nSortedOrder   = NULL;
        return 0;
    }

    /* deallocate */
    if ( at && new_orig_atom->at != at )
        inchi_free( at );
    if ( nCurAtLen && new_orig_atom->nCurAtLen != nCurAtLen )
        inchi_free( nCurAtLen );
    if ( nOldCompNumber && new_orig_atom->nOldCompNumber != nOldCompNumber )
        inchi_free( nOldCompNumber );

    return -1; /* failed */
}
#ifndef TARGET_API_LIB
/*******************************************************************************************/
int GetOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                     INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file, 
                     INCHI_IOSTREAM *prb_file, /*^^^ was: INCHI_IOSTREAM */
                     ORIG_ATOM_DATA *orig_inp_data, long *num_inp, char *pStr, int nStrLen, STRUCT_FPTRS *struct_fptrs )
{
    int nRet, inp_index, out_index, bUseFptr = (NULL != struct_fptrs);

    FreeOrigAtData( orig_inp_data );
    /*
    FreeOrigAtData( orig_inp_data + 1 );
    FreeOrigAtData( orig_inp_data + 2 );
    */

    /* added for TARGET_LIB_FOR_WINCHI early EOF detection */
    inp_index = -1;
    out_index = -1;
    if ( struct_fptrs ) {
        if ( inp_file->f == stdin ) {
            return _IS_FATAL;
        }
        if ( ip->nInputType == INPUT_CMLFILE ) {
            bUseFptr = 0;
        }
        /* initially allocate or increase length of struct_fptrs->fptr array */
        if ( !struct_fptrs->fptr || struct_fptrs->len_fptr <= struct_fptrs->cur_fptr+1 ) {
            INCHI_FPTR *new_fptr = (INCHI_FPTR *)inchi_calloc( struct_fptrs->len_fptr + ADD_LEN_STRUCT_FPTRS, sizeof(new_fptr[0]) );
            if ( new_fptr ) {
                if ( struct_fptrs->fptr ) {
                    if ( struct_fptrs->len_fptr ) {
                        memcpy( new_fptr, struct_fptrs->fptr, struct_fptrs->len_fptr*sizeof(new_fptr[0]));
                    }
                    inchi_free( struct_fptrs->fptr );
                } else {
                    struct_fptrs->len_fptr = 0;
                    struct_fptrs->cur_fptr = 0;
                    struct_fptrs->max_fptr = 0;
                }
                struct_fptrs->len_fptr += ADD_LEN_STRUCT_FPTRS;
                struct_fptrs->fptr = new_fptr;
            } else {
                return _IS_FATAL;  /* new_fptr allocation error */
            }
        }
        if ( struct_fptrs->fptr[struct_fptrs->cur_fptr] == EOF ) {
            return _IS_EOF;
        } else {
            if ( bUseFptr ) {
                if( fseek( inp_file->f, struct_fptrs->fptr[struct_fptrs->cur_fptr], SEEK_SET) ) {
                    return _IS_FATAL;
                }
                if ( struct_fptrs->cur_fptr && struct_fptrs->max_fptr <= struct_fptrs->cur_fptr ) {
                    return _IS_FATAL;
                }
            } else {
                inp_index = struct_fptrs->fptr[struct_fptrs->cur_fptr];
                out_index = EOF;
            }
        }
        *num_inp = struct_fptrs->cur_fptr; /* set structure count */
    }

    nRet = ReadTheStructure( sd, ip, inp_file, orig_inp_data, inp_index, &out_index );

    if ( !nRet ) {
        /*****************************************************
         * In case of no error output structure xml start tag
         * output read the structure errors and warnings
         *****************************************************/
        if ( ip->nInputType == INPUT_INCHI_PLAIN || ip->nInputType == INPUT_INCHI_XML ||
             ip->nInputType == INPUT_MOLFILE     || ip->nInputType == INPUT_SDFILE) {
            if ( ip->lMolfileNumber ) {
                *num_inp = ip->lMolfileNumber;
            } else {
                *num_inp += 1;
            }
        } else {
            *num_inp += 1;
        }
        nRet = TreatReadTheStructureErrors( sd, ip, LOG_MASK_ALL, inp_file, log_file, output_file, prb_file,
                                            orig_inp_data, num_inp, pStr, nStrLen );
    }

    /************************************************************/
    /* added for TARGET_LIB_FOR_WINCHI: look ahead for end of file detection */
    /************************************************************/
    if ( struct_fptrs && struct_fptrs->fptr && struct_fptrs->fptr[struct_fptrs->cur_fptr+1] <= 0 ) {
        int nRet2 = 0;
        INCHI_FPTR next_fptr;
        STRUCT_DATA sd2;

        if ( nRet != _IS_EOF && nRet != _IS_FATAL ) {
            if ( inp_file->f == stdin || struct_fptrs->len_fptr <= struct_fptrs->cur_fptr+1 ) {
                return _IS_FATAL;
            }
            /* get next structure fptr */
            if ( bUseFptr ) {
                next_fptr = ftell( inp_file->f );
            } else {
                inp_index = out_index;
                out_index = EOF;
            }
            /* read the next structure */
            nRet2 = ReadTheStructure( &sd2, ip, inp_file, NULL, inp_index, &out_index );
            /* restore fptr to the next structure */
            if ( bUseFptr ) {
                if ( next_fptr != -1L ) {
                    fseek( inp_file->f, next_fptr, SEEK_SET);
                }
            }
#if ( ADD_CMLPP == 1 )
            else {
                if ( inp_index >= 0 ) {
                    SetCmlStructIndex( inp_index ); /* so far nothing to do */
                }
            }
#endif
        } else {
            /* treat current fatal error as end of file */
            struct_fptrs->fptr[struct_fptrs->cur_fptr] = EOF;
        }
        /* next is end of file or fatal */
        if ( nRet  == _IS_EOF || nRet  == _IS_FATAL ||
             nRet2 == _IS_EOF || nRet2 == _IS_FATAL ) {
            struct_fptrs->fptr[struct_fptrs->cur_fptr+1] = EOF;
        } else {
            struct_fptrs->fptr[struct_fptrs->cur_fptr+1] = bUseFptr? sd->fPtrEnd : inp_index;
        }

        /* update struct_fptrs->max_fptr */
        if ( struct_fptrs->max_fptr <= struct_fptrs->cur_fptr+1  ) {
            struct_fptrs->max_fptr = struct_fptrs->cur_fptr+2;
        }
    }

    switch ( nRet ) {
    case _IS_EOF:
        *num_inp -= 1;
    case _IS_FATAL:
    case _IS_ERROR:
    case _IS_SKIP:
        goto exit_function;
    }

    /*
    if ( !orig_inp_data->num_dimensions ) {
        AddMOLfileError(sd->pStrErrStruct, "0D"); */ /* 0D-structure: no coordinates
    }
    */


exit_function:
    return nRet;
}
#endif
#if ( TEST_RENUMB_ATOMS == 1 ) /* { */
/************************************************************************************************/
int RenumberingTestInit( RENUMB_DATA *pRenumbData, INP_ATOM_DATA *inp_cur_data )
{
    int j;
    pRenumbData->ren_inp_norm_data[0] = &pRenumbData->ren_inp_norm_data1;
    pRenumbData->ren_inp_norm_data[1] = &pRenumbData->ren_inp_norm_data2;
    memset( pRenumbData->ren_INChI2, 0, sizeof( pRenumbData->ren_INChI2 ));
    memset( pRenumbData->ren_INChI_Aux, 0, sizeof( pRenumbData->ren_INChI_Aux ));
    memset( &pRenumbData->orig_inp_cur_data, 0, sizeof( pRenumbData->orig_inp_cur_data ));
    memset( &pRenumbData->saved_inp_cur_data, 0, sizeof( pRenumbData->saved_inp_cur_data ));
    memset( pRenumbData->ren_inp_norm_data[0], 0, sizeof( *pRenumbData->ren_inp_norm_data[0] ));
    memset( pRenumbData->ren_inp_norm_data[1], 0, sizeof( *pRenumbData->ren_inp_norm_data[1] ));
#if ( TEST_RENUMB_ATOMS_SAVE_LONGEST == 1 )
    memset( &pRenumbData->longest_inp_cur_data, 0, sizeof(pRenumbData->longest_inp_cur_data));
#endif
    CopyInpAtomData( &pRenumbData->orig_inp_cur_data, inp_cur_data );
    pRenumbData->ren_counter = pRenumbData->orig_inp_cur_data.num_at * pRenumbData->orig_inp_cur_data.num_at;
    srand(1);  /*  for reproducibility */
    rand();    /*  shift to avoid prev. sequences */
    pRenumbData->nComp = 0;
    /*ren_counter = 29;*/
    pRenumbData->new_ord = (AT_RANK *)inchi_calloc( pRenumbData->orig_inp_cur_data.num_at, sizeof(pRenumbData->new_ord[0]) );
    if ( pRenumbData->new_ord ) {
        for ( j = 0; j < pRenumbData->orig_inp_cur_data.num_at; j ++ ) {
            pRenumbData->new_ord[j] = (AT_RANK)j;
        }
        return 0;
    }
    return -1; /* out of RAM */
}
/************************************************************************************************/
int RenumberingTestUninit( RENUMB_DATA *pRenumbData )
{
    FreeInpAtomData( &pRenumbData->orig_inp_cur_data );
#if ( TEST_RENUMB_ATOMS_SAVE_LONGEST == 1 )
    FreeInpAtomData( &pRenumbData->longest_inp_cur_data );
#endif
    inchi_free( pRenumbData->new_ord );
    return 0;
}
/************************************************************************************************/
int RenumberingTest( PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, ORIG_ATOM_DATA *orig_inp_data, int iINChI,
                     RENUMB_DATA *pRenumbData, INP_ATOM_DATA *inp_cur_data, INP_ATOM_DATA **inp_norm_data,
                     STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *prb_file,
                     int i, long num_inp, NORM_CANON_FLAGS *pncFlags)
{
    int k, bLongerTime;
    CopyInpAtomData( &pRenumbData->saved_inp_cur_data, inp_cur_data );
    pRenumbData->nRet2 = 0;
    pRenumbData->num_taut0     = sd->num_taut[iINChI];
    pRenumbData->num_non_taut0 = sd->num_non_taut[iINChI];
    pRenumbData->ulMaxTime    = 0;
    while ( -- pRenumbData->ren_counter >= 0 && !pRenumbData->nRet2 ) {
        pRenumbData->nComp ++;
        MakeNewOrd( pRenumbData->orig_inp_cur_data.num_at, pRenumbData->new_ord );
        RenumbInpAtomData( inp_cur_data /* output*/, &pRenumbData->orig_inp_cur_data/* input*/, pRenumbData->new_ord/* input*/ );
#if ( TEST_RENUMB_ATOMS_SAVE_LONGEST == 1 )
        CopyInpAtomData( &pRenumbData->longest_inp_cur_data, inp_cur_data );
#endif
        if ( 470 == pRenumbData->nComp ) {
            int stop = 1; /* debug only */
        }

        pRenumbData->nRet2 = CreateOneComponentINChI( sd, ip, inp_cur_data, NULL /*orig_inp_data*/,
                                       pRenumbData->ren_INChI2, pRenumbData->ren_INChI_Aux, iINChI,
                                       0, num_inp, pRenumbData->ren_inp_norm_data, pncFlags, log_file );
        /*
        CreateOneComponentINChI( sd, ip, inp_cur_data, orig_inp_data, pINChI2[iINChI], pINChI_Aux2[iINChI], iINChI,
                                       i, num_inp, inp_norm_data, log_file );
        */
        if ( !pRenumbData->nRet2 ) {
            pRenumbData->c1 = CompareINChI( pINChI[i][TAUT_NON], pRenumbData->ren_INChI2[0][TAUT_NON],
                                           pINChI_Aux[i][TAUT_NON], pRenumbData->ren_INChI_Aux[0][TAUT_NON]);
            pRenumbData->c2 = CompareINChI( pINChI[i][TAUT_YES], pRenumbData->ren_INChI2[0][TAUT_YES],
                                           pINChI_Aux[i][TAUT_YES], pRenumbData->ren_INChI_Aux[0][TAUT_YES]);
            if ( pRenumbData->c1 || pRenumbData->c2 || pRenumbData->nRet2 ) {

                /****** the renumbering result is different ******/

                inchi_ios_eprint( log_file, "Compare (%d,%d) %d (err=%d) %s structure #%d component %d.%s%s%s%s\n",
                                    pRenumbData->c1, pRenumbData->c2,
                                    pRenumbData->nComp, pRenumbData->nRet2, INCHI_NAME,
                                    num_inp, i+1, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
                for ( k = 0; k < pRenumbData->orig_inp_cur_data.num_at; k ++ ) {
                    inchi_ios_eprint( log_file, " %d", (int)pRenumbData->new_ord[k] );
                }
                inchi_ios_eprint( log_file, "\n" );
                pRenumbData->ren_counter = 0; /* force exit */
                pRenumbData->bRenumbErr = 1000*pRenumbData->c2 + pRenumbData->c1;
#if ( TEST_RENUMB_SWITCH == 1 )
                CopyInpAtomData( &pRenumbData->longest_inp_cur_data, inp_cur_data );
                if ( pRenumbData->longest_inp_cur_data.at ) {
                    for ( k = 0; k < pRenumbData->longest_inp_cur_data.num_at; k ++ ) {
                        pRenumbData->longest_inp_cur_data.at[k].orig_at_number = k+1; /* display new atom numbers */
                    }
                }
#endif
            }
#if ( TEST_RENUMB_ATOMS_SAVE_LONGEST == 1 )
            /*  output time per this component */
            inchi_ios_eprint( stderr, "\rComp#%d str#%ld/%d%s%s%s%s Ren %d/%d n(%lu:%lu)c(%lu:%lu)...\r",
                                i+1, num_inp, iINChI, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), pRenumbData->nComp, pRenumbData->ren_counter+pRenumbData->nComp,
                                pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulNormTime, pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulCanonTime,
                                pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulNormTime, pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulCanonTime);
#endif
            /*  make sure the max. time is not overwritten */
            pRenumbData->ulCurTime0 = pRenumbData->ren_INChI_Aux[0][TAUT_NON]?
                                      (pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulNormTime
                                     + pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulCanonTime) : 0;
            pRenumbData->ulCurTime1 = pRenumbData->ren_INChI_Aux[0][TAUT_YES]?
                                      (pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulNormTime
                                      + pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulCanonTime) : 0;
            pRenumbData->ulCurTime  = inchi_max( pRenumbData->ulCurTime0, pRenumbData->ulCurTime1 );

            pRenumbData->ulCurTimeCanon0 = pRenumbData->ren_INChI_Aux[0][TAUT_NON]? pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulCanonTime : 0;
            pRenumbData->ulCurTimeCanon1 = pRenumbData->ren_INChI_Aux[0][TAUT_YES]? pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulCanonTime : 0;
            pRenumbData->ulCurTimeCanon  = inchi_max( pRenumbData->ulCurTimeCanon0, pRenumbData->ulCurTimeCanon1);

            pRenumbData->ulCurTimeNorm0 = pRenumbData->ren_INChI_Aux[0][TAUT_NON]? pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulNormTime:0;
            pRenumbData->ulCurTimeNorm1 = pRenumbData->ren_INChI_Aux[0][TAUT_YES]? pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulNormTime:0;
            pRenumbData->ulCurTimeNorm  = inchi_max( pRenumbData->ulCurTimeNorm0, pRenumbData->ulCurTimeNorm1);


            bLongerTime = 0;
            if ( pRenumbData->ulCurTime > pRenumbData->ulMaxTime ) {
                pRenumbData->ulMaxTime = pRenumbData->ulCurTime;
                bLongerTime = 1;
            }
            if ( pRenumbData->ulMaxTimeCanon > pRenumbData->ulCurTimeCanon ) {
                pRenumbData->ulMaxTimeCanon = pRenumbData->ulCurTimeCanon;
                bLongerTime = 1;
            }
            if ( pRenumbData->ulMaxTimeNorm > pRenumbData->ulCurTimeCanon ) {
                pRenumbData->ulMaxTimeCanon = pRenumbData->ulCurTimeCanon;
                bLongerTime = 1;
            }
#if ( TEST_RENUMB_ATOMS_SAVE_LONGEST == 1 || TEST_RENUMB_SWITCH == 1 )
            if ( bLongerTime || TEST_RENUMB_SWITCH == 1 && (pRenumbData->c1 || pRenumbData->c2 || pRenumbData->nRet2) ) {
                char szLine[512];
                char szValue[512];
                inchi_ios_eprint( stderr, "\n" );
                sprintf( szLine, "Comp#%d str#%ld/%d%s%s%s%s Ren %d/%d n=%lu:%lu c=%lu:%lu",
                                i+1, num_inp, iINChI, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), pRenumbData->nComp, pRenumbData->ren_counter+pRenumbData->nComp,
                                pRenumbData->ren_INChI_Aux[0][TAUT_NON]? pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulNormTime:0,
                                pRenumbData->ren_INChI_Aux[0][TAUT_NON]? pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulCanonTime:0,
                                pRenumbData->ren_INChI_Aux[0][TAUT_YES]? pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulNormTime:0,
                                pRenumbData->ren_INChI_Aux[0][TAUT_YES]? pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulCanonTime:0);
                sprintf( szValue, "%s (c%d/s%ld/i%d, r%d/%d n=%lu:%lu c=%lu:%lu)",
                                  (ip->pSdfValue && ip->pSdfValue[0])? ip->pSdfValue:"unk",
                                i+1, num_inp, iINChI, pRenumbData->nComp, pRenumbData->ren_counter+pRenumbData->nComp,
                                pRenumbData->ren_INChI_Aux[0][TAUT_NON]? pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulNormTime:0,
                                pRenumbData->ren_INChI_Aux[0][TAUT_NON]? pRenumbData->ren_INChI_Aux[0][TAUT_NON]->ulCanonTime:0,
                                pRenumbData->ren_INChI_Aux[0][TAUT_YES]? pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulNormTime:0,
                                pRenumbData->ren_INChI_Aux[0][TAUT_YES]? pRenumbData->ren_INChI_Aux[0][TAUT_YES]->ulCanonTime:0);

                WriteToSDfile( &pRenumbData->longest_inp_cur_data, prb_file, szLine, NULL, ip->pSdfLabel, szValue );
            }
#endif

#if ( TEST_RENUMB_SWITCH == 1 )
            if ( pRenumbData->c1 || pRenumbData->c2 || !pRenumbData->ren_counter ) {
                inchi_swap( (char*)&pINChI[i][TAUT_NON], (char*)&pRenumbData->ren_INChI2[0][TAUT_NON], sizeof(&pRenumbData->ren_INChI2[0][0]) );
                inchi_swap( (char*)&pINChI[i][TAUT_YES], (char*)&pRenumbData->ren_INChI2[0][TAUT_YES], sizeof(&pRenumbData->ren_INChI2[0][0]) );
                inchi_swap( (char*)&pINChI_Aux[i][TAUT_NON], (char*)&pRenumbData->ren_INChI_Aux[0][TAUT_NON], sizeof(&pRenumbData->ren_INChI_Aux[0][0]) );
                inchi_swap( (char*)&pINChI_Aux[i][TAUT_YES], (char*)&pRenumbData->ren_INChI_Aux[0][TAUT_YES], sizeof(&pRenumbData->ren_INChI_Aux[0][0]) );
            }
#endif
        }

        for ( k = 0; k < TAUT_NUM; k ++ ) {
            if ( pRenumbData->ren_INChI2[0][k] ) {
                Free_INChI(&pRenumbData->ren_INChI2[0][k]);
                /*
                inchi_free(pRenumbData->ren_INChI2[0][k]);
                pRenumbData->ren_INChI2[0][k] = NULL;
                */
            }
            if ( pRenumbData->ren_INChI_Aux[0][k] ) {
                Free_INChI_Aux(&pRenumbData->ren_INChI_Aux[0][k]);
                /*
                inchi_free(pRenumbData->ren_INChI_Aux[0][k]);
                pRenumbData->ren_INChI_Aux[0][k] = NULL;
                */
            }
        }
    }
    /*  eliminate overcounting due to multiple renumberings/recalculations */
    pRenumbData->num_taut         = sd->num_taut[iINChI] - pRenumbData->num_taut0;
    pRenumbData->num_non_taut     = sd->num_non_taut[iINChI] - pRenumbData->num_non_taut0;
    sd->num_taut[iINChI]     = pRenumbData->num_taut0;
    sd->num_non_taut[iINChI] = pRenumbData->num_non_taut0;
    if ( pRenumbData->num_taut % pRenumbData->nComp || pRenumbData->num_non_taut % pRenumbData->nComp ) {
        inchi_ios_eprint( log_file, "Compare (%d,%d) %d (err=%d) %s structure #%ld component %d.%s%s%s%s\n",
                    pRenumbData->num_non_taut % pRenumbData->nComp, pRenumbData->num_taut % pRenumbData->nComp,
                    pRenumbData->nComp, 333, INCHI_NAME, num_inp, i+1, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
    }
#if ( TEST_RENUMB_SWITCH == 1 )  /* { */
    CopyInpAtomData( inp_norm_data[TAUT_NON], pRenumbData->ren_inp_norm_data[TAUT_NON] );
    CopyInpAtomData( inp_norm_data[TAUT_YES], pRenumbData->ren_inp_norm_data[TAUT_YES] );
    /*  renumbered input structure */
#ifndef COMPILE_ANSI_ONLY /* { */
    if ( /*ip->bDisplayEachComponentINChI &&*/ !pRenumbData->nRet2 ) {
        int err, len;
        /*
        err = DisplayStructure( inp_cur_data->at, inp_cur_data->num_at, 0, 1, 0, NULL. 1, 0, NULL, NULL,
                                            ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
        */
        err = DisplayStructure( inp_cur_data->at, inp_cur_data->num_at, 0, 1, 0, NULL, 1/*isotopic*/, 0/*taut*/, NULL, NULL,
                                            ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
        if ( pRenumbData->c1 || pRenumbData->c2 ) {
            len = strlen(szTitle);
            strcat( szTitle, " (Renumbered)" );
            err = DisplayStructure( pRenumbData->longest_inp_cur_data.at, pRenumbData->longest_inp_cur_data.num_at,
                                    0, 1, 0, NULL, 1, 0, NULL, NULL, ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
            szTitle[len] = '\0';
        }
        sd->bUserQuitComponentDisplay = (err==ESC_KEY);
        if ( !err ) {
            inchi_ios_eprint( stderr, "Cannot display the structure\n");
        }
    }
#endif  /* } COMPILE_ANSI_ONLY */
#else /* } TEST_RENUMB_SWITCH { */
    CopyInpAtomData( inp_cur_data, &pRenumbData->saved_inp_cur_data );
#endif /* } TEST_RENUMB_SWITCH */
    FreeInpAtomData( &pRenumbData->saved_inp_cur_data );
    FreeInpAtomData( pRenumbData->ren_inp_norm_data[TAUT_NON] );
    FreeInpAtomData( pRenumbData->ren_inp_norm_data[TAUT_YES] );
#if ( TEST_RENUMB_ATOMS_SAVE_LONGEST == 1 || TEST_RENUMB_SWITCH == 1 )
    FreeInpAtomData( &pRenumbData->longest_inp_cur_data );
#endif
    return pRenumbData->nRet2;
}
#endif  /* } TEST_RENUMB_ATOMS  */

/****************************************************************************/
int bCheckUnusualValences( ORIG_ATOM_DATA *orig_at_data, int bAddIsoH, char *pStrErrStruct )
{
    int i, val, num_found = 0;
    char msg[32];
    int len, num_H;
    inp_ATOM *at = ( orig_at_data && orig_at_data->num_inp_atoms > 0 )? orig_at_data->at : NULL;

    if ( at ) {
        for ( i = 0, num_found = 0; i < orig_at_data->num_inp_atoms; i ++ ) {
            num_H = bAddIsoH? NUMH(at,i) : at[i].num_H;
            val = detect_unusual_el_valence( at[i].el_number, at[i].charge, at[i].radical,
                                          at[i].chem_bonds_valence, num_H, at[i].valence );
            if ( val ) {
                num_found ++;
                /* produce message */
                AddMOLfileError(pStrErrStruct, "Accepted unusual valence(s):");
                len = sprintf( msg, "%s", at[i].elname );
                if ( at[i].charge ) {
                    len += sprintf( msg+len, "%+d", at[i].charge );
                }
                if ( at[i].radical ) {
                    len += sprintf( msg + len, ",%s", at[i].radical == RADICAL_SINGLET? "s" :
                                                      at[i].radical == RADICAL_DOUBLET? "d" :
                                                      at[i].radical == RADICAL_TRIPLET? "t" : "?" );
                }
                len += sprintf( msg + len, "(%d)", val ); 
                AddMOLfileError(pStrErrStruct, msg);
            }
        }
    }
    return num_found;
}
/***************************************************************************/
int PreprocessOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data )
{        
        int i;
    INCHI_MODE bTautFlags     = 0;
    INCHI_MODE bTautFlagsDone = 0;
    /*************************************************/
    /* 1. copy orig_inp_data --> prep_inp_data */
    /*************************************************/
    if ( 0 > DuplicateOrigAtom( prep_inp_data, orig_inp_data ) ) {
        AddMOLfileError(sd->pStrErrStruct, "Out of RAM");
        sd->nStructReadError =  99;
        sd->nErrorType = _IS_FATAL;
        goto exit_function;
    }
#if ( bRELEASE_VERSION == 0 && (EXTR_HAS_METAL_ATOM & (EXTR_MASK | EXTR_FLAG) ) )
        if ( bHasMetalAtom( orig_inp_data ) ) {
            sd->bExtract |= EXTR_HAS_METAL_ATOM;
        }
#endif

    /*************************************************/
    /* 2. fix odd things in prep_inp_data            */
    /*************************************************/

    if ( 0 < fix_odd_things( prep_inp_data->num_inp_atoms, prep_inp_data->at, /*0*/ip->bTautFlags & TG_FLAG_FIX_SP3_BUG, ip->bFixNonUniformDraw ) ) { /* changed 2010-03-17 DT */
        AddMOLfileError(sd->pStrErrStruct, "Charges were rearranged");
        if ( sd->nErrorType < _IS_WARNING ) {
            sd->nErrorType = _IS_WARNING;
        }
        sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
    }
#if ( FIX_ADJ_RAD == 1 )
    if ( ip->bTautFlags & TG_FLAG_FIX_ADJ_RADICALS ) {
        if ( 0 < FixAdjacentRadicals( prep_inp_data->num_inp_atoms, prep_inp_data->at ) ) {
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ADJ_RADICALS_DONE;
        }
    }
#endif
#if ( bRELEASE_VERSION == 0 && (EXTR_FLAGS & EXTR_HAS_FEATURE) )
    if ( bFoundFeature( prep_inp_data->at, prep_inp_data->num_inp_atoms ) ) {
        sd->bExtract |= EXTR_HAS_FEATURE;
    }
#endif

    /*******************************************************************
     * Find whether the structure can be disconnected or is a salt
     *******************************************************************/


    /* needs salt disconnection? */
    if ( ip->bTautFlags & TG_FLAG_DISCONNECT_SALTS ) {
        prep_inp_data->bDisconnectSalts = (0 < DisconnectSalts( prep_inp_data, 0 ));
    } else {
        prep_inp_data->bDisconnectSalts = 0;
    }
    /* needs metal disconnection? */
    if ( ip->bTautFlags & TG_FLAG_DISCONNECT_COORD ) {
        i = (0 != (ip->bTautFlags & TG_FLAG_CHECK_VALENCE_COORD));
        bMayDisconnectMetals( prep_inp_data, i, &bTautFlagsDone ); /* changes prep_inp_data->bDisconnectCoord */
        sd->bTautFlagsDone[INCHI_BAS] |= bTautFlagsDone; /* whether any disconnection has been rejected because of the metal proper valence */
#if ( bRELEASE_VERSION == 0 )
        if ( i && (bTautFlagsDone & TG_FLAG_CHECK_VALENCE_COORD_DONE) ) {
            sd->bExtract |= EXTR_METAL_WAS_NOT_DISCONNECTED;
        }
#endif
    } else {
        prep_inp_data->bDisconnectCoord = 0;
    }
    orig_inp_data->bDisconnectSalts = prep_inp_data->bDisconnectSalts;
    orig_inp_data->bDisconnectCoord = prep_inp_data->bDisconnectCoord;

    /*************************************************/
    /* 3. if( orig_inp_data->bDisconnectSalts ) then */
    /*       -- disconnect salts in prep_inp_data */
    /*************************************************/

    if ( ( ip->bTautFlags & TG_FLAG_DISCONNECT_SALTS ) && prep_inp_data->bDisconnectSalts &&
         0 < (i=DisconnectSalts( prep_inp_data, 1 )) ) {
        AddMOLfileError(sd->pStrErrStruct, "Salt was disconnected");
        sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_DISCONNECT_SALTS_DONE;
        if ( sd->nErrorType < _IS_WARNING ) {
            sd->nErrorType = _IS_WARNING;
        }
        if ( (i = ReconcileAllCmlBondParities( prep_inp_data->at, prep_inp_data->num_inp_atoms, 0 )) ) {
            char szErrCode[16];
            sprintf( szErrCode, "%d", i);
            AddMOLfileError( sd->pStrErrStruct, "0D Parities Reconciliation failed:" );
            AddMOLfileError( sd->pStrErrStruct, szErrCode );
        }
#if ( bRELEASE_VERSION == 0 )
        sd->bExtract |= EXTR_SALT_WAS_DISCONNECTED;
#endif
    } else {
        prep_inp_data->bDisconnectSalts = 0;
    }

    /***********************************************************/
    /*  mark the (disconnected) components in prep_inp_data    */
    /***********************************************************/

    prep_inp_data->num_components = MarkDisconnectedComponents( prep_inp_data, 0 );

    if ( prep_inp_data->num_components < 0 ) {
        AddMOLfileError(sd->pStrErrStruct, "Out of RAM");
        sd->nStructReadError =  99;
        sd->nErrorType = _IS_FATAL;
        goto exit_function;
    }

    /***********************************************************/
    /* Detect isotopic H on heteroatoms -- necessary condition */
    /* for global isotopic tautomerism                         */
    /***********************************************************/

    if ( (i = bNumHeterAtomHasIsotopicH( prep_inp_data->at, prep_inp_data->num_inp_atoms )) ) {
        if ( i & 1 ) {
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FOUND_ISOTOPIC_H_DONE;
        }
        if ( i & 2 ) {
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE;
        }
    }


    /****************************************************************************/
    /* 4a. Detect unusual valences                                              */
    /*     should be called before metal disconnection                          */
    /****************************************************************************/


    if ( bCheckUnusualValences( prep_inp_data, 1, sd->pStrErrStruct ) ) {
#if ( bRELEASE_VERSION == 0 )
        sd->bExtract |= EXTR_UNUSUAL_VALENCES;
#else
        ;
#endif
    }
    /***********************************************************/
    /*    5. if( orig_inp_data->bDisconnectCoord ) then        */
    /*        -- copy prep_inp_data --> prep_inp_data+1        */
    /*        -- disconnect metals in prep_inp_data            */
    /***********************************************************/

    if ( prep_inp_data->bDisconnectCoord ) {

        prep_inp_data->num_components = MarkDisconnectedComponents( prep_inp_data, 0 );
        if ( prep_inp_data->num_components < 0 ) {
            AddMOLfileError(sd->pStrErrStruct, "Out of RAM");
            sd->nStructReadError =  99;
            sd->nErrorType = _IS_FATAL;
            goto exit_function;
        }
        /* save Reconnected structure in prep_inp_data+1 if requested */
        if ( 0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD) ) {
            if ( 0 > DuplicateOrigAtom( prep_inp_data+1, prep_inp_data ) ) {
                AddMOLfileError(sd->pStrErrStruct, "Out of RAM");
                sd->nStructReadError =  99;
                sd->nErrorType = _IS_FATAL;
                goto exit_function;
            }
            sd->bTautFlags[INCHI_REC]     = sd->bTautFlags[INCHI_BAS];
            sd->bTautFlagsDone[INCHI_REC] = sd->bTautFlagsDone[INCHI_BAS];
            {   /* remove "parity undefined in disconnected structure" flag from reconnected structure */
                int k, m, p;
                inp_ATOM *at     = (prep_inp_data+1)->at;
                int       num_at = (prep_inp_data+1)->num_inp_atoms;
                for ( k = 0; k < num_at; k ++ ) {
                    for ( m = 0; m < MAX_NUM_STEREO_BONDS && (p=at[k].sb_parity[m]); m ++ ) {
                        at[k].sb_parity[m] &= SB_PARITY_MASK;
                    }
                }
            }
        }

        /* make Disconnected structure in prep_inp_data */
        i = (0 != ( ip->bTautFlags & TG_FLAG_CHECK_VALENCE_COORD ));
        /* prep_inp_data->bDisconnectCoord > 1 means add prep_inp_data->bDisconnectCoord-1 explicit H atoms */
        if ( 0 < (i = DisconnectMetals( prep_inp_data, i, &bTautFlagsDone ) ) ) {
            AddMOLfileError(sd->pStrErrStruct, "Metal was disconnected");
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_DISCONNECT_COORD_DONE;
            if ( sd->nErrorType < _IS_WARNING ) {
                sd->nErrorType = _IS_WARNING;
            }
#if ( bRELEASE_VERSION == 0 )
            sd->bExtract |= EXTR_METAL_WAS_DISCONNECTED;
#endif
            /* last parm=1 means find link between unchanged by Metal Disconnection components */
            prep_inp_data->num_components = MarkDisconnectedComponents( prep_inp_data, 1 );
            if ( prep_inp_data->num_components < 0 ) {
                AddMOLfileError(sd->pStrErrStruct, "Out of RAM");
                sd->nStructReadError =  99;
                sd->nErrorType = _IS_FATAL;
                goto exit_function;
            }

            {   /* set parities for the disconnected structure */
                int k, m, p;
                inp_ATOM *at     = (prep_inp_data)->at;
                int       num_at = (prep_inp_data)->num_inp_atoms;
                for ( k = 0; k < num_at; k ++ ) {
                    for ( m = 0; m < MAX_NUM_STEREO_BONDS && (p=at[k].sb_parity[m]); m ++ ) {
                        if ( p & SB_PARITY_FLAG ) {
                            at[k].sb_parity[m] = (p >> SB_PARITY_SHFT) & SB_PARITY_MASK;
                        }
                    }
                }
            }

            if ( (i = ReconcileAllCmlBondParities( prep_inp_data->at, prep_inp_data->num_inp_atoms, 1 )) ) {
                char szErrCode[16];
                sprintf( szErrCode, "%d", i);
                AddMOLfileError( sd->pStrErrStruct, "0D Parities Reconciliation failed:" );
                AddMOLfileError( sd->pStrErrStruct, szErrCode );
            }

#if ( REMOVE_ION_PAIRS_DISC_STRU == 1 )
            if ( 0 < remove_ion_pairs( prep_inp_data->num_inp_atoms, prep_inp_data->at ) ) {
                AddMOLfileError(sd->pStrErrStruct, "Charges were rearranged");
                if ( sd->nErrorType < _IS_WARNING ) {
                    sd->nErrorType = _IS_WARNING;
                }
                sd->bTautFlagsDone[INCHI_REC] |= TG_FLAG_FIX_ODD_THINGS_DONE;
                sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
            }
#endif

            /*
              if prep_inp_data->nOldCompNumber[i] = iINChI+1 > 0 then
              component #(i+1) in prep_inp_data is identical to component #(iINChI+1) in prep_inp_data+1
            */
        } else
        if ( i < 0 ) {
            AddMOLfileError(sd->pStrErrStruct, "Cannot disconnect metal error");
            sd->nStructReadError =  i;
            sd->nErrorType = _IS_ERROR;
            goto exit_function;
        }
    } else 
    {   /* remove "disconnected structure parities" from the structure */
        int k, m, p;
        inp_ATOM *at     = (prep_inp_data)->at;
        int       num_at = (prep_inp_data)->num_inp_atoms;
        for ( k = 0; k < num_at; k ++ ) {
            for ( m = 0; m < MAX_NUM_STEREO_BONDS && (p=at[k].sb_parity[m]); m ++ ) {
                at[k].sb_parity[m] &= SB_PARITY_MASK;
            }
        }
    }


exit_function:

    if ( sd->nErrorType < _IS_ERROR && prep_inp_data ) {

        if ( 0 < post_fix_odd_things( prep_inp_data->num_inp_atoms, prep_inp_data->at ) ) {
            AddMOLfileError(sd->pStrErrStruct, "Charges were rearranged");
            if ( sd->nErrorType < _IS_WARNING ) {
                sd->nErrorType = _IS_WARNING;
            }
            sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
        }
        if ( (sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE) &&
             (prep_inp_data+1)->at && (prep_inp_data+1)->num_inp_atoms > 0 ) {
            if ( 0 < post_fix_odd_things( (prep_inp_data+1)->num_inp_atoms, (prep_inp_data+1)->at ) ) {
                AddMOLfileError(sd->pStrErrStruct, "Charges were rearranged");
                if ( sd->nErrorType < _IS_WARNING ) {
                    sd->nErrorType = _IS_WARNING;
                }
                sd->bTautFlagsDone[INCHI_REC] |= TG_FLAG_FIX_ODD_THINGS_DONE;
                sd->bTautFlagsDone[INCHI_BAS] |= TG_FLAG_FIX_ODD_THINGS_DONE;
            }
        }
    }
    
    sd->bTautFlags[INCHI_BAS]     |= bTautFlags;  /* TG_FLAG_CHECK_VALENCE_COORD_DONE, TG_FLAG_MOVE_CHARGE_COORD_DONE */
    sd->bTautFlagsDone[INCHI_BAS] |= bTautFlagsDone;  /* TG_FLAG_CHECK_VALENCE_COORD_DONE, TG_FLAG_MOVE_CHARGE_COORD_DONE */
    return sd->nErrorType;
}

#ifndef COMPILE_ANSI_ONLY  /* { */
/************************************************************************************************/
int DisplayTheWholeStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, 
                              INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file,
                              ORIG_ATOM_DATA *orig_inp_data, long num_inp, int iINChI, int bShowStruct, int bINCHI_LIB_Flag )
{

    int bDisplayEqu = 0;
#ifndef TARGET_LIB_FOR_WINCHI
    /* Displaying equivalent input structures when disconnection has been done: */
    /* in case of TARGET_LIB_FOR_WINCHI equivalence info is always unknown here and bOriginalReconnected=0 */
    int bOriginalReconnected = iINChI < 0 && orig_inp_data && orig_inp_data->nEquLabels &&
                               (sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE) &&
                               (ip->bTautFlags & TG_FLAG_RECONNECT_COORD);
    const char *lpszType = bOriginalReconnected? " (Reconnected)"   :
                           (iINChI <  0        )? ""                 :
                           (iINChI == INCHI_BAS )? " (Preprocessed)"  :
                           (iINChI == INCHI_REC )? " (Reconnected)"   : "";
    int err = 0;
    /* Display the original structure */
    bDisplayEqu = bShowStruct && ip->bDisplay &&
                  ip->dp.nEquLabels && 0 < ip->dp.nCurEquLabel && ip->dp.nCurEquLabel <= ip->dp.nNumEquSets;
#else
    if(!DRAWDATA || !DRAWDATA_EXISTS)
        return 0;
#endif
#ifndef TARGET_API_LIB
    /********************************************************************
     * Ask the user whether to process the input structure or quit
     */
    if ( ip->bDisplay && inp_file->f != stdin ) {
        if ( user_quit(bDisplayEqu?"Enter=Display identical components, Esc=Stop ?" : "Enter=Display, Esc=Stop ?", ip->ulDisplTime) ) {
            sd->bUserQuit = 1;
            goto exit_function;
        }
    }
#endif
    /******************************************************
     * Display the whole input structure in console app
     */
/*^^^ #ifndef TARGET_LIB_FOR_WINCHI */
#if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined(TARGET_EXE_STANDALONE) )
    if ( bShowStruct && ip->bDisplay ) {
        if ( bDisplayEqu ) {
            sprintf( szTitle, " Equ Set %d of %d, Input Structure #%ld.%s%s%s%s%s",
                     ip->dp.nCurEquLabel, ip->dp.nNumEquSets,
                     num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), lpszType);
        } else {
            sprintf( szTitle, "Input Structure #%ld.%s%s%s%s%s", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), lpszType);
        }
        err = DisplayStructure( orig_inp_data->at, orig_inp_data->num_inp_atoms, 0, 1, 0, NULL, 1/*isotopic*/, 0/*taut*/, NULL, NULL,
                                            ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
        sd->bUserQuitComponent = (err==ESC_KEY);
        if ( !err ) {
            inchi_fprintf( stderr, "Cannot display the structure\n");
        }
    }
    if( !bDisplayEqu ) {
        /*  console output progress report */
        if ( ip->bDisplay && !sd->bUserQuitComponent ) {
            if ( iINChI == 1 ) {
                if ( ip->bDisplay )
                    inchi_ios_eprint( log_file, "Processing (rec) structure #%ld.%s%s%s%s...\n", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
                else
                    inchi_fprintf( stderr, "Processing (rec) structure #%ld.%s%s%s%s...\r", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            } else {
                if ( ip->bDisplay )
                    inchi_ios_eprint( log_file, "Processing structure #%ld.%s%s%s%s...\n", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
                else
                    inchi_fprintf( stderr, "Processing structure #%ld.%s%s%s%s...\r", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            }
        }
    }
#endif


    /******************************************************
     * Store the whole input structure in GUI application
     */
#ifdef TARGET_LIB_FOR_WINCHI
    if ( ip->bDisplay && bINCHI_LIB_Flag )
#else
    if ( (ip->bDisplay || (ip->bCompareComponents & CMP_COMPONENTS)) && bINCHI_LIB_Flag )
#endif
    {
        int bBit, k, bReconnected, nComponent, bPreprocessed;
        for ( bBit = 1, k = 0; k < 8; k ++, bBit <<= 1 ) {
            /******************************************************************************
             *  bReconnected  = k%2     (0 or 1)
             *  nComponent    = k/4     (0 or 1)
             *  bPreprocessed = (k/2)%2 (0 or 1)
             ******************************************************************************/
            if ( !(bINCHI_LIB_Flag & bBit) ) {
                continue;
            }
            bReconnected  = k%2;
            nComponent    = k/4;
            bPreprocessed = ((k/2)%2);

            sprintf( szTitle, "%s Structure #%ld.%s%s%s%s",
                              bPreprocessed? "Preprocessed" : bReconnected? "Reconnected" : "Input",
                              num_inp,
                              SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue));

#ifdef TARGET_LIB_FOR_WINCHI
            if(DRAWDATA && DRAWDATA_EXISTS)
            {
                struct DrawData vDrawData;
                int    nType = bPreprocessed? COMPONENT_ORIGINAL_PREPROCESSED : COMPONENT_ORIGINAL;
                if ( DRAWDATA_EXISTS( nComponent, bPreprocessed, bReconnected ) ) {
                    sd->nErrorType = _IS_FATAL;
                    sd->nErrorCode = CT_UNKNOWN_ERR;
                    return -1;
                }
                vDrawData.pWindowData = CreateWinData_( orig_inp_data->at, orig_inp_data->num_inp_atoms,
                                                        0, 1 /* bAdd_DT_to_num_H */, 0, NULL, 1, 0, NULL, NULL,
                                                        ip->bAbcNumbers, &ip->dp, ip->nMode );
                if( vDrawData.pWindowData != NULL )
                {
                    vDrawData.nComponent   = nComponent;
                    vDrawData.nType        = nType; /* COMPONENT_ORIGINAL or COMPONENT_ORIGINAL_PREPROCESSED */
                    vDrawData.bReconnected = bReconnected; /* 0=>main; 1=>reconnected */
                    vDrawData.pWindowData->szTitle = _strdup(szTitle);
                    vDrawData.szTitle              = _strdup(szTitle);
                    DRAWDATA(&vDrawData);
                    if ( !nComponent ) {
                        /* keep track of saved INCHI_LIB data */
                        orig_inp_data->bSavedInINCHI_LIB[bReconnected] ++;
                        orig_inp_data->bPreprocessed[bReconnected]    = bPreprocessed;
                    }
                }
            }
#else
            if ( !nComponent ) {
                /* keep track of saved INCHI_LIB data */
                orig_inp_data->bSavedInINCHI_LIB[bReconnected] ++;
                orig_inp_data->bPreprocessed[bReconnected]    = bPreprocessed;
            }
#endif

        }
    }

exit_function:
    return sd->bUserQuit;
}
#endif /* } COMPILE_ANSI_ONLY */
/************************************************************************************************/
int ProcessOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                         PINChI2 *pINChI[INCHI_NUM], PINChI_Aux2 *pINChI_Aux[INCHI_NUM],
                         INCHI_IOSTREAM *inp_file, 
                         INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file, INCHI_IOSTREAM *prb_file, /*^^^ was: INCHI_IOSTREAM */
                         ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                         long num_inp, char *pStr, int nStrLen,
                         unsigned char save_opt_bits)
{
        int nRet = 0, nRet1, i, k, maxINChI=0;
        COMP_ATOM_DATA composite_norm_data[INCHI_NUM][TAUT_NUM+1]; /* [0]:non-taut, [1]:taut, [2]:intermediate taut struct */
        NORM_CANON_FLAGS ncFlags;
        NORM_CANON_FLAGS *pncFlags = &ncFlags;
        ORIG_STRUCT      OrigStruct;
        ORIG_STRUCT      *pOrigStruct = NULL;
        int bSortPrintINChIFlags=0;


#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
        int ret1=0, ret2=0;
#endif
        sd->bUserQuitComponent = 0;
        sd->bUserQuitComponentDisplay = 0;
        memset( composite_norm_data, 0, sizeof(composite_norm_data) );
        memset( pncFlags, 0, sizeof(*pncFlags) );
        /* ip->msec_LeftTime = ip->msec_MaxTime; */ /* start timeout countdown */

        /* for testing only */
#if ( REMOVE_ION_PAIRS_ORIG_STRU == 1 )
        fix_odd_things( orig_inp_data->num_inp_atoms, orig_inp_data->at, 0, ip->bFixNonUniformDraw );
#endif

#if ( UNDERIVATIZE == 1 )  /***** post v.1 feature *****/
        if ( ip->bUnderivatize && 0 > (ret2=underivatize( orig_inp_data )) ) {
            long num_inp2 = num_inp;
            AddMOLfileError(sd->pStrErrStruct, "Underivatization error");
            sd->nStructReadError =  99;
            sd->nErrorType = _IS_ERROR;
            nRet = _IS_ERROR;
            TreatReadTheStructureErrors( sd, ip, LOG_MASK_ALL, inp_file, log_file, output_file, prb_file,
                                        prep_inp_data, &num_inp2, pStr, nStrLen );
            goto exit_function; /* output only if derivatives found */
        }
#endif /* UNDERIVATIZE == 1 */
#if ( RING2CHAIN == 1 )  /***** post v.1 feature *****/
        if ( ip->bRing2Chain && 0 > (ret1 = Ring2Chain( orig_inp_data )) ) {
            long num_inp2 = num_inp;
            AddMOLfileError(sd->pStrErrStruct, "Ring to chain error");
            sd->nStructReadError =  99;
            sd->nErrorType = _IS_ERROR;
            nRet = _IS_ERROR;
            TreatReadTheStructureErrors( sd, ip, LOG_MASK_ALL, inp_file, log_file, output_file, prb_file,
                                        prep_inp_data, &num_inp2, pStr, nStrLen );
            goto exit_function; /* output only if derivatives found */
        }
#endif /* RING2CHAIN == 1 */
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )  /***** post v.1 feature *****/
        if ( ip->bIngnoreUnchanged && !ret1 && !ret2 ) {
            goto exit_function; /* output only if derivatives or ring/chain found */
        }
#endif /* RING2CHAIN == 1 || UNDERIVATIZE == 1 */


        /***** output MOLfile ***************/
        if ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY  ) {
            char szNumber[32];
            int ret1a=0, ret2a=0; /* for derivatives and ring-chain */
/*^^^ #if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined( TARGET_API_LIB ) ) */
#if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined( TARGET_API_LIB ) && !defined(TARGET_EXE_STANDALONE) )
#if ( TEST_RENUMB_ATOMS != 1 )
            /*  log file / console output */
            if ( log_file->f != stderr ) {
                if ( ip->bDisplay )
                    inchi_ios_eprint( log_file, "Writing structure #%ld.%s%s%s%s...\n", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
                else
                    inchi_fprintf( stderr, "Writing structure #%ld.%s%s%s%s...\r", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            }
#endif
#endif
            ret1a = sprintf( szNumber, "Structure #%ld", num_inp );
            ret2a = WriteOrigAtomDataToSDfile( orig_inp_data, output_file, szNumber, NULL,
                (sd->bChiralFlag & FLAG_INP_AT_CHIRAL)? 1:0,
                (ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ATOMS_DT)? 1:0, ip->pSdfLabel, ip->pSdfValue );
            goto exit_function;
        }

        /******* create full reversibility information **************/
        if ( !(ip->bINChIOutputOptions & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO)) ) {
            pOrigStruct = &OrigStruct;
            memset( pOrigStruct, 0, sizeof(*pOrigStruct));
            if ( FillOutOrigStruct( orig_inp_data, pOrigStruct, sd ) ) {
                AddMOLfileError(sd->pStrErrStruct, "Cannot interpret reversibility information");
                sd->nStructReadError =  99;
                sd->nErrorType = _IS_ERROR;
                nRet = _IS_ERROR;
            }
        }
        /* create INChI for each connected component of the structure and optionally display them */
        /* create INChI for the whole disconnected or original structure */
        if ( nRet != _IS_FATAL && nRet != _IS_ERROR ) {
            nRet1 = CreateOneStructureINChI(sd, ip, szTitle, pINChI, pINChI_Aux, INCHI_BAS,
                                            inp_file, log_file, output_file, prb_file,
                                            orig_inp_data, prep_inp_data, 
                                            composite_norm_data, num_inp, 
                                            pStr, nStrLen, pncFlags );
            nRet = inchi_max(nRet, nRet1);
        }
        if ( nRet != _IS_FATAL && nRet != _IS_ERROR ) {
            maxINChI = 1;
        }


        if ( nRet != _IS_FATAL && nRet != _IS_ERROR &&
             (sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE) &&
             (ip->bTautFlags               & TG_FLAG_RECONNECT_COORD)          ) {

            /* create INChI for the whole reconnected structure */
            nRet1 = CreateOneStructureINChI(sd, ip, szTitle, pINChI, pINChI_Aux, INCHI_REC,
                                            inp_file, log_file, output_file, prb_file,
                                            orig_inp_data, prep_inp_data, 
                                            composite_norm_data,num_inp, 
                                            pStr, nStrLen, pncFlags);
            nRet = inchi_max(nRet, nRet1);
            if ( nRet != _IS_FATAL && nRet != _IS_ERROR ) {
                maxINChI = 2;
            }
        }


        if ( nRet != _IS_FATAL && nRet != _IS_ERROR ) {

            if ( (sd->bChiralFlag & FLAG_INP_AT_CHIRAL) &&
                  (ip->nMode & REQ_MODE_STEREO) &&
                 !(ip->nMode & (REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO)) &&
                 !bIsStructChiral( pINChI, sd->num_components ) ) {

                AddMOLfileError(sd->pStrErrStruct, "Not chiral");
            }
            /*************************************/
            /*       Output err/warn messages    */
            /*************************************/
            if ( /*!sd->nErrorCode &&*/ !sd->bUserQuitComponent && !sd->bUserQuit ) {
                /*  if successful then returns 0, otherwise returns _IS_FATAL */
                /*  exctract the structure if requested */
                nRet1 = TreatCreateINChIWarning(sd, ip, prep_inp_data, num_inp,
                                 inp_file, log_file, output_file, prb_file,pStr, nStrLen );
                nRet = inchi_max(nRet, nRet1);
            }
        }


            /************************************************/
            /*  sort and print INChI for the whole structure */
            /************************************************/

        if ( ip->nInputType != INPUT_INCHI ) 
        {
            /* Prepare SaveOpt bits */
            save_opt_bits = 0;
            if ( ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT )
            {
                if ( 0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD) )
                    save_opt_bits |= SAVE_OPT_RECMET;
                if ( 0 != ( ip->nMode & REQ_MODE_BASIC) )
                    save_opt_bits |= SAVE_OPT_FIXEDH;
                if ( 0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO) ) 
                    save_opt_bits |= SAVE_OPT_SLUUD;
                if ( 0 == (ip->nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU)) ) 
                    save_opt_bits |= SAVE_OPT_SUU;
                if ( 0 != (ip->bTautFlags & TG_FLAG_KETO_ENOL_TAUT) )
                    save_opt_bits |= SAVE_OPT_KET;
                if ( 0 != (ip->bTautFlags & TG_FLAG_1_5_TAUT) )
                    save_opt_bits |= SAVE_OPT_15T;
                /* Check if /SNon requested and turn OFF stereo bits if so */
                if ( ! (ip->nMode & REQ_MODE_STEREO) )
                {
                    save_opt_bits &= ~SAVE_OPT_SUU;
                    save_opt_bits &= ~SAVE_OPT_SLUUD;
                }
            }
        }

        if ( nRet != _IS_FATAL && nRet != _IS_ERROR ) {
             
            nRet1 = SortAndPrintINChI(output_file, pStr, nStrLen, log_file, 
                                      ip, orig_inp_data, prep_inp_data,
                                      composite_norm_data, pOrigStruct,
                                      sd->num_components, sd->num_non_taut, sd->num_taut,
                                      sd->bTautFlags, sd->bTautFlagsDone, pncFlags, num_inp,
                                      pINChI, pINChI_Aux, 
                                      &bSortPrintINChIFlags, save_opt_bits);
            nRet = inchi_max(nRet, nRet1);
        }
#ifndef COMPILE_ANSI_ONLY /* { */
        
        /* display equivalent components on original or preprocessed structure(s) */

#ifndef TARGET_LIB_FOR_WINCHI

        if ( nRet != _IS_FATAL && nRet != _IS_ERROR && /*ip->bDisplay &&*/
             (ip->bCompareComponents & CMP_COMPONENTS) && !sd->bUserQuit && !sd->bUserQuitComponent ) 
        {
            int j, ret, ord;
            int bDisplaySaved = ip->bDisplay;
            ORIG_ATOM_DATA *inp_data;
            AT_NUMB         nEquSet;
            for ( ord = -1; ord < INCHI_NUM; ord ++ ) {
                switch( ord ) {
                case -1:
                    j = INCHI_BAS;  /* preprocessed non-tautomeric */
                    break;
                case 0:
                    j = INCHI_REC;  /* preprocessed tautomeric */
                    break;
                case 1:
                    j = -1;        /* original input */
                    break;
                default:
                    continue;
                }
                inp_data   = j < 0? orig_inp_data : prep_inp_data+j;
                if ( inp_data && inp_data->num_inp_atoms && inp_data->at &&
                     inp_data->nEquLabels &&
                     inp_data->nNumEquSets ) {
                    for ( nEquSet = 1; nEquSet <= inp_data->nNumEquSets; nEquSet ++ ) {
                        ip->dp.nEquLabels   = inp_data->nEquLabels;
                        ip->dp.nCurEquLabel = nEquSet;
                        ip->dp.nNumEquSets  = inp_data->nNumEquSets;
                        ip->bDisplay = 1; /* force display if it was not requested */
                        ret = DisplayTheWholeStructure( sd, ip, szTitle, inp_file, log_file, inp_data, num_inp,
                                                       j, 1 /*bShowStructure*/, 0 );
                        ip->dp.nEquLabels   = NULL;
                        ip->dp.nCurEquLabel = 0;
                        ip->dp.nNumEquSets  = 0;
                        ip->bDisplay = bDisplaySaved; /* restore display option */
                        if ( ret ) {
                            /* user pressed Esc */
                            goto exit_loop;
                        }
                    }
                }
            }
exit_loop:;
        }

#endif



        /* display composite results and equivalent components on composite results */
        if ( nRet != _IS_FATAL && nRet != _IS_ERROR && /*ip->bDisplay &&*/
             ip->bDisplayCompositeResults ) {
            int iINChI;
            for ( iINChI = 0; iINChI < maxINChI && !sd->bUserQuitComponentDisplay; iINChI ++ ) {
                DisplayTheWholeCompositeStructure( ip, sd, num_inp, iINChI,
                                               pINChI[iINChI], pINChI_Aux[iINChI],
                                               orig_inp_data, prep_inp_data,
                                               composite_norm_data[iINChI] );
            }
#ifndef TARGET_LIB_FOR_WINCHI
            if( !ip->bDisplay && sd->bUserQuitComponentDisplay ) {
                sd->bUserQuit = 1;
            }
#endif
        }

#endif /* } COMPILE_ANSI_ONLY */


        /* XML struct end tag */
        if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && sd->bXmlStructStarted > 0 ) {
            if ( !OutputINChIXmlStructEndTag( output_file, pStr, nStrLen, 1 ) ) {
                inchi_ios_eprint( log_file, "Cannot create end xml tag for structure #%ld.%s%s%s%s Terminating.\n", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
                sd->bXmlStructStarted = -1; /*  do not repeat same message */
                nRet = _IS_FATAL;
            } else {
                sd->bXmlStructStarted = 0; /*  do not continue xml output for this structure */
            }
        }
        if ( nRet != _IS_FATAL && nRet != _IS_ERROR ) {
            /* Special mode: extract all good MOLfiles into the problem file
             * Do not extract any MOLfile that could not be processed (option /PGO)
             */
            if ( prb_file && prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && ip->bSaveAllGoodStructsAsProblem ) {
                CopyMOLfile(inp_file->f, sd->fPtrStart, sd->fPtrEnd, prb_file->f, 0);
            }
#if ( /*bRELEASE_VERSION != 1 &&*/ EXTR_FLAGS == EXTR_TRANSPOSITION_EXAMPLES && EXTR_MASK == EXTR_FLAGS )
            else
            if ( prb_file->f && (bSortPrintINChIFlags & 
                   ( FLAG_SORT_PRINT_TRANSPOS_BAS | FLAG_SORT_PRINT_TRANSPOS_REC ) )
            ) {
                CopyMOLfile(inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, 0);
            }
#endif
        }
        for ( i = 0; i < INCHI_NUM; i ++ ) {
            for ( k = 0; k < TAUT_NUM+1; k ++ ) {
                FreeCompAtomData( &composite_norm_data[i][k] );
            }
        }
        FreeOrigStruct( pOrigStruct);


/*
        FreeInpAtomData( inp_cur_data     );
        FreeInpAtomData( inp_norm_data[0] );
        FreeInpAtomData( inp_norm_data[1] );
*/
exit_function:


        return nRet;
}
/************************************************************************************************/
int bIsStructChiral( PINChI2 *pINChI2[INCHI_NUM], int num_components[] )
{
    int i, j, k;
    INChI *pINChI;
    INChI_Stereo *Stereo;
    for ( j = 0; j < INCHI_NUM; j ++ ) {  /* disconnected / reconnected */
        if ( !num_components[j] ) {
            continue;
        }
        for ( i = 0; i < num_components[j]; i ++ ) {  /* component */
            for ( k = 0; k < TAUT_NUM; k ++ ) {       /* mobile/immobile H */
                if ( (pINChI = pINChI2[j][i][k]) &&
                      !pINChI->bDeleted &&
                      pINChI->nNumberOfAtoms > 0 ) {
                    
                    if ( (Stereo = pINChI->Stereo) &&
                         Stereo->t_parity &&
                         Stereo->nNumberOfStereoCenters > 0 &&
                         Stereo->nCompInv2Abs ) {
                        return 1; /* inversion changed stereo */
                    }
                    if ( (Stereo = pINChI->StereoIsotopic) &&
                         Stereo->t_parity &&
                         Stereo->nNumberOfStereoCenters > 0 &&
                         Stereo->nCompInv2Abs ) {
                        return 1; /* inversion changed stereo */
                    }
                }
            }
        }
    }
    return 0;
}
/************************************************************************************************/
int CreateOneStructureINChI( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                         PINChI2 *pINChI2[INCHI_NUM], PINChI_Aux2 *pINChI_Aux2[INCHI_NUM], int iINChI,
                         INCHI_IOSTREAM *inp_file, 
                         INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file, INCHI_IOSTREAM *prb_file, /*^^^ was: INCHI_IOSTREAM */
                         ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                         COMP_ATOM_DATA composite_norm_data2[][TAUT_NUM+1],
                         long num_inp, char *pStr, int nStrLen, NORM_CANON_FLAGS *pncFlags )
{
    int i, j, k, /*m,*/ nRet = 0;
#ifndef TARGET_LIB_FOR_WINCHI
    int n;
#ifndef COMPILE_ANSI_ONLY
    int err;
#endif
#endif

    PINChI2     *pINChI     = NULL;
    PINChI_Aux2 *pINChI_Aux = NULL;

    INP_ATOM_DATA InpCurAtData;
    INP_ATOM_DATA *inp_cur_data;

    INP_ATOM_DATA InpNormAtData, InpNormTautData;
    INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
    ORIG_ATOM_DATA *cur_prep_inp_data = prep_inp_data + iINChI;
    inchiTime      ulTStart;
#ifndef COMPILE_ANSI_ONLY
    int            bShowStructure = 0;
    int            bStructurePreprocessed = 0; /* All changes except disconnection */
    int            bStructureDisconnected = 0;
    int            bAlsoOutputReconnected = 0, bINCHI_LIB_Flag = 0;
    COMP_ATOM_DATA *composite_norm_data = composite_norm_data2[iINChI];
    INP_ATOM_DATA2 *all_inp_norm_data = NULL;
#endif

    /*
        if ( orig_inp_data is NOT empty AND
             prep_inp_data[0] IS empty ) then:

            1. copy orig_inp_data --> prep_inp_data[0]
            2. fix odd things in prep_inp_data[0]
            3. if( orig_inp_data->bDisconnectSalts ) then
                  -- disconnect salts in prep_inp_data[0]
            4. move protons to neutralize charges on heteroatoms
            5. if( orig_inp_data->bDisconnectCoord ) then
                  -- copy prep_inp_data[0] --> prep_inp_data[1]
                  -- disconnect metals in prep_inp_data[0]

            [ This all is done in PreprocessOneStructure() ]

        iINChI = 0
        =========
        (normal/disconnected layer)

            1. normalize prep_inp_data[0] in inp_norm_data[0,1]
            2. create INChI[ iINChI ] out of inp_norm_data[0,1]


        iINChI = 1 AND orig_inp_data->bDisconnectCoord > 0
        =================================================
        (reconnected layer)

            1. normalize prep_inp_data[1] in inp_norm_data[0,1]
            2. create INChI[ iINChI ] out of inp_norm_data[0,1]

    */

#if ( TEST_RENUMB_ATOMS == 1 )
    RENUMB_DATA RenumbData;
    RENUMB_DATA *pRenumbData = &RenumbData;
#endif


    ip->msec_LeftTime = ip->msec_MaxTime; /* start timeout countdown for each component */

#if ( TEST_RENUMB_ATOMS == 1 )
    memset( pRenumbData, 0, sizeof(*pRenumbData) );
#endif

    inp_cur_data     = &InpCurAtData;
    inp_norm_data[TAUT_NON] = &InpNormAtData;
    inp_norm_data[TAUT_YES] = &InpNormTautData;

    memset( inp_cur_data      , 0, sizeof( *inp_cur_data     ) );
    memset( inp_norm_data[TAUT_NON], 0, sizeof( *inp_norm_data[0] ) );
    memset( inp_norm_data[TAUT_YES], 0, sizeof( *inp_norm_data[0] ) );

#ifndef COMPILE_ANSI_ONLY
    memset( composite_norm_data+TAUT_NON, 0, sizeof( composite_norm_data[0] ) );
    memset( composite_norm_data+TAUT_YES, 0, sizeof( composite_norm_data[0] ) );
    memset( composite_norm_data+TAUT_INI, 0, sizeof( composite_norm_data[0] ) );
#endif
    if ( ip->bAllowEmptyStructure && !orig_inp_data->at && !orig_inp_data->num_inp_atoms ) {
        ;
    } else
    if ( !orig_inp_data->at || !orig_inp_data->num_inp_atoms ) {
        return 0; /* nothing to do */
    }
    if ( iINChI == 1 && orig_inp_data->bDisconnectCoord <= 0 ) {
        return 0;
    }

   /* m = iINChI; */ /* orig_inp_data index */

    if ( iINChI != INCHI_BAS && iINChI != INCHI_REC ) {
        AddMOLfileError(sd->pStrErrStruct, "Fatal undetermined program error");
        sd->nStructReadError =  97;
        nRet = sd->nErrorType = _IS_FATAL;
        goto exit_function;
    }

    /*******************************************************************
     *                                                                 *
     *                                                                 *
     *  Whole structure preprocessing: 1st step of the normalization   *
     *                                                                 *
     *  Happen only on the first call to CreateOneStructureINChI()      *
     *                                                                 *
     *                                                                 *
     *******************************************************************/

    if ( (!prep_inp_data->at || !prep_inp_data->num_inp_atoms) && orig_inp_data->num_inp_atoms > 0 ) {
        /* the structure has not been preprocessed */
        if ( ip->msec_MaxTime ) {
            InchiTimeGet( &ulTStart );
        }
        PreprocessOneStructure( sd, ip, orig_inp_data, prep_inp_data );
        pncFlags->bTautFlags[iINChI][TAUT_YES] =
        pncFlags->bTautFlags[iINChI][TAUT_NON] = sd->bTautFlags[INCHI_BAS] | ip->bTautFlags;
        pncFlags->bTautFlagsDone[iINChI][TAUT_YES] =
        pncFlags->bTautFlagsDone[iINChI][TAUT_NON] = sd->bTautFlagsDone[INCHI_BAS] | ip->bTautFlagsDone;

#ifndef COMPILE_ANSI_ONLY
        /* in this location the call happens once for each input structure, before preprocessing */
        bStructurePreprocessed = (0 != (sd->bTautFlagsDone[INCHI_BAS] & (
                                        TG_FLAG_MOVE_HPLUS2NEUTR_DONE  |
                                        TG_FLAG_DISCONNECT_SALTS_DONE  |
                                        TG_FLAG_MOVE_POS_CHARGES_DONE  |
                                        TG_FLAG_FIX_ODD_THINGS_DONE    )));
        bStructureDisconnected = (0 != (sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE));

        bShowStructure = ( bStructurePreprocessed ||
                           bStructureDisconnected ||
                           prep_inp_data[0].num_components > 1);

        /* sd->bTautFlags[] contains output flags
           ip->bTautFlags   contains input flags
        */
        bAlsoOutputReconnected = (sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE) &&
                                 (ip->bTautFlags               & TG_FLAG_RECONNECT_COORD);
        bINCHI_LIB_Flag = 0;

        /*************** output structures to TARGET_LIB_FOR_WINCHI conditions *********************
         *
         *  Send to TARGET_LIB_FOR_WINCHI:
         *
         *  type                      component  conditions
         *
         *  COMPONENT_ORIGINAL              #0:  (num_components > 1)
         *  COMPONENT_ORIGINAL_PREPROCESSED #0:  (num_components > 1) && (preprocessed)
         *  COMPONENT_ORIGINAL              #1:  (num_components = 1) && (preprocessed)
         *
         *  Flags explanation:
         *        MAIN => iINChI=0,  RECN => iINChI=1 (Reconnected)
         *        ORIG => Original, PREP => Preprocessed
         *
         *  Possible flags:           k
         *
         *  COMP_ORIG_0_MAIN  0x0001  0  COMPONENT_ORIGINAL, bMain, component #0
         *  COMP_ORIG_0_RECN  0x0002  1  COMPONENT_ORIGINAL, bRecn, component #0
         *
         *  COMP_PREP_0_MAIN  0x0004  2  COMPONENT_ORIGINAL_PREPROCESSED, bMain, component #0
         *  COMP_PREP_0_RECN  0x0008  3  COMPONENT_ORIGINAL_PREPROCESSED, bRecn, component #0
         *
         *  COMP_ORIG_1_MAIN  0x0010  4  COMPONENT_ORIGINAL, bMain, component #1
         *  COMP_ORIG_1_RECN  0x0020  5  COMPONENT_ORIGINAL, bRecn, component #1
         *
         *  bReconnected  = k%2     (0 or 1)
         *  nComponent    = k/4     (0 or 1)
         *  bPreprocessed = (k/2)%2 (0 or 1)
         *
         ******************************************************************************/
        /* Original -> Main, component #0, Original */
        if ( prep_inp_data[INCHI_BAS].num_components > 1 ) {
            bINCHI_LIB_Flag |= COMP_ORIG_0_MAIN;
        } else
        /* Original -> Main, component #1, Original */
        if ( prep_inp_data[INCHI_BAS].num_components == 1 && bStructurePreprocessed ) {
            bINCHI_LIB_Flag |= COMP_ORIG_1_MAIN;
            /* preprocessed will be added when output canonicalization results */
        }
        if ( bAlsoOutputReconnected ) {
            /* Original -> Reconnected, component #0, Original */
            if ( prep_inp_data[INCHI_REC].num_components > 1 ) {
                bINCHI_LIB_Flag |= COMP_ORIG_0_RECN;
            } else
            /* Original -> Reconnected, component #1, Original */
            if ( prep_inp_data[INCHI_BAS].num_components == 1 && bStructurePreprocessed ) {
                bINCHI_LIB_Flag |= COMP_ORIG_1_RECN;
                /* preprocessed will be added when output canonicalization results */
            }
        }
        if ( ip->msec_MaxTime ) {
            ip->msec_LeftTime -= InchiTimeElapsed( &ulTStart );
        }

        /* display the ORIGINAL, UN-PREPROCESSED structure */
        if ( DisplayTheWholeStructure( sd, ip, szTitle, inp_file, log_file, orig_inp_data, num_inp,
                                       -1, bShowStructure, bINCHI_LIB_Flag ) ) {
            goto exit_function;
        }
#endif
        switch (sd->nErrorType) {
        case _IS_ERROR:
        case _IS_FATAL:
            /* error message */
            nRet = TreatReadTheStructureErrors( sd, ip, LOG_MASK_ALL, inp_file, log_file, output_file, prb_file,
                                                prep_inp_data, &num_inp, pStr, nStrLen );
            goto exit_cycle;
        }
    }
    /* tranfer flags from INChI_Aux to sd */





#ifndef COMPILE_ANSI_ONLY /* { */

    /******************************************/
    /*      Displaying the structures         */
    /*          Only under WIN32              */
    /******************************************/
    if ( ip->bDisplayCompositeResults &&
        !sd->bUserQuitComponentDisplay && prep_inp_data[iINChI].num_components > 1) {
        all_inp_norm_data = (INP_ATOM_DATA2 *)inchi_calloc( prep_inp_data[iINChI].num_components, sizeof(all_inp_norm_data[0]));
    }


    /* Display the input structure AFTER PREPROCESSING */
    switch ( iINChI ) {

    case INCHI_BAS:
        /*------------ Possibly disconnected structure -------------------*/
        bStructurePreprocessed = 0 != (sd->bTautFlagsDone[iINChI] & (
                                        TG_FLAG_MOVE_HPLUS2NEUTR_DONE  |
                                        TG_FLAG_DISCONNECT_SALTS_DONE  |
                                        TG_FLAG_MOVE_POS_CHARGES_DONE  |
                                        TG_FLAG_MOVE_CHARGE_COORD_DONE |
                                        TG_FLAG_DISCONNECT_COORD_DONE  |
                                        TG_FLAG_FIX_ODD_THINGS_DONE    ));
        bINCHI_LIB_Flag = 0;

        /* Preprocessed/Main -> Main, component #0, Preprocessed */
        if ( prep_inp_data[iINChI].num_components > 1 && bStructurePreprocessed ) {
            bINCHI_LIB_Flag |= COMP_PREP_0_MAIN;
        }

        bShowStructure = ( bStructurePreprocessed && prep_inp_data[iINChI].num_components > 1);
        break;

    case INCHI_REC:
        /*------------ Reconnected structure ------------------------------*/
        bAlsoOutputReconnected = (sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE) &&
                                 (ip->bTautFlags               & TG_FLAG_RECONNECT_COORD);
        if ( !bAlsoOutputReconnected ) {
            break;
        }
        bStructurePreprocessed = 0 != (sd->bTautFlagsDone[iINChI] & (
                                        TG_FLAG_MOVE_HPLUS2NEUTR_DONE  |
                                        TG_FLAG_DISCONNECT_SALTS_DONE  |
                                        TG_FLAG_MOVE_POS_CHARGES_DONE  |
                                        TG_FLAG_FIX_ODD_THINGS_DONE    ));
        bINCHI_LIB_Flag = 0;

        /* Preprocessed/Reconnected -> Reconnected, component #0, Preprocessed */
        if ( prep_inp_data[iINChI].num_components > 1 && bStructurePreprocessed ) {
            bINCHI_LIB_Flag |= COMP_PREP_0_RECN;
        }

        bShowStructure = ( bStructurePreprocessed && prep_inp_data[iINChI].num_components > 1 );
        break;
    default:
        bShowStructure = 0;
    }
    if ( prep_inp_data[iINChI].num_inp_atoms > 0 ) {
        if ( DisplayTheWholeStructure( sd, ip, szTitle, inp_file, log_file, prep_inp_data+iINChI, num_inp,
                                       iINChI, bShowStructure, bINCHI_LIB_Flag ) ) {
            goto exit_function;
        }
    }
#endif /* } ifndef COMPILE_ANSI_ONLY */



    /* allocate pINChI[iINChI] and pINChI_Aux2[iINChI] -- arrays of pointers to INChI and INChI_Aux */
    /* assign values to sd->num_components[]                                                  */
    MYREALLOC2(PINChI2, PINChI_Aux2, pINChI2[iINChI], pINChI_Aux2[iINChI], sd->num_components[iINChI], cur_prep_inp_data->num_components, k);
    if ( k ) {
        AddMOLfileError(sd->pStrErrStruct, "Cannot allocate output data. Terminating");
        sd->nStructReadError =  99;
        sd->nErrorType = _IS_FATAL;
        goto exit_function;
    }
    pINChI     = pINChI2[iINChI];
    pINChI_Aux = pINChI_Aux2[iINChI];

    /**************************************************************************/
    /*                                                                        */
    /*                                                                        */
    /*   M A I N   C Y C L E:   P R O C E S S    C O M P O N E N T S          */
    /*                                                                        */
    /*                                                                        */
    /*                     O N E   B Y   O N E                                */
    /*                                                                        */
    /*                                                                        */
    /**************************************************************************/

    for ( i = 0, nRet = 0; !sd->bUserQuitComponent && i < cur_prep_inp_data->num_components; i ++ ) {
        if ( ip->msec_MaxTime ) {
            InchiTimeGet( &ulTStart );
        }
#ifndef TARGET_LIB_FOR_WINCHI  /* { */
#if ( bREUSE_INCHI == 1 )
        if ( iINChI == INCHI_REC && ((!ip->bDisplay && !ip->bDisplayCompositeResults && !(ip->bCompareComponents & CMP_COMPONENTS)) ||
                                   sd->bUserQuitComponentDisplay) ) {
            /* reconnected structure (06-20-2005: added "&& !ip->bDisplayCompositeResults" to display composite structure) */
            int m = iINChI-1;
            /* find whether we have already calculated this INChI in basic (disconnected) layer */
            for ( j = n = 0; j < prep_inp_data[m].num_components; j ++ ) {
                if ( i+1 == prep_inp_data[m].nOldCompNumber[j] &&
                     (pINChI2[m][j][TAUT_NON] || pINChI2[m][j][TAUT_YES]) ) {
                    /* yes, we have already done this */
                    if ( !n++ ) {
                        memcpy( pINChI    +i, pINChI2    [m]+j, sizeof(pINChI[0]));
                        memcpy( pINChI_Aux+i, pINChI_Aux2[m]+j, sizeof(pINChI_Aux[0]));
                        for ( k = 0; k < TAUT_NUM; k ++ ) {
                            if ( pINChI[i][k] ) {
                                pINChI[i][k]->nRefCount ++;
                                if ( pINChI[i][k]->nNumberOfAtoms > 0 ) {
                                    switch( k ) {
                                    case TAUT_NON:
                                        sd->num_non_taut[iINChI] ++;
                                        break;
                                    case TAUT_YES:
                                        if ( pINChI[i][k]->lenTautomer > 0 ) {
                                            sd->num_taut[iINChI] ++;
                                        } else
                                        if ( !pINChI[i][TAUT_NON] || !pINChI[i][TAUT_NON]->nNumberOfAtoms ) {
                                            sd->num_non_taut[iINChI] ++;
                                        }
                                        break;
                                    }
                                }
                            }
                            if ( pINChI_Aux[i][k] ) {
                                pINChI_Aux[i][k]->nRefCount ++;
                            }
                        }
                    }
                }
            }
            if ( n == 1 ) {
                continue;
            }
            if ( n > 1 ) { /* ith component is equivalent to more than one another component */
                AddMOLfileError(sd->pStrErrStruct, "Cannot distinguish components");
                sd->nStructReadError =  99;
                sd->nErrorType = _IS_ERROR;
                goto exit_function;
            }
        }
#endif
#endif /* } TARGET_LIB_FOR_WINCHI */

        /*****************************************************/
        /*  a) allocate memory and extract current component */
        /*****************************************************/
        nRet = GetOneComponent( sd, ip, log_file, output_file, inp_cur_data, cur_prep_inp_data, i, num_inp, pStr, nStrLen );
        if ( ip->msec_MaxTime ) {
            ip->msec_LeftTime -= InchiTimeElapsed( &ulTStart );
        }
        switch ( nRet ) {
        case _IS_ERROR:
        case _IS_FATAL:
            goto exit_cycle;
        }
#ifndef TARGET_API_LIB
        /*  console request: Display the component? */
        if ( ip->bDisplay && inp_file->f != stdin ) {
            if ( user_quit("Enter=Display Component, Esc=Stop ?", ip->ulDisplTime) ) {
                sd->bUserQuitComponent = 1;
                break;
            }
        }
#endif
#ifndef COMPILE_ANSI_ONLY  /* { */
        /*  b) Display the extracted original component structure */
        if ( inp_cur_data->at && ip->bDisplay && !sd->bUserQuitComponentDisplay ) {
            if ( cur_prep_inp_data->num_components == 1 ) {
                sprintf( szTitle, "%sInput Structure #%ld.%s%s%s%s%s",
                                  bStructurePreprocessed? "Preprocessed ":"",
                                  num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), iINChI? " (Reconnected)":"");
            } else {
                sprintf( szTitle, "Component #%d of %d, Input Structure #%ld.%s%s%s%s%s",
                                  i+1, cur_prep_inp_data->num_components,
                                  num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), iINChI? " (Reconnected)":"");
            }
#ifndef TARGET_LIB_FOR_WINCHI
            err = DisplayStructure( inp_cur_data->at, inp_cur_data->num_at,
                                    0, 1, 0, NULL, 1/*isotopic*/, 0/*taut*/, NULL, NULL,
                                    ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
            sd->bUserQuitComponentDisplay = (err==ESC_KEY);
            if ( !err ) 
            {
                inchi_fprintf( stderr, "Cannot display the structure\n");
            }
#else
            if(DRAWDATA && DRAWDATA_EXISTS)
            {
                struct DrawData vDrawData;
                int    nType = COMPONENT_ORIGINAL;
                vDrawData.pWindowData = CreateWinData_( inp_cur_data->at, inp_cur_data->num_at,
                                                        0, 1 /* bAdd_DT_to_num_H */, 0, NULL,
                                                        1 /* display isotopic if present */, 0, NULL, NULL,
                                                        ip->bAbcNumbers, &ip->dp, ip->nMode );
                if( vDrawData.pWindowData != NULL )
                {
                    if ( DRAWDATA_EXISTS ( i+1, nType, iINChI ) ) {  /* i = component number */
                        nType = COMPONENT_ORIGINAL_PREPROCESSED;
                    }
                    vDrawData.nComponent   = i+1;
                    vDrawData.nType        = nType;
                       vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                    vDrawData.szTitle              = _strdup(szTitle);
                    vDrawData.pWindowData->szTitle = _strdup(szTitle);
                    DRAWDATA(&vDrawData);
                }
            }
#endif
        }
#endif   /* } COMPILE_ANSI_ONLY */

#if ( TEST_RENUMB_ATOMS == 1 ) /* { */
        /****************************************************************************/
        /*     R E N U M B E R I N G (testing only) Part I  STARTS here             */
        /****************************************************************************/
        RenumberingTestInit( pRenumbData, inp_cur_data );
        if ( log_file != stderr ) {
            if ( ip->bDisplay )
                inchi_ios_eprint( log_file, "Component #%d structure #%ld.%s%s%s%s...\n", i+1, num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            else
                inchi_ios_eprint( stderr, "Component #%d structure #%ld.%s%s%s%s...\r", i+1, num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
        }
        /****************************************************************************/
        /*     R E N U M B E R I N G (testing only) Part I  ENDS here               */
        /****************************************************************************/
#endif  /* } TEST_RENUMB_ATOMS */


        /*******************************************************************************/
        /*                                                                             */
        /*  N O R M A L I Z A T I O N    a n d     C A N O N I C A L I Z A T I O N     */
        /*                                                                             */
        /*         (both tautomeric and non-tautomeric if requested)                   */
        /*                                                                             */
        /*******************************************************************************/
        /*  c) Create the component's INChI ( copies ip->bTautFlags into sd->bTautFlags)*/
        /*******************************************************************************/
        nRet = CreateOneComponentINChI( sd, ip, inp_cur_data, orig_inp_data, pINChI/*2[iINChI]*/, pINChI_Aux/*2[iINChI]*/, iINChI,
                                       i, num_inp, inp_norm_data, pncFlags, log_file );


#if ( TEST_RENUMB_ATOMS == 1 )  /*  { */
        /****************************************************************************/
        /*     R E N U M B E R I N G (testing only) Part II STARTS here             */
        /****************************************************************************/
        if ( !nRet ) {
            nRet = RenumberingTest( pINChI/*2[iINChI]*/, pINChI_Aux/*2[iINChI]*/, orig_inp_data, iINChI, pRenumbData, inp_cur_data, inp_norm_data, sd, ip, szTitle, log_file, prb_file, i, num_inp, pncFlags);
        }
        RenumberingTestUninit( pRenumbData );
        /****************************************************************************/
        /*     R E N U M B E R I N G (testing only)  Part II  ENDS here             */
        /****************************************************************************/
#endif  /*  }  TEST_RENUMB_ATOMS */


        /*  d) Display one component structure and/or INChI results only if there was no error */
#ifndef COMPILE_ANSI_ONLY /* { */
        if ( !nRet ) {
            /*  output one component INChI to the stdout if requested */
            /*
            if ( ip->bDisplayEachComponentINChI ) {
                int cur_num_non_taut = (pINChI[i][TAUT_NON] && pINChI[i][TAUT_NON]->nNumberOfAtoms>0);
                int cur_num_taut     = (pINChI[i][TAUT_YES] && pINChI[i][TAUT_YES]->nNumberOfAtoms>0);
                if ( ip->bDisplayEachComponentINChI && cur_num_non_taut + cur_num_taut ) {
                    SortAndPrintINChI(stdout, pStr, nStrLen, NULL, 
                                      ip, 1, cur_num_non_taut, cur_num_taut,
                                      num_inp, pINChI+i, pINChI_Aux+i, 
                                      save_opt_bits);
                }
            }
            */
            /**************************************************************************
             * display from one up to 4 structure pictures-results for each component *
             * Enable buttons:                                                        *
             * BN (non-tautomeric non-isotopic): inp_norm_data[0]->bExists            *
             * TN (tautomeric non-isotopic):     inp_norm_data[1]->bExists            *
             * BI (non-tautomeric isotopic):     inp_norm_data[0]->bExists &&         *
             *                                   inp_norm_data[0]->bHasIsotopicLayer  *
             * TI (tautomeric isotopic):         inp_norm_data[1]->bExists &&         *
             *                                   inp_norm_data[1]->bHasIsotopicLayer  *
             **************************************************************************/
            int bIsotopic, bTautomeric, bDisplayTaut, bHasIsotopicLayer, bFixedBondsTaut, m_max, m, nNumDisplayedFixedBondTaut=0;
            for ( j = 0; ip->bDisplay && !sd->bUserQuitComponentDisplay && j < TAUT_NUM; j ++ ) {
                if ( inp_norm_data[j]->bExists && !inp_norm_data[j]->bDeleted ) {
                    bTautomeric = (pINChI[i][j]->lenTautomer > 0); /* same as (inp_norm_data[j]->bTautomeric > 0) */
                    /* if requested tautomeric and no tautmerism found then do not say mobile or fixed H. 2004-10-27 */
                    bDisplayTaut = (!(ip->nMode & REQ_MODE_BASIC) && !bTautomeric)? -1 : bTautomeric;
                    bHasIsotopicLayer = (inp_norm_data[j]->bHasIsotopicLayer > 0);
                    for ( k = 0; k <= bHasIsotopicLayer; k ++ ) {
                        bIsotopic = (k > 0);
                        m_max = inp_norm_data[j]->at_fixed_bonds && inp_norm_data[j]->bTautPreprocessed? 1 : 0;
                        for ( m = m_max; 0 <= m; m -- ) {
                            bFixedBondsTaut = (m>0);
                            nNumDisplayedFixedBondTaut += bFixedBondsTaut; /* display only one time */
                            /*  added number of components, added another format for a single component case - DCh */
                            if ( cur_prep_inp_data->num_components > 1 ) {
                                sprintf( szTitle, "%s Component #%d of %d, Structure #%ld%s%s.%s%s%s%s%s",
                                              bFixedBondsTaut? "Preprocessed":"Result for",
                                              i+1, cur_prep_inp_data->num_components, num_inp,
                                              bDisplayTaut==1? ", mobile H": bDisplayTaut==0?", fixed H":"",
                                              bIsotopic? ", isotopic":"",
                                              SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), iINChI? " (Reconnected)":"");
                            } else {
                                sprintf( szTitle, "%s Structure #%ld%s%s.%s%s%s%s%s",
                                              bFixedBondsTaut? "Preprocessed":"Result for",
                                              num_inp,
                                              bDisplayTaut==1? ", mobile H": bDisplayTaut==0?", fixed H":"",
                                              bIsotopic? ", isotopic":"",
                                              SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue), iINChI? " (Reconnected)":"");
                            }
#ifndef TARGET_LIB_FOR_WINCHI
                            if ( bFixedBondsTaut && nNumDisplayedFixedBondTaut != 1 )
                                continue;
                            if ( bFixedBondsTaut ) {
                                err = DisplayStructure( inp_norm_data[j]->at_fixed_bonds, inp_norm_data[j]->num_at,
                                                        inp_norm_data[j]->num_removed_H, 0 /*bAdd_DT_to_num_H*/,
                                                        inp_norm_data[j]->nNumRemovedProtons,
                                                        inp_norm_data[j]->nNumRemovedProtonsIsotopic,
                                                        bHasIsotopicLayer, j, NULL,  NULL,
                                                        ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
                            } else {
                                err = DisplayStructure( inp_norm_data[j]->at, inp_norm_data[j]->num_at,
                                                        0, 0 /*bAdd_DT_to_num_H*/, 0, NULL,
                                                        k, j, pINChI[i], pINChI_Aux[i],
                                                        ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
                            }
                            if ( sd->bUserQuitComponentDisplay = (err==ESC_KEY) ) {
                                break;
                            }
#else
                            if(DRAWDATA && !bFixedBondsTaut)
                            {
                                struct DrawData vDrawData;
                                vDrawData.pWindowData = CreateWinData_( inp_norm_data[j]->at, inp_norm_data[j]->num_at,
                                                        0, 0 /* bAdd_DT_to_num_H */, 0, NULL,
                                                        k, j, pINChI[i], pINChI_Aux[i],
                                                        ip->bAbcNumbers, &ip->dp, ip->nMode );
                                if( vDrawData.pWindowData != NULL )
                                {
                                    int nType;
                                    vDrawData.nComponent = i+1;
                                    if( bTautomeric == 0 )
                                        nType = (bIsotopic == 0) ? COMPONENT_BN: COMPONENT_BI;
                                    else
                                        nType = (bIsotopic == 0) ? COMPONENT_TN: COMPONENT_TI;
                                    vDrawData.nType        = nType;
                                       vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                                    vDrawData.szTitle              = _strdup(szTitle);
                                    vDrawData.pWindowData->szTitle = _strdup(szTitle);
                                    DRAWDATA(&vDrawData);
                                }
                            } else
                            if(DRAWDATA && bFixedBondsTaut)
                            {
                                struct DrawData vDrawData;
                                if ( (ip->bCompareComponents & CMP_COMPONENTS) &&
                                     !(ip->bCompareComponents & CMP_COMPONENTS_NONTAUT) &&
                                     !bIsotopic == !inp_norm_data[j]->bHasIsotopicLayer ) {

                                    vDrawData.pWindowData =
                                         CreateWinData_( inp_norm_data[j]->at_fixed_bonds, inp_norm_data[j]->num_at,
                                                         inp_norm_data[j]->num_removed_H,
                                                         0 /* bAdd_DT_to_num_H */,
                                                         inp_norm_data[j]->nNumRemovedProtons,
                                                         inp_norm_data[j]->nNumRemovedProtonsIsotopic,
                                                         k, j, NULL, NULL,
                                                         ip->bAbcNumbers, &ip->dp, ip->nMode );
                                } else {
                                    continue;
                                }
                                if( vDrawData.pWindowData != NULL )
                                {
                                    vDrawData.nComponent = i+1;
                                    vDrawData.nType        = COMPONENT_ORIGINAL_PREPROCESSED;
                                       vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                                    vDrawData.szTitle              = _strdup(szTitle);
                                    vDrawData.pWindowData->szTitle = _strdup(szTitle);
                                    DRAWDATA(&vDrawData);
                                }
                            }
#endif
                        }
                    }
                }
            }

            /* save normalized components for composite display */
            if ( ip->bDisplayCompositeResults && all_inp_norm_data ) {
                for ( j = 0; j < TAUT_NUM; j ++ ) {
                    if ( inp_norm_data[j]->bExists ) {
                        all_inp_norm_data[i][j] = *inp_norm_data[j];
                        memset( inp_norm_data[j], 0, sizeof(*inp_norm_data[0]) );
                    }
                }
            }

        }
#endif  /* } COMPILE_ANSI_ONLY */
        if ( nRet ) {
            nRet = TreatCreateOneComponentINChIError(sd, ip, cur_prep_inp_data, i, num_inp,
                                 inp_file, log_file, output_file, prb_file,pStr, nStrLen );

            break;
        }
    }
    /**************************************************************************/
    /*                                                                        */
    /*                                                                        */
    /*   E N D   O F   T H E    M A I N   C Y C L E   P R O C E S S I N G     */
    /*                                                                        */
    /*          C O M P O N E N T S    O N E   B Y   O N E                    */
    /*                                                                        */
    /*                                                                        */
    /**************************************************************************/


exit_cycle:

#if ( TEST_RENUMB_ATOMS == 1 )  /* { */
    if ( pRenumbData->bRenumbErr && (!nRet || nRet==_IS_WARNING) ) {
        sd->nErrorCode = pRenumbData->bRenumbErr;
        nRet = TreatCreateOneComponentINChIError(sd, ip, cur_prep_inp_data, -1, num_inp,
                             inp_file, log_file, output_file, prb_file,pStr, nStrLen );
        /* nRet = _IS_ERROR; */
        sd->nErrorCode = 0;
        nRet           = 0;
    }
#endif /* } TEST_RENUMB_ATOMS */
    switch ( nRet ) {

    case _IS_FATAL:
    case _IS_ERROR:
        break;

    default:

#ifndef COMPILE_ANSI_ONLY /* { */
        /* composite results picture(s) */
        if ( all_inp_norm_data ) {
             int res = CreateCompositeNormAtom( composite_norm_data, all_inp_norm_data, pINChI, pINChI_Aux,
                                          prep_inp_data[iINChI].num_components, ip->nMode );
             /*
             for ( i = 0; i < prep_inp_data[iINChI].num_components; i ++ ) {
                 for ( k = 0; k < TAUT_NUM; k ++ ) {
                    FreeInpAtomData( &all_inp_norm_data[i][k] );
                 }
             }
             inchi_free( all_inp_norm_data );
             all_inp_norm_data = NULL;
             */
        }
#endif /* } COMPILE_ANSI_ONLY */

        break;
    }


#ifndef COMPILE_ANSI_ONLY /* { */
        /* avoid memory leaks in case of error */
        if ( all_inp_norm_data ) {
             for ( i = 0; i < prep_inp_data[iINChI].num_components; i ++ ) {
                 for ( k = 0; k < TAUT_NUM; k ++ ) {
                    FreeInpAtomData( &all_inp_norm_data[i][k] );
                 }
             }
             inchi_free( all_inp_norm_data );
             all_inp_norm_data = NULL;
        }
#endif /* } COMPILE_ANSI_ONLY */


    FreeInpAtomData( inp_cur_data     );
    for ( i = 0; i < TAUT_NUM; i ++ ) {
        FreeInpAtomData( inp_norm_data[i] );
    }


exit_function:

    return nRet;
}



#ifndef COMPILE_ANSI_ONLY /* { */
/****************************************************************************/
int CreateCompositeNormAtom( COMP_ATOM_DATA *composite_norm_data, INP_ATOM_DATA2 *all_inp_norm_data,
                             PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, int num_components, INCHI_MODE nMode )
{
    int i, j, jj, k, n, m, tot_num_at, tot_num_H, cur_num_at, cur_num_H, nNumRemovedProtons;
    int num_comp[TAUT_NUM+1], num_taut[TAUT_NUM+1], num_del[TAUT_NUM+1], num_at[TAUT_NUM+1], num_inp_at[TAUT_NUM+1];
    int ret = 0, indicator = 1;
    inp_ATOM *at, *at_from;
    memset( num_comp, 0, sizeof(num_comp) );
    memset( num_taut, 0, sizeof(num_taut) );
    memset( num_del, 0, sizeof(num_taut) );
    /* count taut and non-taut components */
    for ( j = 0; j < TAUT_NUM; j ++ ) {
        num_comp[j] = num_taut[j] = 0;
        for ( i = 0; i < num_components; i ++ ) {
            if ( all_inp_norm_data[i][j].bExists ) {
                num_del[j]  += (0 != all_inp_norm_data[i][j].bDeleted );
                num_comp[j] ++;
                num_taut[j] += (0 != all_inp_norm_data[i][j].bTautomeric);
            }
        }
    }
    /* count intermediate taut structure components */
    if ( num_comp[TAUT_YES] > num_del[TAUT_YES] && num_taut[TAUT_YES] ) {
        /*
        num_comp[TAUT_INI] = num_comp[TAUT_YES] - num_del[TAUT_YES];
        */

        for ( i = 0, j=TAUT_YES; i < num_components; i ++ ) {
            if ( all_inp_norm_data[i][j].bExists &&
                (all_inp_norm_data[i][j].bDeleted ||
                 all_inp_norm_data[i][j].bTautomeric &&
                 all_inp_norm_data[i][j].at_fixed_bonds &&
                 all_inp_norm_data[i][j].bTautPreprocessed) ) {
                num_comp[TAUT_INI] ++;
            }
        }

    }
    /* count atoms and allocate composite atom data */
    for ( jj = 0; jj <= TAUT_INI; jj ++ ) {
        num_at[jj] = num_inp_at[jj] = 0;
        j = inchi_min (jj, TAUT_YES);
        if ( num_comp[jj] ) {
            for ( i = 0; i < num_components; i ++ ) {
                if ( all_inp_norm_data[i][j].bDeleted )
                    continue;
                /* find k = the normaized structure index */
                if ( jj == TAUT_INI ) {
                    if ( all_inp_norm_data[i][j].bExists &&
                         all_inp_norm_data[i][j].at_fixed_bonds ) {
                        k = j;
                    } else
                    if ( all_inp_norm_data[i][ALT_TAUT(j)].bExists && !all_inp_norm_data[i][ALT_TAUT(j)].bDeleted &&
                         !all_inp_norm_data[i][j].bDeleted  ) {
                        k = ALT_TAUT(j);
                    } else
                    if ( all_inp_norm_data[i][j].bExists ) {
                        k = j;
                    } else {
                        continue;
                    }
                } else {
                    if ( all_inp_norm_data[i][j].bExists ) {
                        k = j;
                    } else
                    if ( all_inp_norm_data[i][ALT_TAUT(j)].bExists && !all_inp_norm_data[i][ALT_TAUT(j)].bDeleted) {
                        k = ALT_TAUT(j);
                    } else {
                        continue;
                    }
                }
                num_inp_at[jj] += all_inp_norm_data[i][k].num_at; /* all atoms including terminal H */
                num_at[jj]     += all_inp_norm_data[i][k].num_at - all_inp_norm_data[i][k].num_removed_H;
            }
            if ( num_inp_at[jj] ) {
                if ( !CreateCompAtomData( composite_norm_data+jj, num_inp_at[jj], num_components, jj == TAUT_INI ) )
                    goto exit_error;
                composite_norm_data[jj].num_removed_H = num_inp_at[jj] - num_at[jj];
            }
        }
    }
    /* fill out composite atom */
    for ( jj = 0; jj <= TAUT_INI; jj ++, indicator <<= 1 ) {
        j = inchi_min (jj, TAUT_YES);
        if ( num_comp[jj] ) {
            tot_num_at = 0;
            tot_num_H = 0;
            for ( i = 0; i < num_components; i ++ ) {
                if ( all_inp_norm_data[i][j].bDeleted ) {
                    composite_norm_data[jj].nNumRemovedProtons += all_inp_norm_data[i][j].nNumRemovedProtons;
                    for ( n = 0; n < NUM_H_ISOTOPES; n ++ ) {
                        composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][j].nNumRemovedProtonsIsotopic[n];
                    }
                    continue;
                }
                nNumRemovedProtons = 0;
                k = TAUT_NUM;
                /* find k = the normaized structure index */
                if ( jj == TAUT_INI ) {
                    if ( all_inp_norm_data[i][j].bExists && all_inp_norm_data[i][j].at_fixed_bonds ) {
                        k = j;
                    } else
                    if ( all_inp_norm_data[i][ALT_TAUT(j)].bExists ) {
                        k = ALT_TAUT(j);
                    } else
                    if ( all_inp_norm_data[i][j].bExists && !all_inp_norm_data[i][ALT_TAUT(j)].bDeleted ) {
                        k = j;
                    } else {
                        continue;
                    }
                } else {
                    if ( all_inp_norm_data[i][j].bExists ) {
                        k = j;
                    } else
                    if ( all_inp_norm_data[i][ALT_TAUT(j)].bExists && !all_inp_norm_data[i][ALT_TAUT(j)].bDeleted ) {
                        k = ALT_TAUT(j);
                    } else {
                        continue;
                    }
                }
                /* copy main atoms */
                cur_num_H  = all_inp_norm_data[i][k].num_removed_H;       /* number of terminal H atoms */
                cur_num_at = all_inp_norm_data[i][k].num_at - cur_num_H;  /* number of all but explicit terminal H atoms */

                if ( (tot_num_at + cur_num_at) > num_at[jj] ||
                     (num_at[jj] + tot_num_H + cur_num_H) > num_inp_at[jj] ) {
                    goto exit_error; /* miscount */
                }
                at      = composite_norm_data[jj].at+tot_num_at; /* points to the 1st destination atom */
                at_from = (jj == TAUT_INI && k == TAUT_YES && all_inp_norm_data[i][k].at_fixed_bonds)?
                          all_inp_norm_data[i][k].at_fixed_bonds : all_inp_norm_data[i][k].at;
                memcpy( at, at_from, sizeof(composite_norm_data[0].at[0]) * cur_num_at ); /* copy atoms except terminal H */
                /* shift neighbors of main atoms */
                for ( n = 0; n < cur_num_at; n ++, at ++ ) {
                    for ( m = 0; m < at->valence; m ++ ) {
                        at->neighbor[m] += tot_num_at;
                    }
                }
                /* copy explicit H */
                if ( cur_num_H ) {
                    at = composite_norm_data[jj].at+num_at[jj]+tot_num_H; /* points to the 1st destination atom */
                    memcpy( at, at_from+cur_num_at,
                            sizeof(composite_norm_data[0].at[0]) * cur_num_H );
                    /* shift neighbors of explicit H atoms */
                    for ( n = 0; n < cur_num_H; n ++, at ++ ) {
                        for ( m = 0; m < at->valence; m ++ ) {
                            at->neighbor[m] += tot_num_at;
                        }
                    }
                }
                /* composite counts */
                composite_norm_data[jj].bHasIsotopicLayer   |= all_inp_norm_data[i][k].bHasIsotopicLayer;
                composite_norm_data[jj].num_isotopic        += all_inp_norm_data[i][k].num_isotopic;
                composite_norm_data[jj].num_bonds           += all_inp_norm_data[i][k].num_bonds;
                composite_norm_data[jj].bTautomeric         += (j == jj) && all_inp_norm_data[i][k].bTautomeric;
                composite_norm_data[jj].nNumRemovedProtons  += all_inp_norm_data[i][k].nNumRemovedProtons;
                for ( n = 0; n < NUM_H_ISOTOPES; n ++ ) {
                    composite_norm_data[jj].nNumRemovedProtonsIsotopic[n] += all_inp_norm_data[i][k].nNumRemovedProtonsIsotopic[n];
                    composite_norm_data[jj].num_iso_H[n]                  += all_inp_norm_data[i][k].num_iso_H[n];
                }
                /*
                composite_norm_data[j].num_at            += cur_num_at + cur_num_H;
                composite_norm_data[j].num_removed_H     += cur_num_H;
                */
                /* total count */
                tot_num_at += cur_num_at;
                tot_num_H += cur_num_H;
                /* offset for the next component */
                if ( composite_norm_data[jj].nOffsetAtAndH ) {
                    composite_norm_data[jj].nOffsetAtAndH[2*i]   = tot_num_at;
                    composite_norm_data[jj].nOffsetAtAndH[2*i+1] = num_at[jj]+tot_num_H;
                }
            }
            if ( tot_num_at != num_at[jj] ||
                 num_at[jj] + tot_num_H  != num_inp_at[jj] ) {
                goto exit_error; /* miscount */
            }
            composite_norm_data[jj].bExists       = (tot_num_at>0);
            ret |= indicator;
        }
    }
    return ret;





exit_error:
    return ret;
}
#endif /* } COMPILE_ANSI_ONLY */
