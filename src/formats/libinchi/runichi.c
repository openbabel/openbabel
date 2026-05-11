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


/*
    General processing procedures

*/

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>
#include <locale.h>

#include "mode.h"
#include "ichitime.h"
#ifndef COMPILE_ANSI_ONLY
#include <conio.h>
#endif

#include "ichimain.h"
#include "ichi_io.h"
#include "mol_fmt.h"
#include "inchi_api.h"
#include "readinch.h"
#include "ichicant.h"
#ifdef TARGET_LIB_FOR_WINCHI
#include "../../../IChI_lib/src/ichi_lib.h"
#include "inchi_api.h"
#endif
#include "inchi_gui.h"
#include "readinch.h"
#include "ichirvrs.h"

#include "bcf_s.h"

extern int DisplayTheWholeStructure( struct tagCANON_GLOBALS *pCG,
                                     struct tagINCHI_CLOCK   *ic,
                                     struct tagStructData    *sd,
                                     INPUT_PARMS             *ip,
                                     char                    *szTitle,
                                     INCHI_IOSTREAM          *inp_file,
                                     INCHI_IOSTREAM          *log_file,
                                     ORIG_ATOM_DATA          *orig_inp_data,
                                     long                    num_inp,
                                     int                     iINChI,
                                     int                     bShowStruct,
                                     int                     bINCHI_LIB_Flag );

extern int DisplayStructure( struct tagCANON_GLOBALS *pCG,
                             inp_ATOM   *at,
                             int        num_at,
                             OAD_Polymer *polymer,
                             int        num_removed_H,
                             int        bAdd_DT_to_num_H,
                             int        nNumRemovedProtons,
                             NUM_H      *nNumRemovedProtonsIsotopic,
                             int        bIsotopic,
                             int        j /*bTautomeric*/,
                             INChI      **cur_INChI,
                             INChI_Aux  **cur_INChI_Aux,
                             int        bAbcNumbers,
                             DRAW_PARMS *dp,
                             INCHI_MODE nMode,
                             char       *szTitle );

extern int ReadInChICoord( INCHI_IOSTREAM   *pInp,
                           SEGM_LINE        *pLine,
                           int              *pState,
                           INChI            *pInpInChI[INCHI_NUM][TAUT_NUM],
                           int              nNumComponents[INCHI_NUM][TAUT_NUM] );


/* Local functions */
static int DoOneStructureEarlyPreprocessing( INCHI_CLOCK *ic,
                                             CANON_GLOBALS *pCG,
                                             long num_inp,
                                             STRUCT_DATA *sd,
                                             INPUT_PARMS *ip,
                                             INCHI_IOSTREAM *inp_file,
                                             INCHI_IOSTREAM *log_file,
                                             INCHI_IOSTREAM *out_file,
                                             INCHI_IOSTREAM *prb_file,
                                             ORIG_ATOM_DATA *orig_inp_data,
                                             ORIG_ATOM_DATA *prep_inp_data );

/* Actual worker sitting under ProcessOneStructureEx */
int ProcessOneStructureExCore( struct tagINCHI_CLOCK    *ic,
                               struct tagCANON_GLOBALS  *CG,
                               STRUCT_DATA              *sd,
                               INPUT_PARMS              *ip,
                               char                     *szTitle,
                               PINChI2                  *pINChI2[INCHI_NUM],
                               PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                               INCHI_IOSTREAM           *inp_file,
                               INCHI_IOSTREAM           *log_file,
                               INCHI_IOSTREAM           *out_file,
                               INCHI_IOSTREAM           *prb_file,
                               ORIG_ATOM_DATA           *orig_inp_data,
                               ORIG_ATOM_DATA           *prep_inp_data,
                               long                     num_inp,
                               INCHI_IOS_STRING         *strbuf,
                               unsigned char            save_opt_bits );



static ORIG_STRUCT*
OrigAtData_StoreNativeInput( CANON_GLOBALS  *pCG,
                             int            *nRet,
                             STRUCT_DATA    *sd,
                             INPUT_PARMS    *ip,
                             ORIG_ATOM_DATA *orig_inp_data,
                             ORIG_STRUCT    *pOrigStruct );
static void PrepareSaveOptBits( unsigned char *save_opt_bits,
                                INPUT_PARMS   *ip );

static void DisplayOrigAndResultStructuresAndComponents( int nRet,
                                                         INCHI_CLOCK    *ic,
                                                         CANON_GLOBALS  *pCG,
                                                         STRUCT_DATA    *sd,
                                                         INPUT_PARMS    *ip,
                                                         char           *szTitle,
                                                         PINChI2        *pINChI[INCHI_NUM],
                                                         PINChI_Aux2    *pINChI_Aux[INCHI_NUM],
                                                         INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file,
                                                         INCHI_IOSTREAM *out_file,
                                                         ORIG_ATOM_DATA *orig_inp_data,
                                                         ORIG_ATOM_DATA *prep_inp_data,
                                                         long           num_inp,
                                                         int            maxINChI,
                                                         COMP_ATOM_DATA composite_norm_data[INCHI_NUM][TAUT_NUM + 1] );
static void SaveOkProcessedMolfile( int             nRet,
                                    STRUCT_DATA     *sd,
                                    INPUT_PARMS     *ip,
                                    INCHI_IOSTREAM  *prb_file,
                                    INCHI_IOSTREAM  *inp_file );

static int mark_atoms_to_delete_or_renumber( ORIG_ATOM_DATA *orig_at_data,
                                             OAD_StructureEdits *ed,
                                             int *at_renum);
static int set_renumbered_or_delete( int *number,
                                     int *buf,
                                     int nelems,
                                     int *renum,
                                     int base );

static void swap_atoms_xyz( ORIG_ATOM_DATA *orig_at_data, int ia1, int ia2 );

/* Callbacks */

/*  Console user issued CTRL+C etc. */
int( *ConsoleQuit )( void ) = NULL;
int( *UserAction ) ( void ) = NULL;

#if (BUILD_WITH_ENG_OPTIONS==1)
#if ALLOW_SUBSTRUCTURE_FILTERING==1
static int OrigAtData_CheckForSubstructure(ORIG_ATOM_DATA *orig_inp_data);
#endif
#endif


/**********************************************
 * output " L=V" or " L missing" or ""
 * The fprintf format string must contain %s%s%s%s
 */
const char gsMissing[] = "is missing";
const char gsEmpty[] = "";
const char gsSpace[] = " ";
const char gsEqual[] = "=";



/****************************************************************************
 Process a portion of input data (molecule, InChI string, ...)
 in a relevant way (generate InChI, restore molecule by InChI )
****************************************************************************/
int ProcessOneStructure( INCHI_CLOCK            *ic,
                         CANON_GLOBALS          *pCG,
                         STRUCT_DATA            *sd,
                         INPUT_PARMS            *ip,
                         char                   *szTitle,
                         PINChI2                *pINChI[INCHI_NUM],
                         PINChI_Aux2            *pINChI_Aux[INCHI_NUM],
                         INCHI_IOSTREAM         *inp_file,
                         INCHI_IOSTREAM         *log_file,
                         INCHI_IOSTREAM         *out_file,
                         INCHI_IOSTREAM         *prb_file,
                         ORIG_ATOM_DATA         *orig_inp_data,
                         ORIG_ATOM_DATA         *prep_inp_data,
                         long                   num_inp,
                         INCHI_IOS_STRING       *strbuf,
                         unsigned char          save_opt_bits )
{
    int nRet = 0,
        nRet1, i, k,
        maxINChI = 0,
        bSortPrintINChIFlags = 0;
    COMP_ATOM_DATA
        composite_norm_data[INCHI_NUM][TAUT_NUM + 1];    /*    [0]:non-taut,
                                                        [1]:taut,
                                                        [2]:intermediate taut struct */
    NORM_CANON_FLAGS ncFlags;
    NORM_CANON_FLAGS *pncFlags = &ncFlags;
    ORIG_STRUCT OrigStruct;
    ORIG_STRUCT *pOrigStruct = NULL;
    int err, ret1 = 0;

    /* djb-rwth: removing redundant code */
#ifdef GHI100_FIX
#if ((SPRINTF_FLAG != 1) && (SPRINTF_FLAG != 2))
    setlocale(LC_ALL, "en-US"); /* djb-rwth: setting all locales to "en-US" */
#endif
#endif

    /*    1. Preliminary work */

    /* djb-rwth: fixing coverity ID #499508 */
    if (!orig_inp_data)
    {
        goto exit_function;
    }

    int is_polymer = orig_inp_data->valid_polymer
                     && orig_inp_data->polymer
                     && orig_inp_data->polymer->n ;

    int is_polymer2inchi = is_polymer && ( ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE );

    sd->bUserQuitComponent = 0;
    sd->bUserQuitComponentDisplay = 0;
    memset( composite_norm_data, 0, sizeof( composite_norm_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( pncFlags, 0, sizeof( *pncFlags ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        
    /* For experimental purposes only */
    /*ret1 = DoOneStructureEarlyPreprocessing( num_inp, sd, ip, inp_file,
                                             log_file, out_file, prb_file,
                                             orig_inp_data, prep_inp_data );
    */
    ret1 = DoOneStructureEarlyPreprocessing( ic, pCG, num_inp, sd, ip,
                                             inp_file, log_file, out_file, prb_file,
                                             orig_inp_data, prep_inp_data );
    switch (ret1)
    {
        case _IS_SKIP:
        case _IS_ERROR:
        case _IS_FATAL:
            nRet = ret1;
    }

    if (ret1)
    {
        goto exit_function;
    }


    if (is_polymer)
    {
        /* Polymer house-keeping related to possible CRU frame shift(s) */

        orig_inp_data->polymer->frame_shift_scheme = ip->bFrameShiftScheme;
        orig_inp_data->polymer->treat = ip->bPolymers;

        if (!is_polymer2inchi)
        {
            /* Polymer structure is being restored from InChI string            */
            /* If CRUs were pre-cyclized, re-open them in preferred forms here  */
            if (orig_inp_data->polymer->frame_shift_scheme == FSS_STARS_CYCLED)
            {
                OAD_Polymer_SmartReopenCyclizedUnits( orig_inp_data->polymer,
                                                      orig_inp_data->at,
                                                      orig_inp_data->num_inp_atoms,
                                                      &orig_inp_data->num_inp_bonds );
            }
        }
    }


    ret1 = OrigAtData_SaveMolfile( orig_inp_data, sd, ip, num_inp, out_file );
    if (ret1)
    {
        goto exit_function;
    }


    pOrigStruct = &OrigStruct;
    memset( pOrigStruct, 0, sizeof( *pOrigStruct ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    OrigAtData_StoreNativeInput( pCG, &nRet, sd,  ip,  orig_inp_data, pOrigStruct );

    /*    2. Create INChI for the whole disconnected or original structure */

    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        nRet1 = CreateOneStructureINChI( pCG, ic, sd, ip, szTitle,
                                         pINChI, pINChI_Aux, INCHI_BAS,
                                         inp_file, log_file, out_file, prb_file,
                                         orig_inp_data, prep_inp_data,
                                         composite_norm_data,
                                         num_inp, strbuf, pncFlags );
        nRet = inchi_max( nRet, nRet1 );

        /* If we create InChI from polymer-containing structure */
        if (is_polymer2inchi)
        {
            int polymer_repr_type = OAD_Polymer_GetRepresentation( orig_inp_data->polymer );

#ifdef ALLOW_MIXED_SRU_AND_MON
            if (polymer_repr_type == POLYMER_REPRESENTATION_STRUCTURE_BASED ||
                 polymer_repr_type == POLYMER_REPRESENTATION_MIXED)
#else
            if (polymer_repr_type == POLYMER_REPRESENTATION_STRUCTURE_BASED)
#endif
            {
                /* Temporarily copy ptr to polymer data to prep_inp_data */
                OAD_Polymer *prep_polymer = prep_inp_data->polymer; /* may be NULL */
                prep_inp_data->polymer = orig_inp_data->polymer;

                OAD_Polymer_FindBackbones( prep_inp_data, /* NB: not orig_inp_data! */
                                           &( composite_norm_data[INCHI_BAS][TAUT_YES] ),
                                           &err, sd->pStrErrStruct );
                if (err)
                {
                    ret1 = _IS_ERROR;
                }
                nRet = inchi_max( nRet, ret1 );
                prep_inp_data->polymer = prep_polymer;    /* restore temp copied*/
            }
        }
    }

    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        maxINChI = 1;
    }


    /* 3. Create INChI for the whole metal-reconnected structure */

    if (nRet != _IS_FATAL                                              &&
         nRet != _IS_ERROR &&
         ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
         ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ))
    {

        nRet1 = CreateOneStructureINChI( pCG, ic, sd, ip, szTitle,
                                         pINChI, pINChI_Aux, INCHI_REC, inp_file,
                                         log_file, out_file, prb_file,
                                         orig_inp_data, prep_inp_data,
                                         composite_norm_data,
                                         num_inp, strbuf, pncFlags );
        nRet = inchi_max( nRet, nRet1 );

        if (is_polymer2inchi)
        {
            ret1 = 0;
            /* temporarily copy ptr to polymer data to prep_inp_data */
            prep_inp_data->polymer = orig_inp_data->polymer;

            OAD_Polymer_FindBackbones( prep_inp_data, /* NB: not orig_inp_data! */
                                       &( composite_norm_data[INCHI_REC][TAUT_YES] ),
                                       &err, sd->pStrErrStruct );
            if (err)
            {
                ret1 = _IS_ERROR;
            }
            nRet = inchi_max( nRet, ret1 );
            prep_inp_data->polymer = NULL;    /* remove temp copied */
        }
        if (nRet != _IS_FATAL && nRet != _IS_ERROR)
        {
            maxINChI = 2;
        }
    }

    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        if (( sd->bChiralFlag & FLAG_INP_AT_CHIRAL ) &&
            ( ip->nMode & REQ_MODE_STEREO ) &&
             !( ip->nMode & ( REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO ) ) &&
             !bIsStructChiral( pINChI, sd->num_components ))
        {
            if (!ip->bNoWarnings)
            {
                WarningMessage( sd->pStrErrStruct, "Not chiral" );
            }
        }
        if (!sd->bUserQuitComponent && !sd->bUserQuit)
        {
            nRet1 = TreatCreateINChIWarning( sd, ip, prep_inp_data, num_inp,
                                             inp_file, log_file, out_file, prb_file );
            nRet = inchi_max( nRet, nRet1 );
        }
    }


    /*    4. Sort and print INChI for the whole structure */

    PrepareSaveOptBits( &save_opt_bits, ip );
    if (nRet != _IS_FATAL && nRet != _IS_ERROR)
    {
        nRet1 = SortAndPrintINChI( pCG, out_file, strbuf, log_file, ip,
                                   orig_inp_data, prep_inp_data,
                                   composite_norm_data,
                                   pOrigStruct, sd->num_components,
                                   sd->num_non_taut, sd->num_taut,
                                   sd->bTautFlags, sd->bTautFlagsDone,
                                   pncFlags, num_inp,
                                   pINChI, pINChI_Aux,
                                   &bSortPrintINChIFlags, save_opt_bits ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    }


    /*    5. Post-process */

    DisplayOrigAndResultStructuresAndComponents( nRet, ic, pCG, sd, ip, szTitle,
                                                 pINChI, pINChI_Aux,
                                                 inp_file, log_file, out_file,
                                                 orig_inp_data, prep_inp_data,
                                                 num_inp, maxINChI,
                                                 composite_norm_data );


    SaveOkProcessedMolfile( nRet, sd, ip, prb_file, inp_file );


    /* Cleanup */

    for (i = 0; i < INCHI_NUM; i++)
    {
        for (k = 0; k < TAUT_NUM + 1; k++)
        {
            FreeCompAtomData( &composite_norm_data[i][k] );
        }
    }

    OrigStruct_Free( pOrigStruct );

exit_function:

    return nRet;
}



/****************************************************************************
 Early preprocessing: used if defined
    REMOVE_ION_PAIRS_ORIG_STRU or UNDERIVATIZE or RING2CHAIN
****************************************************************************/


int DoOneStructureEarlyPreprocessing( INCHI_CLOCK *ic,
                                      CANON_GLOBALS *pCG,
                                      long num_inp,
                                      STRUCT_DATA *sd,
                                      INPUT_PARMS *ip,
                                      INCHI_IOSTREAM *inp_file,
                                      INCHI_IOSTREAM *log_file,
                                      INCHI_IOSTREAM *out_file,
                                      INCHI_IOSTREAM *prb_file,
                                      ORIG_ATOM_DATA *orig_inp_data,
                                      ORIG_ATOM_DATA *prep_inp_data )
{

#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )
    int ret1 = 0, ret2 = 0; /* djb-rwth: removing redundant variables */
#endif

#if ( REMOVE_ION_PAIRS_ORIG_STRU == 1 )
    fix_odd_things( orig_inp_data->num_inp_atoms, orig_inp_data->at, 0, ip->bFixNonUniformDraw );
#endif

#if ( UNDERIVATIZE == 1 )
    if (ip->bUnderivatize)
    {
        if (0 > ( ret2 = OAD_Edit_Underivatize( ic, pCG, orig_inp_data, ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ), ip->bUnderivatize & 2, ip->pSdfValue ) ))
        {
            long num_inp2 = num_inp;
            AddErrorMessage( sd->pStrErrStruct, "Underivatization error" );
            sd->nStructReadError = 99;
            sd->nErrorType = _IS_ERROR;
            TreatErrorsInReadTheStructure( sd, ip, LOG_MASK_ALL, inp_file, log_file, out_file, prb_file,
                                        prep_inp_data, &num_inp2 );
            return _IS_ERROR; /* output only if derivatives found */
        }
        else if (0 < ret2)
        {
            if (!ip->bNoWarnings)
            {
                WarningMessage( sd->pStrErrStruct, "Input structure underivatized" );
            }
        }
    }
#endif /* UNDERIVATIZE == 1 */

#if ( RING2CHAIN == 1 )
    if (ip->bRing2Chain && 0 > ( ret1 = Ring2Chain( ic, pCG, orig_inp_data ) ))
    {
        long num_inp2 = num_inp;
        AddErrorMessage( sd->pStrErrStruct, "Ring to chain error" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_ERROR;
        /* djb-rwth: removing redundant code */
        TreatErrorsInReadTheStructure( sd, ip, LOG_MASK_ALL,
                                       inp_file, log_file,
                                       out_file, prb_file,
                                       prep_inp_data, &num_inp2 );
        return _IS_ERROR; /* output only if derivatives found */
    }
#endif /* RING2CHAIN == 1 */
#if ( RING2CHAIN == 1 || UNDERIVATIZE == 1 )  /***** post v.1 feature *****/
    if (ip->bIgnoreUnchanged && !ret1 && !ret2)
    {
        return _IS_SKIP; /* output only if derivatives or ring/chain found */
    }
#endif /* RING2CHAIN == 1 || UNDERIVATIZE == 1 */
    return 0;
}


/* If requested, save input data to a Molfile instead of creating INChI                     */
/* Also used for output in case of combination of options 'InChI2Struct' and 'OutputSDF'    */
int OrigAtData_SaveMolfile( ORIG_ATOM_DATA  *orig_inp_data,
                            STRUCT_DATA     *sd,
                            INPUT_PARMS     *ip,
                            long            num_inp,
                            INCHI_IOSTREAM  *out_file )
{
    int ret = 0;

    if (!( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ))
    {
        return _IS_OKAY;
    }
    else
    {
        char szNumber[256];
        sprintf(szNumber, "Structure #%ld. %s%s%s%s", num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));
        ret = OrigAtData_WriteToSDfile( orig_inp_data, out_file, szNumber, NULL,
                                        ( sd->bChiralFlag&FLAG_INP_AT_CHIRAL ) ? 1 : 0,
                                        ( ip->bINChIOutputOptions&INCHI_OUT_SDFILE_ATOMS_DT ) ? 1 : 0,
                                        ip->pSdfLabel, ip->pSdfValue );
    }

    return ret;
}


/****************************************************************************
 Optionally save native input data as 'OrigStruct' data package
****************************************************************************/
ORIG_STRUCT * OrigAtData_StoreNativeInput( CANON_GLOBALS    *pCG,
                                           int              *nRet,
                                           STRUCT_DATA      *sd,
                                           INPUT_PARMS      *ip,
                                           ORIG_ATOM_DATA   *orig_inp_data,
                                           ORIG_STRUCT      *pOrigStruct )
{

    /*    v. 1.05 always create and fill OrigStruc as it may be used to store e.g. polymer info    */
    /*    If normal AuxInfo is requested, create full reversibility information from native inp data
    if ( ip->bINChIOutputOptions & (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO))
        return NULL; */

    if (OrigStruct_FillOut( pCG, orig_inp_data, pOrigStruct, sd ))
    {
        AddErrorMessage( sd->pStrErrStruct, "Cannot interpret reversibility information" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_ERROR;
        *nRet = _IS_ERROR;
    }

    return pOrigStruct;
}


/****************************************************************************
 Prepare SaveOpt bits
****************************************************************************/
void PrepareSaveOptBits( unsigned char *save_opt_bits, INPUT_PARMS *ip )
{
    if (ip->nInputType != INPUT_INCHI)
    {
        *save_opt_bits = 0;
        if (ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT)
        {
            if (0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ))
            {
                ( *save_opt_bits ) |= SAVE_OPT_RECMET;
            }
            if (0 != ( ip->nMode & REQ_MODE_BASIC ))
            {
                ( *save_opt_bits ) |= SAVE_OPT_FIXEDH;
            }
            if (0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO ))
            {
                ( *save_opt_bits ) |= SAVE_OPT_SLUUD;
            }
            if (0 == ( ip->nMode & ( REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU ) ))
            {
                ( *save_opt_bits ) |= SAVE_OPT_SUU;
            }
            if (0 != ( ip->bTautFlags & TG_FLAG_KETO_ENOL_TAUT ))
            {
                ( *save_opt_bits ) |= SAVE_OPT_KET;
            }
            if (0 != ( ip->bTautFlags & TG_FLAG_1_5_TAUT ))
            {
                ( *save_opt_bits ) |= SAVE_OPT_15T;
            }
            /* djb-rwth: addressing coverity ID #499536 -- despite different bit-sizes, works properly */
            if (0 != (ip->bTautFlags & TG_FLAG_PT_22_00))
                (*save_opt_bits) |= SAVE_OPT_PT_22_00;
            if (0 != (ip->bTautFlags & TG_FLAG_PT_16_00))
                (*save_opt_bits) |= SAVE_OPT_PT_16_00;
            if (0 != (ip->bTautFlags & TG_FLAG_PT_06_00))
                (*save_opt_bits) |= SAVE_OPT_PT_06_00;
            if (0 != (ip->bTautFlags & TG_FLAG_PT_39_00))
                (*save_opt_bits) |= SAVE_OPT_PT_39_00;
            if (0 != (ip->bTautFlags & TG_FLAG_PT_13_00))
                (*save_opt_bits) |= SAVE_OPT_PT_13_00;
            if (0 != (ip->bTautFlags & TG_FLAG_PT_18_00))
                (*save_opt_bits) |= SAVE_OPT_PT_18_00;
            /* Check if /SNon requested and turn OFF stereo bits if so */
            if (!( ip->nMode & REQ_MODE_STEREO ))
            {
                ( *save_opt_bits ) &= ~SAVE_OPT_SUU;
                ( *save_opt_bits ) &= ~SAVE_OPT_SLUUD;
            }
        }
    }
}


/****************************************************************************
 Display structures/components on screen
****************************************************************************/
void DisplayOrigAndResultStructuresAndComponents( int               nRet,
                                                  INCHI_CLOCK       *ic,
                                                  CANON_GLOBALS     *pCG,
                                                  STRUCT_DATA       *sd,
                                                  INPUT_PARMS       *ip,
                                                  char              *szTitle,
                                                  PINChI2           *pINChI[INCHI_NUM],
                                                  PINChI_Aux2       *pINChI_Aux[INCHI_NUM],
                                                  INCHI_IOSTREAM    *inp_file,
                                                  INCHI_IOSTREAM    *log_file,
                                                  INCHI_IOSTREAM    *out_file,
                                                  ORIG_ATOM_DATA    *orig_inp_data,
                                                  ORIG_ATOM_DATA    *prep_inp_data,
                                                  long              num_inp,
                                                  int               maxINChI,
                                                  COMP_ATOM_DATA    composite_norm_data[INCHI_NUM][TAUT_NUM + 1] )
{


    if (ip->bDisplay)    ip->bDisplayCompositeResults = 1;    /* v. 1.05 */

#ifndef COMPILE_ANSI_ONLY /* { */

    /* Display equivalent components on original or preprocessed structure(s) */
#ifndef TARGET_LIB_FOR_WINCHI
    if (nRet != _IS_FATAL && nRet != _IS_ERROR && /*ip->bDisplay &&*/
        ( ip->bCompareComponents & CMP_COMPONENTS ) && !sd->bUserQuit && !sd->bUserQuitComponent)
    {
        int j, ret, ord;
        int bDisplaySaved = ip->bDisplay;
        ORIG_ATOM_DATA *inp_data;
        AT_NUMB         nEquSet;
        for (ord = -1; ord < INCHI_NUM; ord++)
        {
            switch (ord)
            {
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
            inp_data = j < 0 ? orig_inp_data : prep_inp_data + j;
            if (inp_data && inp_data->num_inp_atoms && inp_data->at &&
                 inp_data->nEquLabels &&
                 inp_data->nNumEquSets)
            {
                for (nEquSet = 1; nEquSet <= inp_data->nNumEquSets; nEquSet++)
                {
                    ip->dp.nEquLabels = inp_data->nEquLabels;
                    ip->dp.nCurEquLabel = nEquSet;
                    ip->dp.nNumEquSets = inp_data->nNumEquSets;
                    ip->bDisplay = 1; /* force display if it was not requested */
                    ret = DisplayTheWholeStructure( pCG, ic, sd, ip, szTitle,
                                                    inp_file, log_file, inp_data,
                                                    num_inp, j, 1 /*bShowStructure*/, 0 );
                    ip->dp.nEquLabels = NULL;
                    ip->dp.nCurEquLabel = 0;
                    ip->dp.nNumEquSets = 0;
                    ip->bDisplay = bDisplaySaved; /* restore display option */
                    if (ret)
                    {
            /* user pressed Esc */
                        goto exit_loop;
                    }
                }
            }
        }
    exit_loop:;
    }
#endif

    /* Display composite results and equivalent components on composite results */
    if (nRet != _IS_FATAL && nRet != _IS_ERROR && /*ip->bDisplay &&*/ ip->bDisplayCompositeResults)
    {
        int iINChI;
        for (iINChI = 0; iINChI < maxINChI && !sd->bUserQuitComponentDisplay; iINChI++)
        {
            DisplayTheWholeCompositeStructure( pCG, ic, ip, sd, num_inp,
                                               iINChI, pINChI[iINChI], pINChI_Aux[iINChI],
                                               orig_inp_data, prep_inp_data, composite_norm_data[iINChI] );
        }
#ifndef TARGET_LIB_FOR_WINCHI
        if (!ip->bDisplay && sd->bUserQuitComponentDisplay)
        {
            sd->bUserQuit = 1;
        }
#endif
    }

#endif /* } COMPILE_ANSI_ONLY */

    return;
}


/****************************************************************************
 Special mode (option /PGO) : extract all good MOLfiles into the problem file;
 do not extract any MOLfile that could not be processed.
****************************************************************************/
void SaveOkProcessedMolfile( int            nRet,
                             STRUCT_DATA    *sd,
                             INPUT_PARMS    *ip,
                             INCHI_IOSTREAM *prb_file,
                             INCHI_IOSTREAM *inp_file )
{
    if (ip->bSaveAllGoodStructsAsProblem &&
         nRet != _IS_FATAL                &&
         nRet != _IS_ERROR                &&
         prb_file                         &&
         prb_file->f &&
         0L <= sd->fPtrStart              &&
         sd->fPtrStart < sd->fPtrEnd)
    {
        MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, 0 ); /* djb-rwth: addressing coverity ID #499510 -- return values handled properly */
    }

    return;
}


/****************************************************************************
 Generate InChI for the whole (multi-component) structure
****************************************************************************/
int CreateOneStructureINChI( CANON_GLOBALS          *pCG,
                             INCHI_CLOCK            *ic,
                             STRUCT_DATA            *sd,
                             INPUT_PARMS            *ip,
                             char                   *szTitle,
                             PINChI2                *pINChI2[INCHI_NUM],
                             PINChI_Aux2            *pINChI_Aux2[INCHI_NUM],
                             int                    iINChI,
                             INCHI_IOSTREAM         *inp_file,
                             INCHI_IOSTREAM         *log_file,
                             INCHI_IOSTREAM         *out_file,
                             INCHI_IOSTREAM         *prb_file,
                             ORIG_ATOM_DATA         *orig_inp_data,
                             ORIG_ATOM_DATA         *prep_inp_data,
                             COMP_ATOM_DATA         composite_norm_data2[][TAUT_NUM + 1],
                             long                   num_inp,
                             INCHI_IOS_STRING		*strbuf,
                             NORM_CANON_FLAGS       *pncFlags )
{
    int i, j, k, nRet = 0, n = 0l;
#if defined (TARGET_EXE_STANDALONE) && defined(_WIN32)
    int err_display;
#endif

#ifdef GHI100_FIX
#if ((SPRINTF_FLAG != 1) && (SPRINTF_FLAG != 2))
    setlocale(LC_ALL, "en-US"); /* djb-rwth: setting all locales to "en-US" */
#endif
#endif

    PINChI2     *pINChI = NULL;
    PINChI_Aux2 *pINChI_Aux = NULL;

    INP_ATOM_DATA InpCurAtData;
    INP_ATOM_DATA *inp_cur_data;

    INP_ATOM_DATA InpNormAtData, InpNormTautData;
    INP_ATOM_DATA *inp_norm_data[TAUT_NUM]; /*  = { &InpNormAtData, &InpNormTautData }; */
    ORIG_ATOM_DATA *cur_prep_inp_data = prep_inp_data + iINChI;
    inchiTime      ulTStart;

    /* Always create info data structures (but do not display them always )
    #ifndef COMPILE_ANSI_ONLY
    */
    int            bShowStructure = 0;
    int            bStructurePreprocessed = 0; /* All changes except disconnection */
    int            bStructureDisconnected = 0;
    int            bAlsoOutputReconnected = 0, bINCHI_LIB_Flag = 0;
    COMP_ATOM_DATA *composite_norm_data = composite_norm_data2[iINChI];
    INP_ATOM_DATA2 *all_inp_norm_data = NULL;
    /*#endif*/

    /*        Order of actions:

        if ( orig_inp_data is NOT empty AND
             prep_inp_data[0] IS empty ) then do
             in PreprocessOneStructure()        :

            1. copy orig_inp_data --> prep_inp_data[0]
            2. fix odd things in prep_inp_data[0]
            3. if( orig_inp_data->bDisconnectSalts ) then
                  -- disconnect salts in prep_inp_data[0]
            4. move protons to neutralize charges on heteroatoms
            5. if( orig_inp_data->bDisconnectCoord ) then
                  -- copy prep_inp_data[0] --> prep_inp_data[1]
                  -- disconnect metals in prep_inp_data[0]

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

    ip->msec_LeftTime = ip->msec_MaxTime; /* start timeout countdown for each component */

    inp_cur_data = &InpCurAtData;
    inp_norm_data[TAUT_NON] = &InpNormAtData;
    inp_norm_data[TAUT_YES] = &InpNormTautData;

    memset( inp_cur_data, 0, sizeof( *inp_cur_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( inp_norm_data[TAUT_NON], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( inp_norm_data[TAUT_YES], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    {
        /*#ifndef COMPILE_ANSI_ONLY*/
        memset( composite_norm_data + TAUT_NON, 0, sizeof( composite_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        memset( composite_norm_data + TAUT_YES, 0, sizeof( composite_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        memset( composite_norm_data + TAUT_INI, 0, sizeof( composite_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    }   /*#endif*/

    if (ip->bAllowEmptyStructure && !orig_inp_data->at && !orig_inp_data->num_inp_atoms)
    {
        ;
    }
    else if (!orig_inp_data->at || !orig_inp_data->num_inp_atoms)
    {
        return 0; /* nothing to do */
    }
    if (iINChI == 1 && orig_inp_data->bDisconnectCoord <= 0)
    {
        return 0;
    }

    /* m = iINChI; */ /* orig_inp_data index */
    if (iINChI != INCHI_BAS && iINChI != INCHI_REC)
    {
        AddErrorMessage( sd->pStrErrStruct, "Fatal undetermined program error" );
        sd->nStructReadError = 97;
        nRet = sd->nErrorType = _IS_FATAL;
        inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
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

    if (( !prep_inp_data->at || !prep_inp_data->num_inp_atoms ) &&
         orig_inp_data->num_inp_atoms > 0)
    {
        /* The structure has not been preprocessed */
        if (ip->msec_MaxTime)
        {
            InchiTimeGet( &ulTStart );
        }

        PreprocessOneStructure( ic, sd, ip, orig_inp_data, prep_inp_data );

        pncFlags->bTautFlags[iINChI][TAUT_YES] =
                pncFlags->bTautFlags[iINChI][TAUT_NON] =
                    sd->bTautFlags[INCHI_BAS] | ip->bTautFlags;

        pncFlags->bTautFlagsDone[iINChI][TAUT_YES] =
            pncFlags->bTautFlagsDone[iINChI][TAUT_NON] =
            sd->bTautFlagsDone[INCHI_BAS] | ip->bTautFlagsDone;

        {
            /*#ifndef COMPILE_ANSI_ONLY*/
            /* in this location the call happens once for each input structure, before preprocessing */
            bStructurePreprocessed = ( 0 != ( sd->bTautFlagsDone[INCHI_BAS] & (
                TG_FLAG_MOVE_HPLUS2NEUTR_DONE |
                TG_FLAG_DISCONNECT_SALTS_DONE |
                TG_FLAG_MOVE_POS_CHARGES_DONE |
                TG_FLAG_FIX_ODD_THINGS_DONE ) ) );

            bStructureDisconnected = ( 0 != ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) );

            bShowStructure = ( bStructurePreprocessed ||
                               bStructureDisconnected ||
                               prep_inp_data[0].num_components > 1 );

            /* sd->bTautFlags[] contains output flags
               ip->bTautFlags   contains input flags
            */
            bAlsoOutputReconnected = ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
                ( ip->bTautFlags               & TG_FLAG_RECONNECT_COORD );
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
            if (prep_inp_data[INCHI_BAS].num_components > 1)
            {
                bINCHI_LIB_Flag |= COMP_ORIG_0_MAIN;
            }
            else
            {
                /* Original -> Main, component #1, Original */
                if (prep_inp_data[INCHI_BAS].num_components == 1 && bStructurePreprocessed)
                {
                    bINCHI_LIB_Flag |= COMP_ORIG_1_MAIN;
                    /* preprocessed will be added when output canonicalization results */
                }
            }

            if (bAlsoOutputReconnected)
            {
                /* Original -> Reconnected, component #0, Original */
                if (prep_inp_data[INCHI_REC].num_components > 1)
                {
                    bINCHI_LIB_Flag |= COMP_ORIG_0_RECN;
                }
                else if (prep_inp_data[INCHI_BAS].num_components == 1 && bStructurePreprocessed)
                {
                    /* Original -> Reconnected, component #1, Original */
                    bINCHI_LIB_Flag |= COMP_ORIG_1_RECN;
                    /* preprocessed will be added when output canonicalization results */
                }
            }
            if (ip->msec_MaxTime)
            {
                ip->msec_LeftTime -= InchiTimeElapsed( ic, &ulTStart );
            }

            /* display the ORIGINAL, UN-PREPROCESSED structure */

            if (ip->bDisplay)
            {
                if (DisplayTheWholeStructure( pCG, ic, sd, ip, szTitle,
                    inp_file, log_file, orig_inp_data, num_inp,
                    -1, bShowStructure, bINCHI_LIB_Flag ))
                {
                    inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
                    goto exit_function;
                }
            }
        } /*#endif */

        switch (sd->nErrorType)
        {
            case _IS_ERROR:
            case _IS_FATAL:
                /* error message */
                nRet = TreatErrorsInReadTheStructure( sd, ip,
                                                      LOG_MASK_ALL,
                                                      inp_file, log_file, out_file, prb_file,
                                                      prep_inp_data, &num_inp );
                goto exit_cycle;
        }
    }
    /* tranfer flags from INChI_Aux to sd */


    {
        /*#ifndef COMPILE_ANSI_ONLY */ /* { */

        /******************************************/
        /*      Displaying the structures         */
        /*          Only under WIN32              */
        /******************************************/
        if ( /* ip->bDisplayCompositeResults && !sd->bUserQuitComponentDisplay && */
             prep_inp_data[iINChI].num_components > 1
           )
        {
            all_inp_norm_data = (INP_ATOM_DATA2 *) inchi_calloc( prep_inp_data[iINChI].num_components, sizeof( all_inp_norm_data[0] ) );
        }

        /* Display the input structure AFTER PREPROCESSING */
        switch (iINChI)
        {
            case INCHI_BAS:
                /*------------ Possibly disconnected structure -------------------*/
                bStructurePreprocessed = 0 != ( sd->bTautFlagsDone[iINChI] & (
                    TG_FLAG_MOVE_HPLUS2NEUTR_DONE |
                    TG_FLAG_DISCONNECT_SALTS_DONE |
                    TG_FLAG_MOVE_POS_CHARGES_DONE |
                    TG_FLAG_MOVE_CHARGE_COORD_DONE |
                    TG_FLAG_DISCONNECT_COORD_DONE |
                    TG_FLAG_FIX_ODD_THINGS_DONE ) );
                bINCHI_LIB_Flag = 0;
                /* Preprocessed/Main -> Main, component #0, Preprocessed */
                if (prep_inp_data[iINChI].num_components > 1 &&
                     bStructurePreprocessed)
                {
                    bINCHI_LIB_Flag |= COMP_PREP_0_MAIN;
                }
                bShowStructure = ( bStructurePreprocessed &&
                                   prep_inp_data[iINChI].num_components > 1 );
                break;

            case INCHI_REC:
                /*------------ Reconnected structure ------------------------------*/
                bAlsoOutputReconnected =
                    ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
                    ( ip->bTautFlags               & TG_FLAG_RECONNECT_COORD );

                if (!bAlsoOutputReconnected)
                {
                    break;
                }

                bStructurePreprocessed = 0 != ( sd->bTautFlagsDone[iINChI] & (
                    TG_FLAG_MOVE_HPLUS2NEUTR_DONE |
                    TG_FLAG_DISCONNECT_SALTS_DONE |
                    TG_FLAG_MOVE_POS_CHARGES_DONE |
                    TG_FLAG_FIX_ODD_THINGS_DONE ) );
                bINCHI_LIB_Flag = 0;
                /* Preprocessed/Reconnected -> Reconnected, component #0, Preprocessed */
                if (prep_inp_data[iINChI].num_components > 1 && bStructurePreprocessed)
                {
                    bINCHI_LIB_Flag |= COMP_PREP_0_RECN;
                }
                bShowStructure = ( bStructurePreprocessed &&
                                   prep_inp_data[iINChI].num_components > 1 );
                break;

            default:
                bShowStructure = 0;
        }


        if (ip->bDisplay && prep_inp_data[iINChI].num_inp_atoms > 0)
        {
            if (DisplayTheWholeStructure( pCG, ic, sd, ip, szTitle,
                inp_file, log_file,
                prep_inp_data + iINChI,
                num_inp,
                iINChI,
                bShowStructure,
                bINCHI_LIB_Flag ))
            {
                inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
                goto exit_function;
            }
        }
    } /* #endif */ /* } ifndef COMPILE_ANSI_ONLY */


    /* allocate pINChI[iINChI] and pINChI_Aux2[iINChI] -- arrays of pointers to INChI and INChI_Aux */
    /* assign values to sd->num_components[]                                                  */
    
    /* djb-rwth: MYREALLOC2 has been replaced and the whole block rewritten to address memory leaks and reading from freed memory locations */
    do {
        if( (sd->num_components[iINChI]) <= ((long long)cur_prep_inp_data->num_components) ) {
            PINChI2* newPTR1 = (PINChI2*)inchi_calloc( ((long long)cur_prep_inp_data->num_components)+1, sizeof(PINChI2) );
            PINChI_Aux2* newPTR2 = (PINChI_Aux2*)inchi_calloc( ((long long)cur_prep_inp_data->num_components)+1, sizeof(PINChI_Aux2) );
            if ( newPTR1 && newPTR2 ) {
                if ( (pINChI2[iINChI]) && (sd->num_components[iINChI]) > 0 )
                    memcpy(newPTR1, pINChI2[iINChI], (sd->num_components[iINChI]) * sizeof(PINChI2));
                if ((pINChI_Aux2[iINChI]) && (sd->num_components[iINChI]) > 0)
                    memcpy(newPTR2, pINChI_Aux2[iINChI], (sd->num_components[iINChI]) * sizeof(PINChI_Aux2));
                if (pINChI2[iINChI])
                    inchi_free(pINChI2[iINChI]);
                if (pINChI_Aux2[iINChI])
                    inchi_free(pINChI_Aux2[iINChI]);
                pINChI2[iINChI] = newPTR1;
                pINChI_Aux2[iINChI] = newPTR2;
                sd->num_components[iINChI] = cur_prep_inp_data->num_components;
                k = 0;
            }
            else {
                inchi_free(newPTR1);
                inchi_free(newPTR2);
                k = 1;
            }
        }
        else { k = 0; }
    } while (0);

    if (k)
    {
        AddErrorMessage( sd->pStrErrStruct, "Cannot allocate output data. Terminating" );
        sd->nStructReadError = 99;
        sd->nErrorType = _IS_FATAL;
        inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
        goto exit_function;
    }

    pINChI = pINChI2[iINChI];
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

    for (i = 0, nRet = 0;
            !sd->bUserQuitComponent && i < cur_prep_inp_data->num_components;
                i++)
    {
        if (ip->msec_MaxTime)
        {
            InchiTimeGet( &ulTStart );
        }

#ifndef TARGET_LIB_FOR_WINCHI  /* { */
#if ( bREUSE_INCHI == 1 )

        if ((iINChI == INCHI_REC &&
             /*( !ip->bDisplay &&
               !ip->bDisplayCompositeResults && */
               !( ip->bCompareComponents & CMP_COMPONENTS )) ||
               sd->bUserQuitComponentDisplay) /* djb-rwth: addressing LLVM warning */
        {
            /* Reconnected structure (06-20-2005: added "&& !ip->bDisplayCompositeResults" to display composite structure) */
            int m = iINChI - 1;

            /* Find whether we have already calculated this INChI in basic (disconnected) layer */
            for (j = n = 0; j < prep_inp_data[m].num_components; j++)
            {
                if (i + 1 == prep_inp_data[m].nOldCompNumber[j] &&
                    ( pINChI2[m][j][TAUT_NON] || pINChI2[m][j][TAUT_YES] ))
                {
                    /* Yes, we have already done this */
                    if (!n++)
                    {
                        memcpy(pINChI + i, pINChI2[m] + j, sizeof(pINChI[0]));
                        memcpy(pINChI_Aux + i, pINChI_Aux2[m] + j, sizeof(pINChI_Aux[0]));
                        for (k = 0; k < TAUT_NUM; k++)
                        {
                            if (pINChI[i][k])
                            {
                                pINChI[i][k]->nRefCount++;
                                if (pINChI[i][k]->nNumberOfAtoms > 0)
                                {
                                    switch (k)
                                    {
                                        case TAUT_NON:
                                            sd->num_non_taut[iINChI] ++;
                                            break;
                                        case TAUT_YES:
                                            if (pINChI[i][k]->lenTautomer > 0)
                                            {
                                                sd->num_taut[iINChI] ++;
                                            }
                                            else
                                                if (!pINChI[i][TAUT_NON] ||
                                                     !pINChI[i][TAUT_NON]->nNumberOfAtoms)
                                                {
                                                    sd->num_non_taut[iINChI] ++;
                                                }
                                            break;
                                    }
                                }
                            }
                            if (pINChI_Aux[i][k])
                            {
                                pINChI_Aux[i][k]->nRefCount++;
                            }
                        }
                    }
                }
            }

            if (n == 1)
            {
                continue;
            }
            if (n > 1)
            {
                /* ith component is equivalent to more than one another component */
                AddErrorMessage( sd->pStrErrStruct, "Cannot distinguish components" );
                sd->nStructReadError = 99;
                sd->nErrorType = _IS_ERROR;
                inchi_free(all_inp_norm_data); /* djb-rwth: avoiding memory leak */
                goto exit_function;
            }
        }

#endif
#endif /* } TARGET_LIB_FOR_WINCHI */

        /*****************************************************/
        /*  a) Allocate memory and extract current component */
        /*****************************************************/

        nRet = GetOneComponent( ic, sd, ip,
                                log_file, out_file,
                                inp_cur_data, cur_prep_inp_data,
                                i, num_inp );

        if (ip->msec_MaxTime)
        {
            ip->msec_LeftTime -= InchiTimeElapsed( ic, &ulTStart );
        }

        switch (nRet)
        {
            case _IS_ERROR:
            case _IS_FATAL:
                goto exit_cycle;
        }

#if !defined(TARGET_API_LIB) && !defined(COMPILE_ANSI_ONLY)
        /*  console request: Display the component? */
        if (ip->bDisplay && inp_file->f != stdin)
        {
            if (user_quit( ic, "Enter=Display Component, Esc=Stop ?", ip->ulDisplTime ))
            {
                sd->bUserQuitComponent = 1;
                break;
            }
        }
#endif

        /*#ifndef COMPILE_ANSI_ONLY  
        { */

        /*  b) Display the extracted original component structure */
        if (ip->bDisplay && inp_cur_data->at && !sd->bUserQuitComponentDisplay)
        {
            if (cur_prep_inp_data->num_components == 1)
            {
                sprintf(szTitle, "%sInput Structure #%ld.%s%s%s%s%s",
                    bStructurePreprocessed ? "Preprocessed " : "",
                    num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
            }
            else
            {
                sprintf(szTitle, "Component #%d of %d, Input Structure #%ld.%s%s%s%s%s",
                    i + 1, cur_prep_inp_data->num_components,
                    num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
            }

#if defined (TARGET_EXE_STANDALONE) && defined(_WIN32)
            err_display = DisplayStructure( pCG,
                                    inp_cur_data->at,
                                    inp_cur_data->num_at,
                                    NULL, /* OAD_Polymer *polymer, */
                                    0,
                                    1,
                                    0,
                                    NULL,
                                    1           /*isotopic*/,
                                    0           /*taut*/,
                                    NULL,
                                    NULL,
                                    ip->bAbcNumbers,
                                    &ip->dp,
                                    ip->nMode,
                                    szTitle );

            sd->bUserQuitComponentDisplay = (err_display == ESC_KEY );

            if (!err_display)
            {
                inchi_fprintf( stderr, "Cannot display the structure\n" );
            }
#endif
#ifdef TARGET_LIB_FOR_WINCHI
            if (DRAWDATA && DRAWDATA_EXISTS)
            {
                struct DrawData vDrawData;
                int    nType = COMPONENT_ORIGINAL;
                vDrawData.pWindowData = CreateWinData_( pCG,
                                                        inp_cur_data->at,
                                                        inp_cur_data->num_at,
                                                        0,
                                                        1 /* bAdd_DT_to_num_H */,
                                                        0,
                                                        NULL,
                                                        1 /* display isotopic if present */,
                                                        0,
                                                        NULL,
                                                        NULL,
                                                        ip->bAbcNumbers,
                                                        &ip->dp,
                                                        ip->nMode );
                if (vDrawData.pWindowData != NULL)
                {
                    if (DRAWDATA_EXISTS( i + 1, nType, iINChI ))
                    {
                        /* i = component number */
                        nType = COMPONENT_ORIGINAL_PREPROCESSED;
                    }
                    vDrawData.nComponent = i + 1;
                    vDrawData.nType = nType;
                    vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                    vDrawData.szTitle = inchi__strdup( szTitle );
                    vDrawData.pWindowData->szTitle = inchi__strdup( szTitle );
                    DRAWDATA( &vDrawData );
                }
            }
#endif
        }
        /*#endif */  /* } COMPILE_ANSI_ONLY */


        /*******************************************************************************/
        /*                                                                             */
        /*  N O R M A L I Z A T I O N    a n d     C A N O N I C A L I Z A T I O N     */
        /*                                                                             */
        /*         (both tautomeric and non-tautomeric if requested)                   */
        /*                                                                             */
        /*******************************************************************************/
        /*  c) Create the component's INChI ( copies ip->bTautFlags into sd->bTautFlags)*/
        /*******************************************************************************/

        nRet = CreateOneComponentINChI( pCG, ic, sd, ip,
                                        inp_cur_data, orig_inp_data,
                                        pINChI/*2[iINChI]*/,
                                        pINChI_Aux/*2[iINChI]*/,
                                        iINChI, i, num_inp,
                                        inp_norm_data, pncFlags, log_file );



        /*  d) Display one component structure and/or INChI results only if there was no error */

        /* #ifndef COMPILE_ANSI_ONLY */ /* { */
        if (!nRet)
        {
            /*  output one component INChI to the stdout if requested */
            /*
            if ( ip->bDisplayEachComponentINChI ) {
                int cur_num_non_taut = (pINChI[i][TAUT_NON] && pINChI[i][TAUT_NON]->nNumberOfAtoms>0);
                int cur_num_taut     = (pINChI[i][TAUT_YES] && pINChI[i][TAUT_YES]->nNumberOfAtoms>0);
                if ( ip->bDisplayEachComponentINChI && cur_num_non_taut + cur_num_taut ) {
                    SortAndPrintINChI( pCG, stdout, pStr, nStrLen, NULL,
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

            int bIsotopic, bTautomeric, bDisplayTaut, bHasIsotopicLayer, bFixedBondsTaut, m_max, m, nNumDisplayedFixedBondTaut = 0; /* djb-rwth: ignoring LLVM warning: variable used */

            for ( j = 0;
                  ip->bDisplay && !sd->bUserQuitComponentDisplay && j < TAUT_NUM;
                  j++)
            {
                if (inp_norm_data[j]->bExists && !inp_norm_data[j]->bDeleted)
                {
                    bTautomeric = ( pINChI[i][j]->lenTautomer > 0 );
                     /* same as (inp_norm_data[j]->bTautomeric > 0) */

                    /* If requested tautomeric and no tautmerism found then do not say mobile or fixed H. 2004-10-27 */
                    bDisplayTaut = ( !( ip->nMode & REQ_MODE_BASIC ) && !bTautomeric ) ? -1 : bTautomeric;
                    bHasIsotopicLayer = ( inp_norm_data[j]->bHasIsotopicLayer > 0 );

                    for (k = 0; k <= bHasIsotopicLayer; k++)
                    {
                        bIsotopic = ( k > 0 );
                        m_max = inp_norm_data[j]->at_fixed_bonds && inp_norm_data[j]->bTautPreprocessed ? 1 : 0;
                        for (m = m_max; 0 <= m; m--)
                        {
                            bFixedBondsTaut = ( m > 0 );
                            nNumDisplayedFixedBondTaut += bFixedBondsTaut;
                                /* display only one time */

                            /*  Added number of components, added another format for a single component case - DCh */
                            if (cur_prep_inp_data->num_components > 1)
                            {
                                sprintf(szTitle, "%s Component #%d of %d, Structure #%ld%s%s.%s%s%s%s%s",
                                    bFixedBondsTaut ? "Preprocessed" : "Result for",
                                    i + 1, cur_prep_inp_data->num_components, num_inp,
                                    bDisplayTaut == 1 ? ", mobile H" : bDisplayTaut == 0 ? ", fixed H" : "",
                                    bIsotopic ? ", isotopic" : "",
                                    SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
                            }
                            else
                            {
                                sprintf(szTitle, "%s Structure #%ld%s%s.%s%s%s%s%s",
                                    bFixedBondsTaut ? "Preprocessed" : "Result for",
                                    num_inp,
                                    bDisplayTaut == 1 ? ", mobile H" : bDisplayTaut == 0 ? ", fixed H" : "",
                                    bIsotopic ? ", isotopic" : "",
                                    SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
                            }

#if defined (TARGET_EXE_STANDALONE) && defined(_WIN32)
                            if (bFixedBondsTaut && nNumDisplayedFixedBondTaut != 1)
                            {
                                continue;
                            }
                            if (ip->bDisplay)
                            {
                                if (bFixedBondsTaut)
                                {
                                    err_display = DisplayStructure( pCG,
                                                            inp_norm_data[j]->at_fixed_bonds,
                                                            inp_norm_data[j]->num_at,
                                                            NULL, /* OAD_Polymer *polymer, */
                                                            inp_norm_data[j]->num_removed_H,
                                                            0 /*bAdd_DT_to_num_H*/,
                                                            inp_norm_data[j]->nNumRemovedProtons,
                                                            inp_norm_data[j]->nNumRemovedProtonsIsotopic,
                                                            bHasIsotopicLayer,
                                                            j,
                                                            NULL,
                                                            NULL,
                                                            ip->bAbcNumbers,
                                                            &ip->dp,
                                                            ip->nMode,
                                                            szTitle );
                                }
                                else
                                {
                                    err_display = DisplayStructure( pCG,
                                                            inp_norm_data[j]->at,
                                                            inp_norm_data[j]->num_at,
                                                            NULL, /* OAD_Polymer *polymer, */
                                                            0,
                                                            0 /*bAdd_DT_to_num_H*/,
                                                            0,
                                                            NULL,
                                                            k,
                                                            j,
                                                            pINChI[i],
                                                            pINChI_Aux[i],
                                                            ip->bAbcNumbers,
                                                            &ip->dp,
                                                            ip->nMode, szTitle );
                                }
                                if ((sd->bUserQuitComponentDisplay = (err_display == ESC_KEY ))) /* djb-rwth: addressing LLVM warning */
                                {
                                    break;
                                }
                            }
#endif
#ifdef TARGET_LIB_FOR_WINCHI

                            if (DRAWDATA && !bFixedBondsTaut)
                            {
                                struct DrawData vDrawData;
                                vDrawData.pWindowData =
                                    CreateWinData_( pCG,
                                                    inp_norm_data[j]->at,
                                                    inp_norm_data[j]->num_at,
                                                    0,
                                                    0 /* bAdd_DT_to_num_H */,
                                                    0,
                                                    NULL,
                                                    k,
                                                    j,
                                                    pINChI[i],
                                                    pINChI_Aux[i],
                                                    ip->bAbcNumbers,
                                                    &ip->dp,
                                                    ip->nMode );

                                if (vDrawData.pWindowData != NULL)
                                {
                                    int nType;
                                    vDrawData.nComponent = i + 1;
                                    if (bTautomeric == 0)
                                        nType = ( bIsotopic == 0 ) ? COMPONENT_BN
                                        : COMPONENT_BI;
                                    else
                                        nType = ( bIsotopic == 0 ) ? COMPONENT_TN
                                        : COMPONENT_TI;
                                    vDrawData.nType = nType;

                                    vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                                    vDrawData.szTitle = inchi__strdup( szTitle );
                                    vDrawData.pWindowData->szTitle = inchi__strdup( szTitle );
                                    DRAWDATA( &vDrawData );
                                }
                            }
                            else if (DRAWDATA && bFixedBondsTaut)
                            {
                                struct DrawData vDrawData;
                                if (( ip->bCompareComponents & CMP_COMPONENTS ) &&
                                     !( ip->bCompareComponents & CMP_COMPONENTS_NONTAUT ) &&
                                     !bIsotopic == !inp_norm_data[j]->bHasIsotopicLayer)
                                {

                                    vDrawData.pWindowData =
                                        CreateWinData_( pCG,
                                                        inp_norm_data[j]->at_fixed_bonds,
                                                        inp_norm_data[j]->num_at,
                                                        inp_norm_data[j]->num_removed_H,
                                                        0 /* bAdd_DT_to_num_H */,
                                                        inp_norm_data[j]->nNumRemovedProtons,
                                                        inp_norm_data[j]->nNumRemovedProtonsIsotopic,
                                                        k,
                                                        j,
                                                        NULL,
                                                        NULL,
                                                        ip->bAbcNumbers,
                                                        &ip->dp,
                                                        ip->nMode );
                                }
                                else
                                {
                                    continue;
                                }
                                if (vDrawData.pWindowData != NULL)
                                {
                                    vDrawData.nComponent = i + 1;
                                    vDrawData.nType = COMPONENT_ORIGINAL_PREPROCESSED;
                                    vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                                    vDrawData.szTitle = inchi__strdup( szTitle );
                                    vDrawData.pWindowData->szTitle = inchi__strdup( szTitle );
                                    DRAWDATA( &vDrawData );
                                }
                            }
#endif
                        }
                    }
                }
            }

            /* Save normalized components for composite display */
            if ( /*ip->bDisplayCompositeResults && */
                 all_inp_norm_data
                )
            {
                for (j = 0; j < TAUT_NUM; j++)
                {
                    if (inp_norm_data[j]->bExists)
                    {
                        all_inp_norm_data[i][j] = *inp_norm_data[j];
                        memset( inp_norm_data[j], 0, sizeof( *inp_norm_data[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
                    }
                }
            }
        }

        /* #endif */ /* } COMPILE_ANSI_ONLY */


        if (nRet)
        {
            nRet = TreatErrorsInCreateOneComponentINChI( sd, ip,
                                                         cur_prep_inp_data,
                                                         i, num_inp, inp_file,
                                                         log_file, out_file, prb_file );
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
    switch (nRet)
    {
        case _IS_FATAL:
        case _IS_ERROR:
            break;
        default:

        /* #ifndef COMPILE_ANSI_ONLY *//* { */
            /* composite results picture(s) */
            if (all_inp_norm_data)
            {
                CreateCompositeNormAtom( composite_norm_data,
                                         all_inp_norm_data,
                                         prep_inp_data[iINChI].num_components );
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
            /* #endif */ /* } COMPILE_ANSI_ONLY */

            break;
    }

    /*#ifndef COMPILE_ANSI_ONLY*/ /* { */

        /* avoid memory leaks in case of error */
    if (all_inp_norm_data)
    {
        for (i = 0; i < prep_inp_data[iINChI].num_components; i++)
        {
            for (k = 0; k < TAUT_NUM; k++)
            {
                FreeInpAtomData( &all_inp_norm_data[i][k] );
            }
        }
        inchi_free( all_inp_norm_data );
        all_inp_norm_data = NULL;
    }
/*#endif */ /* } COMPILE_ANSI_ONLY */

    FreeInpAtomData( inp_cur_data );
    for (i = 0; i < TAUT_NUM; i++)
    {
        FreeInpAtomData( inp_norm_data[i] );
    }


exit_function:

    return nRet;
}


/****************************************************************************
 Generate InChI for one connected component
 (of possibly multi-component structure)
****************************************************************************/
int CreateOneComponentINChI( CANON_GLOBALS      *pCG,
                             INCHI_CLOCK        *ic,
                             STRUCT_DATA        *sd,
                             INPUT_PARMS        *ip,
                             INP_ATOM_DATA      *inp_cur_data,
                             ORIG_ATOM_DATA     *orig_inp_data,
                             PINChI2            *pINChI,
                             PINChI_Aux2        *pINChI_Aux,
                             int                iINChI,
                             int                i,
                             long               num_inp,
                             INP_ATOM_DATA      **inp_norm_data,
                             NORM_CANON_FLAGS   *pncFlags,
                             INCHI_IOSTREAM     *log_file )
{
    inchiTime     ulTStart, ulTEnd, *pulTEnd = NULL;
    int           k, num_at, ret = 0;
    int           bOrigCoord;
    INCHI_MODE     bTautFlags = ip->bTautFlags;
    INCHI_MODE     bTautFlagsDone = ( ip->bTautFlagsDone | sd->bTautFlagsDone[INCHI_BAS] );
    INChI       *cur_INChI[TAUT_NUM];
    INChI_Aux   *cur_INChI_Aux[TAUT_NUM];
    long          lElapsedTime;

    int nAllocMode = 0;  /* moved from below 2024-09-01 DT */

#ifdef GHI100_FIX
#if ((SPRINTF_FLAG != 1) && (SPRINTF_FLAG != 2))
    setlocale(LC_ALL, "en-US"); /* djb-rwth: setting all locales to "en-US" */
#endif
#endif

    InchiTimeGet( &ulTStart );
    bOrigCoord =
        !( ip->bINChIOutputOptions & ( INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO ) );

    for (k = 0; k < TAUT_NUM; k++)
    {
        cur_INChI[k] = NULL;
        cur_INChI_Aux[k] = NULL;
    }

    /*  Allocate memory for non-tautomeric (k=0) and tautomeric (k=1) results */
    for (k = 0; k < TAUT_NUM; k++)
    {
        /* djb-rwth: introducing variables for correct nAllocMode expression */
        int nAM1 = 0, nAM2 = 0;

        if (k == TAUT_YES)
            nAM1 = REQ_MODE_TAUT;

        if (bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))
            nAM2 = ip->nMode & REQ_MODE_ISO;

        nAllocMode = nAM1 | nAM2; /* djb-rwth: original sequence of bit-wise operations had to be rewritten */


        if ((k == TAUT_NON && ( ip->nMode & REQ_MODE_BASIC )) ||
             (k == TAUT_YES && ( ip->nMode & REQ_MODE_TAUT ))) /* djb-rwth: addressing LLVM warning */
        {
            /*  alloc INChI and INChI_Aux */
            cur_INChI[k] = Alloc_INChI( inp_cur_data->at,
                                            inp_cur_data->num_at,
                                            &inp_cur_data->num_bonds,
                                            &inp_cur_data->num_isotopic,
                                            nAllocMode );

            cur_INChI_Aux[k] = Alloc_INChI_Aux( inp_cur_data->num_at,
                                                inp_cur_data->num_isotopic,
                                                nAllocMode,
                                                bOrigCoord );
            if (cur_INChI_Aux[k])
            {
                cur_INChI_Aux[k]->bIsIsotopic = inp_cur_data->num_isotopic;
            }
            /*  alloc memory for the output structure: non-tautomeric and tautomeric (for displaying) */

            CreateInpAtomData( inp_norm_data[k], inp_cur_data->num_at, k );
        }
        else
        {
            FreeInpAtomData( inp_norm_data[k] );
        }
    }

    lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
    if (ip->msec_MaxTime)
    {
        ip->msec_LeftTime -= lElapsedTime;
    }
    sd->ulStructTime += lElapsedTime;

    /*^^^#if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined( TARGET_API_LIB ) ) */
    #if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined( TARGET_API_LIB ) && !defined(TARGET_EXE_STANDALONE) )
    #endif

    /******************************************************
     *
     *  Get one component canonical numberings, etc.
     *
     ******************************************************/

    /* Create_INChI() return value:
     * num_at <= 0: error code
     * num_at >  0: number of atoms (excluding terminal hydrogen atoms)
     * inp_norm_data[0] => non-tautomeric, inp_norm_data[1] => tautomeric */

    InchiTimeGet( &ulTStart );

    if (ip->msec_MaxTime)
    {
        ulTEnd = ulTStart;
        pulTEnd = &ulTEnd;
        if (ip->msec_LeftTime > 0)
        {
            InchiTimeAddMsec( ic, pulTEnd, ip->msec_LeftTime );
        }
    }

    num_at = Create_INChI( pCG, ic, ip,
                           cur_INChI, cur_INChI_Aux,
                           orig_inp_data/* not used */,
                           inp_cur_data->at, inp_norm_data, inp_cur_data->num_at,
                           ip->nMode,
                           &bTautFlags, &bTautFlagsDone,
                           pulTEnd, NULL, sd->pStrErrStruct );

    SetConnectedComponentNumber( inp_cur_data->at, inp_cur_data->num_at, i + 1 );
                        /*  NB: normalization alters structure component number */

    for (k = 0; k < TAUT_NUM; k++)
    {
        if (cur_INChI_Aux[k] && cur_INChI_Aux[k]->nNumberOfAtoms > 0)
        {
            pncFlags->bNormalizationFlags[iINChI][k] |=
                cur_INChI_Aux[k]->bNormalizationFlags;
            pncFlags->bTautFlags[iINChI][k] |=
                cur_INChI_Aux[k]->bTautFlags;
            pncFlags->bTautFlagsDone[iINChI][k] |=
                cur_INChI_Aux[k]->bTautFlagsDone;
            pncFlags->nCanonFlags[iINChI][k] |=
                cur_INChI_Aux[k]->nCanonFlags;
        }
    }

    /*  Detect errors */
    if (num_at < 0)
    {
        sd->nErrorCode = num_at;
    }
    else if (num_at == 0)
    {
        sd->nErrorCode = -1;
    }
    else if (cur_INChI[TAUT_NON] && cur_INChI[TAUT_NON]->nErrorCode)
    {   /*  non-tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_NON]->nErrorCode;
    }
    else if (cur_INChI[TAUT_YES] && cur_INChI[TAUT_YES]->nErrorCode)
    {   /*  tautomeric error */
        sd->nErrorCode = cur_INChI[TAUT_YES]->nErrorCode;
    }


#if ( bRELEASE_VERSION == 0 )
    if (cur_INChI[TAUT_NON]) sd->bExtract |= cur_INChI[TAUT_NON]->bExtract;
    if (cur_INChI[TAUT_YES]) sd->bExtract |= cur_INChI[TAUT_YES]->bExtract;
    if (( TG_FLAG_TEST_TAUT3_SALTS_DONE & bTautFlagsDone ))
    {
        sd->bExtract |= EXTR_TEST_TAUT3_SALTS_DONE;
    }
#endif

    /*  Detect and store stereo warnings */
    if (!sd->nErrorCode)
        GetProcessingWarningsOneComponentInChI( cur_INChI, inp_norm_data, sd, ip->bNoWarnings );

    lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
    if (ip->msec_MaxTime)
    {
        ip->msec_LeftTime -= lElapsedTime;
    }
    sd->ulStructTime += lElapsedTime;

#if !defined(TARGET_API_LIB) && !defined(COMPILE_ANSI_ONLY)
    /*  Display the results */
    if (ip->bDisplay)
    {
        eat_keyboard_input( );
    }
#endif

    /*  a) No matter what happened save the allocated INChI pointers */
    /*  save the INChI of the current component */

    InchiTimeGet( &ulTStart );
    for (k = 0; k < TAUT_NUM; k++)
    {
        pINChI[i][k] = cur_INChI[k];
        pINChI_Aux[i][k] = cur_INChI_Aux[k];
        cur_INChI[k] = NULL;
        cur_INChI_Aux[k] = NULL;
    }

    /*  b) Count one component structure and/or INChI results only
           if there was no error
           Set inp_norm_data[j]->num_removed_H = number of removed explicit H
    */

    if (!sd->nErrorCode)
    {
        /*  find where the current processed structure is located */
        int cur_is_in_non_taut = ( pINChI[i][TAUT_NON] &&
                                   pINChI[i][TAUT_NON]->nNumberOfAtoms > 0 );
        int cur_is_in_taut = ( pINChI[i][TAUT_YES] &&
                                   pINChI[i][TAUT_YES]->nNumberOfAtoms > 0 );

        int cur_is_non_taut = (cur_is_in_non_taut && 0 == pINChI[i][TAUT_NON]->lenTautomer) ||
            (cur_is_in_taut && 0 == pINChI[i][TAUT_YES]->lenTautomer); /* djb-rwth: addressing LLVM warnings */
        int cur_is_taut = cur_is_in_taut && 0 < pINChI[i][TAUT_YES]->lenTautomer;

        /*
        sd->bTautFlags[iINChI]     |= bTautFlags;
        sd->bTautFlagsDone[iINChI] |= bTautFlagsDone;
        */

        if (cur_is_non_taut + cur_is_taut)
        {
            /*  count tautomeric and non-tautomeric components of the structures */
            int j1 = cur_is_in_non_taut ? TAUT_NON : TAUT_YES;
            int j2 = cur_is_in_taut ? TAUT_YES : TAUT_NON;
            int j;
            sd->num_non_taut[iINChI] += cur_is_non_taut;
            sd->num_taut[iINChI] += cur_is_taut;
            for (j = j1; j <= j2; j++)
            {
                int bIsotopic = ( pINChI[i][j]->nNumberOfIsotopicAtoms ||
                                  pINChI[i][j]->nNumberOfIsotopicTGroups ||
                                  (pINChI[i][j]->nPossibleLocationsOfIsotopicH && pINChI[i][j]->nPossibleLocationsOfIsotopicH[0] > 1) ); /* djb-rwth: addressing LLVM warning */
                if (pINChI_Aux[i][j] && (j == TAUT_YES)) /* djb-rwth: fixing a NULL pointer dereference */
                {
                    bIsotopic |= ( 0 < pINChI_Aux[i][j]->nNumRemovedIsotopicH[0] +
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[1] +
                                      pINChI_Aux[i][j]->nNumRemovedIsotopicH[2] );
                }

                inp_norm_data[j]->bExists = 1; /*  j=0: non-taut exists, j=1: taut exists */
                inp_norm_data[j]->bHasIsotopicLayer = bIsotopic;
                /*inp_norm_data[j]->num_removed_H = inp_norm_data[j]->num_at - num_at;*/
            }
        }
    }

    /* return (sd->nErrorCode==CT_OUT_OF_RAM || sd->nErrorCode==CT_USER_QUIT_ERR)? _IS_FATAL :
            sd->nErrorCode? _IS_ERROR : 0; */

    if (sd->nErrorCode == CT_OUT_OF_RAM || sd->nErrorCode == CT_USER_QUIT_ERR)
    {
        ret = _IS_FATAL;
    }
    else if (sd->nErrorCode)
    {
        ret = _IS_ERROR;
    }

    lElapsedTime = InchiTimeElapsed( ic, &ulTStart );
    if (ip->msec_MaxTime)
    {
        ip->msec_LeftTime -= lElapsedTime;
    }
    sd->ulStructTime += lElapsedTime;

    return ret;
}


/****************************************************************************
 Extended-functionality version of ProcessOneStructure

 able to handle both polymer related and unrelated pseudoelement atoms
 and to perform advanced polymer treatment (v. 1.06+)
 ****************************************************************************/
int ProcessOneStructureEx( struct tagINCHI_CLOCK    *ic,
                           struct tagCANON_GLOBALS  *CG,
                           STRUCT_DATA              *sd,
                           INPUT_PARMS              *ip,
                           char                     *szTitle,
                           PINChI2                  *pINChI2[INCHI_NUM],
                           PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                           INCHI_IOSTREAM           *inp_file,
                           INCHI_IOSTREAM           *log_file,
                           INCHI_IOSTREAM           *out_file,
                           INCHI_IOSTREAM           *prb_file,
                           ORIG_ATOM_DATA           *orig_inp_data,
                           ORIG_ATOM_DATA           *prep_inp_data,
                           long                     num_inp,
                           INCHI_IOS_STRING         *strbuf,
                           unsigned char            save_opt_bits )
{
    int ret = _IS_OKAY;
    char *sinchi_noedits=NULL, *saux_noedits=NULL;
    

    /* PREPROCESS */

#if (BUILD_WITH_ENG_OPTIONS==1)
#if ALLOW_SUBSTRUCTURE_FILTERING==1
    if (ip->bFilterSS)
    {
        int present, ok = 0;
        
        present = OrigAtData_CheckForSubstructure(orig_inp_data);

        if (ip->bFilterSS == 1 && present)			ok = 1;
        else if (ip->bFilterSS == -1 && !present)	ok = 1;

        if (!ok)
        {
            inchi_ios_eprint(log_file, "Warning (Skip record which does not pass substructure presence/absence filter) structure #%ld.%s%s%s%s\n",
                             num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));
            return _IS_SKIP;
        }
    }
#endif
#endif

    /*  Preprocess Polymer CRUs (collect frame shift info and edit the original input accordingly) */
    ret = PreprocessPolymerCRUData( ic, CG, sd, ip, szTitle,
                                    pINChI2, pINChI_Aux2,
                                    inp_file, log_file, out_file, prb_file,
                                    orig_inp_data, prep_inp_data,
                                    num_inp, strbuf, save_opt_bits,
                                    &sinchi_noedits, &saux_noedits); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    /* CALCULATE INCHI */
    
    /* Perform calculation as usual either for untouched (modes POLYMERS_LEGACY and POLYMERS_LEGACY_PLUS) 
       or just edited, probably (mode POLYMERS_MODERN) structure as passed in orig_inp_data
    */
    ret = ProcessOneStructureExCore( ic, CG, sd, ip,  szTitle, 
                                     pINChI2, pINChI_Aux2,
                                     inp_file, log_file, out_file, prb_file,
                                     orig_inp_data, prep_inp_data,
                                     num_inp, strbuf, save_opt_bits );

    if (ip->bINChIOutputOptions2 & INCHI_OUT_INCHI_GEN_ERROR)
    {
        if (ret == _IS_FATAL || ret == _IS_ERROR)
        {
            if (ip->bINChIOutputOptions & INCHI_OUT_STDINCHI)
            {
                inchi_ios_eprint(out_file, "InChI=1S//\n");
            }
            else
            {
                inchi_ios_eprint(out_file, "InChI=1//\n");
            }
        }
    }

    /* Post-process: add AuxInfo for a very unedited original structure */
    if (ret != _IS_FATAL && ret != _IS_ERROR)
    {
        if (!(ip->bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO) &&
            out_file->s.pStr && strstr(out_file->s.pStr, "AuxInfo=" ) &&
            saux_noedits && strstr(saux_noedits, "AuxInfo=") )
        {
            char *pp = out_file->s.pStr;
            if (saux_noedits[8])
            {
                remove_one_lf(pp);
                out_file->s.nUsedLength = strlen(out_file->s.pStr);
                inchi_ios_eprint(out_file, "/U/%-s\n", saux_noedits +8);
            }
        }
    }

    inchi_free(sinchi_noedits);
    inchi_free(saux_noedits);

    
#ifdef TARGET_LIB_FOR_WINCHI

    push_to_winchi_text_window(out_file);

    if (sd->pStrErrStruct && FWPUSH)
    {
        FWPUSH(sd->pStrErrStruct);
    }
    inchi_ios_flush(out_file);
#endif

    return ret;
}


/****************************************************************************
 Special treatment for polymers: perform CRU frame shift analysis 
 and make related edits in orig_inp_data whenever applicable
****************************************************************************/
int PreprocessPolymerCRUData(	struct tagINCHI_CLOCK    *ic,
                                struct tagCANON_GLOBALS  *CG,
                                STRUCT_DATA              *sd,
                                INPUT_PARMS              *ip,
                                char                     *szTitle,
                                PINChI2                  *pINChI2[INCHI_NUM],
                                PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                                INCHI_IOSTREAM           *inp_file,
                                INCHI_IOSTREAM           *log_file,
                                INCHI_IOSTREAM           *out_file,
                                INCHI_IOSTREAM           *prb_file,
                                ORIG_ATOM_DATA           *orig_inp_data,
                                ORIG_ATOM_DATA           *prep_inp_data,
                                long                     num_inp,
                                INCHI_IOS_STRING         *strbuf,
                                unsigned char            save_opt_bits,
                                char					 **sinchi_noedits,
                                char					 **saux_noedits)
{
    int ret = _IS_OKAY;
    char *sinchi_105p = NULL, *saux_105p = NULL;
    OAD_StructureEdits edits_unit_frame_shift, *ed_fs = &edits_unit_frame_shift;
    OAD_StructureEdits edits_unit_folding, *ed_fold = &edits_unit_folding;
    
    OAD_StructureEdits_Init(ed_fold);
    OAD_StructureEdits_Init(ed_fs);

    if (orig_inp_data)
    {
        orig_inp_data->valid_polymer = 0;
        if (orig_inp_data->polymer)
        {
            orig_inp_data->polymer->treat = ip->bPolymers;
            if (orig_inp_data->polymer->treat != POLYMERS_NO)
            {
                orig_inp_data->valid_polymer = 1;

                if (orig_inp_data->polymer->treat == POLYMERS_MODERN)
                {
                    int n_done, n_todo = 0, n_poly_zz = 0; /* djb-rwth: ignoring LLVM warning: variable used */

                    /*  First, get InChI and AuxInfo for the unedited original input
                    (actually, we are interested in AuxInfo for original structure;
                    we then append it to the final result AuxInfo, in order to
                    preserve total reversibility (restoring the original structure ) */
                    ret = OAD_ProcessOneStructureNoEdits(ic, CG, sd, ip, szTitle,
                                                            pINChI2, pINChI_Aux2,
                                                            inp_file, log_file, out_file, prb_file,
                                                            orig_inp_data, prep_inp_data,
                                                            num_inp, strbuf, save_opt_bits,
                                                            &n_poly_zz,
                                                            sinchi_noedits, saux_noedits);
                    if (ret == _IS_FATAL || ret == _IS_ERROR)
                    {
                        ret = _IS_WARNING; 
                        if (!ip->bNoWarnings)
                        {
                            AddErrorMessage(sd->pStrErrStruct, "CRU folding and frame shift analysis failed");
                        }
                        goto exit_function;
                    }
                    if (n_poly_zz < 2)
                    {
                        /* For now, CRU folding and frame shift analysis are only applicable to */
                        /* CRU having both caps of indefinite nature, Zz                        */
                        goto exit_function;
                    }
                    

                    /* Prepare and perform CRU folding related edits */
                    if (ip->bFoldPolymerSRU != 0) 
                    {
                        /*	Get interim 105+ flavour of InChI and AuxInfo and prepare */
                        int old_bFrameShiftScheme = ip->bFrameShiftScheme;
                        ip->bFrameShiftScheme = FSS_STARS_CYCLED;
                        ret = OAD_ProcessOneStructure105Plus(ic, CG, sd, ip, szTitle,
                            pINChI2, pINChI_Aux2,
                            inp_file, log_file, out_file, prb_file,
                            orig_inp_data, prep_inp_data,
                            num_inp, strbuf, save_opt_bits,
                            &sinchi_105p, &saux_105p);
                        ip->bFrameShiftScheme = old_bFrameShiftScheme;
                        if (ret == _IS_FATAL || ret == _IS_ERROR)
                        {
                            ret = _IS_WARNING;
                            if (!ip->bNoWarnings)
                            {
                                /* AddErrorMessage(sd->pStrErrStruct, "CRU fold analysis failed");*/
                                ;
                            }
                            goto frame_shift; 
                        }

                        ret = OAD_Polymer_PrepareFoldCRUEdits( orig_inp_data, *sinchi_noedits, *saux_noedits, sinchi_105p, saux_105p, ed_fold);
                        if (ret == _IS_FATAL || ret == _IS_ERROR)
                        {
                            if (!ip->bNoWarnings)
                            {
                                /*AddErrorMessage(sd->pStrErrStruct, "CRU fold analysis failed");*/
                                ;
                            }
                            goto frame_shift;
                        }
                        if (ret == _IS_WARNING)
                        {
                           /* inchi_ios_eprint(log_file, "Warning (CRU fold analysis failed) structure #%ld.%s%s%s%s\n",
                                num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));*/                            
                            ;
                        }
                        /* else */
                        {
                            /* Proceed with CRU fold */
                            n_done = 0; /* djb-rwth: ignoring LLVM warning: variable used */
                            n_todo = ed_fold->del_atom->used + ed_fold->del_bond->used + ed_fold->new_bond->used + ed_fold->mod_bond->used;
                            ed_fold->del_side_chains = 1;
                            OAD_StructureEdits_DebugPrint(ed_fold);
                            if (n_todo)
                            {
                                /* Edit the original input data */
                                ed_fold->del_side_chains = 1;
                                n_done = OAD_StructureEdits_Apply(sd, ip, orig_inp_data, ed_fold, &ret); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
                                if (ret == _IS_FATAL || ret == _IS_ERROR)
                                {
                                    ret = _IS_WARNING;
                                    /*inchi_ios_eprint(log_file, "Warning (CRU fold failed) structure #%ld.%s%s%s%s\n",
                                        num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));
                                    */
                                    if (!ip->bNoWarnings)
                                    {
                                        AddErrorMessage(sd->pStrErrStruct, "CRU folding failed");
                                    }
                                }
                                else
                                {
                                    if (!ip->bNoWarnings)
                                    {
                                        WarningMessage(sd->pStrErrStruct, "Atom(s) removed due to CRU folding");
                                    }
                                }
                            }
                        }
                    }

frame_shift:        ;
                    /* Prepare and perform frame shift related edits */
                    if (ip->bFrameShiftScheme != FSS_NONE) 
                    {
                        /* Clear buffers */
                        if (sinchi_105p)
                        {
                            inchi_free(sinchi_105p);
                        }
                        if (saux_105p)
                        {
                            inchi_free(saux_105p);
                        }
                        /*	Get interim 105+ flavour of InChI and AuxInfo (possibly 2nd time) */
                        ret = OAD_ProcessOneStructure105Plus(ic, CG, sd, ip, szTitle,
                                                                pINChI2, pINChI_Aux2,
                                                                inp_file, log_file, out_file, prb_file,
                                                                orig_inp_data, prep_inp_data,
                                                                num_inp, strbuf, save_opt_bits,
                                                                &sinchi_105p, &saux_105p);
                        if (ret == _IS_FATAL || ret == _IS_ERROR)
                        {
                            ret = _IS_WARNING;
                            /*inchi_ios_eprint(log_file, "Warning (Frame shift analysis failed) structure #%ld.%s%s%s%s\n",
                                num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));*/
                            if (!ip->bNoWarnings)
                            {
                                AddErrorMessage(sd->pStrErrStruct, "Frame shift analysis failed");
                            }
                            goto exit_function;
                        }

                        ret = OAD_Polymer_PrepareFrameShiftEdits( orig_inp_data, sinchi_105p, saux_105p, ed_fs);
                     
                        if (ret == _IS_FATAL || ret == _IS_ERROR) /* djb-rwth: logical operator corrected */
                        {
                            ret = _IS_WARNING;
                            /*inchi_ios_eprint(log_file, "Warning (Frame shift analysis failed) structure #%ld.%s%s%s%s\n",
                                             num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));*/
                            if (!ip->bNoWarnings)
                            {
                                AddErrorMessage(sd->pStrErrStruct, "Frame shift analysis failed");
                            }
                            goto exit_function;
                        }
                        else
                        {
                            /* OK, proceed with frame shift */
                            n_done = 0; /* djb-rwth: ignoring LLVM warning: variable used */
                            n_todo = ed_fs->del_atom->used +
                                     ed_fs->del_bond->used + ed_fs->new_bond->used + ed_fs->mod_bond->used +
                                     ed_fs->mod_coord->used;
                            OAD_StructureEdits_DebugPrint(ed_fs);
                            if (n_todo)
                            {
                                /* Edit the original input data according to frame shift info */
                                n_done = OAD_StructureEdits_Apply(sd, ip, orig_inp_data, ed_fs, &ret); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
                                if (ret == _IS_FATAL || ret == _IS_ERROR)
                                {
                                    ret = _IS_WARNING;
                                    /*inchi_ios_eprint(log_file, "Warning (Frame shift failed) structure #%ld.%s%s%s%s\n",
                                                     num_inp, SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));*/
                                    if (!ip->bNoWarnings)
                                    {
                                        AddErrorMessage(sd->pStrErrStruct, "Frame shift failed");
                                    }
                                }
                                else
                                {
                                    if (!ip->bNoWarnings)
                                    {
                                        WarningMessage(sd->pStrErrStruct, "Bond(s) rearranged due to CRU frame shift");
                                    }
                                }
                            }
                        }
                    } /* if (ip->bFrameShiftScheme != FSS_NONE)  */


                } /* if (orig_inp_data->polymer->treat == POLYMERS_MODERN) */
            }

            else /* orig_inp_data->polymer->treat == POLYMERS_NO) */
            {
                /*inchi_ios_eprint(log_file, "Ignore polymer data");*/
                if (!ip->bNoWarnings)
                {
                    AddErrorMessage(sd->pStrErrStruct, "Ignore polymer data");
                }
            }
        }
    }

exit_function:
    if (sinchi_105p)
    {
        inchi_free(sinchi_105p);
    }
    if (saux_105p)
    {
        inchi_free(saux_105p);
    }
    OAD_StructureEdits_Clear(ed_fold);  /* Clear edits collection */
    OAD_StructureEdits_Clear(ed_fs);    /* Clear edits collection */

    
    return ret;
}


/****************************************************************************/
void swap_atoms_xyz( ORIG_ATOM_DATA *orig_at_data, int ia1, int ia2 )
{
    double x, y, z;
    
    if (ia1 != ia2)
    {
        x = orig_at_data->at[ia1].x; 
        y = orig_at_data->at[ia1].y; 
        z = orig_at_data->at[ia1].z;

        orig_at_data->at[ia1].x = orig_at_data->at[ia2].x;
        orig_at_data->at[ia1].y = orig_at_data->at[ia2].y;
        orig_at_data->at[ia1].z = orig_at_data->at[ia2].z;
        
        orig_at_data->at[ia2].x = x;
        orig_at_data->at[ia2].y = y;
        orig_at_data->at[ia2].z = z;
    }

    return;
}


/****************************************************************************
 OAD_StructureEdits_Edit
****************************************************************************/
int OAD_StructureEdits_Apply( STRUCT_DATA *sd, 
                              INPUT_PARMS *ip, 
                              ORIG_ATOM_DATA *orig_at_data, 
                              OAD_StructureEdits *ed, 
                              int *ret)
{
    int ok = 0, fail;
    int i, j, old_a1, old_a2, new_a1, new_a2, n_edits = 0;
    int n_del_atom, n_del_bond, n_new_bond, n_mod_bond, n_mod_coord;
    int a1, a2;
    int bond_type = INCHI_BOND_TYPE_NONE, bond_stereo = INCHI_BOND_STEREO_NONE;
    inp_ATOM *at = orig_at_data->at;
    OAD_Polymer *p = orig_at_data->polymer;
    int *at_renum = NULL;
    int *ibuf = NULL;
    int n_max_stored=-1;

    *ret = _IS_OKAY;

    n_del_atom = ed->del_atom->used;
    n_del_bond = ed->del_bond->used/2;
    n_new_bond = ed->new_bond->used/3;
    n_mod_bond = ed->mod_bond->used/4;
    n_mod_coord = ed->mod_coord->used/2;
    if (n_del_atom + n_del_bond + n_new_bond + n_mod_bond + n_mod_coord < 1)
    {
        return 0;
    }

    n_max_stored = inchi_max(2 * (orig_at_data->num_inp_atoms + 1), 2 * (orig_at_data->num_inp_bonds + 1)); /* set all-purpose buffer */

    ITRACE_("\n***************************\nOrig_at_data BEFORE EDITS:\n");
    OrigAtData_DebugTrace(orig_at_data);
    OAD_Polymer_DebugTrace(orig_at_data->polymer);

    /* Delete bonds */
    if (n_del_bond)
    {
        for (i = 0; i < 2 * n_del_bond; i += 2)
        {
            a1 = ed->del_bond->item[i] - 1;
            a2 = ed->del_bond->item[i + 1] - 1;
            ok = OrigAtData_RemoveBond(a1, a2, at, &bond_type, &bond_stereo, &orig_at_data->num_inp_bonds);
            if (!ok)
            {
                *ret = _IS_ERROR;
                goto exit_function;
            }
            n_edits++;
        }
    }
    
    /* Add bonds */
    if (n_new_bond)
    {
        for (i = 0; i < 2 * n_new_bond; i += 2)
        {
            a1 = ed->new_bond->item[i] - 1;
            a2 = ed->new_bond->item[i + 1] - 1;
            /* TODO: consider real bond_type, bond_stereo */
            /* OrigAtData_AddSingleStereolessBond( a1, a2, at, &dummy ); */
            ok = OrigAtData_AddBond(a1, a2, at, bond_type, bond_stereo, &orig_at_data->num_inp_bonds);
            if (!ok)
            {
                *ret = _IS_ERROR;
                goto exit_function;
            }
            n_edits++;
        }
    }
    
    /* Modify bonds */
    if (n_mod_bond)
    {
        for (j = 0; j < 4 * n_mod_bond; j += 4)
        {
            old_a1 = ed->mod_bond->item[j];
            old_a2 = ed->mod_bond->item[j + 1];
            new_a1 = ed->mod_bond->item[j + 2];
            new_a2 = ed->mod_bond->item[j + 3];
            if ((old_a1 == new_a1&&old_a2 == new_a2) || (old_a2 == new_a1&&old_a1 == new_a2))
            {
                continue;
            }
            ok = OrigAtData_RemoveBond(old_a1 - 1, old_a2 - 1, at, &bond_type, &bond_stereo, &orig_at_data->num_inp_bonds);
            if (!ok)
            {
                *ret = _IS_ERROR;
                goto exit_function;
            }
            ok = OrigAtData_AddBond(new_a1 - 1, new_a2 - 1, at, bond_type, bond_stereo, &orig_at_data->num_inp_bonds);
            if (!ok)
            {
                *ret = _IS_ERROR;
                goto exit_function;
            }
            /* Correct CRU blist lists */
            for (i = 0; i < p->n; i++)
            {
                OAD_PolymerUnit *u = p->units[i];
                if (!u->blist)
                {
                    /* No crossing bonds in the unit */
                    continue;
                }
                if ( bIsSameBond(u->blist[0], u->blist[1], old_a1, old_a2) )
                {
                    u->blist[0] = new_a1;
                    u->blist[1] = new_a2;
                }
                else if ( bIsSameBond(u->blist[2], u->blist[3], old_a1, old_a2) )
                {
                    u->blist[2] = new_a1;
                    u->blist[3] = new_a2;
                }
            }

            n_edits++;
        }
    }

    /* Modify coordinates */
    if (n_mod_coord)
    {
        for (j = 0; j < 2 * n_mod_coord; j += 2)
        {
            old_a1 = ed->mod_coord->item[j];
            new_a1 = ed->mod_bond->item[j + 1];
            swap_atoms_xyz(orig_at_data, old_a1 - 1, new_a1 - 1);
            n_edits++;
        }
    }

    
    /* Delete atoms */
    if (n_del_atom)
    {
        int nat0, nat, nacc;
        inp_ATOM *new_at = NULL, *new_at0=NULL;

        at_renum = (int *)inchi_calloc(n_max_stored, sizeof(int));
        if (!at_renum)
        {
            *ret = _IS_ERROR;
            goto exit_function;
        }
        /* all-purpose buffer */
        ibuf = (int *)inchi_calloc(n_max_stored, sizeof(int));
        if (!ibuf)
        {
            *ret = _IS_ERROR;
            goto exit_function;
        }

        fail = mark_atoms_to_delete_or_renumber(orig_at_data, ed, at_renum); 
        if (fail)
        {
            *ret = _IS_ERROR;
            goto exit_function;
        }
        /* Now remove atom by atom */

        nat0 = orig_at_data->num_inp_atoms;
        nat = nat0 - ed->del_atom->used;
        new_at = (inp_ATOM  *)inchi_calloc(nat, sizeof(new_at[0]));
        if (!new_at)
        {
            *ret = _IS_ERROR;
            goto exit_function;
        }

        for (i = 0, nacc = 0; i < nat0; i++)
        {
            AT_NUMB nbr0[MAXVAL];
            U_CHAR btype0[MAXVAL];
            int m, macc, valen;
            int new_num = at_renum[i];				
            if (new_num == -1)
            {
                /* Skip removed atom */
                continue;
            }

            /* Atom to keep; copy it */ 
            new_at0 = new_at + nacc;
            ++nacc;
            memcpy(new_at0, orig_at_data->at + i, sizeof(new_at[0]));
            /* Correct its own number(s) */
            new_at0->orig_at_number = new_num + 1;
            
            /* Correct its nbr number(s) */
            valen = new_at0->valence;
            memcpy(nbr0, new_at0->neighbor, valen * sizeof(AT_NUMB));
            memcpy(btype0, new_at0->bond_type, valen);
            memset(new_at0->neighbor, 0, valen); /* djb-rwth: memset_s C11/Annex K variant? */
            for (m = 0, macc=0; m < valen; m++)
            {
                int num2 = nbr0[m];
                int renum2 = at_renum[num2];
                if (renum2 == num2)
                {
                    /* keep old */
                    new_at0->neighbor[macc++] = num2; 
                }
                else if (renum2 == -1)
                {
                    /* skip and decrement valences */
                    new_at0->chem_bonds_valence -= btype0[m];
                    new_at0->valence--;
                }
                else 
                {
                    /* set renumbered */
                    new_at0->neighbor[macc++] = renum2;
                }
            }
        }

        if (new_at)
        {
            inchi_free(orig_at_data->at);
            orig_at_data->at = new_at;
        }

        orig_at_data->num_inp_atoms = nacc;
        orig_at_data->num_inp_bonds = 0;
        for (i = 0; i < nacc; i++)
        {
            orig_at_data->num_inp_bonds += new_at[i].valence;
        }
        orig_at_data->num_inp_bonds /= 2;

        /* Correct other data */

        /* Correct polymer data */
        if (p)
        {
            int iu;

            for (iu = 0; iu < p->n; iu++)
            {
                OAD_PolymerUnit* u = p->units[iu];
                int new_na, new_nb, new_bb;
                memset(ibuf, 0, n_max_stored * sizeof(ibuf[0])); /* djb-rwth: memset_s C11/Annex K variant? */
                if (u)
                {
                    if (u->alist)
                    {
                        memcpy(ibuf, u->alist, u->na * sizeof(ibuf[0]));
                        new_na = set_renumbered_or_delete(u->alist, ibuf, u->na, at_renum, 1);
                        if (new_na == -1)
                        {
                            *ret = _IS_ERROR;
                            goto exit_function;
                        }
                        u->na = new_na;
                    }
                    if (u->blist)
                    {
                        memcpy(ibuf, u->blist, 2 * (long long)u->nb * sizeof(int)); /* djb-rwth: cast operator added */
                        new_nb = set_renumbered_or_delete(u->blist, ibuf, 2*u->nb, at_renum, 1);
                        new_nb /= 2;
                        if (new_nb == -1)
                        {
                            *ret = _IS_ERROR;
                            goto exit_function;
                        }
                        u->nb = new_nb;
                    }


                    if (u->bkbonds)
                    {
                        int b, new_nbkbonds;
                        for (b = 0, new_nbkbonds = 0; b < u->nbkbonds; b++)
                        {
                            int bnd[2];
                            bnd[0] = u->bkbonds[b][0];
                            bnd[1] = u->bkbonds[b][1];
                            memcpy(ibuf, bnd, 2 * sizeof(ibuf[0]));
                            memset(u->bkbonds[b], 0, 2 * sizeof(u->bkbonds[0][0])); /* djb-rwth: memset_s C11/Annex K variant? */
                            new_bb = set_renumbered_or_delete(bnd, ibuf, 2, at_renum, 1);
                            if (new_bb == -1)
                            {
                                *ret = _IS_ERROR;
                                goto exit_function;
                            }
                            else if (new_bb == 2)
                            {
                                u->bkbonds[new_nbkbonds][0] = bnd[0];
                                u->bkbonds[new_nbkbonds][1] = bnd[1];
                                new_nbkbonds++;
                                /* OK, settled new nums or kept old */
                            }
                            else
                            {
                                continue;
                            }
                        }
                        u->nbkbonds = new_nbkbonds;
                    }

                    if (u->blist) /* djb-rwth: fixing a NULL pointer dereference */
                    {
                        u->cap1 = u->blist[0];
                        u->end_atom1 = u->blist[1];
                        u->cap2 = u->blist[2];
                        u->end_atom2 = u->blist[3];
                    }
                    if (u->cap1 < 1 || u->cap2 < 1 || u->end_atom1 < 1 || u->end_atom2 < 1)
                    {
                        *ret = _IS_ERROR;
                        goto exit_function;
                    }


                }
            }
        } /* if (p) */

        /* Correct V300 data */
        if (orig_at_data->v3000)
        {
            ;
        }


    } /* if (n_del_atom) */
    

exit_function:
    if (ibuf)
    {
        inchi_free(ibuf);
    }
    if (at_renum)
    {
        inchi_free(at_renum);
    }
    return n_edits;
}
                                

/****************************************************************************
 Set each element of number to renum[element] or delete it if renum==(base -1)
 base is either 0 (0-started numbers) or 1 (1-started)
 Returns new number of elements or -1 at error
****************************************************************************/
int set_renumbered_or_delete( int *number,	/* numbers to correct */
                              int *buf,		/* must be enough size to hold nelems elements */
                              int nelems,	/* initial size of numbers */
                              int *renum,	/* new numbers in order of old numbers; always 0-based */
                              int base)
{
    int i, new_nelems;
    memcpy(buf, number, nelems * sizeof(int));
    memset(number, 0, nelems * sizeof(int)); /* djb-rwth: memset_s C11/Annex K variant? */
    for (i = 0, new_nelems = 0; i < nelems; i++)
    {
        int new_num = renum[ buf[i]-base ] + base;
        if (new_num == (base-1))
        {
            continue;
        }
        else
        {
            number[new_nelems++] = new_num;
        }
    }
    return new_nelems;
}


/****************************************************************************
 Worker placed under ProcessOneStructure wrapper (not the last nested doll)
****************************************************************************/
int ProcessOneStructureExCore( struct tagINCHI_CLOCK *ic,
                               struct tagCANON_GLOBALS  *CG,
                               STRUCT_DATA *sd,
                               INPUT_PARMS *ip,
                               char *szTitle,
                               PINChI2 *pINChI2[INCHI_NUM],
                               PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                               INCHI_IOSTREAM *inp_file,
                               INCHI_IOSTREAM *log_file,
                               INCHI_IOSTREAM *out_file,
                               INCHI_IOSTREAM *prb_file,
                               ORIG_ATOM_DATA *orig_inp_data,
                               ORIG_ATOM_DATA *prep_inp_data,
                               long num_inp,
                               INCHI_IOS_STRING *strbuf,
                               unsigned char save_opt_bits )
{
    int res = _IS_OKAY;
    int mind_polymers;

#ifdef TARGET_LIB_FOR_WINCHI
    inchi_ios_free_str( out_file );
    inchi_ios_print(out_file, "Structure: %d\n", num_inp);
#endif

    /* Polymer and pseudoelement specific */
    res = ValidateAndPreparePolymerAndPseudoatoms( ic, CG, sd, ip, szTitle, pINChI2, pINChI_Aux2,
                                                    inp_file, log_file, out_file, prb_file,
                                                    orig_inp_data, prep_inp_data, num_inp, strbuf,
                                                    save_opt_bits, &mind_polymers);
    if (res== _IS_ERROR || res == _IS_FATAL )
    {
        return res;
    }

    /* Call the very actual worker placed under this (ProcessOneStructureCore) wrapper */
    res = ProcessOneStructure( ic, CG, sd, ip, szTitle, pINChI2, pINChI_Aux2,
                               inp_file, log_file, out_file, prb_file,
                               orig_inp_data, prep_inp_data, num_inp, strbuf,
                               save_opt_bits);

    if ( (res == _IS_OKAY  || res == _IS_WARNING ) && mind_polymers )
    {
        /* Post-edit the polymer layer at older polymer treatment modes (1.05, 1.05+) */
        if (ip->bPolymers == POLYMERS_LEGACY || ip->bPolymers == POLYMERS_LEGACY_PLUS)
        {
            /* Cut and hide "Zz" and related things in InChI (AuxInfo has a specifics). */
            int n_pzz = 0, n_zy = orig_inp_data->n_zy;
            if (orig_inp_data->polymer)
            {
                n_pzz = orig_inp_data->polymer->n_pzz;
            }

            EditINCHI_HidePolymerZz(out_file, n_pzz, n_zy);

        }
    }

#ifdef TARGET_LIB_FOR_WINCHI
/*	if ( res == _IS_ERROR || res == _IS_FATAL )
    {
        inchi_ios_print(out_file, "Error %d (%s)\n", sd->nErrorCode, sd->pStrErrStruct);
    }
*/
    /*push_to_winchi_text_window(out_file); */
    /*inchi_ios_free_str(out_file);*/
    /*inchi_ios_flush(out_file);*/
#endif

    return res;
}


/****************************************************************************
 Treat pseudoelement and polymers: parse, validate and set details 
****************************************************************************/
int ValidateAndPreparePolymerAndPseudoatoms( struct tagINCHI_CLOCK *ic,
                                             struct tagCANON_GLOBALS *CG,
                                             STRUCT_DATA *sd,
                                             INPUT_PARMS *ip,
                                             char *szTitle,
                                             PINChI2 *pINChI2[INCHI_NUM],
                                             PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                                             INCHI_IOSTREAM *inp_file,
                                             INCHI_IOSTREAM *log_file,
                                             INCHI_IOSTREAM *out_file,
                                             INCHI_IOSTREAM *prb_file,
                                             ORIG_ATOM_DATA *orig_inp_data,
                                             ORIG_ATOM_DATA *prep_inp_data,
                                             long num_inp,
                                             INCHI_IOS_STRING *strbuf,
                                             unsigned char save_opt_bits,
                                             int *mind_polymers )
{
    int res = _IS_OKAY;

    int mind_pseudoelements = 0;
    
    /* djb-rwth: fixing coverity ID #499512 */
    if (!orig_inp_data)
    {
        goto exit_function;
    }

    *mind_polymers = orig_inp_data->polymer && orig_inp_data->polymer->n > 0;
    *mind_polymers = *mind_polymers && orig_inp_data->valid_polymer &&
        (ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE);
    mind_pseudoelements = (ip->bNPZz == 1) || (ip->bPolymers != POLYMERS_NO);


    /* Validate the data */
    res = OAD_ValidatePolymerAndPseudoElementData(orig_inp_data,
        ip->bPolymers,
        ip->bNPZz,
        sd->pStrErrStruct,
        ip->bNoWarnings);
    if (res)
    {
        sd->nErrorCode = res;
        inchi_ios_eprint( log_file, "Error %d (%s) structure #%ld.%s%s%s%s\n",
                          sd->nErrorCode, sd->pStrErrStruct, num_inp,
                          SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        res = _IS_ERROR;
        orig_inp_data->num_inp_atoms = -1; /* djb-rwth: fixing coverity ID #499522 */
        goto exit_function;
    }

    if (*mind_polymers || mind_pseudoelements)
    {
        /*OrigAtData_DebugTrace(orig_inp_data);*/
        if (*mind_polymers &&
            ip->bPolymers != POLYMERS_MODERN &&
            (ip->bFrameShiftScheme == FSS_STARS_CYCLED || ip->bFrameShiftScheme == FSS_STARS_CYCLED_SORTED))
        {
            /*  Analyze and cyclize frame-shift eligible CRUs using InChI canonical numbers 
                (do this only at older polymer treatment modes 1.05, 1.05+)					
            */
            res = OAD_Polymer_CyclizeCloseableUnits( orig_inp_data,
                                                     ip->bPolymers,
                                                     sd->pStrErrStruct,
                                                     ip->bNoWarnings );
            if (res)
            {
                sd->nErrorCode = res;
                AddErrorMessage(sd->pStrErrStruct, "Error while processing polymer-related input");
                res = _IS_ERROR;
                orig_inp_data->num_inp_atoms = -1;
                goto exit_function;
            }
            /*OrigAtData_DebugTrace(orig_inp_data);*/
        }
    }

exit_function:

    return res;
}


/****************************************************************************
 Get InChI and AuxInfo of totally unedited original structure.
 The intent is to preserve AuxInfo for the very original structure
 in order to keep a final ability to restore that structure. 
****************************************************************************/
int OAD_ProcessOneStructureNoEdits( struct tagINCHI_CLOCK    *ic,
                                    struct tagCANON_GLOBALS  *CG,
                                    STRUCT_DATA              *sd,
                                    INPUT_PARMS              *ip,
                                    char                     *szTitle,
                                    PINChI2                  *pINChI2[INCHI_NUM],
                                    PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                                    INCHI_IOSTREAM           *inp_file,
                                    INCHI_IOSTREAM           *log_file,
                                    INCHI_IOSTREAM           *out_file,
                                    INCHI_IOSTREAM           *prb_file,
                                    ORIG_ATOM_DATA           *orig_inp_data,
                                    ORIG_ATOM_DATA           *prep_inp_data,
                                    long                     num_inp,
                                    INCHI_IOS_STRING         *strbuf,
                                    unsigned char            save_opt_bits,
                                    int                      *n_pzz,
                                    char					 **sinchi, 
                                    char					 **saux)
{
    size_t slen;
    int ret = _IS_OKAY, dup_fail = 0;
    POSEContext dup_context, *dup = &dup_context;

    *n_pzz = 0;
    *sinchi = NULL;
    *saux = NULL;

    /* Make a working copy of all the native input */
    dup_fail = POSEContext_Init(dup, sd, ip, szTitle, pINChI2, pINChI_Aux2,
        inp_file, log_file, out_file, prb_file,
        orig_inp_data, prep_inp_data,
        num_inp, strbuf, save_opt_bits);
    if (dup_fail)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    /* Set necessary for this specific case options */
    dup->orig_inp_data->polymer->treat = dup->ip.bPolymers = POLYMERS_MODERN;
    dup->ip.bFoldPolymerSRU = 0;
    dup->ip.bFrameShiftScheme = FSS_NONE;
    dup->ip.bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO);
    dup->ip.bDisplay = dup->ip.bDisplayCompositeResults = dup->ip.bDisplayEachComponentINChI = 0;
    dup->ip.bFoldPolymerSRU = 0;
    /*dup->ip.bTautFlags |= TG_FLAG_RECONNECT_COORD;*/

    /* Calculate */
    ret = ProcessOneStructureExCore( ic, CG, &dup->sd, &dup->ip, dup->szTitle,
                                     dup->pINChI2, dup->pINChI_Aux2,
                                     dup->inp_file, dup->log_file,
                                     dup->out_file, dup->prb_file,
                                     dup->orig_inp_data, dup->prep_inp_data,
                                     dup->num_inp, dup->strbuf, dup->save_opt_bits );

    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }
    
    *n_pzz = dup->orig_inp_data->polymer->n_pzz;
    /* Extract InChI */
    slen = dup->out_file->s.nUsedLength;
    extract_inchi_substring(sinchi, dup->out_file->s.pStr, slen);
    if (!*sinchi)
    {
        ret = _IS_ERROR;
    }
    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }
    /* Extract AuxInfo */
    slen = dup->out_file->s.nUsedLength;
    extract_auxinfo_substring(saux, dup->out_file->s.pStr, slen);
    if (!*saux)
    {
        ret = _IS_ERROR;
    }
    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }

exit_function:
    POSEContext_Free(dup);

    return ret;
}


/****************************************************************************
Get InChI and AuxInfo of original structure using 1.05+ InChI version.
****************************************************************************/
int OAD_ProcessOneStructure105Plus( struct tagINCHI_CLOCK    *ic,
                                    struct tagCANON_GLOBALS  *CG,
                                    STRUCT_DATA              *sd,
                                    INPUT_PARMS              *ip,
                                    char                     *szTitle,
                                    PINChI2                  *pINChI2[INCHI_NUM],
                                    PINChI_Aux2              *pINChI_Aux2[INCHI_NUM],
                                    INCHI_IOSTREAM           *inp_file,
                                    INCHI_IOSTREAM           *log_file,
                                    INCHI_IOSTREAM           *out_file,
                                    INCHI_IOSTREAM           *prb_file,
                                    ORIG_ATOM_DATA           *orig_inp_data,
                                    ORIG_ATOM_DATA           *prep_inp_data,
                                    long                     num_inp,
                                    INCHI_IOS_STRING         *strbuf,
                                    unsigned char            save_opt_bits,
                                    char					 **sinchi,
                                    char					 **saux)
{
    size_t slen;
    int ret = _IS_OKAY, dup_fail = 0;
    POSEContext dup_context, *dup = &dup_context;

    *sinchi = NULL;
    *saux = NULL;

    /* Make a working copy of all the native input */
    dup_fail = POSEContext_Init(dup, sd, ip, szTitle, pINChI2, pINChI_Aux2,
                                inp_file, log_file, out_file, prb_file,
                                orig_inp_data, prep_inp_data,
                                num_inp, strbuf, save_opt_bits);
    if (dup_fail)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    /* Set necessary for this specific case options */
    /* 1.05+ polymer treatment mode: hidden Zz, senior bkbond comes the first */
    dup->orig_inp_data->polymer->treat = dup->ip.bPolymers = POLYMERS_LEGACY_PLUS;
    /* include full aux info */
    dup->ip.bINChIOutputOptions &= ~(INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO);
    /* request no /D display in inchi-1 executable */
    dup->ip.bDisplay = dup->ip.bDisplayCompositeResults = dup->ip.bDisplayEachComponentINChI = 0;

    /* Calculate */
    ret = ProcessOneStructureExCore(ic, CG, &dup->sd, &dup->ip, dup->szTitle,
                                    dup->pINChI2, dup->pINChI_Aux2,
                                    dup->inp_file, dup->log_file,
                                    dup->out_file, dup->prb_file,
                                    dup->orig_inp_data, dup->prep_inp_data,
                                    dup->num_inp, dup->strbuf, dup->save_opt_bits);

    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }

    /* Extract InChI */
    slen = dup->out_file->s.nUsedLength;
    extract_inchi_substring(sinchi, dup->out_file->s.pStr, slen);
    if (!*sinchi)
    {
        ret = _IS_ERROR;
    }
    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }
    /* Extract AuxInfo */
    slen = dup->out_file->s.nUsedLength;
    extract_auxinfo_substring(saux, dup->out_file->s.pStr, slen);
    if (!*saux)
    {
        ret = _IS_ERROR;
    }
    if (ret == _IS_FATAL || ret == _IS_ERROR)
    {
        goto exit_function;
    }


exit_function:
    POSEContext_Free(dup);

    return ret;
}


/****************************************************************************
 Mark atoms to be renumbered or removed
 at_renum[initial number] := new_number or -1 (mark of deletion)
****************************************************************************/
int mark_atoms_to_delete_or_renumber( ORIG_ATOM_DATA *orig_at_data,
                                      OAD_StructureEdits *ed,
                                      int *at_renum)
{
    int i, j;
    int fail=0, ret = 0;
    size_t *atnums = NULL; /* djb-rwth: needs to be size_t type */
    size_t max_atoms = orig_at_data->num_inp_atoms;

    /* NB:	new/old ORIG_ATOM_DATA atom numbers are 0-based (==orig_number-1) 
            while those in ed->... are just 1-based orig_numbers */
    
    for (i = 0; (size_t)i < max_atoms; i++)
    {
        at_renum[i] = i;
    }

    if (ed->del_side_chains)
    {
        /* Extend list of atoms to be deleted with those connected to original ones
        (i.e., delete a whole connected component(s) comprising original atoms)
        */
        int natnums = 0;
        atnums = (size_t *)inchi_calloc(max_atoms, sizeof(size_t)); /* djb-rwth: size_t type used for max_atoms to fit the definition of inchi_calloc  */
        if (!atnums)
        {
            return _IS_ERROR;
        }
        for (i = 0; (size_t)i < max_atoms; i++)
        {
            int iatom = ed->del_atom->item[i] - 1;
            subgraf *sg = NULL;
            subgraf_pathfinder *spf = NULL;
            if (i >= ed->del_atom->used) /* NB: ed->del_atom->used may be increased within this loop */
            {
                break;
            }
            for (j = 0; (size_t)j < max_atoms; j++)
            {
                atnums[j] = orig_at_data->at[j].orig_at_number; /*j+1;*/
            }
            sg = subgraf_new(orig_at_data, max_atoms, (int*)atnums);
            memset(atnums, 0, max_atoms * sizeof(int)); /* djb-rwth: memset_s C11/Annex K variant? */
            if (!sg)
            {
                ret = _IS_ERROR;
                goto exit_function;
            }
            spf = subgraf_pathfinder_new(sg, orig_at_data, iatom, iatom);
            if (!spf)
            {
                ret = _IS_ERROR;
                goto exit_function;
            }
            spf->start = iatom;
            spf->nseen = 0;
            natnums = subgraf_pathfinder_collect_all(spf, 0, NULL, (int*)atnums);
            if (natnums)
            {
                for (j = 0; j < natnums && j < max_atoms; j++) /* djb-rwth: fixing buffer overruns */
                {
                    fail = IntArray_AppendIfAbsent(ed->del_atom, atnums[j]);
                    if (fail)
                    {
                        ret = _IS_ERROR;
                        goto exit_function;
                    }
                }
            }
            subgraf_free(sg);
            subgraf_pathfinder_free(spf);
        }

    } /* if (ed->del_side_chains) */

    for (i = max_atoms - 1; i >= 0; i--)
    {
        int orig_num = i + 1;	/* NB: ed->del_atom->item contains orig# which are (OAD# + 1) */
        if (is_in_the_ilist(ed->del_atom->item, orig_num, ed->del_atom->used)) 
        {
            /* mark as deleted atnum */
            at_renum[i] = -1;
            /* shift other atnums */
            for (j = max_atoms - 1; j > i; j--)
            {
                if (at_renum[j] != -1)
                {
                    at_renum[j]--;
                }
            }
        }
    }


exit_function:
    if (atnums)
    {
        inchi_free(atnums);
    }
    return ret;
}


#if (BUILD_WITH_ENG_OPTIONS==1)
#if ALLOW_SUBSTRUCTURE_FILTERING==1
int check_presence_of_the_encoded_substructure(ORIG_ATOM_DATA *orig_inp_data);
int OrigAtData_CheckForSubstructure(ORIG_ATOM_DATA *orig_inp_data)
{
    int ok = 0;

    ok = check_presence_of_the_encoded_substructure(orig_inp_data);

    return ok;
}
int check_presence_of_the_encoded_substructure(ORIG_ATOM_DATA *oad)
{
    int i, ok;

    /* Place sub-structure filtering code below.

       Return 1 if structure matches some hard-coded pattern
       
       In this example pattern is a presence of a pseudo atom.
    */

    ok = 0;
    for (i = 0; i < oad->num_inp_atoms; i++)
    {
        if (oad->at[i].el_number == EL_NUMBER_ZZ || oad->at[i].el_number == EL_NUMBER_ZY)
        {
            ok = 1;
            break;
        }
    }
    return ok;
}
#endif
#endif
