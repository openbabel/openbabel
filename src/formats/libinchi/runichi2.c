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
    Reading input data

    Prepare and make edits

*/

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>

#include "mode.h"
#include "ichitime.h"
#ifndef COMPILE_ANSI_ONLY
#include <conio.h>
#endif

#include "ichister.h"

#include "ichimain.h"
#include "ichi_io.h"
#include "mol_fmt.h"
#include "inchi_api.h"
#include "readinch.h"
#ifdef TARGET_LIB_FOR_WINCHI
#include "../../../IChI_lib/src/ichi_lib.h"
#include "inchi_api.h"
#else
#include "inchi_gui.h"
#endif
#include "readinch.h"

#include "ichirvrs.h"

#include "bcf_s.h"

/* Modified in-CRU stereoventer info */
typedef struct tagModSCenterInfo
{
    int num;			/* atnums are 0-based, internal ones */
    int valence;
    int n_stereo;
    int nbr[MAXVAL];
    int new_nbr[MAXVAL];
} ModSCenterInfo;

typedef struct tagDiylFrag
{
    int na;			/* number of atoms	*/
    int nb;			/* number of bonds	*/
    int end1;		/* "starting" end	*/
    int end2;		/* "ending" end		*/
    int *alist;		/* list of atoms orig numbers	    */
    int *xclist;	/* list of atoms extended classes   */
    INCHI_IOS_STRING sig; /* string signature           */
} DiylFrag;

static DiylFrag* DiylFrag_New( int na, int e1nd, int end2, char *s );
static void DiylFrag_Free( DiylFrag *pfrag );
static void DiylFrag_MakeSignature( DiylFrag *pfrag, int nxc, int *xc, int *tmp );
static int DiylFrag_Diff( DiylFrag *pfrag1, DiylFrag *pfrag2 );
static void DiylFrag_DebugTrace( DiylFrag *pfrag );

static int NDefStereoBonds( inp_ATOM *at, int iatom, int bOnlyPointedEndMatters );

static void ModSCenter_Init( ModSCenterInfo *scinfo, inp_ATOM *at, int iatom );
static void ModSCenter_AddTo( ModSCenterInfo *scinfo, int iadd );
static void ModSCenter_DelFrom( ModSCenterInfo *scinfo, int idel );
static int ModSCenter_IsChanged( ModSCenterInfo *scinfo, inp_ATOM *at );

static int GetFrameShiftInfoFrom105PlusInChI( char *sinchi, int *alist, int max_crossing );
/* static int IsPolymerRequiringEdits(ORIG_ATOM_DATA* orig_inp_data); */ /* djb-rwth: function definition not found*/
static int analyze_CRU_folding( ORIG_ATOM_DATA *orig_at_data,
                                int iunit,
                                int n_all_bkb,
                                int *all_bkb,
                                int nxclasses,
                                int *xc,
                                OAD_StructureEdits *ed );
static int count_colors_in_sequence( int *entries, int n_entries, int max_distinct, int *counts );
static int len_repeating_subsequence( int *color, int *color2, int n );



/****************************************************************************
 Get (the next) one portion of input data of any possible kind
 (Molfile, InChI string, ...) from a sequential input stream
****************************************************************************/
int GetOneStructure( INCHI_CLOCK    *ic,
                     STRUCT_DATA    *sd,
                     INPUT_PARMS    *ip,
                     char           *szTitle,
                     INCHI_IOSTREAM *inp_file,
                     INCHI_IOSTREAM *log_file,
                     INCHI_IOSTREAM *out_file,
                     INCHI_IOSTREAM *prb_file,
                     ORIG_ATOM_DATA *orig_inp_data,
                     long           *num_inp,
                     STRUCT_FPTRS   *struct_fptrs )
{
    int nRet, inp_index, out_index, bUseFptr = ( NULL != struct_fptrs );

    FreeOrigAtData( orig_inp_data );

    /* added for TARGET_LIB_FOR_WINCHI early EOF detection */
    inp_index = -1;
    out_index = -1;

    if (struct_fptrs)
    {
        if (inp_file->f == stdin)
        {
            return _IS_FATAL;
        }
        if (ip->nInputType == INPUT_CMLFILE)
        {
            bUseFptr = 0;
        }

        /* Initially allocate or increase length of struct_fptrs->fptr array */
        if (!struct_fptrs->fptr || struct_fptrs->len_fptr <= struct_fptrs->cur_fptr + 1)
        {

            INCHI_FPTR *new_fptr = (INCHI_FPTR *) 
                                        inchi_calloc( (long long)struct_fptrs->len_fptr + ADD_LEN_STRUCT_FPTRS, sizeof( new_fptr[0] ) ); /* djb-rwth: cast operator added */

            if (new_fptr)
            {
                if (struct_fptrs->fptr)
                {
                    if (struct_fptrs->len_fptr)
                    {
                        memcpy(new_fptr, struct_fptrs->fptr, struct_fptrs->len_fptr * sizeof(new_fptr[0]));
                    }
                    inchi_free( struct_fptrs->fptr );
                }
                else
                {
                    struct_fptrs->len_fptr = 0;
                    struct_fptrs->cur_fptr = 0;
                    struct_fptrs->max_fptr = 0;
                }
                struct_fptrs->len_fptr += ADD_LEN_STRUCT_FPTRS;
                struct_fptrs->fptr = new_fptr;
            }
            else
            {
                return _IS_FATAL;  /* new_fptr allocation error */
            }
        }

        if (struct_fptrs->fptr[struct_fptrs->cur_fptr] == EOF)
        {
            return _IS_EOF;
        }
        else
        {

            if (bUseFptr)
            {
                if (fseek( inp_file->f,
                    struct_fptrs->fptr[struct_fptrs->cur_fptr],
                    SEEK_SET ))
                {
                    return _IS_FATAL;
                }
                if (struct_fptrs->cur_fptr &&
                     struct_fptrs->max_fptr <= struct_fptrs->cur_fptr)
                {
                    return _IS_FATAL;
                }
            }
            else
            {
                inp_index = struct_fptrs->fptr[struct_fptrs->cur_fptr];
                out_index = EOF;
            }
        }

        *num_inp = struct_fptrs->cur_fptr; /* set structure count */
    }


    nRet = ReadTheStructure( ic, sd, ip, inp_file, orig_inp_data, inp_index, &out_index );

    if (!nRet)
    {
        if (ip->nInputType == INPUT_INCHI_PLAIN || ip->nInputType == INPUT_MOLFILE || ip->nInputType == INPUT_SDFILE)
        {
            if (ip->lMolfileNumber)
            {
                *num_inp = ip->lMolfileNumber;
            }
            else
            {
                *num_inp += 1;
            }
        }
        else
        {
            *num_inp += 1;
        }

        nRet = TreatErrorsInReadTheStructure( sd, ip, LOG_MASK_ALL,
                                              inp_file, log_file, out_file, prb_file,
                                              orig_inp_data, num_inp );
    }

    /************************************************************
     Added for TARGET_LIB_FOR_WINCHI:
     look ahead for end of file detection
    ************************************************************/
    if ( inp_file->type == INCHI_IOS_TYPE_FILE              &&
         inp_file->f && struct_fptrs && struct_fptrs->fptr  &&
         struct_fptrs->fptr[struct_fptrs->cur_fptr + 1] <= 0 )
    {

        int nRet2 = 0;
        INCHI_FPTR next_fptr = 0;
        STRUCT_DATA sd2;

        if (nRet != _IS_EOF && nRet != _IS_FATAL)
        {
            if (inp_file->f == stdin || struct_fptrs->len_fptr <= struct_fptrs->cur_fptr + 1)
            {
                return _IS_FATAL;
            }
            /* Get the next structure fptr */
            if (bUseFptr)
            {
                next_fptr = ftell( inp_file->f );
            }
            else
            {
                inp_index = out_index;
                out_index = EOF;
            }

            /* Read the next structure */
            nRet2 = ReadTheStructure( ic, &sd2, ip, inp_file,
                                      NULL, inp_index, &out_index );

            /* Restore fptr to the next structure */
            if (bUseFptr)
            {
                if (next_fptr != -1L)
                {
                    fseek( inp_file->f, next_fptr, SEEK_SET );
                }
            }
        }
        else
        {
            /* Treat current fatal error as end of file */
            struct_fptrs->fptr[struct_fptrs->cur_fptr] = EOF;
        }

        /* Next is end of file or fatal */
        if (nRet == _IS_EOF || nRet == _IS_FATAL ||
             nRet2 == _IS_EOF || nRet2 == _IS_FATAL)
        {
            struct_fptrs->fptr[struct_fptrs->cur_fptr + 1] = EOF;
        }
        else
        {
            struct_fptrs->fptr[struct_fptrs->cur_fptr + 1] = bUseFptr ? sd->fPtrEnd : inp_index;
        }

        /* Update struct_fptrs->max_fptr */
        if (struct_fptrs->max_fptr <= struct_fptrs->cur_fptr + 1)
        {
            struct_fptrs->max_fptr = struct_fptrs->cur_fptr + 2;
        }
    }

    switch (nRet)
    {
        case _IS_EOF:
            *num_inp -= 1;
        case _IS_FATAL:
        case _IS_ERROR:
        case _IS_SKIP:
            goto exit_function;
    }

    /*
    if ( !orig_inp_data->num_dimensions ) {
        WarningMessage(sd->pStrErrStruct, "0D"); */ /* 0D-structure: no coordinates
    }
    */

exit_function:
    return nRet;
}


/****************************************************************************
 Extract one connected component from the input structure
****************************************************************************/
int GetOneComponent( INCHI_CLOCK        *ic,
                     STRUCT_DATA        *sd,
                     INPUT_PARMS        *ip,
                     INCHI_IOSTREAM     *log_file,
                     INCHI_IOSTREAM     *out_file,
                     INP_ATOM_DATA      *inp_cur_data,
                     ORIG_ATOM_DATA     *orig_inp_data,
                     int                i,
                     long               num_inp )
{
    inchiTime ulTStart;

    InchiTimeGet( &ulTStart );

    CreateInpAtomData( inp_cur_data, orig_inp_data->nCurAtLen[i], 0 );

    inp_cur_data->num_at = ExtractConnectedComponent( orig_inp_data->at, orig_inp_data->num_inp_atoms, i + 1, inp_cur_data->at );

    sd->ulStructTime += InchiTimeElapsed( ic, &ulTStart );


    /*  Error processing */
    if (inp_cur_data->num_at <= 0 || orig_inp_data->nCurAtLen[i] != inp_cur_data->num_at)
    {
        /*  Log error message */
        AddErrorMessage( sd->pStrErrStruct, "Cannot extract Component" );
        inchi_ios_eprint( log_file, "%s #%d structure #%ld.%s%s%s%s\n", sd->pStrErrStruct, i + 1, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        sd->nErrorCode = inp_cur_data->num_at < 0 ? inp_cur_data->num_at : ( orig_inp_data->nCurAtLen[i] != inp_cur_data->num_at ) ? CT_ATOMCOUNT_ERR : CT_UNKNOWN_ERR;
        /* num_err ++; */
        sd->nErrorType = _IS_ERROR;
#ifdef TARGET_LIB_FOR_WINCHI
        if (( ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) &&
            ( ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT ))
        {
            sd->nErrorType = ProcessStructError( out_file,
                                                 log_file,
                                                 sd->pStrErrStruct,
                                                 sd->nErrorType,
                                                 num_inp,
                                                 ip );
        }
#endif
    }

    return sd->nErrorType;
}


/****************************************************************************
 Read input data of any possible kind (Molfile, InChI string, ...)
****************************************************************************/
int ReadTheStructure( struct tagINCHI_CLOCK *ic,
                      STRUCT_DATA           *sd,
                      INPUT_PARMS           *ip,
                      INCHI_IOSTREAM        *inp_file,
                      ORIG_ATOM_DATA        *orig_inp_data,
                      /* the further is deprecated (CML support)        */
                      int                   inp_index,
                      int                   *out_index )
{
    inchiTime ulTStart;
    int nRet = 0, nRet2 = 0;
    int bGetOrigCoord = !( ip->bINChIOutputOptions &
        ( INCHI_OUT_NO_AUX_INFO | INCHI_OUT_SHORT_AUX_INFO ) );

    INCHI_MODE InpAtomFlags = 0;
    /* NB: reading Molfile may set FLAG_INP_AT_CHIRAL bit */

    int vABParityUnknown = AB_PARITY_UNDF;
        /* vABParityUnknown holds actual value of an internal constant */
        /* signifying unknown parity: either the same as for undefined */
        /* parity (default==standard)  or a specific one (non-std;       */
        /* requested by SLUUD switch).                                   */

    if (0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO ))
    {
        /* Make labels for unknown and undefined stereo different */
        vABParityUnknown = AB_PARITY_UNKN;
    }

    if (ip->bLargeMolecules)
    {
        InpAtomFlags |= FLAG_SET_INP_LARGE_MOLS;
    }

    memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    switch (ip->nInputType)
    {

        case INPUT_MOLFILE:
        case INPUT_SDFILE:

            /*  Read the original input structure from Molfile */
            if (orig_inp_data)
            {

                if (ip->pSdfValue && ip->pSdfValue[0])
                {
                    /* Added 07-29-2003 to avoid inheriting exact value from prev.
                        structure and to make reference to a (bad) structure
                        with unknown ID Value                                     */
                    char *p, *q;  /* q shadows prev declaration of const char *q */
                    int  n;
                    if (( p = strrchr( ip->pSdfValue, '+' ) ) &&
                         '[' == *( p - 1 ) &&
                         0 < ( n = strtol( p + 1, &q, 10 ) ) &&
                         q[0] &&
                         ']' == q[0] &&
                         !q[1])
                    {
                        sprintf(p + 1, "%d]", n + 1);
                    }
                    else
                    {
                        strcat(ip->pSdfValue, " [+1]");
                    }
                }

                InchiTimeGet( &ulTStart );

                if (inp_file->type == INCHI_IOS_TYPE_FILE && inp_file->f)
                    sd->fPtrStart = ( inp_file->f == stdin ) ? -1 : ftell( inp_file->f );


                nRet2 = CreateOrigInpDataFromMolfile( inp_file,
                                                      orig_inp_data,
                                                      ip->bMergeAllInputStructures,
                                                      bGetOrigCoord,
                                                      ip->bDoNotAddH,
                                                      ip->bPolymers,
                                                      ip->bNPZz,
                                                      ip->pSdfLabel,
                                                      ip->pSdfValue,
                                                      &ip->lSdfId,
                                                      &ip->lMolfileNumber,
                                                      &InpAtomFlags,
                                                      &sd->nStructReadError,
                                                      sd->pStrErrStruct,
                                                      ip->bNoWarnings); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */



                if (!ip->bGetSdfileId || ip->lSdfId == 999999LU)
                {
                    ip->lSdfId = 0LU;
                }

                if (!ip->bGetMolfileNumber || ip->lMolfileNumber < 0)
                {
                    ip->lMolfileNumber = 0;
                }

                if (inp_file->type == INCHI_IOS_TYPE_FILE && inp_file->f)
                {
                    sd->fPtrEnd = ( inp_file->f == stdin ) ? -1 : ftell( inp_file->f );
                }

                sd->ulStructTime += InchiTimeElapsed( ic, &ulTStart );

#if ( bRELEASE_VERSION == 0 )
                sd->bExtract |= orig_inp_data->bExtract;
#endif

            /* 2004-11-16: added Molfile Chiral Flag Mode */
            /* ********************************************
             * Chiral flags are set in:
             * - RunICHI.c -- ReadTheStructure()
             * - e_IchiMain.c -- main()
             * - inchi_dll.c -- ExtractOneStructure
             **********************************************/

            /* 1. Highest precedence: Chiral Flag set by the user */
                if (ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL)
                {
                    InpAtomFlags = FLAG_INP_AT_CHIRAL; /* forced by the user */
                }
                else if (ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL)
                {
                    InpAtomFlags = FLAG_INP_AT_NONCHIRAL; /* forced by the user */
                }
                else if (( InpAtomFlags & FLAG_INP_AT_CHIRAL ) && ( InpAtomFlags & FLAG_INP_AT_NONCHIRAL )) /* djb-rwth: correcting &&->& for bitwise calculation? */
                {
                    InpAtomFlags &= ~FLAG_INP_AT_NONCHIRAL;
                }

                /* Save requested flags in the AuxInfo */
                sd->bChiralFlag &= ~( FLAG_INP_AT_CHIRAL | FLAG_INP_AT_NONCHIRAL );
                sd->bChiralFlag |= InpAtomFlags & ( FLAG_INP_AT_CHIRAL | FLAG_INP_AT_NONCHIRAL );

                /* Quick fix: modify ip->nMode on the fly */

                /* 2. The user requested both Stereo AND Chiral flag */
                if (( ip->nMode & REQ_MODE_CHIR_FLG_STEREO ) && ( ip->nMode & REQ_MODE_STEREO ))
                {
                    if (InpAtomFlags & FLAG_INP_AT_CHIRAL)
                    {
                        /* structure has chiral flag or the user said it is chiral */
                        ip->nMode &= ~( REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO );
                        sd->bChiralFlag |= FLAG_INP_AT_CHIRAL; /* write AuxInfo as chiral */
                    }
                    else
                    {
                        ip->nMode &= ~REQ_MODE_RACEMIC_STEREO;
                        ip->nMode |= REQ_MODE_RELATIVE_STEREO;
                        sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL; /* write AuxInfo as explicitly not chiral */
                    }
                }
            } /* if ( orig_inp_data )  */

            else    /* !orig_inp_data */
            {
                /*  read the next original structure */
                int nStructReadError = 0;

                if (!ip->bMergeAllInputStructures)
                {
                    nRet2 = CreateOrigInpDataFromMolfile( inp_file,
                                                          NULL,     /* *orig_at_data, */
                                                          0,        /* bMergeAllInputStructures */
                                                          0,        /* bGetOrigCoord */
                                                          0,        /* bDoNotAddH */
                                                          0, /* ip->bPolymers */
                                                          0, /* ip->bNPZz */
                                                          NULL,     /* *pSdfLabel */
                                                          NULL,     /* *pSdfValue */
                                                          NULL,     /* *lSdfId */
                                                          NULL,     /* *lMolfileNumber */
                                                          &InpAtomFlags, /*NULL, */
                                                          &nStructReadError,
                                                          NULL,      /* *pStrErr */
                                                          0);

                    if (nRet2 <= 0 && 10 < nStructReadError && nStructReadError < 20)
                    {
                        return _IS_EOF;
                    }
                }
                else
                {
                    return _IS_EOF;
                }
            }

            break;

        case INPUT_INCHI_PLAIN:
            /*  Read the original input data  as text (InChI string ) */
            if (orig_inp_data)
            {
                if (ip->pSdfValue && ip->pSdfValue[0])
                {
                    /* Added 07-29-2003 to avoid inheriting exact value from prev. structure
                       and to make reference to a (bad) structure with unknown ID Value */
                    char *p, *q;
                    int  n;
                    if (( p = strrchr( ip->pSdfValue, '+' ) ) &&
                         '[' == *( p - 1 ) &&
                         0 < ( n = strtol( p + 1, &q, 10 ) ) &&
                         q[0] &&
                         ']' == q[0] &&
                         !q[1])
                    {
                        sprintf(p + 1, "%d]", n + 1);
                    }
                    else
                    {
                        strcat(ip->pSdfValue, " [+1]");
                    }
                }

                InchiTimeGet( &ulTStart );

                if (inp_file->type == INCHI_IOS_TYPE_FILE && inp_file->f)
                {
                    sd->fPtrStart = ( inp_file->f == stdin ) ? -1 : ftell( inp_file->f );
                }


                /*  Read and make internal molecular  data */
                nRet2 = InchiToOrigAtom( inp_file,
                                         orig_inp_data,
                                         ip->bMergeAllInputStructures,
                                         bGetOrigCoord,
                                         ip->bDoNotAddH,
                                         vABParityUnknown,
                                         ip->nInputType,
                                         ip->pSdfLabel,
                                         ip->pSdfValue,
                                         (unsigned long *) &ip->lMolfileNumber,
                                         &InpAtomFlags,
                                         &sd->nStructReadError,
                                         sd->pStrErrStruct ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

                /*if ( !ip->bGetSdfileId || ip->lSdfId == 999999LU) ip->lSdfId = 0;*/
                if (inp_file->type == INCHI_IOS_TYPE_FILE && inp_file->f)
                {
                    sd->fPtrEnd = ( inp_file->f == stdin ) ? -1 : ftell( inp_file->f );
                }

                sd->ulStructTime += InchiTimeElapsed( ic, &ulTStart );

#if ( bRELEASE_VERSION == 0 )
                sd->bExtract |= orig_inp_data->bExtract;
#endif

                /* 2004-11-16: added Molfile Chiral Flag Mode */
                if (ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL)
                {
                    InpAtomFlags = FLAG_INP_AT_CHIRAL; /* forced by the user */
                }
                else if (ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL)
                {
                    InpAtomFlags = FLAG_INP_AT_NONCHIRAL; /* forced by the user */
                }
                else if (( InpAtomFlags & FLAG_INP_AT_CHIRAL ) && ( InpAtomFlags & FLAG_INP_AT_NONCHIRAL )) /* djb-rwth: correcting &&->& for bitwise calculation? */
                {
                    InpAtomFlags &= ~FLAG_INP_AT_NONCHIRAL;
                }

                sd->bChiralFlag |= InpAtomFlags; /* copy chiral flag to AuxInfo */

                /* Quick fix: modify ip->nMode on the fly */
                if (( ip->nMode & REQ_MODE_CHIR_FLG_STEREO ) && ( ip->nMode & REQ_MODE_STEREO ))
                {
                    if (InpAtomFlags & FLAG_INP_AT_CHIRAL)
                    {
                        ip->nMode &= ~( REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO );
                    }
                    else
                    {
                        ip->nMode &= ~REQ_MODE_RACEMIC_STEREO;
                        ip->nMode |= REQ_MODE_RELATIVE_STEREO;
                    }
                }
            }
            else
            {
                /*  Read the next original structure */
                int           nStructReadError = 0;
                if (!ip->bMergeAllInputStructures)
                {
                    nRet2 = InchiToOrigAtom( inp_file,
                                             NULL, 0, 0, 0, 0,
                                             ip->nInputType,
                                             NULL, NULL, NULL, NULL,
                                             &nStructReadError,
                                             NULL );
                    if (nRet2 <= 0 && 10 < nStructReadError && nStructReadError < 20)
                    {
                        return _IS_EOF;
                    }
                }
                else
                {
                    return _IS_EOF;
                }
            }
            break;

        default:
            nRet = _IS_FATAL; /*  wrong file type */
    }

    return nRet;
}


/****************************************************************************
 Interpret and treat input reading errors/warnings
****************************************************************************/
int TreatErrorsInReadTheStructure( STRUCT_DATA      *sd,
                                   INPUT_PARMS      *ip,
                                   int              nLogMask,
                                   INCHI_IOSTREAM   *inp_file,
                                   INCHI_IOSTREAM   *log_file,
                                   INCHI_IOSTREAM   *out_file,
                                   INCHI_IOSTREAM   *prb_file,
                                   ORIG_ATOM_DATA   *orig_inp_data,
                                   long             *num_inp )
{
    int nRet = _IS_OKAY;

    if (10 < sd->nStructReadError && sd->nStructReadError < 20)
    {
        /*  End of file */
        if (sd->pStrErrStruct[0])
        {
            inchi_ios_eprint( log_file, "%s inp structure #%ld: End of file.%s%s%s%s    \n", sd->pStrErrStruct, *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        }

        inchi_ios_eprint( log_file, "End of file detected after structure #%ld.   \n", *num_inp - 1 );

        nRet = _IS_EOF;
        goto exit_function; /*  end of file */
    }

    /*(*num_inp) ++;*/

    if (*num_inp < ip->first_struct_number)
    {
        /*  Skip the structure */

#if ( !defined(TARGET_API_LIB) && !defined(TARGET_EXE_STANDALONE) )
/* #ifndef TARGET_API_LIB */
        if (log_file->f != stderr)
            inchi_fprintf( stderr, "\rSkipping structure #%ld.%s%s%s%s...", *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
#endif

        nRet = sd->nErrorType = _IS_SKIP;
        goto exit_function;
    }

    sd->nErrorType = GetInpStructErrorType( ip, sd->nStructReadError,
                                            sd->pStrErrStruct,
                                            orig_inp_data->num_inp_atoms );

    if (sd->nErrorType == _IS_FATAL)
    {
        /*  Fatal error */
        if (nLogMask & LOG_MASK_FATAL)
        {
            inchi_ios_eprint( log_file, "Fatal Error %d (aborted; %s) inp structure #%ld.%s%s%s%s\n",
                     sd->nStructReadError, sd->pStrErrStruct, *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        }

#if ( bRELEASE_VERSION == 1 || EXTR_FLAGS == 0 )
        /* djb-rwth: fixing oss-fuzz issue #27902 */
        if (prb_file && prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && !ip->bSaveAllGoodStructsAsProblem)
        {
            MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, *num_inp );
        }
#endif
        /* goto exit_function; */
    }

    if (sd->nErrorType == _IS_ERROR)
    {
        /*  Non-fatal errors: do not produce INChI */

        /*  70 => too many atoms */
        if (nLogMask & LOG_MASK_ERR)
        {
            inchi_ios_eprint( log_file, "Error %d (no %s; %s) inp structure #%ld.%s%s%s%s\n",
                 sd->nStructReadError, ( ip->bINChIOutputOptions & INCHI_OUT_SDFILE_ONLY ) ? "Molfile" : INCHI_NAME,
                 sd->pStrErrStruct, *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        }

#if ( bRELEASE_VERSION == 1 || EXTR_FLAGS == 0 )
        if (prb_file && prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd && !ip->bSaveAllGoodStructsAsProblem)
        {
            MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, *num_inp ); /* djb-rwth: addressing coverity ID #499477 -- return values handled properly */
        }
#endif
    }

    if (sd->nErrorType == _IS_WARNING)
    {
        /*  Warnings: try to produce INChI */

        if (nLogMask & LOG_MASK_WARN)
        {
            inchi_ios_eprint( log_file, "Warning: (%s) inp structure #%ld.%s%s%s%s\n",
               sd->pStrErrStruct, *num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        }
    }


#ifdef TARGET_LIB_FOR_WINCHI
    if (( ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) &&
        ( ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT ))
    {
        if (sd->nErrorType != _IS_OKAY && sd->nErrorType != _IS_WARNING)
        {
            sd->nErrorType =
                ProcessStructError( out_file, log_file, sd->pStrErrStruct, sd->nErrorType, *num_inp, ip );
        }
    }
#endif

exit_function:
    if (nRet <= _IS_OKAY && sd->nErrorType > 0)
    {
        nRet = sd->nErrorType;
    }

    return nRet;
}


#ifdef TARGET_EXE_USING_API


/****************************************************************************/
int InchiToInchi_Input( INCHI_IOSTREAM  *inp_file,
                        inchi_Input     *orig_at_data,
                        int             bMergeAllInputStructures,
                        int             bDoNotAddH,
                        int             vABParityUnknown,
                        INPUT_TYPE      nInputType,
                        char            *pSdfLabel,
                        char            *pSdfValue,
                        unsigned long   *lSdfId,
                        INCHI_MODE      *pInpAtomFlags,
                        int             *err,
                        char            *pStrErr )
{
/* inp_ATOM       *at = NULL; */
    int             num_dimensions_new;
    int             num_inp_bonds_new;
    int             num_inp_atoms_new;
    int             num_inp_0D_new;
    inp_ATOM     *at_new = NULL;
    inp_ATOM     *at_old = NULL;
    inchi_Stereo0D *stereo0D_new = NULL;
    inchi_Stereo0D *stereo0D_old = NULL;
    int             nNumAtoms = 0, nNumStereo0D = 0;
    MOL_COORD      *szCoordNew = NULL;
    MOL_COORD      *szCoordOld = NULL;
    int            i, j;
    int               max_num_at;


    max_num_at = MAX_ATOMS;
    if (*pInpAtomFlags  & FLAG_SET_INP_LARGE_MOLS)
        max_num_at = MAX_ATOMS_LARGE_MOL;

    if (pStrErr)
        pStrErr[0] = '\0';


    /*FreeOrigAtData( orig_at_data );*/
    if (lSdfId)
        *lSdfId = 0LU;

    do
    {
        at_old = orig_at_data ? orig_at_data->atom : NULL; /*  save pointer to the previous allocation */

        stereo0D_old = orig_at_data ? orig_at_data->stereo0D : NULL;

        szCoordOld = NULL;

        num_inp_atoms_new =
            InchiToinp_ATOM( inp_file, orig_at_data ? &stereo0D_new : NULL, &num_inp_0D_new,
                          bDoNotAddH, vABParityUnknown, nInputType,
                          orig_at_data ? &at_new : NULL, MAX_ATOMS,
                          &num_dimensions_new, &num_inp_bonds_new,
                          pSdfLabel, pSdfValue, lSdfId, pInpAtomFlags, err, pStrErr );

        if (num_inp_atoms_new <= 0 && !*err)
        {
            TREAT_ERR( *err, 0, "Empty structure" );
            *err = 98;
        }
        else if (orig_at_data &&
                  !num_inp_atoms_new &&
                  10 < *err && *err < 20 &&
                  orig_at_data->num_atoms > 0
                  && bMergeAllInputStructures)
        {
            *err = 0; /* end of file */
            break;
        }
        else if (num_inp_atoms_new > 0 && orig_at_data)
        {
            /*  merge pOrigDataTmp + orig_at_data => pOrigDataTmp; */
            nNumAtoms = num_inp_atoms_new + orig_at_data->num_atoms;
            nNumStereo0D = num_inp_0D_new + orig_at_data->num_stereo0D;

            if (nNumAtoms >= max_num_at) /*MAX_ATOMS ) */
            {
                TREAT_ERR( *err, 0, "Too many atoms [check  'LargeMolecules' switch]" );
                *err = 70;
                orig_at_data->num_atoms = -1;
            }
            else if (!at_old)
            {
                /* the first structure */

                orig_at_data->atom = at_new;            at_new = NULL;
                orig_at_data->num_atoms = num_inp_atoms_new; num_inp_atoms_new = 0;
                orig_at_data->stereo0D = stereo0D_new;      stereo0D_new = NULL;
                orig_at_data->num_stereo0D = num_inp_0D_new;    num_inp_0D_new = 0;
            }
            else if (orig_at_data->atom = Createinp_ATOM( nNumAtoms ))
            {
                /*  switch at_new <--> orig_at_data->at; */

                if (orig_at_data->num_atoms)
                {
                    memcpy( orig_at_data->atom, at_old, orig_at_data->num_atoms * sizeof( orig_at_data->atom[0] ) );
                    /*  adjust numbering in the newly read structure */
                    for (i = 0; i < num_inp_atoms_new; i++)
                    {
                        for (j = 0; j < at_new[i].num_bonds; j++)
                        {
                            at_new[i].neighbor[j] += orig_at_data->num_atoms;
                        }
                    }
                }
                Freeinp_ATOM( &at_old );

                /*  copy newly read structure */
                memcpy( orig_at_data->atom + orig_at_data->num_atoms,
                        at_new,
                        num_inp_atoms_new * sizeof( orig_at_data->atom[0] ) );

                /*  copy newly read 0D stereo */
                if (num_inp_0D_new > 0 && stereo0D_new)
                {
                    if (orig_at_data->stereo0D = CreateInchi_Stereo0D( nNumStereo0D ))
                    {
                        int ncopy = orig_at_data->num_stereo0D * sizeof( orig_at_data->stereo0D[0] );
                        memcpy( orig_at_data->stereo0D, stereo0D_old, ncopy );

                        /*  adjust numbering in the newly read structure */
                        for (i = 0; i < num_inp_0D_new; i++)
                        {
                            if (stereo0D_new[i].central_atom >= 0)
                                stereo0D_new[i].central_atom += orig_at_data->num_atoms;

                            for (j = 0; j < 4; j++)
                                stereo0D_new[i].neighbor[j] += orig_at_data->num_atoms;
                        }

                        FreeInchi_Stereo0D( &stereo0D_old );

                        int ncopy = num_inp_0D_new * sizeof( orig_at_data->stereo0D[0] );
                        memcpy( orig_at_data->stereo0D + orig_at_data->num_stereo0D,
                                stereo0D_new,
                                ncopy );
                    }
                    else
                    {
                        num_inp_0D_new = 0;
                        TREAT_ERR( *err, 0, "Out of RAM" );
                        *err = -1;
                    }
                }
                else
                {
                    num_inp_0D_new = 0;
                }

                /* update lengths */
                orig_at_data->num_atoms += num_inp_atoms_new;
                orig_at_data->num_stereo0D += num_inp_0D_new;
            }
            else
            {
                TREAT_ERR( *err, 0, "Out of RAM" );
                *err = -1;
            }
        }
        else if (num_inp_atoms_new > 0)
        {
            nNumAtoms += num_inp_atoms_new;
        }

        Freeinp_ATOM( &at_new );
        num_inp_atoms_new = 0;
        FreeInchi_Stereo0D( &stereo0D_new );
        num_inp_0D_new = 0;
    }
    while (!*err && bMergeAllInputStructures);

 /*
 if ( !*err )
 {
     orig_at_data->num_components =
         MarkDisconnectedComponents( orig_at_data );
     if ( orig_at_data->num_components == 0 )
     {
         TREAT_ERR (*err, 0, "No components found");
         *err = 99;
     }
     if ( orig_at_data->num_components < 0 )
     {
         TREAT_ERR (*err, 0, "Too many components");
         *err = 99;
     }
 }
 */

    if (szCoordNew)
        inchi_free( szCoordNew );

    if (at_new)
        inchi_free( at_new );

    /*
    if ( !*err )
    {
        if ( ReconcileAllCmlBondParities( orig_at_data->atom, orig_at_data->num_atoms ) )
        {
            TREAT_ERR (*err, 0, "Cannot reconcile stereobond parities");
            if (!orig_at_data->num_atoms)
            {
                *err = 1;
            }
        }
    }
    */

    if (*err)
        FreeInchi_Input( orig_at_data );

    if (*err && !( 10 < *err && *err < 20 ) && pStrErr && !pStrErr[0])
    {
        TREAT_ERR( *err, 0, "Unknown error" );  /*   <BRKPT> */
    }

    return orig_at_data ? orig_at_data->num_atoms : nNumAtoms;
}

#endif /* #ifdef TARGET_EXE_USING_API */


/****************************************************************************
 Read input InChI string and fill internal data structs
****************************************************************************/
int InchiToOrigAtom( INCHI_IOSTREAM *inp_molfile,
                     ORIG_ATOM_DATA *orig_at_data,
                     int            bMergeAllInputStructures,
                     int            bGetOrigCoord,
                     int            bDoNotAddH,
                     int            vABParityUnknown,
                     INPUT_TYPE     nInputType,
                     char           *pSdfLabel,
                     char           *pSdfValue,
                     unsigned long *lSdfId,
                     INCHI_MODE     *pInpAtomFlags,
                     int            *err,
                     char           *pStrErr )
{
    int       num_dimensions_new;
    int       num_inp_bonds_new;
    int       num_inp_atoms_new;
    inp_ATOM  *at_new = NULL;
    inp_ATOM  *at_old = NULL;
    MOL_COORD *szCoordNew = NULL;
    MOL_COORD *szCoordOld = NULL;
    int       nNumAtoms = 0;
    int       i, j, max_num_at;


    max_num_at = MAX_ATOMS;
    if (!( *pInpAtomFlags  & FLAG_SET_INP_LARGE_MOLS ))
    {
        max_num_at = NORMALLY_ALLOWED_INP_MAX_ATOMS;
    }

    if (pStrErr)
    {
        pStrErr[0] = '\0';
    }

    /*FreeOrigAtData( orig_at_data );*/
    if (lSdfId)
    {
        *lSdfId = 0LU;
    }

    do
    {
        at_old = orig_at_data ? orig_at_data->at : NULL; /*  save pointer to the previous allocation */

        szCoordOld = orig_at_data ? orig_at_data->szCoord : NULL;

        num_inp_atoms_new = InchiToInpAtom( inp_molfile,
                                            ( bGetOrigCoord && orig_at_data ) ? &szCoordNew : NULL,
                                            bDoNotAddH,
                                            vABParityUnknown,
                                            nInputType,
                                            orig_at_data ? &at_new : NULL,
                                            MAX_ATOMS,
                                            &num_dimensions_new,
                                            &num_inp_bonds_new,
                                            pSdfLabel,
                                            pSdfValue,
                                            lSdfId,
                                            pInpAtomFlags,
                                            err,
                                            pStrErr );

        if (num_inp_atoms_new <= 0 && !*err)
        {
            TREAT_ERR( *err, 0, "Empty structure" );
            *err = 98;
        }
        else if (orig_at_data && !num_inp_atoms_new &&
                  10 < *err && *err < 20 &&
                  orig_at_data->num_inp_atoms > 0 &&
                  bMergeAllInputStructures)
        {
            *err = 0; /* end of file */
            break;
        }
        else if (num_inp_atoms_new > 0 && orig_at_data)
        {
            /*  merge pOrigDataTmp + orig_at_data => pOrigDataTmp; */
            nNumAtoms = num_inp_atoms_new + orig_at_data->num_inp_atoms;
            if (nNumAtoms >= max_num_at) /*MAX_ATOMS ) */
            {
                TREAT_ERR( *err, 0, "Too many atoms [check  'LargeMolecules' switch]" );
                *err = 70;
                orig_at_data->num_inp_atoms = -1;
            }
            else if (!at_old)
            {
                /* the first structure */
                orig_at_data->at = at_new;
                orig_at_data->szCoord = szCoordNew;
                at_new = NULL;
                szCoordNew = NULL;
                orig_at_data->num_inp_atoms = num_inp_atoms_new;
                orig_at_data->num_inp_bonds = num_inp_bonds_new;
                orig_at_data->num_dimensions = num_dimensions_new;
            }
            else if (( orig_at_data->at = (inp_ATOM*) inchi_calloc( nNumAtoms, sizeof( inp_ATOM ) ) ) &&
                ( !szCoordNew || ( orig_at_data->szCoord = (MOL_COORD *) inchi_calloc( nNumAtoms, sizeof( MOL_COORD ) ) ) ))
            {
                /*  switch at_new <--> orig_at_data->at; */
                if (orig_at_data->num_inp_atoms)
                {
                    memcpy(orig_at_data->at,
                        at_old,
                        orig_at_data->num_inp_atoms * sizeof(orig_at_data->at[0]));
                    /*  adjust numbering in the newly read structure */
                    for (i = 0; i < num_inp_atoms_new; i++)
                    {
                        if (at_new) /* djb-rwth: fixing a NULL pointer dereference */
                        {
                            for (j = 0; j < at_new[i].valence; j++)
                            {
                                at_new[i].neighbor[j] += orig_at_data->num_inp_atoms;
                            }
                            at_new[i].orig_at_number += orig_at_data->num_inp_atoms; /* 12-19-2003 */
                        }
                    }
                    if (orig_at_data->szCoord && szCoordOld)
                    {
                        memcpy(orig_at_data->szCoord,
                            szCoordOld,
                            orig_at_data->num_inp_atoms * sizeof(MOL_COORD));
                    }
                }
                if (at_old)
                {
                    /* inchi_free( at_old ); */ /* djb-rwth: avoiding the use of freed memory */
                    at_old = NULL;
                }
                if (szCoordOld)
                {
                    /* inchi_free( szCoordOld ); */ /* djb-rwth: avoiding the use of freed memory */
                    szCoordOld = NULL;
                }
                /*  copy newly read structure */
                if (at_new) /* djb-rwth: fixing a NULL pointer dereference */
                    memcpy(orig_at_data->at + orig_at_data->num_inp_atoms, at_new, num_inp_atoms_new * sizeof(orig_at_data->at[0]));
                if (orig_at_data->szCoord && szCoordNew)
                {
                    memcpy(orig_at_data->szCoord + orig_at_data->num_inp_atoms,
                        szCoordNew,
                        num_inp_atoms_new * sizeof(MOL_COORD));
                }
                /*  add other things */
                orig_at_data->num_inp_atoms += num_inp_atoms_new;
                orig_at_data->num_inp_bonds += num_inp_bonds_new;
                orig_at_data->num_dimensions = inchi_max( num_dimensions_new, orig_at_data->num_dimensions );
            }
            else
            {
                TREAT_ERR( *err, 0, "Out of RAM" );
                *err = -1;
            }
        }
        else if (num_inp_atoms_new > 0)
        {
            nNumAtoms += num_inp_atoms_new;
        }
        if (at_new)
        {
            /* inchi_free( at_new ); */ /* djb-rwth: avoiding the use of freed memory */
            at_new = NULL;
        }
    }
    while (!*err && bMergeAllInputStructures);

     /*
     if ( !*err ) {
         orig_at_data->num_components =
             MarkDisconnectedComponents( orig_at_data );
         if ( orig_at_data->num_components == 0 ) {
             TREAT_ERR (*err, 0, "No components found");
             *err = 99;
         }
         if ( orig_at_data->num_components < 0 ) {
             TREAT_ERR (*err, 0, "Too many components");
             *err = 99;
         }
     }
     */

    /* djb-rwth: avoiding the use of freed memory */
    /*
    if (szCoordNew)
    {
        inchi_free( szCoordNew );
    }
    if (at_new)
    {
        inchi_free( at_new );
    }
    */

    if (!*err && orig_at_data)
    {
        if (ReconcileAllCmlBondParities( orig_at_data->at,
            orig_at_data->num_inp_atoms, 0 ))
        {
            TREAT_ERR( *err, 0, "Cannot reconcile stereobond parities" );  /* <BRKPT> */
            if (!orig_at_data->num_dimensions)
            {
                *err = 1;
            }
        }
    }

    if (*err)
    {
        FreeOrigAtData( orig_at_data );
    }

    if (*err && !( 10 < *err && *err < 20 ) && pStrErr && !pStrErr[0])
    {
        TREAT_ERR( *err, 0, "Unknown error" );  /*   <BRKPT> */
    }

    return orig_at_data ? orig_at_data->num_inp_atoms : nNumAtoms;
}


/****************************************************************************
Returns 1 if bonds (a1,a2) and (b1,b2) are the same, -1 if atoms swapped,
0 if not the same
****************************************************************************/
int bIsSameBond(int a1, int a2, int b1, int b2)
{
    if (a1 == b1 && a2 == b2) return 1;
    if (a1 == b2 && a2 == b1) return -1;
    return 0;
}


/****************************************************************************
Parse InChI and get a list of crossing bonds in z layer
Return number of frame_shift_info triples or 0
****************************************************************************/
static int GetFrameShiftInfoFrom105PlusInChI(char *sinchi,
    int *frame_shift_info,
    int max_crossing)
{
    int k, c = 0, j, aindex = 0, iunit = 0;
    const char *p, *q;

    p = strstr(sinchi, "/z");  /* must always be there */

                               /* each frame_shift_info triple(iunit, iunit_a1, iunit_a2) contains,
                               for each eligible frame-shiftable unit,
                               iunit - unit_no
                               iunit_a1, iunit_a2 - atom numbers for the senior bkbond
                               note that iunit_a1 is more senior then iunit_a2
                               */

                               /* eligible unit has Z-layer pattern
                               "range-of-numbers(number1,number2,nimbers...)"
                               >=2 bkbonds in CRU; relink may be necessary to shift frame or swap bkbond atoms?
                               OPTIONALLY DO NOT DO SWAP NOW?
                               senior bkbond and right atoms order is (number1,number2)
                               "range-of-numbers(number1.number2)"
                               1-bkbond CRU; relink may still be necessary to swap bkbond atoms so that
                               more senior atom is connected to lesser-numbered Zz
                               OPTIONALLY DO NOT DO SWAP NOW?
                               senior atoms order in bkbond is (number1,number2)
                               "range-of-numbers(number)"  1-atom CRU
                               relink is not applicable, skip it
                               */

    while (p)
    {
        int num[2] = { -1,-1 };
        p = strstr(p + 2, "(");
        if (!p)
        {
            break;
        }
        p++;
        q = p;
        j = 0;
        while ((k = (int)inchi_strtol(p, &q, 10)) && j < 2)
        {
            num[j] = k;
            j++;
            c = UCINT *q;
            if (j==1 && c == '-') /* do not consider pattern "(cap-end, cap-end)" */
            {
                goto find_next_unit;
            }
            else if (c != ')')
            {
                p = q + 1;
            }
            else
            {
                goto find_next_unit;
            }
        }
        if (j < 2)
        {
            goto find_next_unit;
        }
        frame_shift_info[3 * aindex] = iunit;
        frame_shift_info[3 * aindex + 1] = num[0];
        frame_shift_info[3 * aindex + 2] = num[1];
        aindex++;
        if (aindex >= max_crossing)
        {
            break;
        }

    find_next_unit:
        p = strstr(p, ";");
        iunit++;
    }

    return aindex;
}


/****************************************************************************
Parse AuxInfostring and get a list of original atom numbers orig[cano_num]
****************************************************************************/
int extract_orig_nums_from_auxinfo_string(char *saux, int *orig)
{
    int res = _IS_OKAY;
    int k, c = 0, cano_num = 1 /*0*/;
    const char *p, *q;

    p = strstr(saux, "/N:");  /* must always be there */
    if (!p || !p[3] || !isdigit(UCINT p[3]))
    {
        res = _IS_ERROR;
        goto exit_function;
    }

    p += 3;
    q = p;

    while ((k = inchi_strtol(p, &q, 10))) /* djb-rwth: addressing LLVM warning */
    {
        orig[cano_num++] = k/* - 1*/; /* 1-based numbers */
        if ((c = UCINT *q) && c != '/') /* djb-rwth: addressing LLVM warning */
        {
            p = q + 1;
        }
        else
        {
            break;
        }
    }

exit_function:

    return res;
}


/****************************************************************************
 (currently, this function gets only the first list of E: )
****************************************************************************/
int extract_nonstereo_eq_classes_from_auxinfo_string( char *saux,
                                                      int nat,
                                                      int *orig,
                                                      int *nclasses,
                                                      int *eclass,
                                                      int *eclass_by_origs)
{
    int res = _IS_OKAY;
    int k, c = 0, cano_num = 1, orig_num = 1;
    const char *p, *q;

    /* Note that all atom and class numbers here are 1-based */

    *nclasses = 0;
    memset(eclass, -1, ((long long)nat+1) * sizeof(int)); /* djb-rwth: cast operator added; memset_s C11/Annex K variant? */
    memset(eclass_by_origs, -1, ((long long)nat+1) * sizeof(int)); /* djb-rwth: cast operator added; memset_s C11/Annex K variant? */

    p = strstr(saux, "/E:");
    if (!p)
    {
        /* No "/E" means that all atoms are different  */
        return res;
    }

    p += 3;
    q = p;
    while ((k = (AT_NUMB)inchi_strtol(p + 1, &q, 10))) /* djb-rwth: addressing LLVM warning */
    {
        c = UCINT *q;
        if (c == '/')
        {
            break;
        }
        else if (c == ',' || c == ')')
        {
            eclass[k] = *nclasses;
            if (c == ')')
            {
                (*nclasses)++;
                q++;
                c = UCINT *q;
                if (c == '/')
                    break;
                else
                    ;
            }
            p = q;
        }
        else
        {
            return _IS_ERROR;
        }
    }
    /* NB: cano, origs start from 0 */
    for (cano_num = 1; cano_num <= nat; cano_num++)
    {
        if (eclass[cano_num] == -1) /* the atom is unique, add one more eq class for him */
        {
            (*nclasses)++;
            eclass[cano_num] = *nclasses;
        }
    }

    for (cano_num = 1; cano_num <= nat; cano_num++)
    {
        orig_num = orig[cano_num];  /* NB: cano, origs start from 0 */
        eclass_by_origs[orig_num] = eclass[cano_num];
    }

    return res;
}


/****************************************************************************
Make a copy of the context of ProcessOneStructureEx or set a a new one
****************************************************************************/
int  POSEContext_Init(POSEContext *context,
                      STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                      PINChI2 *pINChI2[INCHI_NUM], PINChI_Aux2 *pINChI_Aux2[INCHI_NUM],
                      INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file,
                      INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file,
                      ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                      long num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits)
{
    char *sz = NULL;
    int ret = _IS_OKAY, res = 0, i;

    memset(context, 0, sizeof(*context)); /* djb-rwth: memset_s C11/Annex K variant? */

    if (!sd)
    {
        memset(&context->sd, 0, sizeof(context->sd)); /* djb-rwth: memset_s C11/Annex K variant? */
    }
    else
    {
        memcpy(&context->sd, sd, sizeof(context->sd));
    }

    if (!ip)
    {
        memset(&context->ip, 0, sizeof(context->ip)); /* djb-rwth: memset_s C11/Annex K variant? */
    }
    else
    {
        memcpy(&context->ip, ip, sizeof(context->ip));
        for (i = 0; i < MAX_NUM_PATHS; i++)
        {
            if (ip->path[i])
            {
                sz = (char*)inchi_malloc((strlen(ip->path[i]) + 1) * sizeof(sz[0]));
                if (!sz)
                {
                    ret = _IS_ERROR;
                    goto exit_function;
                }
                strcpy(sz, context->ip.path[i]);
                context->ip.path[i] = sz;
            }
        }
    }

    if (strlen(szTitle))
    {
        strcpy(context->szTitle, szTitle);
    }
    else
    {
        context->szTitle[0] = '\0';
    }

    /* pINChI2, pINChI_Aux2: We do not fill/allocate elements of these structures   */
    /* assuming that NULL's are there. If not just raise an error.                  */

    context->pINChI2[0] = context->pINChI2[1] = NULL;
    if (pINChI2 && (pINChI2[0] || pINChI2[1])) /* djb-rwth: condition corrected */
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    context->pINChI_Aux2[0] = context->pINChI_Aux2[1] = NULL; 
    if (pINChI_Aux2 && (pINChI_Aux2[0] || pINChI_Aux2[1])) /* djb-rwth: condition corrected */
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    context->out_file = context->inchi_file;
    context->log_file = context->inchi_file + 1;
    context->prb_file = context->inchi_file + 2;
    /* Initialize internal for this function output streams as string buffers */
    inchi_ios_init(context->out_file, INCHI_IOS_TYPE_STRING, NULL);
    inchi_ios_init(context->log_file, INCHI_IOS_TYPE_STRING, NULL);
    inchi_ios_init(context->prb_file, INCHI_IOS_TYPE_STRING, NULL);
    context->inp_file = NULL;
    if (inp_file)
    {
        context->inp_file = inp_file;
    }

    context->orig_inp_data = &context->OrigAtData;
    context->prep_inp_data = context->PrepAtData;

    if (orig_inp_data)
    {
        memset(context->orig_inp_data, 0, sizeof(*context->orig_inp_data)); /* djb-rwth: memset_s C11/Annex K variant? */
        res = OrigAtData_Duplicate(context->orig_inp_data, orig_inp_data);
        if (res)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }
    }

    if (prep_inp_data)
    {
        memset(context->prep_inp_data, 0, 2 * sizeof(*context->prep_inp_data)); /* djb-rwth: memset_s C11/Annex K variant? */
        res = OrigAtData_Duplicate(context->prep_inp_data, prep_inp_data);
        if (res)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }
    }

    /* num_inp, strbuf, save_opt_bits */
    context->num_inp = num_inp;
    context->save_opt_bits = save_opt_bits;
    context->strbuf = &context->temp_string_container;
    if (strbuf)
    {
        res = inchi_strbuf_create_copy(context->strbuf, strbuf);
    }
    else
    {
        res = inchi_strbuf_init(context->strbuf, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT);
    }
    if (res == -1)
    {
        ret = _IS_FATAL;
        goto exit_function;
    }

exit_function:

    return ret;
}


/****************************************************************************/
void POSEContext_Free(POSEContext *context)
{
    int i;
    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (context->ip.path[i])
        {
            inchi_free((void*)context->ip.path[i]);
            /*  cast deliberately discards 'const' qualifier */
            context->ip.path[i] = NULL;
        }
    }
    FreeAllINChIArrays(context->pINChI2, context->pINChI_Aux2, context->sd.num_components);
    if (context->inp_file)
    {
        ;
    }
    else
    {
        ;
    }
    inchi_ios_close(context->out_file);
    inchi_ios_close(context->log_file);
    inchi_ios_close(context->prb_file);
    FreeOrigAtData(context->orig_inp_data);
    FreeOrigAtData(context->prep_inp_data);
    FreeOrigAtData( context->prep_inp_data+1);
    context->num_inp = 0;
    context->save_opt_bits = 0;
    inchi_strbuf_close(context->strbuf);

    return;
}


/****************************************************************************/
void POSEContext_DebugPrint(POSEContext *context)
{
    ITRACE_("\nDUMP OF POSEContext OBJECT");
    /* sd */
    ;
    /* ip */
    ;
    /* szTitle */
    ITRACE_("\n\tszTitle = %-s", context->szTitle);
    /* pINChI2, pINChI_Aux2 */
    /* inp_file, log_file, out_file, prb_file */
    if (context->inp_file)
    {
        ;
    }
    else
    {
        ;
    }
    if (context->log_file)
    {
        ;
    }
    else
    {
        ;
    }
    if (context->out_file)
    {
        ;
    }
    else
    {
        ;
    }
    if (context->prb_file)
    {
        ;
    }
    else
    {
        ;
    }
    /* orig_inp_data, prep_inp_data, */
    /* num_inp, strbuf, save_opt_bits */
    ITRACE_("\n\tnum_inp = %-ld, ", context->num_inp);
    ITRACE_("\n\tsave_opt_bits = 0x%x, ", context->save_opt_bits);
    ITRACE_("\n\tsave_opt_bits = 0x%x, ", context->save_opt_bits);
    if (context->strbuf->nUsedLength > 0)
    {
        ITRACE_("\n\tstrbuf = %-s", context->strbuf->pStr);
    }
    else
    {
        ITRACE_("\n\tstrbuf = <empty>", context->strbuf);
    }
    ITRACE_("\n");

    return;
}


/****************************************************************************/
int OAD_StructureEdits_Init(OAD_StructureEdits *ed)
{
    ed->del_side_chains = 0; /* by default, do not delete */

    ed->del_atom = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->del_atom)						goto exitf;
    if (0 != IntArray_Alloc(ed->del_atom, 2))	goto exitf;

    ed->del_bond = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->del_bond)						goto exitf;
    if (0 != IntArray_Alloc(ed->del_bond, 2))	goto exitf;

    ed->new_bond = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->new_bond)						goto exitf;
    if (0 != IntArray_Alloc(ed->new_bond, 2))	goto exitf;

    ed->mod_bond = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->mod_bond)						goto exitf;
    if (0 != IntArray_Alloc(ed->mod_bond, 12))	goto exitf;

    ed->mod_coord = (INT_ARRAY *)inchi_calloc(1, sizeof(INT_ARRAY));
    if (!ed->mod_coord)						goto exitf;
    if (0 != IntArray_Alloc(ed->mod_coord, 4))	goto exitf;


    return 0;

exitf:
    OAD_StructureEdits_Clear(ed);
    return _IS_ERROR;
}


/****************************************************************************/
void OAD_StructureEdits_Clear(OAD_StructureEdits *ed)
{
    if (ed->del_atom)
    {
        IntArray_Free(ed->del_atom);
        inchi_free(ed->del_atom);
        ed->del_atom = NULL;
    }
    if (ed->del_bond)
    {
        IntArray_Free(ed->del_bond);
        inchi_free(ed->del_bond);
        ed->del_bond = NULL;
    }
    if (ed->mod_bond)
    {
        IntArray_Free(ed->mod_bond);
        inchi_free(ed->mod_bond);
        ed->mod_bond = NULL;
    }
    if (ed->new_bond)
    {
        IntArray_Free(ed->new_bond);
        inchi_free(ed->new_bond);
        ed->new_bond = NULL;
    }
    if (ed->mod_coord)
    {
        IntArray_Free(ed->mod_coord);
        inchi_free(ed->mod_coord);
        ed->mod_coord = NULL;
    }

    return;
}


/****************************************************************************/
void OAD_StructureEdits_DebugPrint(OAD_StructureEdits *ed)
{
    ITRACE_("\n*****************************\nOAD_StructureEdits @ %-p\n*****************************", ed); 
    ITRACE_("\nDel_side_chains :\t%-d\n", ed->del_side_chains); 
    ITRACE_("Del_atom:\t%-s", ed->del_atom->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->del_atom);
    ITRACE_("Del_bond:\t%-s", ed->del_bond->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->del_bond);
    ITRACE_("New_bond:\t%-s", ed->new_bond->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->new_bond);
    ITRACE_("Mod_bond:\t%-s", ed->mod_bond->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->mod_bond);
    ITRACE_("Mod_coord:\t%-s", ed->mod_coord->used ? "" : "(empty)\n");
    IntArray_DebugPrint(ed->mod_coord);
    
}


/****************************************************************************
 Prepare CRU fold edits as suggested by the strings with preliminary
 generated interim (1.05+ flavoured) InChI and AuxInfo
****************************************************************************/
/* djb-rwth: placed as global variables to avoid function buffer issues */
int ec_opp[MAX_ATOMS],		/* equivalence classes for atoms, in order of 1-based orig nums	*/
ec_cano_opp[MAX_ATOMS],	/* equivalence classes for atoms, in order of 1-based cano nums	*/
at_stereo_mark_orig_opp[MAX_ATOMS],	/* stereo parities, in order of 1-based orig nums	*/
xc_opp[MAX_ATOMS];      /* Extended (stereo-aware) atom classes.
                            There are 'n_ec' non-stereo atom equivalence classes
                            For ec[i]=k, keep value k for no-stereo atoms while use
                            (k + neclasses)   for '-' parity
                            (k + 2*neclasses) for '+' parity                            */
int  OAD_Polymer_PrepareFoldCRUEdits( ORIG_ATOM_DATA *orig_at_data,
                                      char *sinchi_noedits, 
                                      char *saux_noedits,
                                      char *sinchi,
                                      char *saux,
                                      OAD_StructureEdits *ed)
{
    int ret = _IS_OKAY;
    int i, j, k;
    int err;
    char pStrErr[STR_ERR_LEN];
    int *orig = NULL;
    int nat = orig_at_data->num_inp_atoms;
    int neclasses = 0;		/* No of constitutional equivalence classses for the atoms		*/
    int nxclasses = 0;      /* No of extended (stereo-aware) atom classses == 3*neclasses   */
    
    int *all_bkb_orig = NULL, n_all_bkb_orig = 0;
    OAD_Polymer *p = orig_at_data->polymer;
    int nu = orig_at_data->polymer->n;

    /* Extract cano_nums-->orig_nums mapping from AuxInfo AuxInfo Main Layer */
    orig = (int*)inchi_calloc((long long)nat + 1, sizeof(int)); /* djb-rwth: cast operator added */
    if (!orig)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    ret = extract_orig_nums_from_auxinfo_string(saux, orig);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    /* Extract non-stereo eq. classes data from AuxInfo */
    ret = extract_nonstereo_eq_classes_from_auxinfo_string(saux, nat, orig, &neclasses, ec_cano_opp, ec_opp);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    if (neclasses == 0)
    {
        goto exit_function;
    }
    /* Extract stereocenter data from InChI */

    /*ret = extract_stereo_info_from_inchi_string(sinchi, nat, orig, at_stereo_mark_orig);*/
    ret = extract_stereo_info_from_inchi_string(sinchi_noedits, nat, orig, at_stereo_mark_orig_opp);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    /* Make extended stereo-aware atom classes */
    nxclasses = neclasses * 3;
    for (i = 1; i <= nat; i++)  /* orig # */
    {
        int atom_class = ec_opp[i];

        if (at_stereo_mark_orig_opp[i] == INCHI_PARITY_ODD)
        {
            atom_class += neclasses;
        }
        else if (at_stereo_mark_orig_opp[i] == INCHI_PARITY_EVEN)
        {
            atom_class += 2 * neclasses;
        }
        xc_opp[i] = atom_class;
    }
    /* Extract all backbone bonds, in all units, from InChI (z layer).
        NB: we assume that units are not 'inter-crossing' so
        any particular bkbond belongs to some unique CRU.
    */
    all_bkb_orig = (int*)inchi_calloc(2 * ((long long)orig_at_data->num_inp_bonds + 1), sizeof(int)); /* djb-rwth: cast operator added */
    if (!all_bkb_orig)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    memset(all_bkb_orig, 0, ((long long)orig_at_data->num_inp_bonds + 1) * sizeof(int)); /* djb-rwth: cast operator added; memset_s C11/Annex K variant? */
    ret = extract_all_backbone_bonds_from_inchi_string(sinchi, &n_all_bkb_orig, orig, all_bkb_orig);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    /* just for case, remove those bkbonds which are not single (alternate may be here) */
    for (k = n_all_bkb_orig - 1; k >= 0; k--)
    {
        int orig1 = all_bkb_orig[2 * k];
        int orig2 = all_bkb_orig[2 * k + 1];
        int bond_type = Inp_Atom_GetBondType(orig_at_data->at, orig1 - 1, orig2 - 1);
        if (bond_type > BOND_TYPE_SINGLE) /* not == intentionally, to keep -1 ("no bond") */
        {
            /* remove k-th bond and shift others to start of list */
            int kk;
            for (kk = k; kk < n_all_bkb_orig; kk++)
            {
                all_bkb_orig[2 * kk] = all_bkb_orig[2 * (kk + 1)];
                all_bkb_orig[2 * kk + 1] = all_bkb_orig[2 * (kk + 1) + 1];
            }
            all_bkb_orig[2 * n_all_bkb_orig] = 0;
            all_bkb_orig[2 * n_all_bkb_orig + 1] = 0;
            n_all_bkb_orig--;
        }
    }

    err = OAD_ValidatePolymerAndPseudoElementData(orig_at_data,
        POLYMERS_MODERN,
        1, /* ip->bNPZz,*/
        pStrErr,
        0 /*ip->bNoWarnings*/);
    if (err)
    {
        goto exit_function;
    }

    /* For each unit analyze a possibility of folding (i.e., removal of excess in-CRU repeats) */
    for (j = 0; j < nu; j++)
    {
        OAD_PolymerUnit* u = p->units[j];

        if (u->na < 2)
        {
            goto nextj;
        }
        if (u->nb < 2)
        {
            goto nextj;
        }
        /* this is only for bi-star CRU's */
        if (!u->cap1_is_undef)
        {
            goto nextj;
        }
        if (!u->cap2_is_undef)
        {
            goto nextj;
        }

        err = analyze_CRU_folding(orig_at_data, j,
            n_all_bkb_orig, all_bkb_orig,
            nxclasses, xc_opp,
            ed);
        if (err)
        {
            ret = inchi_max(_IS_WARNING, err);
            goto nextj;
        }

    nextj:;
    }

exit_function:
    if (orig)
    {
        inchi_free(orig);
    }
    if (all_bkb_orig)
    {
        inchi_free(all_bkb_orig);
    }

    return ret;
}



/***************************************************************************/
DiylFrag* DiylFrag_New(int na, int end1, int end2, char *s)
{
    int err = 0;

    DiylFrag *pfrag = NULL;

    pfrag = (DiylFrag *)inchi_calloc(1, sizeof(DiylFrag));
    if (NULL == pfrag)
    {
        err = 1;
        goto exit_function;
    }

    pfrag->na = na; 
    pfrag->end1 = end1;
    pfrag->end2 = end2;
    pfrag->alist = NULL;
    pfrag->xclist = NULL;

    if (na > 0 )
    {
        pfrag->alist = (int *)inchi_calloc(na, sizeof(int));
        pfrag->xclist = (int *)inchi_calloc(na, sizeof(int));
        if (!pfrag->alist || !pfrag->xclist)
        {
            err = 2;
            goto exit_function;
        }
    }

    inchi_strbuf_printf(&pfrag->sig, "%-s", s);

exit_function:
    if (err)
    {
        DiylFrag_Free(pfrag);
        inchi_free(pfrag); /* djb-rwth: addressing coverity ID #499507 */
        return NULL;
    }
    return pfrag;
}
/***************************************************************************/
void DiylFrag_Free(DiylFrag *pfrag)
{
    if (!pfrag)
    {
        return;
    }
    if (pfrag->alist)
    {
        inchi_free(pfrag->alist);
        pfrag->alist = NULL;
    }
    if (pfrag->xclist)
    {
        inchi_free(pfrag->xclist);
        pfrag->xclist = NULL;
    }
    inchi_strbuf_close(&pfrag->sig);
    return;
}
/***************************************************************************/
void DiylFrag_MakeSignature(DiylFrag *pfrag, 
                            int nxc,            /* n xclasses (molecule-wide)       */
                            int *xc,            /* xclasses (molecule-wide)         */
                            int *cnt )          /* temp storage: counts of xclasses */
{
    int i, k, nxc_frag; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    
    inchi_strbuf_printf(&pfrag->sig, "%-d,%-d,%-d{", pfrag->na, xc[pfrag->end1], xc[pfrag->end2]);
    for (i = 0; i < pfrag->na; i++)
    {
        pfrag->xclist[i] = xc[pfrag->alist[i]];
    }  
    nxc_frag = count_colors_in_sequence(pfrag->xclist, pfrag->na, nxc+1, cnt); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    for (k = 0; k < nxc; k++)
    {
        if (cnt[k] > 0)
        {
            /* (xclass:cnt)*/
            inchi_strbuf_printf(&pfrag->sig, "(%-d:%-d)", k, cnt[k]);
        }
    }

    inchi_strbuf_printf(&pfrag->sig, "}");

    return;
}
/***************************************************************************
 Compare two fragments and return 1 if they are different, 0 if equal
***************************************************************************/
int DiylFrag_Diff(DiylFrag *pfrag1, DiylFrag *pfrag2)
{
    if (pfrag1->na != pfrag2->na)
    {
        return 1;
    }
    if (pfrag1->nb != pfrag2->nb)
    {
        return 1;
    }
    if (pfrag1->sig.nUsedLength && pfrag2->sig.nUsedLength)
    {
        int cmp = strcmp(pfrag1->sig.pStr, pfrag2->sig.pStr);
        return cmp;
    }

    return 0;
}
/****************************************************************************
Debug print polymer data for a given SRU
****************************************************************************/
void DiylFrag_DebugTrace(DiylFrag *pfrag)
{
    int k, na;

    if (!pfrag)
    {
        return;
    }
    
    ITRACE_("DiylFrag @ %-p ", pfrag);
    na = pfrag->na;
    ITRACE_("\n\t%-d atoms. List of atoms and their xclasses : { ", na);
    for (k = 0; k < na - 1; k++)
    {
        ITRACE_(" %-d(%-d), ", pfrag->alist[k], pfrag->xclist[k]);
    }
    ITRACE_(" %-d(%-d) }\n", pfrag->alist[na - 1], pfrag->xclist[na - 1]);

    ITRACE_("\tend1 = %-d, end2 = %-d, nb = %-d\n", pfrag->end1, pfrag->end2, pfrag->nb);
    
    ITRACE_("\tSignature = '%-s'\n", pfrag->sig.pStr);

    return;
}


/***************************************************************************/
int analyze_CRU_folding(ORIG_ATOM_DATA *orig_at_data,
                        int iunit,
                        int n_all_bkb,
                        int *all_bkb,
                        int nxclasses, 
                        int *xc,
                        OAD_StructureEdits *ed)
{
    int ret = _IS_OKAY;
    int err, i, j, k, m, fail, a1, a2;
    int n_cuts = 0, n_frags = 0; 
    int n_frags_in_repeating_subunit = 0;
    int n_fold, n_frag_classes = 0;
    int subunit_last_atom, next_subunit_first_atom = 0;
    int *cut = NULL;        /* [ bkbond1at1, bkbond1at2,  bkbond2at1,bkbond2at2, ... ] 
                               these are (atoms of) backbone bonds which are non-cyclic and non-multiple ('breakable')      */
    DiylFrag **frag=NULL;   /* frag is divalent fragment surrounded by 'cut' bonds, so it may be a repeating CRU sub-unit   */
    int *frag_class=NULL;   /* fragments are classified, by their signatures, to produce unique labelling; 
                            if the two fragments have the same class, they have the same signature and whence are equivalent */
    int *frag_xc_counts = NULL; /* counts of xclass atoms in CRU, order of class numbers    */
    char pStrErr[STR_ERR_LEN];

    OAD_PolymerUnit *u = orig_at_data->polymer->units[iunit];
    ITRACE_("\n\n%-s\t\t%-s:%-d", "analyze_CRU_folding()", __FILE__,__LINE__);

    pStrErr[0] = '\0'; /* djb-rwth: fixing coverity ID #499611; pStrErr is a dummy parameter in this function and is never used */

    /* Reserve space for frag-specific xclass counts */
    frag_xc_counts = (int *)inchi_calloc((long long)nxclasses + 1, sizeof(int)); /* djb-rwth: cast operator added */
    if (!frag_xc_counts)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }


    /* Prepare list of cuts - backbone lying on the way from cap1 to cap2 */
    cut = (int *)inchi_calloc(2 * (long long)n_all_bkb, sizeof(int)); /* djb-rwth: cast operator added */
    if (!cut)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    OAD_PolymerUnit_DebugTrace(u);
    OAD_CollectBackboneBonds(orig_at_data,
                            u->na, u->alist,
                            u->end_atom1, u->end_atom2,
                            &(u->nbkbonds), u->bkbonds,
                            &err, pStrErr);
    if (err)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    OAD_PolymerUnit_DebugTrace(u);
    if (u->nbkbonds < 1)
    {
        goto exit_function;
    }

    /* Make 'cut' list from the bonds which are both in all_bkb and u->bkb
       (all_bkb eliminates bonds with order >1 and cyclic ones, 
       but may contain artificial cyclizing bond) 
    */
    for (i = 0; i <u->nbkbonds; i++)
    {
        a1 = u->bkbonds[i][0];
        a2 = u->bkbonds[i][1];
        for (j = 0; j < n_all_bkb; j++)
        {
            if (bIsSameBond(a1, a2, all_bkb[2 * j], all_bkb[2 * j + 1]))
            {
                cut[2 * n_cuts] = a1; /* djb-rwth: buffer overrun implicitly avoided in loop condition */
                cut[2 * n_cuts + 1] = a2;
                n_cuts++;
                break;
            }
        }
    }
    if (n_cuts < 1)
    {
        /* no valid sub-units is available */
        goto exit_function;
    }

    /* Collect fragments */
    n_frags = n_cuts + 1;
    frag = (DiylFrag**) inchi_calloc(n_frags, sizeof(DiylFrag *));
    if (!frag)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    frag_class = (int *) inchi_calloc(n_frags, sizeof(int));
    if (!frag_class)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    n_frag_classes = 0;
    for (i = 0; i < n_frags; i++)
    {
        /* Create fragment */
        int forbidden[4], novel=1;
        DiylFrag *pfrag = NULL; 

        /* Calculate and store signature of the fragment */
        /* 
            end_atom1...cut[i-1])---frag[i]---cut[i]---...end_atom2
        */
        if (i == 0)
        {
            a1 = u->end_atom1;
            forbidden[0] = u->cap1;
        }
        else
        {
            a1 = cut[2 * (i-1) + 1];            /* near end of prev cut */
            forbidden[0] = cut[2 * (i - 1) ];   /* far  end of prev cut */
        }
        forbidden[1] = a1;
        if (i==n_frags-1)
        {
            a2 = u->end_atom2;
            forbidden[2] = u->cap2;
        }
        else
        {
            a2 = cut[2 * i];                    /* near end of next cut */
            forbidden[2] = cut[2 * i + 1];      /* far  end of next cut */
        }
        forbidden[3] = a2;

        pfrag = DiylFrag_New(u->na, a1, a2, "");
        if (!pfrag)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }
        frag[i] = pfrag;

        ret = OAD_CollectReachableAtoms(orig_at_data, a1, 2, forbidden,
                                  &pfrag->na, pfrag->alist, &err, pStrErr);
        if (ret==_IS_ERROR)
        {
            goto exit_function;
        }

        DiylFrag_MakeSignature(pfrag, nxclasses, xc, frag_xc_counts); 

        novel = 1;
        for (j = 0; j < i; j++)
        {
            if ( !DiylFrag_Diff(frag[i], frag[j]) )
            {
                frag_class[i] = frag_class[j];
                novel = 0;
                break;
            }
        }
        if (novel)
        {
            frag_class[i] = n_frag_classes++;
        }

        ITRACE_("\nCANDIDATE CRU SUBUNIT %-d/%-d (CLASS #%-d)\t", i+1, n_frags, frag_class[i]);
        DiylFrag_DebugTrace(pfrag);
    }

    if (n_frag_classes == n_frags)
    {
        /* All classes are distinct ==> no repeats, folding is impossible, skip the CRU */
        goto exit_function;
    }
        
    n_frags_in_repeating_subunit = len_repeating_subsequence(frag_class, NULL, n_frags);
    if (0 == n_frags_in_repeating_subunit)
    {
        /* valid repeating pattern not found */
        goto exit_function;
    }
    n_fold = n_frags / n_frags_in_repeating_subunit;
    if (1==n_fold || (0!=n_frags%n_frags_in_repeating_subunit) )    
    {
        /* valid repeating pattern not found */
        goto exit_function;
    }
    ITRACE_("\n");

    /* {1...2}---{5...6}---{8...9} */

    ITRACE_("\n* Found %-d times foldable unit of %-d fragments\n* First repeating sub-unit formed by %-d-fragment backbone : ",
            n_fold, n_frags, n_frags_in_repeating_subunit);
    
    for (k = 0; k < n_frags_in_repeating_subunit && n_frags_in_repeating_subunit < n_frags && frag[k]; k++) /* djb-rwth: fixing a NULL pointer dereference and buffer overflow */
    {
        if (frag[k]->end1 == frag[k]->end2)
        {
            ITRACE_("-{%-d}-", frag[k]->end1, frag[k]->end2);
        }
        else
        {
            ITRACE_("-{%-d...%-d}-", frag[k]->end1, frag[k]->end2);
        }
    }

    ITRACE_("\n");
    ITRACE_("* Backbone pattern for %-d fragments that may be removed :  ", n_frags - n_frags_in_repeating_subunit);
    for (k = n_frags_in_repeating_subunit; k < n_frags; k++)
    {
        if (frag[k]->end1 == frag[k]->end2)
        {
            ITRACE_("-{%-d}-", frag[k]->end1, frag[k]->end2);
        }
        else
        {
            ITRACE_("-{%-d...%-d}-", frag[k]->end1, frag[k]->end2);
        }
    }
    ITRACE_("\n");
    
    /* Folding is possible, prepare the edits 
        Keep the least in-CRU repeating subunit 
            { frag[0] ... frag[n_frags_in_repeating_subunit-1] }
        and remove 
            { frag[n_frags_in_repeating_subunit]...frag[n_frags-1] } and all side chain attached to that  
    
        NB: which bond is modified and which is broke is important for applying these edits further!			
    */

    /* Break bond from the subunit to the next fragment and replace an original 
       bond to "right" cap with bond from the subunit "right" atom  
    */

    /*djb-rwth: the whole block had to be rewritten to fix NULL pointer dereference */
    if (n_frags_in_repeating_subunit < n_frags && frag[n_frags_in_repeating_subunit] && frag[n_frags_in_repeating_subunit - 1]) /* djb-rwth: fixing a NULL pointer dereference and buffer overflow */
    {
        subunit_last_atom        = frag[n_frags_in_repeating_subunit - 1]->end2;
        next_subunit_first_atom  = frag[n_frags_in_repeating_subunit]->end1;

        fail = 0;
        fail += IntArray_Append(ed->del_bond, subunit_last_atom);
        fail += IntArray_Append(ed->del_bond, next_subunit_first_atom);

        fail += IntArray_Append(ed->mod_bond, u->end_atom2);
        fail += IntArray_Append(ed->mod_bond, u->cap2);
        fail += IntArray_Append(ed->mod_bond, subunit_last_atom);
        fail += IntArray_Append(ed->mod_bond, u->cap2);

        if (fail)
        {
            ret = _IS_ERROR;
            goto exit_function;
        }
    }
            
    /*	Now collect all backbone atoms to be deleted (we will then delete the
    associated side chains also, but no need to reveal them at the moment)	*/

    for (k = n_frags_in_repeating_subunit; k < n_frags; k++)
    {
        if (frag[k]) /* djb-rwth: fixing a NULL pointer dereference */
        {
            for (m = 0; m < frag[k]->na; m++)
            {
                fail = IntArray_AppendIfAbsent(ed->del_atom, frag[k]->alist[m]);
                if (fail)
                {
                    ret = _IS_ERROR;
                    goto exit_function;
                }
            }
        }
    }
    /* Care on atom coordinates: as bond to cap2 changes, 
       we use coordinates of next_subunit_first_atom for cap2 
    */
    fail = 0;
    fail += IntArray_Append(ed->mod_coord, next_subunit_first_atom);
    fail += IntArray_Append(ed->mod_coord, u->cap2);
    if (fail)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

exit_function:
    if (cut)
    {
        inchi_free(cut);
    }
    if (frag)
    {
        for (i = 0; i < n_frags; i++)
        {
            DiylFrag_Free(frag[i]);
            inchi_free(frag[i]);
        }
        inchi_free(frag);
    }
    if (frag_class)
    {
        inchi_free(frag_class);
    }
    if (frag_xc_counts)
    {
        inchi_free(frag_xc_counts);
    }

    return ret;
}

/***************************************************************************
 Return number of colors ncol<=maxcol in the sequence of n colored entries 
 and counts of individiual colors
***************************************************************************/
int count_colors_in_sequence( int *color, int n, int maxcol, int *counts)
{
    int i, ncol=0;
    memset(counts, 0, maxcol * sizeof(int)); /* djb-rwth: memset_s C11/Annex K variant? */
    for (i = 0; i<n; i++) 
    {
        int colori = color[i];
        if (colori < 0) /* removed orig atom (H D etc.) */
        {
            continue;
        }
        if (0==counts[colori ])
        {
            ncol++;
        }
        counts[ color[i] ]++;
    }
    return ncol;
}


/***************************************************************************
 Find repeating starting subsequence in the sequence of n entries
 and return its length m
 each i-th entry, 0<i<m, is characterized by color[i] and optional color2[i]
***************************************************************************/
int len_repeating_subsequence(int *color, int *color2, int n)
{
    int m, k;

    if (n < 2 || !color)
    {
        return 0;
    }

    for (m = 0; m < (n + 1) / 2; m++)
    {
        for (k = m + 1; k < n; k++)
        {
            if (color[k] != color[k - m - 1]) 
            { 
                goto nextm; 
            }
            if (color2 && color2[k] != color2[k - m - 1])
            {
                goto nextm;
            }
        }
        return (m + 1);
nextm:	;
    }

    return 0;
}


/****************************************************************************
 Prepare CRU edits suggested by the string containing  preliminary generated
 interim (1.05+ flavoured) InChI and AuxInfo
****************************************************************************/
int  OAD_Polymer_PrepareFrameShiftEdits( ORIG_ATOM_DATA *orig_at_data,
                                         char *sinchi,
                                         char *saux,
                                         OAD_StructureEdits *ed)
{
    int ret = _IS_OKAY;
    int *orig = NULL, *frame_shift_info = NULL;
    int n_frame_shifts, j;
    ModSCenterInfo *scinfo = NULL;		/* 4 elements; [0]th for old_end1, [1] old_end2, [2] end1, [3] end2	*/
    
    OAD_Polymer *p = orig_at_data->polymer;
    int nu = orig_at_data->polymer->n;
    int nat = orig_at_data->num_inp_atoms;
    
    /* Extract cano_nums-->orig_nums mapping for InChI AuxInfo Main Layer */
    orig = (int *)inchi_calloc((long long)nat + 1, sizeof(int)); /* djb-rwth: cast operator added */
    if (!orig)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    ret = extract_orig_nums_from_auxinfo_string(saux, orig);
    if (ret != _IS_OKAY && ret != _IS_WARNING)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }

    scinfo = (ModSCenterInfo *)inchi_calloc(4, sizeof(scinfo[0]));
    if (!scinfo)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    

    /* Parse InChI and extract, for each 'bistar' CRU, the senior bkbond (to frame-shift brackets to its ends) */
    frame_shift_info = (int *)inchi_calloc(3 * ((long long)nu + 1), sizeof(int)); /* djb-rwth: cast operator added */
    if (!frame_shift_info)
    {
        ret = _IS_ERROR;
        goto exit_function;
    }
    n_frame_shifts = GetFrameShiftInfoFrom105PlusInChI(sinchi, frame_shift_info, nu);
    /* translate atom numbers to orig numbers */
    for (j = 0; j < n_frame_shifts; j++)
    {
        frame_shift_info[3 * j + 1] = orig[frame_shift_info[3 * j + 1]];
        frame_shift_info[3 * j + 2] = orig[frame_shift_info[3 * j + 2]];
    }

    /* Collect OAD edits */
    for (j = 0; j < n_frame_shifts; j++)
    {
        OAD_PolymerUnit *u = NULL;
        int k, iu = -1; /*int iu = frame_shift_info[3 * j];*/
        int end1, cap1, cap1_is_star, end2, cap2, cap2_is_star, old_end1, old_end2, err, fail = 0;

        end1 = frame_shift_info[3 * j + 1];
        end2 = frame_shift_info[3 * j + 2];
        
        /* Find the unit to edit (== that unit whose alist contains the new end atoms) */
        for (k = 0; k < p->n; k++)
        {
            int ak, present=0;
            
            if (NULL == p->units[k]->blist || p->units[k]->nb < 2 )
            {
                /* No crossing bonds in the unit */
                continue;
            }
            /* Find the unit to edit (== that unit whose backbone contains the new end atoms)
            for (bk = 0; bk < p->units[k]->nbkbonds; bk++ )
            {
                if ( bIsSameBond(end1, end2, p->units[k]->bkbonds[bk][0], p->units[k]->bkbonds[bk][1] ) )
                {
                    iu = k;
                    break;
                }
            }*/
            for (ak = 0; ak < p->units[k]->na; ak++ )
            {
                if (p->units[k]->alist[ak] == end1 || p->units[k]->alist[ak] == end2)
                {
                    present++;
                }
                if (present==2)
                {
                    iu = k;
                    break;
                }
            }
        }
        if (iu < 0)
        {
            /* Unit to edit unexpectedly not found, that's an error */
            ret = _IS_ERROR;
            goto exit_function;
        }

        u = p->units[iu];

        OAD_PolymerUnit_FindEndsAndCaps(u, orig_at_data,
                                        &old_end1, &cap1, &cap1_is_star,
                                        &old_end2, &cap2, &cap2_is_star,
                                        &err, NULL);

        if (!err && cap1_is_star && cap2_is_star && end1 && end2 && cap1 && cap2)
        {
            /* find old CRU ends */
            if (cap1 == u->blist[0])		old_end1 = u->blist[1];
            else if (cap1 == u->blist[1])	old_end1 = u->blist[0];
            else if (cap1 == u->blist[2])	old_end1 = u->blist[3];
            else if (cap1 == u->blist[3])	old_end1 = u->blist[2];
            else /* something wrong */
                continue;
            if (cap2 == u->blist[0])		old_end2 = u->blist[1];
            else if (cap2 == u->blist[1])	old_end2 = u->blist[0];
            else if (cap2 == u->blist[2])	old_end2 = u->blist[3];
            else if (cap2 == u->blist[3])	old_end2 = u->blist[2];
            else /* something wrong */
                continue;

            if (!old_end1 || !old_end2 || old_end1 == old_end2)
            {
                continue;
            }
            if (bIsSameBond(old_end1, cap1, end1, cap1) && bIsSameBond(old_end2, cap2, end2, cap2))
            {
                continue; /* ignore swaps for now */
            }

            /*	If applicable, collect bonds to modify */
            
            /* Check if atoms involved in modifications are stereocenters (needs additional care) */
            ModSCenter_Init(&scinfo[0], orig_at_data->at, old_end1 - 1);
            ModSCenter_Init(&scinfo[1], orig_at_data->at, old_end2 - 1);
            ModSCenter_Init(&scinfo[2], orig_at_data->at, end1 - 1);
            ModSCenter_Init(&scinfo[3], orig_at_data->at, end2 - 1);

            /* djb-rwth: removing redundant code */
            if (!bIsSameBond(old_end1, cap1, end1, cap1))
            {
                /* Modify bond: (old_end1-cap1) --> (end1-cap1) */
                fail = 0;
                fail += IntArray_Append(ed->mod_bond, old_end1);
                fail += IntArray_Append(ed->mod_bond, cap1);
                fail += IntArray_Append(ed->mod_bond, end1);
                fail += IntArray_Append(ed->mod_bond, cap1);
                if (fail)
                {
                    ret = _IS_ERROR;
                    goto exit_function;
                }
                ModSCenter_DelFrom(&scinfo[0], cap1 - 1);
                ModSCenter_AddTo(&scinfo[2], cap1-1 );
            }
            if (!bIsSameBond(old_end2, cap2, end2, cap2))
            {
                /* Modify bond: (old_end2-cap2) --> (end2-cap2) */
                fail = 0;
                fail += IntArray_Append(ed->mod_bond, old_end2);
                fail += IntArray_Append(ed->mod_bond, cap2);
                fail += IntArray_Append(ed->mod_bond, end2);
                fail += IntArray_Append(ed->mod_bond, cap2);
                if (fail)
                {
                    ret = _IS_ERROR;
                    goto exit_function;
                }
                ModSCenter_DelFrom(&scinfo[1], cap2 - 1);
                ModSCenter_AddTo(&scinfo[3], cap2 - 1);
            }
            /* Modify bond: (end1-end2) --> (old_end1-old_end2) */
            fail = 0;
            fail += IntArray_Append(ed->mod_bond, end1);
            fail += IntArray_Append(ed->mod_bond, end2);
            fail += IntArray_Append(ed->mod_bond, old_end1);
            fail += IntArray_Append(ed->mod_bond, old_end2);
            if (fail)
            {
                ret = _IS_ERROR;
                goto exit_function;
            }
            ModSCenter_DelFrom(&scinfo[2], end2 - 1);
            ModSCenter_DelFrom(&scinfo[3], end1 - 1);
            ModSCenter_AddTo(&scinfo[0], old_end2 - 1);
            ModSCenter_AddTo(&scinfo[1], old_end1 - 1);

        }

        /* djb-rwth: n_flip and ModSCenter_IsChanged function completely redundant? -- discussion required */
        if (orig_at_data->num_dimensions)
        {
            /* Check if we must flip stereocenter configuration */
            /* (ignore errrors signaled by returning -1)		*/
            int n_flip = 0;
            if (0 < ModSCenter_IsChanged(&scinfo[0], orig_at_data->at))
            {
                n_flip++;
            }
            if (0 < ModSCenter_IsChanged(&scinfo[1], orig_at_data->at))
            {
                n_flip++;
            }
            if (0 < ModSCenter_IsChanged(&scinfo[2], orig_at_data->at))
            {
                n_flip++;
            }
            if (0 < ModSCenter_IsChanged(&scinfo[3], orig_at_data->at))
            {
                n_flip++;
            }
            n_flip = 1;
        }
    }

exit_function:
    if (orig)
    {
        inchi_free(orig);
    }
    if (frame_shift_info)
    {
        inchi_free(frame_shift_info);
    }
    if (scinfo)
    {
        inchi_free(scinfo);
    }

    return ret;
}

/****************************************************************************
 Initialize modifiable stereo center
****************************************************************************/
void ModSCenter_Init(ModSCenterInfo *scinfo, inp_ATOM *at, int iatom)
{
    int i;
    scinfo->num = iatom;
    scinfo->valence = at[iatom].valence;
    scinfo->n_stereo = NDefStereoBonds(at, iatom, 1); /* , bOnlyPointedEndMatters=1 */
    for (i = 0; i < scinfo->valence; i++)
    {
        scinfo->nbr[i] = at[iatom].neighbor[i];
        scinfo->new_nbr[i] = scinfo->nbr[i];
    }

    return;
}
/****************************************************************************/
int NDefStereoBonds(inp_ATOM *at, int iatom, int bOnlyPointedEndMatters)
{
    int i, n_stereo = 0;
    int stereo_value, stereo_type;
    for (i = 0; i < at[iatom].valence; i++)
    {
        stereo_value = at[iatom].bond_stereo[i];
        if (bOnlyPointedEndMatters)
        {
            /* establish the stereo considering only the pointed end of stereo bond */
            stereo_type = stereo_value;
        }
        else
        {
            stereo_type = abs(stereo_value);
        }
        if (stereo_type == STEREO_SNGL_UP || stereo_type == STEREO_SNGL_DOWN)
        {
            n_stereo++;
        }
    }
    return n_stereo;
}


/****************************************************************************
 Add atom to modifiable stereo center
****************************************************************************/
void ModSCenter_AddTo(ModSCenterInfo *scinfo, int iadd)
{
    if (!is_in_the_ilist(scinfo->new_nbr, iadd, scinfo->valence))
    {
        scinfo->new_nbr[scinfo->valence] = iadd;
        scinfo->valence++;
    }
    return;
}
/****************************************************************************
 Delete atom from modifiable stereo center
****************************************************************************/
void ModSCenter_DelFrom(ModSCenterInfo *scinfo, int idel)
{
    int i, j;
    for (i = 0; i < scinfo->valence; i++)
    {
        if (scinfo->nbr[i]==idel )
        {
            for (j=i+1; j < scinfo->valence; j++)
            {
                scinfo->new_nbr[j-1] = scinfo->new_nbr[j];
            }
            scinfo->valence--;
            return;
        }
    }
    return;
}
/****************************************************************************
 Check if stereo configuration of modifiable stereo center changed
****************************************************************************/
/* djb-rwth: n_flip and ModSCenter_IsChanged function completely redundant? -- discussion required */
int ModSCenter_IsChanged(ModSCenterInfo *scinfo, inp_ATOM *at)
{
    int i, ns, base1=-1, base2=-1, new_base2=-1, n_changed=0;
    double a[3], b[3], new_b[3], z[3], new_z[3], zz; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    if (scinfo->n_stereo < 1)
    {
        return 0;
    }
    if (scinfo->valence != at[scinfo->num].valence )
    {
        return -1; /* something went wrong */
    }
    iisort(scinfo->nbr, scinfo->valence);
    iisort(scinfo->new_nbr, scinfo->valence);
    /* Find the kept stereo base atom */
    for (i = 0; i < at[scinfo->num].valence; i++)
    {
        if ( is_in_the_ilist(scinfo->nbr, scinfo->new_nbr[i], scinfo->valence) )
        {
            ns = NDefStereoBonds(at, scinfo->new_nbr[i], 0); /* bOnlyPointedEndMatters=0 */
            if (ns==0)
            {
                base1 = scinfo->new_nbr[i];
                break;
            }
        }
    }
    if (base1==-1)
    {
        return -1; /* something went wrong */
    }
    /* Find the newly appeared stereo base atom */
    for (i = 0; i < at[scinfo->num].valence; i++)
    {
        /*!!! TUT NADO NE TAK
         base2 tot, kogo net v new_nbr
         new_base2 - tot, kogo net v  nbr
        */
        if ( !is_in_the_ilist(scinfo->nbr, scinfo->new_nbr[i], scinfo->valence))
        {
            ns = NDefStereoBonds(at, scinfo->nbr[i], 0);
            if (ns == 0)
            {
                new_base2 = scinfo->new_nbr[i];
                base2 = scinfo->nbr[i];
                n_changed++;
            }
        }
    }
    if (n_changed > 1 || new_base2 == -1 || base2 == -1)
    {
        return -1; /* something went wrong */
    }
    a[0] = at[base1].x - at[scinfo->num].x; a[1] = at[base1].y - at[scinfo->num].y; a[2] = at[base1].z - at[scinfo->num].z;
    b[0] = at[base2].x - at[scinfo->num].x; b[1] = at[base2].y - at[scinfo->num].y; b[2] = at[base2].z - at[scinfo->num].z;
    new_b[0] = at[new_base2].x - at[scinfo->num].x; new_b[1] = at[new_base2].y - at[scinfo->num].y; new_b[2] = at[new_base2].z - at[scinfo->num].z;

    cross_prod3(a, b, z);
    cross_prod3(a, new_b, new_z);
    zz = dot_prod3(z, new_z); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

    return -1;
}
