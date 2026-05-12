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
    SortAndPrintINChI and misc. processing (display, tests, error treatment, etc.)

*/

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>

#include "mode.h"

#ifndef COMPILE_ANSI_ONLY
#include <conio.h>
#endif

#include "ichimain.h"
#include "ichi_io.h"
#include "mol_fmt.h"
#include "ichicant.h"
#include "inchi_api.h"
#include "readinch.h"
#ifdef TARGET_LIB_FOR_WINCHI
#include "../../../IChI_lib/src/ichi_lib.h"
#include "inchi_api.h"
#else
#include "inchi_gui.h"
#endif
#include "readinch.h"

#include "bcf_s.h"

#ifdef TARGET_LIB_FOR_WINCHI

void( *FWPRINT )
( const char * format, va_list argptr ) = NULL;
void( *FWPUSH )
( const char *s );
void( *DRAWDATA )
( struct DrawData * pDrawData ) = NULL;
int( *DRAWDATA_EXISTS )
( int nComponent, int nType, int bReconnected ) = NULL;
struct DrawData * ( *GET_DRAWDATA )
    ( int nComponent, int nType, int bReconnected ) = NULL;
#endif

/* Local */
int GetProcessingWarningsOneInChI( INChI *pINChI,
                                   INP_ATOM_DATA *inp_norm_data,
                                   char *pStrErrStruct,
                                   int bNoWarnings);



/****************************************************************************
 Main InChI serialization procedure
****************************************************************************/
int SortAndPrintINChI( CANON_GLOBALS            *pCG,
                       INCHI_IOSTREAM           *out_file,
                       INCHI_IOS_STRING         *strbuf,
                       INCHI_IOSTREAM           *log_file,
                       INPUT_PARMS              *ip,
                       ORIG_ATOM_DATA           *orig_inp_data,
                       ORIG_ATOM_DATA           *prep_inp_data,
                       COMP_ATOM_DATA           composite_norm_data[INCHI_NUM][TAUT_NUM + 1],
                       ORIG_STRUCT              *pOrigStruct,
                       int                      num_components[INCHI_NUM],
                       int                      num_non_taut[INCHI_NUM],
                       int                      num_taut[INCHI_NUM],
                       INCHI_MODE               bTautFlags[INCHI_NUM],
                       INCHI_MODE               bTautFlagsDone[INCHI_NUM],
                       NORM_CANON_FLAGS         *pncFlags,
                       long                     num_inp,
                       PINChI2                  *pINChI[INCHI_NUM],
                       PINChI_Aux2              *pINChI_Aux[INCHI_NUM],
                       int                      *pSortPrintINChIFlags,
                       unsigned char            save_opt_bits )
{
    INCHI_SORT *pINChISort[INCHI_NUM][TAUT_NUM];
    int j, i, k, k1, ret, ret2, iINChI, max_num_components; /* djb-rwth: ignoring LLVM warning: variable used */
    int INCHI_basic_or_INCHI_reconnected;
    /* djb-rwth: removing redundant variables */
    int bDisconnectedCoord = ( 0 != ( bTautFlagsDone[0] & TG_FLAG_DISCONNECT_COORD_DONE ) );
    int bINChIOutputOptions0, bCurOption, bINChIOutputOptionsCur, bEmbedReconnected;
    static const char szAnnHdr[] = "InChI ANNOTATED CONTENTS";
#ifdef TARGET_LIB_FOR_WINCHI
    int ikflag = 0;
    out_file->type = INCHI_IOS_TYPE_STRING;
#endif

    /*
        Note:

        pINChI[INCHI_BAS] refers to either disconnected or original structure;
            num_components[INCHI_BAS] > 0 if there was input structure
        pINChI[INCHI_REC] refers to the reconnected structure,
            and only if the input structure has been disconnected, that is,
            num_components[INCHI_REC] > 0
    */

    ret = 1;

    for (i = 0; i < INCHI_NUM; i++)
    {
        for (k = 0; k < TAUT_NUM; k++)
        {
            bTautFlags[i] |= pncFlags->bTautFlags[i][k];
            bTautFlagsDone[i] |= pncFlags->bTautFlagsDone[i][k];
        }
    }

    /* djb-rwth: removing redundant code */

    max_num_components = 0;
    for (j = 0; j < INCHI_NUM; j++)
    {
        if (max_num_components < num_components[j])
        {
            max_num_components = num_components[j];
        }
    }
    if (max_num_components <= 0)
    {
        max_num_components = 1;
    }

    for (j = 0, i = 0; j < INCHI_NUM; j++)
    {
        if (num_components[j])
        {
            for (k1 = 0; k1 < TAUT_NUM; k1++)
            {
                pINChISort[j][k1] =
                    (INCHI_SORT *) inchi_calloc( max_num_components,
                                               sizeof( pINChISort[0][0][0] ) );
                i += !pINChISort[j][k1]; /* number of failed allocatons */
            }
        }
        else
        {
            for (k1 = 0; k1 < TAUT_NUM; k1++)
            {
                pINChISort[j][k1] = NULL; /* keep BC happy */
            }
        }
    }

    if (i)
    {
        ret = CT_OUT_OF_RAM;
        goto exit_function;
    }

    for (j = 0; j < INCHI_NUM; j++)
    {
        if (!num_components[j])
        {
            continue;
        }
        iINChI = j;

#if ( OUTPUT_CONNECTED_METAL_ONLY == 1 ) /* test: output connected as the only one INChI */
        if (INCHI_BAS == j && num_components[INCHI_REC])
        {
            j = INCHI_REC;
        }
#endif

        /* j = INCHI_BAS; <- for debug only */
        /* for only normal or disconnected coord compounds */
        /* (j=0=INCHI_BAS => normal or disconnected, j=1=INCHI_REC => reconnected */

        for (k1 = 0; k1 < TAUT_NUM; k1++)
        {
            for (i = 0; i < num_components[j]; i++)
            {
                for (k = 0; k < TAUT_NUM; k++)
                {
                    pINChISort[j][k1][i].pINChI[k] = pINChI[j][i][k];
                    pINChISort[j][k1][i].pINChI_Aux[k] = pINChI_Aux[j][i][k];
                }
                pINChISort[j][k1][i].ord_number = i;
            }
        }

        /* Sort component INChIs */
        for (k1 = 0; k1 < TAUT_NUM; k1++)
        {
            switch (k1)
            {
                case TAUT_NON:
                    qsort( pINChISort[j][k1],
                           num_components[j],
                           sizeof( pINChISort[0][0][0] ),
                           CompINChINonTaut2 );
                    break;
                case TAUT_YES:
                    qsort( pINChISort[j][k1],
                           num_components[j],
                           sizeof( pINChISort[0][0][0] ),
                           CompINChITaut2 );
                    break;
            }
        }

#ifndef COMPILE_ANSI_ONLY
        /* Find equivalent and wINChI display order;
        use requested in ip->bCompareComponents comparison */

        ret = SaveEquComponentsInfoAndSortOrder( iINChI, pINChISort[j],
                                                  num_components,
                                                  orig_inp_data, prep_inp_data,
#if ( FIX_DALKE_BUGS == 1 )
                                                  composite_norm_data ? composite_norm_data[j] : NULL,
#else
                                                  composite_norm_data[j],
#endif
                                                  ip->bCompareComponents );
        if (RETURNED_ERROR( ret ))
        {
            ret = 0;
            goto exit_function;
        }
        else
        {
            ret = 1;
        }
#endif
    } /* j */

    if (!( ip->bINChIOutputOptions & INCHI_OUT_PRINT_OPTIONS ))
    {
        /* Prepare InChI from the structures obtained by
           reversing InChI for returning to the caller */
        for (j = 0; j < INCHI_NUM; j++)
        {
            if (!num_components[j])
            {
                continue;
            }

            /* pINChI[iINCHI][iComponent][bTaut] */
            /* j  = disconnected/connected */
            /* k1 = sort order for Mobile or Fixed H */

            k1 = TAUT_YES; /* in Mobile H order */
            /* store components in Mobile H order */

            for (i = 0; i < num_components[j] && i < max_num_components; i++) /* djb-rwth: fixing undefined index value / buffer overflow */
            {

                if (pINChISort[j][k1][i].pINChI[TAUT_NON] &&
                    !pINChISort[j][k1][i].pINChI[TAUT_YES])
                {
                    /* make sure Mobile-H is always present */
                    for (k = 0; k < TAUT_NUM; k++)
                    {
                        pINChI[j][i][k] = pINChISort[j][k1][i].pINChI[ALT_TAUT( k )];
                        pINChI_Aux[j][i][k] = pINChISort[j][k1][i].pINChI_Aux[ALT_TAUT( k )];
                    }
                }
                else
                {
                    for (k = 0; k < TAUT_NUM; k++)
                    {
                        pINChI[j][i][k] = pINChISort[j][k1][i].pINChI[k];
                        pINChI_Aux[j][i][k] = pINChISort[j][k1][i].pINChI_Aux[k];
                    }
                }
            }
        }
    }
    else
    {
        /* Print InChI string(s) */
        bINChIOutputOptions0 = ip->bINChIOutputOptions & ~INCHI_OUT_PRINT_OPTIONS;
        bEmbedReconnected = ip->bINChIOutputOptions & INCHI_OUT_EMBED_REC;

        for (i = 1; i <= 2; i++)
        {
            bCurOption = INCHI_OUT_PLAIN_TEXT;
            if (i == 2)
            {
                bCurOption = INCHI_OUT_PLAIN_TEXT_COMMENTS;
                /* continue; */
            }
            if (!( ip->bINChIOutputOptions & bCurOption ))
            {
                continue;
            }

            bINChIOutputOptionsCur = bINChIOutputOptions0 | bCurOption;
            if (i == 1)
            {
                /* output INChI */
                bINChIOutputOptionsCur |= bEmbedReconnected;
            }
            else if (i == 2)
            {
                /* output annotation */
                inchi_ios_print( out_file, "\n==== %s ====\n", szAnnHdr );
                ITRACE_( "\n==== %s ====\n", szAnnHdr );
                bINChIOutputOptionsCur |= bEmbedReconnected;
                bINChIOutputOptionsCur &= ~INCHI_OUT_TABBED_OUTPUT;
            }
            else
            {
                continue;
            }

#if 0 /*^^^#ifdef TARGET_LIB_FOR_WINCHI*/
            inchi_ios_free_str( out_file );
#endif

            INCHI_basic_or_INCHI_reconnected = INCHI_BAS;

            ret2 = OutputINChI2( pCG,
                                 strbuf,
                                 pINChISort,
                                 INCHI_basic_or_INCHI_reconnected,
                                 orig_inp_data,
                                 pOrigStruct,
                                 ip,
                                 bDisconnectedCoord,
                                 OUT_TN,
                                 bINChIOutputOptionsCur,
                                 num_components,
                                 num_non_taut,
                                 num_taut,
                                 out_file,
                                 log_file,
                                 num_inp,
                                 pSortPrintINChIFlags,
                                 save_opt_bits );

            ret &= ret2;

#ifdef TARGET_LIB_FOR_WINCHI
            /* always calculate InChIKey */
            winchi_calc_inchikey( ret, &ikflag, ip, out_file, log_file );
            /*push_to_winchi_text_window( out_file );
            inchi_ios_flush( out_file );*/

#endif

            if (!ret)
            {
                break;
            }
        } /* i */
    }


exit_function:

    for (j = 0; j < INCHI_NUM; j++)
    {
        for (k1 = 0; k1 < TAUT_NUM; k1++) /* djb-rwth: removing redundant code */
        {
            if (pINChISort[j][k1])
            {
                inchi_free( pINChISort[j][k1] );
            }
        }
    }


    ret = ret ? 0 : _IS_FATAL;

    return ret;
}


/****************************************************************************
 Specific to winchi-1 calculation of InChIKey
****************************************************************************/
void winchi_calc_inchikey( int            ret,
                           int            *ikflag,
                           INPUT_PARMS    *ip,
                           INCHI_IOSTREAM *out_file,
                           INCHI_IOSTREAM *log_file )
{
    char ik_string[256];    /* Resulting InChIKey string */
    int ik_ret = 0;           /* InChIKey-calc result code */
    int xhash1 = 0, xhash2 = 0;
    char szXtra1[256], szXtra2[256];
    size_t slen = out_file->s.nUsedLength;
    char *buf = NULL;

    ( *ikflag )++;
    if (*ikflag != 1)
    {
        return;
    }

    if (ret == 0)
    {
        inchi_ios_flush( out_file );
        return;
    }

    extract_inchi_substring( &buf, out_file->s.pStr, slen );

    /* Calculate and print InChIKey */
    if (NULL != buf)
    {
        xhash1 = xhash2 = 0;
        if (( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1 ) ||
            ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2 ))
        {
            xhash1 = 1;
        }
        if (( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA2 ) ||
            ( ip->bCalcInChIHash == INCHIHASH_KEY_XTRA1_XTRA2 ))
        {
            xhash2 = 1;
        }
        ik_ret = GetINCHIKeyFromINCHI( buf, xhash1, xhash2, ik_string, szXtra1, szXtra2 );
        inchi_free( buf );
    }
    else
    {
        ik_ret = 3;
    }

    if (ik_ret == INCHIKEY_OK)
    {
        /* NB: correctly treat tabbed output with InChIKey & hash extensions */
        char csep = '\n';
        if (ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT)
        {
            csep = '\t';
        }
        inchi_ios_print( out_file, "InChIKey=%-s", ik_string );
        if (xhash1)
        {
            inchi_ios_print( out_file, "%cXHash1=%-s", csep, szXtra1 );
        }
        if (xhash2)
        {
            inchi_ios_print( out_file, "%cXHash2=%-s", csep, szXtra2 );
        }
        inchi_ios_print( out_file, "\n" );
    }
    else
    {
        inchi_ios_print( log_file, "Warning (Could not compute InChIKey)\n" ); /*: ", num_inp);*/
    }
    /*inchi_ios_flush(out_file);*/

    return;
}


#ifndef COMPILE_ANSI_ONLY /* { */


/****************************************************************************/
int SaveEquComponentsInfoAndSortOrder( int             iINChI,
                                        INCHI_SORT      *pINChISort[TAUT_NUM],
                                        int             *num_components,
                                        ORIG_ATOM_DATA  *orig_inp_data,
                                        ORIG_ATOM_DATA  *prep_inp_data,
                                        COMP_ATOM_DATA  composite_norm_data[TAUT_NUM + 1],
                                        int             bCompareComponents )
{
    int nRet = 0, i, k, nNumDeleted;
    int bCompareIsotopic, bCompareTaut, bCompareAlt;
    ORIG_ATOM_DATA *inp_data = NULL;

    /* Equivalent components and sorting order                  */
    /* bCompareComponents: bit = 1 => compare                   */
    /*                     bit = 2 => compare non-isotopic      */
    /*                     bit = 4 => compare non-tautomeric    */

    if (num_components[iINChI] <= 1)
    {
        return 0;
    }

#ifdef TARGET_LIB_FOR_WINCHI
    if (!DRAWDATA)
        return 0;
#endif

    if (!( bCompareComponents & CMP_COMPONENTS ))
    {
        return 0;
    }

    bCompareIsotopic = !( bCompareComponents & CMP_COMPONENTS_NONISO );

    bCompareTaut = ( bCompareComponents & CMP_COMPONENTS_NONTAUT )
        ? TAUT_NON
        : TAUT_YES;

    bCompareAlt = ALT_TAUT( bCompareTaut );

    if (num_components[iINChI] > 1)
    {
        if (prep_inp_data[iINChI].bSavedInINCHI_LIB[iINChI] &&
             prep_inp_data[iINChI].bPreprocessed[iINChI])
        {
            inp_data = prep_inp_data + iINChI;
        }
        else if (orig_inp_data->bSavedInINCHI_LIB[iINChI] &&
                  !orig_inp_data->bPreprocessed[iINChI])
        {
            inp_data = orig_inp_data;
        }
        else
        {
            inp_data = NULL;
        }

        if (inp_data && !inp_data->nEquLabels &&
             !prep_inp_data[iINChI].nSortedOrder)
        {
            int i1, i2, nSet;
            AT_NUMB nAtNo;
            AT_NUMB nNumAtoms = (AT_NUMB) inp_data->num_inp_atoms;

            if (( prep_inp_data[iINChI].nSortedOrder =
                (AT_NUMB *) inchi_calloc( (long long)num_components[iINChI] + 1,
                    sizeof( prep_inp_data[0].nSortedOrder[0] ) ) )) /* djb-rwth: cast operator added */
            {
                inp_data->nNumEquSets = 0;

                for (i1 = 0, nSet = 0; i1 < num_components[iINChI]; i1 = i2)
                {
                    nNumDeleted =
                        ( pINChISort[bCompareTaut][i1].pINChI[bCompareTaut] &&
                         pINChISort[bCompareTaut][i1].pINChI[bCompareTaut]->bDeleted );

                    for (i2 = i1 + 1; i2 < num_components[iINChI]; i2++)
                    {
                        /* isotopic/non-isotopic comparison does not separate equivalent components */
                        if (CompINChI2( pINChISort[bCompareTaut] + i1,
                            pINChISort[bCompareTaut] + i2,
                            bCompareTaut, bCompareIsotopic ))
                        {
                            break;
                        }
                        else
                        {
                            nNumDeleted +=
                                ( pINChISort[bCompareTaut][i2].pINChI[bCompareTaut] &&
                                 pINChISort[bCompareTaut][i2].pINChI[bCompareTaut]->bDeleted );
                        }
                    }

                    if (i2 - i1 - nNumDeleted > 1)
                    {
                        if (inp_data->nEquLabels ||
                            ( inp_data->nEquLabels =
                            (AT_NUMB *) inchi_calloc( (long long)inp_data->num_inp_atoms + 1,
                                sizeof( inp_data->nEquLabels[0] ) ) )) /* djb-rwth: cast operator added */
                        {
                            nSet++;
                                /* found i2-i1 equivalent components && */
                                /* memory has been allocated */
                            for (i = i1; i < i2; i++)
                            {
                                INChI_Aux *pINChI_Aux;
                                int aux_test;

                                if (pINChISort[bCompareTaut][i].pINChI[bCompareTaut] &&
                                     pINChISort[bCompareTaut][i].pINChI[bCompareTaut]->bDeleted)
                                    continue;

                                aux_test =
                                    ( pINChISort[bCompareTaut][i].pINChI_Aux[bCompareTaut] &&
                                      pINChISort[bCompareTaut][i].pINChI_Aux[bCompareTaut]->nNumberOfAtoms );

                                pINChI_Aux =
                                    aux_test ? pINChISort[bCompareTaut][i].pINChI_Aux[bCompareTaut]
                                    : ( pINChISort[bCompareTaut][i].pINChI_Aux[bCompareAlt] &&
                                        pINChISort[bCompareTaut][i].pINChI_Aux[bCompareAlt]->nNumberOfAtoms )
                                    ? pINChISort[bCompareTaut][i].pINChI_Aux[bCompareAlt]
                                    : (INChI_Aux *) NULL;

                                if (pINChI_Aux && pINChI_Aux->nOrigAtNosInCanonOrd)
                                {
                                    for (k = 0; k < pINChI_Aux->nNumberOfAtoms; k++)
                                    {
                                        if (( nAtNo = pINChI_Aux->nOrigAtNosInCanonOrd[k] ) &&
                                              nAtNo <= nNumAtoms)
                                        {
                                            inp_data->nEquLabels[nAtNo - 1] = nSet;
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            return CT_OUT_OF_RAM;
                        }
                    }
                }

                nRet |= nSet ? 1 : 0;
            }
            else
            {
                return CT_OUT_OF_RAM;
            }

            inp_data->nNumEquSets = nSet;
            /* output order */
            prep_inp_data[iINChI].nSortedOrder[0] = 0;
            for (i1 = 0; i1 < num_components[iINChI]; i1++)
            {
                prep_inp_data[iINChI].nSortedOrder[i1 + 1] =
                    pINChISort[TAUT_YES][i1].ord_number + 1;
            }

#ifdef TARGET_LIB_FOR_WINCHI /* { */
            if (DRAWDATA && GET_DRAWDATA &&
                 inp_data->nNumEquSets > 0 &&
                 inp_data->nEquLabels)
            {
                int    nType = inp_data->bPreprocessed[iINChI]
                    ? COMPONENT_ORIGINAL_PREPROCESSED
                    : COMPONENT_ORIGINAL;

                struct DrawData *pDrawData = GET_DRAWDATA( 0, nType, iINChI );
                if (pDrawData &&
                      pDrawData->pWindowData &&
                     !pDrawData->pWindowData->nEquLabels)
                {
                    /* copy equivalence data from inp_data */
                    /* to pDrawData->pWindowData           */
                    if (inp_data->nEquLabels &&
                        ( pDrawData->pWindowData->nEquLabels =
                        (AT_NUMB *) inchi_calloc( inp_data->num_inp_atoms,
                            sizeof( inp_data->nEquLabels[0] ) ) ))
                    {
                        memcpy( pDrawData->pWindowData->nEquLabels,
                                inp_data->nEquLabels,
                                inp_data->num_inp_atoms *
                                            sizeof( inp_data->nEquLabels[0] ) );
                        pDrawData->pWindowData->nNumEquSets =
                            inp_data->nNumEquSets;
                        pDrawData->pWindowData->nCurEquLabel = 0;
                    }
                }
            }

#endif  /* } TARGET_LIB_FOR_WINCHI */
        }
    }

    return nRet;
}


/****************************************************************************/
int DisplayTheWholeCompositeStructure( struct tagCANON_GLOBALS  *pCG,
                                       struct tagINCHI_CLOCK    *ic,
                                       INPUT_PARMS              *ip,
                                       struct tagStructData     *sd,
                                       long                     num_inp,
                                       int                      iINChI,
                                       PINChI2                  *pINChI2,
                                       PINChI_Aux2              *pINChI_Aux2,
                                       ORIG_ATOM_DATA           *orig_inp_data,
                                       ORIG_ATOM_DATA           *prep_inp_data,
                                       COMP_ATOM_DATA           composite_norm_data[TAUT_NUM + 1] )
{
    ORIG_ATOM_DATA *inp_data = NULL;
    int jj, j, k, err = 0, nNumIntermediateTaut = 0, bDisplayTaut;
    char szTitle[256];
    int nNumTautComponents, m;
    int bCompareIsotopic = !( ip->bCompareComponents & CMP_COMPONENTS_NONISO );
    int bCompareTaut = ( ip->bCompareComponents & CMP_COMPONENTS_NONTAUT ) ? TAUT_NON : TAUT_YES;

    if (ip->bCompareComponents & CMP_COMPONENTS)
    {
        if (prep_inp_data[iINChI].bSavedInINCHI_LIB[iINChI] && prep_inp_data[iINChI].bPreprocessed[iINChI])
        {
            inp_data = prep_inp_data + iINChI;
        }
        else if (orig_inp_data->bSavedInINCHI_LIB[iINChI] && !orig_inp_data->bPreprocessed[iINChI])
        {
            inp_data = orig_inp_data;
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

    for ( jj = 0;
          ip->bDisplayCompositeResults && !sd->bUserQuitComponentDisplay && jj <= TAUT_INI;
          jj++)
    {
        j = ( jj == 0 )
            ? TAUT_NON : ( jj == 1 ) ? TAUT_INI : ( jj == 2 ) ? TAUT_YES : -1;

        if (j < 0)
        {
            continue;
        }

        if (composite_norm_data[j].bExists &&
             composite_norm_data[j].num_components > 1)
        {

            bDisplayTaut = ( !( ip->nMode & REQ_MODE_BASIC ) && !j ) ? -1 : j;
            nNumTautComponents = 0;
            if (bDisplayTaut)
            {
                /* find whether the structure is actually tautomeric */
                for (m = 0; m < composite_norm_data[TAUT_YES].num_components; m++)
                {
                    if (!pINChI2[m][TAUT_YES])
                    {
                        continue;
                    }
                    if (pINChI2[m][TAUT_YES]->bDeleted ||
                         pINChI2[m][TAUT_YES]->lenTautomer > 0)
                    {
                        nNumTautComponents++;
                    }
                }
            }

            for (k = 0;
                 k <= composite_norm_data[j].bHasIsotopicLayer && !sd->bUserQuitComponentDisplay;
                 k++)
            {
                /*  added number of components, added another format */
                /* for a single component case - DCh         */
                int bMobileH = ( bDisplayTaut > 0 && nNumTautComponents );
                sprintf(szTitle, "%s Structure #%ld%s%s.%s%s%s%s%s",
                    j == TAUT_INI ? "Preprocessed" : "Result for", num_inp,
                    bMobileH ? ", mobile H" :
                    bDisplayTaut == 0 ? ", fixed H" : "",
                    /*j? ", mobile H":", fixed H",*/
                    k ? ", isotopic" : "",
                    SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue),
                    iINChI ? " (Reconnected)" : "");
#ifndef TARGET_LIB_FOR_WINCHI
                /****** Display composite Result structure **************/
                nNumIntermediateTaut += ( j == TAUT_INI );
                /* display TAUT_INI (preprocessed) only once */
                if (j != TAUT_INI || nNumIntermediateTaut == 1)
                {
                    err =
                        DisplayCompositeStructure( pCG,
                                                   composite_norm_data,
                                                   orig_inp_data->polymer,
                                                   j == TAUT_INI ? 1 : k,
                                                        /* bIsotopic*/
                                                   j, /*tautomeric*/
                                                   j == TAUT_INI ? NULL
                                                               : pINChI2,
                                                   j == TAUT_INI ? NULL
                                                               : pINChI_Aux2,
                                                   ip->bAbcNumbers,
                                                   &ip->dp,
                                                   ip->nMode,
                                                   szTitle );
                }
                if ((sd->bUserQuitComponentDisplay = ( err == ESC_KEY ))) /* djb-rwth: addressing LLVM warning */
                {
                    break;
                }

                if (inp_data &&
                     inp_data->nEquLabels &&
                     inp_data->nNumEquSets &&
                     !sd->bUserQuitComponentDisplay &&
                     ( ( j == bCompareTaut || (bCompareTaut && j == TAUT_INI) ) ||
                         (bCompareTaut && !composite_norm_data[bCompareTaut].bExists) ) &&
                         ( k == bCompareIsotopic || (bCompareIsotopic && !composite_norm_data[j].bHasIsotopicLayer) ) /* djb-rwth: addressing LLVM warnings */
                      ) /* djb-rwth: addressing LLVM warning */
                {
                    AT_NUMB         nEquSet;
                    int             bDisplaySaved = ip->bDisplay;

                    /****** Display Equ Sets of composite Result structure **************/
                    for (nEquSet = 1;
                                    nEquSet <= inp_data->nNumEquSets;
                                                                    nEquSet++)
                    {
                        sprintf(szTitle, "Equ set %d of %d, %s Structure #%ld%s%s.%s%s%s%s%s",
                            nEquSet, inp_data->nNumEquSets,
                            j == TAUT_INI ? "Preprocessed" : "Result for",
                            num_inp,
                            (bDisplayTaut > 0 && nNumTautComponents) ? ", mobile H" : bDisplayTaut == 0 ? ", fixed H" : "",
                            /*j? ", mobile H":", fixed H",*/
                            k ? ", isotopic" : "",
                            SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue), iINChI ? " (Reconnected)" : "");
                        ip->dp.nEquLabels = inp_data->nEquLabels;
                        ip->dp.nCurEquLabel = nEquSet;
                        ip->dp.nNumEquSets = inp_data->nNumEquSets;
                        ip->bDisplay = 1; /* force display if it was not requested */
                        err = DisplayCompositeStructure( pCG,
                                                         composite_norm_data,
                                                         inp_data->polymer,
                                                         k, j,
                                                         pINChI2, pINChI_Aux2,
                                                         ip->bAbcNumbers,
                                                         &ip->dp,
                                                         ip->nMode,
                                                         szTitle );
                        ip->dp.nEquLabels = NULL;
                        ip->dp.nCurEquLabel = 0;
                        ip->dp.nNumEquSets = 0;
                        ip->bDisplay = bDisplaySaved; /* restore display option */

                        if ((sd->bUserQuitComponentDisplay = ( err == ESC_KEY ))) /* djb-rwth: addressing LLVM warning */
                        {
                            break;
                        }
                    }
                }
#else
                if (DRAWDATA && j <= TAUT_YES)
                {
                    struct DrawData vDrawData;
                    vDrawData.pWindowData =
                        CreateWinDataComposite_( pCG,
                                                 composite_norm_data,
                                                 k,
                                                 j,
                                                 pINChI2,
                                                 pINChI_Aux2,
                                                 ip->bAbcNumbers,
                                                 &ip->dp,
                                                 ip->nMode );

                    /* vDrawData.pWindowData = CreateWinData_( composite_norm_data[j].at, composite_norm_data[j].num_at,
                                         k, j, pINChI[i], pINChI_Aux[i],ip->bAbcNumbers, &ip->dp, ip->nMode ); */

                    if (vDrawData.pWindowData != NULL)
                    {
                        int nType;
                        vDrawData.nComponent = 0;
                        if (j == 0)
                            nType = ( k == 0 ) ? COMPONENT_BN : COMPONENT_BI;
                        else
                            nType = ( k == 0 ) ? COMPONENT_TN : COMPONENT_TI;
                        vDrawData.nType = nType;
                        vDrawData.bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                        vDrawData.szTitle = inchi__strdup( szTitle );
                        vDrawData.pWindowData->szTitle = inchi__strdup( szTitle );
                        if (inp_data && inp_data->nEquLabels &&
                             inp_data->nNumEquSets &&
                             ( j == bCompareTaut || bCompareTaut && !composite_norm_data[bCompareTaut].bExists ) &&
                             ( k == bCompareIsotopic || bCompareIsotopic && !composite_norm_data[j].bHasIsotopicLayer ) &&
                             ( vDrawData.pWindowData->nEquLabels = (AT_NUMB *) inchi_calloc( inp_data->num_inp_atoms,
                                 sizeof( inp_data->nEquLabels[0] ) ) ))
                        {
                            memcpy( vDrawData.pWindowData->nEquLabels, inp_data->nEquLabels,
                                     inp_data->num_inp_atoms * sizeof( inp_data->nEquLabels[0] ) );
                            vDrawData.pWindowData->nNumEquSets = inp_data->nNumEquSets;
                            vDrawData.pWindowData->nCurEquLabel = 0;
                        }
                        DRAWDATA( &vDrawData );
                    }
                }
                else if (DRAWDATA && GET_DRAWDATA && j == TAUT_INI)
                {
                    struct DrawData vDrawData;
                    struct DrawData *pDrawData;

                    if (!( ip->bCompareComponents & CMP_COMPONENTS ) ||
                        ( ip->bCompareComponents & CMP_COMPONENTS_NONTAUT ) ||
                           !k != !composite_norm_data[j].bHasIsotopicLayer)
                    {

                        continue;
                    }
                    vDrawData.pWindowData =
                        CreateWinDataComposite_( pCG,
                                                 composite_norm_data,
                                                 1 /*k*/,
                                                 j,
                                                 NULL,
                                                 NULL,
                                                 ip->bAbcNumbers,
                                                 &ip->dp,
                                                 ip->nMode );



                    if (vDrawData.pWindowData != NULL)
                    {
                        int nType = COMPONENT_ORIGINAL_PREPROCESSED;
                        pDrawData = GET_DRAWDATA( 0, nType, iINChI );
                        if (pDrawData)
                        {
                            FreeDrawData( pDrawData );
                            pDrawData->pWindowData = vDrawData.pWindowData;
                            vDrawData.pWindowData = NULL;
                        }
                        else
                        {
                            pDrawData = &vDrawData;
                        }

                        /* vDrawData.pWindowData = CreateWinData_( composite_norm_data[j].at, composite_norm_data[j].num_at,
                                        k, j, pINChI[i], pINChI_Aux[i],ip->bAbcNumbers, &ip->dp, ip->nMode ); */

                        pDrawData->nComponent = 0;
                        pDrawData->nType = nType;
                        pDrawData->bReconnected = iINChI; /* 0=>main; 1=>reconnected */
                        pDrawData->szTitle = inchi__strdup( szTitle );
                        pDrawData->pWindowData->szTitle = inchi__strdup( szTitle );
                        if (inp_data && inp_data->nEquLabels &&
                             inp_data->nNumEquSets &&
                             /*(j == bCompareTaut     || bCompareTaut     && !composite_norm_data[bCompareTaut].bExists) &&*/
                             /*(k == bCompareIsotopic || bCompareIsotopic && !composite_norm_data[j].bHasIsotopicLayer)  &&*/
                             ( pDrawData->pWindowData->nEquLabels = (AT_NUMB *) inchi_calloc( inp_data->num_inp_atoms,
                                 sizeof( inp_data->nEquLabels[0] ) ) ))
                        {
                            memcpy( pDrawData->pWindowData->nEquLabels,
                                    inp_data->nEquLabels,
                                    inp_data->num_inp_atoms * sizeof( inp_data->nEquLabels[0] ) );
                            pDrawData->pWindowData->nNumEquSets = inp_data->nNumEquSets;
                            pDrawData->pWindowData->nCurEquLabel = 0;
                        }
                        if (pDrawData == &vDrawData)
                        {
                            DRAWDATA( pDrawData );  /* there was no prepocessed structure */
                        }
                    }
                }

#endif /* #ifndef TARGET_LIB_FOR_WINCHI */
            }
        }
    }

    return err;
}


/****************************************************************************/
int DisplayTheWholeStructure( struct tagCANON_GLOBALS *pCG,
                              struct tagINCHI_CLOCK   *ic,
                              STRUCT_DATA             *sd,
                              INPUT_PARMS             *ip,
                              char                    *szTitle,
                              INCHI_IOSTREAM          *inp_file,
                              INCHI_IOSTREAM          *log_file,
                              ORIG_ATOM_DATA          *orig_inp_data,
                              long                    num_inp,
                              int                     iINChI,
                              int                     bShowStruct,
                              int                     bINCHI_LIB_Flag )
{

    int bDisplayEqu = 0;

#ifndef TARGET_LIB_FOR_WINCHI

    /* Displaying equivalent input structures when disconnection has been done: */
    /* in case of TARGET_LIB_FOR_WINCHI equivalence info is always unknown here and bOriginalReconnected=0 */

    int bOriginalReconnected = iINChI < 0 &&
        orig_inp_data && orig_inp_data->nEquLabels &&
        ( sd->bTautFlagsDone[INCHI_BAS] & TG_FLAG_DISCONNECT_COORD_DONE ) &&
        ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD );

    const char *lpszType =
        bOriginalReconnected ? " (Reconnected)"
        : ( iINChI < 0 ) ? ""
        : ( iINChI == INCHI_BAS ) ? " (Preprocessed)"
        : ( iINChI == INCHI_REC ) ? " (Reconnected)"
        : ""; /* djb-rwth: ignoring LLVM warning: variable used */

    int err = 0; /* djb-rwth: ignoring LLVM warning: variable used */

    /* Display the original structure */

    bDisplayEqu = bShowStruct && ip->bDisplay &&
        ip->dp.nEquLabels && 0 < ip->dp.nCurEquLabel && ip->dp.nCurEquLabel <= ip->dp.nNumEquSets;
#else
    if (!DRAWDATA || !DRAWDATA_EXISTS)
        return 0;
#endif


#if !defined(TARGET_API_LIB) && !defined(COMPILE_ANSI_ONLY)
    /* Ask the user whether to process the input structure or quit*/
    if (ip->bDisplay && inp_file->f != stdin)
    {
        if (user_quit( ic, bDisplayEqu ? "Enter=Display identical components, Esc=Stop ?"
            : "Enter=Display, Esc=Stop ?",
            ip->ulDisplTime ))
        {
            sd->bUserQuit = 1;
            goto exit_function;
        }
    }
#endif

    /*
     *    Display the whole input structure in console app
     */
    
/* #ifndef TARGET_LIB_FOR_WINCHI */
    
#if ( !defined( TARGET_LIB_FOR_WINCHI ) && !defined(TARGET_EXE_STANDALONE) )
    if (bShowStruct && ip->bDisplay)
    {
        if (bDisplayEqu)
        {
            sprintf( szTitle, " Equ Set %d of %d, Input Structure #%ld.%s%s%s%s%s",
                     ip->dp.nCurEquLabel, ip->dp.nNumEquSets,
                     num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ), lpszType );
        }
        else
        {
            sprintf( szTitle, "Input Structure #%ld.%s%s%s%s%s", num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ), lpszType );
        }
        err = DisplayStructure( orig_inp_data->at, orig_inp_data->num_inp_atoms, 
                                orig_inp_data->polymer,
                                0, 1, 0, NULL, 1/*isotopic*/, 0/*taut*/, NULL, NULL,
                                ip->bAbcNumbers, &ip->dp, ip->nMode, szTitle );
        sd->bUserQuitComponent = ( err == ESC_KEY );
        if (!err)
        {
            inchi_fprintf( stderr, "Cannot display the structure\n" );
        }
    }
    if (!bDisplayEqu)
    {
        /*  console output progress report */
        if (ip->bDisplay && !sd->bUserQuitComponent)
        {
            if (iINChI == 1)
            {
                if (ip->bDisplay)
                    inchi_ios_eprint( log_file, "Processing (rec) structure #%ld.%s%s%s%s...\n", num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
                else
                    inchi_fprintf( stderr, "Processing (rec) structure #%ld.%s%s%s%s...\r", num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
            }
            else
            {
                if (ip->bDisplay)
                    inchi_ios_eprint( log_file, "Processing structure #%ld.%s%s%s%s...\n", num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
                else
                    inchi_fprintf( stderr, "Processing structure #%ld.%s%s%s%s...\r", num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
            }
        }
    }
#endif


    /*
     *    Store the whole input structure in GUI application
     */

#ifdef TARGET_LIB_FOR_WINCHI
    if (ip->bDisplay && bINCHI_LIB_Flag)
#else
    if (( ip->bDisplay || ( ip->bCompareComponents & CMP_COMPONENTS ) ) && bINCHI_LIB_Flag)
#endif

    {
        int bBit, k, bReconnected, nComponent, bPreprocessed;

        for (bBit = 1, k = 0; k < 8; k++, bBit <<= 1)
        {
            /******************************************************************************
             *  bReconnected  = k%2     (0 or 1)
             *  nComponent    = k/4     (0 or 1)
             *  bPreprocessed = (k/2)%2 (0 or 1)
             ******************************************************************************/
            if (!( bINCHI_LIB_Flag & bBit ))
            {
                continue;
            }

            bReconnected = k % 2;
            nComponent = k / 4;
            bPreprocessed = ( ( k / 2 ) % 2 );

            sprintf(szTitle, "%s Structure #%ld.%s%s%s%s",
                bPreprocessed ? "Preprocessed" : bReconnected ? "Reconnected" : "Input",
                num_inp,
                SDF_LBL_VAL(ip->pSdfLabel, ip->pSdfValue));

#ifdef TARGET_LIB_FOR_WINCHI
            if (DRAWDATA && DRAWDATA_EXISTS)
            {
                struct DrawData vDrawData;
                int    nType = bPreprocessed ? COMPONENT_ORIGINAL_PREPROCESSED : COMPONENT_ORIGINAL;
                if (DRAWDATA_EXISTS( nComponent, bPreprocessed, bReconnected ))
                {
                    /* 
                    TEMPORARILY DISABLED TO MAKE WINCHI TREAT SOME POLYMERS
                    sd->nErrorType = _IS_FATAL;
                    sd->nErrorCode = CT_UNKNOWN_ERR;
                    return -1;
                    */
                }
                vDrawData.pWindowData = CreateWinData_( pCG,
                                                        orig_inp_data->at,
                                                        orig_inp_data->num_inp_atoms,
                                                        0,
                                                        1 /* bAdd_DT_to_num_H */,
                                                        0,
                                                        NULL,
                                                        1,
                                                        0,
                                                        NULL,
                                                        NULL,
                                                        ip->bAbcNumbers,
                                                        &ip->dp,
                                                        ip->nMode );
                if (vDrawData.pWindowData != NULL)
                {
                    vDrawData.nComponent = nComponent;
                    vDrawData.nType = nType; /* COMPONENT_ORIGINAL or COMPONENT_ORIGINAL_PREPROCESSED */
                    vDrawData.bReconnected = bReconnected; /* 0=>main; 1=>reconnected */
                    vDrawData.pWindowData->szTitle = inchi__strdup( szTitle );
                    vDrawData.szTitle = inchi__strdup( szTitle );
                    DRAWDATA( &vDrawData );
                    if (!nComponent)
                    {
                        /* keep track of saved INCHI_LIB data */
                        orig_inp_data->bSavedInINCHI_LIB[bReconnected] ++;
                        orig_inp_data->bPreprocessed[bReconnected] = bPreprocessed;
                    }
                }
            }
#else
            if (!nComponent && orig_inp_data) /* djb-rwth: fixing a NULL pointer dereference */
            {
                /* keep track of saved INCHI_LIB data */
                orig_inp_data->bSavedInINCHI_LIB[bReconnected] ++;
                orig_inp_data->bPreprocessed[bReconnected] = bPreprocessed;
            }
#endif
        }
    }

exit_function:

    return sd->bUserQuit;
}

#else
/* dummies */

/****************************************************************************/
int DisplayTheWholeStructure( struct tagCANON_GLOBALS *pCG,
                              struct tagINCHI_CLOCK   *ic,
                              STRUCT_DATA             *sd,
                              INPUT_PARMS             *ip,
                              char                    *szTitle,
                              INCHI_IOSTREAM          *inp_file,
                              INCHI_IOSTREAM          *log_file,
                              ORIG_ATOM_DATA          *orig_inp_data,
                              long                    num_inp,
                              int                     iINChI,
                              int                     bShowStruct,
                              int                     bINCHI_LIB_Flag )
{
    return 0;
}


/****************************************************************************/
int DisplayStructure( struct tagCANON_GLOBALS *pCG,
                      inp_ATOM                *at,
                      int                     num_at,
                      OAD_Polymer			  *polymer,
                      int                     num_removed_H,
                      int                     bAdd_DT_to_num_H,
                      int                     nNumRemovedProtons,
                      NUM_H                   *nNumRemovedProtonsIsotopic,
                      int                     bIsotopic,
                      int                     j /*bTautomeric*/,
                      INChI                   **cur_INChI,
                      INChI_Aux               **cur_INChI_Aux,
                      int                     bAbcNumbers,
                      DRAW_PARMS              *dp,
                      INCHI_MODE              nMode,
                      char                    *szTitle )
{
    return 0;
}




#endif /*  }COMPILE_ANSI_ONLY */


#ifndef TARGET_API_LIB


/****************************************************************************/
void SplitTime( unsigned long ulTotalTime,
                int *hours, int *minutes,
                int *seconds, int *mseconds )
{

    *mseconds = (int) ( ulTotalTime % 1000 );
    ulTotalTime /= 1000;
    *seconds = (int) ( ulTotalTime % 60 );
    ulTotalTime /= 60;
    *minutes = (int) ( ulTotalTime % 60 );
    ulTotalTime /= 60;
    *hours = (int) ( ulTotalTime );
}
#endif


/****************************************************************************
 Check if structure is chiral
****************************************************************************/
int bIsStructChiral( PINChI2 *pINChI2[INCHI_NUM], int num_components[] )
{
    int i, j, k;
    INChI *pINChI;
    INChI_Stereo *Stereo;

    for (j = 0; j < INCHI_NUM; j++)
    {
        /* disconnected / reconnected */
        if (!num_components[j])
        {
            continue;
        }

        for (i = 0; i < num_components[j]; i++)
        {
            /* i-th component */
            for (k = 0; k < TAUT_NUM; k++)
            {
                /* mobile/immobile H */
                if (( pINChI = pINChI2[j][i][k] ) &&
                      !pINChI->bDeleted                &&
                      pINChI->nNumberOfAtoms > 0)
                {

                    if (( Stereo = pINChI->Stereo ) &&
                         Stereo->t_parity &&
                         Stereo->nNumberOfStereoCenters > 0 &&
                         Stereo->nCompInv2Abs)
                    {
                        return 1; /* inversion changed stereo */
                    }
                    if (( Stereo = pINChI->StereoIsotopic ) &&
                         Stereo->t_parity &&
                         Stereo->nNumberOfStereoCenters > 0 &&
                         Stereo->nCompInv2Abs)
                    {
                        return 1; /* inversion changed stereo */
                    }
                }
            }
        }
    }

    return 0;
}


/****************************************************************************
 Release InChI work memory
****************************************************************************/
void FreeAllINChIArrays( PINChI2 *pINChI[INCHI_NUM],
                         PINChI_Aux2 *pINChI_Aux[INCHI_NUM],
                         int num_components[INCHI_NUM] )
{
    int k;
    for (k = 0; k < INCHI_NUM; k++)
    {
        int nk = num_components[k];

        FreeINChIArrays( pINChI[k], pINChI_Aux[k], num_components[k] );

        num_components[k] = 0;

        if (nk &&                /* added check for nk: 2013-12-15 IPl */
             pINChI[k])
        {
            inchi_free( pINChI[k] );
            pINChI[k] = NULL;
        }

        if (nk &&                /* added check for nk: 2013-12-15 IPl */
             pINChI_Aux[k])
        {
            inchi_free( pINChI_Aux[k] );
            pINChI_Aux[k] = NULL;
        }
    }

    return;
}


/****************************************************************************
 Release InChI work memory
****************************************************************************/
void FreeINChIArrays( PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, int num_components )
{
    int i, k;
    /* Release allocated memory */
    if (pINChI)
    {
        for (i = 0; i < num_components; i++)
        {
            for (k = 0; k < TAUT_NUM; k++)
            {
                Free_INChI( &pINChI[i][k] );
            }
        }
    }
    if (pINChI_Aux)
    {
        for (i = 0; i < num_components; i++)
        {
            for (k = 0; k < TAUT_NUM; k++)
            {
                Free_INChI_Aux( &pINChI_Aux[i][k] );
            }
        }
    }
}


/****************************************************************************
  Treat errors returned by CreateOneComponentINChI( ... )
****************************************************************************/
int TreatErrorsInCreateOneComponentINChI( STRUCT_DATA *sd,
                                          INPUT_PARMS    *ip,
                                          ORIG_ATOM_DATA *orig_inp_data,
                                          int i, long num_inp,
                                          INCHI_IOSTREAM *inp_file,
                                          INCHI_IOSTREAM *log_file,
                                          INCHI_IOSTREAM *out_file,
                                          INCHI_IOSTREAM *prb_file )
{
    if (sd->nErrorCode)
    {
        AddErrorMessage( sd->pStrErrStruct, ErrMsg( sd->nErrorCode ) );
        inchi_ios_eprint( log_file,
                          "Error %d (%s) structure #%ld component %d.%s%s%s%s\n",
                          sd->nErrorCode, sd->pStrErrStruct,
                          num_inp, i + 1, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        sd->nErrorType = ( sd->nErrorCode == CT_OUT_OF_RAM || sd->nErrorCode == CT_USER_QUIT_ERR )
            ? _IS_FATAL
            : _IS_ERROR;

#ifdef TARGET_LIB_FOR_WINCHI
        if (( ip->bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) &&
            ( ip->bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT ))
        {
            sd->nErrorType = ProcessStructError( out_file, log_file,
                                                 sd->pStrErrStruct, sd->nErrorType,
                                                 num_inp, ip );
            /*  Save the problem structure */
            if (prb_file->f &&
                 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd &&
                 !ip->bSaveAllGoodStructsAsProblem)
            {
                MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, num_inp );
            }
        }
        else
        {
            /*  Save the problem structure */
            if (sd->nErrorCode &&
                 prb_file->f &&
                 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd &&
                 !ip->bSaveAllGoodStructsAsProblem)
            {
                MolfileSaveCopy( inp_file, sd->fPtrStart, sd->fPtrEnd, prb_file->f, num_inp );
            }
        }
#endif
    }

/* #ifndef TARGET_API_LIB */
#if ( !defined( TARGET_API_LIB ) && !defined(TARGET_EXE_STANDALONE) )
    /*  print the logfile record */
    if (log_file->f && log_file->f != stderr && ( sd->ulStructTime >= 1000 || sd->nErrorCode ))
    {
        fprintf( log_file->f, "%10lu msec structure #%ld.%s%s%s%s (%d component%s, %d atom%s, error=%d).\n",
                sd->ulStructTime, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ),
                orig_inp_data->num_components, orig_inp_data->num_components == 1 ? "" : "s",
                orig_inp_data->num_inp_atoms, orig_inp_data->num_inp_atoms == 1 ? "" : "s", sd->nErrorCode );
    }
#endif

    return sd->nErrorType;
}


/****************************************************************************/
int TreatCreateINChIWarning( STRUCT_DATA    *sd,
                             INPUT_PARMS    *ip,
                             ORIG_ATOM_DATA *orig_inp_data,
                             long           num_inp,
                             INCHI_IOSTREAM *inp_file,
                             INCHI_IOSTREAM *log_file,
                             INCHI_IOSTREAM *out_file,
                             INCHI_IOSTREAM *prb_file )
{

#if ( bRELEASE_VERSION == 0 && (EXTR_FLAGS || EXTR_MASK) )
    if (EXTR_MASK ? ( ( sd->bExtract & EXTR_MASK ) == EXTR_FLAGS )
                    : ( sd->bExtract & EXTR_FLAGS )
       )
    {
        char szMsg[64];
        sprintf( szMsg, "ExtractStruct.code=0x%X", sd->bExtract );
        if (!ip->bNoWarnings)
        {
            WarningMessage( sd->pStrErrStruct, szMsg );
        }
    }
#endif

    if (!sd->nErrorCode && sd->pStrErrStruct[0])
    {
        inchi_ios_eprint( log_file, "Warning (%s) structure #%ld.%s%s%s%s\n",
            sd->pStrErrStruct, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
        sd->nErrorType = _IS_WARNING;

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
        /*  save the structure as a problem structure if requested */
        if (ip->bSaveWarningStructsAsProblem && !ip->bSaveAllGoodStructsAsProblem &&
             prb_file->f && 0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd)
        {   /* djb-rwth: addressing coverity ID #499545 -- return values handled properly */
            MolfileSaveCopy( inp_file,
                             sd->fPtrStart,
                             sd->fPtrEnd,
                             prb_file->f,
                             num_inp );
        }

#if ( bRELEASE_VERSION == 0 )
        /*  otherwise extract the structure as a problem structure if requested */
        else
            if (( EXTR_MASK ? ( ( sd->bExtract & EXTR_MASK ) == EXTR_FLAGS ) : ( sd->bExtract & EXTR_FLAGS ) )
                            && !ip->bSaveAllGoodStructsAsProblem &&
                               prb_file->f &&
                               0L <= sd->fPtrStart && sd->fPtrStart < sd->fPtrEnd)
            {
                MolfileSaveCopy( inp_file->f,
                                 sd->fPtrStart,
                                 sd->fPtrEnd,
                                 prb_file->f,
                                 num_inp );
            }
#endif
    }

#if ( bRELEASE_VERSION != 1 && bOUTPUT_ONE_STRUCT_TIME == 1 )
#ifndef TARGET_API_LIB
    if (log_file && log_file != stderr)
    {
        fprintf( log_file, "%10lu msec structure %1dD #%ld.%s%s%s%s (%d component%s, %d atom%s, error=%d).\n",
                sd->ulStructTime, orig_inp_data->num_dimensions, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ),
                orig_inp_data->num_components, orig_inp_data->num_components == 1 ? "" : "s",
                orig_inp_data->num_inp_atoms, orig_inp_data->num_inp_atoms == 1 ? "" : "s", sd->nErrorCode );
    }
#else
    if (log_file)
    {
        inchi_ios_eprint( log_file, "%10lu msec structure %1dD #%ld.%s%s%s%s (%d component%s, %d atom%s, error=%d).\n",
                sd->ulStructTime, orig_inp_data->num_dimensions, num_inp, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ),
                orig_inp_data->num_components, orig_inp_data->num_components == 1 ? "" : "s",
                orig_inp_data->num_inp_atoms, orig_inp_data->num_inp_atoms == 1 ? "" : "s", sd->nErrorCode );
    }
#endif
#endif

    return sd->nErrorType;
}


/****************************************************************************
 Treat errors/warnings
****************************************************************************/
int GetProcessingWarningsOneComponentInChI( INChI *cur_INChI[],
                                            INP_ATOM_DATA **inp_norm_data,
                                            STRUCT_DATA *sd,
                                            int bNoWarnings )
{
    int i, ret = 0;
    for (i = 0; i < TAUT_NUM; i++)
    {
        if (cur_INChI[i] && cur_INChI[i]->nNumberOfAtoms > 0)
        {
            ret |= GetProcessingWarningsOneInChI( cur_INChI[i],
                                                  inp_norm_data[i],
                                                  sd->pStrErrStruct,
                                                  bNoWarnings );
        }
    }

    return ret;
}


/****************************************************************************
 Treat errors/warnings
****************************************************************************/
int GetProcessingWarningsOneInChI( INChI *pINChI,
                                   INP_ATOM_DATA *inp_norm_data,
                                   char *pStrErrStruct,
                                   int bNoWarnings )
{
    int j;
    int nAmbiguousStereoAtoms, nAmbiguousStereoBonds;
    nAmbiguousStereoAtoms = 0;
    nAmbiguousStereoBonds = 0;

    if (inp_norm_data->at)
    {
        for (j = 0; j < pINChI->nNumberOfAtoms; j++)
        {
            if (inp_norm_data->at[j].bAmbiguousStereo & ( AMBIGUOUS_STEREO_ATOM | AMBIGUOUS_STEREO_ATOM_ISO ))
            {
                nAmbiguousStereoAtoms++;
            }
            if (inp_norm_data->at[j].bAmbiguousStereo & ( AMBIGUOUS_STEREO_BOND | AMBIGUOUS_STEREO_BOND_ISO ))
            {
                nAmbiguousStereoBonds++;
            }
        }
        if (nAmbiguousStereoAtoms)
        {
            if (!bNoWarnings)
            {
                WarningMessage( pStrErrStruct, "Ambiguous stereo:" );
                WarningMessage( pStrErrStruct, "center(s)" );
            }
        }
        if (nAmbiguousStereoBonds)
        {
            if (!bNoWarnings)
            {
                WarningMessage( pStrErrStruct, "Ambiguous stereo:" );
                WarningMessage( pStrErrStruct, "bond(s)" );
            }
        }
    }

    return ( nAmbiguousStereoAtoms || nAmbiguousStereoBonds );
}
