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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <locale.h>

#include "mode.h"

#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "ichimain.h"
#include "extr_ct.h"
#include "ichi_io.h"
#include "ichicomp.h"
#include "inchi_api.h"
#include "readinch.h"

#include "ichitaut.h"
#include "ichicant.h"
#include "ichitime.h"
#include "bcf_s.h"

#include "inchi_dll.h"

/*************************************************************************
 *
 *   Local prototypes
 *
 *************************************************************************/

int SetAtomProperties( inp_ATOM *at,
                       MOL_COORD *szCoord,
                       inchi_Atom *ati,
                       int a1,
                       int *nDim,
                       char *pStrErr,
                       int *err );
void SetNumImplicitH( inp_ATOM* at, int num_atoms );
int SetBondProperties( inp_ATOM *at,
                       inchi_Atom *ati,
                       int a1,
                       int j,
                       int nNumAtoms,
                       int *nNumBonds,
                       char *pStrErr,
                       int *err );
int SetAtomAndBondProperties( inp_ATOM *at,
                              inchi_Atom *ati,
                              int a1,
                              int bDoNotAddH,
                              char *pStrErr,
                              int *err );
int InpAtom0DToInchiAtom( inp_ATOM *at,
                          int num_inp_atoms,
                          AT_NUM *num_atoms,
                          inchi_Atom **atom,
                          AT_NUM *num_stereo0D,
                          inchi_Stereo0D **stereo0D );
int ExtractOneStructure( STRUCT_DATA *sd,
                         INPUT_PARMS *ip,
                         char *szTitle,
                         inchi_InputEx *inp,
                         INCHI_IOSTREAM *log_file,
                         INCHI_IOSTREAM *out_file,
                         INCHI_IOSTREAM *prb_file,
                         ORIG_ATOM_DATA *orig_inp_data,
                         long *num_inp );

static int GetINCHI1( inchi_InputEx *inp, inchi_Output *out, int enforce_std_format );

int SetExtOrigAtDataByInChIExtInput( OAD_Polymer **ppPolymer,
                                     OAD_V3000 **ppV3000,
                                     inchi_Input_Polymer *polymer,
                                     inchi_Input_V3000 *v3000,
                                     int nat );
int SetInChIExtInputByExtOrigAtData( OAD_Polymer *pPolymer,
                                     OAD_V3000 *pV3000,
                                     inchi_Input_Polymer **ipolymer,
                                     inchi_Input_V3000 **iv3000,
                                     int nat );

/****************************************************************************/

int bInterrupted = 0;



/****************************************************************************
 *
 * INCHI API
 *
 ****************************************************************************/



/****************************************************************************

    FreeINCHI

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
void INCHI_DECL FreeINCHI( inchi_Output *out )
{
    if (!out)
    {
        return;
    }

    if (out->szInChI)
    {
        inchi_free( out->szInChI );
    }
    if (out->szLog)
    {
        inchi_free( out->szLog );
    }
    if (out->szMessage)
    {
        inchi_free( out->szMessage );
    }

    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */
}


/****************************************************************************

    FreeStdINCHI

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
void INCHI_DECL FreeStdINCHI( inchi_Output *out )
{
    FreeINCHI( out );
}



/****************************************************************************

    FreeStructFromStdINCHI
****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
void INCHI_DECL FreeStructFromStdINCHI( inchi_OutputStruct *out )
{
    FreeStructFromINCHI( out );
}



/****************************************************************************

    FreeStructFromINCHI

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
void INCHI_DECL FreeStructFromINCHI( inchi_OutputStruct *out )
{
    if (!out)
    {
        return;
    }

    if (out->atom)
    {
        inchi_free( out->atom );
    }
    if (out->stereo0D)
    {
        inchi_free( out->stereo0D );
    }
    if (out->szLog)
    {
        inchi_free( out->szLog );
    }
    if (out->szMessage)
    {
        inchi_free( out->szMessage );
    }

    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */
}


/****************************************************************************

    GetStdINCHI

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
int INCHI_DECL GetStdINCHI( inchi_Input *inp, inchi_Output *out )
{
    inchi_InputEx extended_input;

    /* No '*' or 'Zz' elements are allowed in the input . */
    if (input_erroneously_contains_pseudoatoms(inp, out))
    {
        return _IS_ERROR;
    }

    extended_input.atom = inp->atom;
    extended_input.num_atoms = inp->num_atoms;
    extended_input.num_stereo0D = inp->num_stereo0D;
    extended_input.stereo0D = inp->stereo0D;
    extended_input.szOptions = inp->szOptions;
    extended_input.polymer = NULL;
    extended_input.v3000 = NULL;

    return GetINCHI1( &extended_input, out, 1 );
}


/****************************************************************************

    GetINCHI

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
int INCHI_DECL GetINCHI( inchi_Input *inp, inchi_Output *out )
{
    inchi_InputEx extended_input;

    /* For back compatibility: no '*' or 'Zz' elements are allowed in the input to GetINCHI() ! */
    if ( input_erroneously_contains_pseudoatoms( inp, out) )
    {
        return _IS_ERROR;
    }

    extended_input.atom = inp->atom;
    extended_input.num_atoms = inp->num_atoms;
    extended_input.stereo0D = inp->stereo0D;
    extended_input.num_stereo0D = inp->num_stereo0D;
    extended_input.szOptions = inp->szOptions;
    extended_input.polymer = NULL;
    extended_input.v3000 = NULL;

    return GetINCHI1( &extended_input, out, 0 );
}


/****************************************************************************/
int input_erroneously_contains_pseudoatoms( inchi_Input *inp,
                                            inchi_Output *out)
{
    char *str_noz = "Unsupported in this mode element \'*\'";
    int i;
    /* Supposed that no '*' or 'Zz' elements are allowed in the input. */
    for (i = 0; i < inp->num_atoms; i++)
    {
        if (!strcmp(inp->atom->elname, "Zz") || !strcmp(inp->atom->elname, "*"))
        {
            if (out)
            {
                memset(out, 0, sizeof(*out)); /* djb-rwth: memset_s C11/Annex K variant? */
                if ((out->szMessage = (char *)inchi_malloc(strlen(str_noz) + 1))) /* djb-rwth: addressing LLVM warning */
                {
                    strcpy(out->szMessage, str_noz);
                }
            }
            return 1;
        }
    }

    return 0;
}


/****************************************************************************

    GetINCHIEx

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
int INCHI_DECL GetINCHIEx( inchi_InputEx *inp, inchi_Output *out )
{
    int i;

    /* Check for star atoms and replace them by Zz atoms */
    for (i = 0; i < inp->num_atoms; i++)
    {
        if (!strcmp( inp->atom[i].elname, "*" ))
        {
            strcpy( inp->atom[i].elname, "Zz" );
        }
    }

    return GetINCHI1( inp, out, 0 );
}


/****************************************************************************
    GetINCHI1 (major worker)
****************************************************************************/
static int GetINCHI1( inchi_InputEx *extended_input,
                      inchi_Output *out,
                      int enforce_std_format )
{
    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;
    char szTitle[MAX_SDF_HEADER + MAX_SDF_VALUE + 256];

    int i;
    long num_inp, num_err; /* djb-rwth: ignoring LLVM warning: variable used */
    char      szSdfDataValue[MAX_SDF_VALUE + 1];
    PINChI2     *pINChI[INCHI_NUM];
    PINChI_Aux2 *pINChI_Aux[INCHI_NUM];

    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */
    unsigned long  ulTotalProcessingTime = 0; /* djb-rwth: ignoring LLVM warning: variable used */

    INPUT_PARMS inp_parms;
    INPUT_PARMS *ip = &inp_parms;

    ORIG_ATOM_DATA OrigAtData; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *orig_inp_data = &OrigAtData;
    ORIG_ATOM_DATA PrepAtData[2]; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *prep_inp_data = PrepAtData;
    int             bReleaseVersion = bRELEASE_VERSION;
    int   nRet = 0, nRet1;

    CANON_GLOBALS CG;
    INCHI_CLOCK ic;

    STRUCT_FPTRS *pStructPtrs = NULL;

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif

    const char *argv[INCHI_MAX_NUM_ARG + 1];
    int   argc;
    char *szOptions = NULL;

    INCHI_IOSTREAM inchi_file[3], *out_file = inchi_file, *log_file = inchi_file + 1;
    INCHI_IOSTREAM prb_file0, *prb_file = &prb_file0;
    INCHI_IOS_STRING temp_string_container;
    INCHI_IOS_STRING *strbuf = &temp_string_container;

    inchi_Input prev_versions_input;
    inchi_Input *pvinp = &prev_versions_input;

#ifdef GHI100_FIX
#if ((SPRINTF_FLAG != 1) && (SPRINTF_FLAG != 2))
    setlocale(LC_ALL, "en-US"); /* djb-rwth: setting all locales to "en-US" */
#endif
#endif

    pvinp->atom = extended_input->atom;
    pvinp->num_atoms = extended_input->num_atoms;
    pvinp->num_stereo0D = extended_input->num_stereo0D;
    pvinp->stereo0D = extended_input->stereo0D;
    pvinp->szOptions = extended_input->szOptions;

#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag( _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF );
/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE
    _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode( _CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG );
#endif

#if ( !defined(__STDC__) || __STDC__ != 1 )
    /* turn on floating point exceptions */
    {
        /* Get the default control word. */
        int cw = _controlfp( 0, 0 );

        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &= ~( EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL );

        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif
#endif

    szTitle[0] = '\0';

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    repeat:
          inchi_ios_close( out_file );
          inchi_ios_close( log_file );
          inchi_ios_close( prb_file );
          pStr = NULL;
#endif

    /* Initialize internal for this function output streams as string buffers */
    inchi_ios_init( out_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( log_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( prb_file, INCHI_IOS_TYPE_STRING, NULL );

    num_inp = 0;
    num_err = 0;
    sd->bUserQuit = 0;

    /* clear original input structure */
    memset( pINChI, 0, sizeof( pINChI ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( pINChI_Aux, 0, sizeof( pINChI_Aux ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( ip, 0, sizeof( *ip ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( orig_inp_data, 0, sizeof( *orig_inp_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( prep_inp_data, 0, 2 * sizeof( *prep_inp_data ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( szSdfDataValue, 0, sizeof( szSdfDataValue ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    memset( &CG, 0, sizeof( CG ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &ic, 0, sizeof( ic ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    if (!out)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }
    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    /* options */
    if (pvinp && pvinp->szOptions)
    {
        szOptions = (char*) inchi_malloc( strlen( pvinp->szOptions ) + 1 );
        if (szOptions)
        {
            strcpy( szOptions, pvinp->szOptions );
            argc = parse_options_string( szOptions, argv, INCHI_MAX_NUM_ARG );
        }
        else
        {
            nRet = _IS_FATAL;
            goto translate_RetVal; /* emergency exit */
        }
    }
    else
    {
        argc = 1;
        argv[0] = "";
        argv[1] = NULL;
    }

    if ((argc == 1
#ifdef TARGET_API_LIB
              && ( !pvinp || pvinp->num_atoms <= 0 || !pvinp->atom ))
#endif
              || (argc == 2 && ( argv[1][0] == INCHI_OPTION_PREFX ) &&
                    ( !strcmp( argv[1] + 1, "?" ) || !inchi_stricmp( argv[1] + 1, "help" )) )) /* djb-rwth: addressing LLVM warnings */
    {
        HelpCommandLineParms( log_file );
        out->szLog = log_file->s.pStr;
        memset( log_file, 0, sizeof( *log_file ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        nRet = _IS_EOF;
        goto translate_RetVal;
    }

    nRet1 = ReadCommandLineParms( argc, argv, ip, szSdfDataValue, &ulDisplTime, bReleaseVersion, log_file );
    if (szOptions)
    {
        inchi_free( szOptions );
        szOptions = NULL;
    }
    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if (0 > nRet1)
    {
        nRet = _IS_FATAL;
        goto exit_function;
    }
    if (ip->bNoStructLabels)
    {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    }
    else
    {
        if (ip->nInputType == INPUT_INCHI_XML || ip->nInputType == INPUT_INCHI_PLAIN || ip->nInputType == INPUT_CMLFILE)
        {
            /* the input may contain both the header and the label of the structure */
            if (!ip->pSdfLabel)
                ip->pSdfLabel = ip->szSdfDataHeader;
            if (!ip->pSdfValue)
                ip->pSdfValue = szSdfDataValue;
        }
    }

    /* Ensure standardness */
    if (enforce_std_format)
    {
        if (ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT)
        {
            ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;
        }
        if (0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD ))
        {
            ip->bTautFlags &= ~TG_FLAG_RECONNECT_COORD;
        }
        if (0 != ( ip->nMode & REQ_MODE_BASIC ))
        {
            ip->nMode &= ~REQ_MODE_BASIC;
        }
        if (0 != ( ip->nMode & REQ_MODE_RELATIVE_STEREO ))
        {
            ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO );
        }
        if (0 != ( ip->nMode & REQ_MODE_RACEMIC_STEREO ))
        {
            ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO );
        }
        if (0 != ( ip->nMode & REQ_MODE_CHIR_FLG_STEREO ))
        {
            ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO );
        }
        if (0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO ))
        {
            ip->nMode &= ~REQ_MODE_DIFF_UU_STEREO;
        }
        if (0 == ( ip->nMode & ( REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU ) ))
        {
            ip->nMode |= REQ_MODE_SB_IGN_ALL_UU;
            ip->nMode |= REQ_MODE_SC_IGN_ALL_UU;
        }
        if (0 != ( ip->bTautFlags & TG_FLAG_KETO_ENOL_TAUT ))
        {
            ip->bTautFlags &= ~TG_FLAG_KETO_ENOL_TAUT;
        }
        if (0 != ( ip->bTautFlags & TG_FLAG_1_5_TAUT ))
        {
            ip->bTautFlags &= ~TG_FLAG_1_5_TAUT;
        }
        /* And anyway... */
        ip->bINChIOutputOptions |= INCHI_OUT_STDINCHI;
        ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;
    }
    /* */

    PrintInputParms( log_file, ip );

    if (0 >= inchi_strbuf_init( strbuf, INCHI_STRBUF_INITIAL_SIZE, INCHI_STRBUF_SIZE_INCREMENT ))
    {
        inchi_ios_eprint( log_file, "Cannot allocate internal string buffer. Terminating\n" );
        nRet = _IS_FATAL;
        goto exit_function;
    }

    /***************************************************
    Main cycle -- read input structures and create their INChI's */ /* djb-rwth: addressing LLVM warning */
    ulTotalProcessingTime = 0;

    if (pStructPtrs)
    {
        memset( pStructPtrs, 0, sizeof( pStructPtrs[0] ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    }

    /* === possible improvement: convert inp to orig_inp_data ==== */
    if (!sd->bUserQuit && !bInterrupted)
    {
        if (ip->last_struct_number && num_inp >= ip->last_struct_number)
        {
            nRet = _IS_EOF; /*  simulate end of file */
            goto exit_function;
        }

        nRet = ExtractOneStructure( sd,ip, szTitle, extended_input,
                                    log_file, out_file, prb_file,
                                    orig_inp_data, &num_inp );

        if (pStructPtrs)
        {
            pStructPtrs->cur_fptr++;
        }

#ifndef TARGET_API_LIB
        if (sd->bUserQuit)
        {
            break;
        }
#endif
        switch (nRet)
        {
            case _IS_FATAL:
                num_err++;
                goto exit_function;
            case _IS_EOF:
                goto exit_function;
            case _IS_ERROR:
                num_err++;
                goto exit_function;
#ifndef TARGET_API_LIB
            case _IS_SKIP:
                continue;
#endif
        }

        /* Create INChI for each connected component of the structure and */
        /* optionally display them ; output INChI for the whole structure */

        nRet1 = ProcessOneStructureEx( &ic, &CG, sd, ip, szTitle,
                                        pINChI, pINChI_Aux,
                                        NULL, /* inp_file is not necessary as all input is already saved in 'ip' */
                                        log_file, out_file, prb_file,
                                        orig_inp_data, prep_inp_data,
                                        num_inp, strbuf, 0 /* save_opt_bits */ );

        /*  Free INChI memory */
        FreeAllINChIArrays( pINChI, pINChI_Aux, sd->num_components );

        /* Free structure data */
        FreeOrigAtData( orig_inp_data );
        FreeOrigAtData( prep_inp_data );
        FreeOrigAtData( prep_inp_data + 1 );

        ulTotalProcessingTime += sd->ulStructTime;
        nRet = inchi_max( nRet, nRet1 );
        switch (nRet)
        {
            case _IS_FATAL:
                /* num_err ++; */
                goto exit_function;
            case _IS_ERROR:
                ; /* num_err ++; */
#ifndef TARGET_API_LIB
                continue;
#endif
        }
    }

exit_function:
    /* Avoid memory leaks in case of fatal error */
    if (pStructPtrs && pStructPtrs->fptr)
    {
        inchi_free( pStructPtrs->fptr );
    }
    /* Free INChI memory */
    FreeAllINChIArrays( pINChI, pINChI_Aux, sd->num_components );
    /*    Free structure data */
    FreeOrigAtData( orig_inp_data );
    FreeOrigAtData( prep_inp_data );
    FreeOrigAtData( prep_inp_data + 1 );

    inchi_strbuf_close( strbuf );

    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (ip->path[i])
        {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( &CG );

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if (num_repeat-- > 0)
    {
        goto repeat;
    }
#endif

    /* output */
    produce_generation_output( out, sd, ip, log_file, out_file );

translate_RetVal:

    /* Close inernal I/O streams */
    inchi_ios_close( log_file );
    inchi_ios_close( out_file );
    inchi_ios_close( prb_file );

    switch (nRet)
    {
        case _IS_SKIP: nRet = inchi_Ret_SKIP; break; /* not used in INChI dll */
        case _IS_EOF: nRet = inchi_Ret_EOF; break; /* no structural data has been provided */
        case _IS_OKAY: nRet = inchi_Ret_OKAY; break; /* Success; break; no errors or warnings */
        case _IS_WARNING: nRet = inchi_Ret_WARNING; break; /* Success; break; warning(s) issued */
        case _IS_ERROR: nRet = inchi_Ret_ERROR; break; /* Error: no INChI has been created */
        case _IS_FATAL: nRet = inchi_Ret_FATAL; break; /* Severe error: no INChI has been created (typically; break; memory allocation failed) */
        case _IS_UNKNOWN:
        default: nRet = inchi_Ret_UNKNOWN; break; /* Unlnown program error */
    }

    return nRet;
}


/****************************************************************************/
void produce_generation_output( inchi_Output *out,
                                STRUCT_DATA *sd,
                                INPUT_PARMS *ip,
                                INCHI_IOSTREAM *log_file,
                                INCHI_IOSTREAM *out_file )

{
    if (sd->pStrErrStruct[0])
    {
        if (out && ( out->szMessage = (char *) inchi_malloc( strlen( sd->pStrErrStruct ) + 1 ) ))
        {
            strcpy(out->szMessage, sd->pStrErrStruct);
        }
    }

    /* Make separate strings with InChI and AuxInfo */
    if (out_file->s.pStr && out_file->s.nUsedLength > 0 && out)
    {
        char *p;
        out->szInChI = out_file->s.pStr;
        out->szAuxInfo = NULL;
        if (!( INCHI_OUT_SDFILE_ONLY & ip->bINChIOutputOptions )) /* do not remove last LF from SDF output - 2008-12-23 DT */
        {
            for (p = strchr( out->szInChI, '\n' ); p; p = strchr( p + 1, '\n' ))
            {
                if (!memcmp( p, "\nAuxInfo", 8 ))
                {
                    *p = '\0';            /* remove LF after INChI */
                    out->szAuxInfo = p + 1; /* save pointer to AuxInfo */
                }
                else if (out->szAuxInfo || !p[1])
                {
                    /* remove LF after aux info or from the last char */
                    *p = '\0';
                    break;
                }
            }
        }
        out_file->s.pStr = NULL;
    }

    copy_corrected_log_tail( out, log_file );
}


/****************************************************************************/
void copy_corrected_log_tail( inchi_Output *out, INCHI_IOSTREAM *log_file )
{
    if (log_file->s.pStr && log_file->s.nUsedLength > 0)
    {
        while (log_file->s.nUsedLength &&
                '\n' == log_file->s.pStr[log_file->s.nUsedLength - 1])
        {
            log_file->s.pStr[--log_file->s.nUsedLength] = '\0';
                                            /* remove last LF */
        }
        if (out)
        {
            char *p;
            out->szLog = log_file->s.pStr;
            log_file->s.pStr = NULL;
            for (p = strchr( out->szLog, ' ' ); p; p = strchr( p + 1, ' ' ))
            {
                if (!memcmp( p, " structure #", 12 ))
                {
                    *p = '\0';
                }
            }
        }
    }
}


/****************************************************************************

    CheckINCHI

    Check if the string represents valid InChI/standard InChI.
    Input:
            szINCHI     source InChI
            strict      if 0, just quickly check for proper layout
                        (prefix, version, etc.)
                        The result may not be strict.
                        If not 0, try to perform InChI2InChI conversion and
                        returns success if a resulting InChI string exactly
                        match source.
                        The result may be 'false alarm' due to imperfect algorithm of
                        conversion.
    Returns:
            success/errors codes

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
int INCHI_DECL CheckINCHI( const char *szINCHI, const int strict )
{
    int ret = INCHI_VALID_NON_STANDARD;
    int ret_i2i;
    inchi_InputINCHI    inchi_inp;
    inchi_Output        inchi_out;
    size_t slen, pos_slash1 = 0;
    char *str = NULL;
    size_t i;
    size_t slen0;
    char pp;

    /* .. non-empty */
    if (szINCHI == NULL)
    {
        return INCHI_INVALID_PREFIX;
    }

    slen = strlen( szINCHI );


    /* .. has valid prefix */
    if (slen < LEN_INCHI_STRING_PREFIX + 3)
    {
        return INCHI_INVALID_PREFIX;
    }
    if (memcmp( szINCHI, INCHI_STRING_PREFIX, LEN_INCHI_STRING_PREFIX ))
    {
        return INCHI_INVALID_PREFIX;
    }

    /* .. has InChI version 1 */
    /* if (!isdigit(szINCHI[LEN_INCHI_STRING_PREFIX]) )  */
    if (szINCHI[LEN_INCHI_STRING_PREFIX] != '1')
    {
        return INCHI_INVALID_VERSION;
    }

    /* .. optionally has a 'standard' flag character */
    pos_slash1 = LEN_INCHI_STRING_PREFIX + 1;
    if (szINCHI[pos_slash1] == 'S')
    {
        /* Standard InChI ==> standard InChIKey */
        ret = INCHI_VALID_STANDARD;
        pos_slash1++;
    }
    else if (szINCHI[pos_slash1] == 'B')
    {
        /* Beta version InChI ==> non-standard */
        ret = INCHI_VALID_BETA;
        pos_slash1++;
    }

    /* .. has trailing slash in the right place */
    if (szINCHI[pos_slash1] != '/')
    {
        return INCHI_INVALID_LAYOUT;
    }

    /* .. the rest of source string contains valid literals */


    /* adjust line len so we not check trailing whitespaces */
    i = slen - 1;
    while (isspace(UCINT szINCHI[i--])) slen--;

    /* Treat possible SaveOpt letters  */
    slen0 = slen;
    if (( szINCHI[slen - 3] == '\\' ) &&
        ( szINCHI[slen - 2] >= 'A' ) && ( szINCHI[slen - 2] <= 'Z' ) &&
        ( szINCHI[slen - 1] >= 'A' ) && ( szINCHI[slen - 1] <= 'Z' )
        )
    {
        slen0 = slen - 3;
    }

    int prev_is_slash = 1;
    for (i = pos_slash1 + 1; i < slen0; i++)
    {
        pp = szINCHI[i];
#if ( FIX_GAF_2020_GENERIC==1 )
        if (prev_is_slash)
        {
            /* After slash: */
            if (pp == '0')
            {
                /* '0' is never allowed */
                return INCHI_INVALID_LAYOUT;
            }
            if (i > pos_slash1 + 1)
            {
                /* Not in main formula layer... */ 
                if (!islower(pp))
                {
                    /* only lowercase letters are allowed */
                    return INCHI_INVALID_LAYOUT;
                }
            }
        }
        prev_is_slash = (pp != '/') ? 0 : 1;
#endif
        if (pp >= 'A' && pp <= 'Z')   continue;
        if (pp >= 'a' && pp <= 'z')   continue;
        if (pp >= '0' && pp <= '9')  continue;
        switch (pp)
        {
            case '(': case ')':
            case '*': case '+':
            case ',': case '-':
            case '.': case '/':
#if ( FIX_GAF_2020_GENERIC==1 )
            case ';': case '?':     continue;
#else
            case ';': case '=':
            case '?': case '@':     continue;
#endif
            default:            return INCHI_INVALID_LAYOUT;
        }
    }

    if (strict)
    {
        char opts[] = "?FixedH ?RecMet ?SUU ?SLUUD";
        extract_inchi_substring( &str, szINCHI, slen );
        if (NULL == str)
        {
            ret = INCHI_FAIL_I2I;
            goto fin;
        }

        inchi_inp.szInChI = str;
        opts[0] = opts[8] = opts[16] = opts[21] = INCHI_OPTION_PREFX;
        inchi_inp.szOptions = opts;

        ret_i2i = GetINCHIfromINCHI( &inchi_inp, &inchi_out );

        if (( ( ret_i2i != inchi_Ret_OKAY ) && ( ret_i2i != inchi_Ret_WARNING ) ) || !inchi_out.szInChI)
        {
            ret = INCHI_FAIL_I2I;
        }
        else
        {
            if (strcmp( inchi_inp.szInChI, inchi_out.szInChI ))
            {
                ret = INCHI_FAIL_I2I;
            }
        }
    }

fin:if (strict)
{
    if (NULL != str)
        inchi_free( str );
}

    return ret;
}


/****************************************************************************/
void SetNumImplicitH( inp_ATOM* at, int num_atoms )
{
    int bNonMetal;
    int a1/*, n1*/;

    /* special valences */
    for (bNonMetal = 0; bNonMetal < 2; bNonMetal++)
    {
        for (a1 = 0; a1 < num_atoms; a1++)
        {
            int bHasMetalNeighbor /*, j*/;
            if (bNonMetal != is_el_a_metal( at[a1].el_number ))
            {
                continue; /* first process all metals, after that all non-metals */
            }

            bHasMetalNeighbor = 0;
            /***********************************************************************
             *  Set number of hydrogen atoms
             */
            at[a1].num_H = get_num_H( at[a1].elname,
                                      at[a1].num_H,
                                      at[a1].num_iso_H,
                                      at[a1].charge,
                                      at[a1].radical,
                                      at[a1].chem_bonds_valence,
                                      0, /* instead of valence entered by the user: it does not exist here*/
                                      ( at[a1].at_type & 1 )  /* bAliased */,
                                      !( at[a1].at_type & 2 ) /* bDoNotAddH */,
                                      bHasMetalNeighbor );
            at[a1].at_type = 0;
        }
    }
}


/****************************************************************************/


#define REPEAT_ALL  0


/****************************************************************************/
int parse_options_string( char *cmd, const char *argv[], int maxargs )
{
    char    *p;
    char    *pArgCurChar;
    int      bInsideQuotes;
    int      bCopyCharToArg;
    int      nNumBackSlashes;
    int      i;

    i = 0;
    argv[i++] = ""; /* zeroth argument is not used */
    p = cmd;
    bInsideQuotes = 0;

    /* arguments, one by one */
    while (i < maxargs - 1)
    {
        /* bypass spaces */
        while (*p == ' ' || *p == '\t')
        {
            p++;
        }
        if (!*p)
        {
            break;
        }

        /* scan an argument */
        argv[i++] = pArgCurChar = p;     /* store preliminary ptr to arg */

        while (1)
        {
            bCopyCharToArg = 1;
            nNumBackSlashes = 0;
            while (*p == '\\')
            {
                ++p;
                ++nNumBackSlashes;
            }

            /* each pair of backslashes => one backslash; one more backslash => literal quote */
            if (*p == '\"')
            {
                /* one " found */
                if (nNumBackSlashes % 2 == 0)
                {
                    if (bInsideQuotes)
                    {
                        if (*( p + 1 ) == '\"')
                        {
                            p++;
                        }
                        else
                        {
                            bCopyCharToArg = 0;
                        }
                    }
                    else
                    {
                        bCopyCharToArg = 0;
                    }
                    bInsideQuotes = !bInsideQuotes;
                }
                nNumBackSlashes /= 2;          /* divide nNumBackSlashes by two */
            }
            while (nNumBackSlashes--)
            {
                *pArgCurChar++ = '\\';
            }
            if (!*p)
            {
                break;
            }
            if (!bInsideQuotes && ( *p == ' ' || *p == '\t' ))
            {
                p++;
                /* move to the next char because this char may become
                 * zero due to  *pArgCurChar++ = '\0'; line below */
                break;
            }
            if (bCopyCharToArg)
            {
                *pArgCurChar++ = *p;
            }
            ++p;
        }
        *pArgCurChar++ = '\0';  /* argument zero termination */
    }

    /* The last argument is NULL */
    argv[i] = NULL;

    return i;
}


/****************************************************************************/

#define MIN_BOND_LENGTH   (1.0e-6)


/****************************************************************************/
int SetAtomProperties( inp_ATOM *at,
                       MOL_COORD *szCoord,
                       inchi_Atom *ati,
                       int a1,
                       int *nDim,
                       char *pStrErr,
                       int *err )
{
    S_CHAR      cRadical;

    /* element, check later */
    strcpy( at[a1].elname, ati[a1].elname );

    /* charge */
    at[a1].charge = ati[a1].charge;

    /* radical */
    switch (ati[a1].radical)
    {
        case   INCHI_RADICAL_NONE:
            cRadical = 0;
            break;
        case   INCHI_RADICAL_SINGLET:
#if( SINGLET_IS_TRIPLET == 1) /* 'singlet' means two electrons make a lone pair instead of 2 bonds*/
                              /* its effect on valence is same as the effect of a triplet */
            cRadical = RADICAL_TRIPLET;
#else
            cRadical = RADICAL_SINGLET;
#endif
            break;
        case   INCHI_RADICAL_DOUBLET:
            cRadical = RADICAL_DOUBLET;
            break;
        case   INCHI_RADICAL_TRIPLET:
            cRadical = RADICAL_TRIPLET;
            break;
        default:
        {
            char szRadicalType[16];
            int nRad = ati[a1].radical;
            while (nRad > RADICAL_TRIPLET)
            {
                nRad -= 2;
            }
            sprintf( szRadicalType, "%d->%d", ati[a1].radical, nRad );
            TREAT_ERR( *err, 0, "Radical center type replaced:" );
            TREAT_ERR( *err, 0, szRadicalType );
            cRadical = nRad;
            if (nRad < 0)
            {
                *err |= 8; /*  Unrecognized Radical replaced with non-radical */
            }
        }
        break;
    }
    at[a1].radical = cRadical;

    /* coordinates */
    at[a1].x = ati[a1].x;
    at[a1].y = ati[a1].y;
    at[a1].z = ati[a1].z;

    if (szCoord)
    {
        /* store text coordinates */
        char str[32];
        MOL_COORD * coord_p = szCoord + a1;
        WriteCoord( str, ati[a1].x );
        memcpy( *coord_p, str, 10 );
        WriteCoord( str, ati[a1].y );
        memcpy( *coord_p + 10, str, 10 );
        WriteCoord( str, ati[a1].z );
        memcpy( *coord_p + 20, str, 10 );
    }

    if (MIN_BOND_LENGTH < fabs( ati[a1].x ) || MIN_BOND_LENGTH < fabs( ati[a1].y ) || MIN_BOND_LENGTH < fabs( ati[a1].z ))
    {
        if (MIN_BOND_LENGTH < fabs( ati[a1].z ))
        {
            *nDim |= 3;
        }
        else
        {
            *nDim |= 2;
        }
    }

    /* orig. at. number */
    at[a1].orig_at_number = a1 + 1;
    return 0;

#undef MIN_BOND_LENGTH
}


/****************************************************************************/
int SetBondProperties( inp_ATOM *at,
                       inchi_Atom *ati,
                       int a1,
                       int j,
                       int nNumAtoms,
                       int *nNumBonds,
                       char *pStrErr,
                       int *err )
{
    int a2;
    S_CHAR     cBondType, cStereoType1, cStereoType2;
    AT_NUMB   *p1, *p2;
    int        n1, n2;

    /* bond type */
    switch (ati[a1].bond_type[j])
    {
        case INCHI_BOND_TYPE_SINGLE:
            cBondType = BOND_TYPE_SINGLE;
            break;
        case INCHI_BOND_TYPE_DOUBLE:
            cBondType = BOND_TYPE_DOUBLE;
            break;
        case INCHI_BOND_TYPE_TRIPLE:
            cBondType = BOND_TYPE_TRIPLE;
            break;
        case INCHI_BOND_TYPE_ALTERN:
            cBondType = BOND_TYPE_ALTERN;
            break;
        default:
        {
            char szBondType[16];
            sprintf( szBondType, "%d", ati[a1].bond_type[j] );
            TREAT_ERR( *err, 0, "Unrecognized bond type:" );
            TREAT_ERR( *err, 0, szBondType );
            *err |= 8; /*  Unrecognized Bond type replaced with single bond */
            cBondType = BOND_TYPE_SINGLE;
        }
        break;
    }

    /* 2D stereo */

    switch (ati[a1].bond_stereo[j])
    {
    /* stereocenter-related; positive: the sharp end points to this atom  */
        case   INCHI_BOND_STEREO_NONE:
            cStereoType1 = 0;
            cStereoType2 = 0;
            break;
        case   INCHI_BOND_STEREO_SINGLE_1UP:
            cStereoType1 = STEREO_SNGL_UP;
            cStereoType2 = -STEREO_SNGL_UP;
            break;
        case   INCHI_BOND_STEREO_SINGLE_1EITHER:
            cStereoType1 = STEREO_SNGL_EITHER;
            cStereoType2 = -STEREO_SNGL_EITHER;
            break;
        case   INCHI_BOND_STEREO_SINGLE_1DOWN:
            cStereoType1 = STEREO_SNGL_DOWN;
            cStereoType2 = -STEREO_SNGL_DOWN;
            break;
        /* stereocenter-related; negative: the sharp end points to the opposite atom  */
        case   INCHI_BOND_STEREO_SINGLE_2UP:
            cStereoType1 = -STEREO_SNGL_UP;
            cStereoType2 = STEREO_SNGL_UP;
            break;
        case   INCHI_BOND_STEREO_SINGLE_2EITHER:
            cStereoType1 = -STEREO_SNGL_EITHER;
            cStereoType2 = STEREO_SNGL_EITHER;
            break;
        case   INCHI_BOND_STEREO_SINGLE_2DOWN:
            cStereoType1 = -STEREO_SNGL_DOWN;
            cStereoType2 = STEREO_SNGL_DOWN;
            break;
        /* stereobond-related */
        case   INCHI_BOND_STEREO_DOUBLE_EITHER:
        case  -INCHI_BOND_STEREO_DOUBLE_EITHER:
            cStereoType1 = STEREO_DBLE_EITHER;
            cStereoType2 = STEREO_DBLE_EITHER;
            break;
        default:
        {
            char szBondType[16];
            sprintf( szBondType, "%d", ati[a1].bond_stereo[j] );
            TREAT_ERR( *err, 0, "Unrecognized bond stereo:" );
            TREAT_ERR( *err, 0, szBondType );
            *err |= 8; /*  Unrecognized Bond stereo replaced with non-stereo bond */
            cStereoType1 = 0;
            cStereoType2 = 0;
        }
        break;
    }

    /* neighbor */
    if (ati[a1].neighbor[j] < 0 || ati[a1].neighbor[j] >= nNumAtoms)
    {
        *err |= 1; /*  bond for impossible atom number(s); ignored */
        TREAT_ERR( *err, 0, "Bond to nonexistent atom" );
        goto err_exit;
    }

    a2 = (AT_NUMB) ati[a1].neighbor[j];
    if (a2 == a1)
    {
        *err |= 1; /*  bond for impossible atom number(s); ignored */
        TREAT_ERR( *err, 0, "Atom has a bond to itself" );
        goto err_exit;
    }

    /* consistency check; locate the bond in the opposite atom */
    p1 = is_in_the_list( at[a1].neighbor, (AT_NUMB) a2, at[a1].valence );
    p2 = is_in_the_list( at[a2].neighbor, (AT_NUMB) a1, at[a2].valence );

    if (p1 && p2)
    {
        n1 = (int) ( p1 - at[a1].neighbor );
        n2 = (int) ( p2 - at[a2].neighbor );
        if ((n1 + 1 < at[a1].valence &&
             is_in_the_list( at[a1].neighbor + n1 + 1, (AT_NUMB) a2, at[a1].valence - n1 - 1 ))
             ||
             (n2 + 1 < at[a2].valence &&
             is_in_the_list( at[a2].neighbor + n2 + 1, (AT_NUMB) a1, at[a2].valence - n2 - 1 ))) /* djb-rwth: addressing LLVM warnings */
        {
            TREAT_ERR( *err, 0, "Multiple bonds between two atoms" );
            *err |= 2; /*  multiple bonds between atoms */
        }
        else if (n1 < at[a1].valence && n2 < at[a2].valence &&
             cBondType == at[a2].bond_type[n2] &&
             cBondType == at[a1].bond_type[n1] &&
             cStereoType1 == at[a1].bond_stereo[n1] &&
             cStereoType2 == at[a2].bond_stereo[n2])
        {
            /*TREAT_ERR (*err, 0, "Duplicated bond(s) between two atoms");*/
        }
        else
        {
            TREAT_ERR( *err, 0, "Multiple bonds between two atoms" );
            *err |= 2; /*  multiple bonds between atoms */
        }
    }
    else if (( p1 || p2 ) &&
        ( p1 || at[a1].valence < MAXVAL ) &&
        ( p2 || at[a2].valence < MAXVAL ))
    {
        n1 = p1 ? (int) ( p1 - at[a1].neighbor ) : at[a1].valence++;
        n2 = p2 ? (int) ( p2 - at[a2].neighbor ) : at[a2].valence++;
        /* the bond is present in one atom only: possibly program error */
        if ((p1 && ( cBondType != at[a1].bond_type[n1] || at[a1].bond_stereo[n1] != cStereoType1 )) ||
             (p2 && ( cBondType != at[a2].bond_type[n2] || at[a2].bond_stereo[n2] != cStereoType2 ))) /* djb-rwth: addressing LLVM warnings */
        {
            TREAT_ERR( *err, 0, "Multiple bonds between two atoms" );
            *err |= 2; /*  multiple bonds between atoms */
        }
        else
        {
            TREAT_ERR( *err, 0, "Duplicated bond(s) between two atoms" );
            /* warning */
        }
    }
    else if (!p1 && !p2 && at[a1].valence < MAXVAL && at[a2].valence < MAXVAL)
    {
        n1 = at[a1].valence++;
        n2 = at[a2].valence++;
        ( *nNumBonds )++;
    }
    else
    {
        char szMsg[64];
        *err |= 4; /*  too large number of bonds. Some bonds ignored. */
        sprintf(szMsg, "Atom '%s' has more than %d bonds",
            at[a1].valence >= MAXVAL ? at[a1].elname : at[a2].elname, MAXVAL);
        TREAT_ERR( *err, 0, szMsg );
        goto err_exit;
    }

    /* store the connection */

    /* bond type */ /* djb-rwth: fixing buffer overruns */
    if ((n1 < MAXVAL) && (n2 < MAXVAL))
    {
        at[a1].bond_type[n1] =
            at[a2].bond_type[n2] = cBondType;
        /* connection */
        at[a1].neighbor[n1] = (AT_NUMB)a2;
        at[a2].neighbor[n2] = (AT_NUMB)a1;
        /* stereo */
        at[a1].bond_stereo[n1] = cStereoType1; /*  >0: the wedge (pointed) end is at this atom */
        at[a2].bond_stereo[n2] = cStereoType2; /*  <0: the wedge (pointed) end is at the opposite atom */
    }
    else
    {
        goto err_exit;
    }

    return 0;

err_exit:

    return 1;
}


/****************************************************************************/
int SetAtomAndBondProperties( inp_ATOM *at,
                              inchi_Atom *ati,
                              int a1,
                              int bDoNotAddH,
                              char *pStrErr,
                              int *err )
{
    int valence, chem_valence, num_alt_bonds, j, n1;
    int nRadical, nCharge;
    static int el_number_H = 0;

    if (!el_number_H)
    {
        el_number_H = EL_NUMBER_H;
    }

    nRadical = nCharge = 0;
    valence = at[a1].valence;
    chem_valence = num_alt_bonds = 0;
    for (j = 0; j < valence; j++)
    {
        if (at[a1].bond_type[j] <= BOND_TYPE_TRIPLE)
        {
            chem_valence += at[a1].bond_type[j];
        }
        else
        {
            num_alt_bonds++;
        }
    }
    switch (num_alt_bonds)
    {
        case 0:
            break;
        case 2:
            chem_valence += 3; /* -C= */
            break;
        case 3:
            chem_valence += 4;  /* >C= */
            break;
        default:
        {
            char szMsg[64];
            *err |= 8; /*  wrong number of alt. bonds */
            sprintf(szMsg, "Atom '%s' has %d alternating bonds",
                at[a1].elname, num_alt_bonds);
            TREAT_ERR( *err, 0, szMsg );
        }
        break;
    }
    at[a1].chem_bonds_valence = chem_valence;

    /* aliased hydrogen atoms */
    if (ERR_ELEM == ( n1 = get_periodic_table_number( at[a1].elname ) ))
    {
        /*  Case when elname contains more than 1 element: extract number of H if possible */
        if (extract_charges_and_radicals( at[a1].elname, &nRadical, &nCharge ))
        {
            if ((nRadical && at[a1].radical && nRadical != at[a1].radical) ||
                 (nCharge && at[a1].charge && nCharge != at[a1].charge)) /* djb-rwth: addressing LLVM warnings */
            {
                TREAT_ERR( *err, 0, "Ignored charge/radical redefinition:" );
                TREAT_ERR( *err, 0, ati[a1].elname );
            }
            else
            {
                if (nRadical)
                {
                    at[a1].radical = nRadical;
                }
                if (nCharge)
                {
                    at[a1].charge = nCharge;
                }
            }
        }

        at[a1].num_H = extract_H_atoms( at[a1].elname, at[a1].num_iso_H );
        if (!at[a1].elname[0] && NUMH( at, a1 ))
        {
            /* alias contains only H. Added 2004-07-21, fixed 2004-07-22
             * move the heaviest isotope to the "central atom"
             * Note: this must be consistent with H-H treatment in remove_terminal_HDT()
             */
            strcpy( at[a1].elname, "H" );
            if (NUM_ISO_H( at, a1 ))
            {
                for (j = NUM_H_ISOTOPES - 1; 0 <= j; j--)
                {
                    if (at[a1].num_iso_H[j])
                    {
                        at[a1].num_iso_H[j] --;
                        at[a1].iso_atw_diff = 1 + j;
                        break;
                    }
                }
            }
            else
            {
                at[a1].num_H--;
            }
        }

        if (ERR_ELEM == ( n1 = get_periodic_table_number( at[a1].elname ) ))
        {
            n1 = 0;
        }
        if (n1)
        {
            at[a1].at_type |= 1; /* "Aliased" atom: data in the element name */
            TREAT_ERR( *err, 0, "Parsed compound atom(s):" );
            TREAT_ERR( *err, 0, ati[a1].elname );
        }
    }

    at[a1].el_number = (U_CHAR) n1;
    if (!n1)
    {
        *err |= 64; /*  Unrecognized aromatic bond(s) replaced with single */
        TREAT_ERR( *err, 0, "Unknown element(s):" );
        TREAT_ERR( *err, 0, at[a1].elname );
    }
    else
    {
        /* replace explicit D or T with isotopic H (added 2003-06-02) */
        if (el_number_H == n1 && !at[a1].iso_atw_diff)
        {
            switch (at[a1].elname[0])
            {
                case 'D':
                    at[a1].iso_atw_diff = 2;
                    mystrncpy( at[a1].elname, "H", sizeof( at->elname ) );
                    break;
                case 'T':
                    at[a1].iso_atw_diff = 3;
                    mystrncpy( at[a1].elname, "H", sizeof( at->elname ) );
                    break;
                case 'H':
                    if (1 <= ati[a1].isotopic_mass)
                    {
                        AT_NUM iso_atw_diff;
                        if (ISOTOPIC_SHIFT_FLAG - ISOTOPIC_SHIFT_MAX <= ati[a1].isotopic_mass &&
                             ISOTOPIC_SHIFT_FLAG + ISOTOPIC_SHIFT_MAX >= ati[a1].isotopic_mass)
                        {
                            /* ati[a1].isotopic_mass is isotopic iso_atw_diff + ISOTOPIC_SHIFT_FLAG */
                            iso_atw_diff = ati[a1].isotopic_mass - ISOTOPIC_SHIFT_FLAG;
                        }
                        else
                        {
                            /* ati[a1].isotopic_mass is isotopic mass */
                            int iso_atw = get_atomic_mass_from_elnum( (int) at[a1].el_number );
                            iso_atw_diff = ati[a1].isotopic_mass - iso_atw;
                        }
                        if (iso_atw_diff >= 0)
                            iso_atw_diff++;
                        /* reproduce Bug04: allowed non-terminal H heavier than T */
                        if (1 <= iso_atw_diff &&
                            ( at[a1].valence != 1 || iso_atw_diff <= NUM_H_ISOTOPES ))
                        {
                            at[a1].iso_atw_diff = (S_CHAR) iso_atw_diff;
                        }
                    }
            }
        }
        else
        {/* isotopic shift */
            if (ati[a1].isotopic_mass)
            {
                AT_NUM iso_atw_diff;
                if (ISOTOPIC_SHIFT_FLAG - ISOTOPIC_SHIFT_MAX <= ati[a1].isotopic_mass &&
                     ISOTOPIC_SHIFT_FLAG + ISOTOPIC_SHIFT_MAX >= ati[a1].isotopic_mass)
                {
                    /* ati[a1].isotopic_mass is isotopic iso_atw_diff + ISOTOPIC_SHIFT_FLAG */
                    iso_atw_diff = ati[a1].isotopic_mass - ISOTOPIC_SHIFT_FLAG;
                }
                else
                {
                    /* ati[a1].isotopic_mass is isotopic mass */
                    iso_atw_diff = get_atomic_mass_from_elnum( (int) at[a1].el_number );
                    iso_atw_diff = ati[a1].isotopic_mass - iso_atw_diff;
                }
                if (iso_atw_diff >= 0)
                    iso_atw_diff++;
                at[a1].iso_atw_diff = (S_CHAR) iso_atw_diff;
            }
        }
    }

    /* add implicit hydrogen atoms flag */
    if (ati[a1].num_iso_H[0] == -1)
    {
        if (!bDoNotAddH)
        {
            at[a1].at_type |= 2; /* user requested to add H */
        }
    }
    else
    {
        at[a1].num_H = ati[a1].num_iso_H[0];
    }

    for (j = 0; j < NUM_H_ISOTOPES; j++)
    {
        at[a1].num_iso_H[j] = ati[a1].num_iso_H[j + 1];
    }

    if (num_alt_bonds)
    {
        /* atom has aromatic bonds AND the chemical valence is not known */
        int num_H = NUMH( at, a1 );
        int chem_valence_alt = at[a1].chem_bonds_valence + num_H;
        int bUnusualValenceArom =
            detect_unusual_el_valence( (int) at[a1].el_number, at[a1].charge,
                                        at[a1].radical, chem_valence_alt,
                                        num_H, at[a1].valence );
        int bUnusualValenceNoArom =
            detect_unusual_el_valence( (int) at[a1].el_number, at[a1].charge,
                                        at[a1].radical, chem_valence_alt - 1,
                                        num_H, at[a1].valence );
        if (bUnusualValenceArom && !bUnusualValenceNoArom && 0 == nBondsValToMetal( at, a1 ))
        {
            /* typically NH in 5-member aromatic ring */
            at[a1].chem_bonds_valence--;
        }
    }

    return 0;
}


/****************************************************************************/
int InpAtom0DToInchiAtom( inp_ATOM *at,
                          int num_inp_atoms,
                          AT_NUM *num_atoms,
                          inchi_Atom **atom,
                          AT_NUM *num_stereo0D,
                          inchi_Stereo0D **stereo0D )
{
    int num_stereo_centers, num_stereo_bonds, num_inp_stereo0D, i, m, m1, m2, n, ret = 0;

    /* count stereobonds, allenes. cumulenes. and stereoatoms */
    num_stereo_centers = num_stereo_bonds = ret = 0;

    *atom = NULL;
    *num_atoms = 0;
    *stereo0D = NULL;
    *num_stereo0D = 0;

    for (i = 0; i < num_inp_atoms; i++)
    {
        if (at[i].p_parity)
        {
            /* stereocenter */
            num_stereo_centers++;
        }
        else
        {
            for (m = 0; m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m]; m++)
            {
                ;
            }
            num_stereo_bonds += m;
        }
    }

    num_stereo_bonds /= 2;
    num_inp_stereo0D = num_stereo_bonds + num_stereo_centers;

    if (num_inp_atoms > 0)
    {
        *atom = (inchi_Atom *) inchi_calloc( num_inp_atoms, sizeof( ( *atom )[0] ) );
    }

    *num_atoms = num_inp_atoms;

    if (num_inp_stereo0D > 0)
    {
        *stereo0D = (inchi_Stereo0D *) inchi_calloc( num_inp_stereo0D, sizeof( ( *stereo0D )[0] ) );
    }

    if ((num_inp_atoms && !( *atom )) || (num_inp_stereo0D > 0 && !( *stereo0D ))) /* djb-rwth: addressing LLVM warnings */
    {
        /* allocation failed */
        ret = -1;
        goto exit_function;
    }

    /* copy atom properties */
    for (i = 0; i < num_inp_atoms; i++)
    {
        ( *atom )[i].num_bonds = at[i].valence;
        for (m = 0; m < at[i].valence; m++)
        {
            ( *atom )[i].bond_type[m] = at[i].bond_type[m];
            ( *atom )[i].neighbor[m] = at[i].neighbor[m];
        }
        ( *atom )[i].charge = at[i].charge;
#if USE_BCF
        memcpy_s( ( *atom )[i].elname, ATOM_EL_LEN, at[i].elname, ATOM_EL_LEN ); /* djb-rwth: function replaced with its safe C11 variant */
#else
        memcpy( ( *atom )[i].elname, at[i].elname, ATOM_EL_LEN );
#endif
        if (at[i].iso_atw_diff)
        {
            ( *atom )[i].isotopic_mass = ISOTOPIC_SHIFT_FLAG + ( at[i].iso_atw_diff > 0 ? at[i].iso_atw_diff - 1 : at[i].iso_atw_diff );
        }
        ( *atom )[i].num_iso_H[0] = at[i].num_H;
        for (m = 0; m < NUM_H_ISOTOPES; m++)
        {
            ( *atom )[i].num_iso_H[m + 1] = at[i].num_iso_H[m];
        }
        ( *atom )[i].radical = at[i].radical;
    }

    /* stereo */
    for (i = n = 0; i < num_inp_atoms; i++)
    {
        if (at[i].p_parity)
        {
            if (n < num_inp_stereo0D)
            {
                ( *stereo0D )[n].central_atom = i;
                ( *stereo0D )[n].parity = at[i].p_parity;
                ( *stereo0D )[n].type = INCHI_StereoType_Tetrahedral;
                for (m = 0; m < MAX_NUM_STEREO_ATOM_NEIGH; m++)
                {
                    ( *stereo0D )[n].neighbor[m] = at[i].p_orig_at_num[m] - 1;
                }
                n++;
            }
            else
            {
                ret |= 1;
                break;
            }
        }
        else
        {
            for (m1 = 0; m1 < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m1]; m1++)
            {
                /* find the opposite atom at the other end of double bond, allene, or cumulene */
                int chain[12], len = 0, nxt_neigh, nxt, cur;
                cur = chain[len++] = i;
                nxt_neigh = at[cur].sb_ord[m1];

                do
                {
                    /* add next atom */
                    chain[len++] = nxt = at[cur].neighbor[nxt_neigh];
                    nxt_neigh = ( at[nxt].neighbor[0] == cur );
                    cur = nxt;
                    /* find nxt_neigh */
                }
                while (!at[cur].sb_parity[0] &&
                       len < 12 &&
                       at[cur].valence == 2);

                if (at[cur].sb_parity[0] && len <= 4 && i < cur /* count bonds only one time */)
                {
                    /* double bond, cumulene, or allene has been found */
                    for (m2 = 0; m2 < MAX_NUM_STEREO_BONDS && at[cur].sb_parity[m2]; m2++)
                    {
                        if (chain[len - 2] == at[cur].neighbor[(int) at[cur].sb_ord[m2]])
                        {
                            if (n < num_inp_stereo0D)
                            {
                                int parity1 = at[i].sb_parity[m1];
                                int parity2 = at[cur].sb_parity[m2];
                                int parity;
                                if (( INCHI_PARITY_ODD == parity1 || INCHI_PARITY_EVEN == parity1 ) &&
                                    ( INCHI_PARITY_ODD == parity2 || INCHI_PARITY_EVEN == parity2 ))
                                {
                                    /* well-defined parity */
                                    parity = ( parity1 == parity2 ) ? INCHI_PARITY_EVEN : INCHI_PARITY_ODD;
                                }
                                else
                                {
                                    parity = inchi_max( parity1, parity2 );
                                }
                                ( *stereo0D )[n].central_atom = ( len == 3 ) ? chain[1] : NO_ATOM;
                                ( *stereo0D )[n].parity = parity;
                                ( *stereo0D )[n].type = len == 3 ? INCHI_StereoType_Allene : INCHI_StereoType_DoubleBond;
                                ( *stereo0D )[n].neighbor[0] = at[i].sn_orig_at_num[m1] - 1;
                                ( *stereo0D )[n].neighbor[1] = i;
                                ( *stereo0D )[n].neighbor[2] = cur;
                                ( *stereo0D )[n].neighbor[3] = at[cur].sn_orig_at_num[m2] - 1;
                                n++;
                            }
                            else
                            {
                                ret |= 1;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    *num_stereo0D = n;

exit_function:
    if (ret < 0)
    {
        if (*atom)
        {
            inchi_free( *atom );
        }
        if (*stereo0D)
        {
            inchi_free( *stereo0D );
        }
        *atom = NULL;
        *stereo0D = NULL;
        *num_atoms = 0;
        *num_stereo0D = 0;
    }

    return ret;
}


/****************************************************************************/
int ExtractOneStructure( STRUCT_DATA *sd,
                         INPUT_PARMS *ip,
                         char *szTitle,
                         inchi_InputEx *inp,
                         INCHI_IOSTREAM *log_file,
                         INCHI_IOSTREAM *out_file,
                         INCHI_IOSTREAM *prb_file,
                         ORIG_ATOM_DATA *orig_inp_data,
                         long *num_inp )
{
    int         *err = &sd->nStructReadError;
    char        *pStrErr = sd->pStrErrStruct;
    inp_ATOM    *at = NULL;
    MOL_COORD   *szCoord = NULL;
    inchi_Atom  *ati = NULL;
    int       nNumAtoms = 0;
    int       a1, j, valence, nDim, nNumBonds, nRet = 0, max_num_at;

    /* vABParityUnknown holds actual value of an internal constant signifying       */
    /* unknown parity: either the same as for undefined parity (default==standard)  */
    /*  or a specific one (non-std; requested by SLUUD switch).                     */
    int vABParityUnknown = AB_PARITY_UNDF;
    if (0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO ))
    {
        /* Make labels for unknown and undefined stereo different */
        vABParityUnknown = AB_PARITY_UNKN;
    }

    /********************************************************
     *
     *   Extract the structure
     *
     ********************************************************/

    FreeOrigAtData( orig_inp_data );
    nDim = 0;
    nNumBonds = 0;

    if (!inp || ( nNumAtoms = inp->num_atoms ) <= 0 || !( ati = inp->atom ))
    {
        TREAT_ERR( *err, 0, "Empty structure" );
        *err = 98;
        goto err_exit;
    }

    max_num_at = ip->bLargeMolecules ? MAX_ATOMS : NORMALLY_ALLOWED_INP_MAX_ATOMS;
    if (nNumAtoms >= max_num_at)
    {
        TREAT_ERR( *err, 0, "Too many atoms [did you forget 'LargeMolecules' switch?]" );
        *err = 70;
        orig_inp_data->num_inp_atoms = -1;
        goto err_exit;
    }

    at = (inp_ATOM  *) inchi_calloc( nNumAtoms, sizeof( at[0] ) );
    szCoord = (MOL_COORD *) inchi_calloc( inchi_max( nNumAtoms, 1 ), sizeof( MOL_COORD ) );

    if (!at || !szCoord)
    {
        TREAT_ERR( *err, 0, "Out of RAM" );
        *err = -1;
        goto err_exit;
    }

    /********************************************************
     *
     *   Extract typical for Molfile structural data
     *
     ********************************************************/
    /* extract atoms and bonds */
    for (a1 = 0; a1 < nNumAtoms; a1++)
    {
        /* extract atoms */
        SetAtomProperties( at, szCoord, ati, a1, &nDim, pStrErr, err );

        if (*err)
        {
            goto err_exit;
        }

        /* extract connections */
        valence = ati[a1].num_bonds;
        for (j = 0; j < valence; j++)
        {
            SetBondProperties( at, ati, a1, j, nNumAtoms, &nNumBonds, pStrErr, err );
        }

        if (*err)
        {
            goto err_exit;
        }
    }

    orig_inp_data->num_inp_atoms = nNumAtoms;
    orig_inp_data->num_inp_bonds = nNumBonds;
    orig_inp_data->num_dimensions = nDim;

    /* extract elements, chemical valences, implicit H, isotopic shifts */
    for (a1 = 0; a1 < nNumAtoms; a1++)
    {
        /* set temp flags in at[a1].at_type */
        /* (1: data in atom name; 2: request to add H) */
        SetAtomAndBondProperties( at,
                                  ati,
                                  a1,
                                  ip->bDoNotAddH,
                                  pStrErr,
                                  err );
        if (*err)
        {
            goto err_exit;
        }
    }

    /* clear temp flags in at[].at_type; add implicit H */
    SetNumImplicitH( at, nNumAtoms );

    if (*err)
    {
        goto err_exit;
    }

    /********************************************************
     *
     *   Extract the 0D parities (typical for CML)
     *
     ********************************************************/
    Extract0DParities( at,
                       nNumAtoms,
                       inp->stereo0D,
                       inp->num_stereo0D,
                       pStrErr,
                       err,
                       vABParityUnknown );

    if (*err)
    {
        goto err_exit;
    }

    orig_inp_data->at = at;
    at = NULL;
    orig_inp_data->num_dimensions = nDim;
    orig_inp_data->num_inp_atoms = nNumAtoms;
    orig_inp_data->num_inp_bonds = nNumBonds;
    orig_inp_data->szCoord = szCoord;
    szCoord = NULL;

    /* chiral flag */

    /* *****************************************************************************
     * Chiral flags are set in:
     * - ReadTheStructure() inchi-1, wInChI
     * - e_IchiMain.c -- main()               -- C example of calling InChI dll
     * - inchi_dll.c  ExtractOneStructure -- InChI dll code (here)
     *******************************************************************************/

    if (( ip->nMode & REQ_MODE_CHIR_FLG_STEREO ) && ( ip->nMode & REQ_MODE_STEREO ))
    {
        if (ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL)
        {
            /* absolute stereo */
            ip->nMode &= ~( REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO );
            sd->bChiralFlag &= ~FLAG_INP_AT_NONCHIRAL;
            sd->bChiralFlag |= FLAG_INP_AT_CHIRAL; /* write AuxInfo as chiral */
        }
        else
        /*if ( ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL )*/
        {
            /* relative stereo */
            ip->nMode &= ~( REQ_MODE_RACEMIC_STEREO );
            ip->nMode |= REQ_MODE_RELATIVE_STEREO;
            sd->bChiralFlag &= ~FLAG_INP_AT_CHIRAL;
            sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL; /* write AuxInfo as non-chiral */
        }
    }
    else if (ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL)
    {
        sd->bChiralFlag &= ~FLAG_INP_AT_NONCHIRAL;
        sd->bChiralFlag |= FLAG_INP_AT_CHIRAL; /* write AuxInfo as chiral */
    }
    else if (ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL)
    {
        sd->bChiralFlag &= ~FLAG_INP_AT_CHIRAL;
        sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL; /* write AuxInfo as non-chiral */
    }

    /* v. 1.05 extensions  */
    {
        int res = SetExtOrigAtDataByInChIExtInput( &orig_inp_data->polymer,
                                                   &orig_inp_data->v3000,
                                                   inp->polymer,
                                                   inp->v3000,
                                                   orig_inp_data->num_inp_atoms );
        if (res)
        {
            TREAT_ERR( res, 0, "General error on treating polymers" );
            *err = -1;
            goto err_exit;
        }
    }
    *num_inp += 1;

err_exit:

    if (at)
    {   /* if not moved to orig_inp_data/then nullified */
        inchi_free( at );
    }
    if (szCoord)
    {
        inchi_free( szCoord );
    }

    nRet = TreatErrorsInReadTheStructure( sd, ip, LOG_MASK_NO_WARN, NULL,
                                          log_file, out_file, prb_file,
                                          orig_inp_data, num_inp );

    return nRet;
}


/****************************************************************************/
int INCHI_DECL GetStringLength( char *p )
{
    if (p)
    {
        return (int) strlen( p );
    }
    else
    {
        return 0;
    }
}

#define MAX_MSG_LEN 512


/****************************************************************************
 GetINCHIfromINCHI does same as -InChI2InChI option: converts InChI into
 InChI for validation purposes
 It may also be used to filter out specific layers. For instance,
 /Snon would remove stereochemical layer
 Omitting /FixedH and/or /RecMet would remove Fixed-H or Reconnected layers
 To keep all InChI layers use options string "/FixedH /RecMet";
 option /InChI2InChI is not needed
 inchi_InputINCHI is created by the user;
 strings in inchi_Output are allocated and deallocated by InChI
 inchi_Output does not need to be initilized out to zeroes;
 see FreeINCHI() on how to deallocate it
****************************************************************************/
int INCHI_DECL GetINCHIfromINCHI( inchi_InputINCHI *inpInChI,
                                  inchi_Output *out )
{
    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;

    static char szMainOption[] = " ?InChI2InChI";

    INCHI_CLOCK ic;
    CANON_GLOBALS CG;

    int i;
    char      szSdfDataValue[MAX_SDF_VALUE + 1];
    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */

    INPUT_PARMS inp_parms;
    INPUT_PARMS *ip = &inp_parms;

    int             bReleaseVersion = bRELEASE_VERSION;
    int   nRet = 0, nRet1;

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif

    const char *argv[INCHI_MAX_NUM_ARG + 1];
    int   argc;
    char *szOptions = NULL;

    INCHI_IOSTREAM inchi_file[3], *out_file = inchi_file, *log_file = inchi_file + 1, *input_file = inchi_file + 2;

#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag( _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF );
/* for execution outside the VC++ debugger uncomment one of the following two */

#ifdef MY_REPORT_FILE
    _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode( _CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG );
#endif

    /* turn on floating point exceptions */
#if ( !defined(__STDC__) || __STDC__ != 1 )
    {
        /* Get the default control word. */
        int cw = _controlfp( 0, 0 );

        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &= ~( EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL );

        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif
#endif

    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */
#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
repeat:
    FreeINCHI( out );
    inchi_ios_close( out_file );
    inchi_ios_close( log_file );
    inchi_ios_reset( input_file );  /* do not close input_file - its string buffer may point to inpInChI->szInChI */
#endif

    /* Initialize internal for this function I/O streams as string buffers */
    inchi_ios_init( input_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( out_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( log_file, INCHI_IOS_TYPE_STRING, NULL );

    sd->bUserQuit = 0;

    /* clear original input structure */
    /* memset( inchi_file, 0, sizeof(inchi_file) ); */
    memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( ip, 0, sizeof( *ip ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( szSdfDataValue, 0, sizeof( szSdfDataValue ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    memset( &ic, 0, sizeof( ic ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &CG, 0, sizeof( CG ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    szMainOption[1] = INCHI_OPTION_PREFX;

    if (!inpInChI)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }

    /* options */
    if (inpInChI)
    {
        int opt_len = (int) ( ( inpInChI->szOptions ? strlen( inpInChI->szOptions ) : 0 ) + sizeof( szMainOption ) + 1 );
        szOptions = (char*) inchi_calloc( (long long)opt_len + 1, sizeof( szOptions[0] ) ); /* djb-rwth: cast operator added */
        if (szOptions)
        {
            if (inpInChI->szOptions)
            {
                strcpy( szOptions, inpInChI->szOptions );
            }
            strcat( szOptions, szMainOption );
            argc = parse_options_string( szOptions, argv, INCHI_MAX_NUM_ARG );
        }
        else
        {
            nRet = _IS_FATAL;
            goto translate_RetVal; /* emergency exit */
        }
    }
    else
    {
        argc = 1;
        argv[0] = "";
        argv[1] = NULL;
    }

    if ((argc == 1
#ifdef TARGET_API_LIB
        && ( !inpInChI || !inpInChI->szInChI ))
#endif

        || (argc == 2 && ( argv[1][0] == INCHI_OPTION_PREFX ) &&
        ( !strcmp( argv[1] + 1, "?" ) || !inchi_stricmp( argv[1] + 1, "help" ) ))) /* djb-rwth: addressing LLVM warnings */
    {
        HelpCommandLineParms( log_file );
        out->szLog = log_file->s.pStr;
        memset( log_file, 0, sizeof( *log_file ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        nRet = _IS_EOF;
        inchi_free(szOptions); /* djb-rwth: avoiding memory leak */
        goto translate_RetVal;
    }

    nRet1 = ReadCommandLineParms( argc, argv, ip, szSdfDataValue,
                                &ulDisplTime, bReleaseVersion, log_file );

    if (szOptions)
    {
        /* argv pointed to strings in szOptions */
        inchi_free( szOptions );
        szOptions = NULL;
    }
    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if (0 > nRet1)
    {
        goto exit_function;
    }

    if (ip->bNoStructLabels)
    {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    }
    else if (ip->nInputType == INPUT_INCHI_XML ||
            ip->nInputType == INPUT_INCHI_PLAIN ||
            ip->nInputType == INPUT_CMLFILE ||
            ip->nInputType == INPUT_INCHI)
    {
        /* the input may contain both the header and the label of the structure */
        if (!ip->pSdfLabel)
        {
            ip->pSdfLabel = ip->szSdfDataHeader;
        }
        if (!ip->pSdfValue)
        {
            ip->pSdfValue = szSdfDataValue;
        }
    }

    if (ip->nInputType && ip->nInputType != INPUT_INCHI)
    {
        inchi_ios_eprint( log_file, "Input type set to INPUT_INCHI\n" );
        ip->nInputType = INPUT_INCHI;
    }

    if (!inpInChI->szInChI)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }
    else
    {
        const int strict = 0;
        nRet = CheckINCHI( inpInChI->szInChI, strict );
        if (nRet != INCHI_VALID_STANDARD     &&
            nRet != INCHI_VALID_NON_STANDARD &&
            nRet != INCHI_VALID_BETA)
        {
            nRet = _IS_ERROR;
            goto exit_function;
        }
    }


    PrintInputParms( log_file, ip );

    /********************************/
    /* InChI -> InChI               */
    /********************************/

    /* input_file simulation */
    input_file->s.pStr = inpInChI->szInChI;
    input_file->s.nUsedLength = (int) strlen( input_file->s.pStr ) + 1;
    input_file->s.nAllocatedLength = input_file->s.nUsedLength;
    input_file->s.nPtr = 0;

    /* buffer for the message */
    out->szMessage = (char *) inchi_calloc( MAX_MSG_LEN, sizeof( out->szMessage[0] ) );
    if (!out->szMessage)
    {
        inchi_ios_eprint( log_file, "Cannot allocate output message buffer.\n" );
        nRet = -1;
    }
    else
    {
        nRet = ReadWriteInChI( &ic, &CG, input_file, out_file, log_file,
                                ip, sd,
                                NULL, 0, NULL,
                                NULL, NULL,
                                out->szMessage, MAX_MSG_LEN,
                                NULL /*out->WarningFlags*/ );
    }

    if (nRet >= 0 && out_file->s.pStr)
    {
        /* success */
        char* p = NULL;
        /* djb-rwth: fixing oss-fuzz issue #40971 */
        int p_len, out_szinchi_len;
        out_szinchi_len = strlen(out_file->s.pStr);
        out->szInChI = out_file->s.pStr;
        out->szAuxInfo = NULL;

        for (p = strchr( out->szInChI, '\n' ); p; p = strchr( p + 1, '\n' ))
        {
            p_len = strlen(p);
            if ((p_len >= 8) && !memcmp( p, "\nAuxInfo", 8 ))
            {
                *p = '\0';            /* remove LF after INChI */
                out->szAuxInfo = p + 1; /* save pointer to AuxInfo */
            }
            else if (out->szAuxInfo || !p[1])
            {
                /* remove LF after aux info or from the last char */
                *p = '\0';
                break;
            }
        }
        out_file->s.pStr = NULL;
    }

    /*
    out->szLog = log_file->pStr;
    log_file->pStr   = NULL;
    */

exit_function:;

    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (ip->path[i])
        {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( &CG );

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if (num_repeat-- > 0)
    {
        goto repeat;
    }
#endif

#ifdef TARGET_API_LIB
    /* output */

        if (log_file->s.pStr && log_file->s.nUsedLength > 0)
        {
            while (log_file->s.nUsedLength && '\n' == log_file->s.pStr[log_file->s.nUsedLength - 1])
            {
                log_file->s.pStr[--log_file->s.nUsedLength] = '\0'; /* remove last LF */
            }
            if (out)
            {
                out->szLog = log_file->s.pStr;
                log_file->s.pStr = NULL;
            }
        }

#endif

translate_RetVal:

    /* close internal output streams */
    inchi_ios_close( out_file );
    inchi_ios_close( log_file );
    inchi_ios_reset( input_file );  /* do not close input_file - its string buffer may point to inpInChI->szInChI */

    switch (nRet)
    {
        case -3: nRet = inchi_Ret_ERROR; break; /* Error: no Structure has been created */
        case -2: nRet = inchi_Ret_ERROR; break; /* Error: no Structure has been created */
        case -1: nRet = inchi_Ret_FATAL; break; /* Severe error: no Structure has been created (typically; break; memory allocation failed) */
        default:
            /*
            if ( !outStruct->atom || !outStruct->num_atoms )
            {
                nRet = inchi_Ret_EOF;
            }
            else
            {
                int m,n,t=0;
                for ( m=0; m < 2; m ++ )
                {
                    for ( n=0; n < 2; n ++ )
                    {
                        if ( outStruct->WarningFlags[m][n] ) {
                            t ++;
                        }
                    }
                }
                nRet = t? inchi_Ret_WARNING : inchi_Ret_OKAY;
            }
            */
            break;
    }

    return nRet;
}


/****************************************************************************

    GetStructFromStdINCHI

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
int INCHI_DECL GetStructFromStdINCHI( inchi_InputINCHI *inpInChI,
                                      inchi_OutputStruct *outStruct )
{
    if (( inpInChI ) &&
        ( inpInChI->szInChI ) &&
        ( strlen( inpInChI->szInChI ) >= LEN_INCHI_STRING_PREFIX + 3 ) &&
        ( inpInChI->szInChI[LEN_INCHI_STRING_PREFIX + 1] == 'S' ))
    {
         /* brief check indicated valid std input (more checks in GetStructFromINCHI) */
        return GetStructFromINCHI( inpInChI, outStruct );
    }
    else
    {
        /* non-std or just invalid input */
        return inchi_Ret_ERROR;
    }
}


/****************************************************************************

    GetStructFromINCHIEx
****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
int INCHI_DECL GetStructFromINCHIEx( inchi_InputINCHI *inpInChI,
                                     inchi_OutputStructEx *outStruct )
{
    INCHI_CLOCK ic;
    CANON_GLOBALS CG;
    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;
    INPUT_PARMS inp_parms;
    INPUT_PARMS *ip = &inp_parms;
    INCHI_IOSTREAM inchi_file[3];
    INCHI_IOSTREAM *out_file = inchi_file, *log_file = inchi_file + 1, *input_file = inchi_file + 2;
    int    i, nRet = 0, nRet1;
    /* djb-rwth: removing redundant variables/code */
    int bReleaseVersion = bRELEASE_VERSION; /* djb-rwth: ignoring LLVM warning: variable used in function call */
    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */
#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif
    static char szMainOption[] = " ?InChI2Struct";
    char szSdfDataValue[MAX_SDF_VALUE + 1];
    const char *argv[INCHI_MAX_NUM_ARG + 1];
    int   argc;
    char *szOptions = NULL;
    /* conversion result */
    inp_ATOM *at = NULL;
    int num_at = 0;
    OAD_Polymer *polymer = NULL;
    OAD_V3000    *v3000 = NULL;

#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag( _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF );

/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE
    _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
    _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
    _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode( _CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG );
#endif

    /* turn on floating point exceptions */
#if ( !defined(__STDC__) || __STDC__ != 1 )
    {
        /* Get the default control word. */
        int cw = _controlfp( 0, 0 );

        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &= ~( EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL );

        /* Set the control word. */
        _controlfp( cw, MCW_EM );
    }
#endif
#endif

    memset( outStruct, 0, sizeof( *outStruct ) ); /* djb-rwth: memset_s C11/Annex K variant? */

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    repeat:
    FreeStructFromINCHI( &outStruct );
    inchi_ios_reset( input_file );  /* do not close input_file - its string buffer may point to inpInChI->szInChI */
    inchi_ios_close( out_file );
    inchi_ios_close( log_file );
#endif

    sd->bUserQuit = 0;

    /* Initialize internal for this function I/O streams as string buffers */
    inchi_ios_init( input_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( out_file, INCHI_IOS_TYPE_STRING, NULL );
    inchi_ios_init( log_file, INCHI_IOS_TYPE_STRING, NULL );

    /* clear original input structure */
    memset( sd, 0, sizeof( *sd ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( ip, 0, sizeof( *ip ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( szSdfDataValue, 0, sizeof( szSdfDataValue ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    memset( &ic, 0, sizeof( ic ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    memset( &CG, 0, sizeof( CG ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    szMainOption[1] = INCHI_OPTION_PREFX;

    if (!inpInChI)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }

    /* options */
    if (inpInChI /*&& inpInChI->szOptions*/)
    {
        /* fix bug discovered by Burt Leland 2008-12-23 */
        int opt_len = ( inpInChI->szOptions ? strlen( inpInChI->szOptions ) : 0 ) + sizeof( szMainOption ) + 1;
        szOptions = (char*)inchi_calloc((long long)opt_len + 1, sizeof(szOptions[0])); /* djb-rwth: cast operator added */
        if (szOptions)
        {
            if (inpInChI->szOptions)
                /* fix bug discovered by Burt Leland 2008-12-23 */
                strcpy(szOptions, inpInChI->szOptions);
            strcat(szOptions, szMainOption);
            argc = parse_options_string( szOptions, argv, INCHI_MAX_NUM_ARG );
        }
        else
        {
            nRet = _IS_FATAL;
            goto translate_RetVal; /* emergency exit */
        }
    }
    else
    {
        argc = 1;
            argv[0] = "";
        argv[1] = NULL;
    }

    if ((argc == 1
#ifdef TARGET_API_LIB
        && ( !inpInChI || !inpInChI->szInChI ))
#endif
        || (argc == 2 && ( argv[1][0] == INCHI_OPTION_PREFX ) &&
        ( !strcmp( argv[1] + 1, "?" ) || !inchi_stricmp( argv[1] + 1, "help" ) ))) /* djb-rwth: addressing LLVM warnings */
    {
        HelpCommandLineParms( log_file );
        outStruct->szLog = log_file->s.pStr;
        nRet = _IS_EOF;
        goto translate_RetVal;
    }

    nRet1 = ReadCommandLineParms( argc, argv, ip, szSdfDataValue,
                                  &ulDisplTime, bReleaseVersion,
                                  log_file );

    if (szOptions)
    {
        /* argv pointed to strings in szOptions */
        inchi_free( szOptions );
        szOptions = NULL;
    }

    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if (0 > nRet1)
    {
        goto exit_function;
    }

    if (ip->bNoStructLabels)
    {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    }
    else if (ip->nInputType == INPUT_INCHI_XML ||
            ip->nInputType == INPUT_INCHI_PLAIN ||
            ip->nInputType == INPUT_CMLFILE ||
            ip->nInputType == INPUT_INCHI)
    {
        /* the input may contain both the header and the label of the structure */
        if (!ip->pSdfLabel)
            ip->pSdfLabel = ip->szSdfDataHeader;
        if (!ip->pSdfValue)
            ip->pSdfValue = szSdfDataValue;
    }

    if (ip->nInputType && ip->nInputType != INPUT_INCHI)
    {
        inchi_ios_eprint( log_file, "Input type set to INPUT_INCHI\n" );
        ip->nInputType = INPUT_INCHI;
    }

    if (!inpInChI->szInChI)
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }
    else
    {
        const int strict = 0;                     /* do not use strict mode, it may be too alarmous */
        nRet = CheckINCHI( inpInChI->szInChI, strict );
        if (nRet == INCHI_VALID_STANDARD || nRet == INCHI_VALID_NON_STANDARD || nRet == INCHI_VALID_BETA) /* djb-rwth: removing redundant code */
        {
            ;
        }
        else
        {
            nRet = _IS_ERROR;
            goto exit_function;
        }
    }

    PrintInputParms( log_file, ip );

    /*********************************/
    /* InChI -> Structure conversion */
    /*********************************/

    /* input_file simulation */

    /*
    that was incorrect simulation, and correct one is much simpler, see below
    input_file->s.pStr = inpInChI->szInChI;
    input_file->s.nUsedLength = (int) strlen(input_file->s.pStr)+1;
    input_file->s.nAllocatedLength = input_file->s.nUsedLength;
    input_file->s.nPtr = 0;
    */
    inchi_ios_print_nodisplay( input_file, inpInChI->szInChI );

    /* buffer for the message */
    /* outStruct->szMessage = (char *)inchi_calloc( MAX_MSG_LEN, sizeof(outStruct->szMessage[0])); */

    outStruct->szMessage = (char *) inchi_calloc( MAX_MSG_LEN, sizeof( char ) );
    if (!outStruct->szMessage)
    {
        inchi_ios_eprint( log_file, "Cannot allocate output message buffer.\n" );
        nRet = -1;
    }
    else
    {
        int num_bonds;
        nRet = ReadWriteInChI( &ic, &CG , input_file, out_file, log_file,
                                ip, sd, &at, &num_at, &num_bonds,
                                &polymer, &v3000,
                                outStruct->szMessage,
                                MAX_MSG_LEN, outStruct->WarningFlags );

        if (nRet >= 0 && polymer && at && (num_at > 0)) /* djb-rwth: fixing oss-fuzz issue #68329, #68286 */
        {
            OAD_Polymer_SmartReopenCyclizedUnits( polymer, at,
                                                 num_at, &num_bonds );
        }
    }

    if (nRet >= 0 && at && num_at)
    {
        /* success */
        nRet = InpAtom0DToInchiAtom( at, num_at,
                                    &outStruct->num_atoms,
                                    &outStruct->atom,
                                    &outStruct->num_stereo0D,
                                    &outStruct->stereo0D );

        if (at)
        {
            inchi_free( at );
            at = NULL;
        }

        if (nRet >= 0 && polymer)
        {
            /* Check for and then replace ZZ for star atoms if Polymer extension is supplied */
            for (i = 0; i < outStruct->num_atoms; i++)
            {
                if (!strcmp( outStruct->atom[i].elname, "Zz" ))
                {
                    strcpy( outStruct->atom[i].elname, "*" );
                }
            }
        }

        if (nRet >= 0)
        {
            if (polymer || v3000)
            {
                nRet = SetInChIExtInputByExtOrigAtData( polymer, v3000,
                                                        &outStruct->polymer,
                                                        &outStruct->v3000,
                                                        outStruct->num_atoms ); /* pair to SetExtOrigAtDataByInChIExtInput */
                FreeExtOrigAtData( polymer, v3000 );
                polymer = NULL;
                v3000 = NULL;
            }
        }
        if (nRet < 0)
        {
            inchi_ios_eprint( log_file, "Final structure conversion failed\n" );
        }
    }
    outStruct->szLog = log_file->s.pStr;

exit_function:;

    for (i = 0; i < MAX_NUM_PATHS; i++)
    {
        if (ip->path[i])
        {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( &CG );

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if (num_repeat-- > 0)
    {
        goto repeat;
    }
#endif

#ifdef TARGET_API_LIB
    /* output */
    if (log_file->s.pStr && log_file->s.nUsedLength > 0)
    {
        while (log_file->s.nUsedLength &&
                '\n' == log_file->s.pStr[log_file->s.nUsedLength - 1])
        {
            log_file->s.pStr[--log_file->s.nUsedLength] = '\0'; /* remove last LF */
        }
        if (outStruct)
        {
            outStruct->szLog = log_file->s.pStr;
            log_file->s.pStr = NULL;
        }
    }
#endif

translate_RetVal:

    /* Close internal I/O streams */
    /* that was incorrect also
    inchi_ios_reset(input_file);  */    /* do not close input_file - its string buffer may point to inpInChI->szInChI */
    inchi_ios_close( input_file );
    inchi_ios_close( out_file );
    inchi_ios_close( log_file );

    switch (nRet)
    {
        case -3: nRet = inchi_Ret_ERROR; break; /* Error: no Structure has been created */
        case -2: nRet = inchi_Ret_ERROR; break; /* Error: no Structure has been created */
        case -1: nRet = inchi_Ret_FATAL; break; /* Severe error: no Structure has been created (typically; break; memory allocation failed) */
        default:
            if (outStruct) /* djb-rwth: fixing a NULL pointer dereference */
            {
                if (!outStruct->atom || !outStruct->num_atoms)
                {
                    nRet = inchi_Ret_EOF;
                }
                else
                {
                    int m, n, t = 0;
                    for (m = 0; m < 2; m++)
                    {
                        for (n = 0; n < 2; n++)
                        {
                            if (outStruct->WarningFlags[m][n])
                            {
                                t++;
                            }
                        }
                    }
                    nRet = t ? inchi_Ret_WARNING : inchi_Ret_OKAY;
                }
                break;
            }
    }

    return nRet;
}


/****************************************************************************

    GetStructFromINCHI

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
int INCHI_DECL GetStructFromINCHI( inchi_InputINCHI *inpInChI,
                                   inchi_OutputStruct *out )
{
    int ret = 0;

    inchi_OutputStructEx outex;
    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */

    ret = GetStructFromINCHIEx( inpInChI, &outex );

    out->szLog = outex.szLog;
    out->szMessage = outex.szMessage;
    out->WarningFlags[0][0] = outex.WarningFlags[0][0];
    out->WarningFlags[0][1] = outex.WarningFlags[0][1];
    out->WarningFlags[1][0] = outex.WarningFlags[1][0];
    out->WarningFlags[1][1] = outex.WarningFlags[1][1];

    if (ret == inchi_Ret_OKAY || ret == inchi_Ret_WARNING)
    {
        out->num_atoms = outex.num_atoms;
        out->atom = outex.atom;
        out->num_stereo0D = outex.num_stereo0D;
        out->stereo0D = outex.stereo0D;
    }

    return ret;
}



/****************************************************************************

    FreeStructFromINCHIEx

****************************************************************************/
EXPIMP_TEMPLATE INCHI_API
void INCHI_DECL FreeStructFromINCHIEx( inchi_OutputStructEx *out )
{
    if (!out)
        return;

    if (out->atom)
    {
        inchi_free( out->atom );
    }
    if (out->stereo0D)
    {
        inchi_free( out->stereo0D );
    }
    if (out->szLog)
    {
        inchi_free( out->szLog );
    }
    if (out->szMessage)
    {
        inchi_free( out->szMessage );
    }
    if (out->polymer || out->v3000)
    {
        FreeInChIExtInput( out->polymer, out->v3000 );
    }

    memset( out, 0, sizeof( *out ) ); /* djb-rwth: memset_s C11/Annex K variant? */
}


/****************************************************************************/
void FreeInChIExtInput( inchi_Input_Polymer *polymer, inchi_Input_V3000 *v3000 )
{
    int k;
    if (polymer)
    {
        if (polymer->n && polymer->units)
        {
            for (k = 0; k < polymer->n; k++)
            {
                if (polymer->units[k])
                {
                    if (polymer->units[k]->alist)
                    {
                        inchi_free( polymer->units[k]->alist );  polymer->units[k]->alist = NULL;
                    }
                    if (polymer->units[k]->blist)
                    {
                        inchi_free( polymer->units[k]->blist );  polymer->units[k]->blist = NULL;
                    }
                }
                inchi_free( polymer->units[k] );
            }
            inchi_free( polymer->units );
            polymer->units = NULL;
            inchi_free( polymer );
        }
    }
    if (v3000)
    {
        if (v3000->atom_index_orig)
        {
            inchi_free( v3000->atom_index_orig );
            v3000->atom_index_orig = NULL;
        }
        if (v3000->atom_index_fin)
        {
            inchi_free( v3000->atom_index_fin );
            v3000->atom_index_fin = NULL;
        }
        if (v3000->n_haptic_bonds && v3000->lists_haptic_bonds)
        {
            for (k = 0; k < v3000->n_haptic_bonds; k++)
            {
                if (v3000->lists_haptic_bonds[k])
                {
                    inchi_free( v3000->lists_haptic_bonds[k] );
                    v3000->lists_haptic_bonds[k] = NULL;
                }
            }
            inchi_free( v3000->lists_haptic_bonds );
            v3000->lists_haptic_bonds = NULL;
        }
        if (v3000->n_steabs && v3000->lists_steabs)
        {
            for (k = 0; k < v3000->n_steabs; k++)
            {
                if (v3000->lists_steabs[k])
                {
                    inchi_free( v3000->lists_steabs[k] );
                    v3000->lists_steabs[k] = NULL;
                }
            }
            inchi_free( v3000->lists_steabs );
            v3000->lists_steabs = NULL;
        }
        if (v3000->n_sterel && v3000->lists_sterel)
        {
            for (k = 0; k < v3000->n_sterel; k++)
            {
                if (v3000->lists_sterel[k])
                {
                    inchi_free( v3000->lists_sterel[k] );
                    v3000->lists_sterel[k] = NULL;
                }
            }
            inchi_free( v3000->lists_sterel );
            v3000->lists_sterel = NULL;
        }
        if (v3000->n_sterac && v3000->lists_sterac)
        {
            for (k = 0; k < v3000->n_sterac; k++)
            {
                if (v3000->lists_sterac[k])
                {
                    inchi_free( v3000->lists_sterac[k] );
                    v3000->lists_sterac[k] = NULL;
                }
            }
            inchi_free( v3000->lists_sterac );
            v3000->lists_sterac = NULL;
        }
        inchi_free( v3000 );
        /*memset( v3000, 0, sizeof( *v3000 ) );*/
    }
}


/****************************************************************************/
int SetExtOrigAtDataByInChIExtInput( OAD_Polymer **ppPolymer,
                                     OAD_V3000 **ppV3000,
                                     inchi_Input_Polymer *iep,
                                     inchi_Input_V3000 *iev,
                                     int nat )
{
    int    k, m, err = 0;
    OAD_V3000 *pv = NULL;

    /* Polymers */
    if (iep && iep->n)
    {
        /* Prepare OAD_Polymer container */
        *ppPolymer = (OAD_Polymer *) inchi_calloc( 1, sizeof( OAD_Polymer ) );
        if (!*ppPolymer)
        {
            err = 9001;
            goto exitf;
        }

        /* Convert Molfile's Sgroup's to OAD_PolymerUnit's */
        ( *ppPolymer )->units = (OAD_PolymerUnit**) inchi_calloc( iep->n, sizeof( ( *ppPolymer )->units[0] ) );
        if (!( *ppPolymer )->units)
        {
            err = 9001;
            goto exitf;
        }
        memset( ( *ppPolymer )->units, 0, sizeof( *( *ppPolymer )->units ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        ( *ppPolymer )->n = iep->n;
        /*( *ppPolymer )->valid = -1;*/
        ( *ppPolymer )->really_do_frame_shift = 0;

        for (k = 0; k < iep->n; k++)
        {
            int q = 0;
            OAD_PolymerUnit *unitk;

            inchi_Input_PolymerUnit *groupk = iep->units[k];
            ( *ppPolymer )->units[k] = (OAD_PolymerUnit*) inchi_calloc( 1, sizeof( OAD_PolymerUnit ) );
            unitk = ( *ppPolymer )->units[k];
            if (!unitk)
            {
                err = 9001;
                goto exitf;
            }

            memset( unitk, 0, sizeof( *unitk ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            unitk->id = groupk->id;
            unitk->type = groupk->type;
            unitk->subtype = groupk->subtype;
            unitk->conn = groupk->conn;
            unitk->label = groupk->label;

            for (q = 0; q < 4; q++)
            {
                unitk->xbr1[q] = groupk->xbr1[q];
                unitk->xbr2[q] = groupk->xbr2[q];
            }
            strcpy( unitk->smt, groupk->smt );
            unitk->na = groupk->na;
            unitk->alist = (int *) inchi_calloc( unitk->na, sizeof( int ) );
            if (!unitk->alist )
            {
                err = 9001;
                goto exitf;
            }
            for (m = 0; m < unitk->na; m++)
            {
                unitk->alist[m] = groupk->alist[m];
            }
            unitk->nb = groupk->nb;
            if (unitk->nb > 0)
            {
                unitk->blist = (int *) inchi_calloc( 2 * (long long)unitk->nb, sizeof( int ) ); /* djb-rwth: cast operator added */
                if (!unitk->blist )
                {
                    err = 9001;
                    goto exitf;
                }
                for (m = 0; m < 2 * groupk->nb; m++)
                {
                    unitk->blist[m] = groupk->blist[m];
                }
            }
            else
            {
                unitk->blist = NULL;
            }
        }
    }

    /* V3000 Extensions */
    if (iev)
    {
        int nn;
        *ppV3000 = (OAD_V3000 *) inchi_calloc( 1, sizeof( OAD_V3000 ) );
        pv = *ppV3000;
        if (!pv)
        {
            err = 9001;
            goto exitf;
        }
        memset( pv, 0, sizeof( *pv ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        pv->n_collections = iev->n_collections;
        pv->n_haptic_bonds = iev->n_haptic_bonds;
        pv->n_non_haptic_bonds = iev->n_non_haptic_bonds;
        pv->n_sgroups = iev->n_sgroups;
        pv->n_non_star_atoms = iev->n_non_star_atoms;
        pv->n_star_atoms = iev->n_star_atoms;
        pv->n_steabs = iev->n_steabs;
        pv->n_sterac = iev->n_sterac;
        pv->n_sterel = iev->n_sterel;
        pv->n_3d_constraints = iev->n_3d_constraints;

        if (iev->atom_index_orig)
        {
            pv->atom_index_orig = (int *) inchi_calloc( nat, sizeof( int ) );
            if (NULL == pv->atom_index_orig)
            {
                err = 9001;
                goto exitf;
            }
            memcpy( pv->atom_index_orig, iev->atom_index_orig, nat );
        }
        if (iev->atom_index_fin)
        {
            pv->atom_index_fin = (int *) inchi_calloc( nat, sizeof( int ) );
            if (NULL == pv->atom_index_fin)
            {
                err = 9001;
                goto exitf;
            }
            memcpy( pv->atom_index_fin, iev->atom_index_fin, nat );
        }
        if (iev->n_haptic_bonds && iev->lists_haptic_bonds)
        {
            pv->lists_haptic_bonds = (int **) inchi_calloc( iev->n_haptic_bonds, sizeof( int* ) );
            if (NULL == pv->lists_haptic_bonds)
            {
                err = 9001;
                goto exitf;
            }
            for (m = 0; m < iev->n_haptic_bonds; m++)
            {
                int *lst = NULL;
                int *mol_lst = iev->lists_haptic_bonds[m];
                nn = mol_lst[2] + 3;
                lst = pv->lists_haptic_bonds[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (iev->n_steabs && iev->lists_steabs)
        {
            pv->lists_steabs = (int **) inchi_calloc( iev->n_steabs, sizeof( int* ) );
            if (NULL == pv->lists_steabs) { err = 9001; goto exitf; }
            for (m = 0; m < iev->n_steabs; m++)
            {
                int *lst = NULL;
                int *mol_lst = iev->lists_steabs[m];
                nn = mol_lst[1] + 2;
                lst = pv->lists_steabs[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (iev->n_sterac && iev->lists_sterac)
        {
            pv->lists_sterac = (int **) inchi_calloc( iev->n_sterac, sizeof( int* ) );
            if (NULL == pv->lists_sterac) { err = 9001; goto exitf; }
            for (m = 0; m < iev->n_sterac; m++)
            {
                int *lst = NULL;
                int *mol_lst = iev->lists_sterac[m];
                nn = mol_lst[1] + 2;
                lst = pv->lists_sterac[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (iev->n_sterel && iev->lists_sterel)
        {
            pv->lists_sterel = (int **) inchi_calloc( iev->n_sterel, sizeof( int* ) );
            if (NULL == pv->lists_sterel) { err = 9001; goto exitf; }
            for (m = 0; m < iev->n_sterel; m++)
            {
                int *lst = NULL;
                int *mol_lst = iev->lists_sterel[m];
                nn = mol_lst[1] + 2;
                lst = pv->lists_sterel[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
    }

exitf:
    if (err)
    {
        FreeExtOrigAtData( *ppPolymer, pv );
    }

    return err;
}


/****************************************************************************/
int SetInChIExtInputByExtOrigAtData( OAD_Polymer     *orp,
                                     OAD_V3000     *orv,
                                     inchi_Input_Polymer **iip,
                                     inchi_Input_V3000     **iiv,
                                     int nat )
{
    int    k, m, err = 0;

        /* Polymers */
    if (orp && orp->n > 0)
    {
        /* djb-rwth: fixing oss-fuzz issue #67695, #66748 */
        inchi_Input_Polymer* iip_tmp = (inchi_Input_Polymer*) inchi_calloc( 1, sizeof( inchi_Input_Polymer ) );
        inchi_Input_PolymerUnit** units_tmp = (inchi_Input_PolymerUnit**) inchi_calloc( orp->n, sizeof( ( *iip )->units[0] ) );
        int** uk_al_tmp = (int**)inchi_malloc((orp->n) * sizeof(int*));
        inchi_Input_PolymerUnit** unitk = (inchi_Input_PolymerUnit**)inchi_malloc((orp->n) * sizeof(inchi_Input_PolymerUnit*));

        if (!iip_tmp || !units_tmp || !uk_al_tmp || !unitk)
        {
            err = 9001;
            goto preexitf;
        }
        /* *iip = iip_tmp; */
        iip_tmp->n = orp->n;
        iip_tmp->units = units_tmp;
        memset(units_tmp, 0, sizeof( *units_tmp) ); /* djb-rwth: memset_s C11/Annex K variant? */
        for (k = 0; k < orp->n; k++)
        {
            int q = 0;
            unitk[k] = (inchi_Input_PolymerUnit*)inchi_calloc(1, sizeof(inchi_Input_PolymerUnit));
            OAD_PolymerUnit    *groupk = orp->units[k];
            /* unitk = ( *iip )->units[k]; */
            if (!unitk[k])
            {
                err = 9001; 
                goto preexitf;
            }
            iip_tmp->units[k] = unitk[k];
            memset( unitk[k], 0, sizeof(*unitk[k])); /* djb-rwth: memset_s C11/Annex K variant? */
            unitk[k]->id = groupk->id;
            unitk[k]->type = groupk->type;
            unitk[k]->subtype = groupk->subtype;
            unitk[k]->conn = groupk->conn;
            unitk[k]->label = groupk->label;
            for (q = 0; q < 4; q++)
            {
                unitk[k]->xbr1[q] = groupk->xbr1[q];
                unitk[k]->xbr2[q] = groupk->xbr2[q];
            }
            strcpy( unitk[k]->smt, groupk->smt);
            unitk[k]->na = groupk->na;
            uk_al_tmp[k] = (int*)inchi_calloc(unitk[k]->na, sizeof(int));
            if (!uk_al_tmp[k])
            {
                err = 9001; 
                goto preexitf;
            }
            unitk[k]->alist = uk_al_tmp[k];
            for (m = 0; m < unitk[k]->na; m++)
            {
                uk_al_tmp[k][m] = groupk->alist[m];
            }
            unitk[k]->nb = groupk->nb;
            if (unitk[k]->nb > 0)
            {
                unitk[k]->blist = (int*)inchi_calloc(2 * (long long)unitk[k]->nb, sizeof(int)); /* djb-rwth: cast operator added */
                if (!unitk[k]->blist)
                {
                    err = 9001;
                    goto preexitf;
                }
                for (m = 0; m < 2 * groupk->nb; m++)
                {
                    unitk[k]->blist[m] = groupk->blist[m];
                }
            }
            else
            {
                unitk[k]->blist = NULL;
            }

            inchi_free(unitk[k]);
            inchi_free(uk_al_tmp[k]);
        }
    /* djb-rwth: avoiding memory leak */
    preexitf:
        /* djb-rwth: fixing GHI #165 */
        if (iip_tmp && *iip) /* djb-rwth: fixing oss-fuzz issue #455987437 */
        {
            memcpy(*iip, iip_tmp, sizeof(inchi_Input_Polymer));
        }
        inchi_free(iip_tmp);
        inchi_free(units_tmp);
        inchi_free(uk_al_tmp);
        inchi_free(unitk);
        if (err)
        {
            goto exitf;
        }
    }

    if (orv)
    {
        int nn;
        *iiv = (inchi_Input_V3000 *) inchi_calloc( 1, sizeof(inchi_Input_V3000) ); /* djb-rwth: fixing the incorrect type of variable */
        if (!*iiv)
        {
            err = 9001;
            goto exitf;
        }
        memset( *iiv, 0, sizeof( **iiv ) ); /* djb-rwth: memset_s C11/Annex K variant? */

        ( *iiv )->n_collections = orv->n_collections;
        ( *iiv )->n_haptic_bonds = orv->n_haptic_bonds;
        ( *iiv )->n_non_haptic_bonds = orv->n_non_haptic_bonds;
        ( *iiv )->n_sgroups = orv->n_sgroups;
        ( *iiv )->n_non_star_atoms = orv->n_non_star_atoms;
        ( *iiv )->n_star_atoms = orv->n_star_atoms;
        ( *iiv )->n_steabs = orv->n_steabs;
        ( *iiv )->n_sterac = orv->n_sterac;
        ( *iiv )->n_sterel = orv->n_sterel;
        ( *iiv )->n_3d_constraints = orv->n_3d_constraints;

        if (orv->atom_index_orig)
        {
            ( *iiv )->atom_index_orig = (int *) inchi_calloc( nat, sizeof( int ) );
            if (NULL == ( *iiv )->atom_index_orig)
            {
                err = 9001;
                goto exitf;
            }
            memcpy( ( *iiv )->atom_index_orig, orv->atom_index_orig, nat );
        }
        if (orv->atom_index_fin)
        {
            ( *iiv )->atom_index_fin = (int *) inchi_calloc( nat, sizeof( int ) );
            if (NULL == ( *iiv )->atom_index_fin)
            {
                err = 9001;
                goto exitf;
            }
            memcpy( ( *iiv )->atom_index_fin, orv->atom_index_fin, nat );
        }
        if (orv->n_haptic_bonds && orv->lists_haptic_bonds)
        {
            ( *iiv )->lists_haptic_bonds = (int **) inchi_calloc( orv->n_haptic_bonds, sizeof( int* ) );
            if (NULL == ( *iiv )->lists_haptic_bonds)
            {
                err = 9001;
                goto exitf;
            }
            for (m = 0; m < orv->n_haptic_bonds; m++)
            {
                int *lst = NULL;
                int *mol_lst = orv->lists_haptic_bonds[m];
                nn = mol_lst[2] + 3;
                lst = ( *iiv )->lists_haptic_bonds[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (orv->n_steabs && orv->lists_steabs)
        {
            ( *iiv )->lists_steabs = (int **) inchi_calloc( orv->n_steabs, sizeof( int* ) );
            if (NULL == ( *iiv )->lists_steabs) { err = 9001; goto exitf; }
            for (m = 0; m < orv->n_steabs; m++)
            {
                int *lst = NULL;
                int *mol_lst = orv->lists_steabs[m];
                nn = mol_lst[1] + 2;
                lst = ( *iiv )->lists_steabs[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (orv->n_sterac && orv->lists_sterac)
        {
            ( *iiv )->lists_sterac = (int **) inchi_calloc( orv->n_sterac, sizeof( int* ) );
            if (NULL == ( *iiv )->lists_sterac) { err = 9001; goto exitf; }
            for (m = 0; m < orv->n_sterac; m++)
            {
                int *lst = NULL;
                int *mol_lst = orv->lists_sterac[m];
                nn = mol_lst[1] + 2;
                lst = ( *iiv )->lists_sterac[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
        if (orv->n_sterel && orv->lists_sterel)
        {
            ( *iiv )->lists_sterel = (int **) inchi_calloc( orv->n_sterel, sizeof( int* ) );
            if (NULL == ( *iiv )->lists_sterel) { err = 9001; goto exitf; }
            for (m = 0; m < orv->n_sterel; m++)
            {
                int *lst = NULL;
                int *mol_lst = orv->lists_sterel[m];
                nn = mol_lst[1] + 2;
                lst = ( *iiv )->lists_sterel[m] = (int *) inchi_calloc( nn, sizeof( int ) );
                if (NULL == lst)
                {
                    err = 9001;
                    goto exitf;
                }
                for (k = 0; k < nn; k++)
                {
                    lst[k] = mol_lst[k];
                }
            }
        }
    }

exitf:
    if (err)
    {
        FreeInChIExtInput( *iip, *iiv );
    }

    return err;
}


#if( defined( _WIN32 ) && defined( _MSC_VER ) && _MSC_VER >= 800 && defined(_USRDLL) && defined(BUILD_LINK_AS_DLL) )
/* Win32 & MS VC ++, compile and link as a DLL */

/*********************************************************/
/*   C calling conventions export from Win32 dll         */
/*********************************************************/

/* prototypes */
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

    int cdecl_GetINCHI( inchi_Input *inp, inchi_Output *out );
    int cdecl_GetStdINCHI( inchi_Input *inp, inchi_Output *out );
    void cdecl_FreeINCHI( inchi_Output *out );
    void cdecl_FreeStdINCHI( inchi_Output *out );
    int  cdecl_GetStringLength( char *p );
    int  cdecl_Get_inchi_Input_FromAuxInfo( char *szInchiAuxInfo,
                                            int bDoNotAddH,
                                            int bDiffUnkUndfStereo,
                                            InchiInpData *pInchiInp );
    int  cdecl_Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo,
                                                int bDoNotAddH,
                                                InchiInpData *pInchiInp );
    /*void cdecl_Free_inchi_Input( inchi_Input *pInp );*/
    void cdecl_Free_std_inchi_Input( inchi_Input *pInp );
    int cdecl_GetStructFromINCHI( inchi_InputINCHI *inpInChI,
                                  inchi_OutputStruct *outStruct );
    int cdecl_GetStructFromStdINCHI( inchi_InputINCHI *inpInChI,
                                     inchi_OutputStruct *outStruct );
    int cdecl_GetINCHIfromINCHI( inchi_InputINCHI *inpInChI,
                                 inchi_Output *out );
    void cdecl_FreeStructFromINCHI( inchi_OutputStruct *outStruct );
    void cdecl_FreeStructFromStdINCHI( inchi_OutputStruct *outStruct );
    int cdecl_CheckINCHI( const char *szINCHI, const int strict );

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

/* implementation */
/* libinchi.def provides export without cdecl_ prefixes */


/****************************************************************************/
int cdecl_GetINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetINCHI( inp, out );
}


/****************************************************************************/
int cdecl_GetStdINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetStdINCHI( inp, out );
}


/****************************************************************************/
void cdecl_FreeINCHI( inchi_Output *out )
{
    FreeINCHI( out );
}


/****************************************************************************/
void cdecl_FreeStdINCHI( inchi_Output *out )
{
    FreeStdINCHI( out );
}


/****************************************************************************/
int cdecl_GetStringLength( char *p )
{
    return GetStringLength( p );
}


/****************************************************************************/
int cdecl_Get_inchi_Input_FromAuxInfo( char *szInchiAuxInfo,
                                      int bDoNotAddH,
                                      int bDiffUnkUndfStereo,
                                      InchiInpData *pInchiInp )
{
    return Get_inchi_Input_FromAuxInfo( szInchiAuxInfo,
                                        bDoNotAddH,
                                        bDiffUnkUndfStereo,
                                        pInchiInp );
}


/****************************************************************************/
int cdecl_Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo,
                                          int bDoNotAddH,
                                          InchiInpData *pInchiInp )
{
    return Get_std_inchi_Input_FromAuxInfo( szInchiAuxInfo,
                                            bDoNotAddH,
                                            pInchiInp );
}


/****************************************************************************/
void cdecl_Free_std_inchi_Input( inchi_Input *pInp )
{
    Free_std_inchi_Input( pInp );
}


/****************************************************************************/
void cdecl_Free_inchi_Input( inchi_Input *pInp )
{
    Free_inchi_Input( pInp );
}


/****************************************************************************/
int cdecl_GetStructFromINCHI( inchi_InputINCHI *inpInChI,
                              inchi_OutputStruct *outStruct )
{
    return GetStructFromINCHI( inpInChI, outStruct );
}


/****************************************************************************/
int cdecl_GetStructFromStdINCHI( inchi_InputINCHI *inpInChI,
                                 inchi_OutputStruct *outStruct )
{
    return GetStructFromStdINCHI( inpInChI, outStruct );
}

/********************************************************/
void cdecl_FreeStructFromINCHI( inchi_OutputStruct *outStruct )
{
    FreeStructFromINCHI( outStruct );
}


/****************************************************************************/
int cdecl_GetINCHIfromINCHI( inchi_InputINCHI *inpInChI,
                             inchi_Output *out )
{
    return GetINCHIfromINCHI( inpInChI, out );
}


/****************************************************************************/
void cdecl_FreeStructFromStdINCHI( inchi_OutputStruct *outStruct )
{
    FreeStructFromStdINCHI( outStruct );
}


/****************************************************************************/
int cdecl_CheckINCHI( const char *szINCHI, const int strict )
{
    return CheckINCHI( szINCHI, strict );
}
#endif

#if( defined(__GNUC__) && __GNUC__ >= 3 && defined(__MINGW32__) && defined(_WIN32) )
#include <windows.h>
/*********************************************************/
/*   Pacal calling conventions export from Win32 dll     */
/*********************************************************/
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif
/* prototypes */

    int  PASCAL pasc_GetINCHI( inchi_Input *inp, inchi_Output *out );
    int  PASCAL pasc_GetStdINCHI( inchi_Input *inp, inchi_Output *out );
    void PASCAL pasc_FreeINCHI( inchi_Output *out );
    void PASCAL pasc_FreeStdINCHI( inchi_Output *out );
    int  PASCAL pasc_GetStringLength( char *p );
    int  PASCAL pasc_Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo,
                                                      int bDoNotAddH,
                                                      InchiInpData *pInchiInp );
    int  PASCAL pasc_Get_inchi_Input_FromAuxInfo( char *szInchiAuxInfo,
                                                      int bDoNotAddH,
                                                      int bDiffUnkUndfStereo,
                                                      InchiInpData *pInchiInp );
    void PASCAL pasc_Free_inchi_Input( inchi_Input *pInp );
    void PASCAL pasc_Free_std_inchi_Input( inchi_Input *pInp );
    void PASCAL pasc_FreeStructFromINCHI( inchi_OutputStruct *out );
    void PASCAL pasc_FreeStructFromStdINCHI( inchi_OutputStruct *out );
    int PASCAL pasc_GetStructFromINCHI( inchi_InputINCHI *inp, inchi_OutputStruct *out );
    int PASCAL pasc_GetStructFromStdINCHI( inchi_InputINCHI *inp, inchi_OutputStruct *out );
    int PASCAL pasc_CheckINCHI( const char *szINCHI, const int strict );

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

/* implementation */
/* libinchi.def provides export without PASCAL pasc_ prefixes */


/****************************************************************************/
int PASCAL pasc_GetINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetINCHI( inp, out );
}


/****************************************************************************/
int PASCAL pasc_GetStdINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetStdINCHI( inp, out );
}


/****************************************************************************/
void PASCAL pasc_FreeINCHI( inchi_Output *out )
{
    FreeINCHI( out );
}


/****************************************************************************/
void PASCAL pasc_FreeStdINCHI( inchi_Output *out )
{
    FreeStdINCHI( out );
}


/****************************************************************************/
int PASCAL pasc_GetStringLength( char *p )
{
    return GetStringLength( p );
}


/****************************************************************************/
int PASCAL pasc_Get_inchi_Input_FromAuxInfo( char *szInchiAuxInfo,
                                                int bDoNotAddH,
                                                int bDiffUnkUndfStereo,
                                                InchiInpData *pInchiInp )
{
    return Get_inchi_Input_FromAuxInfo( szInchiAuxInfo, bDoNotAddH,
                                            bDiffUnkUndfStereo, pInchiInp );
}


/****************************************************************************/
int PASCAL pasc_Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo,
                                                int bDoNotAddH,
                                                InchiInpData *pInchiInp )
{
    return Get_std_inchi_Input_FromAuxInfo( szInchiAuxInfo, bDoNotAddH, pInchiInp );
}


/****************************************************************************/
void PASCAL pasc_Free_inchi_Input( inchi_Input *pInp )
{
    Free_inchi_Input( pInp );
}


/****************************************************************************/
void PASCAL pasc_Free_std_inchi_Input( inchi_Input *pInp )
{
    Free_std_inchi_Input( pInp );
}


/****************************************************************************/
void PASCAL pasc_FreeStructFromINCHI( inchi_OutputStruct *out )
{
    FreeStructFromINCHI( out );
}


/****************************************************************************/
void PASCAL pasc_FreeStructFromStdINCHI( inchi_OutputStruct *out )
{
    FreeStructFromStdINCHI( out );
}


/****************************************************************************/
int PASCAL pasc_GetStructFromINCHI( inchi_InputINCHI *inp, inchi_OutputStruct *out )
{
    return GetStructFromINCHI( inp, out );
}


/****************************************************************************/
int PASCAL pasc_GetStructFromStdINCHI( inchi_InputINCHI *inp, inchi_OutputStruct *out )
{
    return GetStructFromStdINCHI( inp, out );
}


/****************************************************************************/
int PASCAL pasc_CheckINCHI( const char *szINCHI, const int strict )
{
    return CheckINCHI( szINCHI, strict );
}

#endif
