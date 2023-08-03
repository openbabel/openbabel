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

#include "mode.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <float.h>
#include <math.h>
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








/*************************************************************************
 *
 *   Local prototypes
 *
 *************************************************************************/


int SetAtomProperties( inp_ATOM *at, MOL_COORD *szCoord, inchi_Atom *ati,
                       int a1, int *nDim, char *pStrErr, int *err );
int SetBondProperties( inp_ATOM *at, inchi_Atom *ati, int a1, int j,
                       int nNumAtoms, int *nNumBonds, char *pStrErr, int *err );
int SetAtomAndBondProperties( inp_ATOM *at, inchi_Atom *ati, int a1,
                              int bDoNotAddH, char *pStrErr, int *err );
void SetNumImplicitH(inp_ATOM* at, int num_atoms);
int Extract0DParities(inp_ATOM *at, int nNumAtoms, inchi_Stereo0D *stereo0D,
                       int num_stereo0D, char *pStrErr, int *err, int vABParityUnknown);
int parse_options_string ( char *cmd, const char *argv[], int maxargs );

int InpAtom0DToInchiAtom( inp_ATOM *at, int num_atoms, inchi_OutputStruct *outStruct );


int ExtractOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
         inchi_Input *inp, 
         INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file, INCHI_IOSTREAM *prb_file,
         ORIG_ATOM_DATA *orig_inp_data, long *num_inp, char *pStr, int nStrLen );


static int GetINCHI1(inchi_Input *inp, inchi_Output *out, int bStdFormat);

/*************************************************************************/

int bInterrupted = 0;

/********************************************************************
 *
 * INCHI API: DEALLOCATE INCHI OUTPUT
 *
 ********************************************************************/


EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeINCHI( inchi_Output *out )
{
    if ( out->szInChI ) {
        inchi_free( out->szInChI );
    }
    if ( out->szLog ) {
        inchi_free( out->szLog );
    }
    if ( out->szMessage ) {
        inchi_free( out->szMessage );
    }
    memset( out, 0, sizeof(*out) );
}


EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeStdINCHI( inchi_Output *out )
{
    FreeINCHI( out );
}

EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeStructFromStdINCHI( inchi_OutputStruct *out )
{
    FreeStructFromINCHI( out );
}


/*******************************************************************/
EXPIMP_TEMPLATE INCHI_API void INCHI_DECL FreeStructFromINCHI( inchi_OutputStruct *out )
{
    if ( out->atom ) {
        inchi_free( out->atom );
    }
    if ( out->stereo0D ) {
        inchi_free( out->stereo0D );
    }
    if ( out->szLog ) {
        inchi_free( out->szLog );
    }
    if ( out->szMessage ) {
        inchi_free( out->szMessage );
    }
    memset( out, 0, sizeof(*out) );
}
/********************************************************************/
#define INCHI_MAX_NUM_ARG 32
/********************************************************************
 *
 *    INCHI API: MAIN ENTRY POINT
 *
 ********************************************************************/

int bLibInchiSemaphore = 0;


EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStdINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetINCHI1( inp, out, 1 );
}




EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetINCHI1( inp, out, 0 );
}


static int GetINCHI1(inchi_Input *inp, inchi_Output *out, int bStdFormat)
{
    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;       
    char szTitle[MAX_SDF_HEADER+MAX_SDF_VALUE+256];

    int i;
    long num_inp, num_err;
    char      szSdfDataValue[MAX_SDF_VALUE+1];
    PINChI2     *pINChI[INCHI_NUM];
    PINChI_Aux2 *pINChI_Aux[INCHI_NUM];

    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */
    unsigned long  ulTotalProcessingTime = 0;

    INPUT_PARMS inp_parms;
    INPUT_PARMS *ip = &inp_parms;

    ORIG_ATOM_DATA OrigAtData; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *orig_inp_data = &OrigAtData;
    ORIG_ATOM_DATA PrepAtData[2]; /* 0=> disconnected, 1=> original */
    ORIG_ATOM_DATA *prep_inp_data = PrepAtData;
    int             bReleaseVersion = bRELEASE_VERSION;
    const int nStrLen = 64000;
    char *pStr = NULL;
    int   nRet = 0, nRet1;

    STRUCT_FPTRS *pStructPtrs = NULL;

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif

    const char *argv[INCHI_MAX_NUM_ARG+1];
    int   argc;
    char *szOptions = NULL;

    INCHI_IOSTREAM inchi_file[3], *output_file = inchi_file, *log_file = inchi_file+1;
    INCHI_IOSTREAM prb_file0, *prb_file = &prb_file0;


    if ( bLibInchiSemaphore ) {  /* does not work properly under sufficient stress */
        return inchi_Ret_BUSY;
    }
    bLibInchiSemaphore = 1;

#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);
/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE 
   _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode(_CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG);
#endif

#if ( !defined(__STDC__) || __STDC__ != 1 )
    /* turn on floating point exceptions */
    {
        /* Get the default control word. */
        int cw = _controlfp( 0,0 );

        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_ZERODIVIDE|EM_DENORMAL);

        /* Set the control word. */
        _controlfp( cw, MCW_EM );
 
    }
#endif
#endif




    szTitle[0] = '\0';

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
repeat:
    inchi_ios_close(output_file);
    inchi_ios_close(log_file);
    inchi_ios_close(prb_file);
    pStr = NULL;
#endif

    /*^^^ Initialize internal for this function output streams as string buffers */
    inchi_ios_init(output_file, INCHI_IOSTREAM_STRING, NULL);
    inchi_ios_init(log_file, INCHI_IOSTREAM_STRING, NULL);
    inchi_ios_init(prb_file, INCHI_IOSTREAM_STRING, NULL);


    num_inp    = 0;
    num_err    = 0;
    sd->bUserQuit  = 0;

    /* clear original input structure */
    memset( pINChI,     0, sizeof(pINChI    ) );
    memset( pINChI_Aux, 0, sizeof(pINChI_Aux) );
    memset( sd,         0, sizeof(*sd) );
    memset( ip,         0, sizeof(*ip) );
    memset( orig_inp_data     , 0,   sizeof( *orig_inp_data  ) );
    memset( prep_inp_data     , 0, 2*sizeof( *prep_inp_data  ) );
    memset( szSdfDataValue    , 0, sizeof( szSdfDataValue    ) );

    if ( !out ) {
        nRet = _IS_ERROR;
        goto exit_function;
    }
    memset( out, 0, sizeof(*out) );

    /* options */
    if ( inp && inp->szOptions ) {
        szOptions = (char*)inchi_malloc( strlen(inp->szOptions) + 1 );
        if ( szOptions ) {
            strcpy( szOptions, inp->szOptions );
            argc = parse_options_string ( szOptions, argv, INCHI_MAX_NUM_ARG );
        } else {
            nRet = _IS_FATAL;
            goto translate_RetVal; /* emergency exit */
        }
    } else {
        argc = 1;
        argv[0] = "";
        argv[1] = NULL;
    }

    if ( (argc == 1
#ifdef TARGET_API_LIB
        && (!inp || inp->num_atoms <= 0 || !inp->atom)
#endif   
        )
        || (argc==2 && ( argv[1][0]==INCHI_OPTION_PREFX ) &&
        (!strcmp(argv[1]+1, "?") || !stricmp(argv[1]+1, "help") )) ) {
        HelpCommandLineParms(log_file);
        out->szLog = log_file->s.pStr;
        memset( log_file, 0, sizeof(*log_file) );
        nRet = _IS_EOF;
        goto translate_RetVal;
    }

    nRet1 = ReadCommandLineParms( argc, argv, ip, szSdfDataValue, &ulDisplTime, bReleaseVersion, log_file );
    if ( szOptions ) {
        inchi_free( szOptions );
        szOptions = NULL;
    }
    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if ( 0 > nRet1 ) {
        nRet = _IS_FATAL;
        goto exit_function;
    }
    if ( ip->bNoStructLabels ) {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    } else
    if ( ip->nInputType == INPUT_INCHI_XML || ip->nInputType == INPUT_INCHI_PLAIN  || ip->nInputType == INPUT_CMLFILE ) {
        /* the input may contain both the header and the label of the structure */
        if ( !ip->pSdfLabel ) 
            ip->pSdfLabel  = ip->szSdfDataHeader;
        if ( !ip->pSdfValue )
            ip->pSdfValue  = szSdfDataValue;
    }

    /* Ensure standardness */
    if ( bStdFormat )
    {
        if ( ip->bINChIOutputOptions & INCHI_OUT_SAVEOPT )
        {
            ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;
        }
        if ( 0 != ( ip->bTautFlags & TG_FLAG_RECONNECT_COORD) )
        {
            ip->bTautFlags &= ~TG_FLAG_RECONNECT_COORD;
        }
        if ( 0 != (ip->nMode & REQ_MODE_BASIC) )
        {
            ip->nMode &= ~REQ_MODE_BASIC;
        }
        if ( 0 != ( ip->nMode & REQ_MODE_RELATIVE_STEREO) ) 
        {
            ip->nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
        }
        if ( 0 != ( ip->nMode & REQ_MODE_RACEMIC_STEREO) ) 
        {
            ip->nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
        }
        if ( 0 != ( ip->nMode & REQ_MODE_CHIR_FLG_STEREO) ) 
        {
            ip->nMode &= ~(REQ_MODE_RACEMIC_STEREO | REQ_MODE_RELATIVE_STEREO | REQ_MODE_CHIR_FLG_STEREO);
        }
        if ( 0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO) ) 
        {
            ip->nMode &= ~REQ_MODE_DIFF_UU_STEREO;
        }
        if ( 0 == (ip->nMode & (REQ_MODE_SB_IGN_ALL_UU | REQ_MODE_SC_IGN_ALL_UU))  ) 
        {
            ip->nMode |= REQ_MODE_SB_IGN_ALL_UU;
            ip->nMode |= REQ_MODE_SC_IGN_ALL_UU;
        }	
        if ( 0 != (ip->bTautFlags & TG_FLAG_KETO_ENOL_TAUT) )
        {
            ip->bTautFlags  &= ~TG_FLAG_KETO_ENOL_TAUT;
        }
        if ( 0 != (ip->bTautFlags & TG_FLAG_1_5_TAUT) )
        {
            ip->bTautFlags  &= ~TG_FLAG_1_5_TAUT;
        }
        /* And anyway... */
        ip->bINChIOutputOptions |= INCHI_OUT_STDINCHI;
        ip->bINChIOutputOptions &= ~INCHI_OUT_SAVEOPT;
    }
    /* */


    PrintInputParms( log_file, ip );
    if ( !(pStr = (char*)inchi_malloc(nStrLen))) {
        inchi_ios_eprint( log_file, "Cannot allocate output buffer. Terminating\n");
        goto exit_function;
    }
    pStr[0] = '\0';



    /**********************************************************************************************/
    /*  Main cycle */
    /*  read input structures and create their INChI */
    ulTotalProcessingTime = 0;

    if ( pStructPtrs ) {
        memset ( pStructPtrs, 0, sizeof(pStructPtrs[0]) );
    }

    /* === possible improvement: convert inp to orig_inp_data ==== */
    if ( !sd->bUserQuit && !bInterrupted )
    {
        if ( ip->last_struct_number && num_inp >= ip->last_struct_number ) {
            nRet = _IS_EOF; /*  simulate end of file */
            goto exit_function;
        }

        nRet = ExtractOneStructure( sd, ip, szTitle, inp, log_file, output_file, prb_file,
                                orig_inp_data, &num_inp, pStr, nStrLen );

        if ( pStructPtrs ) {
            pStructPtrs->cur_fptr ++;
        }

#ifndef TARGET_API_LIB
        if ( sd->bUserQuit ) {
            break;
        }
#endif
        switch ( nRet ) {
        case _IS_FATAL:
            num_err ++;
            goto exit_function;
        case _IS_EOF:
            goto exit_function;
        case _IS_ERROR:
            num_err ++;
            goto exit_function;
#ifndef TARGET_API_LIB
        case _IS_SKIP:
            continue;
#endif
        }

        /* create INChI for each connected component of the structure and optionally display them */
        /* output INChI for the whole structure */
        nRet1 = ProcessOneStructure( sd, ip, szTitle, pINChI, pINChI_Aux,
                                     NULL, /* inp_file is not necessary as all input is already saved in 'ip' */
                                     log_file, output_file, prb_file,
                                     orig_inp_data, prep_inp_data,
                                     num_inp, pStr, nStrLen,
                                     0 /* save_opt_bits */);

        /*  free INChI memory */
        FreeAllINChIArrays( pINChI, pINChI_Aux, sd->num_components );
        
        /* free structure data */
        FreeOrigAtData( orig_inp_data );
        FreeOrigAtData( prep_inp_data );
        FreeOrigAtData( prep_inp_data+1 );

        ulTotalProcessingTime += sd->ulStructTime;
        nRet = inchi_max(nRet, nRet1);
        switch ( nRet ) {
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
    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && sd->bXmlStructStarted > 0 ) {
        if ( !OutputINChIXmlStructEndTag( output_file, pStr, nStrLen, 1 ) ) {
            inchi_ios_eprint( log_file, "Cannot create end xml tag for structure #%d.%s%s%s%s Terminating.\n", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
            sd->bXmlStructStarted = -1; /*  do not repeat same message */
        }
    }


    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && ip->bXmlStarted ) {
        OutputINChIXmlRootEndTag( output_file );
        ip->bXmlStarted = 0;
    }

    
    /* avoid memory leaks in case of fatal error */
    if ( pStructPtrs && pStructPtrs->fptr ) {
        inchi_free( pStructPtrs->fptr );
    }

    /*  free INChI memory */
    FreeAllINChIArrays( pINChI, pINChI_Aux, sd->num_components );
    /* free structure data */
    FreeOrigAtData( orig_inp_data );
    FreeOrigAtData( prep_inp_data );
    FreeOrigAtData( prep_inp_data+1 );

#if( ADD_CMLPP == 1 )
        /* BILLY 8/6/04 */
        /* free CML memory */
        FreeCml ();
        FreeCmlDoc( 1 );
#endif

    if ( pStr ) {
        inchi_free( pStr );
    }

    for ( i = 0; i < MAX_NUM_PATHS; i ++ ) {
        if ( ip->path[i] ) {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( );


#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if ( num_repeat-- > 0 ) {
        goto repeat;
    }
#endif


    /* output */
    if ( sd->pStrErrStruct[0] ) {
        if ( out && (out->szMessage = (char *)inchi_malloc( strlen(sd->pStrErrStruct) + 1 )) ) {
            strcpy( out->szMessage, sd->pStrErrStruct );
        }
    }
    if ( output_file->s.pStr && output_file->s.nUsedLength > 0 && out ) {
        char *p;
        out->szInChI   = output_file->s.pStr;
        out->szAuxInfo = NULL;
        if ( !(INCHI_OUT_SDFILE_ONLY & ip->bINChIOutputOptions ) ) /* do not remove last LF from SDF output - 2008-12-23 DT */
        for ( p = strchr(out->szInChI, '\n'); p; p = strchr(p+1, '\n') ) {
            if ( !memcmp( p, "\nAuxInfo", 8 ) ) {
                *p = '\0';            /* remove LF after INChI */
                out->szAuxInfo = p+1; /* save pointer to AuxInfo */
            } else
            if ( out->szAuxInfo || !p[1]) {   /* remove LF after aux info or from the last char */
                *p = '\0';
                break;
            }
        }
       output_file->s.pStr = NULL;
    }


    if ( log_file->s.pStr && log_file->s.nUsedLength > 0 ) {
        while ( log_file->s.nUsedLength && '\n' == log_file->s.pStr[log_file->s.nUsedLength-1] ) {
            log_file->s.pStr[-- log_file->s.nUsedLength]  = '\0'; /* remove last LF */
        }
        if ( out ) {
            out->szLog = log_file->s.pStr;
            log_file->s.pStr = NULL;
        }
    }

    
    
translate_RetVal:

    /* Close inernal I/O streams */
    inchi_ios_close(log_file);
    inchi_ios_close(output_file);
    inchi_ios_close(prb_file);
    
    switch (nRet) { 
    case _IS_SKIP   : nRet = inchi_Ret_SKIP   ; break; /* not used in INChI dll */
    case _IS_EOF    : nRet = inchi_Ret_EOF    ; break; /* no structural data has been provided */
    case _IS_OKAY   : nRet = inchi_Ret_OKAY   ; break; /* Success; break; no errors or warnings */
    case _IS_WARNING: nRet = inchi_Ret_WARNING; break; /* Success; break; warning(s) issued */
    case _IS_ERROR  : nRet = inchi_Ret_ERROR  ; break; /* Error: no INChI has been created */
    case _IS_FATAL  : nRet = inchi_Ret_FATAL  ; break; /* Severe error: no INChI has been created (typically; break; memory allocation failed) */
    case _IS_UNKNOWN:
    default         : nRet = inchi_Ret_UNKNOWN; break; /* Unlnown program error */
    }
    bLibInchiSemaphore = 0;

    return nRet;
}


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Check if the string represents valid InChI/standard InChI.          
Input:
        szINCHI     source InChI
        strict      if 0, just quickly check for proper layout (prefix, version, etc.)
                    The result may not be strict.
                    If not 0, try to perform InChI2InChI conversion and 
                    returns success if a resulting InChI string exactly match source.
                    The result may be 'false alarm' due to imperfect algorithm of
                    conversion.
Returns:
        success/errors codes

*/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL CheckINCHI(const char *szINCHI, const int strict)
{
int ret=INCHI_VALID_NON_STANDARD;
int ret_i2i;
inchi_InputINCHI    inchi_inp;
inchi_Output        inchi_out;
size_t slen, pos_slash1=0;
char *str = NULL;
size_t i;
size_t slen0;
char pp;

    /* .. non-empty */
    if (szINCHI==NULL)
        return INCHI_INVALID_PREFIX;
    
    slen = strlen(szINCHI);
    
    /* .. has valid prefix */
    if (slen<LEN_INCHI_STRING_PREFIX+3)                         
        return INCHI_INVALID_PREFIX;    
    if (memcmp(szINCHI,INCHI_STRING_PREFIX,LEN_INCHI_STRING_PREFIX))    
        return INCHI_INVALID_PREFIX;

    /* .. has InChI version 1 */
    /* if (!isdigit(szINCHI[LEN_INCHI_STRING_PREFIX]) )  */
    if ( szINCHI[LEN_INCHI_STRING_PREFIX] != '1' )  
        return INCHI_INVALID_VERSION;

    /* .. optionally has a 'standard' flag character */
    pos_slash1 = LEN_INCHI_STRING_PREFIX+1;
    if (szINCHI[pos_slash1]=='S')            
    {
        /* Standard InChI ==> standard InChIKey */
        ret = INCHI_VALID_STANDARD;
        pos_slash1++;
    }

    /* .. has trailing slash in the right place */
    if (szINCHI[pos_slash1]!='/')            
        return INCHI_INVALID_LAYOUT;

    /* .. the rest of source string contains valid literals */
#if 0
    if (!isalnum(szINCHI[pos_slash1+1] ) &&
       ( szINCHI[pos_slash1+1]!='/' )     ) 
        return INCHI_INVALID_LAYOUT;    
#endif
    
    /* Treat possible SaveOpt letters  */
    slen0 = slen;
    if ( (szINCHI[slen-3]=='\\') &&
         (szINCHI[slen-2] >= 'A') && (szINCHI[slen-2] <='Z') &&
         (szINCHI[slen-1] >= 'A') && (szINCHI[slen-1] <='Z')		 
        )
        slen0 = slen -3;

    for (i=pos_slash1+1; i<slen0; i++)
    {
        pp = szINCHI[i];
        if (pp >= 'A' && pp <='Z')   continue; 
        if (pp >= 'a' && pp <='z')   continue; 
        if (pp >= '0' && pp <='9')   continue;             
        switch ( pp ) 
        { 
            case '(': case ')': 
            case '*': case '+': 
            case ',': case '-': 
            case '.': case '/': 
            case ';': case '=': 
            case '?': case '@': continue;             
            
            default:            break; 
        }
        return INCHI_INVALID_LAYOUT;
    } 

    if ( strict )
    {
        char opts[]="?FixedH ?RecMet ?SUU ?SLUUD";
        extract_inchi_substring(&str, szINCHI, slen);
        if (NULL==str)
        {
            ret = INCHI_FAIL_I2I;
            goto fin; 
        }

        inchi_inp.szInChI = str;
        opts[0] = opts[8] = opts[16] = opts[21] = INCHI_OPTION_PREFX;
        inchi_inp.szOptions  = opts;

        ret_i2i = GetINCHIfromINCHI(&inchi_inp, &inchi_out);
        
        if ( ((ret_i2i!=inchi_Ret_OKAY) && (ret_i2i!=inchi_Ret_WARNING)) || !inchi_out.szInChI )
        {
            ret = INCHI_FAIL_I2I;
        }
        else
        {
            if (strcmp(inchi_inp.szInChI, inchi_out.szInChI))
            {
                ret = INCHI_FAIL_I2I;
            }
        }
    }

fin:if ( strict )
    {
        if (NULL!=str)      
            inchi_free(str);
    }
    return ret;
}



/*************************************************************************/
/******************************** from readmol.c *************************/
/*************************************************************************/
int AddMOLfileError( char *pStrErr, const char *szMsg )
{
    if ( pStrErr && szMsg && szMsg[0] ) {
        int lenStrErr = strlen( pStrErr );
        int lenMsg    = strlen( szMsg );
        char *p = strstr( pStrErr, szMsg );
        if ( p && (p==pStrErr || (*(p-1) == ' ' && (*(p-2) == ';' || *(p-2) == ':' ))) &&
                  (p+lenMsg == pStrErr+lenStrErr || 
                  (p[lenMsg] == ';' && p[lenMsg+1] == ' ') ||
                  (p[lenMsg-1]==':' && p[lenMsg]==' ')) ) {
            return 1; /*  reject duplicates */
        }
        if ( lenStrErr + lenMsg + 2*(lenStrErr > 0) < STR_ERR_LEN ) {
            /*  enough room to add */
            if (lenStrErr > 0) {
                if ( pStrErr[lenStrErr-1] != ':' ) {
                    strcat( pStrErr, ";" );
                }
                strcat( pStrErr, " " );
            }
            strcat( pStrErr, szMsg );
            return 1;
        }
        /*  no room */
        if ( strstr( pStrErr, "..." ) ) {
            return 0; /*  no room mark has already been set */
        }
        if ( lenStrErr + 3 < STR_ERR_LEN ) {
            strcat( pStrErr, "..." );
        }
    }
    return 0;
}
/****************************************************************/
int CopyMOLfile(FILE *inp_file, long fPtrStart, long fPtrEnd, 
                FILE *prb_file,
                long lNumb)
{
    return 0; /* dummy */
}
/****************************************************************/
/************************** from mol2atom.c *********************/
/****************************************************************/
void SetNumImplicitH(inp_ATOM* at, int num_atoms)
{
    int bNonMetal;
    int a1/*, n1*/;

    /* special valences */
    for ( bNonMetal = 0; bNonMetal < 2; bNonMetal ++ ) {
        for ( a1 = 0; a1 < num_atoms; a1 ++ ) {
            int bHasMetalNeighbor /*, j*/;
            if ( bNonMetal != is_el_a_metal( at[a1].el_number ) ) {
                continue; /* first process all metals, after that all non-metals */
            }

            bHasMetalNeighbor = 0;
            /***********************************************************************
             *  Set number of hydrogen atoms
             */
            at[a1].num_H = get_num_H( at[a1].elname, at[a1].num_H, at[a1].num_iso_H,
                                    at[a1].charge, at[a1].radical,
                                    at[a1].chem_bonds_valence,
                                    0, /* instead of valence entered by the user: it does not exist here*/
                                    (at[a1].at_type & 1)  /* bAliased */,
                                    !(at[a1].at_type & 2) /* bDoNotAddH */,
                                    bHasMetalNeighbor );
            at[a1].at_type = 0;
        }
    }
}

/******************************************************************************************************/
void FreeInpAtom( inp_ATOM **at )
{
    if ( at && *at ) {
        inchi_free( *at );
        *at = NULL;
    }
}
/******************************************************************************************************/
inp_ATOM *CreateInpAtom( int num_atoms )
{
    /*
    void *p = inchi_calloc(num_atoms, sizeof(inp_ATOM) );
    if ( p == (void*)0x009143A8 ) {
        int stop = 1;
    }
    return (inp_ATOM* )p;
    */
   return (inp_ATOM* ) inchi_calloc(num_atoms, sizeof(inp_ATOM) );
}
/******************************************************************************************************/
void FreeInpAtomData( INP_ATOM_DATA *inp_at_data )
{
    if ( inp_at_data ) {
        if ( inp_at_data->at ) {
            FreeInpAtom( &inp_at_data->at );
        }
        if ( inp_at_data->at_fixed_bonds ) {
            FreeInpAtom( &inp_at_data->at_fixed_bonds );
        }
        memset( inp_at_data, 0, sizeof(*inp_at_data) );
    }
}
/******************************************************************************************************/
int CreateInpAtomData( INP_ATOM_DATA *inp_at_data, int num_atoms, int create_at_fixed_bonds )
{
    FreeInpAtomData( inp_at_data );
    if ( (inp_at_data->at = CreateInpAtom( num_atoms )) &&
         (!create_at_fixed_bonds || (inp_at_data->at_fixed_bonds = CreateInpAtom( num_atoms) ) ) ) {
        inp_at_data->num_at = num_atoms;
        return 1;
    }
    FreeInpAtomData( inp_at_data );
    return 0;
}
/******************************************************************************************************/
void FreeCompAtomData( COMP_ATOM_DATA *inp_at_data )
{
    FreeInpAtom( &inp_at_data->at );
    if ( inp_at_data->nOffsetAtAndH )
        inchi_free( inp_at_data->nOffsetAtAndH );
    memset( inp_at_data, 0, sizeof(*inp_at_data) );
}
/******************************************************************************************************/
#if( TEST_RENUMB_ATOMS == 1 )  /*  { */
/******************************************************************************************************/
int CopyInpAtomData( INP_ATOM_DATA *dest_inp_at_data, INP_ATOM_DATA *src_inp_at_data )
{
    int ret = 1;
    if ( !dest_inp_at_data->at  || dest_inp_at_data->num_at != src_inp_at_data->num_at ) {
        ret = CreateInpAtomData( dest_inp_at_data, src_inp_at_data->num_at, (NULL != src_inp_at_data->at_fixed_bonds) );
    } else {
        inp_ATOM *at  = dest_inp_at_data->at;  /*  save ptr to already allocated memory */
        inp_ATOM *at2 = dest_inp_at_data->at_fixed_bonds;
        *dest_inp_at_data = *src_inp_at_data; /*  copy all other (scalar) data */
        dest_inp_at_data->at = at;            /*  restore ptr to already allocated memory */
        dest_inp_at_data->at_fixed_bonds = at2;
    }
    if ( ret ) {
        memcpy( dest_inp_at_data->at, src_inp_at_data->at,
                src_inp_at_data->num_at*sizeof(dest_inp_at_data->at[0]) );
        if ( dest_inp_at_data->at_fixed_bonds && src_inp_at_data->at_fixed_bonds ) {
            memcpy( dest_inp_at_data->at_fixed_bonds, src_inp_at_data->at_fixed_bonds,
                src_inp_at_data->num_at*sizeof(dest_inp_at_data->at_fixed_bonds[0]) );
        }
    }
    return ret;
}
/******************************************************************************************************/
void RenumbInpAtomData( INP_ATOM_DATA *dest_inp_at_data, INP_ATOM_DATA *src_inp_at_data, AT_RANK *new_ord )
{
    int j, n, m, val;
#if( TEST_RENUMB_NEIGH == 1 )
    int i, k;
#endif
    int       num_atoms = src_inp_at_data->num_at;
    inp_ATOM *dest_at   = dest_inp_at_data->at;
    for ( n = 0; n < num_atoms; n ++ ) {
        m = new_ord[n];
        dest_at[m] = src_inp_at_data->at[n];
        dest_at[m].orig_compt_at_numb = (AT_NUMB)(m+1);  /*  new ordering number within the component */
        val = dest_at[m].valence;
        for ( j = 0; j < val; j ++ ) {
            dest_at[m].neighbor[j] = new_ord[dest_at[m].neighbor[j]];
        }
#if( TEST_RENUMB_NEIGH == 1 )
        for ( i = 0; i < 3*val; i ++ ) {
            j = (rand() * val) / (RAND_MAX+1);
            k = (rand() * val) / (RAND_MAX+1);
            if ( j >= val || k >= val || j == k ) {
                continue;
            }
            inchi_swap( (char*)&dest_at[m].neighbor[j],    (char*)&dest_at[m].neighbor[k],    sizeof(dest_at[0].neighbor[0]) );
            inchi_swap( (char*)&dest_at[m].bond_stereo[j], (char*)&dest_at[m].bond_stereo[k], sizeof(dest_at[0].bond_stereo[0]) );
            inchi_swap( (char*)&dest_at[m].bond_type[j],   (char*)&dest_at[m].bond_type[k],   sizeof(dest_at[0].bond_type[0]) );
            /* adjust stereo bond links */
            if ( dest_at[m].sb_parity[0] ) {
                int a;
                for ( a = 0; a < MAX_NUM_STEREO_BONDS && dest_at[m].sb_parity[a]; a ++ ) {
                    
                    if ( k == (int)dest_at[m].sb_ord[a] ) {
                        dest_at[m].sb_ord[a] = j;
                    } else
                    if ( j == (int)dest_at[m].sb_ord[a] ) {
                        dest_at[m].sb_ord[a] = k;
                    }

                    if ( k == (int)dest_at[m].sn_ord[a] ) {
                        dest_at[m].sn_ord[a] = j;
                    } else
                    if ( j == (int)dest_at[m].sn_ord[a] ) {
                        dest_at[m].sn_ord[a] = k;
                    }
                }
            }
        }
#endif
    }

}
/******************************************************************************************************/
void MakeNewOrd( int num_atoms, AT_RANK *new_ord )
{
    int i, j, k;
    for ( i = 0; i < 3*num_atoms; i ++ ) {
        j = (rand() * num_atoms) / (RAND_MAX+1);
        k = (rand() * num_atoms) / (RAND_MAX+1);
        if ( j >= num_atoms || k >= num_atoms || j == k ) {
            continue;
        }
        inchi_swap( (char*)&new_ord[j], (char*)&new_ord[k], sizeof(new_ord[0]) );
    }
}
#endif /*  } TEST_RENUMB_ATOMS == 1  */
/**********************************************************************************/
void FreeOrigAtData( ORIG_ATOM_DATA *orig_at_data )
{
    if ( !orig_at_data )
        return;
    FreeInpAtom( &orig_at_data->at );
    if ( NULL != orig_at_data->nCurAtLen ) {
        inchi_free( orig_at_data->nCurAtLen );
    }
    if ( NULL != orig_at_data->nOldCompNumber ) {
        inchi_free( orig_at_data->nOldCompNumber );
    }
    if ( NULL != orig_at_data->szCoord ) {
        inchi_free( orig_at_data->szCoord );
    }
    if ( NULL != orig_at_data->nEquLabels ) {
        inchi_free( orig_at_data->nEquLabels );
    }
    if ( NULL != orig_at_data->nSortedOrder ) {
        inchi_free( orig_at_data->nSortedOrder );
    }
    memset( orig_at_data, 0, sizeof(*orig_at_data) );
}
/********************************************************************/

#define REPEAT_ALL  0
/********************************************************************/
int parse_options_string ( char *cmd, const char *argv[], int maxargs )
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
    while( i < maxargs-1 ) {
        /* bypass spaces */
        while ( *p == ' ' || *p == '\t' )
            p ++;
        if ( !*p )
            break;
        /* scan an argument */
        argv[i++] = pArgCurChar = p;     /* store preliminary ptr to arg */
        while ( 1 ) {
            bCopyCharToArg = 1;
            nNumBackSlashes = 0;
            while (*p == '\\') {
                ++p;
                ++nNumBackSlashes;
            }
            /* each pair of backslashes => one backslash; one more backslash => literal quote */
            if ( *p == '\"' ) {
                /* one " found */
                if ( nNumBackSlashes % 2 == 0 ) {
                    if (bInsideQuotes) {
                        if (*(p+1) == '\"') {
                            p++;
                        } else {
                            bCopyCharToArg = 0;
                        }
                    } else {
                        bCopyCharToArg = 0;
                    }
                    bInsideQuotes = !bInsideQuotes;
                }
                nNumBackSlashes /= 2;          /* divide nNumBackSlashes by two */
            }
            while (nNumBackSlashes--) {
                *pArgCurChar++ = '\\';
            }
            if (!*p) {
                break;
            }
            if (!bInsideQuotes && (*p == ' ' || *p == '\t')) {
                p ++; 
                /* move to the next char because this char may become
                 * zero due to  *pArgCurChar++ = '\0'; line below */
                break;
            }
            if (bCopyCharToArg) {
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
/*****************************************************************/
#define MIN_BOND_LENGTH   (1.0e-6)
int SetAtomProperties( inp_ATOM *at, MOL_COORD *szCoord, inchi_Atom *ati, int a1, int *nDim, char *pStrErr, int *err )
{        
    S_CHAR      cRadical;
    /* element, check later */

    strcpy( at[a1].elname, ati[a1].elname );

    /* charge */

    at[a1].charge = ati[a1].charge;

    /* radical */

    switch ( ati[a1].radical ) {
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
            while ( nRad > RADICAL_TRIPLET ) {
                nRad -= 2;
            }
            sprintf( szRadicalType, "%d->%d", ati[a1].radical, nRad );
            MOLFILE_ERR_SET (*err, 0, "Radical center type replaced:");
            MOLFILE_ERR_SET (*err, 0, szRadicalType);
            cRadical = nRad;
            if ( nRad < 0 ) {
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

    if ( szCoord ) {
        /* store text coordinates */
        char str[32];
        MOL_COORD * coord_p = szCoord + a1;
        WriteCoord( str, ati[a1].x );
        memcpy( *coord_p, str, 10 );
        WriteCoord( str, ati[a1].y );
        memcpy( *coord_p+10, str, 10 );
        WriteCoord( str, ati[a1].z );
        memcpy( *coord_p+20, str, 10 );
    }

    if ( MIN_BOND_LENGTH < fabs(ati[a1].x) || MIN_BOND_LENGTH < fabs(ati[a1].y) || MIN_BOND_LENGTH < fabs(ati[a1].z) ) {
        if ( MIN_BOND_LENGTH < fabs(ati[a1].z) ) {
            *nDim |= 3;
        } else {
            *nDim |= 2;
        }
    }

    /* orig. at. number */
    at[a1].orig_at_number = a1+1;
    return 0;
#undef MIN_BOND_LENGTH
}
/*********************************************************************/
int SetBondProperties( inp_ATOM *at, inchi_Atom *ati, int a1, int j,
                       int nNumAtoms, int *nNumBonds, char *pStrErr, int *err )
{ 
    int a2;
    S_CHAR     cBondType, cStereoType1, cStereoType2;
    AT_NUMB   *p1, *p2;
    int        n1, n2;

    /* bond type */
    switch( ati[a1].bond_type[j] ) {
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
        MOLFILE_ERR_SET (*err, 0, "Unrecognized bond type:");
        MOLFILE_ERR_SET (*err, 0, szBondType);
        *err |= 8; /*  Unrecognized Bond type replaced with single bond */
        cBondType = BOND_TYPE_SINGLE;
        }
        break;
    }

    /* 2D stereo */

    switch( ati[a1].bond_stereo[j] ) {
    /* stereocenter-related; positive: the sharp end points to this atom  */
    case   INCHI_BOND_STEREO_NONE:
        cStereoType1 = 0;
        cStereoType2 = 0;
        break;
    case   INCHI_BOND_STEREO_SINGLE_1UP:
        cStereoType1 =  STEREO_SNGL_UP;
        cStereoType2 = -STEREO_SNGL_UP;
        break;
    case   INCHI_BOND_STEREO_SINGLE_1EITHER:
        cStereoType1 =  STEREO_SNGL_EITHER;
        cStereoType2 = -STEREO_SNGL_EITHER;
        break;
    case   INCHI_BOND_STEREO_SINGLE_1DOWN:
        cStereoType1 =  STEREO_SNGL_DOWN;
        cStereoType2 = -STEREO_SNGL_DOWN;
        break;
    /* stereocenter-related; negative: the sharp end points to the opposite atom  */
    case   INCHI_BOND_STEREO_SINGLE_2UP:
        cStereoType1 = -STEREO_SNGL_UP;
        cStereoType2 =  STEREO_SNGL_UP;
        break;
    case   INCHI_BOND_STEREO_SINGLE_2EITHER:
        cStereoType1 = -STEREO_SNGL_EITHER;
        cStereoType2 =  STEREO_SNGL_EITHER;
        break;
    case   INCHI_BOND_STEREO_SINGLE_2DOWN:
        cStereoType1 = -STEREO_SNGL_DOWN;
        cStereoType2 =  STEREO_SNGL_DOWN;
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
        MOLFILE_ERR_SET (*err, 0, "Unrecognized bond stereo:");
        MOLFILE_ERR_SET (*err, 0, szBondType);
        *err |= 8; /*  Unrecognized Bond stereo replaced with non-stereo bond */
        cStereoType1 = 0;
        cStereoType2 = 0;
        }
        break;
    }

    /* neighbor */

    if ( ati[a1].neighbor[j] < 0 || ati[a1].neighbor[j] >= nNumAtoms ) {
        *err |= 1; /*  bond for impossible atom number(s); ignored */
        MOLFILE_ERR_SET (*err, 0, "Bond to nonexistent atom");
        goto err_exit;
    }
    a2 = (AT_NUMB) ati[a1].neighbor[j];
    if ( a2 == a1 ) {
        *err |= 1; /*  bond for impossible atom number(s); ignored */
        MOLFILE_ERR_SET (*err, 0, "Atom has a bond to itself");
        goto err_exit;
    }

    /* consistency check; locate the bond in the opposite atom */

    p1 = is_in_the_list( at[a1].neighbor, (AT_NUMB)a2, at[a1].valence );
    p2 = is_in_the_list( at[a2].neighbor, (AT_NUMB)a1, at[a2].valence );
    if ( p1 && p2 ) {
        n1 = (p1 - at[a1].neighbor);
        n2 = (p2 - at[a2].neighbor);
        if ( (n1+1 < at[a1].valence && is_in_the_list( at[a1].neighbor+n1+1, (AT_NUMB)a2, at[a1].valence-n1-1 )) ||
             (n2+1 < at[a2].valence && is_in_the_list( at[a2].neighbor+n2+1, (AT_NUMB)a1, at[a2].valence-n2-1 )) ) {
            MOLFILE_ERR_SET (*err, 0, "Multiple bonds between two atoms");
            *err |= 2; /*  multiple bonds between atoms */
        } else
        if ( n1 < at[a1].valence && n2 < at[a2].valence &&
             cBondType == at[a2].bond_type[n2] &&
             cBondType == at[a1].bond_type[n1] &&
             cStereoType1 == at[a1].bond_stereo[n1] &&
             cStereoType2 == at[a2].bond_stereo[n2] ) {
            /*MOLFILE_ERR_SET (*err, 0, "Duplicated bond(s) between two atoms");*/
        } else {
            MOLFILE_ERR_SET (*err, 0, "Multiple bonds between two atoms");
            *err |= 2; /*  multiple bonds between atoms */
        }
    } else
    if ( (p1 || p2) && (p1 || at[a1].valence < MAXVAL) && (p2 || at[a2].valence < MAXVAL) ) {
        n1 = p1? (p1 - at[a1].neighbor) : at[a1].valence ++;
        n2 = p2? (p2 - at[a2].neighbor) : at[a2].valence ++;
        /* the bond is present in one atom only: possibly program error */
        if ( (p1 && (cBondType != at[a1].bond_type[n1] || at[a1].bond_stereo[n1] != cStereoType1)) ||
             (p2 && (cBondType != at[a2].bond_type[n2] || at[a2].bond_stereo[n2] != cStereoType2)) ) {
            MOLFILE_ERR_SET (*err, 0, "Multiple bonds between two atoms");
            *err |= 2; /*  multiple bonds between atoms */
        } else {
            MOLFILE_ERR_SET (*err, 0, "Duplicated bond(s) between two atoms");
            /* warning */
        }
    } else
    if ( !p1 && !p2 && at[a1].valence < MAXVAL && at[a2].valence < MAXVAL ) {
        n1 = at[a1].valence ++;
        n2 = at[a2].valence ++;
        (*nNumBonds) ++;
    } else {
        char szMsg[64];
        *err |= 4; /*  too large number of bonds. Some bonds ignored. */
        sprintf( szMsg, "Atom '%s' has more than %d bonds",
                        at[a1].valence>= MAXVAL? at[a1].elname:at[a2].elname, MAXVAL );
        MOLFILE_ERR_SET (*err, 0, szMsg);
        goto err_exit;
    }
    
    /* store the connection */

    /* bond type */
    at[a1].bond_type[n1] =
    at[a2].bond_type[n2] = cBondType;
    /* connection */
    at[a1].neighbor[n1] = (AT_NUMB)a2;
    at[a2].neighbor[n2] = (AT_NUMB)a1;
    /* stereo */
    at[a1].bond_stereo[n1] =  cStereoType1; /*  >0: the wedge (pointed) end is at this atom */
    at[a2].bond_stereo[n2] =  cStereoType2; /*  <0: the wedge (pointed) end is at the opposite atom */
    return 0;
err_exit:
    return 1;
}
/******************************************************************/
int SetAtomAndBondProperties( inp_ATOM *at, inchi_Atom *ati, int a1,
                              int bDoNotAddH, char *pStrErr, int *err )
{
    int valence, chem_valence, num_alt_bonds, j, n1;
    int nRadical, nCharge;
    static int el_number_H = 0;
    
    if ( !el_number_H ) {
        el_number_H = get_periodic_table_number( "H" );
    }

    nRadical = nCharge = 0;
    valence = at[a1].valence;
    chem_valence = num_alt_bonds = 0;
    for ( j = 0; j < valence; j ++ ) {
        if ( at[a1].bond_type[j] <= BOND_TYPE_TRIPLE ) {
            chem_valence += at[a1].bond_type[j];
        } else {
            num_alt_bonds ++;
        }
    }
    switch( num_alt_bonds ) {
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
        sprintf( szMsg, "Atom '%s' has %d alternating bonds",
                        at[a1].elname, num_alt_bonds );
        MOLFILE_ERR_SET (*err, 0, szMsg);
        }
        break;
    }
    at[a1].chem_bonds_valence = chem_valence;

    /* aliased hydrogen atoms */
    if ( ERR_ELEM == (n1 = get_periodic_table_number( at[a1].elname ) ) ) {
        /*  Case when elname contains more than 1 element: extract number of H if possible */
        if ( extract_ChargeRadical( at[a1].elname, &nRadical, &nCharge ) ) {
            if ( (nRadical && at[a1].radical && nRadical != at[a1].radical) ||
                 (nCharge  && at[a1].charge  && nCharge  != at[a1].charge) ) {
                MOLFILE_ERR_SET (*err, 0, "Ignored charge/radical redefinition:");
                MOLFILE_ERR_SET (*err, 0, ati[a1].elname);
            } else {
                if ( nRadical )
                    at[a1].radical = nRadical;
                if ( nCharge )
                    at[a1].charge  = nCharge;
            }
        }
        at[a1].num_H = extract_H_atoms( at[a1].elname, at[a1].num_iso_H );
        if ( !at[a1].elname[0] && NUMH(at, a1) ) {
            /* alias contains only H. Added 2004-07-21, fixed 2004-07-22
             * move the heaviest isotope to the "central atom"
             * Note: this must be consistent with H-H treatment in remove_terminal_HDT()
             */
            strcpy( at[a1].elname, "H" );
            if ( NUM_ISO_H(at,a1) ) {
                for ( j = NUM_H_ISOTOPES-1; 0 <= j; j -- ) {
                    if ( at[a1].num_iso_H[j] ) {
                        at[a1].num_iso_H[j] --;
                        at[a1].iso_atw_diff = 1 + j;
                        break;
                    }
                }
            } else {
                at[a1].num_H --;
            }
        }
        if ( ERR_ELEM == (n1 = get_periodic_table_number( at[a1].elname ) ) ) {
            n1 = 0;
        }
        if ( n1 ) {
            at[a1].at_type |= 1; /* "Aliased" atom: data in the element name */
            MOLFILE_ERR_SET (*err, 0, "Parsed compound atom(s):");
            MOLFILE_ERR_SET (*err, 0, ati[a1].elname);
        }
    }

    at[a1].el_number = (U_CHAR) n1;
    if ( !n1 ) {
        *err |= 64; /*  Unrecognized aromatic bond(s) replaced with single */
        MOLFILE_ERR_SET (*err, 0, "Unknown element(s):");
        MOLFILE_ERR_SET (*err, 0, at[a1].elname);
    } else
    /* replace explicit D or T with isotopic H (added 2003-06-02) */
    if ( el_number_H == n1 && !at[a1].iso_atw_diff ) {
        switch( at[a1].elname[0] ) {
        case 'D':
            at[a1].iso_atw_diff = 2;
            mystrncpy( at[a1].elname, "H", sizeof(at->elname) );
            break;
        case 'T':
            at[a1].iso_atw_diff = 3;
            mystrncpy( at[a1].elname, "H", sizeof(at->elname) );
            break;
        case 'H':
            if ( 1 <= ati[a1].isotopic_mass ) {
                AT_NUM iso_atw_diff;
                if ( ISOTOPIC_SHIFT_FLAG - ISOTOPIC_SHIFT_MAX <=  ati[a1].isotopic_mass &&
                     ISOTOPIC_SHIFT_FLAG + ISOTOPIC_SHIFT_MAX >=  ati[a1].isotopic_mass ) {
                    /* ati[a1].isotopic_mass is isotopic iso_atw_diff + ISOTOPIC_SHIFT_FLAG */
                    iso_atw_diff = ati[a1].isotopic_mass - ISOTOPIC_SHIFT_FLAG;
                } else {
                    /* ati[a1].isotopic_mass is isotopic mass */
                    iso_atw_diff = get_atw_from_elnum( (int) at[a1].el_number );
                    iso_atw_diff = ati[a1].isotopic_mass - iso_atw_diff;
                }
                if ( iso_atw_diff >= 0 )
                    iso_atw_diff ++;
                /* reproduce Bug04: allowed non-terminal H heavier than T */
                if ( 1 <= iso_atw_diff &&
                     (at[a1].valence != 1 || iso_atw_diff <= NUM_H_ISOTOPES) ) {
                    at[a1].iso_atw_diff = (S_CHAR)iso_atw_diff;
                }
            }
        }
    } else
    /* isotopic shift */
    if ( ati[a1].isotopic_mass ) {
        AT_NUM iso_atw_diff;
        if ( ISOTOPIC_SHIFT_FLAG - ISOTOPIC_SHIFT_MAX <=  ati[a1].isotopic_mass &&
             ISOTOPIC_SHIFT_FLAG + ISOTOPIC_SHIFT_MAX >=  ati[a1].isotopic_mass ) {
            /* ati[a1].isotopic_mass is isotopic iso_atw_diff + ISOTOPIC_SHIFT_FLAG */
            iso_atw_diff = ati[a1].isotopic_mass - ISOTOPIC_SHIFT_FLAG;
        } else {
            /* ati[a1].isotopic_mass is isotopic mass */
            iso_atw_diff = get_atw_from_elnum( (int) at[a1].el_number );
            iso_atw_diff = ati[a1].isotopic_mass - iso_atw_diff;
        }
        if ( iso_atw_diff >= 0 )
            iso_atw_diff ++;
        at[a1].iso_atw_diff = (S_CHAR)iso_atw_diff;
    }

    /* add implicit hydrogen atoms flag */

    if ( ati[a1].num_iso_H[0] == -1 ) {
        if ( !bDoNotAddH ) {
            at[a1].at_type |= 2; /* user requested to add H */
        }
    } else {
        at[a1].num_H = ati[a1].num_iso_H[0];
    }
    for ( j = 0; j < NUM_H_ISOTOPES; j ++ ) {
        at[a1].num_iso_H[j] = ati[a1].num_iso_H[j+1];
    }
    if ( num_alt_bonds ) {
        /* atom has aromatic bonds AND the chemical valence is not known */
        int num_H = NUMH(at, a1);
        int chem_valence_alt = at[a1].chem_bonds_valence + num_H;
        int bUnusualValenceArom = 
            detect_unusual_el_valence( (int)at[a1].el_number, at[a1].charge,
                                        at[a1].radical, chem_valence_alt,
                                        num_H, at[a1].valence );
        int bUnusualValenceNoArom = 
            detect_unusual_el_valence( (int)at[a1].el_number, at[a1].charge,
                                        at[a1].radical, chem_valence_alt-1,
                                        num_H, at[a1].valence );
        if ( bUnusualValenceArom && !bUnusualValenceNoArom && 0 == nBondsValToMetal( at, a1) ) {
            /* typically NH in 5-member aromatic ring */
            at[a1].chem_bonds_valence --;
        }
    }

    return 0;
}
/****************************************************************************************/
int InpAtom0DToInchiAtom( inp_ATOM *at, int num_atoms, inchi_OutputStruct *outStruct )
{
    int num_stereo_centers, num_stereo_bonds, num_stereo0D, i, m, m1, m2, n, ret=0;
    /* count stereobonds, allenes. cumulenes. and stereoatoms */
    num_stereo_centers = num_stereo_bonds = ret = 0;
    
    outStruct->atom = NULL;
    outStruct->num_atoms = 0;
    outStruct->stereo0D = NULL;
    outStruct->num_stereo0D = 0;

    for ( i = 0; i < num_atoms; i ++ ) {
        if ( at[i].p_parity ) {
            /* stereocenter */
            num_stereo_centers ++;
        } else {
            for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m]; m ++ )
                ;
            num_stereo_bonds += m;
        }
    }
    num_stereo_bonds /= 2;
    num_stereo0D = num_stereo_bonds + num_stereo_centers;

    if ( num_atoms > 0 ) {
        outStruct->atom = (inchi_Atom *)inchi_calloc( num_atoms, sizeof( outStruct->atom[0] ) );
    }
    outStruct->num_atoms = num_atoms;
    if ( num_stereo0D > 0 ) {
        outStruct->stereo0D = (inchi_Stereo0D *)inchi_calloc( num_stereo0D, sizeof(outStruct->stereo0D[0]));
    }
    if ( (num_atoms && !outStruct->atom) || (num_stereo0D > 0 && !outStruct->stereo0D) ) {
        /* allocation failed */
        ret = -1;
        goto exit_function;
    }

    /* copy atom properties */
    for ( i = 0; i < num_atoms; i ++ ) {
        outStruct->atom[i].num_bonds = at[i].valence;
        for ( m = 0; m < at[i].valence; m ++ ) {
            outStruct->atom[i].bond_type[m] = at[i].bond_type[m];
            outStruct->atom[i].neighbor[m]  = at[i].neighbor[m];
        }
        outStruct->atom[i].charge = at[i].charge;
        memcpy( outStruct->atom[i].elname, at[i].elname, ATOM_EL_LEN );
        if ( at[i].iso_atw_diff ) {
            outStruct->atom[i].isotopic_mass = ISOTOPIC_SHIFT_FLAG + (at[i].iso_atw_diff > 0? at[i].iso_atw_diff-1 : at[i].iso_atw_diff);
        }
        outStruct->atom[i].num_iso_H[0] = at[i].num_H;
        for ( m = 0; m < NUM_H_ISOTOPES; m ++ ) {
            outStruct->atom[i].num_iso_H[m+1] = at[i].num_iso_H[m];
        }
        outStruct->atom[i].radical = at[i].radical;
    }
    /* stereo */
    for ( i = n = 0; i < num_atoms; i ++ ) {
        if ( at[i].p_parity ) {
            if ( n < num_stereo0D ) {
                outStruct->stereo0D[n].central_atom = i;
                outStruct->stereo0D[n].parity       = at[i].p_parity;
                outStruct->stereo0D[n].type         = INCHI_StereoType_Tetrahedral;
                for ( m = 0; m < MAX_NUM_STEREO_ATOM_NEIGH; m ++ ) {
                    outStruct->stereo0D[n].neighbor[m] = at[i].p_orig_at_num[m] - 1;
                }
                n ++;
            } else {
                ret |= 1;
                break;
            }
        } else {
            for ( m1 = 0; m1 < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m1]; m1 ++ ) {
                
                /* find the opposite atom at the other end of double bond, allene, or cumulene */
                int chain[12], len = 0, nxt_neigh, nxt, cur;
                cur = chain[len++] = i;
                nxt_neigh = at[cur].sb_ord[m1];
                
                do {
                    /* add next atom */
                    chain[len ++] = nxt = at[cur].neighbor[nxt_neigh];
                    nxt_neigh = (at[nxt].neighbor[0] == cur);
                    cur = nxt;
                    /* find nxt_neigh */
                } while ( !at[cur].sb_parity[0] && len < 12 && at[cur].valence == 2 );

                if ( at[cur].sb_parity[0] && len <= 4 && i < cur /* count bonds only one time */ ) {
                    /* double bond, cumulene, or allene has been found */
                    for ( m2 = 0; m2 < MAX_NUM_STEREO_BONDS && at[cur].sb_parity[m2]; m2 ++ ) {
                        if ( chain[len-2] == at[cur].neighbor[(int)at[cur].sb_ord[m2]] ) {
                            if ( n < num_stereo0D ) {
                                int parity1 = at[i].sb_parity[m1];
                                int parity2 = at[cur].sb_parity[m2];
                                int parity;
                                if ( (INCHI_PARITY_ODD == parity1 || INCHI_PARITY_EVEN == parity1) &&
                                    (INCHI_PARITY_ODD == parity2 || INCHI_PARITY_EVEN == parity2) ) {
                                    /* well-defined parity */
                                    parity = (parity1==parity2)? INCHI_PARITY_EVEN : INCHI_PARITY_ODD;
                                } else {
                                    parity = inchi_max(parity1, parity2);
                                }
                                outStruct->stereo0D[n].central_atom = (len==3)? chain[1] : NO_ATOM;
                                outStruct->stereo0D[n].parity       = parity;
                                outStruct->stereo0D[n].type         = len == 3? INCHI_StereoType_Allene : INCHI_StereoType_DoubleBond;
                                outStruct->stereo0D[n].neighbor[0]  = at[i].sn_orig_at_num[m1]-1;
                                outStruct->stereo0D[n].neighbor[1]  = i;
                                outStruct->stereo0D[n].neighbor[2]  = cur;
                                outStruct->stereo0D[n].neighbor[3]  = at[cur].sn_orig_at_num[m2] - 1;
                                n ++;
                            } else {
                                ret |= 1;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    outStruct->num_stereo0D = n;
exit_function:
    if ( ret < 0 ) {
        if ( outStruct->atom ) inchi_free( outStruct->atom );
        if ( outStruct->stereo0D ) inchi_free( outStruct->stereo0D );
        outStruct->atom = NULL;
        outStruct->stereo0D = NULL;
        outStruct->num_atoms = 0;
        outStruct->num_stereo0D = 0;
    }
    return ret;
}
/****************************************************************************************/
int ExtractOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
         inchi_Input *inp, 
         INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file, INCHI_IOSTREAM *prb_file,
         ORIG_ATOM_DATA *orig_inp_data, long *num_inp,
         char *pStr, int nStrLen )
{
    int         *err           = &sd->nStructReadError;
    char        *pStrErr       = sd->pStrErrStruct;
    inp_ATOM    *at            = NULL;
    MOL_COORD   *szCoord       = NULL; 
    inchi_Atom  *ati           = NULL;  
    int       nNumAtoms = 0;
    int       a1, j, valence, nDim, nNumBonds, nRet = 0;

    /* vABParityUnknown holds actual value of an internal constant signifying       */
    /* unknown parity: either the same as for undefined parity (default==standard)  */
    /*  or a specific one (non-std; requested by SLUUD switch).                     */
    int vABParityUnknown = AB_PARITY_UNDF;
    if ( 0 != ( ip->nMode & REQ_MODE_DIFF_UU_STEREO) ) 
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
    nDim      = 0;
    nNumBonds = 0;

    if ( !inp || (nNumAtoms = inp->num_atoms) <= 0 || !(ati = inp->atom) ) {
        MOLFILE_ERR_SET (*err, 0, "Empty structure");
        *err = 98;
        goto err_exit;
    }
    if ( nNumAtoms >= MAX_ATOMS ) {
        MOLFILE_ERR_SET (*err, 0, "Too many atoms");
        *err = 70;
        orig_inp_data->num_inp_atoms = -1;
        goto err_exit;
    }

    at      = (inp_ATOM  *) inchi_calloc( nNumAtoms, sizeof(at[0]) );
    szCoord = (MOL_COORD *) inchi_calloc (inchi_max(nNumAtoms, 1), sizeof (MOL_COORD));

    if ( !at || !szCoord ) {
        MOLFILE_ERR_SET (*err, 0, "Out of RAM");
        *err = -1;
        goto err_exit;
    }


    /********************************************************
     *
     *   Extract typical for Molfile structural data
     *
     ********************************************************/
    /* extract atoms and bonds */
    for ( a1 = 0; a1 < nNumAtoms; a1 ++ ) {
        /* extract atoms */
        SetAtomProperties( at, szCoord, ati, a1, &nDim, pStrErr, err );
        if ( *err ) {
            goto err_exit;
        }
        /* extract connections */
        valence = ati[a1].num_bonds;
        for ( j = 0; j < valence; j ++ ) {
            SetBondProperties( at, ati, a1, j, nNumAtoms, &nNumBonds, pStrErr, err );
        }
        if ( *err ) {
            goto err_exit;
        }
    }

    orig_inp_data->num_inp_atoms = nNumAtoms;
    orig_inp_data->num_inp_bonds = nNumBonds;
    orig_inp_data->num_dimensions = nDim;

    /* extract elements, chemical valences, implicit H, isotopic shifts */
    for ( a1 = 0; a1 < nNumAtoms; a1 ++ ) {
        /* set temp flags in at[a1].at_type (1: data in atom name; 2: request to add H) */
        SetAtomAndBondProperties( at, ati, a1, ip->bDoNotAddH, pStrErr, err );
        if ( *err ) {
            goto err_exit;
        }
    }
    /* clear temp flags in at[].at_type; add implicit H */
    SetNumImplicitH( at, nNumAtoms );
    if ( *err ) {
        goto err_exit;
    }
            
    /********************************************************
     *
     *   Extract the 0D parities (typical for CML)
     *
     ********************************************************/
    Extract0DParities(at, nNumAtoms, inp->stereo0D, inp->num_stereo0D, 
                       pStrErr, err, vABParityUnknown);

    if ( *err ) {
        goto err_exit;
    }
    orig_inp_data->at             = at;          at     = NULL;
    orig_inp_data->num_dimensions = nDim;
    orig_inp_data->num_inp_atoms  = nNumAtoms;
    orig_inp_data->num_inp_bonds  = nNumBonds;
    orig_inp_data->szCoord        = szCoord;     szCoord = NULL;

    /* chiral flag */
    /* *****************************************************************************
     * Chiral flags are set in: 
     * - RunICHI.c #1610 -- ReadTheStructure()     -- cInChI, wInChI
     * - e_IchiMain.c #273 -- main()               -- C example of calling InChI dll  
     * - inchi_dll.c  #1662 -- ExtractOneStructure -- InChI dll code (here)
     *******************************************************************************/   
    if ( (ip->nMode & REQ_MODE_CHIR_FLG_STEREO) && (ip->nMode & REQ_MODE_STEREO) ) {
        if ( ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL ) {
            /* absolute stereo */
            ip->nMode &= ~(REQ_MODE_RELATIVE_STEREO | REQ_MODE_RACEMIC_STEREO);
            sd->bChiralFlag &= ~FLAG_INP_AT_NONCHIRAL;
            sd->bChiralFlag |= FLAG_INP_AT_CHIRAL; /* write AuxInfo as chiral */
        } else
        /*if ( ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL )*/ {
            /* relative stereo */
            ip->nMode &= ~(REQ_MODE_RACEMIC_STEREO);
            ip->nMode |=   REQ_MODE_RELATIVE_STEREO;
            sd->bChiralFlag &= ~FLAG_INP_AT_CHIRAL;
            sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL; /* write AuxInfo as non-chiral */
        }
    } else
    if ( ip->bChiralFlag & FLAG_SET_INP_AT_CHIRAL ) {
        sd->bChiralFlag &= ~FLAG_INP_AT_NONCHIRAL;
        sd->bChiralFlag |= FLAG_INP_AT_CHIRAL; /* write AuxInfo as chiral */
    } else
    if ( ip->bChiralFlag & FLAG_SET_INP_AT_NONCHIRAL ) {
        sd->bChiralFlag &= ~FLAG_INP_AT_CHIRAL;
        sd->bChiralFlag |= FLAG_INP_AT_NONCHIRAL; /* write AuxInfo as non-chiral */
    }

    *num_inp += 1;

err_exit:

    if ( at )
        inchi_free( at );
    if ( szCoord )
        inchi_free( szCoord );

    nRet = TreatReadTheStructureErrors( sd, ip, LOG_MASK_NO_WARN, NULL, log_file, output_file, prb_file,
                                                orig_inp_data, num_inp, pStr, nStrLen );

    return nRet;

}
/********************************************************/
int INCHI_DECL GetStringLength( char *p )
{
    if ( p ) {
        return strlen(p);
    } else {
        return 0;
    }
}
#define MAX_MSG_LEN 512



/* GetINCHIfromINCHI does same as -InChI2InChI option: converts InChI into InChI for validation purposes */
/* It may also be used to filter out specific layers. For instance, /Snon would remove stereochemical layer */
/* Omitting /FixedH and/or /RecMet would remove Fixed-H or Reconnected layers */
/* To keep all InChI layers use options string "/FixedH /RecMet"; option /InChI2InChI is not needed */
/* inchi_InputINCHI is created by the user; strings in inchi_Output are allocated and deallocated by InChI */
/* inchi_Output does not need to be initilized out to zeroes; see FreeINCHI() on how to deallocate it */
/*************************************************************/
int INCHI_DECL GetINCHIfromINCHI( inchi_InputINCHI *inpInChI, inchi_Output *out )
{
    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;    
    
    static char szMainOption[] = " ?InChI2InChI";

    int i;
    char      szSdfDataValue[MAX_SDF_VALUE+1];
    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */

    INPUT_PARMS inp_parms;
    INPUT_PARMS *ip = &inp_parms;

    int             bReleaseVersion = bRELEASE_VERSION;
    int   nRet = 0, nRet1;

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif

    const char *argv[INCHI_MAX_NUM_ARG+1];
    int   argc;
    char *szOptions = NULL;

    INCHI_IOSTREAM inchi_file[3], *output_file = inchi_file, *log_file = inchi_file+1, *input_file = inchi_file+2;




    if ( bLibInchiSemaphore ) {  /* does not work properly under sufficient stress */
        return inchi_Ret_BUSY;
    }
    bLibInchiSemaphore = 1;

#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);
/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE 
   _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode(_CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG);
#endif
    
    /* turn on floating point exceptions */
#if ( !defined(__STDC__) || __STDC__ != 1 )
    {
        /* Get the default control word. */
        int cw = _controlfp( 0,0 );

        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_ZERODIVIDE|EM_DENORMAL);

        /* Set the control word. */
        _controlfp( cw, MCW_EM );
 
    }
#endif
#endif



    memset( out, 0, sizeof(*out) );   
#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
repeat:
    FreeINCHI( out );
    inchi_ios_close(output_file);
    inchi_ios_close(log_file);
    inchi_ios_reset(input_file);  /* do not close input_file - its string buffer may point to inpInChI->szInChI */
#endif

    /*^^^ Initialize internal for this function I/O streams as string buffers */
    inchi_ios_init(input_file, INCHI_IOSTREAM_STRING, NULL);
    inchi_ios_init(output_file, INCHI_IOSTREAM_STRING, NULL);
    inchi_ios_init(log_file, INCHI_IOSTREAM_STRING, NULL);


    sd->bUserQuit  = 0;

    /* clear original input structure */
    /*^^^ memset( inchi_file, 0, sizeof(inchi_file) ); */
    memset( sd,         0, sizeof(*sd) );
    memset( ip,         0, sizeof(*ip) );
    memset( szSdfDataValue    , 0, sizeof( szSdfDataValue    ) );
    szMainOption[1] = INCHI_OPTION_PREFX;

    if ( !inpInChI ) {
        nRet = _IS_ERROR;
        goto exit_function;
    }

    /* options */
    if ( inpInChI ) {
        int opt_len = (inpInChI->szOptions? strlen(inpInChI->szOptions) : 0) + sizeof(szMainOption) + 1;
        szOptions = (char*)inchi_calloc( opt_len+1, sizeof(szOptions[0]) );
        if ( szOptions ) {
            if ( inpInChI->szOptions ) {
                strcpy( szOptions, inpInChI->szOptions );
            }
            strcat( szOptions, szMainOption );
            argc = parse_options_string ( szOptions, argv, INCHI_MAX_NUM_ARG );
        } else {
            nRet = _IS_FATAL;
            goto translate_RetVal; /* emergency exit */
        }
    } else {
        argc = 1;
        argv[0] = "";
        argv[1] = NULL;
    }

    if ( (argc == 1
#ifdef TARGET_API_LIB
        && (!inpInChI || !inpInChI->szInChI)
#endif   
        )
        || (argc==2 && ( argv[1][0]==INCHI_OPTION_PREFX ) &&
        (!strcmp(argv[1]+1, "?") || !stricmp(argv[1]+1, "help") )) ) {
        HelpCommandLineParms(log_file);
        out->szLog = log_file->s.pStr;
        memset( log_file, 0, sizeof(*log_file) );
        nRet = _IS_EOF;
        goto translate_RetVal;
    }

    nRet1 = ReadCommandLineParms( argc, argv, ip, szSdfDataValue, &ulDisplTime, bReleaseVersion, log_file );
    if ( szOptions ) {
        /* argv pointed to strings in szOptions */
        inchi_free( szOptions );
        szOptions = NULL;
    }
    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if ( 0 > nRet1 ) {
        goto exit_function;
    }
    if ( ip->bNoStructLabels ) {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    } else
    if ( ip->nInputType == INPUT_INCHI_XML || ip->nInputType == INPUT_INCHI_PLAIN  ||
         ip->nInputType == INPUT_CMLFILE || ip->nInputType == INPUT_INCHI ) {
        /* the input may contain both the header and the label of the structure */
        if ( !ip->pSdfLabel ) 
            ip->pSdfLabel  = ip->szSdfDataHeader;
        if ( !ip->pSdfValue )
            ip->pSdfValue  = szSdfDataValue;
    }
    if ( ip->nInputType && ip->nInputType != INPUT_INCHI ) {
        inchi_ios_eprint( log_file, "Input type set to INPUT_INCHI\n" );
        ip->nInputType = INPUT_INCHI;
    }

    PrintInputParms( log_file, ip );
    /*********************************/
    /* InChI -> Structure conversion */
    /*********************************/

    /* input_file simulation */
    input_file->s.pStr = inpInChI->szInChI;
    input_file->s.nUsedLength = strlen(input_file->s.pStr)+1;
    input_file->s.nAllocatedLength = input_file->s.nUsedLength;
    input_file->s.nPtr = 0;

    /* buffer for the message */
    out->szMessage = (char *)inchi_calloc( MAX_MSG_LEN, sizeof(out->szMessage[0]));
    if ( !out->szMessage ) {
         inchi_ios_eprint( log_file, "Cannot allocate output message buffer.\n");
        nRet = -1;
    } else {
        nRet = ReadWriteInChI( input_file, output_file, log_file,
                               ip,  sd, NULL,  NULL, out->szMessage, MAX_MSG_LEN, NULL /*out->WarningFlags*/ );
    }
    
    if ( nRet >= 0 && output_file->s.pStr ) 
    {
        /* success */
        char *p;
        out->szInChI = output_file->s.pStr;
        out->szAuxInfo = NULL;
        for ( p = strchr(out->szInChI, '\n'); p; p = strchr(p+1, '\n') ) {
            if ( !memcmp( p, "\nAuxInfo", 8 ) ) {
                *p = '\0';            /* remove LF after INChI */
                out->szAuxInfo = p+1; /* save pointer to AuxInfo */
            } else
            if ( out->szAuxInfo || !p[1]) {   /* remove LF after aux info or from the last char */
                *p = '\0';
                break;
            }
        }
        output_file->s.pStr = NULL;
    }
    /*
    out->szLog = log_file->pStr;
    log_file->pStr   = NULL;
    */

exit_function:;

#if( ADD_CMLPP == 1 )
        /* BILLY 8/6/04 */
        /* free CML memory */
        FreeCml ();
        FreeCmlDoc( 1 );
#endif



    for ( i = 0; i < MAX_NUM_PATHS; i ++ ) {
        if ( ip->path[i] ) {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( );


#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if ( num_repeat-- > 0 ) {
        goto repeat;
    }
#endif


#ifdef TARGET_API_LIB
    /* output */

    if ( log_file->s.pStr && log_file->s.nUsedLength > 0 ) {
        while ( log_file->s.nUsedLength && '\n' == log_file->s.pStr[log_file->s.nUsedLength-1] ) {
            log_file->s.pStr[-- log_file->s.nUsedLength]  = '\0'; /* remove last LF */
        }
        if ( out ) {
            out->szLog = log_file->s.pStr;
            log_file->s.pStr = NULL;
        }
    }


#endif

    
translate_RetVal:

    /* Close internal output streams */
    inchi_ios_close(output_file);
    inchi_ios_close(log_file);
    inchi_ios_reset(input_file);  /* do not close input_file - its string buffer may point to inpInChI->szInChI */


    switch (nRet) { 
    case -3         : nRet = inchi_Ret_ERROR  ; break; /* Error: no Structure has been created */
    case -2         : nRet = inchi_Ret_ERROR  ; break; /* Error: no Structure has been created */
    case -1         : nRet = inchi_Ret_FATAL  ; break; /* Severe error: no Structure has been created (typically; break; memory allocation failed) */
    default         :
        /*
        if ( !outStruct->atom || !outStruct->num_atoms ) {
            nRet = inchi_Ret_EOF;
        } else {
            int m,n,t=0;
            for ( m=0; m < 2; m ++ ) {
                for ( n=0; n < 2; n ++ ) {
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

    bLibInchiSemaphore = 0;
    return nRet;
}


EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStructFromStdINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct )
{
    if ( ( inpInChI ) && 
         ( inpInChI->szInChI ) &&
         ( strlen(inpInChI->szInChI) >= LEN_INCHI_STRING_PREFIX+3 ) &&
         ( inpInChI->szInChI[LEN_INCHI_STRING_PREFIX+1] == 'S' ) 
       ) 
        /* brief check indicated valid std input (more checks in GetStructFromINCHI) */
        return GetStructFromINCHI( inpInChI, outStruct );
    else
        /* non-std or just invalid input */
        return inchi_Ret_ERROR;
}



/*************************************************************/
EXPIMP_TEMPLATE INCHI_API int INCHI_DECL GetStructFromINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct )
{
    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;
    

    INCHI_IOSTREAM inchi_file[3];
    INCHI_IOSTREAM *output_file = inchi_file, *log_file = inchi_file+1, *input_file = inchi_file+2;
    

    static char szMainOption[] = " ?InChI2Struct";


    int i;
    char      szSdfDataValue[MAX_SDF_VALUE+1];
    unsigned long  ulDisplTime = 0;    /*  infinite, milliseconds */

    INPUT_PARMS inp_parms;
    INPUT_PARMS *ip = &inp_parms;

    int             bReleaseVersion = bRELEASE_VERSION;
    int   nRet = 0, nRet1;
    int bStdFormat          = 0;

    /* conversion result */
    inp_ATOM *at=NULL;
    int num_at = 0;

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif

    const char *argv[INCHI_MAX_NUM_ARG+1];
    int   argc;
    char *szOptions = NULL;

    if ( bLibInchiSemaphore ) {  /* does not work properly under sufficient stress */
        return inchi_Ret_BUSY;
    }
#if 0
    /* moved to after call to CheckINCHI - Marc 2010 */
    bLibInchiSemaphore = 1;
#endif

#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);
/* for execution outside the VC++ debugger uncomment one of the following two */
#ifdef MY_REPORT_FILE 
   _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_WARN, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ERROR, MY_REPORT_FILE );
   _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ASSERT, MY_REPORT_FILE );
#else
    _CrtSetReportMode(_CRT_WARN | _CRT_ERROR, _CRTDBG_MODE_DEBUG);
#endif
    
    /* turn on floating point exceptions */
#if ( !defined(__STDC__) || __STDC__ != 1 )    
    {
        /* Get the default control word. */
        int cw = _controlfp( 0,0 );

        /* Set the exception masks OFF, turn exceptions on. */
        /*cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_INEXACT|EM_ZERODIVIDE|EM_DENORMAL);*/
        cw &=~(EM_OVERFLOW|EM_UNDERFLOW|EM_ZERODIVIDE|EM_DENORMAL);

        /* Set the control word. */
        _controlfp( cw, MCW_EM );
 
    }
#endif
#endif

    memset( outStruct, 0, sizeof(*outStruct) );

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
repeat:
    FreeStructFromINCHI( &outStruct );
    inchi_ios_reset(input_file);  /* do not close input_file - its string buffer may point to inpInChI->szInChI */
    inchi_ios_close(output_file);
    inchi_ios_close(log_file);
#endif

    sd->bUserQuit  = 0;



    /*^^^ Initialize internal for this function I/O streams as string buffers */
    inchi_ios_init(input_file, INCHI_IOSTREAM_STRING, NULL);
    inchi_ios_init(output_file, INCHI_IOSTREAM_STRING, NULL);
    inchi_ios_init(log_file, INCHI_IOSTREAM_STRING, NULL);


    /* clear original input structure */

    memset( sd,         0, sizeof(*sd) );
    memset( ip,         0, sizeof(*ip) );
    memset( szSdfDataValue    , 0, sizeof( szSdfDataValue    ) );
    szMainOption[1] = INCHI_OPTION_PREFX;

    if ( !inpInChI ) 
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }

    /* options */
    if ( inpInChI /*&& inpInChI->szOptions*/ ) {  /* fix bug discovered by Burt Leland 2008-12-23 */
        int opt_len = (inpInChI->szOptions? strlen(inpInChI->szOptions) : 0) + sizeof(szMainOption) + 1;
        szOptions = (char*)inchi_calloc( opt_len+1, sizeof(szOptions[0]) );
        if ( szOptions ) {
            if ( inpInChI->szOptions ) /* fix bug discovered by Burt Leland 2008-12-23 */
                strcpy( szOptions, inpInChI->szOptions );
            strcat( szOptions, szMainOption );
            argc = parse_options_string ( szOptions, argv, INCHI_MAX_NUM_ARG );
        } else {
            nRet = _IS_FATAL;
            goto translate_RetVal; /* emergency exit */
        }
    } else {
        argc = 1;
        argv[0] = "";
        argv[1] = NULL;
    }

    if ( (argc == 1
#ifdef TARGET_API_LIB
        && (!inpInChI || !inpInChI->szInChI)
#endif   
        )
        || (argc==2 && ( argv[1][0]==INCHI_OPTION_PREFX ) &&
        (!strcmp(argv[1]+1, "?") || !stricmp(argv[1]+1, "help") )) ) {
        HelpCommandLineParms(log_file);
        outStruct->szLog = log_file->s.pStr;
        nRet = _IS_EOF;
        goto translate_RetVal;
    }

    nRet1 = ReadCommandLineParms( argc, argv, ip, szSdfDataValue, &ulDisplTime, bReleaseVersion, log_file );
    if ( szOptions ) {
        /* argv pointed to strings in szOptions */
        inchi_free( szOptions );
        szOptions = NULL;
    }
    /* INChI DLL specific */
    ip->bNoStructLabels = 1;

    if ( 0 > nRet1 ) {
        goto exit_function;
    }
    if ( ip->bNoStructLabels ) {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    } else
    if ( ip->nInputType == INPUT_INCHI_XML || ip->nInputType == INPUT_INCHI_PLAIN  ||
         ip->nInputType == INPUT_CMLFILE || ip->nInputType == INPUT_INCHI ) {
        /* the input may contain both the header and the label of the structure */
        if ( !ip->pSdfLabel ) 
            ip->pSdfLabel  = ip->szSdfDataHeader;
        if ( !ip->pSdfValue )
            ip->pSdfValue  = szSdfDataValue;
    }
    if ( ip->nInputType && ip->nInputType != INPUT_INCHI ) {
        inchi_ios_eprint( log_file, "Input type set to INPUT_INCHI\n" );
        ip->nInputType = INPUT_INCHI;
    }

    if ( !inpInChI->szInChI ) 
    {
        nRet = _IS_ERROR;
        goto exit_function;
    }
    else
    {
        const int strict=0; /* do not use strict mode, it may be too alarmous */
        nRet = CheckINCHI(inpInChI->szInChI, strict);
        if (nRet == INCHI_VALID_STANDARD)
        {
            bStdFormat = 1;
        }
        else if (nRet == INCHI_VALID_NON_STANDARD)
        {
            ;
        }
        else
        {
            nRet = _IS_ERROR;
            goto exit_function;
        }
    }


    if ( bLibInchiSemaphore ) {  /* does not work properly under sufficient stress */
        return inchi_Ret_BUSY;
    }
    bLibInchiSemaphore = 1;


    PrintInputParms( log_file, ip );
    /*********************************/
    /* InChI -> Structure conversion */
    /*********************************/

    /* input_file simulation */
    input_file->s.pStr = inpInChI->szInChI;
    input_file->s.nUsedLength = strlen(input_file->s.pStr)+1;
    input_file->s.nAllocatedLength = input_file->s.nUsedLength;
    input_file->s.nPtr = 0;
    /* buffer for the message */
    outStruct->szMessage = (char *)inchi_calloc( MAX_MSG_LEN, sizeof(outStruct->szMessage[0]));
    if ( !outStruct->szMessage ) {
         inchi_ios_eprint( log_file, "Cannot allocate output message buffer.\n");
        nRet = -1;
    } else {
        nRet = ReadWriteInChI( input_file, output_file, log_file,
                               ip,  sd, &at,  &num_at, outStruct->szMessage, MAX_MSG_LEN, outStruct->WarningFlags );
    }
    if ( nRet >= 0 && at && num_at ) {
        /* success */
        nRet = InpAtom0DToInchiAtom( at, num_at, outStruct );
        if ( at ) {
            inchi_free( at );
            at = NULL;
        }
        if ( nRet < 0 ) {
            inchi_ios_eprint( log_file, "Final structure conversion failed\n" );
        }
    }
    outStruct->szLog = log_file->s.pStr;


exit_function:;

#if( ADD_CMLPP == 1 )
        /* BILLY 8/6/04 */
        /* free CML memory */
        FreeCml ();
        FreeCmlDoc( 1 );
#endif



    for ( i = 0; i < MAX_NUM_PATHS; i ++ ) {
        if ( ip->path[i] ) {
            inchi_free( (char*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( );


#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if ( num_repeat-- > 0 ) {
        goto repeat;
    }
#endif


#ifdef TARGET_API_LIB
    /* output */

    if ( log_file->s.pStr && log_file->s.nUsedLength > 0 ) {
        while ( log_file->s.nUsedLength && '\n' == log_file->s.pStr[log_file->s.nUsedLength-1] ) {
            log_file->s.pStr[-- log_file->s.nUsedLength]  = '\0'; /* remove last LF */
        }
        if ( outStruct ) {
            outStruct->szLog = log_file->s.pStr;
            log_file->s.pStr = NULL;
        }
    }
#endif
    
translate_RetVal:

    /* Close internal I/O streams */
    inchi_ios_reset(input_file);  /* do not close input_file - its string buffer may point to inpInChI->szInChI */
    inchi_ios_close(output_file);
    inchi_ios_close(log_file);

    switch (nRet) { 
    case -3         : nRet = inchi_Ret_ERROR  ; break; /* Error: no Structure has been created */
    case -2         : nRet = inchi_Ret_ERROR  ; break; /* Error: no Structure has been created */
    case -1         : nRet = inchi_Ret_FATAL  ; break; /* Severe error: no Structure has been created (typically; break; memory allocation failed) */
    default         :
        if ( !outStruct->atom || !outStruct->num_atoms ) {
            nRet = inchi_Ret_EOF;
        } else {
            int m,n,t=0;
            for ( m=0; m < 2; m ++ ) {
                for ( n=0; n < 2; n ++ ) {
                    if ( outStruct->WarningFlags[m][n] ) {
                        t ++;
                    }
                }
            }
            nRet = t? inchi_Ret_WARNING : inchi_Ret_OKAY;
        }
            break;
    }

    bLibInchiSemaphore = 0;
    return nRet;
}

/********************************************************************/

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
                                       int bDoNotAddH, int bDiffUnkUndfStereo,
                                       InchiInpData *pInchiInp );
int  cdecl_Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, 
                                           int bDoNotAddH,
                                           InchiInpData *pInchiInp );
void cdecl_Free_inchi_Input( inchi_Input *pInp );
void cdecl_Free_std_inchi_Input( inchi_Input *pInp );
int cdecl_GetStructFromINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct );
int cdecl_GetStructFromStdINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct );
int cdecl_GetINCHIfromINCHI( inchi_InputINCHI *inpInChI, inchi_Output *out );
void cdecl_FreeStructFromINCHI( inchi_OutputStruct *outStruct );
void cdecl_FreeStructFromStdINCHI( inchi_OutputStruct *outStruct );
int cdecl_CheckINCHI(const char *szINCHI, const int strict);
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

/* implementation */
/* libinchi.def provides export without cdecl_ prefixes */

/********************************************************/
int cdecl_GetINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetINCHI( inp, out );
}
/********************************************************/
int cdecl_GetStdINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetStdINCHI( inp, out );
}
/********************************************************/
void cdecl_FreeINCHI( inchi_Output *out )
{
    FreeINCHI( out );
}
/********************************************************/
void cdecl_FreeStdINCHI( inchi_Output *out )
{
    FreeStdINCHI( out );
}
/********************************************************/
int cdecl_GetStringLength( char *p )
{
    return GetStringLength( p );
}
/********************************************************/
int cdecl_Get_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, 
                                      int bDoNotAddH, int bDiffUnkUndfStereo,
                                      InchiInpData *pInchiInp )
{
    return Get_inchi_Input_FromAuxInfo( szInchiAuxInfo, bDoNotAddH, bDiffUnkUndfStereo,
                                        pInchiInp );
}
/********************************************************/
/********************************************************/
int cdecl_Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, 
                                          int bDoNotAddH,
                                          InchiInpData *pInchiInp )
{
    return Get_std_inchi_Input_FromAuxInfo( szInchiAuxInfo, bDoNotAddH, pInchiInp );
}
/********************************************************/
void cdecl_Free_std_inchi_Input( inchi_Input *pInp )
{
    Free_std_inchi_Input( pInp );
}
/********************************************************/
void cdecl_Free_inchi_Input( inchi_Input *pInp )
{
    Free_inchi_Input( pInp );
}
/********************************************************/
int cdecl_GetStructFromINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct )
{
    return GetStructFromINCHI( inpInChI, outStruct );
}
/********************************************************//********************************************************/
int cdecl_GetStructFromStdINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct )
{
    return GetStructFromStdINCHI( inpInChI, outStruct );
}
/********************************************************/
void cdecl_FreeStructFromINCHI( inchi_OutputStruct *outStruct )
{
    FreeStructFromINCHI( outStruct );
}
/********************************************************/
int cdecl_GetINCHIfromINCHI( inchi_InputINCHI *inpInChI, inchi_Output *out )
{
    return GetINCHIfromINCHI( inpInChI, out );
}
/********************************************************/
void cdecl_FreeStructFromStdINCHI( inchi_OutputStruct *outStruct )
{
    FreeStructFromStdINCHI( outStruct );
}
/********************************************************/
int cdecl_CheckINCHI(const char *szINCHI, const int strict)
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
int PASCAL pasc_CheckINCHI(const char *szINCHI, const int strict);

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

/* implementation */
/* libinchi.def provides export without PASCAL pasc_ prefixes */
/********************************************************/
int PASCAL pasc_GetINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetINCHI( inp, out );
}
/********************************************************/
int PASCAL pasc_GetStdINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetStdINCHI( inp, out );
}
/********************************************************/
void PASCAL pasc_FreeINCHI( inchi_Output *out )
{
    FreeINCHI( out );
}
/********************************************************/
void PASCAL pasc_FreeStdINCHI( inchi_Output *out )
{
    FreeStdINCHI( out );
}
/********************************************************/
int PASCAL pasc_GetStringLength( char *p )
{
    return GetStringLength( p );
}
/********************************************************/
int PASCAL pasc_Get_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, 
                                                int bDoNotAddH, 
                                                int bDiffUnkUndfStereo,
                                                InchiInpData *pInchiInp )
{
    return Get_inchi_Input_FromAuxInfo( szInchiAuxInfo, bDoNotAddH, 
                                            bDiffUnkUndfStereo, pInchiInp );
}
/********************************************************/
int PASCAL pasc_Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, 
                                                int bDoNotAddH, 
                                                InchiInpData *pInchiInp )
{
    return Get_std_inchi_Input_FromAuxInfo( szInchiAuxInfo, bDoNotAddH, pInchiInp );
}
/********************************************************/
void PASCAL pasc_Free_inchi_Input( inchi_Input *pInp )
{
    Free_inchi_Input( pInp );
}
/********************************************************/
void PASCAL pasc_Free_std_inchi_Input( inchi_Input *pInp )
{
    Free_std_inchi_Input( pInp );
}
/********************************************************/
void PASCAL pasc_FreeStructFromINCHI( inchi_OutputStruct *out )
{
    FreeStructFromINCHI( out );
}
/********************************************************/
void PASCAL pasc_FreeStructFromStdINCHI( inchi_OutputStruct *out )
{
    FreeStructFromStdINCHI( out );
}
/********************************************************//********************************************************/
int PASCAL pasc_GetStructFromINCHI( inchi_InputINCHI *inp, inchi_OutputStruct *out )
{
    return GetStructFromINCHI( inp, out );
}
/********************************************************//********************************************************/
int PASCAL pasc_GetStructFromStdINCHI( inchi_InputINCHI *inp, inchi_OutputStruct *out )
{
    return GetStructFromStdINCHI( inp, out );
}
/********************************************************/
int PASCAL pasc_CheckINCHI(const char *szINCHI, const int strict)
{
    return CheckINCHI( szINCHI, strict );
}

#endif 


