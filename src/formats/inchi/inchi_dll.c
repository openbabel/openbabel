/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.00
 * April 13, 2005
 * Developed at NIST
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

#include "ichicomp.h"
#include "inchi_api.h"
#include "inchi_dll.h"

/*************************************************************************
 *
 *   Local protopypes
 *
 *************************************************************************/


int SetAtomProperties( inp_ATOM *at, MOL_COORD *szCoord, inchi_Atom *ati,
                       int a1, int *nDim, char *pStrErr, int *err );
int SetBondProperties( inp_ATOM *at, inchi_Atom *ati, int a1, int j,
                       int nNumAtoms, int *nNumBonds, char *pStrErr, int *err );
int SetAtomAndBondProperties( inp_ATOM *at, inchi_Atom *ati, int a1,
                              int bDoNotAddH, char *pStrErr, int *err );
void SetNumImplicitH(inp_ATOM* at, int num_atoms);
int Extract0DParities( inp_ATOM *at, int nNumAtoms, inchi_Stereo0D *stereo0D,
                       int num_stereo0D, char *pStrErr, int *err );
int parse_options_string ( char *cmd, const char *argv[], int maxargs );

/*************************************************************************/

int bInterrupted = 0;

#if( defined( _WIN32 ) && defined( _CONSOLE ) )

#ifndef INCHI_ANSI_ONLY
BOOL WINAPI MyHandlerRoutine(
  DWORD dwCtrlType   /*   control signal type */
  ) {
    if ( dwCtrlType == CTRL_C_EVENT     ||
         dwCtrlType == CTRL_BREAK_EVENT ||
         dwCtrlType == CTRL_CLOSE_EVENT ||
         dwCtrlType == CTRL_LOGOFF_EVENT ) {
        bInterrupted = 1;
        return TRUE;
    }
    return FALSE;
}
#endif
int WasInterrupted(void) {
#ifdef _DEBUG            
    if ( bInterrupted ) {
        int stop=1;  /*  for debug only <BRKPT> */
    }
#endif
    return bInterrupted;
}

#endif


/********************************************************************
 *
 * INCHI API: DEALLOCATE INCHI OUTPUT
 *
 ********************************************************************/
void INCHI_DECL FreeINCHI( inchi_Output *out )
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
/********************************************************************/
#define INCHI_MAX_NUM_ARG 32
/********************************************************************
 *
 *    INCHI API: MAIN ENTRY POINT
 *
 ********************************************************************/

int bLibInchiSemaphore = 0;

int INCHI_DECL GetINCHI( inchi_Input *inp, inchi_Output *out )
{

    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;
    FILE *inp_file = NULL;
    INCHI_FILE inchi_file[3];
    INCHI_FILE *output_file = inchi_file, *log_file = inchi_file+1, *prb_file = inchi_file+2;
    char szTitle[MAX_SDF_HEADER+MAX_SDF_VALUE+256];

    int i, num_inp, num_err, num_output;
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


#if( defined( _WIN32 ) && defined( _CONSOLE ) && !defined( INCHI_ANSI_ONLY ) )
    if ( SetConsoleCtrlHandler( MyHandlerRoutine, 1 ) ) {
        ConsoleQuit = WasInterrupted;
    }
#endif
    
    memset( inchi_file, 0, sizeof(inchi_file) );

#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
repeat:
    inp_file = output_file = log_file = prb_file = NULL;
    pStr = NULL;
#endif

    num_inp    = 0;
    num_err    = 0;
    num_output = 0;
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

    if ( argc == 1
#ifdef INCHI_LIBRARY
        && (!inp || inp->num_atoms <= 0 || !inp->atom)
#endif        
        || argc==2 && ( argv[1][0]==INCHI_OPTION_PREFX ) &&
        (!strcmp(argv[1]+1, "?") || !stricmp(argv[1]+1, "help") ) ) {
        HelpCommandLineParms(log_file);
        out->szLog = log_file->pStr;
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
        goto exit_function;
    }
#ifndef INCHI_LIBRARY
    if ( !OpenFiles( &inp_file, &output_file, &log_file, &prb_file, ip ) ) {
        goto exit_function;
    }
#endif
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
    PrintInputParms( log_file, ip );
    if ( !(pStr = (char*)inchi_malloc(nStrLen))) {
        my_fprintf( log_file, "Cannot allocate output buffer. Terminating\n");
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

    /* === TOMORROW: remove while cycling and convert inp to orig_inp_data ==== */
#ifdef INCHI_LIBRARY
    if ( !sd->bUserQuit && !bInterrupted )
#else
    while ( !sd->bUserQuit && !bInterrupted )
#endif
    {
        if ( ip->last_struct_number && num_inp >= ip->last_struct_number ) {
            nRet = _IS_EOF; /*  simulate end of file */
            goto exit_function;
        }

#ifndef INCHI_LIBRARY
        /*  read one structure from input and display optionally it */
        nRet = GetOneStructure( sd, ip, szTitle, inp_file, log_file, output_file, prb_file,
                                orig_inp_data, &num_inp, pStr, nStrLen, pStructPtrs );
#else
        nRet = ExtractOneStructure( sd, ip, szTitle, inp, log_file, output_file, prb_file,
                                orig_inp_data, &num_inp, pStr, nStrLen );
#endif


        if ( pStructPtrs ) {
            pStructPtrs->cur_fptr ++;
        }

#ifndef INCHI_LIBRARY
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
#ifndef INCHI_LIBRARY
        case _IS_SKIP:
            continue;
#endif
        }

        /* create INChI for each connected component of the structure and optionally display them */
        /* output INChI for the whole structure */
        nRet1 = ProcessOneStructure( sd, ip, szTitle, pINChI, pINChI_Aux,
                                     inp_file, log_file, output_file, prb_file,
                                     orig_inp_data, prep_inp_data,
                                     num_inp, pStr, nStrLen );
        
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
            num_err ++;
            goto exit_function;
        case _IS_ERROR:
            num_err ++;
#ifndef INCHI_LIBRARY
            continue;
#endif
        }

    }

exit_function:
    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && sd->bXmlStructStarted > 0 ) {
        if ( !OutputINChIXmlStructEndTag( output_file, pStr, nStrLen, 1 ) ) {
            my_fprintf( log_file, "Cannot create end xml tag for structure #%d.%s%s%s%s Terminating.\n", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
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

#ifndef INCHI_LIBRARY
    if ( inp_file && inp_file != stdin) {
        fclose ( inp_file );
    }
    if ( prb_file ) {
        fclose ( prb_file );
    }
    if ( output_file && output_file != stdout ) {
        fclose( output_file );
    }
    
    if ( log_file )
    {
        int hours, minutes, seconds, mseconds;
        SplitTime( ulTotalProcessingTime, &hours, &minutes, &seconds, &mseconds );
        my_fprintf( log_file, "Finished processing %d structure%s: %d error%s, processing time %d:%02d:%02d.%02d\n",
                                num_inp, num_inp==1?"":"s",
                                num_err, num_err==1?"":"s",
                                hours, minutes, seconds,mseconds/10);
    }

    if ( log_file && log_file != stderr ) {
        fclose( log_file );
    }
#endif

    if ( pStr ) {
        inchi_free( pStr );
    }

    for ( i = 0; i < MAX_NUM_PATHS; i ++ ) {
        if ( ip->path[i] ) {
            inchi_free( (void*) ip->path[i] ); /*  cast deliberately discards 'const' qualifier */
            ip->path[i] = NULL;
        }
    }

    SetBitFree( );


#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    if ( num_repeat-- > 0 ) {
        goto repeat;
    }
#endif

#ifndef INCHI_LIBRARY
#if( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    if ( inp_file && inp_file != stdin ) {
        user_quit("Press Enter to exit ?", ulDisplTime);
    }
#endif
#endif


#ifdef INCHI_LIBRARY
    /* output */
    if ( sd->pStrErrStruct[0] ) {
        if ( out && (out->szMessage = (char *)inchi_malloc( strlen(sd->pStrErrStruct) + 1 )) ) {
            strcpy( out->szMessage, sd->pStrErrStruct );
        }
    }
    if ( output_file->pStr && output_file->nUsedLength > 0 && out ) {
        char *p;
        out->szInChI   = output_file->pStr;
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
        output_file->pStr = NULL;
    }
    if ( log_file->pStr && log_file->nUsedLength > 0 ) {
        while ( log_file->nUsedLength && '\n' == log_file->pStr[log_file->nUsedLength-1] ) {
            log_file->pStr[-- log_file->nUsedLength]  = '\0'; /* remove last LF */
        }
        if ( out ) {
            out->szLog = log_file->pStr;
            log_file->pStr = NULL;
        }
    }
    if ( output_file->pStr )
        inchi_free( output_file->pStr );
    if ( log_file->pStr    )
        inchi_free( log_file->pStr );
    

#endif
    
translate_RetVal:

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


/*************************************************************************/
/******************************** from readmol.c *************************/
/*************************************************************************/
int AddMOLfileError( char *pStrErr, const char *szMsg )
{
    if ( pStrErr && szMsg && szMsg[0] ) {
        int lenStrErr = strlen( pStrErr );
        int lenMsg    = strlen( szMsg );
        char *p = strstr( pStrErr, szMsg );
        if ( p && (p==pStrErr || *(p-1) == ' ' && (*(p-2) == ';' || *(p-2) == ':' )) &&
                  (p+lenMsg == pStrErr+lenStrErr || 
                  p[lenMsg] == ';' && p[lenMsg+1] == ' ' ||
                  p[lenMsg-1]==':' && p[lenMsg]==' ') ) {
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
int CopyMOLfile(FILE *inp_file, long fPtrStart, long fPtrEnd, INCHI_FILE *prb_file, long lNumb)
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
    FreeInpAtom( &inp_at_data->at );
    FreeInpAtom( &inp_at_data->at_fixed_bonds );
    memset( inp_at_data, 0, sizeof(*inp_at_data) );
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
            swap( (char*)&dest_at[m].neighbor[j],    (char*)&dest_at[m].neighbor[k],    sizeof(dest_at[0].neighbor[0]) );
            swap( (char*)&dest_at[m].bond_stereo[j], (char*)&dest_at[m].bond_stereo[k], sizeof(dest_at[0].bond_stereo[0]) );
            swap( (char*)&dest_at[m].bond_type[j],   (char*)&dest_at[m].bond_type[k],   sizeof(dest_at[0].bond_type[0]) );
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
        swap( (char*)&new_ord[j], (char*)&new_ord[k], sizeof(new_ord[0]) );
    }
}
#endif /*  } TEST_RENUMB_ATOMS == 1  */
/**********************************************************************************/
void FreeOrigAtData( ORIG_ATOM_DATA *orig_at_data )
{
    if ( !orig_at_data )
        return;
    FreeInpAtom( &orig_at_data->at );
    if ( orig_at_data->nCurAtLen ) {
        inchi_free( orig_at_data->nCurAtLen );
    }
    if ( orig_at_data->nOldCompNumber ) {
        inchi_free( orig_at_data->nOldCompNumber );
    }
    if ( orig_at_data->szCoord ) {
        inchi_free( orig_at_data->szCoord );
    }
    if ( orig_at_data->nEquLabels ) {
        inchi_free( orig_at_data->nEquLabels );
    }
    if ( orig_at_data->nSortedOrder ) {
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
        sprintf( szRadicalType, "%d", ati[a1].radical );
        MOLFILE_ERR_SET (*err, 0, "Radical center type ignored:");
        MOLFILE_ERR_SET (*err, 0, szRadicalType);
        *err |= 8; /*  Unrecognized Radical replaced with non-radical */
        cRadical = 0;
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
        if ( n1+1 < at[a1].valence && is_in_the_list( at[a1].neighbor+n1+1, (AT_NUMB)a2, at[a1].valence-n1-1 ) ||
             n2+1 < at[a2].valence && is_in_the_list( at[a2].neighbor+n2+1, (AT_NUMB)a1, at[a2].valence-n2-1 ) ) {
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
        if ( p1 && (cBondType != at[a1].bond_type[n1] || at[a1].bond_stereo[n1] != cStereoType1 )||
             p2 && (cBondType != at[a2].bond_type[n2] || at[a2].bond_stereo[n2] != cStereoType2 ) ) {
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
            if ( nRadical && at[a1].radical && nRadical != at[a1].radical ||
                 nCharge  && at[a1].charge  && nCharge  != at[a1].charge ) {
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
        if ( bUnusualValenceArom && !bUnusualValenceNoArom ) {
            /* typically NH in 5-member aromatic ring */
            at[a1].chem_bonds_valence --;
        }
    }

    return 0;
}
/****************************************************************************************/
int ExtractOneStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
         inchi_Input *inp, INCHI_FILE *log_file, INCHI_FILE *output_file,
         INCHI_FILE *prb_file, ORIG_ATOM_DATA *orig_inp_data, int *num_inp,
         char *pStr, int nStrLen )
{
    int         *err           = &sd->nStructReadError;
    char        *pStrErr       = sd->pStrErrStruct;
    inp_ATOM    *at            = NULL;
    MOL_COORD   *szCoord       = NULL; 
    inchi_Atom  *ati           = NULL;  
    int       nNumAtoms = 0;
    int       a1, j, valence, nDim, nNumBonds, nRet = 0;

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
    Extract0DParities( at, nNumAtoms, inp->stereo0D, inp->num_stereo0D, pStrErr, err );
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
#if( defined( _WIN32 ) && defined( _MSC_VER ) && _MSC_VER >= 800 && defined(_USRDLL) && defined(INCHI_LINK_AS_DLL) )
    /* Win32 & MS VC ++, compile and link as a DLL */
/*********************************************************/
/*   C calling conventions export from Win32 dll         */
/*********************************************************/
/* prototypes */
int  cdecl_GetINCHI( inchi_Input *inp, inchi_Output *out );
void cdecl_FreeINCHI( inchi_Output *out );
int  cdecl_GetStringLength( char *p );
int  cdecl_Get_inchi_Input_FromAuxInfo
             ( char *szInchiAuxInfo, int bDoNotAddH, InchiInpData *pInchiInp );
void cdecl_Free_inchi_Input( inchi_Input *pInp );

/* implementation */
/* vc6_libinchi.def provides export withou cdecl_ prefixes */
/********************************************************/
int cdecl_GetINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetINCHI( inp, out );
}
/********************************************************/
void cdecl_FreeINCHI( inchi_Output *out )
{
    FreeINCHI( out );
}
/********************************************************/
int cdecl_GetStringLength( char *p )
{
    return GetStringLength( p );
}
/********************************************************/
int cdecl_Get_inchi_Input_FromAuxInfo
             ( char *szInchiAuxInfo, int bDoNotAddH, InchiInpData *pInchiInp )
{
    return Get_inchi_Input_FromAuxInfo
             ( szInchiAuxInfo, bDoNotAddH, pInchiInp );
}
/********************************************************/
void cdecl_Free_inchi_Input( inchi_Input *pInp )
{
    Free_inchi_Input( pInp );
}
#endif

#if( defined(__GNUC__) && __GNUC__ >= 3 && defined(__MINGW32__) && defined(_WIN32) )
#include <windows.h>
/*********************************************************/
/*   Pacal calling conventions export from Win32 dll     */
/*********************************************************/
/* prototypes */
int  PASCAL pasc_GetINCHI( inchi_Input *inp, inchi_Output *out );
void PASCAL pasc_FreeINCHI( inchi_Output *out );
int  PASCAL pasc_GetStringLength( char *p );
int  PASCAL pasc_Get_inchi_Input_FromAuxInfo
             ( char *szInchiAuxInfo, int bDoNotAddH, InchiInpData *pInchiInp );
void PASCAL pasc_Free_inchi_Input( inchi_Input *pInp );

/* implementation */
/* vc6_libinchi.def provides export withou PASCAL pasc_ prefixes */
/********************************************************/
int PASCAL pasc_GetINCHI( inchi_Input *inp, inchi_Output *out )
{
    return GetINCHI( inp, out );
}
/********************************************************/
void PASCAL pasc_FreeINCHI( inchi_Output *out )
{
    FreeINCHI( out );
}
/********************************************************/
int PASCAL pasc_GetStringLength( char *p )
{
    return GetStringLength( p );
}
/********************************************************/
int PASCAL pasc_Get_inchi_Input_FromAuxInfo
             ( char *szInchiAuxInfo, int bDoNotAddH, InchiInpData *pInchiInp )
{
    return Get_inchi_Input_FromAuxInfo
             ( szInchiAuxInfo, bDoNotAddH, pInchiInp );
}
/********************************************************/
void PASCAL pasc_Free_inchi_Input( inchi_Input *pInp )
{
    Free_inchi_Input( pInp );
}
#endif 


