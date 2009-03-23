/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02
 * November 30, 2008
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


#include "mode.h"

#ifndef INCHI_ANSI_ONLY
#ifndef INCHI_LIB
#include <windows.h>
#endif
#endif

#include <stdio.h>
#include <stdlib.h>

#ifndef INCHI_ANSI_ONLY
#include <conio.h>
#endif

#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include <float.h>

#include "ichitime.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichidrp.h"
#include "ichierr.h"
#include "ichimain.h"

#include "ichicomp.h"

#include "ichi_io.h"

#if( ADD_CMLPP == 1 )
#include "readcml.hpp"
#endif

/*^^^ */
#ifdef BUILD_CINCHI_WITH_INCHIKEY
#include "inchi_api.h"
void extract_inchi_substring(char ** buf, char *str, size_t slen);
#endif
/*^^^ */

/*  console-specific */
#ifndef INCHI_ANSI_ONLY
/********************************************************************/
int user_quit( const char *msg, unsigned long ulMaxTime )
{
#ifdef INCHI_LIB
    return 0;
#endif
#if( defined(_WIN32) && !defined(INCHI_LIB) )
    int quit, enter, ret;
    printf(msg);
    if ( ulMaxTime ) {
        inchiTime  ulEndTime;
        InchiTimeGet( &ulEndTime );
        InchiTimeAddMsec( &ulEndTime, ulMaxTime );
        while ( !_kbhit() ) {
            if ( bInchiTimeIsOver( &ulEndTime ) ) {
                printf("\n");
                return 0;
            }
            MySleep( 100 );
        }
    }
    while( 1 ) {
        quit  = ( 'q' == (ret = _getch()) || 'Q'==ret || /*Esc*/ 27==ret );
        enter = ( '\r' == ret );
        if ( ret == 0xE0 )
            ret = _getch();
        else {
            _putch(ret); /* echo */
        }
        if ( quit || enter )
            break;
        printf( "\r" );
        printf( msg );
    }

    _putch('\n');
    return quit;
#else
    return 0;
#endif
}
/*****************************************************************/
void eat_keyboard_input( void )
{
#ifndef INCHI_LIB
    while ( _kbhit() ) {
        if ( 0xE0 == _getch() )
            _getch();
    }
#endif
}
#endif /* ifndef INCHI_ANSI_ONLY */

#ifdef INCHI_ANSI_ONLY
/*****************************************************************/
int user_quit( const char *msg, unsigned long ulMaxTime )
{
    return 0;
}
/*****************************************************************/
void eat_keyboard_input( void )
{
}
#endif


/*****************************************************************/
#ifndef INCHI_LIB

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
int WasInterrupted(void) {
#ifdef _DEBUG            
    if ( bInterrupted ) {
        int stop=1;  /*  for debug only <BRKPT> */
    }
#endif
    return bInterrupted;
}
#endif

#endif

#define REPEAT_ALL  0
/*#define TEST_FPTRS*/ /* uncomment for INCHI_LIB testing only */
/********************************************************************/
int main( int argc, char *argv[ ] )
{

    STRUCT_DATA struct_data;
    STRUCT_DATA *sd = &struct_data;

    INCHI_IOSTREAM outputstr, logstr, prbstr, instr;
    INCHI_IOSTREAM *pout=&outputstr, *plog = &logstr, *pprb = &prbstr, *inp_file = &instr;

    char szTitle[MAX_SDF_HEADER+MAX_SDF_VALUE+256];

    /* long rcPict[4] = {0,0,0,0}; */
    int   i;
    long  num_err, num_output;
    long num_inp;
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
    const int nStrLen = INCHI_SEGM_BUFLEN;
    char *pStr = NULL;
    int   nRet = 0, nRet1;

#ifndef TEST_FPTRS
    STRUCT_FPTRS *pStructPtrs = NULL;
#else
    STRUCT_FPTRS struct_fptrs, *pStructPtrs =&struct_fptrs; /* INCHI_LIB debug only */
#endif
#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
    int  num_repeat = REPEAT_ALL;
#endif

/*^^^ */
    int inchi_ios_type;
#ifdef BUILD_CINCHI_WITH_INCHIKEY
    char ik_string[256];    /*^^^ Resulting InChIKey string */
    int ik_ret=0;           /*^^^ InChIKey-calc result code */
    inchi_ios_type = INCHI_IOSTREAM_STRING;
#else
    inchi_ios_type = INCHI_IOSTREAM_FILE;
#endif



#if( TRACE_MEMORY_LEAKS == 1 )
    _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_ALLOC_MEM_DF);
/* for execution outside the VC++ debugger uncomment one of the following two */
/*#define MY_REPORT_FILE  _CRTDBG_FILE_STDERR */
/*#define MY_REPORT_FILE  _CRTDBG_FILE_STDOUT */
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


#if ( defined(REPEAT_ALL) && REPEAT_ALL > 0 )
repeat:
    inchi_ios_close(inp_file);
    inchi_ios_close(pout);
    inchi_ios_close(plog);
    inchi_ios_close(pprb);
    pStr = NULL;
#endif

    num_inp    = 0;
    num_err    = 0;
    num_output = 0;
    sd->bUserQuit  = 0;
    
#if( defined( _WIN32 ) && defined( _CONSOLE ) && !defined( INCHI_ANSI_ONLY ) )
    if ( SetConsoleCtrlHandler( MyHandlerRoutine, 1 ) ) {
        ConsoleQuit = WasInterrupted;
    }
#endif

    inchi_ios_init(inp_file, INCHI_IOSTREAM_FILE, NULL);
    inchi_ios_init(pout, inchi_ios_type, NULL);
    inchi_ios_init(plog, inchi_ios_type, stdout);
    inchi_ios_init(pprb, inchi_ios_type, NULL);



    if ( argc == 1 || argc==2 && ( argv[1][0]==INCHI_OPTION_PREFX ) &&
        (!strcmp(argv[1]+1, "?") || !stricmp(argv[1]+1, "help") ) ) 
    {
        HelpCommandLineParms(plog); 
        inchi_ios_flush(plog);
        return 0;
    }

    /*  original input structure */
    memset( orig_inp_data     , 0,   sizeof( *orig_inp_data  ) );
    memset( prep_inp_data     , 0, 2*sizeof( *prep_inp_data  ) );
    memset( pINChI,     0, sizeof(pINChI    ) );
    memset( pINChI_Aux, 0, sizeof(pINChI_Aux) );
    memset( szSdfDataValue    , 0, sizeof( szSdfDataValue    ) );

    /* explicitly cast to (const char **) to avoid a warning about "incompatible pointer type":*/
    
    plog->f = stderr;
    if ( 0 > ReadCommandLineParms( argc, (const char **)argv, ip, szSdfDataValue, &ulDisplTime, bReleaseVersion, plog) ) 
        goto exit_function;

    if ( !OpenFiles( &(inp_file->f), &(pout->f), &(plog->f), &(pprb->f), ip ) ) 
        goto exit_function;

    if ( ip->bNoStructLabels ) 
    {
        ip->pSdfLabel = NULL;
        ip->pSdfValue = NULL;
    } else
    if ( ip->nInputType == INPUT_INCHI_XML || ip->nInputType == INPUT_INCHI_PLAIN  ||
         ip->nInputType == INPUT_CMLFILE   || ip->nInputType == INPUT_INCHI  ) {
        /* the input may contain both the header and the label of the structure */
        if ( !ip->pSdfLabel ) 
            ip->pSdfLabel  = ip->szSdfDataHeader;
        if ( !ip->pSdfValue )
            ip->pSdfValue  = szSdfDataValue;
    }

    PrintInputParms(plog,ip);
    inchi_ios_flush2(plog, stderr);

    if ( !(pStr = (char*)inchi_malloc(nStrLen))) 
    {
        inchi_ios_eprint( plog, "Cannot allocate output buffer. Terminating\n"); 
        inchi_ios_flush2(plog, stderr);
        goto exit_function;
    }
    pStr[0] = '\0';

#if( READ_INCHI_STRING == 1 )
    if ( ip->nInputType == INPUT_INCHI ) {
        memset( sd, 0, sizeof(*sd) );
        ReadWriteInChI( inp_file, pout, plog, ip,  sd, NULL, NULL, NULL, 0, NULL);
        inchi_ios_flush2(plog, stderr);
        ulTotalProcessingTime = sd->ulStructTime;
        num_inp               = sd->fPtrStart;
        num_err               = sd->fPtrEnd;
        goto exit_function;
    }
#endif
    /**********************************************************************************************/
    /*  Main cycle */
    /*  read input structures and create their INChI */
    ulTotalProcessingTime = 0;

    if ( pStructPtrs ) {
        memset ( pStructPtrs, 0, sizeof(pStructPtrs[0]) );
        /* debug: set CML reading sequence
        pStructPtrs->fptr = (INCHI_FPTR *)inchi_calloc(16, sizeof(INCHI_FPTR));
        for ( i = 0; i < 15; i ++ )
            pStructPtrs->fptr[i] = 15-i;
        pStructPtrs->cur_fptr = 7;
        pStructPtrs->len_fptr = 16;
        pStructPtrs->max_fptr = 14;
        */
    }

    while ( !sd->bUserQuit && !bInterrupted ) {
        if ( ip->last_struct_number && num_inp >= ip->last_struct_number ) {
            nRet = _IS_EOF; /*  simulate end of file */
            goto exit_function;
        }

        /*  read one structure from input and display optionally it */

        nRet = GetOneStructure( sd, ip, szTitle, inp_file, plog, pout, pprb,
                                orig_inp_data, &num_inp, pStr, nStrLen, pStructPtrs );
        inchi_ios_flush2(plog, stderr);

        if ( pStructPtrs ) {
            pStructPtrs->cur_fptr ++;
        }

        if ( sd->bUserQuit ) {
            break;
        }
        switch ( nRet ) {
        case _IS_FATAL:
            num_err ++;
        case _IS_EOF:
            goto exit_function;
        case _IS_ERROR:
            num_err ++;
        case _IS_SKIP:
            continue;
        }

        /* create INChI for each connected component of the structure and optionally display them */
        /* output INChI for the whole structure */
        
        nRet1 = ProcessOneStructure( sd, ip, szTitle, pINChI, pINChI_Aux,
                                     inp_file, plog, pout, pprb,
                                     orig_inp_data, prep_inp_data,
                                     num_inp, pStr, nStrLen );        
        inchi_ios_flush2(plog, stderr);

#ifdef BUILD_CINCHI_WITH_INCHIKEY
        /*^^^ post-1.02b addition - correctly treat tabbed output with InChIKey */
        if ( ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT ) 
            if ( ip->bCalcInChIKey )
                if (pout->s.pStr)
                    if (pout->s.nUsedLength>0)
                        if (pout->s.pStr[pout->s.nUsedLength-1]=='\n')
                            /* replace LF with TAB */
                            pout->s.pStr[pout->s.nUsedLength-1] = '\t';


#endif

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
            continue;
        }

#ifdef BUILD_CINCHI_WITH_INCHIKEY
        if (ip->bCalcInChIKey)
        {
            /*^^^ post-1.02b */
            char *buf = NULL;
            size_t slen = pout->s.nUsedLength;

            extract_inchi_substring(&buf, pout->s.pStr, slen);
            
            if (NULL!=buf)
            {
#ifdef ONLY_STDINCHI_STDKEY
                ik_ret = GetStdINCHIKeyFromStdINCHI(buf, ik_string);                     
#else
                /* to be filled */
                ik_ret = INCHIKEY_UNKNOWN_ERROR;
#endif
                inchi_free(buf);
            }
            else
                ik_ret = INCHIKEY_NOT_INCHI_INPUT;                     
        

            if (ik_ret==INCHIKEY_OK)   
            {
                inchi_ios_print(pout, "InChIKey=%-s\n",ik_string);
            }
            else    
            {
                inchi_ios_print(plog, "Warning (Could not compute InChIKey: ", num_inp);
                switch(ik_ret)
                {
                case INCHIKEY_UNKNOWN_ERROR:
                        inchi_ios_print(plog, "unresolved error)");
                        break;
                case INCHIKEY_EMPTY_INPUT:
                        inchi_ios_print(plog,  "got an empty string)");
                        break;
                case INCHIKEY_NOT_INCHI_INPUT:
                        inchi_ios_print(plog, "got non-InChI string)");
                        break;
                case INCHIKEY_NOT_ENOUGH_MEMORY:
                        inchi_ios_print(plog, "not enough memory to treat the string)");
                        break;
                case INCHIKEY_ERROR_IN_FLAG_CHAR:
                        inchi_ios_print(plog, "detected error in flag character)");
                        break;
                default:inchi_ios_print(plog, "internal program error)");
                        break;
                }
                inchi_ios_print(plog, " structure #%-lu.\n", num_inp);
                if ( ip->bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT )
                    inchi_ios_print(pout, "\n");
            }
            inchi_ios_flush(pout);
            inchi_ios_flush2(plog, stderr);
        }

        else
            inchi_ios_flush(pout);


#endif

   /*   --- debug only ---
        if ( pStructPtrs->cur_fptr > 5 ) {
            pStructPtrs->cur_fptr = 5;
        }
   */


            
    }

exit_function:
    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && sd->bXmlStructStarted > 0 ) 
    {
            if ( !OutputINChIXmlStructEndTag( pout,  pStr, nStrLen, 1 ) ) 
            {
                inchi_ios_eprint( plog, "Cannot create end xml tag for structure #%ld.%s%s%s%s Terminating.\n", num_inp, SDF_LBL_VAL(ip->pSdfLabel,ip->pSdfValue) );
                inchi_ios_flush2(plog, stderr);
                sd->bXmlStructStarted = -1; /*  do not repeat same message */
        }
    }


    if ( (ip->bINChIOutputOptions & INCHI_OUT_XML) && ip->bXmlStarted ) 
    {
        OutputINChIXmlRootEndTag( pout ); 
        inchi_ios_flush(pout);
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


    /*^^^ close output(s) */
    inchi_ios_close(inp_file);    
    inchi_ios_close(pout);    
    inchi_ios_close(pprb);
    {
        int hours, minutes, seconds, mseconds;
        SplitTime( ulTotalProcessingTime, &hours, &minutes, &seconds, &mseconds );
        inchi_ios_eprint( plog, "Finished processing %ld structure%s: %ld error%s, processing time %d:%02d:%02d.%02d\n",
                                num_inp, num_inp==1?"":"s",
                                num_err, num_err==1?"":"s",
                                hours, minutes, seconds,mseconds/10);
        inchi_ios_flush2(plog, stderr);
    }
    inchi_ios_close(plog);
    if ( pStr )                             
        inchi_free( pStr );


    /* frees */
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

#if( bRELEASE_VERSION != 1 && defined(_DEBUG) )
    if ( inp_file->f && inp_file->f != stdin ) {
        user_quit("Press Enter to exit ?", ulDisplTime);
    }
#endif

    
    return 0;
}




#endif  /* ifndef INCHI_LIB */


