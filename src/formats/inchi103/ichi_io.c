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
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include "mode.h"
#include "ichi_io.h"
#include "ichicomp.h"
#include "util.h"

#ifndef INCHI_ADD_STR_LEN
#define INCHI_ADD_STR_LEN   32768
#endif


#ifdef INCHI_LIB
extern void (*FWPRINT) (const char * format, va_list argptr );
#endif


/*^^^ Internal functions */

int inchi_ios_str_getc( INCHI_IOSTREAM *ios );
char *inchi_ios_str_gets( char *szLine, int len, INCHI_IOSTREAM *ios );
char *inchi_ios_str_getsTab( char *szLine, int len, INCHI_IOSTREAM *ios );
int GetMaxPrintfLength( const char *lpszFormat, va_list argList);
char *inchi_fgetsTab( char *szLine, int len, FILE *f ); 
int inchi_vfprintf( FILE* f, const char* lpszFormat, va_list argList ); 

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


                                    INCHI_IOSTREAM OPERATIONS 


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Init INCHI_IOSTREAM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void inchi_ios_init(INCHI_IOSTREAM* ios, int io_type, FILE *f)
{
    memset( ios, 0, sizeof(*ios) );
    switch (io_type)
    {
        case INCHI_IOSTREAM_FILE:	ios->type = INCHI_IOSTREAM_FILE; 
                                    break;
        case INCHI_IOSTREAM_STRING:
        default:                    ios->type = INCHI_IOSTREAM_STRING;
                                    break;
    }
    ios->f = f; 
    return;
}



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    If INCHI_IOSTREAM type is INCHI_IOSTREAM_STRING,
        flush INCHI_IOSTREAM string buffer to file (if non-NULL); then free buffer.
    If INCHI_IOSTREAM type is INCHI_IOSTREAM_FILE, just flush the file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void inchi_ios_flush(INCHI_IOSTREAM* ios)
{
    
    if (ios->type == INCHI_IOSTREAM_STRING) 
    {
        if (ios->s.pStr)
        {
            if (ios->s.nUsedLength > 0)
            {       
                if (ios->f)
                {
                    fprintf(ios->f,"%-s", ios->s.pStr); 
                    fflush(ios->f);
                }
                inchi_free(ios->s.pStr );
                ios->s.pStr = NULL; 
                ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;
            }       
        }
    }
    
    else if (ios->type == INCHI_IOSTREAM_FILE)
    {
        /* output to plain file: just flush it. */
        fflush(ios->f);
    }
    
    return;
}



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    If INCHI_IOSTREAM type is INCHI_IOSTREAM_STRING,
        flush INCHI_IOSTREAM string buffer to file (if non-NULL) and
        another file f2 supplied as parameter (typically, it will be stderr); then free buffer.
    If INCHI_IOSTREAM type is INCHI_IOSTREAM_FILE, just flush the both files.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void inchi_ios_flush2(INCHI_IOSTREAM* ios, FILE *f2)
{

    if (ios->type == INCHI_IOSTREAM_STRING) 
    {
        if (ios->s.pStr)
        {
            if (ios->s.nUsedLength > 0)
            {       
                if (ios->f)
                {
                    fprintf(ios->f,"%-s", ios->s.pStr); 
                    fflush(ios->f);
                }
                if (f2!=ios->f)
                    fprintf(f2,"%-s", ios->s.pStr); 
                
                inchi_free(ios->s.pStr );
                ios->s.pStr = NULL; 
                ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;
            }       
        }
    }

    else if (ios->type == INCHI_IOSTREAM_FILE)
    {
        /* output to plain file: just flush it. */
        if ( (ios->f) && (ios->f!=stderr) && (ios->f!=stdout) )
            fflush(ios->f);
        if ( f2 && f2!=stderr && f2!=stdout)
            fflush(f2);


    }

    return;
}



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Close INCHI_IOSTREAM: free string buffer and close the file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void inchi_ios_close(INCHI_IOSTREAM* ios)
{
    if (ios->s.pStr)
        inchi_free(ios->s.pStr);
    ios->s.pStr = NULL; 
    ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;
    if ( (ios->f) && (ios->f!=stderr) && (ios->f!=stdout) && (ios->f!=stdin))
        fclose(ios->f);
    return;
}



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Reset INCHI_IOSTREAM: set string buffer ptr to NULL (but do _not_ free memory)and close the file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
void inchi_ios_reset(INCHI_IOSTREAM* ios)
{
    ios->s.pStr = NULL; 
    ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;
    if ( (ios->f) && (ios->f!=stderr) && (ios->f!=stdout) && (ios->f!=stdin))
        fclose(ios->f);
    return;
}





/*******************************************************************/
int inchi_ios_str_getc(INCHI_IOSTREAM *ios)
{

    if (ios->type==INCHI_IOSTREAM_STRING)
    {
        if ( ios->s.nPtr < ios->s.nUsedLength ) 
        {
            return (int)ios->s.pStr[ios->s.nPtr++];
        }
        return EOF;
    }
    
    else if (ios->type==INCHI_IOSTREAM_FILE)
    {
        return fgetc( ios->f );
    }
    
    /* error */
    return EOF;
}



/*******************************************************************/
char *inchi_ios_str_gets(char *szLine, int len, INCHI_IOSTREAM *f)
{
int  length=0, c=0;
    if ( -- len < 0 ) 
    {
        return NULL;
    }
    while ( length < len && EOF != (c = inchi_ios_str_getc( f )) ) 
    {
        szLine[length++] = (char)c;
        if ( c == '\n' )
            break;
    }
    if ( !length && EOF == c ) 
    {
        return NULL;
    }
    szLine[length] = '\0';
    return szLine;
}



/********************************************************************************/
/* read up to len or tab or LF; if empty read next until finds non-empty line   */
/* remove leading and trailing white spaces; keep zero termination              */
/********************************************************************************/
char *inchi_ios_str_getsTab( char *szLine, int len, INCHI_IOSTREAM *f )
{
int  length=0, c=0;
    if ( --len < 0 ) 
    {
        return NULL;
    }
    while ( length < len && EOF != (c = inchi_ios_str_getc(f)) ) 
    {
        if ( c == '\t' )
            c = '\n';
        szLine[length++] = (char)c;
        if ( c == '\n' )
            break;
    }
    if ( !length && EOF == c ) 
    {
        return NULL;
    }
    szLine[length] = '\0';
    return szLine;
}


/*******************************************************************/
int inchi_ios_gets( char *szLine, int len, INCHI_IOSTREAM *f, int *bTooLongLine )
{
int  length;
char *p;
    do 
    {
        p = inchi_ios_str_gets( szLine, len-1, f );
        if ( !p ) 
        {
            *bTooLongLine = 0;
            return -1; /* end of file or cannot read */
        }
        szLine[len-1] = '\0';
        /*
        *bTooLongLine = !strchr( szLine, '\n' );
        */
        p = strchr( szLine, '\n' );
        *bTooLongLine = ( !p && ((int)strlen(szLine)) == len-2 );
        LtrimRtrim( szLine, &length );
    } while ( !length );

    return length;
}


/*******************************************************************/
/* read up to len or tab or LF; if empty read next until finds non-empty line   */
/* remove leading and trailing white spaces; keep zero termination */
/*******************************************************************/
int inchi_ios_getsTab( char *szLine, int len, INCHI_IOSTREAM *f, int *bTooLongLine )
{
int  length;
char *p;
    do 
    {
        p = inchi_ios_str_getsTab( szLine, len-1, f );
        if ( !p ) 
        {
            *bTooLongLine = 0;
            return -1; /* end of file or cannot read */
        }
        szLine[len-1] = '\0';
        /*
        *bTooLongLine = !strchr( szLine, '\n' );
        */
        p = strchr( szLine, '\n' );
        *bTooLongLine = ( !p && ((int)strlen(szLine)) == len-2 );
        LtrimRtrim( szLine, &length );
    } while ( !length );
    return length;
}

/*******************************************************************/
int inchi_ios_getsTab1( char *szLine, int len, INCHI_IOSTREAM *f, int *bTooLongLine )
{
int  length;
char *p;
    /*do {*/
        p = inchi_ios_str_getsTab( szLine, len-1, f );
        if ( !p ) 
        {
            *bTooLongLine = 0;
            return -1; /* end of file or cannot read */
        }
        szLine[len-1] = '\0';
        /*
        *bTooLongLine = !strchr( szLine, '\n' );
        */
        p = strchr( szLine, '\n' );
        *bTooLongLine = ( !p && ((int)strlen(szLine)) == len-2 );
        LtrimRtrim( szLine, &length );
    /*} while ( !length );*/
    return length;
}






/*****************************************************************/
int inchi_ios_print( INCHI_IOSTREAM * ios, const char* lpszFormat, ... )
{
int ret=0, ret2=0;
va_list argList;

    if (!ios) 
        return -1;
    
    if (ios->type == INCHI_IOSTREAM_STRING) 
    {    
        /* output to string buffer */        
        int max_len;
        my_va_start( argList, lpszFormat );
        max_len = GetMaxPrintfLength( lpszFormat, argList);
        va_end( argList );
        if ( max_len >= 0 ) 
        {
            if ( ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len ) 
            {
                /* enlarge output string */
                int  nAddLength = inchi_max( INCHI_ADD_STR_LEN, max_len );
                char *new_str = 
                    (char *)inchi_calloc( ios->s.nAllocatedLength + nAddLength, sizeof(new_str[0]) );
                if ( new_str ) 
                {
                    if ( ios->s.pStr ) 
                    {
                        if ( ios->s.nUsedLength > 0 ) 
                            memcpy( new_str, ios->s.pStr, sizeof(new_str[0])* ios->s.nUsedLength );
                        inchi_free( ios->s.pStr );
                    }
                    ios->s.pStr              = new_str;
                    ios->s.nAllocatedLength += nAddLength;
                } 
                else return -1; /* failed */
            }
            /* output */
            my_va_start( argList, lpszFormat );
            ret = vsprintf( ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList );
            va_end(argList);
            if ( ret >= 0 ) 
                ios->s.nUsedLength += ret;
#ifdef INCHI_LIB
            if( FWPRINT )
            {
                my_va_start( argList, lpszFormat );
                FWPRINT( lpszFormat, argList );
                va_end( argList );
            }
#endif
            return ret;
        }
        return -1;
    }
    
    else if (ios->type == INCHI_IOSTREAM_FILE)
    {
        /* output to file */
        if (ios->f) 
        {
            my_va_start( argList, lpszFormat );
            ret = vfprintf( ios->f, lpszFormat, argList );
            va_end( argList );
        } 
        else 
        {
            my_va_start( argList, lpszFormat );
            ret2 = vfprintf( stdout, lpszFormat, argList );
            va_end( argList );
        }
#ifdef INCHI_LIB
        if( FWPRINT )
        {
            my_va_start( argList, lpszFormat );
            FWPRINT( lpszFormat, argList );
            va_end( argList );
        }
#endif
        return ret? ret : ret2;
    }


    /* no output */
    return 0;
}




/**********************************************************************/
/* This function's output should not be displayed in the output pane  */
/**********************************************************************/
int inchi_ios_print_nodisplay( INCHI_IOSTREAM * ios, const char* lpszFormat, ... )
{
va_list argList;

    if (!ios) 
        return -1;
    
    if (ios->type == INCHI_IOSTREAM_STRING) 
    {
        /* output to string buffer */
        int ret=0, max_len;
        my_va_start( argList, lpszFormat );
        max_len = GetMaxPrintfLength( lpszFormat, argList);
        va_end( argList );
        if ( max_len >= 0 ) 
        {
            if ( ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len ) 
            {
                /* enlarge output string */
                int  nAddLength = inchi_max( INCHI_ADD_STR_LEN, max_len );
                char *new_str = (char *)inchi_calloc( ios->s.nAllocatedLength + nAddLength, sizeof(new_str[0]) );
                if ( new_str ) 
                {
                    if ( ios->s.pStr ) 
                    {
                        if ( ios->s.nUsedLength > 0 ) 
                        {
                            memcpy( new_str, ios->s.pStr, sizeof(new_str[0])*ios->s.nUsedLength );
                        }
                        inchi_free( ios->s.pStr );
                    }
                    ios->s.pStr              = new_str;
                    ios->s.nAllocatedLength += nAddLength;
                } else 
                {
                    return -1; /* failed */
                }
            }
            /* output */
            my_va_start( argList, lpszFormat );
            ret = vsprintf( ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList );
            va_end(argList);
            if ( ret >= 0 ) 
            {
                ios->s.nUsedLength += ret;
            }
            return ret;
        }
        return -1;
    }

    else if (ios->type == INCHI_IOSTREAM_FILE)
    {
        my_va_start( argList, lpszFormat );
        inchi_print_nodisplay( ios->f, lpszFormat, argList);
        va_end(argList);        
    }

    /* no output */
    return 0;
}




/*****************************************************************/
/* Print to string buffer or to file+stderr */
int inchi_ios_eprint( INCHI_IOSTREAM * ios, const char* lpszFormat, ... )
{
int ret=0, ret2=0;
va_list argList;

    if (!ios) 
        return -1;

    if (ios->type == INCHI_IOSTREAM_STRING) /* was #if ( defined(INCHI_LIBRARY) || defined(INCHI_STANDALONE_EXE) ) */
    {
        /* output to string buffer */
        int max_len, nAddLength = 0;
        char *new_str = NULL;

        my_va_start( argList, lpszFormat );
        max_len = GetMaxPrintfLength( lpszFormat, argList);
        va_end( argList );

        if ( max_len >= 0 ) 
        {
            if ( ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len ) 
            {
                /* enlarge output string */
                nAddLength = inchi_max( INCHI_ADD_STR_LEN, max_len );
                new_str = (char *)inchi_calloc( ios->s.nAllocatedLength + nAddLength, sizeof(new_str[0]) );
                if ( new_str ) 
                {
                    if ( ios->s.pStr ) 
                    {
                        if ( ios->s.nUsedLength > 0 ) 
                        {
                            memcpy( new_str, ios->s.pStr, sizeof(new_str[0])* ios->s.nUsedLength );
                        }
                        inchi_free( ios->s.pStr );
                    }
                    ios->s.pStr              = new_str;
                    ios->s.nAllocatedLength += nAddLength;
                } 
                else 
                {
                    return -1; /* failed */
                }
            }
         
            /* output */
            my_va_start( argList, lpszFormat );
            ret = vsprintf( ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList );
            va_end(argList);
            if ( ret >= 0 ) 
            {
                ios->s.nUsedLength += ret;
            }
            return ret;
        }
        return -1;
    }

    else if (ios->type == INCHI_IOSTREAM_FILE)
    {
        if ( ios->f) 
        {
            /* output to plain file */
            my_va_start( argList, lpszFormat );
            ret = inchi_vfprintf( ios->f, lpszFormat, argList ); 
            va_end( argList );
/*^^^  No output to stderr from within dll or GUI program */
#if ( !defined(INCHI_LIBRARY) && !defined(INCHI_LIB) )
            if ( ios->f != stderr ) 
            { 
                my_va_start( argList, lpszFormat );
                ret2 = vfprintf( stderr, lpszFormat, argList );
                va_end( argList );
            }
#endif
            return ret? ret : ret2;
        }
    } 

    /* no output */
    return 0;
}




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


                                    PLAIN FILE OPERATIONS 


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/* Print to file, echoing to stderr */
int inchi_fprintf( FILE* f, const char* lpszFormat, ... )
{
int ret=0, ret2=0;
va_list argList;
    if (f) 
    {
        my_va_start( argList, lpszFormat );
        ret = inchi_vfprintf( f, lpszFormat, argList ); 
        va_end( argList );
/*^^^  No output to stderr from within dll or GUI program */
#if ( !defined(INCHI_LIBRARY) && !defined(INCHI_LIB) )
        if ( f != stderr ) 
        { 
            my_va_start( argList, lpszFormat );
            ret2 = vfprintf( stderr, lpszFormat, argList );
            va_end( argList );
        }
#endif
        return ret? ret : ret2;
    }
    return 0;
}



/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/* Print to file */
int inchi_vfprintf( FILE* f, const char* lpszFormat, va_list argList )
{
int ret=0;
    if ( f == stderr && lpszFormat && lpszFormat[0] && '\r' == lpszFormat[strlen(lpszFormat)-1] ) 
    {
#define CONSOLE_LINE_LEN 80
#ifndef INCHI_ANSI_ONLY
        char szLine[CONSOLE_LINE_LEN];
        ret = _vsnprintf( szLine, CONSOLE_LINE_LEN-1, lpszFormat, argList );
        if ( ret < 0 ) 
        {
            /*  output is longer than the console line */
            /*^^^ Fixed bug: (CONSOLE_LINE_LEN-4) --> (CONSOLE_LINE_LEN-4-1) 11-22-08 IPl */
            strcpy(szLine+CONSOLE_LINE_LEN-5, "...\r");
        }
        fputs( szLine, f );
#else
        ret = vfprintf( f, lpszFormat, argList );
#endif
#undef CONSOLE_LINE_LEN
    } 
    else 
    {
        ret = vfprintf( f, lpszFormat, argList );
    }
    return ret;
}



/**********************************************************************/
/* This function's output should not be displayed in the output pane  */
/**********************************************************************/
int inchi_print_nodisplay( FILE* f, const char* lpszFormat, ... )
{
int ret=0;
va_list argList;
FILE* fi;
    if (f)	
        fi = f;
    else	
        fi = stdout;
    my_va_start( argList, lpszFormat );
    ret = vfprintf( fi, lpszFormat, argList );
    return ret;
}



#if ( FIX_READ_LONG_LINE_BUG == 1 )
/********************************************************************/
int inchi_fgetsLfTab( char *szLine, int len, FILE *f )
{
    int  length;
    char *p;
    char szSkip[256];
    int  bTooLongLine = 0;
    do {
        p = inchi_fgetsTab( szLine, len, f );
        if ( !p ) {
            return -1; /* end of file or cannot read */
        }
        bTooLongLine = ( (int)strlen(szLine) == len-1 && szLine[len-2] != '\n' );
        LtrimRtrim( szLine, &length );
    } while ( !length );
    if ( bTooLongLine ) {
        while ( p = inchi_fgetsTab( szSkip, sizeof(szSkip)-1, f ) ) {
            if ( strchr( szSkip, '\n' ) )
                break;
        }
    }
    return length;
}
#else
/********************************************************************/
int inchi_fgetsLfTab( char *szLine, int len, FILE *f )
{
    int  length;
    char *p;
    char szSkip[256];
    int  bTooLongLine = 0;
    do {
        p = inchi_fgetsTab( szLine, len-1, f );
        if ( !p ) {
            return -1; /* end of file or cannot read */
        }
        szLine[len-1] = '\0';
        /*
        bTooLongLine = !strchr( szLine, '\n' );
        */
        bTooLongLine = ( !p && ((int)strlen(szLine)) == len-2 );
        LtrimRtrim( szLine, &length );
    } while ( !length );
    if ( bTooLongLine ) {
        while ( p = inchi_fgetsTab( szSkip, sizeof(szSkip)-1, f ) ) {
            szSkip[sizeof(szSkip)-1] = '\0';
            if ( strchr( szSkip, '\n' ) )
                break;
        }
    }
    return length;
}
#endif


/*******************************************************************/
/* read up to len or tab or LF; if empty read next until finds non-empty line   */
/* remove leading and trailing white spaces; keep zero termination */
/*******************************************************************/
char *inchi_fgetsTab( char *szLine, int len, FILE *f )
{
    int  length=0, c=0;
    len --;
    while ( length < len && EOF != (c = fgetc( f )) ) {
        if ( c == '\t' )
            c = '\n';
        szLine[length++] = (char)c;
        if ( c == '\n' )
            break;
    }
    if ( !length && EOF == c ) {
        return NULL;
    }
    szLine[length] = '\0';
    return szLine;
}



/******************************************************************/
/* read not more than line_len bytes from an lf-terminated line   */
/* if input line is too long quietly ignore the rest of the line  */
char* inchi_fgetsLf( char* line, int line_len, FILE* inp )
{
    char *p, *q;
    memset( line, 0, line_len );
    if ( NULL != (p = fgets( line, line_len, inp ) ) && NULL == strchr(p, '\n' ) ){
        char temp[64]; /* bypass up to '\n' or up to end of file whichever comes first*/
        while ( NULL != fgets( temp, sizeof(temp), inp ) && NULL == strchr(temp,'\n') )
            ;
    }
    if ( p && (q = strchr(line, '\r')) ) { /*  fix CR CR LF line terminator. */
        q[0] = '\n';
        q[1] = '\0';
    }    
    return p;
}










/*****************************************************************
 *
 *  Estimate printf string length
 *  
 *  The code is based on Microsoft Knowledge Base article Q127038:
 *  "FIX: CString::Format Gives Assertion Failed, Access Violation"
 *  (Related to Microsoft Visual C++, 32-bit Editions, versions 2.0, 2.1)
 *
 *****************************************************************/

#define FORCE_ANSI      0x10000
#define FORCE_UNICODE   0x20000

/* formatting (using wsprintf style formatting)*/
int GetMaxPrintfLength( const char *lpszFormat, va_list argList)
{
     /*ASSERT(AfxIsValidString(lpszFormat, FALSE));*/
     const char * lpsz;
     int nMaxLen, nWidth, nPrecision, nModifier, nItemLen;

     nMaxLen = 0;
     /* make a guess at the maximum length of the resulting string */
     for ( lpsz = lpszFormat; *lpsz; lpsz ++ )
     {
          /* handle '%' character, but watch out for '%%' */
          if (*lpsz != '%' || *( ++ lpsz ) == '%')
          {
               nMaxLen += 1;
               continue;
          }

          nItemLen = 0;

          /*  handle '%' character with format */
          nWidth = 0;
          for (; *lpsz; lpsz ++ )
          {
               /* check for valid flags */
               if (*lpsz == '#')
                    nMaxLen += 2;   /* for '0x' */
               else if (*lpsz == '*')
                    nWidth = va_arg(argList, int);
               else if (*lpsz == '-' || *lpsz == '+' || *lpsz == '0'
                        || *lpsz == ' ')
                           ;
               else /* hit non-flag character */
                          break;
          }
          /* get width and skip it */
          if (nWidth == 0)
          {
               /* width indicated by */
               nWidth = atoi(lpsz);
               for (; *lpsz && isdigit(*lpsz); lpsz ++ )
                     ;
          }
          /*ASSERT(nWidth >= 0);*/
          if ( nWidth < 0 )
              goto exit_error; /* instead of exception */

          nPrecision = 0;
          if (*lpsz == '.')
          {
               /* skip past '.' separator (width.precision)*/
               lpsz ++;

               /* get precision and skip it*/
               if (*lpsz == '*')
               {
                    nPrecision = va_arg(argList, int);
                    lpsz ++;
               }
               else
               {
                    nPrecision = atoi(lpsz);
                    for (; *lpsz && isdigit(*lpsz); lpsz ++)
                            ;
               }
              if ( nPrecision < 0 )
                  goto exit_error; /* instead of exception */
          }

          /* should be on type modifier or specifier */
          nModifier = 0;
          switch (*lpsz)
          {
          /* modifiers that affect size */
          case 'h':
               switch ( lpsz[1] ) {
               case 'd':
               case 'i':
               case 'o':
               case 'x':
               case 'X':
               case 'u':
                   /* short unsigned, short double, etc. -- added to the original MS example */
                   /* ignore the fact that these modifiers do affect size */
                   lpsz ++;
                   break;
               default:
                   nModifier = FORCE_ANSI;
                   lpsz ++;
                   break;
               }
               break;
          case 'l':
               switch ( lpsz[1] ) {
               case 'd':
               case 'i':
               case 'o':
               case 'x':
               case 'X':
               case 'u':
               case 'f': /* long float -- post ANSI C */
                   /* long unsigned, long double, etc. -- added to the original MS example */
                   /* ignore the fact that these modifiers do affect size */
                   lpsz ++;
                   break;
               default:
                   /*
                   nModifier = FORCE_UNICODE;
                   lpsz ++;
                   break;
                   */
                   goto exit_error;  /* no UNICODE, please */
               }
               break;
          /* modifiers that do not affect size */
          case 'F':
          case 'N':
          case 'L':
               lpsz ++;
               break;
          }

          /* now should be on specifier */
          switch (*lpsz | nModifier)
          {
          /* single characters*/
          case 'c':
          case 'C':
               nItemLen = 2;
               va_arg(argList, int);
               break;
          case 'c'|FORCE_ANSI:
          case 'C'|FORCE_ANSI:
               nItemLen = 2;
               va_arg(argList, int);
               break;
          case 'c'|FORCE_UNICODE:
          case 'C'|FORCE_UNICODE:
               goto exit_error;  /* no UNICODE, please */
               /*
               nItemLen = 2;
               va_arg(argList, wchar_t);
               break;
               */

          /* strings*/
          case 's':
          case 'S':
               nItemLen = strlen(va_arg(argList, char*));
               nItemLen = inchi_max(1, nItemLen);
               break;
          case 's'|FORCE_ANSI:
          case 'S'|FORCE_ANSI:
               nItemLen = strlen(va_arg(argList, char*));
               nItemLen = inchi_max(1, nItemLen);
               break;

          case 's'|FORCE_UNICODE:
          case 'S'|FORCE_UNICODE:
               goto exit_error;  /* no UNICODE, please */
               /*
               nItemLen = wcslen(va_arg(argList, wchar_t*));
               nItemLen = inchi_max(1, nItemLen);
               break;
               */

          }

          /* adjust nItemLen for strings */
          if (nItemLen != 0)
          {
               nItemLen = inchi_max(nItemLen, nWidth);
               if (nPrecision != 0)
                    nItemLen = inchi_min(nItemLen, nPrecision);
          }
          else
          {
               switch (*lpsz)
               {
               /* integers */
               case 'd':
               case 'i':
               case 'u':
               case 'x':
               case 'X':
               case 'o':
                    va_arg(argList, int);
                    nItemLen = 32;
                    nItemLen = inchi_max(nItemLen, nWidth+nPrecision);
                    break;

               case 'e':
               case 'f':
               case 'g':
               case 'G':
                    va_arg(argList, double);
                    nItemLen = 32;
                    nItemLen = inchi_max(nItemLen, nWidth+nPrecision);
                    break;

               case 'p':
                    va_arg(argList, void*);
                    nItemLen = 32;
                    nItemLen = inchi_max(nItemLen, nWidth+nPrecision);
                    break;

               /* no output */
               case 'n':
                    va_arg(argList, int*);
                    break;

               default:
                   /*ASSERT(FALSE);*/  /* unknown formatting option*/
                   goto exit_error; /* instead of exception */
               }
          }

          /* adjust nMaxLen for output nItemLen */
          nMaxLen += nItemLen;
     }
     return nMaxLen;

exit_error:
     return -1; /* wrong format */
} 

