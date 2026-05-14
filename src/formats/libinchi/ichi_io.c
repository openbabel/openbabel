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

#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <time.h>

#include "mode.h"
#include "ichi_io.h"
#include "ichicomp.h"
#include "ichidrp.h"
#include "strutil.h"
#include "util.h"

#include "bcf_s.h"

#ifndef INCHI_ADD_STR_LEN
#define INCHI_ADD_STR_LEN   32768
#endif

#ifdef TARGET_LIB_FOR_WINCHI
extern void(*FWPRINT) (const char* format, va_list argptr);
extern void(*FWPUSH) (const char* s);

#endif


/* Internal functions */

int inchi_ios_str_getc(INCHI_IOSTREAM* ios);
char* inchi_ios_str_gets(char* szLine, int len, INCHI_IOSTREAM* ios);
char* inchi_ios_str_getsTab(char* szLine, int len, INCHI_IOSTREAM* ios);
int GetMaxPrintfLength(const char* lpszFormat, va_list argList);
char* inchi_fgetsTab(char* szLine, int len, FILE* f);
int inchi_vfprintf(FILE* f, const char* lpszFormat, va_list argList);


/*
    INCHI_IOSTREAM OPERATIONS
*/


/****************************************************************************
 Init INCHI_IOSTREAM
****************************************************************************/
void inchi_ios_init(INCHI_IOSTREAM* ios, int io_type, FILE* f)
{
    memset(ios, 0, sizeof(*ios)); /* djb-rwth: memset_s C11/Annex K variant? */
    switch (io_type)
    {
    case INCHI_IOS_TYPE_FILE:
        ios->type = INCHI_IOS_TYPE_FILE;
        break;
    case INCHI_IOS_TYPE_STRING:
    default:
        ios->type = INCHI_IOS_TYPE_STRING;
        break;
    }
    ios->f = f;
    return;
}


/****************************************************************************
 Make a copy of INCHI_IOSTREAM
****************************************************************************/
int inchi_ios_create_copy(INCHI_IOSTREAM* ios, INCHI_IOSTREAM* ios0)
{
    if (ios) /* djb-rwth: fixing a NULL pointer dereference */
    {
        memset(ios, 0, sizeof(*ios)); /* djb-rwth: memset_s C11/Annex K variant? */
        ios->type = ios0->type;

        if (ios->type == INCHI_IOS_TYPE_STRING)
        {
            if (ios->s.pStr)
            {
                inchi_free(ios->s.pStr);
            }
            ios->s.pStr = (char*)inchi_calloc(ios0->s.nAllocatedLength, sizeof(char));
            if (ios->s.pStr)
            {
                ios->s.nUsedLength = ios0->s.nUsedLength;
                ios->s.nPtr = ios0->s.nPtr;
            }
            else
            {
                return -1; /* no memory */
            }
        }
        ios->f = ios0->f;
    }

    return 0;
}


/****************************************************************************
    If INCHI_IOSTREAM type is INCHI_IOS_TYPE_STRING
    and associated file exists, flush the string buffer
    to that file, then free the buffer.
    If INCHI_IOSTREAM type is INCHI_IOS_TYPE_FILE,
    just flush the file.
****************************************************************************/
void inchi_ios_flush(INCHI_IOSTREAM* ios)
{

    if (ios->type == INCHI_IOS_TYPE_STRING)
    {
        if (ios->s.pStr)
        {
            if (ios->s.nUsedLength > 0)
            {
                if (ios->f)
                {
                    fprintf(ios->f, "%-s", ios->s.pStr);
                    fflush(ios->f);
                }
                inchi_free(ios->s.pStr);
                ios->s.pStr = NULL;
                ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;
            }
        }
    }

    else if (ios->type == INCHI_IOS_TYPE_FILE)
    {
        /* output to plain file: just flush it. */
        fflush(ios->f);
    }

    return;
}


/****************************************************************************
    If INCHI_IOSTREAM type is INCHI_IOS_TYPE_STRING,
    flush INCHI_IOSTREAM string buffer to associated file
    (if non-NULL) echoing to another file (supplied as
    parameter; typically, stderr); then free buffer.
    If INCHI_IOSTREAM type is INCHI_IOS_TYPE_FILE,
    just flush the both files.
****************************************************************************/
void inchi_ios_flush2(INCHI_IOSTREAM* ios, FILE* f2)
{

    if (ios->type == INCHI_IOS_TYPE_STRING)
    {
        if (ios->s.pStr)
        {
            if (ios->s.nUsedLength > 0)
            {
                if (ios->f)
                {
                    fprintf(ios->f, "%-s", ios->s.pStr);
                    fflush(ios->f);
                }
                if (f2 != ios->f)
                {
                    fprintf(f2, "%-s", ios->s.pStr);
                }

                inchi_free(ios->s.pStr);
                ios->s.pStr = NULL;
                ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;
            }
        }
    }

    else if (ios->type == INCHI_IOS_TYPE_FILE)
    {
        /* Output to plain file: just flush it. */
        if ((ios->f) && (ios->f != stderr) && (ios->f != stdout))
        {
            fflush(ios->f);
        }
        if (f2 && f2 != stderr && f2 != stdout)
        {
            fflush(f2);
        }
    }

    return;
}


/****************************************************************************
    Close INCHI_IOSTREAM: free string buffer and close associated file.
****************************************************************************/
void inchi_ios_close(INCHI_IOSTREAM* ios)
{
    if (NULL == ios)
    {
        return;
    }
    if (ios->s.pStr)
    {
        inchi_free(ios->s.pStr);
    }
    ios->s.pStr = NULL;
    ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;

    if (NULL != ios->f && stdout != ios->f && stderr != ios->f && stdin != ios->f)
    {
        fclose(ios->f);
    }

    return;
}


/****************************************************************************
    Reset INCHI_IOSTREAM: set string buffer ptr to NULL
    (but do _not_ free memory)and close associated file.
****************************************************************************/
void inchi_ios_reset(INCHI_IOSTREAM* ios)
{
    ios->s.pStr = NULL;
    ios->s.nUsedLength = ios->s.nAllocatedLength = ios->s.nPtr = 0;
    if (NULL != ios->f && stdout != ios->f && stderr != ios->f && stdin != ios->f)
    {
        fclose(ios->f);
    }

    return;
}


/****************************************************************************
    Reset INCHI_IOSTREAM: set string buffer ptr to NULL
    (after freeing memory) but do not close associated file.
****************************************************************************/
void inchi_ios_free_str(INCHI_IOSTREAM* ios)
{
    if (NULL == ios)
    {
        return;
    }
    if (ios->s.pStr) /* && ios->s.nAllocatedLength)*/
    {
        inchi_free(ios->s.pStr);
    }
    ios->s.pStr = NULL;
    ios->s.nUsedLength = 0;
    ios->s.nAllocatedLength = 0;
    ios->s.nPtr = 0;

    return;
}


/****************************************************************************
    [str] getc()
****************************************************************************/
int inchi_ios_str_getc(INCHI_IOSTREAM* ios)
{
    int c;
    if (ios->type == INCHI_IOS_TYPE_STRING)
    {
        if (ios->s.nPtr < ios->s.nUsedLength)
        {
            return (int)ios->s.pStr[ios->s.nPtr++];
        }
        return EOF;
    }

    else if (ios->type == INCHI_IOS_TYPE_FILE)
    {
        c = fgetc(ios->f);
        if (ferror(ios->f))
        {
            c = EOF;
        }
        return c;
    }

    /* error */
    return EOF;
}


/****************************************************************************
    [str] gets()
****************************************************************************/
char* inchi_ios_str_gets(char* szLine, int len, INCHI_IOSTREAM* f)
{
    int  length = 0, c = 0;
    if (--len < 0)
    {
        return NULL;
    }
    while (length < len && EOF != (c = inchi_ios_str_getc(f)))
    {
        szLine[length++] = (char)c;
        if (c == '\n')
        {
            break;
        }
    }
    if (!length && EOF == c)
    {
        return NULL;
    }
    szLine[length] = '\0';

    return szLine;
}

/****************************************************************************
    Read up to tab or LF but not more than len chars;
    if empty line, read further until a non-empty line;
    remove leading and trailing white spaces;
    keep zero termination
****************************************************************************/
char* inchi_ios_str_getsTab(char* szLine, int len, INCHI_IOSTREAM* f)
{
    int  length = 0, c = 0;
    if (--len < 0)
    {
        return NULL;
    }
    while (length < len && EOF != (c = inchi_ios_str_getc(f)))
    {
        if (c == '\t')
        {
            c = '\n';
        }
        szLine[length++] = (char)c;
        if (c == '\n')
        {
            break;
        }
    }
    if (!length && EOF == c)
    {
        return NULL;
    }
    szLine[length] = '\0';

    return szLine;
}


/****************************************************************************
    gets()
****************************************************************************/
int inchi_ios_gets(char* szLine,
    int len,
    INCHI_IOSTREAM* f,
    int* bTooLongLine)
{
    int  length;
    char* p;
    do
    {
        p = inchi_ios_str_gets(szLine, len - 1, f);
        if (!p)
        {
            *bTooLongLine = 0;
            return -1; /* end of file or cannot read */
        }
        szLine[len - 1] = '\0';
        /*
        *bTooLongLine = !strchr( szLine, '\n' );
        */
        p = strchr(szLine, '\n');
        *bTooLongLine = (!p && ((int)strlen(szLine)) == len - 2);
        lrtrim(szLine, &length);
    } while (!length);

    return length;
}


/****************************************************************************
    Read up to tab or LF but not more than len chars;
    if empty line, read further until a non-empty line;
    remove leading and trailing white spaces;
    keep zero termination
****************************************************************************/
int inchi_ios_getsTab(char* szLine,
    int len,
    INCHI_IOSTREAM* f,
    int* bTooLongLine)
{
    int  length;
    char* p;
    do
    {

        p = inchi_ios_str_getsTab(szLine, len - 1, f);
        if (!p)
        {
            *bTooLongLine = 0;
            return -1; /* end of file or cannot read */
        }
        szLine[len - 1] = '\0';
        /*
        *bTooLongLine = !strchr( szLine, '\n' );
        */
        p = strchr(szLine, '\n');
        *bTooLongLine = (!p && ((int)strlen(szLine)) == len - 2);
        lrtrim(szLine, &length);

    } while (!length);

    return length;
}


/****************************************************************************/
int inchi_ios_getsTab1(char* szLine,
    int len,
    INCHI_IOSTREAM* f,
    int* bTooLongLine)
{
    int  length;
    char* p;

    p = inchi_ios_str_getsTab(szLine, len - 1, f);
    if (!p)
    {
        *bTooLongLine = 0;
        return -1; /* end of file or cannot read */
    }
    szLine[len - 1] = '\0';
    p = strchr(szLine, '\n');
    *bTooLongLine = (!p && ((int)strlen(szLine)) == len - 2);
    lrtrim(szLine, &length);

    return length;
}


/****************************************************************************
 General procedure for printing to INCHI_IOSTREAM
****************************************************************************/
int inchi_ios_print(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
{
    int ret = 0, ret2 = 0;
    va_list argList;

    if (!ios)
    {
        return -1;
    }

    if (ios->type == INCHI_IOS_TYPE_STRING)
    {
        /* output to string buffer */
        int max_len;
        my_va_start(argList, lpszFormat);
        max_len = GetMaxPrintfLength(lpszFormat, argList);
        va_end(argList);
        if (max_len >= 0)
        {
            /* djb-rwth: fixing oss-fuzz issue #30152 */
            int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
            long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
            if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
            {
                /* enlarge output string */
                char* new_str = (char*)inchi_calloc(new_str_len, sizeof(char)); /* djb-rwth: cast operators added */
                if (new_str)
                {
                    if (ios->s.pStr)
                    {
                        if (ios->s.nUsedLength > 0)
                        {
                            memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
                        }
                        inchi_free(ios->s.pStr);
                    }
                    ios->s.pStr = new_str;
                    ios->s.nAllocatedLength += nAddLength;
                }
                else
                {
                    return -1; /* failed */
                }
            }
            /* output */
            my_va_start(argList, lpszFormat);
            ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument */
            va_end(argList);
            if (ret >= 0)
            {
                ios->s.nUsedLength += ret;
            }
#ifdef TARGET_LIB_FOR_WINCHI
#if 0
            if (FWPRINT)
            {
                my_va_start(argList, lpszFormat);
                FWPRINT(lpszFormat, argList);
                va_end(argList);
            }
#endif
#endif
            return ret;
        }
        return -1;
    }

    else if (ios->type == INCHI_IOS_TYPE_FILE)
    {
        /* output to file */
        if (ios->f)
        {
            my_va_start(argList, lpszFormat);
            ret = vfprintf(ios->f, lpszFormat, argList);
            va_end(argList);
        }
        else
        {
            my_va_start(argList, lpszFormat);
            ret2 = vfprintf(stdout, lpszFormat, argList);
            va_end(argList);
        }
#ifdef TARGET_LIB_FOR_WINCHI
        if (FWPRINT)
        {
            my_va_start(argList, lpszFormat);
            FWPRINT(lpszFormat, argList);
            va_end(argList);
        }
#endif
        return ret ? ret : ret2;
    }

    /* no output */
    return 0;
}


/****************************************************************************/
int push_to_winchi_text_window(INCHI_IOSTREAM* ios)
/*, const char* lpszFormat, ... ) */
{
#ifndef TARGET_LIB_FOR_WINCHI
    return -1;
#else
    if (ios->type != INCHI_IOS_TYPE_STRING)
    {
        return -1;
    }
    if (!ios || !ios->s.pStr)
    {
        return -1;
    }
    if (!FWPUSH)
    {
        return -1;
    }

    FWPUSH(ios->s.pStr);

    return 0;
#endif
}

/****************************************************************************
    This function's output should not be displayed in the output pane
****************************************************************************/
int inchi_ios_print_nodisplay( INCHI_IOSTREAM * ios,
                               const char* lpszFormat,
                               ... )
{
    va_list argList;

    if (!ios)
    {
        return -1;
    }

    if (ios->type == INCHI_IOS_TYPE_STRING)
    {
        /* output to string buffer */
        int ret = 0, max_len;
        my_va_start( argList, lpszFormat );
        max_len = GetMaxPrintfLength( lpszFormat, argList );
        va_end( argList );
        if (max_len >= 0)
        {
            /* djb-rwth: fixing oss-fuzz issue #30163 */
            int  nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
            long long new_str_len = (long long)ios->s.nAllocatedLength + (long long)nAddLength;
            if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
            {
                /* enlarge output string */
                char* new_str = (char*)inchi_calloc(new_str_len, sizeof(new_str[0])); /* djb-rwth: cast operators added */
                if (new_str)
                {
                    if (ios->s.pStr)
                    {
                        if (ios->s.nUsedLength > 0)
                        {
                            memcpy( new_str, ios->s.pStr, sizeof( new_str[0] )*ios->s.nUsedLength );
                        }
                        inchi_free( ios->s.pStr );
                    }
                    ios->s.pStr = new_str;
                    ios->s.nAllocatedLength += nAddLength;
                }
                else
                {
                    return -1; /* failed */
                }
            }
            /* output */
            /* djb-rwth: fixing oss-fuzz issue #67676 */
            my_va_start(argList, lpszFormat);
            ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList); /* djb-rwth: not using vsnprintf due to variable length argument; fixing GHI #71 */
            va_end(argList);
            if (ret >= 0)
            {
                ios->s.nUsedLength += ret;
            }
            return ret;
        }
        return -1;
    }

    else if (ios->type == INCHI_IOS_TYPE_FILE)
    {
        my_va_start( argList, lpszFormat );
        inchi_print_nodisplay( ios->f, lpszFormat, argList );
        va_end( argList );
    }

    /* no output */
    return 0;
}


/****************************************************************************
    This function's flushes previously hidden output and resets string stream
    returns n chars on success, otherwise -1
****************************************************************************/
int inchi_ios_flush_not_displayed(INCHI_IOSTREAM* ios)
{
    char* obuf = NULL;
    int ret;

    if (!ios)
    {
        return -1;
    }

    obuf = (char*)inchi_calloc((long long)ios->s.nUsedLength + 1, sizeof(char)); /* djb-rwth: cast operator added */

    if (!obuf)
    {
        return -1;
    }

    strcpy(obuf, ios->s.pStr);
    ios->s.nUsedLength = 0;
    ret = inchi_ios_print(ios, "%s", obuf);
    inchi_free(obuf);

    return ret;
}


/****************************************************************************
 Print to string buffer or to file+stderr
****************************************************************************/
int inchi_ios_eprint(INCHI_IOSTREAM* ios, const char* lpszFormat, ...)
{
    int ret = 0, ret2 = 0;
    va_list argList;

    if (!ios)
    {
        return -1;
    }

    if (ios->type == INCHI_IOS_TYPE_STRING)
        /* was #if ( defined(TARGET_API_LIB) || defined(INCHI_STANDALONE_EXE) ) */
    {
        /* output to string buffer */
        int max_len, nAddLength = 0;
        char* new_str = NULL;

        my_va_start(argList, lpszFormat);
        max_len = GetMaxPrintfLength(lpszFormat, argList);
        va_end(argList);

        if (max_len >= 0)
        {
            if (ios->s.nAllocatedLength - ios->s.nUsedLength <= max_len)
            {
                /* enlarge output string */
                nAddLength = inchi_max(INCHI_ADD_STR_LEN, max_len);
                new_str = (char*)inchi_calloc((long long)ios->s.nAllocatedLength + (long long)nAddLength, sizeof(new_str[0])); /* djb-rwth: cast operators added */
                if (new_str)
                {
                    if (ios->s.pStr)
                    {
                        if (ios->s.nUsedLength > 0)
                        {
                            memcpy(new_str, ios->s.pStr, sizeof(new_str[0]) * ios->s.nUsedLength);
                        }
                        inchi_free(ios->s.pStr);
                    }
                    ios->s.pStr = new_str;
                    ios->s.nAllocatedLength += nAddLength;
                }
                else
                {
                    return -1; /* failed */
                }
            }

            /* output */
            my_va_start(argList, lpszFormat);
            ret = vsprintf(ios->s.pStr + ios->s.nUsedLength, lpszFormat, argList);
            va_end(argList);
            if (ret >= 0)
            {
                ios->s.nUsedLength += ret;
            }
            return ret;
        }
        return -1;
    }

    else if (ios->type == INCHI_IOS_TYPE_FILE)
    {
        if (ios->f)
        {
            /* output to plain file */
            my_va_start(argList, lpszFormat);
            ret = inchi_vfprintf(ios->f, lpszFormat, argList);
            va_end(argList);
            /*  No output to stderr from within dll or GUI program */
#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
            if (ios->f != stderr)
            {
                my_va_start(argList, lpszFormat);
                ret2 = vfprintf(stderr, lpszFormat, argList);
                va_end(argList);
            }
#endif
            return ret ? ret : ret2;
        }
    }

    /* no output */
    return 0;
}


/*    PLAIN FILE OPERATIONS */


/****************************************************************************
 Print to file, echoing to stderr
****************************************************************************/
int inchi_fprintf(FILE* f, const char* lpszFormat, ...)
{
    int ret = 0, ret2 = 0;
    va_list argList;
    if (f)
    {
        my_va_start(argList, lpszFormat);
        ret = inchi_vfprintf(f, lpszFormat, argList);
        va_end(argList);
        /*  No output to stderr from within dll or GUI program */
#if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
        if (f != stderr)
        {
            my_va_start(argList, lpszFormat);
            ret2 = vfprintf(stderr, lpszFormat, argList);
            va_end(argList);
        }
#endif
        return ret ? ret : ret2;
    }

    return 0;
}


/****************************************************************************
 Print to file
****************************************************************************/
int inchi_vfprintf(FILE* f, const char* lpszFormat, va_list argList)
{
    int ret = 0;
    if (lpszFormat && lpszFormat[0]) /* djb-rwth: condition added as lpszFormat == 0 may lead to undefined ret value */
    {
        if (f == stderr && '\r' == lpszFormat[strlen(lpszFormat) - 1])
        {
#define CONSOLE_LINE_LEN 80
#ifndef COMPILE_ANSI_ONLY
            char szLine[CONSOLE_LINE_LEN];

            ret = _vsnprintf(szLine, CONSOLE_LINE_LEN - 1, lpszFormat, argList);

            if (ret < 0)
            {
                /*  output is longer than the console line */
                /* Fixed bug: (CONSOLE_LINE_LEN-4) --> (CONSOLE_LINE_LEN-4-1) 11-22-08 IPl */
                strcpy(szLine + CONSOLE_LINE_LEN - 5, "...\r");
            }

            fputs(szLine, f);
#else
            ret = vfprintf(f, lpszFormat, argList);
#endif
#undef CONSOLE_LINE_LEN
        }
        else
        {
            ret = vfprintf(f, lpszFormat, argList);
        }
    }

    return ret;
}


/****************************************************************************
    This function's output should not be displayed in the output pane
****************************************************************************/
int inchi_print_nodisplay(FILE* f, const char* lpszFormat, ...)
{
    int ret = 0;
    va_list argList;
    FILE* fi;

    if (f)
    {
        fi = f;
    }
    else
    {
        fi = stdout;
    }

    my_va_start(argList, lpszFormat);
    ret = vfprintf(fi, lpszFormat, argList);
    va_end(argList);
    return ret;
}


#if ( FIX_READ_LONG_LINE_BUG == 1 )


/****************************************************************************/
int inchi_fgetsLfTab(char* szLine, int len, FILE* f)
{
    int  length;
    char* p;
    char szSkip[256];
    int  bTooLongLine = 0;
    do
    {
        p = inchi_fgetsTab(szLine, len, f);
        if (!p)
        {
            return -1; /* end of file or cannot read */
        }
        bTooLongLine = ((int)strlen(szLine) == len - 1 && szLine[len - 2] != '\n');
        lrtrim(szLine, &length);
    } while (!length);
    if (bTooLongLine)
    {
        while ((p = inchi_fgetsTab(szSkip, sizeof(szSkip) - 1, f))) /* djb-rwth: ignoring LLVM warning: function returning value */
        {
            if (strchr(szSkip, '\n'))
                break;
        }
    }

    return length;
}


#else


/****************************************************************************/
int inchi_fgetsLfTab(char* szLine, int len, FILE* f)
{
    int  length;
    char* p;
    char szSkip[256];
    int  bTooLongLine = 0;
    do
    {
        p = inchi_fgetsTab(szLine, len - 1, f);
        if (!p)
        {
            return -1; /* end of file or cannot read */
        }
        szLine[len - 1] = '\0';
        /*
        bTooLongLine = !strchr( szLine, '\n' );
        */
        bTooLongLine = (!p && ((int)strlen(szLine)) == len - 2);
        lrtrim(szLine, &length);
    } while (!length);
    if (bTooLongLine)
    {
        while (p = inchi_fgetsTab(szSkip, sizeof(szSkip) - 1, f))
        {
            szSkip[sizeof(szSkip) - 1] = '\0';
            if (strchr(szSkip, '\n'))
                break;
        }
    }
    return length;
}
#endif


/****************************************************************************
    Read up to tab or LF but not more than len chars;
    if empty line, read further until a non-empty line;
    remove leading and trailing white spaces;
    keep zero termination
****************************************************************************/
char* inchi_fgetsTab(char* szLine, int len, FILE* f)
{
    int  length = 0, c = 0;
    len--;
    while (length < len && EOF != (c = fgetc(f)))
    {
        if (c == '\t')
        {
            c = '\n';
        }
        szLine[length++] = (char)c;
        if (c == '\n')
        {
            break;
        }
    }
    if (!length && EOF == c)
    {
        return NULL;
    }
    szLine[length] = '\0';

    return szLine;
}


/****************************************************************************
    Read up to LF but not more than line_len bytes;
    if input line is too long, quietly ignore the rest of the line
****************************************************************************/
char* inchi_fgetsLf(char* line, int line_len, INCHI_IOSTREAM* inp_stream)
{
    char* p = NULL, * q;
    FILE* finp = NULL;

    if (inp_stream->type == INCHI_IOS_TYPE_FILE)
    {
        /* Read from file */
        finp = inp_stream->f;
        memset(line, 0, line_len); /* djb-rwth: memset_s C11/Annex K variant? */
        if (NULL != (p = fgets(line, line_len, finp)) &&
            NULL == strchr(p, '\n'))
        {
            char temp[64];
            /* bypass up to '\n' or up to end of file whichever comes first */
            while (NULL != fgets(temp, sizeof(temp), finp) && NULL == strchr(temp, '\n'))
            {
                ;
            }
        }
    }
    else if (inp_stream->type == INCHI_IOS_TYPE_STRING)
    {
        /* Read from supplied string representing Molfile */
        memset(line, 0, line_len); /* djb-rwth: memset_s C11/Annex K variant? */
        if (NULL != (p = inchi_sgets(line, line_len, inp_stream)) &&
            NULL == strchr(p, '\n'))
        {
            char temp[64];
            /* bypass up to '\n' or up to end of file whichever comes first */
            while (NULL != inchi_sgets(temp, sizeof(temp), inp_stream) && NULL == strchr(temp, '\n'))
            {
                ;
            }
        }
    }
    else
    {
        ;
    }

    if (p)
    {
        if ((q = strchr(line, '\r'))) /* djb-rwth: addressing LLVM warning */
        {
            /*  fix CR CR LF line terminator. */
            q[0] = '\n';
            q[1] = '\0';
        }
    }

    return p;
}


/****************************************************************************
    Estimate printf string length.

    The code is based on Microsoft Knowledge Base article Q127038:
    "FIX: CString::Format Gives Assertion Failed, Access Violation"
    (Related to Microsoft Visual C++, 32-bit Editions, versions 2.0, 2.1)
****************************************************************************/

#define FORCE_ANSI      0x10000
#define FORCE_UNICODE   0x20000

/****************************************************************************
 Formatting (using wsprintf style formatting)
****************************************************************************/
int GetMaxPrintfLength(const char* lpszFormat, va_list argList)
{
    /*ASSERT(AfxIsValidString(lpszFormat, FALSE));*/
    const char* lpsz;
    int nMaxLen, nWidth, nPrecision, nModifier, nItemLen;

    nMaxLen = 0;
    /* make a guess at the maximum length of the resulting string */
    for (lpsz = lpszFormat; *lpsz; lpsz++)
    {
        /* moved from below for C syntax reason - 2024-09-01 DT */
        /* djb-rwth: return values needed for va_arg; djb-rwth: ignoring LLVM warnings: function returning value */
        int ivarg;
        double dvarg;
        void* ivvarg;
        int* ipvarg;

        /* handle '%' character, but watch out for '%%' */
        if (*lpsz != '%' || *(++lpsz) == '%')
        {
            nMaxLen += 1;
            continue;
        }

        nItemLen = 0;

        /*  handle '%' character with format */
        nWidth = 0;
        for (; *lpsz; lpsz++)
        {
            /* check for valid flags */
            if (*lpsz == '#')
            {
                nMaxLen += 2;   /* for '0x' */
            }
            else if (*lpsz == '*')
            {
                nWidth = va_arg(argList, int);
            }
            else if (*lpsz == '-' || *lpsz == '+' || *lpsz == '0'
                || *lpsz == ' ')
            {
                ;
            }
            else /* hit non-flag character */
            {
                break;
            }
        }
        /* get width and skip it */
        if (nWidth == 0)
        {
            /* width indicated by */
            nWidth = atoi(lpsz);
            for (; *lpsz && isdigit(*lpsz); lpsz++)
            {
                ;
            }
        }
        /*ASSERT(nWidth >= 0);*/
        if (nWidth < 0)
        {
            goto exit_error; /* instead of exception */
        }

        nPrecision = 0;
        if (*lpsz == '.')
        {
            /* skip past '.' separator (width.precision)*/
            lpsz++;

            /* get precision and skip it*/
            if (*lpsz == '*')
            {
                nPrecision = va_arg(argList, int);
                lpsz++;
            }
            else
            {
                nPrecision = atoi(lpsz);
                for (; *lpsz && isdigit(*lpsz); lpsz++)
                {
                    ;
                }
            }
            if (nPrecision < 0)
            {
                goto exit_error; /* instead of exception */
            }
        }

        /* should be on type modifier or specifier */
        nModifier = 0;
        switch (*lpsz)
        {
            /* modifiers that affect size */
        case 'h':
            switch (lpsz[1])
            {
            case 'd':
            case 'i':
            case 'o':
            case 'x':
            case 'X':
            case 'u':
                /* short unsigned, short double, etc. -- added to the original MS example */
                /* ignore the fact that these modifiers do affect size */
                lpsz++;
                break;
            default:
                nModifier = FORCE_ANSI;
                lpsz++;
                break;
            }
            break;
        case 'l':
            switch (lpsz[1])
            {
            case 'd':
            case 'i':
            case 'o':
            case 'x':
            case 'X':
            case 'u':
            case 'f': /* long float -- post ANSI C */
                /* long unsigned, long double, etc. -- added to the original MS example */
                /* ignore the fact that these modifiers do affect size */
                lpsz++;
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
            lpsz++;
            break;
        }

        /* now should be on specifier */

        switch (*lpsz | nModifier)
        {
            /* single characters*/
        case 'c':
        case 'C':
            nItemLen = 2;
            ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
            break;
        case 'c' | FORCE_ANSI:
        case 'C' | FORCE_ANSI:
            nItemLen = 2;
            ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
            break;
        case 'c' | FORCE_UNICODE:
        case 'C' | FORCE_UNICODE:
            goto exit_error;  /* no UNICODE, please */
            /*
            nItemLen = 2;
            va_arg(argList, wchar_t);
            break;
            */

            /* strings*/
        case 's':
        case 'S':
            nItemLen = (int)strlen(va_arg(argList, char*));
            nItemLen = inchi_max(1, nItemLen);
            break;
        case 's' | FORCE_ANSI:
        case 'S' | FORCE_ANSI:
            nItemLen = (int)strlen(va_arg(argList, char*));
            nItemLen = inchi_max(1, nItemLen);
            break;

        case 's' | FORCE_UNICODE:
        case 'S' | FORCE_UNICODE:
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
            {
                nItemLen = inchi_min(nItemLen, nPrecision);
            }
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
                ivarg = va_arg(argList, int); /* djb-rwth: int return value; ignoring LLVM warning */
                nItemLen = 32;
                nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
                break;

            case 'e':
            case 'f':
            case 'g':
            case 'G':
                dvarg = va_arg(argList, double); /* djb-rwth: double return value; ignoring LLVM warning */
                nItemLen = 32;
                nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
                break;

            case 'p':
                ivvarg = va_arg(argList, void*); /* djb-rwth: void* return value; ignoring LLVM warning */
                nItemLen = 32;
                nItemLen = inchi_max(nItemLen, nWidth + nPrecision);
                break;

                /* no output */
            case 'n':
                ipvarg = va_arg(argList, int*); /* djb-rwth: int* return value; ignoring LLVM warning */
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


/****************************************************************************
    Get at most n-1 chars, plus a null, then advance input's start.
    Return emulates fgets()
****************************************************************************/
char* inchi_sgets(char* s, int n, INCHI_IOSTREAM* ios)
{
    int c = 0;
    char* p;
    char* inp;

    inp = ios->s.pStr + ios->s.nPtr;

    if (n <= 0)
    {
        return NULL;
    }

    if (NULL == inp)
    {
        /* like EOF */
        return NULL; /* djb-rwth: addressing coverity ID #499480 -- inp can be NULL */
    }

    p = s;

    /*
    if ( *inp == '\0' )
        s = '\0';
    else
    */

    while (--n > 0 && (c = *inp++))
    {
        ios->s.nPtr++;
        if ((*p++ = c) == '\n')
        {
            break;
        }
    }
    *p = '\0';

    /* printf("\n*** {%-s}",s); */

    return (c == '\0' && p == s)
        ? NULL /* like EOF reached */
        : s;
}


/****************************************************************************
 Init expandable buffer of type INCHI_IOS_STRING
****************************************************************************/
int inchi_strbuf_init(INCHI_IOS_STRING* buf, int start_size, int incr_size)
{
    char* new_str = NULL;
    memset(buf, 0, sizeof(*buf)); /* djb-rwth: memset_s C11/Annex K variant? */

    if (start_size <= 0)
    {
        start_size = INCHI_STRBUF_INITIAL_SIZE;
    }
    if (incr_size <= 0)
    {
        incr_size = INCHI_STRBUF_SIZE_INCREMENT;
    }

    new_str = (char*)inchi_calloc(start_size, sizeof(char));

    if (!new_str)
    {
        return -1;
    }

    buf->pStr = new_str;
    buf->nAllocatedLength = start_size;
    buf->nPtr = incr_size;

    return start_size;
}


/****************************************************************************
    Reset INCHI_IOS_STRING object holding an expandable buffer string
    (place '\0' at the start and do _not_ free memory).
****************************************************************************/
void inchi_strbuf_reset(INCHI_IOS_STRING* buf)
{
    if (!buf)
    {
        return;
    }
    if (buf->pStr)
    {
        buf->pStr[0] = '\0';
    }

    buf->nUsedLength = buf->nPtr = 0;
}


/****************************************************************************
Close INCHI_IOS_STRING object holding an expandable buffer string,
    free previously allocated sring memory
****************************************************************************/
void inchi_strbuf_close(INCHI_IOS_STRING* buf)
{
    if (!buf)
    {
        return;
    }
    if (buf->pStr)
    {
        inchi_free(buf->pStr);
    }

    memset(buf, 0, sizeof(*buf)); /* djb-rwth: memset_s C11/Annex K variant? */
}


/****************************************************************************/
int inchi_strbuf_create_copy(INCHI_IOS_STRING* buf2, INCHI_IOS_STRING* buf)
{
    char* new_str = NULL;

    new_str = (char*)inchi_calloc(buf->nAllocatedLength, sizeof(char));
    buf2->pStr = new_str;
    if (!new_str)
    {
        return -1;
    }
    buf2->nAllocatedLength = buf->nAllocatedLength;
    buf2->nUsedLength = buf->nUsedLength;
    buf2->nPtr = buf->nPtr;

    return 0;
}


/****************************************************************************
 Check size and if necessary expand string buffer in INCHI_IOS_STRING
****************************************************************************/
int inchi_strbuf_update(INCHI_IOS_STRING* buf, int new_addition_size)
{
    int requsted_len;

    if (!buf)
    {
        return -1;
    }

    if (new_addition_size <= 0)
    {
        return buf->nAllocatedLength;
    }

    requsted_len = buf->nUsedLength + new_addition_size;

    if (requsted_len >= buf->nAllocatedLength)
    {
        /* Expand */
        int  nAddLength = inchi_max(buf->nPtr, new_addition_size);
        /* buf->nPtr stores size increment for this buffer */
        char* new_str =
            (char*)inchi_calloc((long long)buf->nAllocatedLength + (long long)nAddLength,
                sizeof(new_str[0])); /* djb-rwth: cast operators added */
        if (!new_str)
        {
            return -1; /* failed */
        }
        if (buf->pStr)
        {
            if (buf->nUsedLength > 0)
            {
                memcpy(new_str, buf->pStr, sizeof(new_str[0]) * buf->nUsedLength);
            }
            inchi_free(buf->pStr);
        }
        buf->pStr = new_str;
        buf->nAllocatedLength += nAddLength;
    }

    return buf->nAllocatedLength;
}


/****************************************************************************
 Add to the end of string in INCHI_IOS_STRING object,
    expanding buffer if necessary
****************************************************************************/
int inchi_strbuf_printf(INCHI_IOS_STRING* buf, const char* lpszFormat, ...)
{
    int ret = 0, max_len;
    va_list argList;

    if (!buf)
    {
        return -1;
    }

    my_va_start(argList, lpszFormat);
    max_len = GetMaxPrintfLength(lpszFormat, argList);
    va_end(argList);
    if (max_len < 0)
    {
        return 0;
    }

    inchi_strbuf_update(buf, max_len);

    my_va_start(argList, lpszFormat);
    ret = vsprintf(buf->pStr + buf->nUsedLength, lpszFormat, argList);
    va_end(argList);
    if (ret >= 0)
    {
        buf->nUsedLength += ret;
    }

    return ret;
}


/****************************************************************************
    Print to string in INCHI_IOS_STRING object
    from specified position 'npos', expanding buffer if necessary.
    NB: be careful, intentionally no checks on where is 'npos'!
****************************************************************************/
int inchi_strbuf_printf_from(INCHI_IOS_STRING* buf,
    int npos,
    const char* lpszFormat, ...)
{
    int ret = 0, max_len;
    va_list argList;

    if (!buf)
    {
        return -1;
    }

    my_va_start(argList, lpszFormat);
    max_len = GetMaxPrintfLength(lpszFormat, argList);
    va_end(argList);
    if (max_len < 0)
    {
        return 0;
    }

    max_len += npos;

    inchi_strbuf_update(buf, max_len);

    my_va_start(argList, lpszFormat);
    ret = vsprintf(buf->pStr + npos, lpszFormat, argList);
    va_end(argList);
    if (ret >= 0)
    {
        buf->nUsedLength = npos + ret;
    }

    return ret;
}


/****************************************************************************
    Reads the next line to growing str buf.
    Returns n of read chars, -1 at end of file or at error.
*****************************************************************************/
int inchi_strbuf_getline(INCHI_IOS_STRING* buf,
    FILE* f,
    int crlf2lf,
    int preserve_lf)
{
    int c;
    inchi_strbuf_reset(buf);

    while (1)
    {
        c = fgetc(f);
        if (ferror(f))
        {
            return -1;
        }
        if (c == EOF)
        {
            return -1;
        }
        inchi_strbuf_printf(buf, "%c", c);
        if (c == '\n')
        {
            break;
        }
    }

    if (crlf2lf)
    {
        if (buf->nUsedLength > 2)
        {
            if (buf->pStr[buf->nUsedLength - 2] == '\r')
            {
                buf->pStr[buf->nUsedLength - 2] = '\n';
                buf->pStr[--buf->nUsedLength] = '\0';
            }
        }
    }
    if (!preserve_lf)
    {
        buf->pStr[--buf->nUsedLength] = '\0';
    }

    return buf->nUsedLength;
}



/****************************************************************************
    Adds the next line to growing str buf (does not reset buf before adding).
    Returns n of read chars, -1 at end of file or at error.
****************************************************************************/
int inchi_strbuf_addline(INCHI_IOS_STRING* buf,
    INCHI_IOSTREAM* inp_stream,
    int crlf2lf,
    int preserve_lf)
{
    int c;

    while (1)
    {
        c = inchi_ios_str_getc(inp_stream);
        if (c == EOF)
        {
            return -1;
        }
        inchi_strbuf_printf(buf, "%c", c);
        if (c == '\n')
        {
            break;
        }
    }
    if (crlf2lf)
    {
        if (buf->nUsedLength > 2)
        {
            if (buf->pStr[buf->nUsedLength - 2] == '\r')
            {
                buf->pStr[buf->nUsedLength - 2] = '\n';
                buf->pStr[--buf->nUsedLength] = '\0';
            }
        }
    }
    if (!preserve_lf)
    {
        buf->pStr[--buf->nUsedLength] = '\0';
    }

    return buf->nUsedLength;
}


#if ( defined(_WIN32) && defined(_DEBUG) && defined(_MSC_VER) )
#include <windows.h>


/****************************************************************************/
/* djb-rwth: placed as a global variable to avoid function buffer issues */
char it_buffer[32767];
int _inchi_trace(char* format, ...)
{
    /*
    TCHAR buffer[32767];
    char buffer[32767];
    */
    int ret;

    va_list argptr;
    va_start(argptr, format);
    /*wvsprintf(buffer, format, argptr);*/
    ret = vsprintf(it_buffer, format, argptr);
    va_end(argptr);
    OutputDebugString(it_buffer);
    return 1;
}
#else
int _inchi_trace(char* format, ...)
{
    return 1;
}
#endif


/****************************************************************************
Output structure (compound) header for current record
****************************************************************************/
int Output_RecordInfo(INCHI_IOSTREAM* out_file,
    int              num_input_struct,
    int              bNoStructLabels,
    const char* szSdfLabel,
    const char* szSdfValue,
    unsigned long    lSdfId,
    char* pLF,
    char* pTAB)
{
    if (bNoStructLabels)
    {
        return 0;
    }

#ifdef TARGET_LIB_FOR_WINCHI
    /*out_file->s.nUsedLength = 0;
    inchi_ios_reset( out_file );*/
#endif
    if (!(szSdfLabel && szSdfLabel[0]) && !(szSdfValue && szSdfValue[0]))
    {
        inchi_ios_print_nodisplay(out_file, "%sStructure: %d", pLF, num_input_struct);
        inchi_ios_print_nodisplay(out_file, "%s", pTAB);
    }
    else
    {
        inchi_ios_print_nodisplay(out_file, "%sStructure: %d.%s%s%s%s",
            pLF, num_input_struct,
            SDF_LBL_VAL(szSdfLabel, szSdfValue));

        if (lSdfId)
        {
            (out_file->s.nUsedLength)--;
            inchi_ios_print_nodisplay(out_file, ":%lu", lSdfId);
        }
        inchi_ios_print_nodisplay(out_file, "%s", pTAB);
    }

    return 0;
}
