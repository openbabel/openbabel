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


#ifndef _ICHI_IO_H_
#define _ICHI_IO_H_


/*
    INCHI I/O (IOSTREAM MAINLY) OPERATIONS
*/


#include "mode.h"
#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


    void inchi_ios_init( INCHI_IOSTREAM *ios, int io_type, FILE *f );
    void inchi_ios_flush( INCHI_IOSTREAM *ios );
    void inchi_ios_flush2( INCHI_IOSTREAM *ios, FILE *f2 );
    void inchi_ios_close( INCHI_IOSTREAM *ios );
    void inchi_ios_reset( INCHI_IOSTREAM *ios );
    void inchi_ios_free_str( INCHI_IOSTREAM *ios );
    int inchi_ios_create_copy( INCHI_IOSTREAM* ios, INCHI_IOSTREAM* ios0 );

    int inchi_ios_gets( char *szLine, int len, INCHI_IOSTREAM *ios, int *bTooLongLine );
    int inchi_ios_getsTab( char *szLine, int len, INCHI_IOSTREAM *ios, int *bTooLongLine );
    int inchi_ios_getsTab1( char *szLine, int len, INCHI_IOSTREAM *ios, int *bTooLongLine );

    int inchi_ios_print( INCHI_IOSTREAM *ios, const char* lpszFormat, ... );
    int inchi_ios_print_nodisplay( INCHI_IOSTREAM *ios, const char* lpszFormat, ... );
    int inchi_ios_flush_not_displayed( INCHI_IOSTREAM * ios );
    int push_to_winchi_text_window( INCHI_IOSTREAM * ios ); /*, const char* lpszFormat, ... ) */

    /* Print to string buffer or to file+stderr */
    int inchi_ios_eprint( INCHI_IOSTREAM *ios, const char* lpszFormat, ... );


    /*
        PLAIN FILE OPERATIONS
    */

    /* Print to file, echoing to stderr */
    int inchi_fprintf( FILE* f, const char* lpszFormat, ... );
    int inchi_print_nodisplay( FILE* f, const char* lpszFormat, ... );

    char* inchi_fgetsLf( char* line, int line_len, INCHI_IOSTREAM* inp_stream );
    int inchi_fgetsLfTab( char *szLine, int len, FILE *f );

    char *inchi_sgets( char *s, int n, INCHI_IOSTREAM* ios );




    /*
        Support of simplistic print buffer (growing string)
    */

#define INCHI_STRBUF_INITIAL_SIZE 262144
#define INCHI_STRBUF_SIZE_INCREMENT 262144
#define INCHI_STRBUF_SMALLER_INITIAL_SIZE 1024
#define INCHI_STRBUF_SMALLER_SIZE_INCREMENT 4096

int inchi_strbuf_init( INCHI_IOS_STRING *buf, int start_size, int incr_size );
void inchi_strbuf_reset( INCHI_IOS_STRING *buf );
void inchi_strbuf_close( INCHI_IOS_STRING *buf );
int  inchi_strbuf_printf( INCHI_IOS_STRING *buf, const char* lpszFormat, ... );
int  inchi_strbuf_printf_from( INCHI_IOS_STRING *buf, int npos, const char* lpszFormat, ... );
int inchi_strbuf_create_copy( INCHI_IOS_STRING *buf2, INCHI_IOS_STRING *buf );
int  inchi_strbuf_update( INCHI_IOS_STRING *buf, int new_addition_size );
int inchi_strbuf_getline( INCHI_IOS_STRING *buf, FILE *f, int crlf2lf, int preserve_lf );

int inchi_strbuf_addline( INCHI_IOS_STRING *buf, INCHI_IOSTREAM *inp_stream, int crlf2lf, int preserve_lf );

void print_sequence_of_nums_compressing_ranges( int n, int *num, INCHI_IOS_STRING *strbuf );

int _inchi_trace( char *format, ... );

int Output_RecordInfo( INCHI_IOSTREAM *out_file,
                        int num_input_struct,
                        int bNoStructLabels,
                        const char *szSdfLabel,
                        const char *szSdfValue,
                        unsigned long lSdfId,
                        char *pLF,
                        char *pTAB);


#if ( defined(_WIN32) && defined(_DEBUG) && defined(_MSC_VER) )
#define ITRACE_ _inchi_trace
#else
#define ITRACE_ 0 && _inchi_trace
#endif

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif    /* _ICHI_IO_H_ */
