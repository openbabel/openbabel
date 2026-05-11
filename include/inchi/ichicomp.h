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


#ifndef _INCHI_COMP_H_
#define _INCHI_COMP_H_

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

    /* compatibility */

#if ( defined(__GNUC__) && defined(__MINGW32__) && __GNUC__ == 3 && __GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ == 0 && defined(_WIN32) )
/* replace with the proper definition for GNU gcc & MinGW-2.0.0-3 (mingw special 20020817-1), gcc 3.2 core */
#define my_va_start(A,B) ((A)=(va_list)__builtin_next_arg(lpszFormat))
#else
#define my_va_start va_start
#endif



/*  ANSI redefinitions */
#ifdef COMPILE_ANSI_ONLY  /* { */
#ifndef __isascii
#define __isascii(val)  ((unsigned)(val) <= 0x7F)
#endif

/* #ifndef __GNUC__ */
/* these non-ANSI functions are implemented in gcc */
#include <stdio.h>
 /* this #include provides size_t definition */
 /* implementation is located in util.c */
/*#if ( !defined(_MSC_VER) || defined(__STDC__) && __STDC__ == 1 )*/

#if ( defined(COMPILE_ADD_NON_ANSI_FUNCTIONS) || defined(__STDC__) && __STDC__ == 1 )

/* support (VC++ Language extensions) = OFF && defined(COMPILE_ANSI_ONLY) */
    int   inchi_memicmp( const void*, const void*, size_t );
    int   inchi_stricmp( const char *s1, const char *s2 );
    char *inchi__strnset( char *string, int c, size_t count );
    char *inchi__strdup( const char *string );
#endif
/* #endif */

#endif /* } */


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif    /* _INCHI_COMP_H_ */
