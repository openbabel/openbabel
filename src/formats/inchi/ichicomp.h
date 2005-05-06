/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.00
 * April 13, 2005
 * Developed at NIST
 */

#ifndef __INCHI_COMPAT_H__
#define __INCHI_COMPAT_H__

/* compatibility */

#if( defined(__GNUC__) && defined(__MINGW32__) && __GNUC__ == 3 && __GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ == 0 && defined(_WIN32) ) 
/* replace with the proper definition for GNU gcc & MinGW-2.0.0-3 (mingw special 20020817-1), gcc 3.2 core */
#define my_va_start(A,B) ((A)=(va_list)__builtin_next_arg(lpszFormat))
#else
#define my_va_start va_start
#endif



/*  ANSI redefinitions */
#ifdef INCHI_ANSI_ONLY  /* { */
#ifndef __isascii
#define __isascii(val)  ((unsigned)(val) <= 0x7F)
#endif

/* #ifndef __GNUC__ */
/* these non-ANSI functions are implemented in gcc */
#include <stdio.h>
 /* this #include provides size_t definition */
 /* implementation is located in util.c */
#if ( !defined(_MSC_VER) || defined(__STDC__) && __STDC__ == 1 )
/* support (VC++ Language extensions) = OFF && defined(INCHI_ANSI_ONLY) */
int   memicmp (const void*, const void*, size_t);
int   stricmp( const char *s1, const char *s2 );
char *_strnset( char *string, int c, size_t count );
char *_strdup( const char *string );
#endif
/* #endif */

#endif /* } */

#define inchi_max(a,b)  (((a)>(b))?(a):(b))
#define inchi_min(a,b)  (((a)<(b))?(a):(b))


#endif /* __INCHI_COMPAT_H__ */
