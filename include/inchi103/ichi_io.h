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




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


                                    INCHI_IOSTREAM OPERATIONS 


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


void inchi_ios_init(INCHI_IOSTREAM *ios, int io_type, FILE *f);
void inchi_ios_flush(INCHI_IOSTREAM *ios);
void inchi_ios_flush2(INCHI_IOSTREAM *ios, FILE *f2);
void inchi_ios_close(INCHI_IOSTREAM *ios);
void inchi_ios_reset(INCHI_IOSTREAM *ios);

int inchi_ios_gets( char *szLine, int len, INCHI_IOSTREAM *ios, int *bTooLongLine );
int inchi_ios_getsTab( char *szLine, int len, INCHI_IOSTREAM *ios, int *bTooLongLine );
int inchi_ios_getsTab1( char *szLine, int len, INCHI_IOSTREAM *ios, int *bTooLongLine );

int inchi_ios_print( INCHI_IOSTREAM *ios, const char* lpszFormat, ... );
int inchi_ios_print_nodisplay( INCHI_IOSTREAM *ios, const char* lpszFormat, ... );

/* Print to string buffer or to file+stderr */
int inchi_ios_eprint( INCHI_IOSTREAM *ios, const char* lpszFormat, ... );


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


                                    PLAIN FILE OPERATIONS 


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/* Print to file, echoing to stderr */
int inchi_fprintf( FILE* f, const char* lpszFormat, ... );
int inchi_print_nodisplay( FILE* f, const char* lpszFormat, ... );

char* inchi_fgetsLf( char* line, int line_len, FILE* inp );
int inchi_fgetsLfTab( char *szLine, int len, FILE *f );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif




