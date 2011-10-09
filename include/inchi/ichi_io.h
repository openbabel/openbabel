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




/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


                                    INCHI_IOSTREAM OPERATIONS 


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
#ifndef COMPILE_ALL_CPP
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

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif




