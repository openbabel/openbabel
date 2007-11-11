/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02-beta
 * August 23, 2007
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Former COMDEF.H
Renamed 06/12/07 to avoid occassional conflict with Microsoft's COMDEF.H
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/* common definitions -- do not change */
#ifndef __INCOMDEF_H__
#define __INCOMDEF_H__

#include "ichisize.h"


/* SDF treatment */
#define MAX_SDF_HEADER        64 /* max length of the SDFile data header */
#define MAX_SDF_VALUE        255 /* max lenght of the SDFile data value */

/* size resrictions */
#define ATOM_EL_LEN           6  /* length of atom name string including zero termination */
#define ATOM_INFO_LEN         36 /* inf_ATOM output string ^123Al^+2H12..(+)/999/999/999/999: 32 chars */
#define MAXVAL                20 /* max number of bonds per atom */
#define MAX_STEREO_BONDS      3  /* max number of stereogenic bonds per atom */
#define NUM_H_ISOTOPES        3  /* number of hydrogen isotopes: protium, deuterium, tritium */
#define ATW_H                 1  /* hydrogen atomic weight */

/* input bond type definition */
#define MIN_INPUT_BOND_TYPE 1
#define MAX_INPUT_BOND_TYPE 4

#define BOND_TYPE_SINGLE    1
#define BOND_TYPE_DOUBLE    2
#define BOND_TYPE_TRIPLE    3
#define BOND_TYPE_ALTERN    4

#define STEREO_SNGL_UP       1
#define STEREO_SNGL_EITHER   4
#define STEREO_SNGL_DOWN     6
#define STEREO_DBLE_EITHER   3


/* MOlfile */
#define INPUT_STEREO_SNGL_UP       1
#define INPUT_STEREO_SNGL_EITHER   4
#define INPUT_STEREO_SNGL_DOWN     6
#define INPUT_STEREO_DBLE_EITHER   3

/*
#define BOND_MARK_ODD      0x10
#define BOND_MARK_EVEN     0x20
*/
#define BOND_MARK_PARITY    0x30
#define BOND_MARK_HIGHLIGHT 0x40  /* highlight equivalent components */

#define BOND_MARK_ODD      '-'
#define BOND_MARK_EVEN     '+'
#define BOND_MARK_UNDF     '?'
#define BOND_MARK_UNKN     'u'
#define BOND_MARK_ERR      '*'

#define SALT_DONOR_H      1
#define SALT_DONOR_Neg    2
#define SALT_ACCEPTOR     4
#define SALT_p_DONOR      8  /* >C-SH   */
#define SALT_p_ACCEPTOR  16  /* >C-S(-) */
#define SALT_DONOR_ALL     (SALT_DONOR_Neg | SALT_DONOR_H | SALT_p_ACCEPTOR | SALT_p_DONOR)
#define SALT_DONOR_Neg2    (SALT_DONOR_Neg | SALT_p_ACCEPTOR)
#define SALT_DONOR_H2      (SALT_DONOR_H   | SALT_p_DONOR)
#define SALT_DONOR         (SALT_DONOR_Neg | SALT_DONOR_H)

#define SALT_SELECTED    32

/* radical definitions */
#define RADICAL_SINGLET 1
#define RADICAL_DOUBLET 2
#define RADICAL_TRIPLET 3

/* metal definition */
#define METAL           1          /* definition of an element: lowest valence */
#define METAL2          3          /* definition of an element: lowest and next to it valence */
#define IS_METAL        3          /* metal bitmap */
/* isotopic shift */
#define ZERO_ATW_DIFF        127   /* mark mass of the most abundant isotope */

/* other types */

#define UCINT  (int)(unsigned char)

#ifndef INCHI_US_CHAR_DEF
typedef signed char   S_CHAR;
typedef unsigned char U_CHAR;
#define INCHI_US_CHAR_DEF
#endif

#ifndef INCHI_US_SHORT_DEF
typedef signed short   S_SHORT;
typedef unsigned short U_SHORT;
#define INCHI_US_SHORT_DEF
#endif

/* BILLY 8/6/04 */
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

#define STR_ERR_LEN 256
int AddMOLfileError( char *pStrErr, const char *szMsg );

/* allocator */
#ifndef inchi_malloc
void *inchi_malloc(size_t c);
#endif
#ifndef inchi_calloc
void *inchi_calloc(size_t c, size_t n);
#endif
#ifndef inchi_free
void inchi_free(void *p);
#endif

/* output */
int my_fprintf( INCHI_FILE* f, const char* lpszFormat, ... );
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
int my_fileprintf( FILE* f, const char* lpszFormat, ... );


int inchi_print( INCHI_FILE* f, const char* lpszFormat, ... );
int inchi_print_nodisplay( INCHI_FILE* f, const char* lpszFormat, ... );
/* sorting etc */
void swap ( char *a, char *b, size_t width );
int insertions_sort( void *base, size_t num, size_t width, int ( *compare )(const void *e1, const void *e2 ) );
int insertions_sort_AT_NUMBERS( AT_NUMB *base, int num, int ( *compare )(const void *e1, const void *e2 ) );


#define MOLFILE_ERR_FIN(err, new_err, err_fin, msg) \
        if ( !(err) && (new_err) ) { (err) = (new_err);} AddMOLfileError(pStrErr, (msg)); goto err_fin
#define MOLFILE_ERR_SET(err, new_err, msg) \
        if ( !(err) && (new_err) ) { (err) = (new_err);} AddMOLfileError(pStrErr, (msg))


/* BILLY 8/6/04 */
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif  /* __INCOMDEF_H__ */

