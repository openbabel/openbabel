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


#ifndef _INCOMDEF_H_
#define _INCOMDEF_H_


#include "ichisize.h"
#include "mode.h"

/* Common definitions -- do not change */

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

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


/* allocator */
#ifndef inchi_malloc
    void *inchi_malloc( size_t c );
#endif
#ifndef inchi_calloc
    void *inchi_calloc( size_t c, size_t n );
#endif
#ifndef inchi_free
    void inchi_free( void *p );
#endif



/* sorting etc */

    void inchi_swap( char *a, char *b, size_t width );

    int insertions_sort( void *pCG,
                         void *base, size_t num, size_t width, int( *compare )( const void *e1, const void *e2, void * ) );
    int insertions_sort_AT_NUMBERS( void *pCG,
                                    AT_NUMB *base, int num, int( *compare )( const void *e1, const void *e2, void * ) );
    /*
    int insertions_sort( void *base, size_t num, size_t width, int ( *compare )(const void *e1, const void *e2 ) );
    int insertions_sort_AT_NUMBERS( AT_NUMB *base, int num, int ( *compare )(const void *e1, const void *e2 ) );
    */


    /* min-max */

#define inchi_max(a,b)  (((a)>(b))?(a):(b))
#define inchi_min(a,b)  (((a)<(b))?(a):(b))

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif    /* _INCOMDEF_H_ */
