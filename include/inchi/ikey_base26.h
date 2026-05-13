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


#ifndef _IKEY_BASE26_H_
#define _IKEY_BASE26_H_

/*
    Base-26 encoding procedures.

    'Base26' characters here are considered to be uppercase English
    letters 'A..Z'
*/


/* Uncomment the next line to fix base-26 encoding bug */
/*#define FIX_BASE26_ENC_BUG 1*/

typedef unsigned int UINT32;
typedef unsigned short int UINT16;

#ifdef __cplusplus
extern "C" {
#endif

/*    Get a character representing 1st 14-bit triplet
    (bits 0..13 of contiguous array of octets)        */
    const char* base26_triplet_1( const unsigned char *a );
    /*    Get a character representing 2nd 14-bit triplet (bits 14..27)    */
    const char* base26_triplet_2( const unsigned char *a );
    /*    Get a character representing 3rd 14-bit triplet (bits 28..41)    */
    const char* base26_triplet_3( const unsigned char *a );
    /*    Get a character representing 4th 14-bit triplet (bits 42..55)    */
    const char* base26_triplet_4( const unsigned char *a );

    /*
        Tail dublets
    */

    /*    Get dublet (bits 28..36)    */
    const char* base26_dublet_for_bits_28_to_36( unsigned char *a );
    /*    Get dublet (bits 56..64)    */
    const char* base26_dublet_for_bits_56_to_64( unsigned char *a );
    /*    Get hash extension in hexadecimal representation for the major block.
        Len(extension) = 256 - 65 = 191 bit.                                */
    void get_xtra_hash_major_hex( const unsigned char *a, char* szXtra );
    /*    Get hash extension in hexadecimal representation for the minor block.
        Len(extension) = 256 - 37 = 219 bit.                                */
    void get_xtra_hash_minor_hex( const unsigned char *a, char* szXtra );

    /*    Used instead of isupper() to avoid locale interference.    */
#define isbase26(_c)    ( ((unsigned)(_c) >= 'A') && ((unsigned)(_c) <= 'Z') )

#ifdef __cplusplus
}
#endif


#endif    /* _IKEY_BASE26_H_ */
