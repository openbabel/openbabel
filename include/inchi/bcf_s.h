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

/* djb-rwth: bounds-checking functions / C11 Annex K */

/* djb-rwth: adding guard idiom instead of #pragma once for code portabilty */
#ifndef BCF_S_H 
#define BCF_S_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define CURRENT_VER "1.07.5" /* djb-rwth: full version number in the output, required for fixing GHI #61 */

#ifdef __STDC_LIB_EXT1__
#if (__STDC_LIB_EXT1__ >= 201112L)
#define USE_BCF_GCC 1 /* djb-rwth: use bounds-checking functions / C11 Annex K */
#define __STDC_WANT_LIB_EXT1__ 1 /* Want the ext1 functions */
#endif
#endif

#ifdef __STDC_SECURE_LIB__
#define USE_BCF_MS 1 /* djb-rwth: use bounds-checking functions / C11 Annex K */
#endif

#if (USE_BCF_GCC || USE_BCF_MS)
/* djb-rwth: change the following parameter to 1 (or other unsigned integer) to use bounds-checking functions */
#define USE_BCF 0
#endif

/* djb-rwth: control macros */
#define RINCHI_TEST 0 /* djb-rwrh: RInChI testing */
/* #define GHI100_FIX */ /* djb-rwth: enabling any of the three GHI #100 fixes defined by SPRINTF_FLAG macro */
#define SPRINTF_FLAG 2 /* djb-rwth: 1 - stb_(v)s(n)pritnf, 2 - custom double to 2*int; any other value (different from 1 or 2) - standard sprintf w/ setlocale */

/* djb-rwth: custom functions */
int max_3(int a, int b, int c);
int memcpy_custom(char** dst, char* src, unsigned long long len);
int dbl2int(char* str, int fwidth, int ndecpl, char dbl_flag, double dblinp);

#endif