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

#include "mode.h"

#if defined(COMPILE_ANSI_ONLY) && defined(__APPLE__)
/*    For build under OSX, advice from Burt Leland */
#include "ichicomp.h"    /* Needed for __isascii define */
#endif

/* djb-rwth: defining __isascii */
#if defined(__isascii)
#define is_ascii __isascii
#elif defined(isascii)
#define is_ascii isascii
#else
#define is_ascii(c)   ((unsigned)(c) < 0x80)
#endif

#include "util.h"
#include "extr_ct.h"

#include "bcf_s.h"

#define MIN_ATOM_CHARGE        (-2)
#define MAX_ATOM_CHARGE         2
#define NEUTRAL_STATE          (-MIN_ATOM_CHARGE)
#define NUM_ATOM_CHARGES       (MAX_ATOM_CHARGE - MIN_ATOM_CHARGE + 1)
#define MAX_NUM_VALENCES        5                /* max. number + 1 to provide zero termination */



/*
    Local
*/

/*    CHEMICAL ELEMENTS & ATOMIC VALENCE MODEL
      FOR VARIOUS OXIDATION STATES
*/
typedef struct tagElData
{
    /*    Element chemical symbol */
    const char *szElName;
    /*    Average atomic mass from the Periodic Chart of the Elements
        (Fisher cat. no. 05-702-10) */
    int     nAtMass;
    /*    (not used currently) Atomic mass of the most abundant isotope */
    int     nNormAtMass;
    /*    (not used currently) Exact mw of the most abundant isotope (not used)    */
    double  dAtMass;
    /*    METAL or METAL2 */
    int     nType;
    /*    (not used currently) Pauling electronegativity x 10; 0 means unknown    */
    int     nElNegPauling10;
    /*    InChI does not add implicit H to atoms that have non-zero bSkipAddingH */
    /*    NB: was called bDoNotAddH, renamed to avoid confusion with other procedures */
    int     bSkipAddingH;
    S_CHAR  cValence[NUM_ATOM_CHARGES][MAX_NUM_VALENCES];
} ELDATA;


const ELDATA ElData[] =
{
/*       avg  norm                      El    No  -------- Valence(s) of an ion or neutral atom -------------*/
/*        mw  mass  exact mw     type   neg   H   -2          -1          0           +1         +2          */
{ "H",    1,   1,   1.007825035,     0 , 21,  0, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "D",    2,   2,   2.014101778,     0 , 21,  0, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "T",    3,   3,   3.016049268,     0 , 21,  0, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "He",   4,   4,   4.002600000,     0 ,  0,  0, {{0,},       {0,},       {0,},       {0,},       {0,}       }},
{ "Li",   7,   7,   7.016000000, METAL , 10,  0, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Be",   9,   9,   9.012180000, METAL , 15,  0, {{0,},       {0,},       {2,},       {1,},       {0,}       }},
{ "B",   11,  11,  11.009300000,     0 , 20,  0, {{3,},       {4,},       {3,},       {2,},       {1,}       }},
{ "C",   12,  12,  12.000000000,     0 , 25,  0, {{2,},       {3,},       {4,},       {3,},       {2,}       }},
{ "N",   14,  14,  14.003074000,     0 , 30,  0, {{1,},       {2,},       {3,5},      {4,},       {3,}       }},
{ "O",   16,  16,  15.994914630,     0 , 35,  0, {{0,},       {1,},       {2,},       {3,5,},     {4,}       }},
{ "F",   19,  19,  18.998403220,     0 , 40,  0, {{0,},       {0,},       {1,},       {2,},       {3,5}      }},
{ "Ne",  20,  20,  19.992440000,     0 ,  0,  0, {{0,},       {0,},       {0,},       {0,},       {0,}       }},
{ "Na",  23,  23,  22.989770000, METAL ,  9,  0, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Mg",  24,  24,  23.985000000, METAL , 12,  0, {{0,},       {0,},       {2,},       {1,},       {0,}       }},
{ "Al",  27,  27,  26.981540000, METAL , 15,  0, {{3,5,},     {4,},       {3,},       {2,},       {1,}       }},
{ "Si",  28,  28,  27.976927100,     0 , 18,  0, {{2,},       {3,5},      {4,},       {3,},       {2,}       }},
{ "P",   31,  31,  30.973762000,     0 , 21,  0, {{1,3,5,7,}, {2,4,6,},   {3,5,},     {4,},       {3,}       }},
{ "S",   32,  32,  31.972070700,     0 , 25,  0, {{0,},       {1,3,5,7,}, {2,4,6},    {3,5,},     {4,}       }},
{ "Cl",  35,  35,  34.968852730,     0 , 30,  0, {{0,},       {0,},       {1,3,5,7},  {2,4,6},    {3,5,}     }},
{ "Ar",  40,  40,  39.962400000,     0 ,  0,  0, {{0,},       {0,},       {0,},       {0,},       {0,}       }},
{ "K",   39,  39,  38.963700000, METAL ,  8,  0, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Ca",  40,  40,  39.962600000, METAL , 10,  0, {{0,},       {0,},       {2,},       {1,},       {0,}       }},
{ "Sc",  45,  45,  44.955910000, METAL , 13,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Ti",  48,  48,  47.947950000, METAL , 15,  1, {{0,},       {0,},       {3,4},      {0,},       {0,}       }},
{ "V",   51,  51,  50.943960000, METAL , 16,  1, {{0,},       {0,},       {2,3,4,5,}, {0,},       {0,}       }},
{ "Cr",  52,  52,  51.940500000, METAL , 16,  1, {{0,},       {0,},       {2,3,6,},   {0,},       {0,}       }},
{ "Mn",  55,  55,  54.938050000, METAL2, 15,  1, {{0,},       {0,},       {2,3,4,6,}, {0,},       {0,}       }},
{ "Fe",  56,  56,  55.934900000, METAL2, 18,  1, {{0,},       {0,},       {2,3,4,6,}, {0,},       {0,}       }},
{ "Co",  59,  59,  58.933200000, METAL2, 18,  1, {{0,},       {0,},       {2,3,},     {0,},       {0,}       }},
{ "Ni",  59,  58,  57.935300000, METAL2, 18,  1, {{0,},       {0,},       {2,3,},     {0,},       {0,}       }},
{ "Cu",  64,  63,  62.929600000, METAL , 19,  1, {{0,},       {0,},       {1,2,},     {0,},       {0,}       }},
{ "Zn",  65,  64,  63.929147000, METAL , 16,  1, {{0,},       {0,},       {2,},       {0,},       {0,}       }},
{ "Ga",  70,  69,  68.925600000, METAL , 18,  0, {{3,5,},     {4,},       {3,},       {0,},       {1,}       }},
{ "Ge",  73,  74,  73.921177400,     0 , 18,  0, {{2,4,6,},   {3,5,},     {4,},       {3,},       {0,}       }},
{ "As",  75,  75,  74.921594200,     0 , 20,  0, {{1,3,5,7,}, {2,4,6,},   {3,5,},     {4,},       {3,}       }},
{ "Se",  79,  80,  79.916519600,     0 , 24,  0, {{0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     {4,}       }},
{ "Br",  80,  79,  78.918336100,     0 , 28,  0, {{0,},       {0,},       {1,3,5,7,}, {2,4,6,},   {3,5,}     }},
{ "Kr",  84,  84,  83.911500000,     0 ,  0,  0, {{0,},       {0,},       {0,},       {0,},       {0,}       }},
{ "Rb",  85,  85,  84.911800000, METAL ,  8,  0, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Sr",  88,  88,  87.905600000, METAL , 10,  0, {{0,},       {0,},       {2,},       {1,},       {0,}       }},
{ "Y",   89,  89,  88.905860000, METAL , 12,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Zr",  91,  90,  89.904700000, METAL , 14,  1, {{0,},       {0,},       {4,},       {0,},       {0,}       }},
{ "Nb",  93,  93,  92.906400000, METAL , 16,  1, {{0,},       {0,},       {3,5,},     {0,},       {0,}       }},
{ "Mo",  96,  98,  97.905400000, METAL , 18,  1, {{0,},       {0,},       {3,4,5,6,}, {0,},       {0,}       }},
{ "Tc",  98,  98,  97.907200000, METAL , 19,  1, {{0,},       {0,},       {7,},       {0,},       {0,}       }},
{ "Ru", 101, 102, 101.904300000, METAL , 22,  1, {{0,},       {0,},       {2,3,4,6,}, {0,},       {0,}       }},
{ "Rh", 103, 103, 102.905500000, METAL , 22,  1, {{0,},       {0,},       {2,3,4,},   {0,},       {0,}       }},
{ "Pd", 106, 106, 105.903500000, METAL , 22,  1, {{0,},       {0,},       {2,4,},     {0,},       {0,}       }},
{ "Ag", 108, 107, 106.905100000, METAL , 19,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Cd", 112, 114, 113.903400000, METAL , 17,  1, {{0,},       {0,},       {2,},       {0,},       {0,}       }},
{ "In", 115, 115, 114.903900000, METAL , 17,  0, {{3,5,},     {2,4,},     {3,},       {0,},       {1,}       }},
{ "Sn", 119, 120, 119.902200000, METAL2, 18,  0, {{2,4,6,},   {3,5},      {2,4,},     {3,},       {0,}       }},
{ "Sb", 122, 121, 120.903800000, METAL,  19,  0, {{1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,},     {3,}       }},
{ "Te", 128, 130, 129.906200000,     0 , 21,  0, {{0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,}     }},
{ "I",  127, 127, 126.904500000,     0 , 25,  0, {{0,},       {0,},       {1,3,5,7,}, {2,4,6},    {3,5,}     }},
{ "Xe", 131, 132, 131.904100000,     0 ,  0,  0, {{0,},       {0,},       {0,},       {0,},       {0,}       }},
{ "Cs", 133, 133, 132.905430000, METAL ,  7,  0, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Ba", 137, 138, 137.905200000, METAL ,  9,  0, {{0,},       {0,},       {2,},       {1,},       {0,}       }},
{ "La", 139, 139, 138.906360000, METAL , 11,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Ce", 140, 140, 139.905400000, METAL2,  0,  1, {{0,},       {0,},       {3,4,},     {0,},       {0,}       }},
{ "Pr", 141, 141, 140.907660000, METAL2,  0,  1, {{0,},       {0,},       {3,4,},     {0,},       {0,}       }},
{ "Nd", 144, 142, 141.907719000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Pm", 145, 145, 144.912800000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Sm", 150, 152, 151.919700000, METAL2,  0,  1, {{0,},       {0,},       {2,3,},     {0,},       {0,}       }},
{ "Eu", 152, 153, 152.921200000, METAL2,  0,  1, {{0,},       {0,},       {2,3,},     {0,},       {0,}       }},
{ "Gd", 157, 158, 157.924099000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Tb", 159, 159, 158.925350000, METAL2,  0,  1, {{0,},       {0,},       {3,4,},     {0,},       {0,}       }},
{ "Dy", 163, 164, 163.929200000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }}, /*  mw rounding uncertain */
{ "Ho", 165, 165, 164.930300000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Er", 167, 166, 165.930300000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Tm", 169, 169, 168.934230000, METAL2,  0,  1, {{0,},       {0,},       {2,3,},     {0,},       {0,}       }},
{ "Yb", 173, 174, 173.938900000, METAL2,  0,  1, {{0,},       {0,},       {2,3,},     {0,},       {0,}       }},
{ "Lu", 175, 175, 174.940800000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Hf", 178, 180, 179.946600000, METAL , 13,  1, {{0,},       {0,},       {4,},       {0,},       {0,}       }},
{ "Ta", 181, 181, 180.948010000, METAL , 15,  1, {{0,},       {0,},       {5,},       {0,},       {0,}       }},
{ "W",  184, 184, 183.951000000, METAL2, 17,  1, {{0,},       {0,},       {3,4,5,6,}, {0,},       {0,}       }},
{ "Re", 186, 187, 186.955800000, METAL2, 19,  1, {{0,},       {0,},       {2,4,6,7,}, {0,},       {0,}       }},
{ "Os", 190, 192, 191.961500000, METAL2, 22,  1, {{0,},       {0,},       {2,3,4,6,}, {0,},       {0,}       }},
{ "Ir", 192, 193, 192.962900000, METAL2, 22,  1, {{0,},       {0,},       {2,3,4,6,}, {0,},       {0,}       }},
{ "Pt", 195, 195, 194.964800000, METAL2, 22,  1, {{0,},       {0,},       {2,4,},     {0,},       {0,}       }},
{ "Au", 197, 197, 196.966560000, METAL , 24,  1, {{0,},       {0,},       {1,3,},     {0,},       {0,}       }},
{ "Hg", 201, 202, 201.970617000, METAL2, 19,  1, {{0,},       {0,},       {1,2,},     {0,},       {0,}       }},
{ "Tl", 204, 205, 204.974400000, METAL2, 18,  0, {{3,5,},     {2,4,},     {1,3,},     {0,},       {0,}       }},
{ "Pb", 207, 208, 207.976627000, METAL2, 18,  0, {{2,4,6,},   {3,5},      {2,4,},     {3,},       {0,}       }},
{ "Bi", 209, 209, 208.980390000, METAL , 19,  0, {{1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,},     {3,}       }},
{ "Po", 209, 209, 208.982400000, METAL2, 20,  0, {{0,},       {1,3,5,7,}, {2,4,6,},   {3,5,},     {2,4,}     }},
{ "At", 210, 210, 209.987100000,     0 , 22,  0, {{0,},       {0,},       {1,3,5,7,}, {2,4,6},    {3,5,}     }},
{ "Rn", 222, 222, 222.017500000,     0 ,  0,  0, {{0,},       {0,},       {0,},       {0,},       {0,}       }},
{ "Fr", 223, 223, 223.019700000, METAL ,  0,  0, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Ra", 226, 226, 226.025410000, METAL ,  0,  0, {{0,},       {0,},       {2,},       {1,},       {0,}       }},
{ "Ac", 227, 227, 227.027750000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Th", 232, 232, 232.038050000, METAL2,  0,  1, {{0,},       {0,},       {3,4,},     {0,},       {0,}       }},
{ "Pa", 231, 231, 231.035880000, METAL2,  0,  1, {{0,},       {0,},       {3,4,5,},   {0,},       {0,}       }},
{ "U",  238, 238, 238.050790000, METAL2,  0,  1, {{0,},       {0,},       {3,4,5,6,}, {0,},       {0,}       }},
{ "Np", 237, 237, 237.048170000, METAL2,  0,  1, {{0,},       {0,},       {3,4,5,6,}, {0,},       {0,}       }},
{ "Pu", 244, 244, 244.064200000, METAL2,  0,  1, {{0,},       {0,},       {3,4,5,6,}, {0,},       {0,}       }},
{ "Am", 243, 243, 243.061370000, METAL2,  0,  1, {{0,},       {0,},       {3,4,5,6,}, {0,},       {0,}       }},
{ "Cm", 247, 247, 247.070300000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Bk", 247, 247, 247.070300000, METAL ,  0,  1, {{0,},       {0,},       {3,4,},     {0,},       {0,}       }},
{ "Cf", 251, 251, 251.079600000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Es", 252, 252, 252.082800000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Fm", 257, 257, 257.095100000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "Md", 258, 258, 258.098600000, METAL ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }},
{ "No", 259, 259, 259.100900000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Lr", 260, 260, 260.105400000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Rf", 261, 261, 261.108700000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},

/*
    The elements below were added after v. 1.03.
    When available, the mass is given for isotope with the longest half-life.
    Standard valences given here are just placeholders.
    v. 1.04: added elements 105-112.
                Ref.: M. E. WIESER AND T. B. COPLEN.
                Atomic weights of the elements 2009 (IUPAC Technical Report).
                Pure Appl. Chem., Vol. 83, No. 2, pp. 359-396, 2011.
    v. 1.05: added elements 114 and 116;
                updated data for elements 105-112.
                Ref.: J. Meija, T.B. Coplen, M.Berglund et al.
                Atomic weights of the elements 2013 (IUPAC Technical Report).
                Pure Appl. Chem., Vol. 88, No. 3, pp. 265-291, 2016.
            added elements 113, 115, 117, and 118, according to IUPAC provisional recommendations:
                Ref.: L. Ohrstrom, J. Reedijk.
                Names and Symbols of the Elements with Atomic Numbers 113, 115, 117 and 118.
                Pure Appl. Chem., May 1, 2016, Manuscript ID PAC-REC-16-05-01
                http://iupac.org/cms/wp-content/uploads/2016/06/names-and-symbols-of-elements.pdf
*/

/* 105 dubnium Db                ? Like: Ta */
{ "Db", 270, 270, 270.131000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 106 seaborgium Sg            ? Like: W */
{ "Sg", 269, 269, 269.129000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 107 bohrium Bh                ? Like: Re */
{ "Bh", 270, 270, 270.133000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 108 hassium Hs                ? Like: Os */
{ "Hs", 270, 270, 270.134000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 109 meitnerium Mt            ? Like: Ir */
{ "Mt", 278, 278, 278.156000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 110 darmstadtium Ds            ? Like: Pt */
{ "Ds", 281, 281, 281.165000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 111 roentgenium Rg            ? Like: Au */
{ "Rg", 281, 281, 281.166000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 112 copernicium Cn            ? Like: Hg */
{ "Cn", 285, 285, 285.177000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 113 nihonium Nh                ? Like: ? */
{ "Nh", 278, 278, 278.000000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 114 flerovium Fl                ? Like: Pb */
{ "Fl", 289, 289, 289.190000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 115 moscovium Mc                ? Like: ? */
{ "Mc", 289, 289, 289.000000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 116 livermorium Lv            ? Like: Po */
{ "Lv", 293, 293, 293.204000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 117 tennessine Ts            ? Like: ? */
{ "Ts", 297, 297, 297.000000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 118 oganesson Og            ? Like: ? */
{ "Og", 294, 294, 294.000000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* End of added in v. 1.04 - 1.05 */
{ "Zy",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
{ "Zz",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
#ifdef INCHI_ZFRAG
{ "Zu",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }}, 
{ "Zv",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {2,},       {0,},       {0,}       }}, 
{ "Zw",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }}, 
{ "Zx",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {1,2,},     {0,},       {0,}       }}, 
#endif

{ "",     0,   0,   0.000000000,     0 ,  0,  0, {{0,},       {0,},       {0,},       {0,},       {0,}       }},
};


const int ERR_ELEM = 255;
const int nElDataLen = sizeof( ElData ) / sizeof( ElData[0] ) - 1;


/*
    MISC. CHEMICAL-STRUCTURE RELATED UTILITIES AND HELPERS
*/


/****************************************************************************
  Finds chemical symbol for element of given number.
  Returns 0 if OK and -1 if element was not found.
****************************************************************************/
int get_element_chemical_symbol( int nAtNum, char *szElement )
{
    nAtNum -= 1;

    if (0 < nAtNum)
    {
        nAtNum += 2; /*  bypass D, T */
    }

    if (0 <= nAtNum && nAtNum < nElDataLen)
    {
        /* valid element symbol found */
        strcpy(szElement, ElData[nAtNum].szElName);
        return 0;
    }

    /* not found */
    strcpy(szElement, "??");
    return -1;
}


/****************************************************************************
Finds symbol for element of given number.
Accounts for (translates)pseudoelements.
Returns 0 if OK and -1 if element was not found.
****************************************************************************/
int get_element_or_pseudoelement_symbol( int nAtNum,
                                         char *szElement )
{
    nAtNum -= 1;

    if (0 < nAtNum)
    {
        nAtNum += 2; /*  bypass D, T */
    }

    if (0 <= nAtNum && nAtNum < nElDataLen)
    {
        /* valid element symbol found */
        strcpy(szElement, ElData[nAtNum].szElName);

        if (!strcmp( szElement, "Zy" ))
        {
            strcpy(szElement, "Zz");
        }

        return 0;
    }

    /* not found */
    strcpy(szElement, "??");

    return -1;
}


/****************************************************************************/
int el_number_in_internal_ref_table( const char* elname )
{
    int i;
    const char *p;

    for (i = 0; ( p = ElData[i].szElName )[0] && strcmp( p, elname ); i++)
    {
        ;
    }

    return p[0] ? i : ERR_ELEM;
}


/****************************************************************************
  Get element number by symbol
****************************************************************************/
int get_periodic_table_number( const char* elname )
{
    int num;

    if (elname == NULL)
    {
        return ERR_ELEM;
    }

    if (strlen(elname) == 0)
    {
        return ERR_ELEM;
    }

    /* the single letter (common) elements */
    if (!elname[1]) 
    {
        switch (elname[0]) 
        {
            case 'H': 
                return EL_NUMBER_H; 
                break;
            case 'B': 
                return EL_NUMBER_B; 
                break;
            case 'C': 
                return EL_NUMBER_C; 
                break;
            case 'N': 
                return EL_NUMBER_N; 
                break;
            case 'O': 
                return EL_NUMBER_O; 
                break;
            case 'P': 
                return EL_NUMBER_P; 
                break;
            case 'S': 
                return EL_NUMBER_S; 
                break;
            case 'F': 
                return EL_NUMBER_F; 
                break;
            case 'I': 
                return EL_NUMBER_I; 
                break;
        }
    }

    num = el_number_in_internal_ref_table( elname );

    if (num < ERR_ELEM)
    {
        /* account for D,T in internal table (but not Mendeleev's table) */
        num = inchi_max( 1, num - 1 );
    }

    return num;
}


/****************************************************************************
  Check if no H addition allowed
****************************************************************************/
int if_skip_add_H( int nPeriodicNum )
/* was called if_skip_add_H(, renamed to avoid confusion with other procedures   */
{
    return
        ElData[nPeriodicNum > 1 ? nPeriodicNum + 1 : 0].bSkipAddingH;
}


/****************************************************************************
  Get reference value of atom valence at given charge
****************************************************************************/
int get_el_valence( int nPeriodicNum, int charge, int val_num )
{
    if (charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE || val_num >= MAX_NUM_VALENCES)
    {
        return 0;
    }

    return
        ElData[nPeriodicNum > 1 ? nPeriodicNum + 1 : 0].cValence[NEUTRAL_STATE + charge][val_num];
}


/****************************************************************************
    Output valence needed to unambiguosly reconstruct bonds
****************************************************************************/
int get_unusual_el_valence( int nPeriodicNum,
                            int charge,
                            int radical,
                            int bonds_valence,
                            int num_H,
                            int num_bonds )
{
    int i, num_found, chem_valence, rad_adj, known_chem_valence, exact_found;

    if (!num_bonds && !num_H)
    {
        return 0;
    }

    if (charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE)
    {
        if (bonds_valence == num_bonds)
        {
            return 0; /* all single bonds */
        }
        return bonds_valence;
    }

    if (!get_el_valence( nPeriodicNum, charge, 0 ) && bonds_valence == num_bonds)
    {
        return 0;
    }

    chem_valence = bonds_valence + num_H;
    rad_adj = 0;
    num_found = 0;
    exact_found = 0;

    /* Take into account a radical */
    if (radical == RADICAL_DOUBLET)
    {
        rad_adj = 1;
    }
    else if (radical == RADICAL_TRIPLET)
    {
        rad_adj = 2;
    }

    for (i = 0; i < MAX_NUM_VALENCES; i++)
    {
        if (0 < ( known_chem_valence = get_el_valence( nPeriodicNum, charge, i ) - rad_adj ) &&
             num_bonds <= known_chem_valence && known_chem_valence <= chem_valence)
        {
            num_found++;
            if (known_chem_valence == chem_valence)
            {
                exact_found = 1;
                break;
            }
        }
    }

    return ( exact_found && 1 == num_found ) ? 0 : chem_valence;
}


/****************************************************************************
    Output valence needed to unambiguosly reconstruct number of H
****************************************************************************/
int needed_unusual_el_valence( int nPeriodicNum,
                               int charge,
                               int radical,
                               int bonds_valence,
                               int actual_bonds_valence,
                               int num_H, int
                               num_bonds )
{
    int chem_valence, num_H_expected; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    char szElement[4];

    /*
    if ( !num_bonds && !num_H )
        return 0;
    */

    if (num_bonds && get_element_chemical_symbol( nPeriodicNum, szElement ) != -1)
    {
        num_H_expected = get_num_H( szElement, 0, NULL, charge, radical, actual_bonds_valence, 0, 0, 0, 0 );
    }
    else
    {
        num_H_expected = num_H;
    }

    chem_valence = bonds_valence + num_H;

#if ( (BUILD_WITH_ENG_OPTIONS==1) && (SDF_OUTPUT_HETERO_VALENCE==1) )
    if ((nPeriodicNum == 1 && chem_valence != 1) /* H */ || (nPeriodicNum == 6 && chem_valence != 4) /* C */ ||
        (nPeriodicNum != 1 && nPeriodicNum != 6) || charge || radical) /* djb-rwth: addressing LLVM warning */
    {
        return chem_valence ? chem_valence : -1;
    }
    else
    {
        return 0;
    }
#else
    {
        int i, num_found, num_found_known, rad_adj, known_chem_valence, exact_found;

        if (charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE ||
             !get_el_valence( nPeriodicNum, charge, 0 ) ||
             if_skip_add_H( nPeriodicNum ) || bonds_valence != actual_bonds_valence ||
             num_H_expected != num_H)
        {
            if (!num_H && !num_H_expected && bonds_valence == actual_bonds_valence)
            {
                return 0; /* no H */
            }
            return chem_valence; /* needs to add H-atoms */
        }

        /* take into account radical */
        if (radical == RADICAL_DOUBLET)
        {
            rad_adj = 1;
        }
        else if (radical == RADICAL_TRIPLET)
        {
            rad_adj = 2;
        }
        else
        {
            rad_adj = 0;
        }

        num_found_known = 0;
        num_found = 0;
        exact_found = 0;

        for (i = 0; i < MAX_NUM_VALENCES; i++)
        {
            if (0 < ( known_chem_valence = get_el_valence( nPeriodicNum, charge, i ) ) &&
                 bonds_valence <= ( known_chem_valence -= rad_adj ))
            {
                /* found known valence that fits without H */
                num_found_known++;
                if (known_chem_valence <= chem_valence)
                {
                    /* known valence is large enough to accommodate (implicit) H */
                    num_found++;
                }
                if (known_chem_valence == chem_valence)
                {
                    exact_found = 1;
                    break;
                }
            }
        }

        return ( exact_found && 1 == num_found && 1 == num_found_known )
            ? 0
            : chem_valence ? chem_valence : -1;    /* needs zero */
    }
#endif
}


/****************************************************************************
    Output valence that does not fit any known valences
****************************************************************************/
int detect_unusual_el_valence( int nPeriodicNum,
                               int charge,
                               int radical,
                               int bonds_valence,
                               int num_H,
                               int num_bonds )
{
    int i, chem_valence, rad_adj, known_chem_valence;

    if (!num_bonds && !num_H)
    {
        return 0;
    }

    if (charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE)
    {
        if (bonds_valence == num_bonds)
        {
            return 0; /* all single bonds */
        }
        return bonds_valence;
    }

    if (!get_el_valence( nPeriodicNum, charge, 0 ) && bonds_valence == num_bonds)
    {
        return 0;
    }

    chem_valence = bonds_valence + num_H;
    rad_adj = 0;

    /* take into account radical */
    if (radical == RADICAL_DOUBLET)
    {
        rad_adj = 1;
    }
    else if (radical == RADICAL_TRIPLET || radical == RADICAL_SINGLET)
    {
        rad_adj = 2;
    }

    for (i = 0; i < MAX_NUM_VALENCES; i++)
    {
        if (0 < ( known_chem_valence = get_el_valence( nPeriodicNum, charge, i ) - rad_adj ))
        {
            if (known_chem_valence == chem_valence)
            {
                return 0;
            }
        }
    }

    return chem_valence;
}


/****************************************************************************
 Return element type
****************************************************************************/
int get_el_type( int nPeriodicNum )
{
    return ElData[nPeriodicNum + 1].nType;
}


/****************************************************************************
 Check if element is metal
****************************************************************************/
int is_el_a_metal( int nPeriodicNum )
{
    return 0 != ( ElData[nPeriodicNum + 1].nType & IS_METAL );
}


/*#ifndef TARGET_API_LIB*/


/****************************************************************************
 Extract radicals and charges
****************************************************************************/
int extract_charges_and_radicals( char *elname, int *pnRadical, int *pnCharge )
{
    char *q, *r, *p;
    int  nCharge = 0, nRad = 0, charge_len = 0, k, nVal, nSign, nLastSign = 1; /* djb-rwth: removing redundant variables */

    p = elname;

    /*  extract radicals & charges */
    while ((q = strpbrk(p, "+-^"))) /* djb-rwth: addressing LLVM warning */
    {
        switch (*q)
        {
            case '+':
            case '-':
                for (k = 0, nVal = 0; ( nSign = ( '+' == q[k] ) ) || ( nSign = -( '-' == q[k] ) ); k++)
                {
                    nVal += ( nLastSign = nSign );
                    charge_len++;
                }
                if ((nSign = (int)strtol(q + k, &r, 10))) /* djb-rwth: addressing LLVM warning */
                {
                    /*  fixed 12-5-2001 */
                    nVal += nLastSign * ( nSign - 1 );
                }
                charge_len = (int) ( r - q );
                nCharge += nVal;
                break;
            /* case '.': */ /*  singlet '.' may be confused with '.' in formulas like CaO.H2O */
            case '^':
                nRad = 1; /* doublet here is 1. See below */
                charge_len = 1;
                for (k = 1; q[0] == q[k]; k++)
                {
                    nRad++;
                    charge_len++;
                }
                break;
        }
        memmove(q, q + charge_len, strlen(q + charge_len) + 1);
    }

    /* djb-rwth: removing redundant code */

    /*  radical */
    if (( q = strrchr( p, ':' ) ) && !q[1])
    {
        nRad = RADICAL_SINGLET;
        q[0] = '\0';
        /* djb-rwth: removing redundant code */
    }
    else
    {
        while (( q = strrchr( p, '.' ) ) && !q[1])
        {
            nRad++;
            q[0] = '\0';
            /* djb-rwth: removing redundant code */
        }

        nRad = nRad == 1 ? RADICAL_DOUBLET :
            nRad == 2 ? RADICAL_TRIPLET : 0;
    }

    *pnRadical = nRad;
    *pnCharge = nCharge;

    return ( nRad || nCharge );
}


/*#endif*/


/****************************************************************************/
int extract_H_atoms( char *elname, S_CHAR num_iso_H[] )
{
    int i, len, c, k, num_H, val;
    char *q;
    char elname1 = '\0';

    i = 0;
    num_H = 0;
    len = (int) strlen( elname );
    c = UCINT elname[0];

    if (len > 1)
    {
        elname1 = elname[1];
    }

    while (i < len)
    {
        switch (c)
        {
            case 'H':
                k = 0;
                break;
            case 'D':
                k = 1;
                break;
            case 'T':
                k = 2;
                break;
            default:
                k = -1;
                break;
        }

        q = elname + i + 1; /*  pointer to the next to elname[i] character */
        c = UCINT q[0];

        if (k >= 0 && !islower( c ))
        {
            /*  found a hydrogen */
            if (isdigit( c ))
            {
                val = (int) strtol( q, &q, 10 );
                /*  q = pointer to the next to number of hydrogen atom(s) character */
            }
            else
            {
                val = 1;
            }
            if (k)
            {
                num_iso_H[k] += val;
            }
            else
            {
                num_H += val;
            }

            /*  remove the hydrogen atom from the string */
            len -= (int) ( q - elname ) - i;
            memmove(elname + i, q, (long long)len + 1); /* djb-rwth: cast operator added */
            /*  c =  UCINT elname[i]; */
        }
        else
        {
            i++;
        }

        c = UCINT elname[i]; /*  moved here 11-04-2002 */
    }

    len = (int) strlen( elname );
    if (len == 2)
    {
        if (elname[1] != elname1)
            /* Error, incorrect 2nd char of elname appears after 'subtracting' {H,D,T}  */
            /* See a bug reported to inchi-discuss by A. Dalke for alias atom "pH4d"    */
            /*^^^ 2017-01-06                                                            */
            elname[1] = '?';
    }

    return num_H;
}


/****************************************************************************
 Return number of attached hydrogens
****************************************************************************/
int get_num_H( const char* elname,
                int inp_num_H,
                S_CHAR inp_num_iso_H[],
                int charge,
                int radical,
                int chem_bonds_valence,
                int atom_input_valence,
                int bAliased,
                int bDoNotAddH,
                int bHasMetalNeighbor )
{
    int val, i, el_number, num_H = 0, num_iso_H;
    static int intl_el_number_N = 0, intl_el_number_S=0, intl_el_number_O=0, intl_el_number_C=0;

    if (!intl_el_number_N)
    {
        intl_el_number_N = el_number_in_internal_ref_table( "N" );
    }
    if (!intl_el_number_S)
    {
        intl_el_number_S = el_number_in_internal_ref_table( "S" );
    }
    if (!intl_el_number_O)
    {
        intl_el_number_O = el_number_in_internal_ref_table( "O" );
    }
    if (!intl_el_number_C)
    {
        intl_el_number_C = el_number_in_internal_ref_table( "C" );
    }


    /*  atom_input_valence (cValence) cannot be specified in case of */
    /*  aliased MOLFile atom with known inp_num_H or inp_num_iso_H[] */

    if (bAliased)
    {
        num_H = inp_num_H;
    }
    else if (atom_input_valence && ( atom_input_valence != 15 || chem_bonds_valence ))
    {
        num_H = inchi_max( 0, atom_input_valence - chem_bonds_valence );
    }
    else if (atom_input_valence == 15 && !chem_bonds_valence)
    {
        num_H = 0;
    }
    else if (MIN_ATOM_CHARGE <= charge &&
              MAX_ATOM_CHARGE >= charge &&
              ERR_ELEM != ( el_number = el_number_in_internal_ref_table( elname ) ) &&
              !ElData[el_number].bSkipAddingH && !bDoNotAddH)
    {
        /* add hydrogen atoms according to standard element valence */
        if (radical && radical != RADICAL_SINGLET)
        {
            if ((val = ElData[el_number].cValence[NEUTRAL_STATE + charge][0])) /* djb-rwth: addressing LLVM warning */
            {
                val -= ( radical == RADICAL_DOUBLET ) ? 1
                    : ( radical == RADICAL_SINGLET || radical == RADICAL_TRIPLET ) ? 2 : val;
                /* if unknown radical then do not add H */
                num_H = inchi_max( 0, val - chem_bonds_valence );
            }
        }
        else
        {
            /* find the smallest valence that is greater than the sum of the chemical bond valences */
            for (i = 0;
                 ( val = ElData[el_number].cValence[NEUTRAL_STATE + charge][i] ) &&
                 val < chem_bonds_valence;
                 i++)
            {
                ;
            }

            /* special case: do not add H to N(IV), S(III), S+(II), S-(II) */ /* S ions added 2004-05-10 */
            if (el_number == intl_el_number_N && !charge && !radical && val == 5)
            {
                val = 3;
            }
            /*else if ( el_number == el_number_N && !charge && !radical && val == 3 &&
                 chem_bonds_valence == 2 && bHasMetalNeighbor )
              {
              val = 2;
              }
            */
            else if (el_number == intl_el_number_S && !charge && !radical && val == 4 && chem_bonds_valence == 3)
            {
                val = 3;
            }
            else if (bHasMetalNeighbor && el_number != intl_el_number_C && val > 0)
            {
                val--;
            }
            /*
            if ( (el_number == el_number_S || el_number == el_number_O) &&
                 abs(charge)==1 && !radical && val == 3 && chem_bonds_valence == 2 && bHasMetalNeighbor )
              {
              val = 2;
              }
            else
            */

            num_H = inchi_max( 0, val - chem_bonds_valence );
        }

        num_iso_H = 0;
        if (inp_num_iso_H)
        {
            for (i = 0; i < NUM_H_ISOTOPES; i++)
            {
                num_iso_H += inp_num_iso_H[i];
            }
        }

        /*  should not happen because atom here is not aliased */
        if (num_iso_H)
        {
            if (num_H >= num_iso_H)
            {
                num_H -= num_iso_H;
            }
            else
            {
                num_H = inp_num_H; /*  as requested in the alias */
                /* num_H = (num_iso_H - num_H) % 2; */ /*  keep unchanged parity of the total number of H atoms */
            }
        }

        /*  should not happen because atom here is not aliased */
        if (inp_num_H > num_H)
        {
            num_H = inp_num_H;  /*  as requested in the alias */
            /* num_H = inp_num_H + (inp_num_H - num_H)%2; */ /*  keep unchanged parity of the number of non-isotopic H atoms */
        }
    }
    else
    {
        num_H = inp_num_H;
    }

    return num_H;
}


/****************************************************************************/
int get_atomic_mass_from_elnum( int nAtNum )
{
    nAtNum -= 1;

    if (0 < nAtNum)
    {
        nAtNum += 2; /*  bypass D, T */
    }

    if (0 <= nAtNum && nAtNum < nElDataLen)
    {
        return (int) ElData[nAtNum].nAtMass;
    }

    return 0;
}


/****************************************************************************
int get_mw(char elname[])
{
    int i;

    for (i=0; i<NUMEL; i++)
        if (strcmp(elname,elements[i])==0)
            return(atomic_wt[i]);
    return(0);
}
****************************************************************************/



/****************************************************************************/
int get_atomic_mass( const char *elname )
{
    int el_number, atw;
    if (ERR_ELEM != ( el_number = el_number_in_internal_ref_table( elname ) ))
    {
        atw = ElData[el_number].nAtMass;
    }
    else
    {
        atw = 0;
    }

    return atw;
}


/****************************************************************************
  Check if atom is in the list
****************************************************************************/
AT_NUMB *is_in_the_list( AT_NUMB *pathAtom, AT_NUMB nNextAtom, int nPathLen )
{
    for (; nPathLen && *pathAtom != nNextAtom; nPathLen--, pathAtom++)
    {
        ;
    }
    return nPathLen ? pathAtom : NULL;
}


/****************************************************************************
  Check if integer is in the list
****************************************************************************/
int *is_in_the_ilist( int *pathAtom, int nNextAtom, int nPathLen )
{
    for (; nPathLen && *pathAtom != nNextAtom; nPathLen--, pathAtom++)
    {
        ;
    }
    return nPathLen ? pathAtom : NULL;
}


/****************************************************************************
  Check if list (of integers) is within the embedding list2
****************************************************************************/
int is_ilist_inside( int *ilist, int nlist, int *ilist2, int nlist2 )
{
    int k;
    for (k = 0; k < nlist; k++)
    {
        if (!is_in_the_ilist( ilist2, ilist[k], nlist2 ))
        {
            return 0;
        }
    }
    return 1;
}


/****************************************************************************/
int nBondsValToMetal( inp_ATOM* at, int iat )
{
    int i, bond_type, nVal2Metal = 0; /* djb-rwth: removing redundant variables */
    inp_ATOM* a = at + iat;

    for (i = 0; i < a->valence; i++)
    {
        /* djb-rwth: removing redundant code */

        if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
        {
            bond_type = a->bond_type[i];

            if (bond_type <= BOND_TYPE_TRIPLE)
            {
                nVal2Metal += bond_type;
            }
            else
            {
                return -1;  /* bond to metal order is not well defined */
            }
        }
    }

    return nVal2Metal;
}


/****************************************************************************/
int num_of_H( inp_ATOM *at, int iat )
{
    static int el_number_H = (int)EL_NUMBER_H;
    int    i, n, num_explicit_H = 0;
    inp_ATOM *a = at + iat;

    for (i = 0; i < a->valence; i++)
    {
        n = a->neighbor[i];
        num_explicit_H += ( 1 == at[n].valence && el_number_H == at[n].el_number );
    }

    return num_explicit_H + NUMH( at, iat );
}

/****************************************************************************/
/* Get the element group of an element. The base element rather than the    */
/* periodic group is used to aid readability.                               */
/* - NitrogenGroup = 7 (EL_NUMBER_N)                                        */
/* - OxygenGroup   = 8 (EL_NUMBER_O)                                        */
/* - Cargbon       = 6 (EL_NUMBER_C)                                        */
/****************************************************************************/
U_CHAR ion_el_group( int el )
{
    switch ( el ) 
    {
        case EL_NUMBER_C: /* fallthrough */
#if ( FIX_REM_ION_PAIRS_Si_BUG == 1 )        
        case EL_NUMBER_SI:
#endif        
            return EL_NUMBER_C;
        case EL_NUMBER_N: /* fallthrough */
        case EL_NUMBER_P:
        case EL_NUMBER_AS:
        case EL_NUMBER_SB:
            return EL_NUMBER_N;
        case EL_NUMBER_O: /* fallthrough */
        case EL_NUMBER_S:
        case EL_NUMBER_SE:
        case EL_NUMBER_TE:
            return EL_NUMBER_O;
        default:
            return 0;
    }
}

int has_other_ion_neigh( inp_ATOM *at,
                         int iat,
                         int iat_ion_neigh)
{
    int charge = at[iat_ion_neigh].charge;
    int i, neigh;

    for (i = 0; i < at[iat].valence; i++)
    {
        neigh = at[iat].neighbor[i];

        if (neigh != iat_ion_neigh && at[neigh].charge == charge &&
            ion_el_group( at[neigh].el_number ))
        {
            return 1;
        }
    }

    return 0;
}


/****************************************************************************
  Check if has_other_ion_in_sphere_2
    BFS r=2
****************************************************************************/
int has_other_ion_in_sphere_2( inp_ATOM *at, int iat,
                               int iat_ion_neigh )
{
#define MAXQ 16
    AT_NUMB q[MAXQ];
    int lenq = 0, lenq2, dist = 0, i = 0, iq, neigh, j, nRet = 0;
    q[lenq++] = iat;
    at[iat].cFlags = 1;

    iq = 0;
    dist = 1;
    /* use at->cFlags as an indicator */

    while (dist <= 2)
    {
        for (lenq2 = lenq; iq < lenq2; iq++)
        {
            i = q[iq];

            for (j = 0; j < at[i].valence; j++)
            {
                neigh = at[i].neighbor[j];

                if (!at[neigh].cFlags &&
                     at[neigh].valence <= 3 &&
                     ion_el_group( at[neigh].el_number ))
                {
                    q[lenq++] = neigh;
                    at[neigh].cFlags = 1;
                    if (neigh != iat_ion_neigh &&
                         at[iat_ion_neigh].charge == at[neigh].charge)
                    {
                        nRet++;
                    }
                }
            }
        }

        dist++;
    }

    for (iq = 0; iq < lenq; iq++)
    {
        i = q[iq];
        at[i].cFlags = 0;
    }

    return nRet;
}


/****************************************************************************/
int nNoMetalNumBonds( inp_ATOM *at, int at_no )
{
    int i;

    inp_ATOM *a = at + at_no;
    int num_H = NUMH( a, 0 );
    int std_chem_bonds_valence = get_el_valence( a->el_number, a->charge, 0 );

    if (a->chem_bonds_valence + num_H > std_chem_bonds_valence)
    {
        int valence_to_metal = 0;
        int num_bonds_to_metal = 0;

        for (i = 0; i < a->valence; i++)
        {
            if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
            {
                if (( a->bond_type[i] & BOND_TYPE_MASK ) >= BOND_TYPE_ALTERN)
                {
                    return a->valence; /* fall back */
                }
                num_bonds_to_metal++;
                valence_to_metal += ( a->bond_type[i] & BOND_TYPE_MASK );
            }
        }

        if (a->chem_bonds_valence + num_H - valence_to_metal == std_chem_bonds_valence)
        {
            /* removing bonds to metal produces standard valence */
            return a->valence - num_bonds_to_metal;
        }
    }

#if ( S_VI_O_PLUS_METAL_FIX_BOND == 1 )
    else
    {
        if (1 == a->charge && 2 == get_endpoint_valence( a->el_number ) &&
                a->chem_bonds_valence + num_H == std_chem_bonds_valence)
        {
            int valence_to_metal = 0;
            int num_bonds_to_metal = 0;
            for (i = 0; i < a->valence; i++)
            {
                if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
                {
                    if (( a->bond_type[i] & BOND_TYPE_MASK ) >= BOND_TYPE_ALTERN)
                    {
                        return a->valence; /* fall back */
                    }
                    num_bonds_to_metal++;
                    valence_to_metal += ( a->bond_type[i] & BOND_TYPE_MASK );
                }
            }
            if (1 == valence_to_metal)
            {
                /* removing bonds to metal produces standard valence */
                return a->valence - num_bonds_to_metal;
            }
        }
    }
#endif

    return a->valence;
}


/****************************************************************************/
int nNoMetalBondsValence( inp_ATOM *at, int at_no )
{
    int i;

    inp_ATOM *a = at + at_no;
    int num_H = NUMH( a, 0 );
    int std_chem_bonds_valence = get_el_valence( a->el_number, a->charge, 0 );

    if (a->chem_bonds_valence + num_H > std_chem_bonds_valence)
    {
        int valence_to_metal = 0;
        /*int num_bonds_to_metal = 0;*/

        for (i = 0; i < a->valence; i++)
        {
            if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
            {
                if (( a->bond_type[i] & BOND_TYPE_MASK ) >= BOND_TYPE_ALTERN)
                {
                    return a->valence; /* fall back */
                }
                /* num_bonds_to_metal ++;*/
                valence_to_metal += ( a->bond_type[i] & BOND_TYPE_MASK );
            }
        }

        if (a->chem_bonds_valence + num_H - valence_to_metal == std_chem_bonds_valence)
        {
            /* removing bonds to metal produces standard valence */
            return a->chem_bonds_valence - valence_to_metal;
        }
    }

#if ( S_VI_O_PLUS_METAL_FIX_BOND == 1 )
    else if (1 == a->charge && 2 == get_endpoint_valence( a->el_number ) &&
         a->chem_bonds_valence + num_H == std_chem_bonds_valence)
    {

        int valence_to_metal = 0;
        /* int num_bonds_to_metal = 0;*/

        for (i = 0; i < a->valence; i++)
        {
            if (is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
            {
                if (( a->bond_type[i] & BOND_TYPE_MASK ) >= BOND_TYPE_ALTERN)
                {
                    return a->valence; /* fall back */
                }
                /* num_bonds_to_metal ++;*/
                valence_to_metal += ( a->bond_type[i] & BOND_TYPE_MASK );
            }
        }

        if (1 == valence_to_metal)
        {/* removing bonds to metal produces standard valence */
            return a->chem_bonds_valence - valence_to_metal;
        }
    }
#endif

    return a->chem_bonds_valence;
}


/****************************************************************************/
int nNoMetalNeighIndex( inp_ATOM *at, int at_no )
{
    int i;

    inp_ATOM *a = at + at_no;

    for (i = 0; i < a->valence; i++)
    {
        if (!is_el_a_metal( at[(int) a->neighbor[i]].el_number ))
        {
            return i;
        }
    }

    return -1;
}


/****************************************************************************/
int nNoMetalOtherNeighIndex( inp_ATOM *at, int at_no, int cur_neigh )
{
    int i, neigh;

    inp_ATOM *a = at + at_no;

    for (i = 0; i < a->valence; i++)
    {
        neigh = (int) a->neighbor[i];

        if (neigh != cur_neigh && !is_el_a_metal( at[neigh].el_number ))
        {
            return i;
        }
    }

    return -1;
}


/****************************************************************************/
int nNoMetalOtherNeighIndex2( inp_ATOM *at,
                              int at_no,
                              int cur_neigh,
                              int cur_neigh2 )
{
    int i, neigh;

    inp_ATOM *a = at + at_no;

    for (i = 0; i < a->valence; i++)
    {
        neigh = (int) a->neighbor[i];

        if (neigh != cur_neigh && neigh != cur_neigh2 && !is_el_a_metal( at[neigh].el_number ))
        {
            return i;
        }
    }

    return -1;
}



#ifndef COMPILE_ANSI_ONLY


/****************************************************************************/
int MakeRemovedProtonsString( int nNumRemovedProtons,
                              NUM_H *nNumExchgIsotopicH,
                              NUM_H *nNumRemovedProtonsIsotopic,
                              int bIsotopic,
                              char *szRemovedProtons,
                              int *num_removed_iso_H )
{
    int i, j, len, num;

    len = 0;

    if (nNumRemovedProtons)
    {
        len = sprintf(szRemovedProtons, "Proton balance: %c %d H+",
            nNumRemovedProtons >= 0 ? '+' : '-', abs(nNumRemovedProtons));
    }

    if (bIsotopic && ( nNumRemovedProtonsIsotopic || nNumExchgIsotopicH ))
    {

        for (i = 0, j = 0; i < NUM_H_ISOTOPES; i++)
        {

            num = ( nNumExchgIsotopicH ? nNumExchgIsotopicH[i] : 0 ) +
                ( nNumRemovedProtonsIsotopic ? nNumRemovedProtonsIsotopic[i] : 0 );

            if (num)
            {
                len += sprintf(szRemovedProtons + len, "%s %d^%dH",
                    j ? ", " : "  [ removed ", num, i + 1);
                j++;
            }
        }

        if (j)
        {
            len += sprintf(szRemovedProtons + len, " ]");
            if (num_removed_iso_H)
            {
                *num_removed_iso_H = j;
            }
        }
    }

    if (!len)
    {
        szRemovedProtons[0] = '\0';
    }

    return len;
}
#endif


/****************************************************************************/
int get_endpoint_valence( U_CHAR el_number )
{   
    switch (el_number) 
    {
        case EL_NUMBER_O:  /* fallthrough */
        case EL_NUMBER_S:  
        case EL_NUMBER_SE: 
        case EL_NUMBER_TE: 
            return 2;
        case EL_NUMBER_N:  
            return 3;
        default: 
            return 0;
    }
}


#if ( KETO_ENOL_TAUT == 1 )  /* post v.1 feature */


/****************************************************************************/
int get_endpoint_valence_KET( U_CHAR el_number )
{
    switch (el_number) 
    {
        case EL_NUMBER_C: 
            return 4;
        case EL_NUMBER_O: 
            return 2;
        default: 
            return 0;
    }
}
#endif


/*
    MEMORY MANAGEMENT
*/


#ifndef inchi_malloc
/****************************************************************************/
void *inchi_malloc( size_t c )
{
    return  malloc( c );
}
#endif


#ifndef inchi_calloc
/****************************************************************************/
void *inchi_calloc( size_t c, size_t n )
{
    return calloc( c, n );
}
#endif


#ifndef inchi_free
/****************************************************************************/
void inchi_free( void *p )
{
    if (p)
    {
        inchi_free( p ); /*added check if zero*/
    }
}
#endif


/*
    STRINGS AND TEXT HANDLING
*/


/****************************************************************************
  Remove leading & trailing spaces, replace consecutive
  spaces with a single space.
****************************************************************************/
int normalize_string( char* name )
{
    int i, len, n;

    len = (int) strlen( name );

    for (i = 0, n = 0; i < len; i++)
    {
        if (isspace( UCINT name[i] ) /*|| !isprint( UCINT name[i] )*/)
        {
            name[i] = ' '; /* exterminate tabs !!! */
            n++;
        }
        else
        {
            if (n > 0)
            {
                memmove((void*)&name[i - n], (void*)&name[i], (long long)len - (long long)i + 1); /* djb-rwth: cast operators added */
                i -= n;
                len -= n;
            }
            n = -1;
        }
    }
    if (n == len) /* empty line */
    {
        name[len = 0] = '\0';
    }
    else if (++n && n <= len)
    {
        len -= n;
        name[len] = '\0';
    }

    return len;
}


/****************************************************************************
 Replace non-ASCII characters with '.', return number of replacements
****************************************************************************/
int dotify_non_printable_chars( char *line )
{
    int i, c, num = 0;

    if (line)
    {
        for (i = 0; (c = UCINT line[i]); i++) /* djb-rwth: addressing LLVM warning */
        {
            /* assuming ASCII charset */
            if (c < ' ' || c >= 0x7F)
            {
                line[i] = '.';
                num++;
            }
        }
    }

    return num;
}


/****************************************************************************
    Reads char sequence pointed to by *pstring ( char *p = *ppstring) for
    not more than maxlen bytes) to 'field' up to first occurrence of any of
    delimiters in 'delims' or end of line, whichever occurs first.
    Sets *pstring to point to character which matches delimiter.
    Returns number of bytes copied, -1 on error.
****************************************************************************/
int read_upto_delim( char **pstring, char *field, int maxlen, char* delims )
{
    int i, n;
    char *p = *pstring;

    if (!p)
    {
        return -1;
    }

    /* skip leading spaces */
    for (i = 0; p[i] && isspace( UCINT p[i] ); i++)
    {
        ;
    }
    p += i;

    /* read up to next delim or eol */
    n = 0;
    while (p[n] && !is_matching_any_delim( p[n], delims ))
    {
        n++;
    }

    if (n + 1 > maxlen)
    {
        return -1;
    }

    mystrncpy( field, p, n + 1 );
    field[n + 1] = '\0';

    if (!p[n])
    {
        /* reached EOL */
        *pstring = NULL;
    }
    else
    {
        /* advance reading pos */
        *pstring = *pstring + i + n;
    }

    return n;
}


/****************************************************************************
    Check if a character is in the list of possible delimiters
    NB: same as isspace if delims is " \t\n\v\f\r"
    (0x20 and 0x09-0x0D)
****************************************************************************/
int is_matching_any_delim( char c, char* delims )
{
    int ic = UCINT c;
    while (*delims)
    {
        if (ic == *delims)
        {
            return 1;
        }
        delims++;
    }
    return 0;
}


/****************************************************************************
  Remove trailing spaces
****************************************************************************/
void remove_trailing_spaces( char* p )
{
    int   len;
    for (len = (int) strlen( p ) - 1; len >= 0 && isspace( UCINT p[len] ); len--)
    {
        ;
    }
    p[++len] = '\0';
}


/****************************************************************************/
void remove_one_lf( char* p )
{
    size_t len;
    if (p && 0 < ( len = strlen( p ) ) && p[len - 1] == '\n')
    {
        p[len - 1] = '\0';
        if (len >= 2 && p[len - 2] == '\r')
        {
            p[len - 2] = '\0';
        }
    }
}


/****************************************************************************
    Copies up to maxlen characters INCLUDING end null from source to target
    Fills out the rest of the target with null bytes

    protected from non-zero-terminated source and overlapped target/source.
****************************************************************************/
int mystrncpy( char *target, const char *source, unsigned maxlen )
{
    const char *p;
    unsigned len, source_len;

    if (target == NULL || maxlen == 0 || source == NULL)
    {
        return 0;
    }

    /* giallu: PR #163 */
    /* Find actual source length first to limit memchr search */ 
    source_len = (unsigned)strlen(source);

    if (source_len < maxlen)
    {
        /* Source is shorter than maxlen, use actual source length */
        len = source_len;
    }
    else if ((p = (const char*)memchr(source, 0, maxlen))) /* djb-rwth: addressing LLVM warning */
    {    
        /* maxlen does not include the found zero termination */
        len = (int) ( p - source );
    }
    else
    {    
        /*  reduced length does not include one more byte for zero termination */
        len = maxlen - 1;
    }

    if (len)
    {
        memmove(target, source, len);
    }

    memset(target + len, 0, maxlen - len); /*  zero termination */ /* djb-rwth: memset_s C11/Annex K variant? */

    return 1;
}


/****************************************************************************
 Remove leading and trailing white spaces
****************************************************************************/
char* lrtrim( char *p, int* nLen )
{
    int i, len = 0;

    if (p && ( len = (int) strlen( p ) ))
    {
        for (i = 0; i < len && is_ascii( p[i] ) && isspace( p[i] ); i++)
        {
            ;
        }
        if (i)
        {
            len -= i; /* djb-rwth: variable has to be decreased before memmove */
            (memmove)(p, p + i, ((long long)len + 1)); /* djb-rwth: now cast operator can be added */
        }
            
        for (; 0 < len && is_ascii( p[len - 1] ) && isspace( p[len - 1] ); len--)
        {
            ;
        }
        p[len] = '\0';
    }

    if (nLen)
    {
        *nLen = len;
    }

    return p;
}


 /****************************************************************************
   Extract InChI substring embedded into a longer string.

    InChI should start from "InChI=".

    As for the end of InChI,
    consider that according to
        http://info-uri.info/registry/OAIHandler?verb=GetRecord&metadataPrefix=reg&identifier=info:inchi/
    an InChI identifier may contain the following characters:
        A-Z
        a-z
        0-9
        ()*+,-./;=?@

    Here we treat any character not conforming this specification as a whitespace
    which marks the end of the InChI string.
    For example:
        "InChI=1/Ar%"
        "InChI=1/Ar\n"
        "InChI=1/Ar\r\t"
    all will be trimmed to
        "InChI=1/Ar"

****************************************************************************/
void extract_inchi_substring( char ** buf, const char *str, size_t slen )
{
    size_t i;
    const char *p;
    char* bufp;
    char pp;

    bufp = *buf;
    *buf = NULL;

    if (str == NULL)
    {
        return;
    }
    if (strlen( str ) < 1)
    {
        return;
    }

    p = strstr( str, "InChI=" );
    if (NULL == p)
        return;

    for (i = 0; i < slen; i++)
    {
        pp = p[i];

        if (pp >= 'A' && pp <= 'Z')   continue;
        if (pp >= 'a' && pp <= 'z')   continue;
        if (pp >= '0' && pp <= '9')   continue;
        switch (pp)
        {
            case '(':
            case ')':
            case '*':
            case '+':
            case ',':
            case '-':
            case '.':
            case '/':
            case ';':
            case '=':
            case '?':
            case '@':    continue;

            default:    break;
        }

        break;
    }

    *buf = (char*) inchi_calloc( i + 1, sizeof( char ) );
    memcpy(*buf, p, i);
    if (*buf)
        (*buf)[i] = '\0';

    return;
}


/****************************************************************************/
void extract_auxinfo_substring( char ** buf, const char *str, size_t slen )
{
    size_t i;
    const char *p;
    char* bufp;
    char pp;

    bufp = *buf;
    *buf = NULL;

    if (str == NULL)
    {
        return;
    }
    if (strlen( str ) < 1)
    {
        return;
    }

    p = strstr( str, "AuxInfo=" );
    if (NULL == p)
    {
        return;
    }

    for (i = 0; i < slen; i++)
    {
        pp = p[i];
        if (isspace( UCINT pp ))    break;
    }

    *buf = (char*) inchi_calloc( i + 1, sizeof( char ) );
    memcpy(*buf, p, i);
    if (*buf)
        (*buf)[i] = '\0';

    return;
}


/****************************************************************************
    For compatibility: local implementation of non-ANSI
    (MS-specific) functions, prefixed with "inchi_"
****************************************************************************/


#define __MYTOLOWER(c) ( ((c) >= 'A') && ((c) <= 'Z') ? ((c) - 'A' + 'a') : (c) )


/****************************************************************************/
int inchi_memicmp( const void * p1, const void * p2, size_t length )
{
    const U_CHAR *s1 = (const U_CHAR*) p1;
    const U_CHAR *s2 = (const U_CHAR*) p2;
    while (length--)
    {
        if (*s1 == *s2 ||
              __MYTOLOWER( (int) *s1 ) == __MYTOLOWER( (int) *s2 ))
        {
            s1++;
            s2++;
        }
        else
        {
            return
                __MYTOLOWER( (int) *s1 ) - __MYTOLOWER( (int) *s2 );
        }
    }

    return 0;
}


/****************************************************************************/
int inchi_stricmp( const char *s1, const char *s2 )
{
    while (*s1)
    {
        if (*s1 == *s2 ||
              __MYTOLOWER( (int) *s1 ) == __MYTOLOWER( (int) *s2 ))
        {
            s1++;
            s2++;
        }
        else
        {
            return
                __MYTOLOWER( (int) *s1 ) - __MYTOLOWER( (int) *s2 );
        }
    }

    if (*s2)
    {
        return -1;
    }

    return 0;
}


/****************************************************************************/
char *inchi__strnset( char *s, int val, size_t length )
{
    char *ps = s;
    while (length-- && *ps)
    {
        *ps++ = (char) val;
    }

    return s;
}


/****************************************************************************/
char *inchi__strdup( const char *string )
{
    char *p = NULL;
    if (string)
    {
        size_t length = strlen( string );
        p = (char *) inchi_malloc( length + 1 );
        if (p)
        {
            strcpy(p, string);
        }
    }

    return p;
}


#undef __MYTOLOWER
