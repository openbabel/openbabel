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


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "mode.h"

#include "incomdef.h"
#include "inpdef.h"
#include "util.h"
#include "extr_ct.h"


#include "ichicomp.h"

#define MIN_ATOM_CHARGE        (-2)
#define MAX_ATOM_CHARGE         2
#define NEUTRAL_STATE          (-MIN_ATOM_CHARGE)
#define NUM_ATOM_CHARGES       (MAX_ATOM_CHARGE - MIN_ATOM_CHARGE + 1)
#define MAX_NUM_VALENCES        5                /* max. number + 1 to provide zero termination */

typedef struct tagElData {
     const char *szElName;
     int     nAtMass;      /* Avg. atomic mass from the Periodic Chart of the Elements (Fisher cat. no. 05-702-10) */
     int     nNormAtMass;  /* Atomic mass of the most abundant isotope */
     double  dAtMass;      /* exact mw of the most abundant isotope */
     int     nType;        /* METAL or METAL2 */
     int     nElNegPauling10; /* Pauling electronegativity x 10; 0 => unknown */
     int     bDoNotAddH;   /* InChI does not add implicit H to atoms that have bDoNotAddH != 0 */
     S_CHAR  cValence[NUM_ATOM_CHARGES][MAX_NUM_VALENCES];
} ELDATA;

/* 2004=05-10: Added valences {1,3,5,7,} for As(2-) */

const ELDATA ElData[] = {
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

/*^^^ Added in v. 1.04 */

/*
    Reference:
        M. E. WIESER AND T. B. COPLEN.
        Atomic weights of the elements 2009 (IUPAC Technical Report).
        Pure Appl. Chem., Vol. 83, No. 2, pp. 359�396, 2011.
    When available, the mass is given for isotope with the longest half-life.
*/
/* 105 dubnium Db */		/* ? Like: Ta */
{ "Db", 268, 268, 268.125000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 106 seaborgium Sg */		/* ? Like: W */
{ "Sg", 271, 271, 271.133000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 107 bohrium Bh */		/* ? Like: Re */
{ "Bh", 267, 267, 267.127700000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 108 hassium Hs */		/* ? Like: Os */
{ "Hs", 277, 277, 277.150000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 109 meitnerium Mt */		/* ? Like: Ir */
{ "Mt", 276, 276, 276.151000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 110 darmstadtium Ds */	/* ? Like: Pt */
{ "Ds", 281, 281, 281.162000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 111 roentgenium Rg */	/* ? Like: Au */
{ "Rg", 280, 280, 280.164000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/* 112 copernicium Cn */	/* ? Like: Hg */
{ "Cn", 285, 285, 285.174000000, METAL ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }},
/*^^^ End of added in v. 1.04 */

#ifdef INCHI_ZFRAG
{ "Zu",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {1,},       {0,},       {0,}       }}, //single bond fragment
{ "Zv",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {2,},       {0,},       {0,}       }}, //double bond fragment
{ "Zw",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {3,},       {0,},       {0,}       }}, //triple bond fragment
{ "Zx",   0,   0,   0.000000000,     0 ,  0,  1, {{0,},       {0,},       {1,2,},     {0,},       {0,}       }}, //aromatic bond fragment
#endif
{ "",     0,   0,   0.000000000,     0 ,  0,  0, {{0,},       {0,},       {0,},       {0,},       {0,}       }},

};
 /*
#ifdef __cplusplus
}
#endif
*/

int ERR_ELEM = 255;
int nElDataLen = sizeof(ElData)/sizeof(ElData[0])-1;

/***********************************************************************************/
int GetElementFormulaFromAtNum(int nAtNum, char *szElement )
{
    nAtNum -= 1;
    if ( 0 < nAtNum )
        nAtNum += 2; /*  bypass D, T */
    if ( 0 <= nAtNum && nAtNum < nElDataLen ) {
        strcpy( szElement, ElData[nAtNum].szElName );
        return 0;
    }
    strcpy( szElement, "??" );
    return -1;
}
/***********************************************************************************/
int get_el_number( const char* elname )
{
    int i;
    const char *p;
    for ( i = 0; (p=ElData[i].szElName)[0] && strcmp( p, elname ); i++ )
        ;
    return p[0]? i : ERR_ELEM;
}
/***********************************************************************************/
int get_periodic_table_number( const char* elname )
{
    int num;
    num = get_el_number( elname );
    if ( num < ERR_ELEM )
        num = inchi_max(1, num-1);
    return num;
}
/***********************************************************************************/
int do_not_add_H( int nPeriodicNum )
{
    return ElData[nPeriodicNum>1? nPeriodicNum+1:0].bDoNotAddH;
}
/***********************************************************************************/
int get_el_valence( int nPeriodicNum, int charge, int val_num )
{
    if ( charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE || val_num >= MAX_NUM_VALENCES )
        return 0;
    return ElData[nPeriodicNum>1? nPeriodicNum+1:0].cValence[NEUTRAL_STATE+charge][val_num];
}
/***********************************************************************************
 *  output valence needed to unumbiguosly reconstruct bonds
 ***********************************************************************************/
int get_unusual_el_valence( int nPeriodicNum, int charge, int radical, int bonds_valence, int num_H, int num_bonds )
{
    int i, num_found, chem_valence, rad_adj, known_chem_valence, exact_found;
    if ( !num_bonds && !num_H )
        return 0;
    if ( charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE ) {
        if ( bonds_valence == num_bonds )
            return 0; /* all single bonds */
        return bonds_valence;
    }
    if ( !get_el_valence( nPeriodicNum, charge, 0 ) && bonds_valence == num_bonds )
        return 0;

    chem_valence = bonds_valence + num_H;
    rad_adj     = 0;
    num_found   = 0;
    exact_found = 0;

    /* take into account radical */
    if (radical==RADICAL_DOUBLET)
        rad_adj = 1;
    else
    if (radical==RADICAL_TRIPLET )
        rad_adj = 2;

    for ( i = 0; i < MAX_NUM_VALENCES; i ++ ) {
        if ( 0 < (known_chem_valence = get_el_valence( nPeriodicNum, charge, i )-rad_adj) &&
             num_bonds <= known_chem_valence && known_chem_valence <= chem_valence ) {
                num_found ++;
            if ( known_chem_valence == chem_valence ) {
                exact_found = 1;
                break;
            }
        }
    }
    return (exact_found && 1 == num_found)? 0 : chem_valence;
}
/***********************************************************************************
 *  output valence needed to unumbiguosly reconstruct number of H
 ***********************************************************************************/
int needed_unusual_el_valence( int nPeriodicNum, int charge, int radical, int bonds_valence,
                               int actual_bonds_valence, int num_H, int num_bonds )
{
    int i, num_found, num_found_known, chem_valence, rad_adj, known_chem_valence, exact_found;
    int num_H_expected;
    char szElement[4];
    /*
    if ( !num_bonds && !num_H )
        return 0;
    */
    if ( num_bonds && !GetElementFormulaFromAtNum(nPeriodicNum, szElement ) ) {
        num_H_expected = get_num_H( szElement, 0, NULL, charge, radical, actual_bonds_valence, 0,0,0,0 );
    } else {
        num_H_expected = num_H;
    }

    chem_valence = bonds_valence + num_H;
    if ( charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE ||
         !get_el_valence( nPeriodicNum, charge, 0 ) ||
         do_not_add_H( nPeriodicNum ) || bonds_valence != actual_bonds_valence ||
         num_H_expected != num_H ) {
        if ( !num_H && !num_H_expected && bonds_valence == actual_bonds_valence )
            return 0; /* no H */
        return chem_valence; /* needs to add H-atoms */
    }

    /* take into account radical */
    if (radical==RADICAL_DOUBLET)
        rad_adj = 1;
    else
    if (radical==RADICAL_TRIPLET )
        rad_adj = 2;
    else
        rad_adj = 0;

    num_found_known = 0;
    num_found       = 0;
    exact_found     = 0;

    for ( i = 0; i < MAX_NUM_VALENCES; i ++ ) {
        if ( 0 <  (known_chem_valence = get_el_valence( nPeriodicNum, charge, i )) &&
             bonds_valence <= (known_chem_valence -= rad_adj) ) {
            /* found known valence that fits without H */
            num_found_known ++;
            if ( known_chem_valence <= chem_valence ) {
                /* known valence is large enough to accommodate (implicit) H */
                num_found ++;
            }
            if ( known_chem_valence == chem_valence ) {
                exact_found = 1;
                break;
            }
        }
    }
    return (exact_found && 1 == num_found && 1 == num_found_known)? 0 : chem_valence? chem_valence : -1 /* needs zero */;
}
/***********************************************************************************
 *  output valence that does not fit any known valences
 ***********************************************************************************/
int detect_unusual_el_valence( int nPeriodicNum, int charge, int radical, int bonds_valence, int num_H, int num_bonds )
{
    int i, chem_valence, rad_adj, known_chem_valence;

    if ( !num_bonds && !num_H )
        return 0;

    if ( charge < MIN_ATOM_CHARGE || charge > MAX_ATOM_CHARGE ) {
        if ( bonds_valence == num_bonds )
            return 0; /* all single bonds */
        return bonds_valence;
    }
    if ( !get_el_valence( nPeriodicNum, charge, 0 ) && bonds_valence == num_bonds )
        return 0;

    chem_valence = bonds_valence + num_H;
    rad_adj     = 0;

    /* take into account radical */
    if (radical==RADICAL_DOUBLET)
        rad_adj = 1;
    else
    if (radical==RADICAL_TRIPLET || radical==RADICAL_SINGLET )
        rad_adj = 2;

    for ( i = 0; i < MAX_NUM_VALENCES; i ++ ) {
        if ( 0 < (known_chem_valence = get_el_valence( nPeriodicNum, charge, i )-rad_adj) ) {
            if ( known_chem_valence == chem_valence ) {
                return 0;
            }
        }
    }
    return chem_valence;
}
/***********************************************************************************/
int get_el_type( int nPeriodicNum )
{
    return ElData[nPeriodicNum+1].nType;
}
/***********************************************************************************/
int is_el_a_metal( int nPeriodicNum )
{
    return 0!=(ElData[nPeriodicNum+1].nType & IS_METAL);
}
/******************************************************************************************************/
/*#ifndef TARGET_API_LIB*/
int extract_ChargeRadical( char *elname, int *pnRadical, int *pnCharge )
{
    char *q, *r, *p;
    int  nCharge=0, nRad = 0, charge_len = 0, k, nVal, nSign, nLastSign=1, len;

    p = elname;

    /*  extract radicals & charges */
    while ( (q = strpbrk( p, "+-^" )) ) {
        switch ( *q ) {
        case '+':
        case '-':
            for ( k = 0, nVal=0; (nSign = ('+' == q[k])) || (nSign = -('-' == q[k])); k++ ) {
                nVal += (nLastSign = nSign);
                charge_len ++;
            }
            if ( (nSign = (int)strtol( q+k, &r, 10 )) ) { /*  fixed 12-5-2001 */
                nVal += nLastSign * (nSign-1);
            }
            charge_len = r - q;
            nCharge += nVal;
            break;
        /* case '.': */ /*  singlet '.' may be confused with '.' in formulas like CaO.H2O */
        case '^':
            nRad = 1; /* doublet here is 1. See below */
            charge_len = 1;
            for ( k = 1; q[0] == q[k]; k++ ) {
                nRad ++;
                charge_len ++;
            }
            break;
        }
        memmove( q, q+charge_len, strlen(q+charge_len)+1 );
    }
    len = (int) strlen(p);
    /*  radical */
    if ( (q = strrchr( p, ':' )) && !q[1]) {
        nRad = RADICAL_SINGLET;
        q[0] = '\0';
        len --;
    } else {
        while( (q = strrchr( p, '.' )) && !q[1] ) {
            nRad ++;
            q[0] = '\0';
            len --;
        }

        nRad = nRad == 1? RADICAL_DOUBLET :
               nRad == 2? RADICAL_TRIPLET : 0;
    }
    *pnRadical = nRad;
    *pnCharge  = nCharge;
    return ( nRad || nCharge );

}
/*#endif*/
/****************************************************************/
int extract_H_atoms( char *elname, S_CHAR num_iso_H[] )
{
    int i, len, c, k, num_H, val;
    char *q;
    i = 0;
    num_H = 0;
    len = (int)strlen(elname);
    c =  UCINT elname[0];
    while ( i < len ) {
        switch ( c ) {
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
        q = elname+i+1; /*  pointer to the next to elname[i] character */
        c =  UCINT q[0];
        if ( k >= 0 && !islower( c ) ) {
            /*  found a hydrogen */
            if ( isdigit( c ) ) {
                val = (int)strtol( q, &q, 10 );
                /*  q = pointer to the next to number of hydrogen atom(s) character */
            } else {
                val = 1;
            }
            if ( k ) {
                num_iso_H[k] += val;
            } else {
                num_H += val;
            }
            /*  remove the hydrogen atom from the string */
            len -= (q-elname)-i;
            memmove( elname+i, q, len + 1 );
            /*  c =  UCINT elname[i]; */
        } else {
            i ++;
        }
        c =  UCINT elname[i]; /*  moved here 11-04-2002 */
    }
    return num_H;
}
/***********************************************************************************/
int get_num_H (const char* elname, int inp_num_H, S_CHAR inp_num_iso_H[],
               int charge, int radical, int chem_bonds_valence, int atom_input_valence,
               int bAliased, int bDoNotAddH, int bHasMetalNeighbor )
{
    int val, i, el_number, num_H = 0, num_iso_H;
    static int el_number_N = 0, el_number_S, el_number_O, el_number_C;
    if ( !el_number_N ) {
        el_number_N = get_el_number( "N" );
        el_number_S = get_el_number( "S" );
        el_number_O = get_el_number( "O" );
        el_number_C = get_el_number( "C" );
    }
    /*  atom_input_valence (cValence) cannot be specified in case of */
    /*  aliased MOLFile atom with known inp_num_H or inp_num_iso_H[] */
    if ( bAliased ) {
        num_H = inp_num_H;
    } else
    if ( atom_input_valence && (atom_input_valence !=15 || chem_bonds_valence) ) {
        num_H = inchi_max( 0, atom_input_valence - chem_bonds_valence );
    } else
    if ( atom_input_valence == 15 && !chem_bonds_valence ) {
        num_H = 0;
    } else
    if ( MIN_ATOM_CHARGE <= charge &&
         MAX_ATOM_CHARGE >= charge &&
         ERR_ELEM != (el_number = get_el_number( elname ) ) &&
         !ElData[el_number].bDoNotAddH && !bDoNotAddH ) {
        /* add hydrogen atoms according to standard element valence */
        if ( radical && radical != RADICAL_SINGLET ) {
            if ( (val = ElData[el_number].cValence[NEUTRAL_STATE+charge][0]) ) {
                val -= (radical==RADICAL_DOUBLET)? 1 :
                       (radical==RADICAL_SINGLET || radical==RADICAL_TRIPLET )? 2 : val;
                /* if unknown radical then do not add H */
                num_H = inchi_max( 0, val - chem_bonds_valence );
            }
        } else {
            /* find the smallest valence that is greater than the sum of the chemical bond valences */
            for ( i = 0;
                  (val=ElData[el_number].cValence[NEUTRAL_STATE+charge][i]) &&
                   val < chem_bonds_valence;
                   i++ )
                ;
            /* special case: do not add H to N(IV), S(III), S+(II), S-(II) */ /* S ions added 2004-05-10 */
            if ( el_number == el_number_N && !charge && !radical && val == 5 )
                val = 3;
            else
            /*
            if ( el_number == el_number_N && !charge && !radical && val == 3 &&
                 chem_bonds_valence == 2 && bHasMetalNeighbor )
                val = 2;
            else
            */
            if ( el_number == el_number_S && !charge && !radical && val == 4 && chem_bonds_valence == 3 )
                val = 3;
            else
            if ( bHasMetalNeighbor && el_number != el_number_C && val > 0 ) {
                val --;
            }
            /*
            if ( (el_number == el_number_S || el_number == el_number_O) &&
                 abs(charge)==1 && !radical && val == 3 && chem_bonds_valence == 2 && bHasMetalNeighbor )
                val = 2;
            else
            */
            num_H = inchi_max( 0, val - chem_bonds_valence );
        }
        num_iso_H = 0;
        if ( inp_num_iso_H ) {
            for ( i = 0; i < NUM_H_ISOTOPES; i ++ ) {
                num_iso_H += inp_num_iso_H[i];
            }
        }
        /*  should not happen because atom here is not aliased */
        if ( num_iso_H ) {
            if ( num_H >= num_iso_H ) {
                num_H -= num_iso_H;
            } else {
                num_H = inp_num_H; /*  as requested in the alias */
                /* num_H = (num_iso_H - num_H) % 2; */ /*  keep unchanged parity of the total number of H atoms */
            }
        }
        /*  should not happen because atom here is not aliased */
        if ( inp_num_H > num_H ) {
            num_H = inp_num_H;  /*  as requested in the alias */
            /* num_H = inp_num_H + (inp_num_H - num_H)%2; */ /*  keep unchanged parity of the number of non-isotopic H atoms */
        }
    } else {
        num_H = inp_num_H;
    }
    return num_H;
}
/***********************************************************************************/
int get_atw_from_elnum( int nAtNum )
{
    nAtNum -= 1;
    if ( 0 < nAtNum )
        nAtNum += 2; /*  bypass D, T */
    if ( 0 <= nAtNum && nAtNum < nElDataLen ) {
        return (int)ElData[nAtNum].nAtMass;
    }
    return 0;
}
/***********************************************************************************/
/*
int get_mw(char elname[])
{
    int i;

    for (i=0; i<NUMEL; i++)
        if (strcmp(elname,elements[i])==0)
            return(atomic_wt[i]);
    return(0);
}
*/
/***********************************************************************************/
#ifndef TARGET_API_LIB
/***********************************************************************************/
int get_atw(const char *elname)
{
    int el_number, atw;
    if ( ERR_ELEM != (el_number = get_el_number( elname )) ) {
        atw = ElData[el_number].nAtMass;
    } else {
        atw = 0;
    }
    return atw;
}
/***********************************************************************************/
int normalize_name( char* name )
{
    /* remove leading & trailing spaces; replace consecutive spaces with a single space */
    /* Treat non-printable characters (Greeks) as spaces. 11-23-99 DCh. */
    int i, len, n;
    len = (int)strlen(name);
    for ( i = 0, n = 0; i < len; i++ ) {
        if ( isspace( UCINT name[i] ) /*|| !isprint( UCINT name[i] )*/ ) {
            name[i] = ' '; /* exterminate tabs !!! */
            n++;
        } else {
            if ( n > 0 ) {
                memmove( (void*) &name[i-n], (void*) &name[i], len-i+1 );
                i   -= n;
                len -= n;
            }
            n = -1;
        }
    }
    if ( n == len ) /* empty line */
        name[len=0] = '\0';
    else
    if ( ++n && n <= len ) {
        len -= n;
        name[len] = '\0';
    }
    return len;
}
#endif /* ifndef TARGET_API_LIB */
/************************************************/
#ifndef inchi_malloc
void *inchi_malloc(size_t c)
{
    return  malloc(c);
}
#endif
#ifndef inchi_calloc
void *inchi_calloc(size_t c, size_t n)
{
    return calloc(c,n);
}
#endif
#ifndef inchi_free
void inchi_free(void *p)
{
    if(p) {
        free(p); /*added check if zero*/
    }
}
#endif





#ifndef TARGET_API_LIB
/*************************************************************************/
void remove_trailing_spaces( char* p )
{
    int   len;
    for( len = (int)strlen( p ) - 1; len >= 0 && isspace( UCINT p[len] ); len-- )
        ;
    p[++len] = '\0';
}
/*************************************************************************/
void remove_one_lf( char* p)
{
    size_t len;
    if ( p && 0 < (len = strlen(p)) && p[len-1] == '\n' ){
        p[len-1] = '\0';
        if ( len >= 2 && p[len-2] == '\r' )
            p[len-2] = '\0';
    }
}
#endif /* ifndef TARGET_API_LIB */



/***************************************************************************/
/* Copies up to maxlen characters INCLUDING end null from source to target */
/* Fills out the rest of the target with null bytes */
int mystrncpy(char *target,const char *source,unsigned maxlen)
{   /*  protected from non-zero-terminated source and overlapped target/source. 7-9-99 DCh. */
    const char  *p;
    unsigned    len;

    if (target==NULL || maxlen == 0 || source == NULL)
        return 0;
    if ( (p = (const char*)memchr(source, 0, maxlen)) ) {
        len = p-source; /*  maxlen does not include the found zero termination */
    } else {
        len = maxlen-1; /*  reduced length does not include one more byte for zero termination */
    }
    if ( len )
        memmove( target, source, len );
    /* target[len] = '\0'; */
    memset( target+len, 0, maxlen-len); /*  zero termination */
    return 1;
}
/************************************************************************/
/* Remove leading and trailing white spaces                             */
char* LtrimRtrim( char *p, int* nLen )
{
    int i, len=0;
    if ( p &&  (len = (int) strlen( p )) ) {
        for ( i = 0; i < len && __isascii( p[i] ) && isspace( p[i] ); i++ )
            ;
        if ( i )
            (memmove)( p, p+i, (len -= i)+1 );
        for ( ; 0 < len && __isascii( p[len-1] ) && isspace( p[len-1] ); len-- )
            ;
        p[len] = '\0';
    }
    if ( nLen )
        *nLen = len;
    return p;
}
/*************************************************************************/
AT_NUMB *is_in_the_list( AT_NUMB *pathAtom, AT_NUMB nNextAtom, int nPathLen )
{
    for ( ; nPathLen && *pathAtom != nNextAtom; nPathLen--,  pathAtom++ )
        ;
    return nPathLen? pathAtom : NULL;
}
/******************************************************************************************************/
int nBondsValToMetal( inp_ATOM* at, int iat )
{
    int i, neigh, bond_type, nVal2Metal = 0;
    inp_ATOM* a  = at + iat;
    for ( i = 0; i < a->valence; i ++ ) {
        neigh = a->neighbor[i];
        if ( is_el_a_metal( at[(int)a->neighbor[i]].el_number ) ) {
            bond_type = a->bond_type[i];
            if ( bond_type <= BOND_TYPE_TRIPLE ) {
                nVal2Metal += bond_type;
            } else {
                return -1;  /* bond to metal order is not well defined */
            }
        }
    }
    return nVal2Metal;
}
/************************************************************************/
int num_of_H( inp_ATOM *at, int iat )
{
    static int el_number_H;
    int    i, n, num_explicit_H = 0;
    inp_ATOM *a = at + iat;
    if ( !el_number_H )
        el_number_H = get_periodic_table_number( "H" );
    for ( i = 0; i < a->valence; i ++ ) {
        n = a->neighbor[i];
        num_explicit_H += ( 1 == at[n].valence && el_number_H == at[n].el_number );
    }
    return num_explicit_H+NUMH(at,iat);
}
/************************************************************************/
int has_other_ion_neigh( inp_ATOM *at, int iat, int iat_ion_neigh, const char *el, int el_len )
{
    int charge = at[iat_ion_neigh].charge;
    int i, neigh;
    for ( i = 0; i < at[iat].valence; i ++ ) {
        neigh = at[iat].neighbor[i];
        if ( neigh != iat_ion_neigh && at[neigh].charge == charge &&
             NULL != memchr( el, at[neigh].el_number, el_len ) ) {
            return 1;
        }
    }
    return 0;
}
/************************************************************************/
/* BFS r=2 */
int has_other_ion_in_sphere_2(inp_ATOM *at, int iat, int iat_ion_neigh, const char *el, int el_len )
{
#define MAXQ 16
    AT_NUMB q[MAXQ];
    int lenq=0, lenq2, dist = 0, i = 0, iq, neigh, j, nRet=0;
    q[lenq++] = iat;
    at[iat].cFlags = 1;

    iq  = 0;
    dist = 1;
    /* use at->cFlags as an indicator */
    while ( dist <= 2 ) {
        for ( lenq2 = lenq; iq < lenq2; iq ++ ) {
            i = q[iq];
            for ( j = 0; j < at[i].valence; j ++ ) {
                neigh = at[i].neighbor[j];
                if ( !at[neigh].cFlags &&
                     at[neigh].valence <= 3 &&
                     NULL != memchr( el, at[neigh].el_number, el_len ) ) {
                    q[lenq ++] = neigh;
                    at[neigh].cFlags = 1;
                    if ( neigh != iat_ion_neigh &&
                         at[iat_ion_neigh].charge == at[neigh].charge ) {
                        nRet ++;
                    }
                }
            }
        }
        dist ++;
    }
    for ( iq = 0; iq < lenq; iq ++ ) {
        i = q[iq];
        at[i].cFlags = 0;
    }
    return nRet;
}
/************************************************************************/
int nNoMetalNumBonds( inp_ATOM *at, int at_no )
{
    inp_ATOM *a = at + at_no;
    int num_H = NUMH(a, 0);
    int std_chem_bonds_valence = get_el_valence( a->el_number, a->charge, 0 );
    int i;
    if ( a->chem_bonds_valence + num_H > std_chem_bonds_valence ) {
        int valence_to_metal = 0;
        int num_bonds_to_metal = 0;
        for ( i = 0; i < a->valence; i ++ ) {
            if ( is_el_a_metal( at[(int)a->neighbor[i]].el_number ) ) {
                if ( (a->bond_type[i] & BOND_TYPE_MASK) >= BOND_TYPE_ALTERN ) {
                    return a->valence; /* fall back */
                }
                num_bonds_to_metal ++;
                valence_to_metal += (a->bond_type[i] & BOND_TYPE_MASK);
            }
        }
        if ( a->chem_bonds_valence + num_H - valence_to_metal == std_chem_bonds_valence ) {
            /* removing bonds to metal produces standard valence */
            return a->valence - num_bonds_to_metal;
        }
    }
#if ( S_VI_O_PLUS_METAL_FIX_BOND == 1 )
    else
    if ( 1 == a->charge && 2 == get_endpoint_valence(a->el_number) &&
         a->chem_bonds_valence + num_H == std_chem_bonds_valence ) {
        int valence_to_metal = 0;
        int num_bonds_to_metal = 0;
        for ( i = 0; i < a->valence; i ++ ) {
            if ( is_el_a_metal( at[(int)a->neighbor[i]].el_number ) ) {
                if ( (a->bond_type[i] & BOND_TYPE_MASK) >= BOND_TYPE_ALTERN ) {
                    return a->valence; /* fall back */
                }
                num_bonds_to_metal ++;
                valence_to_metal += (a->bond_type[i] & BOND_TYPE_MASK);
            }
        }
        if ( 1 == valence_to_metal ) {
            /* removing bonds to metal produces standard valence */
            return a->valence - num_bonds_to_metal;
        }
    }
#endif

    return a->valence;
}
/************************************************************************/
int nNoMetalBondsValence( inp_ATOM *at, int at_no )
{
    inp_ATOM *a = at + at_no;
    int num_H = NUMH(a, 0);
    int std_chem_bonds_valence = get_el_valence( a->el_number, a->charge, 0 );
    int i;
    if ( a->chem_bonds_valence + num_H > std_chem_bonds_valence ) {
        int valence_to_metal = 0;
        /*int num_bonds_to_metal = 0;*/
        for ( i = 0; i < a->valence; i ++ ) {
            if ( is_el_a_metal( at[(int)a->neighbor[i]].el_number ) ) {
                if ( (a->bond_type[i] & BOND_TYPE_MASK) >= BOND_TYPE_ALTERN ) {
                    return a->valence; /* fall back */
                }
                /*num_bonds_to_metal ++;*/
                valence_to_metal += (a->bond_type[i] & BOND_TYPE_MASK);
            }
        }
        if ( a->chem_bonds_valence + num_H - valence_to_metal == std_chem_bonds_valence ) {
            /* removing bonds to metal produces standard valence */
            return a->chem_bonds_valence - valence_to_metal;
        }
    }
#if ( S_VI_O_PLUS_METAL_FIX_BOND == 1 )
    else
    if ( 1 == a->charge && 2 == get_endpoint_valence(a->el_number) &&
         a->chem_bonds_valence + num_H == std_chem_bonds_valence ) {
        int valence_to_metal = 0;
        /*int num_bonds_to_metal = 0;*/
        for ( i = 0; i < a->valence; i ++ ) {
            if ( is_el_a_metal( at[(int)a->neighbor[i]].el_number ) ) {
                if ( (a->bond_type[i] & BOND_TYPE_MASK) >= BOND_TYPE_ALTERN ) {
                    return a->valence; /* fall back */
                }
                /*num_bonds_to_metal ++;*/
                valence_to_metal += (a->bond_type[i] & BOND_TYPE_MASK);
            }
        }
        if ( 1 == valence_to_metal ) {
            /* removing bonds to metal produces standard valence */
            return a->chem_bonds_valence - valence_to_metal;
        }
    }
#endif
    return a->chem_bonds_valence;
}
/************************************************************************/
int nNoMetalNeighIndex( inp_ATOM *at, int at_no )
{
    inp_ATOM *a = at + at_no;
    int i;
    for ( i = 0; i < a->valence; i ++ ) {
        if ( !is_el_a_metal( at[(int)a->neighbor[i]].el_number ) ) {
            return i;
        }
    }
    return -1;
}
/************************************************************************/
int nNoMetalOtherNeighIndex( inp_ATOM *at, int at_no, int cur_neigh )
{
    inp_ATOM *a = at + at_no;
    int i, neigh;
    for ( i = 0; i < a->valence; i ++ ) {
        neigh = (int)a->neighbor[i];
        if ( neigh != cur_neigh && !is_el_a_metal( at[neigh].el_number ) ) {
            return i;
        }
    }
    return -1;
}
/************************************************************************/
int nNoMetalOtherNeighIndex2( inp_ATOM *at, int at_no, int cur_neigh, int cur_neigh2 )
{
    inp_ATOM *a = at + at_no;
    int i, neigh;
    for ( i = 0; i < a->valence; i ++ ) {
        neigh = (int)a->neighbor[i];
        if ( neigh != cur_neigh && neigh != cur_neigh2 && !is_el_a_metal( at[neigh].el_number ) ) {
            return i;
        }
    }
    return -1;
}


#ifndef COMPILE_ANSI_ONLY
/**************************************************************************/
int MakeRemovedProtonsString( int nNumRemovedProtons, NUM_H *nNumExchgIsotopicH, NUM_H *nNumRemovedProtonsIsotopic,
                              int bIsotopic, char *szRemovedProtons, int *num_removed_iso_H )
{
    int i, j, len, num;
    len = 0;
    if ( nNumRemovedProtons ) {
        len = sprintf ( szRemovedProtons, "Proton balance: %c %d H+",
                        nNumRemovedProtons>=0? '+':'-', abs(nNumRemovedProtons) );
    }
    if ( bIsotopic && (nNumRemovedProtonsIsotopic || nNumExchgIsotopicH) ) {
        for ( i = 0, j = 0; i < NUM_H_ISOTOPES; i ++ ) {
            num = (nNumExchgIsotopicH? nNumExchgIsotopicH[i]:0) +
                  (nNumRemovedProtonsIsotopic? nNumRemovedProtonsIsotopic[i]:0);
            if ( num ) {
                len += sprintf( szRemovedProtons+len, "%s %d^%dH",
                                j? ", ":"  [ removed ", num, i+1);
                j ++;
            }
        }
        if ( j ) {
            len += sprintf( szRemovedProtons+len, " ]" );
            if ( num_removed_iso_H )
                *num_removed_iso_H = j;
        }
    }
    if ( !len ) {
        szRemovedProtons[0] = '\0';
    }
    return len;
}
#endif

/*
    According to
    http://info-uri.info/registry/OAIHandler?verb=GetRecord&metadataPrefix=reg&identifier=info:inchi/

    An InChI identifier may contain the following characters:

    A-Z
    a-z
    0-9
    ()*+,-./;=?@


    Here we consider any character not conforming this specification as a whitespace
    which marks the end of the InChI string.
    For example:
    "InChI=1/Ar%"
    "InChI=1/Ar\n"
    "InChI=1/Ar\r\t"
    all will be trimmed to
    "InChI=1/Ar"

*/

/**************************************************************************/
void extract_inchi_substring(char ** buf, const char *str, size_t slen)
{
size_t i;
char *p, pp;

    *buf = NULL;
    if (str==NULL)
        return;
    if (strlen(str)<1)
        return;

    p = strstr(str, "InChI=");
    if (NULL==p)
        return;

    for (i=0; i<slen; i++)
    {
        pp = p[i];

        if (pp >= 'A' && pp <='Z')   continue;
        if (pp >= 'a' && pp <='z')   continue;
        if (pp >= '0' && pp <='9')   continue;
        switch ( pp )
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
            case '@': continue;

            default:            break;
        }
        break;
    }

    *buf = (char*) inchi_calloc(i+1, sizeof(char));
    memcpy(*buf, p, i);
    (*buf)[i] = '\0';

    return;
}





#ifdef COMPILE_ANSI_ONLY
/*************************************************************************/
/*************          non-ANSI functions                ****************/
/*************************************************************************/
#define __MYTOLOWER(c) ( ((c) >= 'A') && ((c) <= 'Z') ? ((c) - 'A' + 'a') : (c) )

#if ( defined(COMPILE_ADD_NON_ANSI_FUNCTIONS) || defined(__STDC__) && __STDC__ == 1 )
/* support (VC++ Language extensions) = OFF && defined(COMPILE_ANSI_ONLY) */
int memicmp ( const void * p1, const void * p2, size_t length )
{
    const U_CHAR *s1 = (const U_CHAR*)p1;
    const U_CHAR *s2  = (const U_CHAR*)p2;
    while ( length-- ) {
        if ( *s1 == *s2 ||
              __MYTOLOWER( (int)*s1 ) == __MYTOLOWER( (int)*s2 )) {
            s1 ++;
            s2  ++;
        } else {
            return __MYTOLOWER( (int)*s1 ) - __MYTOLOWER( (int)*s2 );
        }
    }
    return 0;
}
/*************************************************************************/
int stricmp( const char *s1, const char *s2 )
{
    while ( *s1 ) {
        if ( *s1 == *s2 ||
              __MYTOLOWER( (int)*s1 ) == __MYTOLOWER( (int)*s2 )) {
            s1 ++;
            s2  ++;
        } else {
            return __MYTOLOWER( (int)*s1 ) - __MYTOLOWER( (int)*s2 );
        }
    }
    if ( *s2 )
        return -1;
    return 0;
}
/*************************************************************************/
char *_strnset( char *s, int val, size_t length )
{
    char *ps = s;
    while (length-- && *ps)
        *ps++ = (char)val;
    return s;
}
/*************************************************************************/
char *_strdup( const char *string )
{
    char *p = NULL;
    if ( string ) {
        size_t length = strlen( string );
        p = (char *) inchi_malloc( length + 1 );
        if ( p ) {
            strcpy( p, string );
        }
    }
    return p;
}
#endif
#endif
