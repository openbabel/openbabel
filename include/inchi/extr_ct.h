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


#ifndef _EXTR_CT_H_
#define _EXTR_CT_H_



#include "mode.h"
#include "ichisize.h"
#include "util.h"



struct AtData
{
    char element[3];
    int maxvalence;
};



#define NUM_CHEM_ELEMENTS 127               /* well above number of known chem. elements    */



#define         AT_ISO_SORT_KEY_MULT  32    /* up to 32 identical hydrogen isotopes         */
                                            /* (similar to T_GROUP_ISOWT_MULT)              */
                                            /* changed from 16 9-12-2003                    */

typedef long    AT_ISO_SORT_KEY;            /* signed, should hold up to 4096*max_iso_diff  */
                                            /* (similar to T_GROUP_ISOWT)                   */
/*
   = num_1H + AT_ISO_SORT_KEY_MULT*(num_D + AT_ISO_SORT_KEY_MULT*(num_T+AT_ISO_SORT_KEY_MULT*iso_atw_diff))
*/

/* typedef signed   char AT_ISOTOPIC; */ /* + or - */
typedef struct tagStereoCarb
{
    AT_NUMB at_num;
    U_CHAR  parity;
} AT_STEREO_CARB;

typedef struct tagStereoDble
{
    AT_NUMB at_num1;
    AT_NUMB at_num2;
    U_CHAR  parity;
} AT_STEREO_DBLE;

typedef struct tagIsotopicAtom
{
    AT_NUMB at_num;
    NUM_H   num_1H;
    NUM_H   num_D;
    NUM_H   num_T;
    NUM_H   iso_atw_diff;
} AT_ISOTOPIC;

typedef AT_NUMB AT_STEREO;



#define BYTE_BITS 8                     /* number of bits in one byte                               */
#define BOND_MASK 0xf                   /* 4 bits                                                   */
#define BOND_BITS 4                     /* 3 or 4 does not matter; 2 is too small for BOND_TAUTOM   */
#define BOND_ADD  (BOND_BITS==2?-1:0)   /* subtract 1 from bonds stored in CT                       */



typedef struct tagAtom
{
    char elname[ATOM_EL_LEN];
    AT_NUMB    neighbor[MAXVAL];    /* changed to unsigned 2-2-95. D.Ch.                            */
    AT_NUMB    init_rank;           /* also used in remove_terminal_HDT() to save                   */
                                    /* orig. at. number                                             */
    AT_NUMB    orig_at_number;
    AT_NUMB    orig_compt_at_numb;
    /* low 3 bits=bond type;                                                                        */
    /* high 5 bits (in case of cut-vertex atom) = an attached part number                           */
    U_CHAR bond_type[MAXVAL];
    U_CHAR el_number;               /* periodic table number = charge of the nucleus =              */
                                    /* number of the protons                                        */
    /* U_CHAR hill_type; */         /* number in pseudo Hill order                                  */
    S_CHAR valence;
    S_CHAR chem_bonds_valence;      /* 8-24-00 to treat tautomer centerpoints, etc.                 */
    S_CHAR num_H;                   /* first not including D, T; add_DT_to_num_H() includes.        */
    S_CHAR num_iso_H[NUM_H_ISOTOPES];   /* num 1H, 2H(D), 3H(T)                                     */
    S_CHAR cFlags;
    S_CHAR iso_atw_diff;            /* abs(iso_atw_diff) < 127 or 31 - ???                          */
    AT_ISO_SORT_KEY iso_sort_key;   /* = num_1H + AT_ISO_SORT_KEY_MULT^1*num_D                      */
                                    /*         + AT_ISO_SORT_KEY_MULT^2*num_T                       */
                                    /*         + AT_ISO_SORT_KEY_MULT^3*iso_atw_diff                */
    S_CHAR charge;
    S_CHAR radical;                 /* 1=>doublet(.), 2=> triplet as singlet (:)                    */
                                    /* ???? why are they same ????                                  */
    S_CHAR marked;

    AT_NUMB endpoint;               /* tautomer analysis. If != 0 then the hydrogens & (-)charge    */
                                    /* are in the tautomer group.                                   */

    /*
       Pairs stereo_bond_neighbor[] and stereo_bond_neighbor2[], etc
       initially refer to non-isotopic and isotopic cases, respectively.
       To use same stereo processing code these arrays are swapped when
       switching from non-isotopic to isotopic processing and back.
    */
    AT_NUMB stereo_bond_neighbor[MAX_NUM_STEREO_BONDS]; /* Original number of an opposite atom      */
    AT_NUMB stereo_bond_neighbor2[MAX_NUM_STEREO_BONDS];/*     (stereo bond neighbor) +1;           */
    S_CHAR  stereo_bond_ord[MAX_NUM_STEREO_BONDS];      /* Ordering number of a bond/neighbor       */
                                                        /* in the direction to the                  */
    S_CHAR  stereo_bond_ord2[MAX_NUM_STEREO_BONDS];     /* stereo bond opposite atom (important     */
                                                        /* for cumulenes);                          */
    S_CHAR  stereo_bond_z_prod[MAX_NUM_STEREO_BONDS];   /* Relative  atom-neighbors                 */
    S_CHAR  stereo_bond_z_prod2[MAX_NUM_STEREO_BONDS];  /* double bond planes orientation;          */
    S_CHAR  stereo_bond_parity[MAX_NUM_STEREO_BONDS];   /* parity + MULT_STEREOBOND*chain_length,   */
    S_CHAR  stereo_bond_parity2[MAX_NUM_STEREO_BONDS];  /* where:                                   */
                    /*
                        parity (Mask 0x07=BITS_PARITY):

                         0   = AB_PARITY_NONE = not a stereo bond
                         1/2 = AB_PARITY_ODD/EVEN = bond parity defined from initial ranks
                         3   = AB_PARITY_UNKN = geometry is unknown to the user
                         4   = AB_PARITY_UNDF = not enough geometry info to find the parity
                         6   = AB_PARITY_CALC = calculate later from the neighbor ranks;
                               some ot them can be replaced with AB_PARITY_ODD/EVEN
                               after equivalence ranks have been determined

                         length (Mask 0x38=MASK_CUMULENE_LEN, length=stereo_bond_parity[i]/MULT_STEREOBOND):

                         0   => double or alternating stereogenic bond
                         1   => cumulene with 2 double bonds (stereogenic center)
                         2   => cumulene with 3 double bonds (stereogenic bond)
                         length <= (MAX_CUMULENE_LEN=2)
                         bit KNOWN_PARITIES_EQL =  0x40: all pairs of const. equ. atoms are connected
                         by stereo bonds and these bonds have identical parities
                    */

    S_CHAR parity;  /* -- Mask 0x07=BITS_PARITY: --
                        0 = AB_PARITY_NONE => no parity; also parity&0x38 = 0
                        1 = AB_PARITY_ODD  => odd parity
                        2 = AB_PARITY_EVEN => even parity
                        3 = AB_PARITY_UNKN => user marked as unknown parity
                        4 = AB_PARITY_UNDF => parity cannot be defined because of symmetry
                            or not well defined geometry
                    */
    S_CHAR parity2; /* parity including parity due to isotopic terminal H */
    /* bit msks: 0x07 => known parity (1,2,3,4) or AB_PARITY_CALC=6, AB_PARITY_IISO = 6             */
    /*           0x40 => KNOWN_PARITIES_EQL                                                         */
    S_CHAR stereo_atom_parity;  /* similar to stereo_bond_parity[]: known in
                                /( advance AB_PARITY_* value + KNOWN_PARITIES_EQL bit               */
    S_CHAR stereo_atom_parity2;
    S_CHAR final_parity;        /* defined by equivalence ranks                                          */
    S_CHAR final_parity2;       /* defined by equivalence ranks, incl. due to terminal isotopic H        */
    S_CHAR bAmbiguousStereo;
    S_CHAR bHasStereoOrEquToStereo;
    S_CHAR bHasStereoOrEquToStereo2;

#if ( FIND_RING_SYSTEMS == 1 )
    S_CHAR    bCutVertex;
    AT_NUMB   nRingSystem;
    AT_NUMB   nNumAtInRingSystem;
    AT_NUMB   nBlockSystem;

#if ( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    AT_NUMB   nDistanceFromTerminal;
#endif

#endif

    S_CHAR z_dir[3];
} sp_ATOM;

#define BOND_SINGLE BOND_TYPE_SINGLE  /* 1                                  */
#define BOND_DOUBLE BOND_TYPE_DOUBLE  /* 2                                  */
#define BOND_TRIPLE BOND_TYPE_TRIPLE  /* 3                                  */
#define BOND_ALTERN BOND_TYPE_ALTERN  /* 4 single/double                    */

#define BOND_ALT_123             5   /* single/double/triple                */
#define BOND_ALT_13              6   /* single/triple                       */
#define BOND_ALT_23              7   /* double/triple                       */
#define BOND_TAUTOM              8
#define BOND_ALT12NS             9
#define BOND_NUMDIF              9  /* number of different kinds of bonds  */

#define BOND_TYPE_MASK        0x0f

#define BOND_MARK_ALL         0xf0  /* complement to BOND_TYPE_MASK         */

#define BOND_MARK_ALT12       0x10
#define BOND_MARK_ALT123      0x20
#define BOND_MARK_ALT13       0x30
#define BOND_MARK_ALT23       0x40
#define BOND_MARK_ALT12NS     0x50  /* 1 or 2, non-stereo                   */
#define BOND_MARK_MASK        0x70

#define ACTUAL_ORDER(PBNS, IAT,IBOND, BTYPE)  ( ((PBNS) && (PBNS)->edge && (PBNS)->vert &&\
    ((BTYPE)==BOND_ALT_123 || (BTYPE)==BOND_ALT_13 || (BTYPE)==BOND_ALT_23))? (PBNS)->edge[(PBNS)->vert[IAT].iedge[IBOND]].flow+BOND_TYPE_SINGLE:(BTYPE))


#define BITS_PARITY        0x07  /* mask to retrieve half-bond parity */
#define MASK_CUMULENE_LEN  0x38  /* mask to retrieve (cumulene chain length - 1)*MULT_STEREOBOND    */
#define KNOWN_PARITIES_EQL 0x40  /* parity is same for all pairs of constit. equivalent atoms       */
#define MAX_CUMULENE_LEN   2     /* max number of bonds in a cumulene chain - 1                     */

#define MULT_STEREOBOND    0x08  /* multiplier for cumulene chain length
                                    odd length => chiral, even length => stereogenic bond           */

#define MAKE_BITS_CUMULENE_LEN(X)   ((X)*MULT_STEREOBOND)
#define GET_BITS_CUMULENE_LEN(X)    ((X)&MASK_CUMULENE_LEN)
#define BOND_CHAIN_LEN(X)           (GET_BITS_CUMULENE_LEN(X)/MULT_STEREOBOND) /* 0 => double bond, 1 => allene, 2 => cumulene,..*/
#define IS_ALLENE_CHAIN(X)          ((GET_BITS_CUMULENE_LEN(X)/MULT_STEREOBOND)%2)

/* atom or bond parity value definitions */
#define AB_PARITY_NONE   0  /* 0 => no parity; also parity&0x38 = 0 */
#define AB_PARITY_ODD    1  /* 1 => odd parity */
#define AB_PARITY_EVEN   2  /* 2 => even parity */
#define AB_PARITY_UNKN   3  /* 3 => user marked as unknown parity */
#define AB_PARITY_UNDF   4  /* 4 => parity cannot be defined because of symmetry or not well defined geometry */
#define AB_PARITY_IISO   5  /* 5 => no parity because of identical atoms */
#define AB_PARITY_CALC   6  /* 6 => calculate parity later */
#define AB_PARITY_0D     8  /* 8 => bit signifies 0D case -- not used */

#define AB_INV_PARITY_BITS (AB_PARITY_ODD ^ AB_PARITY_EVEN)


#define AB_MAX_KNOWN_PARITY        4 /* precalculated from const. equivalence parities */
#define AB_MIN_KNOWN_PARITY        1

#define AB_MAX_PART_DEFINED_PARITY 3 /* 1, 2, 3 => defined parities, uncluding 'unknown' */
#define AB_MIN_PART_DEFINED_PARITY 1 /* min(AB_PARITY_ODD, AB_PARITY_EVEN, AB_PARITY_UNKN) */

#define AB_MAX_WELL_DEFINED_PARITY 2 /* 1, 2 => well defined parities, uncluding 'unknown' */
#define AB_MIN_WELL_DEFINED_PARITY 1 /* min(AB_PARITY_ODD, AB_PARITY_EVEN) */

#define AB_MIN_ILL_DEFINED_PARITY  3
#define AB_MAX_ILL_DEFINED_PARITY  4

#define AB_MAX_ANY_PARITY          4
#define AB_MIN_ANY_PARITY          1

#define AMBIGUOUS_STEREO           1
#define AMBIGUOUS_STEREO_ATOM      2
#define AMBIGUOUS_STEREO_BOND      4
#define AMBIGUOUS_STEREO_ATOM_ISO  8
#define AMBIGUOUS_STEREO_BOND_ISO  16
#define AMBIGUOUS_STEREO_ERROR     32


#define MIN_DOT_PROD 50          /* min value of at->stereo_bond_z_prod[i] to define parity */

#define ATOM_PARITY_VAL(X)          (X)
#define ATOM_PARITY_PART_DEF(X)     (AB_MIN_PART_DEFINED_PARITY <= (X) && (X) <= AB_MAX_PART_DEFINED_PARITY)
#define ATOM_PARITY_ILL_DEF(X)      (AB_MIN_ILL_DEFINED_PARITY <= (X) && (X) <= AB_MAX_ILL_DEFINED_PARITY)
#define ATOM_PARITY_KNOWN(X)        (AB_MIN_KNOWN_PARITY <= (X) && (X) <= AB_MAX_KNOWN_PARITY)
#define ATOM_PARITY_WELL_DEF(X)     (AB_MIN_WELL_DEFINED_PARITY <= (X) && (X) <= AB_MAX_WELL_DEFINED_PARITY)
#define ATOM_PARITY_NOT_UNKN(X)     (ATOM_PARITY_KNOWN(X) && (X) != AB_PARITY_UNKN)

#define PARITY_VAL(X)               ((X) & BITS_PARITY)
#define PARITY_PART_DEF(X)          (AB_MIN_PART_DEFINED_PARITY <= PARITY_VAL(X) && PARITY_VAL(X) <= AB_MAX_PART_DEFINED_PARITY)
#define PARITY_ILL_DEF(X)           (AB_MIN_ILL_DEFINED_PARITY <= PARITY_VAL(X) && PARITY_VAL(X) <= AB_MAX_ILL_DEFINED_PARITY)
#define PARITY_KNOWN(X)             (AB_MIN_KNOWN_PARITY <= PARITY_VAL(X) && PARITY_VAL(X) <= AB_MAX_KNOWN_PARITY)
#define PARITY_WELL_DEF(X)          (AB_MIN_WELL_DEFINED_PARITY <= PARITY_VAL(X) && PARITY_VAL(X) <= AB_MAX_WELL_DEFINED_PARITY)
#define PARITY_CALCULATE(X)         (AB_PARITY_CALC == PARITY_VAL(X))
#define BOND_PARITY_PART_DEFINED(X) (PARITY_PART_DEF(X) || PARITY_CALCULATE(X))
#define BOND_PARITY_PART_KNOWN(X)   (PARITY_KNOWN(X) || PARITY_CALCULATE(X))
#define ALL_BUT_PARITY(X)           ((X)&~BITS_PARITY)

#define ALWAYS_SET_STEREO_PARITY           0
#define NO_ISOLATED_NON_6RING_AROM_BOND    0  /* for Yuri */
#define SAVE_6_AROM_CENTERS                0  /* for Yuri */



#endif    /* _EXTR_CT_H_ */
