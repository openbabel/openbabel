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


#ifndef _ICHI_H_
#define _ICHI_H_


#include "incomdef.h"
#include "mode.h"


#define REQ_MODE_BASIC              0x000001    /* B include Fixed-H layer */
#define REQ_MODE_TAUT               0x000002    /* T include Mobile-H layer */
#define REQ_MODE_ISO                0x000004    /* I    */
#define REQ_MODE_NON_ISO            0x000008    /* NI   */
#define REQ_MODE_STEREO             0x000010    /* S    */
#define REQ_MODE_ISO_STEREO         0x000020    /* IS   */
#define REQ_MODE_NOEQ_STEREO        0x000040    /* SS   */
#define REQ_MODE_REDNDNT_STEREO     0x000080    /* RS   */
#define REQ_MODE_NO_ALT_SBONDS      0x000100    /* NASB */

/* new 10-10-2003 */
#define REQ_MODE_RELATIVE_STEREO    0x000200    /* REL All Relative Stereo */

#define REQ_MODE_RACEMIC_STEREO     0x000400    /* RAC All Racemic Stereo */

#define REQ_MODE_SC_IGN_ALL_UU      0x000800    /* IAUSC Ignore stereocenters if All Undef/Unknown */

#define REQ_MODE_SB_IGN_ALL_UU      0x001000    /* IAUSC Ignore stereobonds if All Undef/Unknown */

#define REQ_MODE_CHIR_FLG_STEREO    0x002000    /* SUCF  If Chiral flag then Abs otherwise Rel stereo */
/* end of 10-10-2003 */



/* 2009-12-05 */
#define REQ_MODE_DIFF_UU_STEREO     0x004000    /* SLUUD Make labels for unknown and undefined stereo different */
/* 2009-12-05 */

#define REQ_MODE_MIN_SB_RING_MASK   0x0F0000    /* RSB  */
#define REQ_MODE_MIN_SB_RING_SHFT      16

#define REQ_MODE_DEFAULT  (REQ_MODE_BASIC | REQ_MODE_TAUT | REQ_MODE_ISO | REQ_MODE_NON_ISO | REQ_MODE_STEREO)

#define WARN_FAILED_STEREO          0x0001
#define WARN_FAILED_ISOTOPIC        0x0002
#define WARN_FAILED_ISOTOPIC_STEREO 0x0004
#define ERR_NO_CANON_RESULTS        0x0008

/*********** compare components flags **********************************/
#define CMP_COMPONENTS              0x0001     /* perform compare components */
#define CMP_COMPONENTS_NONISO       0x0002     /* ignore isotopic */
#define CMP_COMPONENTS_NONTAUT      0x0004     /* compare non-tautomeric */

/****************** chemical identifier member definitions *************/
typedef struct tagINChI_IsotopicAtom
{

    AT_NUMB   nAtomNumber;      /* Canonical atom number                        */

    NUM_H     nIsoDifference;   /* 0=non-isotopic; 1=rounded avg. atomic mass   */

    NUM_H     nNum_H;           /* number of 1H isotopic atoms attached         */

    NUM_H     nNum_D;           /* number of 2H isotopic atoms attached         */

    NUM_H     nNum_T;           /* number of 3H isotopic atoms attached         */
} INChI_IsotopicAtom;



typedef struct tagINChI_IsotopicTGroup
{
    AT_NUMB   nTGroupNumber;    /* Tautomeric group number          */
    AT_NUMB   nNum_H;           /* number of 1H isotopic atoms      */
    AT_NUMB   nNum_D;           /* number of 2H isotopic atoms      */
    AT_NUMB   nNum_T;           /* number of 3H isotopic atoms      */
} INChI_IsotopicTGroup;




typedef struct tagINChI_Stereo
{
    /* [N] = allocated length */
    /* ---- possibly tetrahedral stereogenic atoms */

    int         nNumberOfStereoCenters;

    AT_NUMB    *nNumber;      /* Canonical number of a possibly tetrahedral
                               * stereogenic atom or allenes [nNumberOfStereoCenters]       */

    S_CHAR     *t_parity;     /* tetrahedral (relative, see nCompInv2Abs) atom parities [nNumberOfStereoCenters] */
    /* ---- possibly tetrahedral stereogenic atoms of the inverted structure         */

    AT_NUMB    *nNumberInv;   /* Canonical number of a possibly tetrahedral
                               * stereogenic atom or allene [nNumberOfStereoCenters]        */

    S_CHAR     *t_parityInv;  /* tetrahedral inverted atom parities [nNumberOfStereoCenters]*/

    /* bFlagAbsStereoIsInverted = nCompInv2Abs==-1: Abs stereo = Inverted           */
    int         nCompInv2Abs; /* 0=>Inv = Abs stereo; -1=> Inv < Abs stereo,
                              +1=> Inv > Abs stereo; +2=> in reading InChI:
                              no /m was found and in /sN N>0                        */

    int         bTrivialInv;  /* 1=> nCompInv2Abs!= 0 && Inverted = Abs stereo with inverted parities 1<-->2 */

    /* ---- possibly stereogenic bonds and tetrahedral cumulenes */
    int         nNumberOfStereoBonds;
    AT_NUMB    *nBondAtom1; /* Canonical number of a first atom
                             * [number of bonds]                                    */
    AT_NUMB    *nBondAtom2; /* Canonical number of a second atom
                             * [number of bonds]                                    */
    S_CHAR     *b_parity;   /* possibly stereogenic bond parities
                             * [number of bonds]                                    */
} INChI_Stereo;
#define INCHI_FLAG_ACID_TAUT            0x0001   /* tautomerism of dissociated acid invoked */
#define INCHI_FLAG_REL_STEREO           0x0002   /* requested relative stereo               */
#define INCHI_FLAG_RAC_STEREO           0x0004   /* requested racemic stereo                */
#define INCHI_FLAG_SC_IGN_ALL_UU        0x0008   /* ignored all undefined/unknown stereocenters, non-isotopic */
#define INCHI_FLAG_SB_IGN_ALL_UU        0x0010   /* ignored all undefined/unknown stereocenters, non-isotopic */
#define INCHI_FLAG_SC_IGN_ALL_ISO_UU    0x0020   /* ignored all undefined/unknown stereocenters, isotopic */
#define INCHI_FLAG_SB_IGN_ALL_ISO_UU    0x0040   /* ignored all undefined/unknown stereocenters, isotopic */
#define INCHI_FLAG_HARD_ADD_REM_PROTON  0x0080   /* in normalization a proton has been added or removed along alt path */

#define INCHI_OUT_NO_AUX_INFO           0x0001   /* do not output Aux Info                  */
#define INCHI_OUT_SHORT_AUX_INFO        0x0002   /* output short version of Aux Info        */
#define INCHI_OUT_ONLY_AUX_INFO         0x0004   /* output only Aux Info                    */
#define INCHI_OUT_EMBED_REC             0x0008   /* embed reconnected INChI into disconnected INChI */
#define INCHI_OUT_SDFILE_ONLY           0x0010   /* save input data in a Molfile instead of creating INChI */
#define INCHI_OUT_XML                   0x0020   /* *obsolete* output xml INChI             */
#define INCHI_OUT_PLAIN_TEXT            0x0040   /* output plain text INChI                 */
#define INCHI_OUT_PLAIN_TEXT_COMMENTS   0x0080   /* output plain text annotation            */
#define INCHI_OUT_XML_TEXT_COMMENTS     0x0100   /* *obsolete* output xml text annotation   */
#define INCHI_OUT_WINCHI_WINDOW         0x0200   /* output into wINChI text window          */
#define INCHI_OUT_TABBED_OUTPUT         0x0400   /* tab-delimited (only for plain text)     */
#define INCHI_OUT_SDFILE_ATOMS_DT       0x0800   /* SDfile output H isotopes as D and T     */
#define INCHI_OUT_SDFILE_SPLIT          0x1000   /* Split SDfile into components            */

#define INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG 0x2000
                                        /* used to accomodate FIX_TRANSPOSITION_CHARGE_BUG */
#define INCHI_OUT_STDINCHI 0x4000
#define INCHI_OUT_SAVEOPT  0x8000

/* bINChIOutputOptions2+ */
#define INCHI_OUT_INCHI_GEN_ERROR        0x0001 /* v. 1.05 */
#define INCHI_OUT_MISMATCH_AS_ERROR      0x0002 /* v. 1.05 */
#define INCHI_OUT_REQ_LARGE_MOL             0x0004 /* v. 1.06 */
#define INCHI_OUT_REQ_POLYMERS             0x0008 /* v. 1.06 */

/* Bits encoding InChI creation options to be saved */
#define SAVE_OPT_SLUUD  0x0001
#define SAVE_OPT_SUU    0x0002
#define SAVE_OPT_FIXEDH 0x0004
#define SAVE_OPT_RECMET 0x0008
#define SAVE_OPT_KET    0x0010
#define SAVE_OPT_15T    0x0020
#define SAVE_OPT_PT_22_00 0x0040
#define SAVE_OPT_PT_16_00 0x0080
#define SAVE_OPT_PT_06_00 0x0100
#define SAVE_OPT_PT_39_00 0x0200
#define SAVE_OPT_PT_13_00 0x0400
#define SAVE_OPT_PT_18_00 0x0800

#define INCHI_OUT_PRINT_OPTIONS       (INCHI_OUT_EMBED_REC |           \
                                       INCHI_OUT_PLAIN_TEXT |          \
                                       INCHI_OUT_PLAIN_TEXT_COMMENTS)


/*******REQ_MODE_SB_IGN_ALL_UU*************** chemical identifier definition *****************/
typedef struct tagINChI {  /* [N] = allocated length */

    int        nErrorCode;  /* 0 = success */
    INCHI_MODE  nFlags;      /* INCHI_FLAG_ACID_TAUT            tautomerism of dissociated acid invoked
                                INCHI_FLAG_REL_STEREO           requested relative stereo
                                INCHI_FLAG_RAC_STEREO           requested racemic stereo
                                INCHI_FLAG_SC_IGN_ALL_UU        ignored all undefined/unknown stereocenters, non-isotopic
                                INCHI_FLAG_SB_IGN_ALL_UU        ignored all undefined/unknown stereocenters, non-isotopic
                                INCHI_FLAG_SC_IGN_ALL_ISO_UU    ignored all undefined/unknown stereocenters, isotopic
                                INCHI_FLAG_SB_IGN_ALL_ISO_UU    ignored all undefined/unknown stereocenters, isotopic
                                INCHI_FLAG_HARD_ADD_REM_PROTON  in normalization a proton has been added or removed along alt path
                              */
    /* ---- basic & tautomer layer */
    int        nTotalCharge;
    int        nNumberOfAtoms;
    char      *szHillFormula;
    U_CHAR    *nAtom;       /* atomic numbers [nNumberOfAtoms] from the Periodic Table */
    int        lenConnTable;
    AT_NUMB   *nConnTable;  /* Connection table [nNumberOfAtoms+NumberOfBonds] */
    int        lenTautomer;
    AT_NUMB   *nTautomer;   /* NumGroups; ((NumAt+1, NumH,         At1..AtNumAt),...); {INCHI_T_NUM_MOVABLE = 1} - old
                             * NumGroups; ((NumAt+2, NumH, Num(-), At1..AtNumAt),...); {INCHI_T_NUM_MOVABLE = 2} - new
                             * Allocated length: [5*nNumberOfAtoms/2+1], see Alloc_INChI(...) */
    S_CHAR    *nNum_H;      /* number of terminal hydrogen atoms on each atom; in tautomeric
                             * representation these H on tautomeric atoms are not included [nNumberOfAtoms] */
    S_CHAR    *nNum_H_fixed;/* number of terminal hydrogen atoms on tautomeric atoms,
                             * in non-atautomeric representation only [nNumberOfAtoms] */
    /* ---- isotopic & isotopic tautomeric layer */
    int                  nNumberOfIsotopicAtoms;
    INChI_IsotopicAtom   *IsotopicAtom;              /* [nNumberOfIsotopicAtoms] */
    int                  nNumberOfIsotopicTGroups;
    /* in reversing InChI keeps a pointer to stolen from AuxInfo coordinates */
    INChI_IsotopicTGroup *IsotopicTGroup;             /* [nNumberOfIsotopicAtoms] */
    /* ---- stereo layer */
    INChI_Stereo *Stereo;
    INChI_Stereo *StereoIsotopic;
                                                /* not including mobile H groups */
    AT_NUMB *nPossibleLocationsOfIsotopicH;     /* [0]=> length including 0th element, location1,...*/
    int      bDeleted;
#if ( bREUSE_INCHI == 1 )
    int nRefCount;
#endif
#if ( bRELEASE_VERSION == 0 )
    int bExtract;
#endif
#if ( READ_INCHI_STRING == 1 )
    int nLink;  /* negative: ignore InChI; positive: index of (Reconnected component) + 1 linked to it */
#endif
} INChI;

typedef INChI *PINChI2[TAUT_NUM];

typedef struct tagOrigInfo {
    S_CHAR cCharge;
    S_CHAR cRadical;       /* 0=none, 1=doublet, 2=triplet, 3=unknown */
    S_CHAR cUnusualValence; /* see get_unusual_el_valence() */
} ORIG_INFO;
/******************** auxiliary chemical identifier info **************/
typedef struct tagINChI_Aux { /* [N] = allocated length */

    int			nErrorCode;  /* 0 = success */
    int			nNumberOfAtoms;
    int			nNumberOfTGroups;   /* non-zero only in tautomeric representation */
    int			bIsIsotopic;        /* filled out even though isotopic has not been requested */
    int			bIsTautomeric;      /* filled out even though tautomeric has not been requested; non-zero if taut exists */
    /* canonical numbers of the atoms: nOrigAtNosInCanonOrd[i-1]+1 =       */
    /*                       input atom number for the canonical number i  */
    AT_NUMB		*nOrigAtNosInCanonOrd;             /* [nNumberOfInputAtoms*1.5]; */
    AT_NUMB		*nIsotopicOrigAtNosInCanonOrd;     /* [nNumberOfInputAtoms*1.5]; */
    /* same for the inverted structure */
    AT_NUMB		*nOrigAtNosInCanonOrdInv;          /* inverted stereo [nNumberOfInputAtoms*1.5]; */
    AT_NUMB		*nIsotopicOrigAtNosInCanonOrdInv;  /* [nNumberOfInputAtoms*1.5]; */
    AT_NUMB		*nConstitEquNumbers;               /* [nNumberOfAtoms*1.5] */
    AT_NUMB		*nConstitEquTGroupNumbers;         /* [nNumberOfAtoms/2] */
    AT_NUMB		*nConstitEquIsotopicNumbers;       /* [nNumberOfAtoms*1.5] */
    AT_NUMB		*nConstitEquIsotopicTGroupNumbers; /* [nNumberOfAtoms/2] */
#if ( bREUSE_INCHI == 1 )
    int			nRefCount;
#endif

    ORIG_INFO	*OrigInfo;
    MOL_COORD	*szOrigCoord;
    NUM_H		nNumRemovedProtons;
    NUM_H		nNumRemovedIsotopicH[NUM_H_ISOTOPES]; /* isotopic H that may be exchanged and considered
                                                     randomly distributed, including removed protons;
                                                     order: 0=>1H, 1=>D, 2=>T */
    int         bDeleted;
    INCHI_MODE  bTautFlags;             /* t_group_info->bTautFlags */
    INCHI_MODE  bTautFlagsDone;         /* t_group_info->bTautFlagsDone */
    INCHI_MODE  bNormalizationFlags;    /* t_group_info->tni.bNormalizationFlags */
    int         nCanonFlags;
} INChI_Aux;

typedef INChI_Aux *PINChI_Aux2[TAUT_NUM];

/********************* array of pointers for sorting components and INChI output *********/
typedef struct tagINChIforSort {
    INChI     *pINChI[TAUT_NUM];
    INChI_Aux *pINChI_Aux[TAUT_NUM];
    short      ord_number;    /* for stable sort */
    short      n1; /* points to the original; used in structure reconstruction only */
    short      n2; /* points to the original; used in structure reconstruction only */
    short      n3; /* points to the original; used in structure reconstruction only */
}INCHI_SORT;



#endif    /* _ICHI_H_ */
