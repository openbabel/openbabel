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

#ifndef _STRUTIL_H_
#define _STRUTIL_H_

/* Mol structure utlities */

#include "inpdef.h"
#include "ichi.h"

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif

    /* forward declaration */
    struct tagTautomerGroupsInfo;
    struct tagCANON_GLOBALS;

    /**
     * @brief Extract one (connected) component
     *
     * @param at Pointer to atom array
     * @param num_at Number of atoms
     * @param component_number Component number
     * @param component_at Pointer to component atom array
     * @return int Number of extracted atoms
     */
    int ExtractConnectedComponent(inp_ATOM *at, int num_at,
                                  int component_number,
                                  inp_ATOM *component_at);
    /**
     * @brief Set the Connected Component Number object
     *
     * @param at Pointer to atom array
     * @param num_at Number of atoms
     * @param component_number Component number
     * @return int Returns 0
     */
    int SetConnectedComponentNumber(inp_ATOM *at, int num_at, int component_number);

    /**
     * @brief Allocate INChI data structure
     *
     * @param at Pointer to atom array
     * @param num_at Number of atoms
     * @param found_num_bonds Pointer to number of bonds
     * @param found_num_isotopic Pointer to number of isotopic atoms
     * @param nAllocMode Allocation mode
     * @return INChI* Pointer to allocated INChI
     */
    INChI *Alloc_INChI(inp_ATOM *at, int num_at, int *found_num_bonds,
                       int *found_num_isotopic, int nAllocMode);

    /**
     * @brief Free INChI data structure
     *
     * @param ppINChI Pointer to INChI data structure
     * @return int Returns 0
     */
    int Free_INChI(INChI **ppINChI);

    /**
     * @brief Free INChI members
     *
     * @param pINChI Pointer to INChI data structure
     * @return int Returns 0
     */
    int Free_INChI_Members(INChI *pINChI);

    /**
     * @brief Free INChI_Stereo data structure
     *
     * @param pINChI_Stereo Pointer to INChI_Stereo
     * @return int Returns 0
     */
    int Free_INChI_Stereo(INChI_Stereo *pINChI_Stereo);

    /**
     * @brief Allocate AuxInfo data structure
     *
     * @param num_at Number of atoms
     * @param num_isotopic_atoms Number of isotopic atoms
     * @param nAllocMode Allocation mode
     * @param bOrigData Flag indicating if the atom data is original
     * @return INChI_Aux* Pointer to allocated INChI_Aux or NULL if failed
     */
    INChI_Aux *Alloc_INChI_Aux(int num_at, int num_isotopic_atoms,
                               int nAllocMode, int bOrigData);

    /**
     * @brief Free INChI_Aux data structure
     *
     * @param ppINChI_Aux Pointer to INChI_Aux data structure
     * @return int Returns 0
     */
    int Free_INChI_Aux(INChI_Aux **ppINChI_Aux);

    /**
     * @brief Create INChI
     *
     * @param pCG Pointer to canonicalization globals data structure
     * @param ic Pointer to InChI clock data structure
     * @param ip Pointer to input parameters
     * @param ppINChI Pointer to INChI data structure
     * @param ppINChI_Aux Pointer to AuxInfo data structure
     * @param orig_inp_data Pointer to original atom data
     * @param inp_at Pointer to atom array
     * @param inp_norm_data Pointer to normalized atom data
     * @param num_inp_at Number of atoms
     * @param nUserMode User mode
     * @param pbTautFlags Tautomer flags
     * @param pbTautFlagsDone Tautomer flags completed
     * @param ulMaxTime Maximum time
     * @param ti_out Tautomer info
     * @param pStrErrStruct Store error or warning messages
     * @return int Status code
     */
    int Create_INChI(struct tagCANON_GLOBALS *pCG,
                     struct tagINCHI_CLOCK *ic,
                     INPUT_PARMS *ip,
                     INChI **ppINChI,
                     INChI_Aux **ppINChI_Aux,
                     ORIG_ATOM_DATA *orig_inp_data,
                     inp_ATOM *inp_at,
                     INP_ATOM_DATA *inp_norm_data[2],
                     int num_inp_at,
                     INCHI_MODE nUserMode,
                     INCHI_MODE *pbTautFlags,
                     INCHI_MODE *pbTautFlagsDone,
                     struct tagInchiTime *ulMaxTime,
                     struct tagTautomerGroupsInfo *ti_out,
                     char *pStrErrStruct);

    /**
     * @brief Fill out
     *
     * @param pCG Pointer to InChI clock data structure
     * @param norm_at Pointer to normalized atoms
     * @param inf_norm_at_data Pointer to normalized info atom data structure
     * @param init_num_at Number of atoms
     * @param num_removed_H Number of removed hydrogens
     * @param bAdd_DT_to_num_H Add isotopic hydrogens to the number of hydrogens
     * @param nNumRemovedProtons Number of removed protons
     * @param nNumRemovedProtonsIsotopic Pointer to array of removed protons (isotopic)
     * @param bIsotopic Flag indicating if the atom data is isotopic
     * @param pINChI Pointer to INChI data structure
     * @param pINChI_Aux Pointer to AuxInfo data structure
     * @param bAbcNumbers Flag indicating whether to use alphanumeric values (?)
     * @param nMode Stereochemistry mode
     * @return int Returns 0
     */
    int FillOutInfAtom(struct tagCANON_GLOBALS *pCG,
                       inp_ATOM *norm_at,
                       INF_ATOM_DATA *inf_norm_at_data,
                       int init_num_at,
                       int num_removed_H,
                       int bAdd_DT_to_num_H,
                       int nNumRemovedProtons,
                       NUM_H *nNumRemovedProtonsIsotopic,
                       int bIsotopic,
                       INChI *pINChI,
                       INChI_Aux *pINChI_Aux,
                       int bAbcNumbers,
                       INCHI_MODE nMode);

    /**
     * @brief Fill out composite canonical info atom data structure
     *
     * @param pCG Pointer to canonicalization globals data structure
     * @param composite_norm_data Pointer to composite normalized atom data
     * @param inf_norm_at_data Pointer to normalized info atom data structure
     * @param bIsotopic Flag indicating if the atom data is isotopic
     * @param bTautomeric Flag indicating if the atom data is tautomeric
     * @param pINChI2 Pointer to INChI2 data structure (?)
     * @param pINChI_Aux2 Pointer to AuxInfo2 data structure (?)
     * @param bAbcNumbers Flag indicating whether to use alphanumeric values (?)
     * @param nMode Stereochemistry mode
     * @return int Returns 1
     */
    int FillOutCompositeCanonInfAtom(struct tagCANON_GLOBALS *pCG,
                                     COMP_ATOM_DATA *composite_norm_data,
                                     INF_ATOM_DATA *inf_norm_at_data,
                                     int bIsotopic,
                                     int bTautomeric,
                                     PINChI2 *pINChI2,
                                     PINChI_Aux2 *pINChI_Aux2,
                                     int bAbcNumbers,
                                     INCHI_MODE nMode);

#if (FIX_DALKE_BUGS == 1)
    /**
     * @brief Allocate and fill hill formula
     *
     * @param pINChI Pointer to INChI data structure
     * @return char* Pointer to allocated hill formula
     */
    char *AllocateAndFillHillFormula(INChI *pINChI);
#endif

    typedef enum tagInchiDiffBits
    {
        IDIF_PROBLEM = 0x00000001,           /* severe: at least one InChI does not exist */
        IDIF_NUM_AT = 0x00000001,            /* severe: different number of atoms in the skeleton */
        IDIF_ATOMS = 0x00000001,             /* severe: diiferent types of skeleton atoms */
        IDIF_NUM_EL = 0x00000001,            /* severe: formulas differ in another element */
        IDIF_CON_LEN = 0x00000001,           /* severe: different connection table lengths */
        IDIF_CON_TBL = 0x00000001,           /* severe: different connection tables */
        IDIF_POSITION_H = 0x00000002,        /* difference in non-taut (Mobile-H) or all H (Fixed-H) location/number */
        IDIF_MORE_FH = 0x00000004,           /* extra fixed H */
        IDIF_LESS_FH = 0x00000008,           /* missing fixed H */
        IDIF_MORE_H = 0x00000010,            /* formulas differ in number of H */
        IDIF_LESS_H = 0x00000020,            /* formulas differ in number of H */
        /*IDIF_TAUT_LEN      = 0x00000008,*/ /* different lengths of tautomer lists */
        IDIF_NO_TAUT = 0x00000040,           /* restored structure has no taut groups while the original InChI has some */
        IDIF_WRONG_TAUT = 0x00000080,        /* restored has tautomerism while the original does not have it */
        IDIF_SINGLE_TG = 0x00000100,         /* restored has 1 taut. group while the original InChI has multiple tg */
        IDIF_MULTIPLE_TG = 0x00000200,       /* restored has multiple tg while the original InChI has only one tg */
        IDIF_NUM_TG = 0x00000400,            /* different number of tautomeric groups */
        /*IDIF_LESS_TG_ENDP  = 0x00000200,*/ /* restores structure has less taut. endpoints */
        /*IDIF_MORE_TG_ENDP  = 0x00000400,*/ /* restores structure has more taut. endpoints */
        IDIF_EXTRA_TG_ENDP = 0x00000800,     /* extra tautomeric endpoint(s) in restored structure */
        IDIF_MISS_TG_ENDP = 0x00001000,      /* one or more tg endpoint is not in the restored structure */
        IDIF_DIFF_TG_ENDP = 0x00002000,      /* lists of tg endpoints are different */
        IDIF_TG = 0x00004000,                /* different tautomeric groups */
        IDIF_NUM_ISO_AT = 0x00008000,        /* ?severe: restored struct. has different number of isotopic atoms */
        IDIF_ISO_AT = 0x00010000,            /* ?severe: restored struct. has different locations/isotopes of isotopic atoms */
        IDIF_CHARGE = 0x00020000,            /* restored structure has different charge */
        IDIF_REM_PROT = 0x00040000,          /* proton(s) removed/added from the restored structure */
        IDIF_REM_ISO_H = 0x00080000,         /* isotopic H removed */
        IDIF_SC_INV = 0x00100000,            /* restores structure has different inversion stereocenter mark */
        IDIF_SC_PARITY = 0x00200000,         /* restored structure has stereoatoms or allenes with different parity */
        IDIF_SC_EXTRA_UNDF = 0x00400000,     /* restored structure has extra undefined stereocenter(s) */
        IDIF_SC_EXTRA = 0x00800000,          /* restored structure has extra stereocenter(s) */
        IDIF_SC_MISS_UNDF = 0x01000000,      /* restored structure has not some undefined stereocenter(s) */
        IDIF_SC_MISS = 0x02000000,           /* restored structure has not some stereocenters that are not undefined */
        IDIF_SB_PARITY = 0x04000000,         /* restored structure has stereobonds or cumulenes with different parity */
        IDIF_SB_EXTRA_UNDF = 0x08000000,     /* restored structure has extra undefined stereobond(s) */
        IDIF_SB_EXTRA = 0x10000000,          /* restored structure has extra stereobond(s) */
        IDIF_SB_MISS_UNDF = 0x20000000,      /* restored structure has not some undefined stereocenters */
        IDIF_SB_MISS = 0x40000000            /* restored structure has not some stereobonds that are not undefined */
    } IDIF;

#define IDIFF_SB (IDIF_SB_PARITY | IDIF_SB_EXTRA_UNDF | IDIF_SB_EXTRA | IDIF_SB_MISS_UNDF | IDIF_SB_MISS)
#define IDIFF_SC (IDIF_SC_PARITY | IDIF_SC_EXTRA_UNDF | IDIF_SC_EXTRA | IDIF_SC_MISS_UNDF | IDIF_SC_MISS)

#define IDIFF_CONSTIT (IDIF_POSITION_H | IDIF_MORE_FH | IDIF_LESS_FH | IDIF_MORE_H | IDIF_LESS_H |    \
                       IDIF_NO_TAUT | IDIF_WRONG_TAUT | IDIF_SINGLE_TG | IDIF_MULTIPLE_TG |           \
                       IDIF_NUM_TG | IDIF_EXTRA_TG_ENDP | IDIF_MISS_TG_ENDP | IDIF_TG |               \
                       IDIF_NUM_ISO_AT | IDIF_ISO_AT | IDIF_CHARGE | IDIF_REM_PROT | IDIF_REM_ISO_H | \
                       IDIF_DIFF_TG_ENDP)
#define IDIFF_STEREO (IDIF_SC_INV | IDIF_SC_PARITY | IDIF_SC_EXTRA_UNDF | IDIF_SC_EXTRA |      \
                      IDIF_SC_MISS_UNDF | IDIF_SC_MISS | IDIF_SB_PARITY | IDIF_SB_EXTRA_UNDF | \
                      IDIF_SB_EXTRA | IDIF_SB_MISS_UNDF | IDIF_SB_MISS)

/*************************************************************************************/
#define ICR_MAX_ENDP_IN1_ONLY 32
#define ICR_MAX_ENDP_IN2_ONLY 32
#define ICR_MAX_DIFF_FIXED_H 32
#define ICR_MAX_SB_IN1_ONLY 32
#define ICR_MAX_SB_IN2_ONLY 32
#define ICR_MAX_SC_IN1_ONLY 32
#define ICR_MAX_SC_IN2_ONLY 32
#define ICR_MAX_SB_UNDF 32
#define ICR_MAX_SC_UNDF 32

    typedef struct tagInChICompareResults
    {
        INCHI_MODE flags;

        int tot_num_H1;
        int tot_num_H2;
        int num_taut_H1;
        int num_taut_H2;
        int num_taut_M1;
        int num_taut_M2;

        /* 1 => InChI from reversed struct. 2 => input InChI */

        AT_NUMB endp_in1_only[ICR_MAX_ENDP_IN1_ONLY]; /* endpoint canonical number = index+1 */
        int num_endp_in1_only;

        AT_NUMB endp_in2_only[ICR_MAX_ENDP_IN2_ONLY]; /* endpoint canonical number = index+1 */
        int num_endp_in2_only;

        AT_NUMB diff_pos_H_at[ICR_MAX_DIFF_FIXED_H]; /* non-tautomeric H */
        S_CHAR diff_pos_H_nH[ICR_MAX_DIFF_FIXED_H];
        int num_diff_pos_H;

        AT_NUMB fixed_H_at1_more[ICR_MAX_DIFF_FIXED_H]; /* extra fixed_H */
        S_CHAR fixed_H_nH1_more[ICR_MAX_DIFF_FIXED_H];
        int num_fixed_H1_more;

        AT_NUMB fixed_H_at2_more[ICR_MAX_DIFF_FIXED_H]; /* missed fixed_H */
        S_CHAR fixed_H_nH2_more[ICR_MAX_DIFF_FIXED_H];
        int num_fixed_H2_more;

        AT_NUMB sc_in1_only[ICR_MAX_SC_IN1_ONLY];
        int num_sc_in1_only;
        AT_NUMB sc_in2_only[ICR_MAX_SC_IN2_ONLY];
        int num_sc_in2_only;

        AT_NUMB sb_in1_only[ICR_MAX_SB_IN1_ONLY];
        int num_sb_in1_only;
        AT_NUMB sb_in2_only[ICR_MAX_SB_IN2_ONLY];
        int num_sb_in2_only;

        AT_NUMB sb_undef_in1_only[ICR_MAX_SC_UNDF];
        int num_sb_undef_in1_only;
        AT_NUMB sb_undef_in2_only[ICR_MAX_SC_UNDF];
        int num_sb_undef_in2_only;

        AT_NUMB sc_undef_in1_only[ICR_MAX_SB_UNDF];
        int num_sc_undef_in1_only;
        AT_NUMB sc_undef_in2_only[ICR_MAX_SB_UNDF];
        int num_sc_undef_in2_only;
    } ICR; /* tagInChICompareResults */

    /**
     * @brief Compares two InChIs (2)
     *
     * @param i1 Pointer to first InChI (InChI from reversed struct)
     * @param i2 Pointer to second InChI
     * @param a1 Pointer to first InChI AuxInfo
     * @param a2 Pointer to second InChI AuxInfo
     * @param picr Pointer to comparison results
     * @param err Pointer to error code
     * @return INCHI_MODE Returns code from IDIF
     */
    INCHI_MODE CompareReversedINChI2(INChI *i1 /* InChI from reversed struct */,
                                     INChI *i2 /* input InChI */,
                                     INChI_Aux *a1, INChI_Aux *a2,
                                     ICR *picr, int *err);

    /**
     * @brief Compares results from InChI comparison
     *
     * @param picr1 Pointer to first comparison results
     * @param picr2 Pointer to second comparison results
     * @param pin1 Pointer to first InChI mode (never used?)
     * @param pin2 Pointer to second InChI mode (never used?)
     * @param mask Either IDIFF_CONSTIT or IDIFF_STEREO
     * @return * int 2 if undefined results, 1 if unequal, 0 if equal, -1 if unequal
     */
    int CompareIcr(ICR *picr1, ICR *picr2, INCHI_MODE *pin1, INCHI_MODE *pin2, INCHI_MODE mask);

    /**
     * @brief Compares two InChIs (1)
     *
     * @param i1 Pointer to first InChI (InChI from reversed struct)
     * @param i2 Pointer to second InChI
     * @param a1 Pointer to first InChI AuxInfo
     * @param a2 Pointer to second InChI AuxInfo
     * @return int 0 if NULL, 1 if one is NULL, >1 if difference in structure
     */
    int CompareReversedINChI(INChI *i1, INChI *i2, INChI_Aux *a1, INChI_Aux *a2);

    // Where is this function?
    const char *CompareReversedInchiMsg(int code);

#define EQL_EXISTS 1
#define EQL_SP3 2
#define EQL_SP3_INV 4
#define EQL_SP2 8

    /**
     * @brief Compares stereo information of two structures
     *
     * @param s1 Pointer to stereo information of first structure
     * @param eql1 Flag for stereo information check
     * @param s2 Pointer to stereo information of second structure
     * @param eql2 Flag for stereo information check
     * @param bRelRac Flag to compare racemic stereo information
     * @return int 0 if unequal, 1 if equal
     */
    int Eql_INChI_Stereo(INChI_Stereo *s1, int eql1, INChI_Stereo *s2, int eql2, int bRelRac);

    /**
     * @brief Compares isotopic information of two InChIs
     *
     * @param i1 Pointer to first InChI
     * @param i2 Pointer to second InChI
     * @return int 0 if unequal, 1 if equal
     */
    int Eql_INChI_Isotopic(INChI *i1, INChI *i2);

#define EQL_EQU 0
#define EQL_EQU_TG 1
#define EQL_EQU_ISO 2

    /**
     * @brief Compares two InChI AuxInfo objects
     *
     * @param a1 Pointer to first AuxInfo
     * @param eql1 Flag to compare type 1
     * @param a2 Pointer to second AuxInfo
     * @param eql2 Flag to compare type 2
     * @return int 0 if unequal or empty, 1 if equal
     */
    int Eql_INChI_Aux_Equ(INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2);

#define EQL_NUM 0
#define EQL_NUM_INV 1
#define EQL_NUM_ISO 2

    /**
     * @brief Compares two InChI AuxInfo objects in terms of numbering
     *
     * @param a1 Pointer to first AuxInfo
     * @param eql1 Flag to compare type 1
     * @param a2 Pointer to second AuxInfo
     * @param eql2 Flag to compare type 2
     * @return int 0 if unequal or empty, 1 if numbering matches
     */
    int Eql_INChI_Aux_Num(INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2);

    /**
     * @brief Checks if a given array of equivalence numbers (LinearCT) contains any repetitions
     *
     * @param LinearCT Pointer to equivalence numbers
     * @param nLenCT Length of the array
     * @return int 1 if the array contains repetitions, 0 otherwise
     */
    int bHasEquString(AT_NUMB *LinearCT, int nLenCT);

    /**
     * @brief Compare function for qsort (normal structure)
     *
     * @param p1 Pointer to first INCHI_SORT structure
     * @param p2 Pointer to second INCHI_SORT structure
     * @return int -1 if *p1 > *p2, return +1 if *p1 < *p2
     */
    int CompINChINonTaut2(const void *p1, const void *p2);

    /**
     * @brief Compare function for qsort (tautomeric structures)
     *
     * @param p1 Pointer to first INCHI_SORT structure
     * @param p2 Pointer to second INCHI_SORT structure
     * @return int -1 if *p1 > *p2, return +1 if *p1 < *p2
     */
    int CompINChITaut2(const void *p1, const void *p2);

    /**
     * @brief Compares two INChI data structures (qsort)
     *
     * @param p1 Pointer to the first INCHI_SORT data structure
     * @param p2 Pointer to the second INCHI_SORT data structure
     * @param bTaut Flag for tautomeric comparison
     * @param bCompareIsotopic Flag for isotopic comparison
     * @return int -1 if *p1 > *p2, return +1 if *p1 < *p2
     */
    int CompINChI2(const INCHI_SORT *p1, const INCHI_SORT *p2, int bTaut, int bCompareIsotopic);

    /**
     * @brief Compare tautomeric vs non-tautomeric information
     *
     * @param p1 Pointer to the first INCHI_SORT data structure
     * @param p2 Pointer to the second INCHI_SORT data structure
     * @param bCompareIsotopic Flag for comparing isotopic non-tautomeric part
     * @return int -1 if *p1 > *p2, return +1 if *p1 < *p2
     */
    int CompINChITautVsNonTaut(const INCHI_SORT *p1, const INCHI_SORT *p2, int bCompareIsotopic);

    typedef enum tagDiffINChISegments
    {                  /* r = repetitive, n = non-repetitive */
      DIFS_f_FORMULA,  /*  0 r; fixed-H <-> mobile-H */
      DIFS_c_CONNECT,  /*  1 n; connection table; mobile-H only */
      DIFS_h_H_ATOMS,  /*  2 n; hydrogen atoms: mobile-H and Fixed-H; have different meanings */
      DIFS_q_CHARGE,   /*  3 r; charge; fixed-H <-> mobile-H */
      DIFS_p_PROTONS,  /*  4 n; protons; mobile-H only */
      DIFS_b_SBONDS,   /*  5 r: stereobonds: fixed-H <-> mobile-H * isotopic <-> non-isotopic */
      DIFS_t_SATOMS,   /*  6 r: stereoatoms: fixed-H <-> mobile-H * isotopic <-> non-isotopic */
      DIFS_m_SP3INV,   /*  7 r: stereo-abs-inv: fixed-H <-> mobile-H * isotopic <-> non-isotopic */
      DIFS_s_STYPE,    /*  8 r: stereo-type: fixed-H <-> mobile-H * isotopic <-> non-isotopic */
      DIFS_i_IATOMS,   /*  9 r: isotopic atoms: fixed-H <-> mobile-H * isotopic <-> non-isotopic */
      DIFS_o_TRANSP,   /* 10 n: Fixed-H transposition */
      DIFS_idf_LENGTH, /* 11    length of the array relevant to the INChI Identifier */
      /* later elements referring to AuxInfo may be added */
      DIFS_LENGTH = DIFS_idf_LENGTH /* length of the array */
    } DIF_SEGMENTS;

    typedef enum tagDiffINChILayers
    {
        DIFL_M,     /* 0 main layer */
        DIFL_MI,    /* 1 main isotopic */
        DIFL_F,     /* 2 fixed-H */
        DIFL_FI,    /* 3 fixed-H isotopic */
        DIFL_LENGTH /* number of layers */
    } DIF_LAYERS;

    /* Value meaning */
    typedef enum tagMarkDiff
    {
        DIFV_BOTH_EMPTY = 0, /* both this and the component in the preceding namesake segment are empty  */
        DIFV_EQL2PRECED = 1, /* equal to the component in the preceding namesake segment                 */
        DIFV_NEQ2PRECED = 2, /* different from the component in the preceding namesake segment           */
        DIFV_IS_EMPTY = 4,   /* is empty while the preceding namesake segment is not empty               */
        DIFV_FI_EQ_MI = 8,   /* FI stereo component is equal to the component in the MI namesake segment */
                             /* while M and F components are empty                                       */
        /* decision_F = bitmask: bits that should not be present */
        /* decision_T = bitmask: at least one of the bits should be present */
        /* decision = true if( !( BITS & decision_F ) && ( BITS & decision_F ) ) */
        DIFV_OUTPUT_EMPTY_T = (DIFV_IS_EMPTY),                                     /* bits present for empty segment output */
        DIFV_OUTPUT_EMPTY_F = (DIFV_EQL2PRECED | DIFV_NEQ2PRECED | DIFV_FI_EQ_MI), /* bits NOT present */

        DIFV_OUTPUT_OMIT_F = (DIFV_NEQ2PRECED | DIFV_IS_EMPTY), /* bits NOT present for omitting */

        DIFV_OUTPUT_FILL_T = (DIFV_EQL2PRECED | DIFV_NEQ2PRECED | DIFV_FI_EQ_MI)
    } DIF_VALUES;

    typedef enum tagINChISegmAction
    {
        INCHI_SEGM_OMIT = 0,
        INCHI_SEGM_FILL = 1, /* the value is used in str_LineEnd() */
        INCHI_SEGM_EMPTY = 2 /* the value is used in str_LineEnd() */
    } INCHI_SEGM_ACTION;

    /**
     * @brief Compare InChI layers
     *
     * @param p1 Pointer to InChI component 1
     * @param p2 Pointer to InChI component 2
     * @param sDifSegs Array of different layers
     * @param bFixTranspChargeBug Flag to fix charge bug (?)
     * @return int -1 if *p1 > *p2, return +1 if *p1 < *p2
     */
    int CompINChILayers(const INCHI_SORT *p1, const INCHI_SORT *p2,
                        char sDifSegs[][DIFS_LENGTH], int bFixTranspChargeBug);

    /**
     * @brief Mark unused and empty layers
     *
     * @param sDifSegs Array of different layers
     * @return int returns 0
     */
    int MarkUnusedAndEmptyLayers(char sDifSegs[][DIFS_LENGTH]);

    /**
     * @brief Action to take per segment (?)
     *
     * @param cDifSegs Segment bits
     * @return int Type of segment action (INCHI_SEGM_OMIT, INCHI_SEGM_EMPTY, INCHI_SEGM_FILL)
     */
    int INChI_SegmentAction(char cDifSegs);

#define FLAG_SORT_PRINT_TRANSPOS_BAS 1   /* transposition in the main InChI layer */
#define FLAG_SORT_PRINT_TRANSPOS_REC 2   /* transposition in the reconnected InChI layer */
#define FLAG_SORT_PRINT_NO_NFIX_H_BAS 4  /* no fixed H non-isotopic in the main InChI layer */
#define FLAG_SORT_PRINT_NO_NFIX_H_REC 8  /* no fixed H non-isotopic in the reconnected InChI layer */
#define FLAG_SORT_PRINT_NO_IFIX_H_BAS 16 /* no fixed H isotopic in the main InChI layer */
#define FLAG_SORT_PRINT_NO_IFIX_H_REC 32 /* no fixed H isotopic in the the reconnected InChI layer */
#define FLAG_SORT_PRINT_ReChI_PREFIX 64  /* Output ReChI instead of InChI */

    struct tagCANON_GLOBALS;

    /**
     * @brief Main actual worker which serializes InChI to string (called from OutputINChI2( ... ) and from itself)
     *
     * @param pCG Pointer to canonicalisation parameters
     * @param strbuf Point to string buffer
     * @param pINChISortTautAndNonTaut2
     * @param iINChI InChI basic or InChI reconnect information (?)
     * @param orig_inp_data Pointer to original atom data structure
     * @param pOrigStruct Pointer to original stucture
     * @param ip Pointer to input parameters
     * @param bDisconnectedCoord Flag for disconnected coordinates
     * @param bOutputType Flag for output type
     * @param bINChIOutputOptions Flag for InChI output options
     * @param num_components2 Pointer to number of components
     * @param num_non_taut2 Pointer to number of non tautomeric units (?)
     * @param num_taut2 Pointer to number of tautomeric units
     * @param out_file Pointer to output file
     * @param log_file Pointer to log file
     * @param num_input_struct Number of input structures
     * @param pSortPrintINChIFlags Pointer to sort and print InChI flags (?)
     * @param save_opt_bits Flag to save optional bits (?)
     * @return int 0 if error, 1 otherwise
     */
    int OutputINChI1(struct tagCANON_GLOBALS *pCG,
                     INCHI_IOS_STRING *strbuf,
                     INCHI_SORT *pINChISortTautAndNonTaut2[][TAUT_NUM],
                     int iINChI,
                     ORIG_ATOM_DATA *orig_inp_data,
                     ORIG_STRUCT *pOrigStruct,
                     INPUT_PARMS *ip,
                     int bDisconnectedCoord,
                     int bOutputType,
                     int bINChIOutputOptions,
                     int num_components2[],
                     int num_non_taut2[],
                     int num_taut2[],
                     INCHI_IOSTREAM *out_file,
                     INCHI_IOSTREAM *log_file,
                     int num_input_struct,
                     int *pSortPrintINChIFlags,
                     unsigned char save_opt_bits);

    int OutputINChI2(struct tagCANON_GLOBALS *pCG,
                     INCHI_IOS_STRING *strbuf,
                     INCHI_SORT *pINChISortTautAndNonTaut2[][TAUT_NUM],
                     int INCHI_basic_or_INCHI_reconnected,
                     ORIG_ATOM_DATA *orig_inp_data,
                     ORIG_STRUCT *pOrigStruct,
                     INPUT_PARMS *ip,
                     int bDisconnectedCoord,
                     int bOutputType,
                     int bINChIOutputOptions,
                     int num_components2[],
                     int num_non_taut2[],
                     int num_taut2[],
                     INCHI_IOSTREAM *out_file,
                     INCHI_IOSTREAM *log_file,
                     int num_input_struct,
                     int *pSortPrintINChIFlags,
                     unsigned char save_opt_bits);

    /**
     * @brief Save equivalent components information and sort order (not used ?)
     *
     * @param iINChI Index to component
     * @param pINChISort Pointer to InChI sort structure
     * @param num_components Pointer to number of components
     * @param orig_inp_data Pointer to original atom data
     * @param prep_inp_data Pointer to prepared original atom data (?)
     * @param composite_norm_data Pointer to composite data (?)
     * @param bCompareComponents Type of comparisson: 1 default, 2 non-isotopic, 4 -non-tautomeric
     * @return * int 0 if error, 1 if otherwise
     */
    int SaveEquComponentsInfoAndSortOrder(int iINChI,
                                          INCHI_SORT *pINChISort[TAUT_NUM],
                                          int *num_components,
                                          ORIG_ATOM_DATA *orig_inp_data,
                                          ORIG_ATOM_DATA *prep_inp_data,
                                          COMP_ATOM_DATA composite_norm_data[TAUT_NUM + 1],
                                          int bCompareComponents); /* djb-rwth: matching composite_norm_data bounds */

    /**
     * @brief Print error message (plain text)
     *
     * @param out_file Pointer to output file
     * @param pErrorText Pointer to error text
     * @param bError Error type
     * @return int Returns 1
     */
    int OutputINChIPlainError(INCHI_IOSTREAM *out_file,
                              char *pErrorText,
                              int bError);

    /**
     * @brief Get the input struct error type
     *
     * @param ip Pointer to input parameters
     * @param err Error code
     * @param pStrErrStruct Pointer to error string
     * @param num_inp_atoms Number of input atoms
     * @return int Returns type of error
     */
    int GetInpStructErrorType(INPUT_PARMS *ip, int err, char *pStrErrStruct, int num_inp_atoms);

    /**
     * @brief
     *
     * @param out_file Pointer to output file
     * @param log_file Pointer to log file
     * @param pStrErrStruct Pointer to error string
     * @param nErrorType Error type
     * @param num_inp Number of input structure (?)
     * @param ip Pointer to input parameters
     * @return int Retunrs error type
     */
    int ProcessStructError(INCHI_IOSTREAM *out_file,
                           INCHI_IOSTREAM *log_file,
                           char *pStrErrStruct,
                           int nErrorType,
                           long num_inp,
                           INPUT_PARMS *ip);

    /**
     * @brief Check if hetero atoms has isotopic hydrogens
     *
     * @param atom Pointer to atom array
     * @param num_atoms Number of atoms
     * @return int 2 if isotopic atoms, 1 if isotopic hydrogens, 0 if not
     */
    int bNumHeterAtomHasIsotopicH(inp_ATOM *atom, int num_atoms);

    // does this function exist?
    int WriteToSDfile(const INP_ATOM_DATA *inp_at_data, INCHI_IOSTREAM *fcb, const char *name, const char *comment,
                      const char *szLabel, const char *szValue);

    /**
     * @brief Check if atom is a metal salt
     *
     * @param at Pointer to atom array
     * @param i atom index
     * @return int 1 if metal salt, 0 if not
     */
    int bIsMetalSalt(inp_ATOM *at, int i);

    /**
     * @brief Check if two bonds are the same
     *
     * @param a1 Atom a1 index
     * @param a2 Atom a2 index
     * @param b1 Atom b1 index
     * @param b2 Atom b2 index
     * @return int Returns 1 if bonds (a1,a2) and (b1,b2) are the same, -1 if atoms swapped, 0 if not the same
     */
    int bIsSameBond(int a1, int a2, int b1, int b2);

    extern const char gsMissing[];
    extern const char gsEmpty[];
    extern const char gsSpace[];
    extern const char gsEqual[];

#define SDF_LBL_VAL(L, V) ((L) && (L)[0]) ? gsSpace : gsEmpty, ((L) && (L)[0]) ? L : gsEmpty, ((L) && (L)[0]) ? (((V) && (V)[0]) ? gsEqual : gsSpace) : gsEmpty, ((V) && (V)[0]) ? V : ((L) && (L)[0]) ? gsMissing \
                                                                                                                                                                                                       : gsEmpty

    /* Handle integer matrix [mxn] */
    /**
     * @brief Allocate integer matrix
     *
     * @param m Size of matrix
     * @param n Size of matrix
     * @param a Pointer to matrix
     * @return int 1 if success, 0 if failure
     */
    int imat_new(int m, int n, int ***a);

    /**
     * @brief Free integer matrix
     *
     * @param m Size of matrix
     * @param a Pointer to matrix
     */
    void imat_free(int m, int **a);

    /* Light-weight {nodes, connections} subgraph of molecule, subgraf
       refers to parent structure via orig atom numbers stored internally
            nodes and edges: are living in 0-based node numbers domain
            atoms and bonds: are living in 1-based orig atom numbers domain
    */

    /* subgraf_edge is node connections element */
    typedef struct subgraf_edge
    {
        int nbr;   /* node neighbor									*/
        int etype; /* edge type echoing corresponding bond type		*/
    } subgraf_edge;
    /* subgraf is set of nodes and their connections */
    typedef struct subgraf
    {
        int nnodes;         /* number of nodes (atoms)							*/
        int *nodes;         /* nodes[i]  is 1-based orig atom number of node #i	*/
                            /* (i = 0 to nnodes)								*/
        int *degrees;       /* degrees of nodes									*/
        int *orig2node;     /* orig2node[k] is 0-based node number for orig #k	*/
        subgraf_edge **adj; /* node connections representing adjacency relation	*/
                            /* adj[i] is vector[0-degree] of edges at node #i	*/
    } subgraf;
    /* helper structure for finding paths in subgraf*/
    typedef struct subgraf_pathfinder
    {
        subgraf *sg;  /* base subgraf										*/
        int start;    /* start node of path, 0-based						*/
        int end;      /* end node of path									*/
        int maxbonds; /* max number of bonds in path						*/
        int nbonds;   /* actual number of bonds in path					*/
        int nseen;    /* number of nodes of found path(s)					*/
        int *seen;    /* nodes of found path(s)							*/
    } subgraf_pathfinder;

    /**
     * @brief Create graph from atom data
     *
     * @param orig_inp_data Pointer to original atom data
     * @param nnodes Number of nodes
     * @param nodes Pointer to nodes
     * @return subgraf* NULL if error, a new graph otherwise
     */
    subgraf *subgraf_new(ORIG_ATOM_DATA *orig_inp_data,
                         int nnodes,
                         int *nodes);

    /**
     * @brief Free graph data structure
     *
     * @param sg Pointer to graph
     * @return void Returns NULL (?)
     */
    void subgraf_free(subgraf *sg);
    /**
     * @brief Debug graph data structure
     *
     * @param sg Pointer to graph
     */
    void subgraf_debug_trace(subgraf *sg);

    /**
     * @brief Allocate new graph pathfinder data structure
     *
     * @param sg Pointer to graph
     * @param orig_inp_data Pointer to original atom data
     * @param start Starting node
     * @param end End node
     * @return subgraf_pathfinder* Returns pointer to graph pathfinder data structure
     */
    subgraf_pathfinder *subgraf_pathfinder_new(subgraf *sg,
                                               ORIG_ATOM_DATA *orig_inp_data,
                                               int start,
                                               int end);

    /**
     * @brief Frees subgraph pathfinder data structure
     *
     * @param spf Pointer to graph pathfinder data structure
     * @return * void
     */
    void subgraf_pathfinder_free(subgraf_pathfinder *spf);

    /**
     * @brief Find path(s) from subgraf node spf->start to spf->end and fill bonds[nbonds] and atoms[natoms]. Does not traverse through supplied forbidden edges (if not zero/NULL)
     *
     * @param spf Pointer to graph pathfinder data structure
     * @param nforbidden Number of forbidden edges
     * @param forbidden_orig Pointer forbidden edges
     * @param nbonds Number of bonds
     * @param bonds Pointer to bonds
     * @param natoms Number of atoms
     * @param atoms Pointer to atoms
     */
    void subgraf_pathfinder_run(subgraf_pathfinder *spf,
                                int nforbidden,
                                int *forbidden_orig,
                                int *nbonds, int **bonds,
                                int *natoms, int *atoms);

    /**
     * @brief Collects atom numbers along path
     *
     * @param spf Pointer to graph pathfinder data structure
     * @param nforbidden Number of edges forbidden for traversal
     * @param forbidden Pointer to forbidden edges
     * @param atnums Pointer to atom numbers
     * @return * int
     */
    int subgraf_pathfinder_collect_all(subgraf_pathfinder *spf,
                                       int nforbidden,
                                       int *forbidden,
                                       int *atnums);

    void CompAtomData_GetNumMapping(COMP_ATOM_DATA *adata, int *orig_num, int *curr_num);

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* _STRUTIL_H_ */
