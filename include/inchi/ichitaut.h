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


#ifndef _ICHITAUT_H_
#define _ICHITAUT_H_

#include "ichi_bns.h"
#include "extr_ct.h"

/*******************************************************
  ---     Header of tautomers groups      ---
  ---  Each entry is AT_TAUTOMER_HDR type ---
          number of tautomer groups (nNumTautGroups)
          index of the first tautomer group (#1)
          ...
          index of the last tautomer group  (#nNumTautGroups)
  --- end of the Header of tautomers groups description ---

  ---  One endpoint group description  ---
  ---  Each entry has AT_TAUTOMER type members ---
    <fixed portion (6 entries)>
          number of endpoints  (0=end of list)
          number of mobile groups, including number of negative charges (=num(H)+num(-))
          number of negative charges
          number of 1H atoms
          number of 2H (deuterium) atoms
          number of 3H (tritium) atoms
     <variable portion, sorted in ascending order>
          atom rank #1 (ascending order)
          ...
          atom rank #endpoints
  --- end of the endpoint group description ----

  ----------------------------------------------
     Note:
     In the following Linear CT Tautomer descriptions
     we assume the tautomeric groups and the endpoints
     within them have been properly sorted

  --------- Linear CT Tautomer description -----
    <for each tautomeric group>
        -- fixed length part, non-isotopic --
        number of endpoints        = t_group->nNumEndpoints
        number of mobile atoms     = t_group->num[0]
        ...
        number of negative charges = t_group->num[T_NUM_NO_ISOTOPIC-1]
        -- fixed length part, isotopic --
        number of T (3H)           = t_group->num[T_NUM_NO_ISOTOPIC]
        ...
        number of 1H               = t_group->num[T_NUM_NO_ISOTOPIC+T_NUM_ISOTOPIC-1]
        -- variable length part --
        rank of the first endpoint = nRank[t_group_info->nEndpointAtomNumber[t_group->nFirstEndpointAtNoPos]];
        ...
        rank of the last endpoint  = nRank[t_group_info->nEndpointAtomNumber[t_group->nFirstEndpointAtNoPos+t_group->nNumEndpoints-1]];

  --------- Linear CT Isotopic Tautomer description -----
    <for each isotopic tautomeric group>
        number of T (3H)           = t_group->num[T_NUM_NO_ISOTOPIC]
        ...
        number of 1H               = t_group->num[T_NUM_NO_ISOTOPIC+T_NUM_ISOTOPIC-1]
        t-group ordering number in the Linear CT Tautomer, starts from 1

***************************************************************/


#define           T_NUM_NO_ISOTOPIC 2
#define           T_NUM_ISOTOPIC    NUM_H_ISOTOPES /* was 2, now 3 */

#define           T_GROUP_HDR_LEN   (1+T_NUM_NO_ISOTOPIC /*+T_NUM_ISOTOPIC*/) /* LinearCTTautomer */

typedef AT_NUMB   AT_TAUTOMER;  /* LinearCTTautomer */

typedef AT_ISO_SORT_KEY  T_GROUP_ISOWT; /* must hold value up to T_GROUP_ISOWT_MULT^3-1 */
                                     /* similar to AT_ISO_SORT_KEY */
/* = num_1H + T_GROUP_ISOWT_MULT*(num_D + T_GROUP_ISOWT_MULT*num_T) */


#define T_GROUP_ISOWT_MULT   1024 /* (max. number of identical isotopic hydrogens in a taut. group) + 1 */
                                   /* changed from 256U 9-12-2003 */
                                   /* (similar to AT_ISO_SORT_KEY_MULT ) */
/* note: (long)T_GROUP_ISOWT should always be positive (have zero sign bit) */

typedef struct tagIsotopicTautomerGroup {
    AT_NUMB  tgroup_num;  /* ordering number of a tautomer group with isotopes > 0 */
    /*
    union {
        struct {
            AT_NUMB num_T;
            AT_NUMB num_D;
            AT_NUMB num_1H;
        };
        AT_NUMB num[T_NUM_ISOTOPIC];
    };
    */
    AT_NUMB  num[T_NUM_ISOTOPIC]; /* inverted order: num_T, num_D, num_1H */
} AT_ISO_TGROUP;

typedef enum tagTG_NumDA {  /* 2004-02-26 */
    TG_Num_dH, /* number of H donors that have only H (all single bonds) */
    TG_Num_dM, /* number of H donors that have (-)  (all single bonds) */
    TG_Num_aH, /* number of H acceptors that have H and no (-) (+a double bond) */
    TG_Num_aM, /* number of H acceptors that have (-) and possibly H (+ one double bond) */
    TG_Num_dO, /* number of H donors =C-OH or =C-O(-) */
    TG_Num_aO, /* number of H acceptors -C=O */
    TG_NUM_DA  /* number of elements in an array */
} TGNUMDA;

typedef struct tagTautomerGroup {
#if 0
union {
    struct {
           /*T_NUM_NO_ISOTOPIC = 2 elements:*/
        AT_RANK num_Mobile; /*Num_H+num_D+num_T+num_NegCharges*/
        AT_RANK num_NegCharges;
           /* T_NUM_ISOTOPIC = 3 elements*/
        AT_RANK num_T;      /*here the isotopic part (num+T_NUM_NO_ISOTOPIC) starts*/
        AT_RANK num_D;
        AT_RANK num_1H;
    };
    AT_RANK num[T_NUM_NO_ISOTOPIC+T_NUM_ISOTOPIC];  /*same size and meaning as num[] in T_ENDPOINT*/
};
#endif /* 0 */

    AT_RANK num[T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC];  /* same size and meaning as num[] in T_ENDPOINT */
                                                    /* isotopic inv. order: num_T, num_D, num_1H */
    AT_RANK num_DA[TG_NUM_DA];
    T_GROUP_ISOWT iWeight;   /* isotopic "weight" = T_GROUP_ISOWT_MULT*(T_GROUP_ISOWT_MULT*num_T + num_D)+num_1H; */
    AT_NUMB   nGroupNumber;  /* positive tautomer group ID = atom->endpoint */
    AT_NUMB   nNumEndpoints; /* number of the atom numbers in T_GROUP_INFO::nEndpointAtomNumber[] */
    AT_NUMB   nFirstEndpointAtNoPos; /* the first index of the atom number in T_GROUP_INFO::nEndpointAtomNumber[] */
} T_GROUP;

/* offsets/num_t_groups within T_GROUP_INFO::tGroupNumber  */
#define TGSO_CURR_ORDER  0  /* tGroupNumber:     current sorting order */
#define TGSO_SYMM_RANK   1  /* tSymmRank:     symmetry ranks (no isotopes) = min. ordering number > 0. */
#define TGSO_SYMM_IORDER 2  /* tiGroupNumber: isotopic symmetry rank sorting order */
#define TGSO_SYMM_IRANK  3  /* tiSymmRank:    isotopic symmetry ranks */
#define TGSO_TOTAL_LEN   4

/***************************************************/
/* flags for t_group_info->tni.bNormalizationFlags */
/***************************************************/
#define FLAG_PROTON_NPO_SIMPLE_REMOVED 0x0001
#define FLAG_PROTON_NP_HARD_REMOVED    0x0002
#define FLAG_PROTON_AC_SIMPLE_ADDED    0x0004
#define FLAG_PROTON_AC_SIMPLE_REMOVED  0x0008
#define FLAG_PROTON_AC_HARD_REMOVED    0x0010
#define FLAG_PROTON_AC_HARD_ADDED      0x0020
#define FLAG_PROTON_CHARGE_CANCEL      0x0040
#define FLAG_PROTON_SINGLE_REMOVED     0x0080

/* signifies tautomeric structure even though no t-group discovered */
#define FLAG_NORM_CONSIDER_TAUT      ( FLAG_PROTON_NPO_SIMPLE_REMOVED | \
                                       FLAG_PROTON_NP_HARD_REMOVED    | \
                                       FLAG_PROTON_AC_SIMPLE_ADDED    | \
                                       FLAG_PROTON_AC_SIMPLE_REMOVED  | \
                                       FLAG_PROTON_AC_HARD_REMOVED    | \
                                       FLAG_PROTON_AC_HARD_ADDED      | \
                                       FLAG_PROTON_SINGLE_REMOVED     | \
                                       FLAG_PROTON_CHARGE_CANCEL    )

#if ( FIX_N_MINUS_NORN_BUG == 1 )
#define FLAG_FORCE_SALT_TAUT         ( FLAG_PROTON_NP_HARD_REMOVED  | \
                                       FLAG_PROTON_AC_HARD_REMOVED  | \
                                       FLAG_PROTON_AC_HARD_ADDED    | \
                                       FLAG_PROTON_CHARGE_CANCEL    )
#else
/* force salt tautomerism exploration */
#define FLAG_FORCE_SALT_TAUT         ( FLAG_PROTON_NP_HARD_REMOVED  | \
                                       FLAG_PROTON_AC_HARD_REMOVED  | \
                                       FLAG_PROTON_AC_HARD_ADDED    )
#endif

typedef struct tagTautomerNormInfo {
    NUM_H       nNumRemovedExplicitH; /* keeps track of explicit H */
    NUM_H       nNumRemovedProtons;
    NUM_H       nNumRemovedProtonsIsotopic[NUM_H_ISOTOPES];
    INCHI_MODE   bNormalizationFlags;
} TNI;

/***************************************************/
/*      t_group_info definition                    */
/***************************************************/
typedef struct tagTautomerGroupsInfo {
    T_GROUP   *t_group;  /* max_num_t_groups elements */
    AT_NUMB   *nEndpointAtomNumber; /* nNumEndpoints elements; also see comments to T_GROUP */
    AT_NUMB   *tGroupNumber;
    int       nNumEndpoints;
    int       num_t_groups;
    int       max_num_t_groups;
    int       bIgnoreIsotopic;

    AT_NUMB   *nIsotopicEndpointAtomNumber; /* [0]: number of the following atoms; [1...]: non-tautomeric atoms that may have isotopic H */
    int       nNumIsotopicEndpoints;     /* allocated length of nIsotopicEndpointAtomNumber */
    NUM_H     num_iso_H[NUM_H_ISOTOPES]; /* isotopic H on tautomeric atoms and those in nIsotopicEndpointAtomNumber */

    TNI       tni;

    INCHI_MODE bTautFlags;
    INCHI_MODE bTautFlagsDone;
} T_GROUP_INFO;

#define CANON_FLAG_NO_H_RECANON           0x0001  /* iOther: second canonicalization of the no H structure */
#define CANON_FLAG_NO_TAUT_H_DIFF         0x0002  /* iOther: NoTautH eq. partition differs from NoH */
#define CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF 0x0004  /* iOther: eq. partition in isotopic only non-taut differs from non-isotopic */
#define CANON_FLAG_ISO_TAUT_DIFF          0x0008  /* iBase:  isotopic eq. partition in isotopic taut differs from non-isotopic taut */
#define CANON_FLAG_ISO_FIXED_H_DIFF       0x0010  /* iOther: isotopic eq. partition in fixed H non-taut differs from non-isotopic fixed H */

/* Note: rank of tautomer atom #i = Rank[nEndpointAtomNumber[i]] */
/*       for each tautomer atom group (t_group) t_group.nFirstEndpointAtNoPos */
/*       is the first index of the atom number in nEndpointAtomNumber[] */

typedef struct tagTautomerEndpoint {
    /*
    union {
        struct {
            AT_RANK num_Mobile; // Num_H+num_D+num_T+num_NegCharges
            AT_RANK num_NegCharges;
            AT_RANK num_T;
            AT_RANK num_D;
        };
        AT_RANK num[T_NUM_NO_ISOTOPIC+T_NUM_ISOTOPIC];    // same size and meaning as num[] in T_GROUP
    };
    */
    AT_RANK num[T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC];    /* same size and meaning as num[] in T_GROUP */
    AT_RANK num_DA[TG_NUM_DA];
    AT_NUMB nGroupNumber;
    AT_NUMB nEquNumber;  /* same for endpoints connected by alt paths */
    AT_NUMB nAtomNumber;
    /*AT_NUMB neighbor_index; */
} T_ENDPOINT;

typedef struct tagTautomerBondLocation {
    AT_NUMB nAtomNumber;
    AT_NUMB neighbor_index;
} T_BONDPOS;

typedef struct tagEndpointInfo {
    S_CHAR cMoveableCharge;
    S_CHAR cNeutralBondsValence;
    S_CHAR cMobile;
    S_CHAR cDonor;
    S_CHAR cAcceptor;
    S_CHAR cKetoEnolCode; /* 1 => carbon, 2 => oxygen */ /* post v.1 feature */
} ENDPOINT_INFO;


/* positive charge group (extended onium) */

#define CHARGED_CPOINT(X,i) ((X)[i].charge==1)

typedef struct tagChargeCandidate {
    AT_NUMB   atnumber;
    S_CHAR    type;
    S_CHAR    subtype;
} C_CANDIDATE;

typedef struct tagChargeGroup {
    AT_RANK   num[2]; /* [0]: number of (+), [1]: number atoms that have H, including H accessible through tautomerism */
    AT_RANK   num_CPoints;
    AT_NUMB   nGroupNumber;
    U_CHAR    cGroupType;
} C_GROUP;

typedef struct tagChargeGroupsInfo {
    C_GROUP *c_group;
    int     num_c_groups;
    int     max_num_c_groups;

    C_CANDIDATE *c_candidate;
    int          max_num_candidates;
    int          num_candidates; /* 0=>unimitialized, -1=>no candidates found */
} C_GROUP_INFO;

/* salts */
typedef struct tagSaltChargeCandidate {
    AT_NUMB   atnumber;
    S_CHAR    type;
    S_CHAR    subtype;
    AT_NUMB   endpoint; /* MAX_ATOMS+1 => found alt path to the candidate */
} S_CANDIDATE;

typedef struct tagSaltGroupInfo {
    S_CANDIDATE *s_candidate;
    int          max_num_candidates;
    int          num_candidates; /* 0=>unimitialized, -1=>no candidates found */
    int          num_other_candidates; /* num. non-"acidic O" candidates */
    int          num_p_only_candidates; /* num. non-tautomeric p-donor/acceptor candidates like -CH2-SH */
} S_GROUP_INFO;

/********************* ATOM_SIZES *******************************/
/* sizes of a component */
typedef struct tagAtomSizes {
    /* for tautomeric and non-tautomeric structures */
    int nMaxNumStereoAtoms; /* max. number of stereo atoms in isotopic case */
    int nMaxNumStereoBonds; /* max. number of stereo bonds in isotopic case */
    int num_isotopic_atoms;  /* includes atoms that have isotopic tautomeric H */
    int nLenCT;
    int nLenBonds;
    int nLenIsotopic;
    int nLenCTAtOnly;
    int nLenLinearCTStereoDble; /* max. number of stereo bonds in non-isotopic case */
    int nLenLinearCTStereoCarb; /* max. number of stereo atoms in non-isotopic case */
    /* int bHasIsotopicAtoms; */
    int bMayHaveStereo;

    int bIgnoreIsotopic;

    /* tautomeric structure only; zeroes in non-tautomeric */
    int nLenLinearCTTautomer;
    int nLenLinearCTIsotopicTautomer;
    int bHasIsotopicTautGroups;
    int nLenIsotopicEndpoints;
} ATOM_SIZES;


typedef struct tagDfsPath {
    AT_RANK       at_no;
    /*AT_RANK       nDfsLevel;*/
    U_CHAR        bond_type;
    S_CHAR        bond_pos;
} DFS_PATH;


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

    int is_centerpoint_elem( U_CHAR el_number );
    int is_centerpoint_elem_strict( U_CHAR el_number );
#if ( KETO_ENOL_TAUT == 1 )
    int is_centerpoint_elem_KET( U_CHAR el_number );
#endif
    int bIsCenterPointStrict( inp_ATOM *atom, int iat );
    int nGetEndpointInfo( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif );

#if ( TAUT_PT_22_00 == 1 )
    int nGetEndpointInfo_PT_22_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
#endif  
#if ( TAUT_PT_16_00 == 1 )  
    int nGetEndpointInfo_PT_16_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
#endif
#if ( TAUT_PT_06_00 == 1 )  
    int nGetEndpointInfo_PT_06_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
#endif
#if ( TAUT_PT_39_00 == 1 )  
    int nGetEndpointInfo_PT_39_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
#endif
#if ( TAUT_PT_13_00 == 1 )  
    int nGetEndpointInfo_PT_13_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
#endif
#if ( TAUT_PT_18_00 == 1 )  
    int nGetEndpointInfo_PT_18_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
#endif   

#if ( KETO_ENOL_TAUT == 1 )
    int nGetEndpointInfo_KET( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif );
#endif
    void AddAtom2DA( AT_RANK num_DA[], inp_ATOM *atom, int at_no, int bSubtract );
    int AddAtom2num( AT_RANK num[], inp_ATOM *atom, int at_no, int bSubtract );
    int AddEndPoint( T_ENDPOINT *pEndPoint, inp_ATOM *at, int iat );
    int bHasAcidicHydrogen( inp_ATOM *at, int i );
    int bHasOtherExchangableH( inp_ATOM *at, int i );
    int bHasAcidicMinus( inp_ATOM *at, int i );


    int nGet15TautIn6MembAltRing( struct tagCANON_GLOBALS *pCG,
                                  inp_ATOM *atom,
                                  int nStartAtom,
                                  AT_RANK  *nBfsTreePos,
                                  DFS_PATH *DfsPath,
                                  int nMaxLenBfsTree,
                                  T_ENDPOINT *EndPoint,
                                  int nMaxNumEndPoint,
                                  T_BONDPOS  *BondPos,
                                  int nMaxNumBondPos,
                                  int *pnNumEndPoint,
                                  int *pnNumBondPos,
                                  struct BalancedNetworkStructure *pBNS,
                                  struct BalancedNetworkData *pBD,
                                  int num_atoms );

    int nGet12TautIn5MembAltRing( struct tagCANON_GLOBALS *pCG,
                                  inp_ATOM *atom, int nStartAtom,
                                  int nStartAtomNeighbor,
                                  AT_RANK  *nBfsTreePos,
                                  DFS_PATH *DfsPath,
                                  int nMaxLenBfsTree,
                                  T_ENDPOINT *EndPoint,
                                  int nMaxNumEndPoint,
                                  T_BONDPOS  *BondPos,
                                  int nMaxNumBondPos,
                                  int *pnNumEndPoint,
                                  int *pnNumBondPos,
                                  struct BalancedNetworkStructure *pBNS,
                                  struct BalancedNetworkData *pBD,
                                  int num_atoms );

    int nGet14TautIn7MembAltRing( struct tagCANON_GLOBALS *pCG,
                                  inp_ATOM *atom,
                                  int nStartAtom,
                                  int nStartAtomNeighbor,
                                  int nStartAtomNeighborEndpoint,
                                  int nStartAtomNeighborNeighborEndpoint,
                                  AT_RANK  *nDfsPathPos,
                                  DFS_PATH *DfsPath,
                                  int nMaxLenDfsPath,
                                  T_ENDPOINT *EndPoint,
                                  int nMaxNumEndPoint,
                                  T_BONDPOS  *BondPos,
                                  int nMaxNumBondPos,
                                  int *pnNumEndPoint,
                                  int *pnNumBondPos,
                                  struct BalancedNetworkStructure *pBNS,
                                  struct BalancedNetworkData *pBD,
                                  int num_atoms );

    int nGet14TautIn5MembAltRing( struct tagCANON_GLOBALS *pCG,
                                  inp_ATOM *atom,
                                  int nStartAtom,
                                  int nStartAtomNeighbor,
                                  int nStartAtomNeighborEndpoint,
                                  int nStartAtomNeighborNeighborEndpoint,
                                  AT_RANK  *nDfsPathPos,
                                  DFS_PATH *DfsPath,
                                  int nMaxLenDfsPath,
                                  T_ENDPOINT *EndPoint,
                                  int nMaxNumEndPoint,
                                  T_BONDPOS  *BondPos,
                                  int nMaxNumBondPos,
                                  int *pnNumEndPoint,
                                  int *pnNumBondPos,
                                  struct BalancedNetworkStructure *pBNS,
                                  struct BalancedNetworkData *pBD,
                                  int num_atoms );

    int nGet15TautInAltPath( struct tagCANON_GLOBALS *pCG,
                                  inp_ATOM *atom,
                                  int nStartAtom, AT_RANK  *nDfsPathPos,
                                  DFS_PATH *DfsPath,
                                  int nMaxLenDfsPath,
                                  T_ENDPOINT *EndPoint,
                                  int nMaxNumEndPoint,
                                  T_BONDPOS  *BondPos,
                                  int nMaxNumBondPos,
                                  int *pnNumEndPoint,
                                  int *pnNumBondPos,
                                  struct BalancedNetworkStructure *pBNS,
                                  struct BalancedNetworkData *pBD,
                                  int num_atoms );


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif    /* _ICHITAUT_H_ */
