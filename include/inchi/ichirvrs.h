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


#ifndef _ICHIRVRS_H_
#define _ICHIRVRS_H_

#include "ichimain.h"
#include "ichiring.h"

#define ICHICONST      const

#define RI_ERR_ALLOC   (-1)
#define RI_ERR_SYNTAX  (-2)
#define RI_ERR_PROGR   (-3)
#define RI_ERR_EOL     (-4)
#define RI_ERR_EOF     (0)

#define RI_ERR_MISMATCH (-9)

#define NO_VALUE_INT    9999
#define NOT_READ_INT    9998 /* has not been read yet */

#define VALUE_OCTET    8     /* number of electrons in a full shell */

#define INC_EDGE_LIST_DEFAULT  64

typedef struct tagXYZCoord {
    double xyz[3];
} XYZ_COORD;

typedef struct tagStructRestoreMode {
    int bMetalAddFlower;       /* 1 => allow adjustable metal valence and charge; 0=> use std charge/valence */
    /* the following three apply only if bMetalAddFlower = 1 */
    int nMetalMinBondOrder;    /* edge_flow=f means bond order=cMetalMinBondOrder+f */
    int nMetalInitEdgeFlow;    /* one bond contribution to metal's (st-cap - st-flow) = */
                               /*   (nMetalInitBondOrder-nMetalMinBondOrder) - nMetalInitEdgeFlow */
    int nMetalInitBondOrder;   /* >= nMetalMinBondOrder + nMetalInitEdgeFlow */
    /* same for metal-endpoint bonds */
    int nMetal2EndpointMinBondOrder;
    int nMetal2EndpointInitBondOrder;
    int nMetal2EndpointInitEdgeFlow;
    int nMetalFlowerParam_D;      /* additional edge capacity for canceling radicals */
    int nMetalMaxCharge_D;        /* cap and/or flow for metal charge group */
    int bStereoRemovesMetalFlag;  /* 1=> treat stereogenic atoms and atoms connected by a stereo bond as non-metals */

    int bFixStereoBonds;       /* 1=> forbid stereogenic double bonds from changing */
} SRM;

typedef struct tagReversedInChI {
    PINChI2     *pINChI[INCHI_NUM];
    PINChI_Aux2 *pINChI_Aux[INCHI_NUM];
    int          num_components[INCHI_NUM];
    int          nRetVal;
} REV_INCHI;

/**************************************/
#define BFS_Q_CLEAR  (-1)
#define BFS_Q_FREE   (-2)
typedef struct tagBfsQueue {
    QUEUE   *q;
    AT_RANK *nAtomLevel;
    S_CHAR  *cSource;
    int      num_at;
    AT_RANK  min_ring_size;  /* 8 => detect 7-member and smaller rings */
} BFS_Q;
/**************************************/

#define EXTRACT_STRUCT_NUMBER 1

/* additional Mobile-H parities to be added to Fixed-H parities  */
/* This allows to set parities that exist in Mobile-H layer only */
typedef struct tagInpAtomAddParities {
    /* cml 0D parities */
    S_CHAR        bUsed0DParity;          /* bit=1 => stereobond; bit=2 => stereocenter */
    /* cml tetrahedral parity */
    S_CHAR        p_parity;
    AT_NUMB       p_orig_at_num[MAX_NUM_STEREO_ATOM_NEIGH];
    /* cml bond parities */
    S_CHAR        sb_ord[MAX_NUM_STEREO_BONDS];  /* stereo bond/neighbor ordering number, starts from 0 */
    /* neighbors on both sides of stereobond have same sign=> trans/T/E, diff. signs => cis/C/Z */
    S_CHAR        sn_ord[MAX_NUM_STEREO_BONDS]; /* ord. num. of the neighbor adjacent to the SB; starts from 0;
                                                   -1 means removed explicit H */
    /* neighbors on both sides of stereobond have same parity => trans/T/E/2, diff. parities => cis/C/Z/1 */
    S_CHAR        sb_parity[MAX_NUM_STEREO_BONDS];
    AT_NUMB       sn_orig_at_num[MAX_NUM_STEREO_BONDS]; /* orig. at number of sn_ord[] neighbors */
} inp_ATOM_STEREO;

#define FIX_STEREO_BOND_ORDER  0   /* 1=> fix stereobonds; treat metal as non-metal if it is stereogenic or has a stereobond */

#define METAL_FREE_CHARGE_VAL  1  /* 1=> allow free changing charges/valences of metals; initial bond order=0 or 1 */
#define ALLOW_METAL_BOND_ZERO  1  /* 1=> allow zero flow (bobd order) to metals */

#if ( ALLOW_METAL_BOND_ZERO == 1 )
/* INIT_METAL_BOND_ZERO=1 => INIT_METAL_BOND_FLOW=0 */
#define INIT_METAL_BOND_ZERO   0  /* 1=> initialize zero order bond to metals */
#define INIT_METAL_BOND_FLOW   1  /* 1=> init flow=1, 0 => init. flow = 0 */
#else
#define INIT_METAL_BOND_ZERO   0
#define INIT_METAL_BOND_FLOW   0  /* always 0 */
#endif

#define I2A_FLAG_FIXEDH  0x0001
#define I2A_FLAG_RECMET  0x0002


#define ATYPE_H   1
#define ATYPE_Na  2
#define ATYPE_Mg  3
#define ATYPE_B   4
#define ATYPE_C   5
#define ATYPE_N   6
#define ATYPE_O   7
#define ATYPE_Cl  8


/* first  bonds valence for charge = c is cValence[0][c+VAL_BASE]; VAL_MIN_CHARGE <= c <= VAL_MAX_CHARGE */
/* second bonds valence for charge = c is cValence[1][c+VAL_BASE]; VAL_MIN_CHARGE <= c <= VAL_MAX_CHARGE */
/* total number of valences is 2 = VAL_NUMBER */
/* neutral bond valence orders are cValence[0][VAL_NEUTR_ORDER], cValence[1][VAL_NEUTR_ORDER] */
#define VAL_BASE         ( 1)
#define VAL_MIN_CHARGE   (-1)
#define VAL_MAX_CHARGE   ( 1)
#define VAL_NUMBER       ( 2)
#define VAL_NEUTR_ORDER  (VAL_MAX_CHARGE-VAL_MIN_CHARGE+1)
#define VAL_LENGTH       (VAL_NEUTR_ORDER+1)
#define VAL_NEGAT_CHARGE    0
#define VAL_NEUTR_CHARGE    1
#define VAL_POSIT_CHARGE    2

typedef struct tagAtomIonPrperies {
    /* char cValence[VAL_NUMBER][VAL_LENGTH]; */ /* ordering numbers of minimal valence, 0-based */
    char   cDoNotAddH;         /* InChI does not add H to this element */
    char   cMetal;             /* the element is a metal */
    char   cNumBondsToMetal;   /* number of bonds to metal */
    char   cInitFlowToMetal;    /* sum of init flow to metal atoms */
    char   cInitValenceToMetal; /* sum of init adjusted bond orders to metal atoms */
    char   cInitOrigValenceToMetal; /* sum of init bond orders to metal atoms */
    char   cMaxFlowToMetal;    /* max total edge flow to metal atoms */
    char   cInitFreeValences;  /* number of 'dots' to connect; charges are marked separately */
    S_CHAR cInitCharge;        /* initial charge on the atom (not included in ChargeStruct */
    char   cNumValenceElectrons;
    char   cPeriodicRowNumber;
    char   cMinRingSize;       /* min ring size for atoms that have 2 bonds only */
    U_CHAR cPeriodicNumber;  /* number in Periodic Table of elements */
    S_CHAR cnListIndex;        /* (index in the cnList) + 1; 0 => none */
    int    nCMinusGroupEdge;   /* (index of the edge to the atom's (-) group) + 1 */
    int    nCPlusGroupEdge;    /* (index of the edge to the atom's (+) group) + 1 */
    int    nMetalGroupEdge;    /* index of the edge to the atom's M-group + 1 */
    int    nTautGroupEdge;     /* index of the edge from the atom to the t-group + 1 */
} VAL_AT;


/****************************************************************************/
#define INI_NUM_TCGROUPS 16
#define INC_NUM_TCGROUPS 16

typedef enum tagTgRestoreFlags {
    TGRF_MINUS_FIRST = 1
} TGRF;
typedef struct tagTCGroup {
    int type;    /* group type */
    int ord_num; /* ordering number within the type, typically t-group number */
    int num_edges;
    /* charge group specific */
    int st_cap;
    int st_flow;
    int edges_cap;
    int edges_flow;
    int nVertexNumber; /* group vertex number; 0 = unassigned */
    int nForwardEdge;  /* edge index: from c-group to central Y-connecting vertex
                          or from supergroup to (+/-) vertex; 0 => unassigned */
    int nBackwardEdge; /* edge index: from central Y-connecting vertex
                          to supergroup; 0 => unassigned */
    /* tautomeric group specific */
    short tg_num_H;      /* number of H in a tautomeric group */
    short tg_num_Minus;  /* negative charge on t-group */
    Vertex tg_set_Minus;  /* the vertex+1 that has to have (-) */
    short tg_RestoreFlags; /* Set (-) to first memberst of a t-group (usually, N) */
} TC_GROUP;

typedef enum tagTCGroupTypes {
    TCG_None = -1, /* so far only ord=0 is used */
                   /* group type      ord */
                   TCG_Plus0 = 0, /* BNS_VT_C_POS      0 */
                   TCG_Plus1,     /* BNS_VT_C_POS      1 */
                   TCG_Minus0,    /* BNS_VT_C_NEG      0 */
                   TCG_Minus1,    /* BNS_VT_C_NEG      1 */
                   TCG_Plus_C0,   /* BNS_VT_C_POS_C    0 */
                   TCG_Plus_C1,   /* BNS_VT_C_POS_C    1 */
                   TCG_Minus_C0,  /* BNS_VT_C_NEG_C    0 */
                   TCG_Minus_C1,  /* BNS_VT_C_NEG_C    1 */
                   TCG_Plus_M0,   /* BNS_VT_C_POS_M    0 */
                   TCG_Plus_M1,   /* BNS_VT_C_POS_M    1 */
                   TCG_Minus_M0,  /* BNS_VT_C_NEG_M    0 */
                   TCG_Minus_M1,  /* BNS_VT_C_NEG_M    1 */
                   TCG_MeFlower0, /* BNS_VT_M_GROUP    0 */  /* base */
                   TCG_MeFlower1, /* BNS_VT_M_GROUP    1 */
                   TCG_MeFlower2, /* BNS_VT_M_GROUP    2 */
                   TCG_MeFlower3, /* BNS_VT_M_GROUP    3 */

                   TCG_Plus,      /* BNS_VT_C_POS_ALL  0 */
                   TCG_Minus,     /* BNS_VT_C_NEG_ALL  0 */

                   NUM_TCGROUP_TYPES /* number of group types */
}TCGR_TYPE;

typedef struct tagAllTCGroups {
    TC_GROUP *pTCG;
    int      num_tc_groups; /* number of charge groups and metal-flower vertices */
    int      max_tc_groups; /* number of allocated of pTCG[] elements */
    int      nGroup[NUM_TCGROUP_TYPES];  /* tagTCGroupTypes */
    int      nVertices;  /* total number of vertices */
    int      nEdges;     /* total number of edges */
    int      nAddIedges; /* additional increase of number of iedges for edge switching to another group */
    int      num_atoms;  /* number of atoms */
    int      num_bonds;  /* number of bonds */
    int      num_tgroups;       /* number t-groups */
    int      num_tgroup_edges;  /* number of edges to t-groups */
    /* charges */
    int      tgroup_charge;     /* total charge of all t-groups */
    int      charge_on_atoms;   /* charge permanently sitting on atoms */
    int      added_charge;      /* charge added to the c-groups */
    int      total_charge;      /* total charge of the component */
    int      total_electrons;   /* total number of electrons on all atoms */
    int      total_electrons_metals; /* total number of electrons on unbonded metals */

    int      num_metal_atoms;   /* number of metal atoms */
    int      num_metal_bonds;   /* number of atom-metal bonds */
    /* excess_charge = total_charge - added_charge - tgroup_charge: add to metals etc. */

    int      nEdge4charge;      /* edge used to add charges; neighbor1=supercharge, another = (+/-) vertex */
    int      nEdgePlus;         /* edge to (+) supergroup; 0 means none */
    int      nEdgeMinus;        /* edge to (-) supergroup; 0 means none */
    int      iComponent;        /* component number */
    int      iAtNoOffset;       /* first atom number -- always 0 for now */
} ALL_TC_GROUPS;

/**************************************/
#define EDGE_LIST_CLEAR  (-1)
#define EDGE_LIST_FREE   (-2)

typedef struct tagEdgeList {
    int num_alloc;
    int num_edges;
    EdgeIndex *pnEdges;
} EDGE_LIST;
/**************************************/

#define BOND_MARK_STEREO 0x10
#define BOND_TYPE_STEREO (BOND_TYPE_SINGLE | BOND_MARK_STEREO)

/* local */
#define RESET_EDGE_FORBIDDEN_MASK 0

#define TREAT_ATOM_AS_METAL 99


/************************************************************************************/
typedef struct tagChargeValence {
    int nValence;
    int nCharge;
    int nValenceOrderingNumber;
} CHARGE_VAL;

#define MY_CONST  const
/*************************************************************************************/
typedef struct tagChargeChangeCandidate {
    Vertex iat;
    char   num_bonds;
    char   chem_valence;
    char   cMetal;
    char   cNumBondsToMetal;
    char   cNumValenceElectrons;
    char   cPeriodicRowNumber;
    char   cNumChargeStates;
    U_CHAR el_number;
} CC_CAND;

typedef struct tagOneComponentRemovedAndExchangeableH {
    NUM_H   nNumRemovedProtons;
    NUM_H   nNumRemovedIsotopicH[NUM_H_ISOTOPES]; /* isotopic H that may be exchanged and considered
                                                     randomly distributed, including removed protons */
} COMPONENT_REM_PROTONS;

typedef struct tagRemovedAndExchangeableH {
    /* totals for Mobile-H layer */
    NUM_H   nNumRemovedProtons;
    NUM_H   nNumRemovedIsotopicH[NUM_H_ISOTOPES]; /* isotopic H that may be exchanged and considered
                                                     randomly distributed, including removed protons */
    /* for individual components from comparing Fixed-H vs Mobile-H formulas; NULL if not available */
    COMPONENT_REM_PROTONS *pNumProtons;
} REM_PROTONS;

typedef struct tagInputInChI {
    INChI      *pInpInChI[INCHI_NUM][TAUT_NUM];
    int         nNumComponents[INCHI_NUM][TAUT_NUM];
    REM_PROTONS nNumProtons[INCHI_NUM][TAUT_NUM];
    int         s[INCHI_NUM][TAUT_NUM][2];
                                /* s[0=non-iso, 1=iso] = 0,1,2,3 <= regular /s; -1=> "/s" (empty) */
    long        num_inp;
    inp_ATOM *atom;                /* the whole restored structure made out of all components */
    int         num_atoms;        /* number of atoms including explicit H */
    int         num_explicit_H; /* number of explicit H in the atom */
    INCHI_MODE  CompareInchiFlags[INCHI_NUM][TAUT_NUM];
    /* v. 1.05 extensions */
    OAD_Polymer *polymer;
    OAD_V3000    *v3000;
    int valid_polymer;
} InpInChI;

typedef struct tagStructFromInChI {
    /* InChI component -> Structure result */
    inp_ATOM        *at;  /* length = num_atoms + num_deleted_H, zero pint struct for BNS->struct conversion */
    inp_ATOM_STEREO *st;  /* additional stereo that exists only in Mobile-H layer */
    inp_ATOM        *at2; /* length = num_atoms + num_deleted_H, the conversion result */

    /* information from InChI only */
    T_GROUP_INFO  ti; /* from original InChI[0] if Mobile-H from the beginning or later from InChI[1] if Fixed-H  */
    AT_NUMB      *endpoint; /* from original InChI[1] in case of Fixed-H only */
    S_CHAR       *fixed_H;  /* from original InChI[0] in case of Fixed-H only */
    XYZ_COORD    *pXYZ;
    int           num_atoms;
    int           num_deleted_H; /* if requested and Fixed-H InChI is available */
    int           nNumRemovedProtonsMobHInChI; /* number of protons removed from Mobile-H struct in original InChI */
    S_CHAR        charge;
    char          bIsotopic;

    /* InChI -> Structure conversion parms and intermediate data */
    BN_STRUCT *pBNS;
    BN_DATA   *pBD;
    ICHICONST SRM       *pSrm;

    /* InChI layer to reverse */
    char        bMobileH;
    char        iINCHI;
    char        bFixedHExists; /* fixed-H InChI exists or not */

    /* InChI -> Struct component -> Full InChI result (both disconnected and connected if exist) */
    REV_INCHI  RevInChI;
    int        nRemovedProtonsByNormFromRevrs; /* number of H(+) removed by normalization after
                                                 Struct Restore and before Add/Remove Protons */
    int        nNumRemovedProtonsByRevrs;       /* number of H(+) removed by the reconstruction,
                                                   before Add/Remove Protons, only from TAUT_YES */

    int        bExtract; /* for debugging */

    /* single component InChI calculation */
    INChI         *pOneINChI[TAUT_NUM];       /* InChI of restored structure */
    INChI_Aux     *pOneINChI_Aux[TAUT_NUM];
    INP_ATOM_DATA *pOne_norm_data[TAUT_NUM];  /* normalized restored structure */
    S_CHAR        *pOne_fixed_H;              /* !!! from normalized restored structure in case of Fixed-H only */
    T_GROUP_INFO   One_ti;                    /* t-groups of normalized canonicalized restored structure */
    int            nOneINChI_bMobileH;        /* type of restored structure InChI */
    int            nNumRemovedProtons;        /* =0 for Fixed-H, = num. removed protons in case of Mobile-H InChI */
                                              /* in case of Fixed-H processing see pStruct->One_ti.tni.nNumRemovedProtons */
    AT_NUMB       *nAtno2Canon[TAUT_NUM];     /* nAtno2Canon[restored_at_no][*] = (atom canon number in restored struct)-1*/
    AT_NUMB       *nCanon2Atno[TAUT_NUM];     /* nCanon2Atno[(atom canon number in restored struct)-1][*] = restored_at_no; */

    int            nError;
    /* other parms */
    char        iInchiRec; /* index in the original InChI array */
    char        iMobileH;  /* index in the original InChI array */
    char        bDeleted;  /* InChI component marked as Deleted, means a proton in Mobile-H layer */
    /* struct. ordering number "Structure: nnn" if present */
    long        num_inp_actual;

    /* utility data */
    BFS_Q   *pbfsq;
    VAL_AT  *pVA;

    int   nLink; /* same as in INChI */
    int   bPostProcessed; /* recalculate after add/remove protons */

    /* TAUT_YES layer charges */
    int nChargeRevrs;  /* component charge of the reconstructed structure, TAUT_YES layer */
    int nChargeInChI;  /* component charge from the original InChI, TAUT_YES layer */

    int n_zy;
    int n_pzz;
} StrFromINChI;


#define EL_TYPE_O     0x0001
#define EL_TYPE_S     0x0002
#define EL_TYPE_N     0x0004
#define EL_TYPE_P     0x0008
#define EL_TYPE_C     0x0010
#define EL_TYPE_X     0x0020 /* any not metal */
#define EL_TYPE_MASK  0x003f
#define EL_TYPE_OSt   0x0100 /* terminal -OH, -O(-), -SH, -S(-), ... from fix_special_bonds(...) */
#define EL_TYPE_PT    0x0200 /* may be a tautomeric endpoint */

/* the atom to which the node is attached has number 1; added atoms have numbers 2,3,... */
#define MAX_CN_VAL 3
typedef struct tagVertCapFlow {
    S_SHORT type;
    S_CHAR  cap;
    S_CHAR  flow;
    S_CHAR  valence;
} VCF;
typedef struct tagEdgeCapFlow {
    S_SHORT neigh;
    S_CHAR  cap;
    S_CHAR  bForbiddenEdge;
    S_CHAR  flow;
} ECF;
typedef struct tagChargeNodes {
    VCF v;
    ECF e[MAX_CN_VAL];
} C_NODE;


#define cn_bits_N  1  /* Neutral: charge =  0 */
#define cn_bits_P  2  /* Plus 1:  charge = +1 */
#define cn_bits_M  4  /* Minus 1: charge = -1 */
#define cn_bits_shift 3
#define MAX_NUM_CN_BITS 4
#define MAKE_CN_BITS(A, B, C, D ) (( (( ((D) << cn_bits_shift | (C)) << cn_bits_shift ) | (B)) << cn_bits_shift ) | (A))

#define cn_bits_PNPN MAKE_CN_BITS(cn_bits_P, cn_bits_N, cn_bits_P, cn_bits_N)
#define cn_bits_NPNP MAKE_CN_BITS(cn_bits_N, cn_bits_P, cn_bits_N, cn_bits_P)
#define cn_bits_NPN MAKE_CN_BITS(cn_bits_N, cn_bits_P, cn_bits_N, 0)
#define cn_bits_PNP MAKE_CN_BITS(cn_bits_P, cn_bits_N, cn_bits_P, 0)
#define cn_bits_MNP  MAKE_CN_BITS(cn_bits_M, cn_bits_N, cn_bits_P, 0)
#define cn_bits_PNM  MAKE_CN_BITS(cn_bits_P, cn_bits_N, cn_bits_M, 0)
#define cn_bits_EN  MAKE_CN_BITS(cn_bits_P | cn_bits_M, cn_bits_N, 0, 0)
#define cn_bits_NMN  MAKE_CN_BITS(cn_bits_N, cn_bits_M, cn_bits_N, 0)
#define cn_bits_NE  MAKE_CN_BITS(cn_bits_N, cn_bits_P | cn_bits_M, 0, 0)
#define cn_bits_NEN  MAKE_CN_BITS(cn_bits_N, cn_bits_M | cn_bits_N, cn_bits_N, 0)
#define cn_bits_NP  MAKE_CN_BITS(cn_bits_N, cn_bits_P, 0, 0)
#define cn_bits_PN  MAKE_CN_BITS(cn_bits_P, cn_bits_N, 0, 0)
#define cn_bits_NM  MAKE_CN_BITS(cn_bits_N, cn_bits_M, 0, 0)
#define cn_bits_MN  MAKE_CN_BITS(cn_bits_M, cn_bits_N, 0, 0)
#define cn_bits_P_  MAKE_CN_BITS(cn_bits_P, 0, 0, 0)
#define cn_bits_M_  MAKE_CN_BITS(cn_bits_M, 0, 0, 0)
#define cn_bits_N_  MAKE_CN_BITS(cn_bits_N, 0, 0, 0)
#define cn_bits_Me (-1)

#define cnListIndexMe (17)  /* index of {cnMe, cn_bits_Me,... } element of cnList[] */

extern const int cnListNumEl;  /* number of elements in cnList[] */

typedef struct tagChargeNodeList {
    MY_CONST C_NODE *pCN;
    int     bits;
    int     nInitialCharge;
    int     len;
} CN_LIST;

extern MY_CONST CN_LIST cnList[];

/************************ fixed H comparison ******************************************************/
#define MAX_DIFF_FIXH 256
#define MAX_DIFF_MOBH 256
typedef struct tagAtomsCmpTwoFixedH {
    AT_NUMB endptInChI;
    AT_NUMB endptRevrs;
    AT_NUMB atomNumber;
    U_CHAR nValElectr;
    U_CHAR nPeriodNum;
    S_CHAR nFixHInChI;
    S_CHAR nFixHRevrs;
    S_CHAR nMobHInChI;
    S_CHAR nMobHRevrs;
    S_CHAR nNumHRevrs;
    S_CHAR nAtChargeRevrs;
    S_CHAR nValue; /* flag(s) */
} CMP2FHATOMS;

typedef struct tagStructCmpTwoFixedH {
    CMP2FHATOMS c2at[MAX_DIFF_FIXH];
    short  len_c2at;
    short  nNumRemHInChI;
    short  nNumRemHRevrs;
    short  nNumTgInChI;
    short  nNumTgRevrs;
    short  nNumEndpInChI;
    short  nNumEndpRevrs;
    short  nNumTgDiffMinus; /* number of would-be-identical t-groups that have different number of (-) */
    short  nNumTgDiffH;     /* number of would-be-identical t-groups that have different number of H */
    short  nNumTgMInChI;    /* number of (-) in orig. InChI t-groups */
    short  nNumTgHInChI;    /* number of H in orig. InChI t-groups */
    short  nNumTgMRevrs;    /* number of (-) in reversed structure t-groups */
    short  nNumTgHRevrs;    /* number of H in reversed structure t-groups */
    S_CHAR nChargeFixHInChI;
    S_CHAR nChargeMobHInChI;
    S_CHAR nChargeFixHRevrs;
    S_CHAR nChargeMobHRevrs;
    S_CHAR nChargeFixHRevrsNonMetal; /* charge does not include charges on metals */
    S_CHAR nChargeMobHRevrsNonMetal; /* charge does not include charges on metals */
    char   bFixedHLayerExistsRevrs;
    char   bHasDifference;
    U_CHAR nNumDiffMobH;
} CMP2FHINCHI;

/************************ Mobile H comparison *********************************************/
typedef struct tagAtomsCmpTwoMobileH {
    AT_NUMB endptInChI;
    AT_NUMB endptRevrs;
    AT_NUMB atomNumber;
    U_CHAR nValElectr;
    U_CHAR nPeriodNum;
    S_CHAR nMobHInChI; /* number of H on the atom in the orig. InChI */
    S_CHAR nMobHRevrs; /* number of H on the atom in InChI from the reconstructed structure */
    S_CHAR nNumHRevrs; /* number of H on the atom in the being reconstructed structure */
    S_CHAR nAtChargeRevrs;
    S_CHAR nValue;         /* flag(s) */
} CMP2MHATOMS;

typedef struct tagStructCmpTwoMobileH {
    CMP2MHATOMS c2at[MAX_DIFF_FIXH];
    short len_c2at;
    short nNumRemHInChI;
    short nNumRemHRevrs;
    short nNumTgInChI;
    short nNumTgRevrs;
    short nNumEndpInChI;
    short nNumEndpRevrs;
    short nNumTgDiffMinus;   /* number of would-be-identical t-groups that have different number of (-) */
    short nNumTgDiffH;       /* number of would-be-identical t-groups that have different number of H */

    short nNumTgMInChI;      /* number of (-) in orig. InChI t-groups */
    short nNumTgHInChI;      /* number of H in orig. InChI t-groups */
    short nNumTgOInChI;      /* number of tautomeric O,S,Se in orig. InChI t-groups */
    short nNumTgNInChI;      /* number of tautomeric N in orig. InChI t-groups */

    short nNumTgMRevrs;      /* number of (-) in reversed structure t-groups */
    short nNumTgHRevrs;      /* number of H in reversed structure t-groups */
    short nNumTgORevrs;      /* number of tautomeric O,S,Se in reversed structure t-groups */
    short nNumTgNRevrs;       /* number of tautomeric N in reversed structure t-groups */

    short nNumTgOMinusRevrs;   /* number of -O(-)  on endpoints found in restored structure */
    short nNumTgOHRevrs;       /* number of -OH    on endpoints found in restored structure */
    short nNumTgDBORevrs;      /* number of =O     on endpoints found in restored structure */
    short nNumTgNMinusRevrs;   /* number of -N(-)- on endpoints found in restored structure */
    short nNumTgNHMinusRevrs;  /* number of -NH(-) on endpoints found in restored structure */
    short nNumTgNHRevrs;       /* number of -NH-   on endpoints found in restored structure */
    short nNumTgNH2Revrs;      /* number of -NH2   on endpoints found in restored structure */
    short nNumTgDBNHRevrs;     /* number of =NH    on endpoints found in restored structure */
    short nNumTgDBNMinusRevrs; /* number of =N(-)  on endpoints found in restored structure */
    short nNumTgDBNRevrs;      /* number of =N-    on endpoints found in restored structure */

    short nNumTgOMinusInChI;   /* number of -O(-)  on endpoints according to original InChI */
    short nNumTgOHInChI;       /* number of -OH    on endpoints according to original InChI */
    short nNumTgDBOInChI;      /* number of =O     on endpoints according to original InChI */
    short nNumTgNMinusInChI;   /* number of -N(-)- on endpoints according to original InChI */
    short nNumTgNHMinusInChI;  /* number of -NH(-) on endpoints according to original InChI */
    short nNumTgNHInChI;       /* number of -NH-   on endpoints according to original InChI */
    short nNumTgNH2InChI;      /* number of -NH2   on endpoints according to original InChI */
    short nNumTgDBNHInChI;     /* number of =NH    on endpoints according to original InChI */
    short nNumTgDBNMinusInChI; /* number of =N(-)  on endpoints according to original InChI */
    short nNumTgDBNInChI;      /* number of =N-    on endpoints according to original InChI */

    S_CHAR nChargeMobHInChI;
    S_CHAR nChargeMobHRevrs;
    S_CHAR nChargeMobHRevrsNonMetal; /* charge does not include charges on metals; later add ion pairs rejection */
    char   bFixedHLayerExistsRevrs;
    char   bHasDifference;
    U_CHAR nNumDiffMobH;
} CMP2MHINCHI;


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


struct tagCANON_GLOBALS;
struct tagINCHI_CLOCK;
int OneInChI2Atom( struct tagINCHI_CLOCK *ic,
                   struct tagCANON_GLOBALS *pCG,
                   ICHICONST INPUT_PARMS *ip,
                   STRUCT_DATA *sd,
                   const char *szCurHdr,
                   long num_inp,
                   StrFromINChI *pStruct,
                   int iComponent,
                   int iAtNoOffset,
                   int bHasSomeFixedH,
                   INChI *pInChI[] );

int get_sp_element_type( int nPeriodicNumber, int *nRow );

int get_bonds_valences( int nPeriodicNum,
                        int bonds_valence,
                        int num_H,
                        VAL_AT *pVA );


/* local prototypes */
int AddExplicitDeletedH( inp_ATOM *at,
                         int jv,
                         int num_at,
                         int *iDeletedH,
                         int *iH,
                         int nNumDeletedH,
                         int bTwoStereo );
int bFindCumuleneChain( inp_ATOM *at,
                        AT_NUMB i1,
                        AT_NUMB i2,
                        AT_NUMB nCumulene[],
                        int nMaxLen );
int set_bond_type( inp_ATOM *at,
                   AT_NUMB i1,
                   AT_NUMB i2,
                   int bType );
int set_cumulene_0D_parity( inp_ATOM *at,
                            inp_ATOM_STEREO *st,
                            int num_at,
                            int idelH1,
                            int i1,
                            int i2,
                            int idelH2,
                            int parity,
                            int len );
int set_atom_0D_parity( inp_ATOM *at,
                        inp_ATOM_STEREO *st,
                        int num_at,
                        int num_deleted_H,
                        int i1,
                        int parity );
int GetTgroupInfoFromInChI( T_GROUP_INFO *ti,
                            inp_ATOM *at,
                            AT_NUMB *endpoint,
                            INChI *pInChI );
int FillOutpStructEndpointFromInChI( INChI *pInChI, AT_NUMB **pEndpoint );
int SetStereoBondTypeFor0DParity( inp_ATOM *at, int i1, int m1 );
int SetStereoBondTypesFrom0DStereo( StrFromINChI *pStruct, INChI *pInChI );
void CopyAt2St( inp_ATOM *at, inp_ATOM_STEREO *st, int num_atoms );
void CopySt2At( inp_ATOM *at, inp_ATOM_STEREO *st, int num_atoms );
int RestoreAtomConnectionsSetStereo( StrFromINChI *pStruct,
                                     int iComponent,
                                     int iAtNoOffset,
                                     INChI *pInChI,
                                     INChI *pInChIMobH );

int RestoreAtomMakeBNS( struct tagINCHI_CLOCK *ic,
                        struct tagCANON_GLOBALS *pCG,
                        ICHICONST INPUT_PARMS *ip,
                        STRUCT_DATA *sd,
                        StrFromINChI *pStruct,
                        int iComponent,
                        int iAtNoOffset,
                        INChI *pInChI[],
                        const char *szCurHdr,
                        long num_inp,
                        int bHasSomeFixedH );

int nAddSuperCGroups( ALL_TC_GROUPS *pTCGroups );

int AddCGroups2TCGBnStruct( BN_STRUCT *pBNS,
                            StrFromINChI *pStruct,
                            VAL_AT *pVA,
                            ALL_TC_GROUPS *pTCGroups,
                            int nMaxAddEdges );
int AddTGroups2TCGBnStruct( BN_STRUCT *pBNS,
                            StrFromINChI *pStruct,
                            VAL_AT *pVA,
                            ALL_TC_GROUPS *pTCGroups,
                            int nMaxAddEdges );
BN_STRUCT* AllocateAndInitTCGBnStruct( StrFromINChI *pStruct,
                                       VAL_AT *pVA,
                                       ALL_TC_GROUPS *pTCGroups,
                                       int nMaxAddAtoms,
                                       int nMaxAddEdges,
                                       int max_altp,
                                       int *pNum_changed_bonds );
int nCountBnsSizes( inp_ATOM *at,
                    int num_at,
                    int nAddEdges2eachAtom,
                    int nAddVertices,
                    T_GROUP_INFO *ti,
                    VAL_AT *pVA,
                    ICHICONST SRM *pSrm,
                    ALL_TC_GROUPS *pTCGroups );

int GetAtomRestoreInfo( struct tagCANON_GLOBALS *pCG,
                        inp_ATOM *atom, int iat, VAL_AT *pVArray,
                        ICHICONST SRM *pSrm, int bMobileH, AT_NUMB *endpoint );

int AddEdgeFlow( int edge_cap,
                 int edge_flow,
                 BNS_EDGE *e01,
                 BNS_VERTEX *pv0 /*src*/,
                 BNS_VERTEX *pv1/*dest*/,
                 int *tot_st_cap,
                 int *tot_st_flow );

void SetEdgeCapFlow( BNS_EDGE *e, int edge_cap, int edge_flow );

void AddStCapFlow( BNS_VERTEX *vert_ficpoint,
                   int *tot_st_flow,
                   int *tot_st_cap,
                   int cap,
                   int flow );
void SetStCapFlow( BNS_VERTEX *vert_ficpoint,
                   int *tot_st_flow,
                   int *tot_st_cap,
                   int cap,
                   int flow );

int ConnectSuperCGroup( int nTCG_Plus,
                        int nAddGroups[],
                        int num_add,
                        int *pcur_num_vertices,
                        int *pcur_num_edges,
                        int *tot_st_cap,
                        int *tot_st_flow,
                        BN_STRUCT *pBNS,
                        ALL_TC_GROUPS *pTCGroups );
int ConnectTwoVertices( BNS_VERTEX *p1,
                        BNS_VERTEX *p2,
                        BNS_EDGE *e,
                        BN_STRUCT *pBNS,
                        int bClearEdge );

int nTautEndpointEdgeCap( inp_ATOM *at, VAL_AT *pVA, int i );

/*
int GetAtomBondFlow( inp_ATOM *atom, VAL_AT *pVA, ICHICONST SRM *pSrm, int iat, int ineigh );
void GetAtomStCapFlow( inp_ATOM *atom, VAL_AT *pVA, ICHICONST SRM *pSrm, int iat, int *pCap, int *pFlow );
int GetAtomToMCGroupInitEdgeCapFlow( EdgeFlow *nEdgeCap, EdgeFlow *nEdgeFlow, ICHICONST SRM *pSrm, inp_ATOM *at, VAL_AT *pVA, int iat );
void GetAtomToMetalInitEdgeCapFlow( EdgeFlow *nEdgeCap, EdgeFlow *nEdgeFlow );
*/

int AtomStcapStflow( inp_ATOM *atom,
                     VAL_AT *pVA,
                     ICHICONST SRM *pSrm,
                     int iat,
                     int *pnStcap,
                     int *pnStflow,
                     EdgeFlow *pnMGroupEdgeCap,
                     EdgeFlow *pnMGroupEdgeFlow );

int BondFlowMaxcapMinorder( inp_ATOM *atom,
                            VAL_AT *pVA,
                            ICHICONST SRM *pSrm,
                            int iat,
                            int ineigh,
                            int *pnMaxcap,
                            int *pnMinorder,
                            int *pbNeedsFlower );


int clean_charge_val( struct tagCANON_GLOBALS *pCG,
                      CHARGE_VAL *pChargeVal, int len,
                      inp_ATOM *atom, VAL_AT *pVA, int iat,
                      int bIsMetal, int bMobileH, AT_NUMB *endpoint );

int ReallocTCGroups( ALL_TC_GROUPS *pTCGroups, int nAdd );
int RegisterTCGroup( ALL_TC_GROUPS *pTCGroups, int nGroupType, int nGroupOrdNum,
                     int nVertexCap, int nVertexFlow, int nEdgeCap, int nEdgeFlow, int nNumEdges );
int CopyBnsToAtom( StrFromINChI *pStruct, BN_STRUCT *pBNS, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                   int bAllowZeroBondOrder );
int CheckBnsConsistency( StrFromINChI *pStruct, BN_STRUCT *pBNS, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int bNoRad );

int RunBnsRestore1( struct tagCANON_GLOBALS *pCG,
                    struct tagINCHI_CLOCK *ic, ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS,
                    BN_DATA *pBD, StrFromINChI *pStruct,
                    VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, INChI *pInChI[], long num_inp, int bHasSomeFixedH );

int RunBnsRestoreOnce( BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups );
int nNumEdgesToCnVertex( MY_CONST C_NODE *pCN, int len, int v );
int ConnectMetalFlower( int *pcur_num_vertices, int *pcur_num_edges,
                        int *tot_st_cap, int *tot_st_flow, ICHICONST SRM *pSrm,
                        BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups );
int AddRadicalToMetal( int *tot_st_cap, int *tot_st_flow, ICHICONST SRM *pSrm, BN_STRUCT *pBNS,
                       ALL_TC_GROUPS *pTCGroups );
int bMayBeACationInMobileHLayer( inp_ATOM *at, VAL_AT *pVA, int iat, int bMobileH );

int MakeOneInChIOutOfStrFromINChI( struct tagCANON_GLOBALS *pCG, struct tagINCHI_CLOCK *ic, ICHICONST INPUT_PARMS *ip,
                                   STRUCT_DATA *sd, StrFromINChI *pStruct, inp_ATOM *at2, inp_ATOM *at3,
                                   ALL_TC_GROUPS *pTCGroups );
void IncrZeroBondsAndClearEndpts( inp_ATOM *at, int num_at, int iComponent );
void IncrZeroBonds( inp_ATOM *at, int num_at, int iComponent );
void ClearEndpts( inp_ATOM *at, int num_at );
int DisplayRestoredComponent( struct tagCANON_GLOBALS *pCG, StrFromINChI *pStruct, int iComponent, int iAtNoOffset,
                              INChI *pInChI, const char *szCurHdr );
int cmp_charge_val( const void *a1, const void *a2, void * );

int EvaluateChargeChanges( BN_STRUCT *pBNS, VAL_AT *pVA, int *pnDeltaH, int *pnDeltaCharge, int *pnNumVisitedAtoms );
int RunBnsTestOnce( BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA, Vertex *pvFirst, Vertex *pvLast,
                    int *pPathLen, int *pnDeltaH, int *pnDeltaCharge, int *pnNumVisitedAtoms );
int comp_cc_cand( const void *a1, const void *a2 );
int MoveRadToAtomsAddCharges( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                              inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int forbidden_mask );
void RemoveForbiddenBondFlowBits( BN_STRUCT *pBNS, int forbidden_edge_mask_int );

int PlusFromDB_N_DB_O_to_Metal( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );

int AdjustTgroupsToForbiddenEdges2( BN_STRUCT *pBNS, inp_ATOM *at, VAL_AT *pVA, int num_atoms, int forbidden_mask );
int AllocEdgeList( EDGE_LIST *pEdges, int nLen );
int AddToEdgeList( EDGE_LIST *pEdges, int iedge, int nAddLen );
int RemoveFromEdgeListByIndex( EDGE_LIST *pEdges, int index );
int RemoveFromEdgeListByValue( EDGE_LIST *pEdges, int iedge );
int FindInEdgeList( EDGE_LIST *pEdges, int iedge );
int AllocBfsQueue( BFS_Q *pQ, int num_at, int min_ring_size );
void RemoveForbiddenEdgeMask( BN_STRUCT *pBNS, EDGE_LIST *pEdges, int forbidden_edge_mask );
void SetForbiddenEdgeMask( BN_STRUCT *pBNS, EDGE_LIST *pEdges, int forbidden_edge_mask );
int ForbidCarbonChargeEdges( BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups, EDGE_LIST *pCarbonChargeEdges,
                             int forbidden_edge_mask );
int ForbidMetalCarbonEdges( BN_STRUCT *pBNS, inp_ATOM *at, int num_at, VAL_AT *pVA,
                            ALL_TC_GROUPS *pTCGroups, EDGE_LIST *pMetalCarbonEdges, int forbidden_edge_mask );
int ForbidNintrogenPlus2BondsInSmallRings( BN_STRUCT *pBNS, inp_ATOM *at, int num_at,
                                           VAL_AT *pVA, int min_ring_size, ALL_TC_GROUPS *pTCGroups,
                                           EDGE_LIST *pNplus2BondsEdges, int forbidden_edge_mask );
int RearrangePlusMinusEdgesFlow( BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA,
                                 ALL_TC_GROUPS *pTCGroups, int forbidden_edge_mask );
int IncrementZeroOrderBondsToHeteroat( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                       inp_ATOM *at, inp_ATOM *at2,
                                       VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                       int *pnNumRunBNS, int *pnTotalDelta,
                                       int forbidden_edge_mask );
int MoveChargeFromHeteroatomsToMetals( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                       inp_ATOM *at, inp_ATOM *at2,
                                       VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                       int *pnNumRunBNS, int *pnTotalDelta,
                                       int forbidden_edge_mask );
int EliminateChargeSeparationOnHeteroatoms( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                            inp_ATOM *at, inp_ATOM *at2,
                                            VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                            int *pnNumRunBNS, int *pnTotalDelta,
                                            int forbidden_edge_mask, int forbidden_stereo_edge_mask );
int MovePlusFromS2DiaminoCarbon( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                 inp_ATOM *at, inp_ATOM *at2,
                                 VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                 int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int RestoreCyanoGroup( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                       inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                       int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int RestoreIsoCyanoGroup( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                          inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                          int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int RestoreNNNgroup( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                     inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                     int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int FixMetal_Nminus_Ominus( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                            inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                            int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int EliminateNitrogen5Val3Bonds( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                 inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                 int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int Convert_SIV_to_SVI( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                        inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                        int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int MoveMobileHToAvoidFixedBonds( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                  inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                  int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int RemoveRadFromMobileHEndpoint( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                  inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                  int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int RemoveRadFromMobileHEndpointFixH( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                      inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                      int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int MoveChargeToMakeCenerpoints( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                 inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                 int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int CheckAndRefixStereobonds( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                              inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int MoveChargeToRemoveCenerpoints( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                   inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                   int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int MakeSingleBondsMetal2ChargedHeteroat( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                                          inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                          int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int SaltBondsToCoordBonds( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct,
                           inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                           int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int FixLessHydrogenInFormula( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
                              inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int FixMoreHydrogenInFormula( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
                              inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                              int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int FixAddProtonForADP( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
                        inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, ICR *picr,
                        int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int ConnectDisconnectedH( inp_ATOM *at, int num_atoms, int num_deleted_H );
int DisconnectedConnectedH( inp_ATOM *at, int num_atoms, int num_deleted_H );

int MakeInChIOutOfStrFromINChI2( struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, ICHICONST INPUT_PARMS *ip_inp,
                                 STRUCT_DATA *sd_inp, StrFromINChI *pStruct,
                                 int iComponent, int iAtNoOffset, long num_inp );

int GetChargeFlowerUpperEdge( BN_STRUCT *pBNS, VAL_AT *pVA, int nChargeEdge );
int get_pVA_atom_type( VAL_AT *pVA, inp_ATOM *at, int iat, int bond_type );

int NormalizeAndCompare( struct tagCANON_GLOBALS *pCG, struct tagINCHI_CLOCK *ic, ICHICONST INPUT_PARMS *ip,
                         STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
                         StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA,
                         ALL_TC_GROUPS *pTCGroups, INChI *pInChI[], long num_inp, int bHasSomeFixedH,
                         int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask, int forbidden_stereo_edge_mask );

/* call InChI normalization only */
int NormalizeStructure( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS,
                        StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2,
                        VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                        inp_ATOM *at_norm, inp_ATOM *at_fixed_bonds_out, T_GROUP_INFO *t_group_info );
/* create one InChI */
int MakeOneInChIOutOfStrFromINChI2( struct tagCANON_GLOBALS *pCG, struct tagINCHI_CLOCK *ic, ICHICONST INPUT_PARMS *ip,
                                    STRUCT_DATA *sd, BN_STRUCT *pBNS, StrFromINChI *pStruct,
                                    inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                    T_GROUP_INFO **t_group_info, inp_ATOM **at_norm, inp_ATOM **at_prep );

/* fixed-H */
int FillOutExtraFixedHDataRestr( StrFromINChI *pStruct );
int FillOutExtraFixedHDataInChI( StrFromINChI *pStruct, INChI *pInChI[] );

int FixFixedHRestoredStructure( struct tagCANON_GLOBALS *pCG,
                                struct tagINCHI_CLOCK *ic,
                                ICHICONST INPUT_PARMS *ip,
                                STRUCT_DATA *sd,
                                BN_STRUCT *pBNS, BN_DATA *pBD,
                                StrFromINChI *pStruct,
                                inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3,
                                VAL_AT *pVA,
                                ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ti,
                                inp_ATOM **at_norm, inp_ATOM **at_prep,
                                INChI *pInChI[], long num_inp,
                                int bHasSomeFixedH, int *pnNumRunBNS,
                                int *pnTotalDelta,
                                int forbidden_edge_mask,
                                int forbidden_stereo_edge_mask );

int FixRemoveExtraTautEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at,
                                 inp_ATOM *at2, inp_ATOM *atf, inp_ATOM *atn, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups,
                                 ICR *picr,
                                 int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask );
int FillOutCMP2FHINCHI( StrFromINChI *pStruct, inp_ATOM *at2, VAL_AT *pVA, INChI *pInChI[], CMP2FHINCHI *pc2i );
int FillOutCMP2MHINCHI( StrFromINChI *pStruct, ALL_TC_GROUPS *pTCGroups, inp_ATOM *at2,
                        VAL_AT *pVA, INChI *pInChI[], CMP2MHINCHI *pc2i );

int bHas_N_V( inp_ATOM *at2, int num_atoms );

int GetPlusMinusVertex( BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups, int bCheckForbiddenPlus, int bCheckForbiddenMinus );

int FixMobileHRestoredStructure( struct tagCANON_GLOBALS *pCG,
                                 struct tagINCHI_CLOCK *ic,
                                 ICHICONST INPUT_PARMS *ip,
                                 STRUCT_DATA *sd,
                                 BN_STRUCT *pBNS,
                                 BN_DATA *pBD,
                                 StrFromINChI *pStruct,
                                 inp_ATOM *at,
                                 inp_ATOM *at2,
                                 inp_ATOM *at3,
                                 VAL_AT *pVA,
                                 ALL_TC_GROUPS *pTCGroups,
                                 T_GROUP_INFO **ppt_group_info,
                                 inp_ATOM **ppat_norm,
                                 inp_ATOM **ppat_prep,
                                 INChI *pInChI[],
                                 long num_inp,
                                 int bHasSomeFixedH,
                                 int *pnNumRunBNS,
                                 int *pnTotalDelta,
                                 int forbidden_edge_mask,
                                 int forbidden_stereo_edge_mask );

int FixRestoredStructureStereo( struct tagCANON_GLOBALS *pCG, struct tagINCHI_CLOCK *ic, INCHI_MODE cmpInChI, ICR *icr,
                                INCHI_MODE cmpInChI2, ICR *icr2,
                                ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD,
                                StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA,
                                ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ppt_group_info, inp_ATOM **ppat_norm,
                                inp_ATOM **ppat_prep, INChI *pInChI[], long num_inp,
                                int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask,
                                int forbidden_stereo_edge_mask );

int AddRemProtonsInRestrStruct( struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, ICHICONST INPUT_PARMS *ip,
                                STRUCT_DATA *sd, long num_inp,
                                int bHasSomeFixedH,
                                StrFromINChI *pStruct, int num_components,
                                StrFromINChI *pStructR, int num_componentsR,
                                NUM_H *pProtonBalance, int *recmet_change_balance );
int AllInchiToStructure( struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, ICHICONST INPUT_PARMS *ip,
                         STRUCT_DATA *sd, long num_inp, char *szCurHdr,
                         ICHICONST SRM *pSrm, int bReqNonTaut, StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM],
                         InpInChI *pOneInput );
int AddProtonAndIsoHBalanceToMobHStruct( struct tagINCHI_CLOCK *ic,
                                         struct tagCANON_GLOBALS *pCG,
                                         ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd,
                                         long num_inp, int bHasSomeFixedH,
                                         char *szCurHdr,
                                         StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM],
                                         InpInChI *pOneInput );
int InChI2Atom( struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG,
                ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd,
                const char *szCurHdr, long num_inp,
                StrFromINChI *pStruct, int iComponent,
                int iAtNoOffset, int bI2A_Flag, int bHasSomeFixedH,
                InpInChI *pOneInput );


int MarkDisconectedIdenticalToReconnected( InpInChI *pOneInput );
void RemoveFixHInChIIdentical2MobH( InpInChI *pOneInput );
void SetUpSrm( SRM *pSrm );
void FreeInputInChI2Struct( InpInChI *pOneInput );
void FreeStrFromINChI( StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM], int nNumComponents[INCHI_NUM][TAUT_NUM] );
int OldPrintCompareOneOrigInchiToRevInChI( StrFromINChI *pStruct, INChI *pInChI[TAUT_NUM], int bMobileH,
                                           int iComponent, long num_inp, char *szCurHdr );
int CompareOneOrigInchiToRevInChI( StrFromINChI *pStruct, INChI *pInChI[TAUT_NUM], int bMobileH, int iComponent,
                                   long num_inp, char *szCurHdr,
                                   COMPONENT_REM_PROTONS *nCurRemovedProtons, INCHI_MODE CompareInchiFlags[] );
int CompareTwoPairsOfInChI( INChI *pInChI1[TAUT_NUM], INChI *pInChI2[TAUT_NUM],
                            int bMobileH, INCHI_MODE CompareInchiFlags[] );
INCHI_MODE CompareReversedINChI3( INChI *i1 /* InChI from reversed struct */, INChI *i2 /* input InChI */,
                                  INChI_Aux *a1, INChI_Aux *a2, int *err );
INCHI_MODE CompareReversedStereoINChI3( INChI_Stereo *s1/* InChI from reversed struct */,
                                        INChI_Stereo *s2 /* input InChI */, ICR *picr );
int CompareAllOrigInchiToRevInChI( StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM], InpInChI *pOneInput, int bReqNonTaut,
                                   long num_inp, char *szCurHdr );
int CompareAllDisconnectedOrigInchiToRevInChI( StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM],
                                               InpInChI *pOneInput, int bHasSomeFixedH,
                                               long num_inp, char *szCurHdr );
int insertions_sort_AT_NUMB( AT_NUMB *base, int num );

int AddRemIsoProtonsInRestrStruct( struct tagINCHI_CLOCK *ic,
                                   struct tagCANON_GLOBALS *pCG,
                                   ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd,
                                   long num_inp, int bHasSomeFixedH,
                                   StrFromINChI *pStruct, int num_components,
                                   StrFromINChI *pStructR, int num_componentsR,
                                   NUM_H pProtonBalance[],
                                   NUM_H recmet_change_balance[] );

int OutputInChIOutOfStrFromINChI( struct tagINCHI_CLOCK *ic,
                                  struct tagCANON_GLOBALS *pCG,
                                  ICHICONST INPUT_PARMS *ip_inp,
                                  STRUCT_DATA *sd_inp,
                                  long num_inp, int bINChIOutputOptions,
                                  INCHI_IOSTREAM *pout, INCHI_IOSTREAM *plog,
                                  InpInChI *pOneInput, int bHasSomeFixedH,
                                  unsigned char save_opt_bits );

int MergeStructureComponents( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, long num_inp, char *szCurHdr,
                              ICHICONST SRM *pSrm, int bReqNonTaut, StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM],
                              InpInChI *pOneInput );
int AddOneMsg( char *szMsg, int used_len, int tot_len, const char *szAddMsg, const char *szDelim );
int FillOutCompareMessage( char *szMsg, int nLenMsg, INCHI_MODE bits[] );
void clear_t_group_info( T_GROUP_INFO *ti );
int bInpInchiComponentExists( InpInChI *pOneInput, int iINCHI, int bMobileH, int k );
int bInpInchiComponentDeleted( InpInChI *pOneInput, int iInChI, int bMobileH, int k );
int bRevInchiComponentExists( StrFromINChI *pStruct, int iInChI, int bMobileH, int k );
int bRevInchiComponentDeleted( StrFromINChI *pStruct, int iInChI, int bMobileH, int k );
int DetectInpInchiCreationOptions( InpInChI *pOneInput, int *bHasReconnected, int *bHasMetal,
                                   int *bHasFixedH, int *sFlag, int *bTautFlag );

int DisplayStructureComponents( struct tagCANON_GLOBALS *pCG,
                                ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd,
                                long num_inp, char *szCurHdr,
                                ICHICONST SRM *pSrm,
                                int bReqNonTaut,
                                StrFromINChI *pStruct[INCHI_NUM][TAUT_NUM],
                                InpInChI *pOneInput );
int DisplayOneRestoredComponent( struct tagCANON_GLOBALS *pCG,
                                 StrFromINChI *pStruct, inp_ATOM *at,
                                 int iComponent, int nNumComponents, int bMobileH,
                                 const char *szCurHdr );
int DisplayAllRestoredComponents( struct tagCANON_GLOBALS *pCG,
                                  inp_ATOM *at, int num_at, const char *szCurHdr );

int CountStereoTypes( INChI *pInChI, int *num_known_SB, int *num_known_SC,
                      int *num_unk_und_SB, int *num_unk_und_SC,
                      int *num_SC_PIII, int *num_SC_AsIII );
int GetNumNeighborsFromInchi( INChI *pInChI, AT_NUMB nAtNumber );
int bIsUnsatCarbonInASmallRing( inp_ATOM *at, VAL_AT *pVA, int iat, BFS_Q *pbfsq, int min_ring_size );

int MakeProtonComponent( StrFromINChI *pStruct, int iComponent, int num_prot );

void FreeInpInChI( InpInChI *pOneInput );


/* extra configurarion */
#define KEEP_METAL_EDGE_FLOW                 0  /* counterexample: mdb0-1738.sdf.txt */
#define MOVE_CHARGES_FROM_HETEREO_TO_METAL   0  /* disabled */
#define FIX_ADD_PROTON_FOR_ADP               0  /* not used */

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif    /* _ICHIRVRS_H_ */
