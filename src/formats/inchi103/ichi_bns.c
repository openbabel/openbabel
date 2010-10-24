/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.03
 * May 9, 2010
 *
 * Originally developed at NIST
 * Modifications and additions by IUPAC and the InChI Trust
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-2.1.php
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "mode.h"

#include "ichierr.h"
#include "incomdef.h" 
#include "inpdef.h"
#include "extr_ct.h"
#include "ichitaut.h"
#include "ichinorm.h"
#include "util.h"

#include "ichicomp.h"
#include "ichister.h"

#include "ichi_bns.h"


#define BNS_MARK_ONLY_BLOCKS        1  /* 1 => find only blocks, do not search for ring systems */
#define ALLOW_ONLY_SIMPLE_ALT_PATH  0  /* 0 => allow alt. path to contain same bond 2 times (in opposite directions) */

#define CHECK_TG_ALT_PATH           0  /* 1=> when chacking alt path of a tautomeric atom modify
                                              t-group, not the atom */
                                       /* 0=> old mode */

#define FIX_CPOINT_BOND_CAP         1  /* 1=> fix bug in case of double bond from neutral cpoint */

#define RESET_EDGE_FORBIDDEN_MASK   1  /* 1: previous; 0: do not apply "edge->forbidden &= pBNS->edge_forbidden_mask" */
#if ( RESET_EDGE_FORBIDDEN_MASK == 1 )
#define IS_FORBIDDEN(EDGE_FORBIDDEN, PBNS)     (EDGE_FORBIDDEN)
#else
#define IS_FORBIDDEN(EDGE_FORBIDDEN, PBNS)     (EDGE_FORBIDDEN & PBNS->edge_forbidden_mask)
#endif
                                       

typedef enum tagAtTypeTotals {
    /* counts do not include:
       charged atom adjacent to another charged atom
       atom in an unusual valence state or adjacent to an atom in an unusual valence state
       radicals different from singlet
    */
    /*ATTOT_NUM_Plus. */    /* number of positive charges, +1, is (ATTOT_NUM_CHARGES + ATTOT_TOT_CHARGE)/2 */
    /*ATTOT_NUM_Minus.*/    /* number of negative charges, -1, is (ATTOT_NUM_CHARGES - ATTOT_TOT_CHARGE)/2 */
    ATTOT_NUM_NP_Plus,    /*  0 no H: =N(+)=, #N(+)-, =N(+)<, does not include onium cations >P(+)<, >N(+)< */
    ATTOT_NUM_NP_Proton,  /*  1 H(+): -NH3(+), =NH2(+), >NH2(+), =NH(+)-, >NH(+)-, #NH(+), N=N,P */
    ATTOT_NUM_NP_H,       /*  2 H:    -NH2, =NH, >NH -NH(-) */
    ATTOT_NUM_N_Minus,    /*  3 (-):  -NH(-), >N(-), =N(-) */
    ATTOT_NUM_NP,         /*  4 no H: >N- =N-, #N  */
    ATTOT_NUM_ON,         /*  5 -N=O: do not allow -N=O => -NH-OH during H(+) add/removal */
    ATTOT_NUM_COH,        /*  6 =C-OH, #C-OH; O=O,S,Se,Te */
    ATTOT_NUM_CSH,        /*  7 -C-SH, -C-SeH -C-TeH  */
    ATTOT_NUM_ZOH,        /*  8 =Z-OH, #Z-OH; O=O,S,Se,Te; Z may have charge, Z != C */
    ATTOT_NUM_OOH,        /*  9  -O-OH, O=O,S,Se,Te */
    ATTOT_NUM_ZOOH,       /* 10  O=Z-OH, O=O,S,Se,Te */
    ATTOT_NUM_NOH,        /* 11  =N-OH, -N(-)-OH */
    ATTOT_NUM_N_OH,       /* 12 >N-OH, -NH-OH, >NH(+)-OH, -N(-)-OH */
    ATTOT_NUM_CO,         /* 13 -C=O, =C=O; O=O,S,Se,Te  */
    ATTOT_NUM_ZO,         /* 14 -Z=O, =Z=O; O=O,S,Se,Te; Z may have charge  */
    ATTOT_NUM_NO,         /* 15 -N=O, =N(+)=O */
    ATTOT_NUM_N_O,        /* 16  >N(+)=O, =N(+)=O */
    ATTOT_NUM_CO_Minus,   /* 17 =C-O(-), #C-O(-); O=O,S,Se,Te  */
    ATTOT_NUM_CS_Minus,   /* 18 -C-S(-); S = S, Se, Te */
    ATTOT_NUM_ZO_Minus,   /* 19 =Z-O(-), #Z-O(-); O = O, S, Se, Te */
    ATTOT_NUM_OO_Minus,   /* 20 -O-O(-), O=O,S,Se,Te */
    ATTOT_NUM_ZOO_Minus,  /* 21 O=Z-O(-), O=O,S,Se,Te */
    ATTOT_NUM_NO_Minus,   /* 22 >N-O(-), -NH-O(-) */
    ATTOT_NUM_N_O_Minus,  /* 23 -NH-O(-), >N-O(-); O = O, S, Se, Te */
    ATTOT_NUM_O_Minus,    /* 24 -Z-O(-); O=O,S,Se,Te */
    ATTOT_NUM_OH_Plus,    /* 25 any OH(+) */
    ATTOT_NUM_O_Plus,     /* 26 any O(+) without H */
    ATTOT_NUM_Proton,     /* 27 proton */
    ATTOT_NUM_HalAnion,   /* 28 Halogen anion */
    ATTOT_NUM_HalAcid,    /* 29 Halogen acid */
    ATTOT_NUM_Errors,     /* 30 for debugging */
    ATTOT_TOT_CHARGE,     /* 31 total of positive and negative single charges, +1 and -1 */
    ATTOT_NUM_CHARGES,    /* 32 number of positive and negative single charges, +1 and -1 */
    ATTOT_ARRAY_LEN       /* 33 array length */
} AT_TYPE_TOTALS;
    
#define ATBIT_NP_Plus    (1 << ATTOT_NUM_NP_Plus)  
#define ATBIT_NP_Proton  (1 << ATTOT_NUM_NP_Proton)
#define ATBIT_NP_H       (1 << ATTOT_NUM_NP_H)     
#define ATBIT_N_Minus    (1 << ATTOT_NUM_N_Minus)  
#define ATBIT_NP         (1 << ATTOT_NUM_NP)
#define ATBIT_ON         (1 << ATTOT_NUM_ON)
#define ATBIT_COH        (1 << ATTOT_NUM_COH)      
#define ATBIT_CSH        (1 << ATTOT_NUM_CSH)      
#define ATBIT_ZOH        (1 << ATTOT_NUM_ZOH)      
#define ATBIT_OOH        (1 << ATTOT_NUM_OOH)
#define ATBIT_ZOOH       (1 << ATTOT_NUM_ZOOH)
#define ATBIT_NOH        (1 << ATTOT_NUM_NOH)
#define ATBIT_N_OH       (1 << ATTOT_NUM_N_OH)
#define ATBIT_CO         (1 << ATTOT_NUM_CO)       
#define ATBIT_ZO         (1 << ATTOT_NUM_ZO)
#define ATBIT_NO         (1 << ATTOT_NUM_NO)
#define ATBIT_N_O        (1 << ATTOT_NUM_N_O)
#define ATBIT_CO_Minus   (1 << ATTOT_NUM_CO_Minus) 
#define ATBIT_CS_Minus   (1 << ATTOT_NUM_CS_Minus) 
#define ATBIT_ZO_Minus   (1 << ATTOT_NUM_ZO_Minus)
#define ATBIT_OO_Minus   (1 << ATTOT_NUM_OO_Minus)
#define ATBIT_ZOO_Minus  (1 << ATTOT_NUM_ZOO_Minus)
#define ATBIT_NO_Minus   (1 << ATTOT_NUM_NO_Minus)
#define ATBIT_N_O_Minus  (1 << ATTOT_NUM_N_O_Minus)
#define ATBIT_O_Minus    (1 << ATTOT_NUM_O_Minus)
#define ATBIT_OH_Plus    (1 << ATTOT_NUM_OH_Plus)
#define ATBIT_O_Plus     (1 << ATTOT_NUM_O_Plus) 
#define ATBIT_Proton     (1 << ATTOT_NUM_Proton)
#define ATBIT_HalAnion   (1 << ATTOT_NUM_HalAnion) 
#define ATBIT_HalAcid    (1 << ATTOT_NUM_HalAcid) 


#define ATBIT_Errors     (1 << ATTOT_NUM_Errors)   

typedef struct tagProtonRemovalMaskAndType {
    int typePos;  /* atoms accessible to positive charges */
    int maskPos; 
    int typeNeg; /* atoms accessible to negative charges */
    int maskNeg; 
    int typeH;     /* atoms accessible to hydrogen atoms */
    int maskH; 
} PRMAT;

#define PR_SIMPLE_MSK  (ATBIT_NP_Proton | ATBIT_OH_Plus)
#define PR_SIMPLE_TYP  (ATT_ATOM_N | ATT_ATOM_P | ATT_O_PLUS)

#define ATBIT_MSK_NP   (ATBIT_NP_Plus | ATBIT_NP_Proton | ATBIT_NP_H | ATBIT_N_Minus | ATBIT_NP)

#define KNOWN_ACIDIC_TYPE   (ATT_ACIDIC_CO | ATT_ACIDIC_S | ATT_OO | ATT_ZOO | ATT_NO)

#define ATBIT_MSK_OS   (ATBIT_COH | ATBIT_CSH | ATBIT_ZOH | ATBIT_OOH | ATBIT_ZOOH | ATBIT_NOH | ATBIT_N_OH |\
                        ATBIT_CO | ATBIT_ZO  | ATBIT_NO | ATBIT_N_O |\
                        ATBIT_CO_Minus | ATBIT_CS_Minus | ATBIT_ZO_Minus | ATBIT_OO_Minus |\
                        ATBIT_ZOO_Minus | ATBIT_NO_Minus | ATBIT_N_O_Minus /*| ATBIT_O_Minus*/ )
#define ATBIT_MSK_H    (ATBIT_NP_Proton | ATBIT_NP_H | ATBIT_COH | ATBIT_CSH | ATBIT_ZOH | ATBIT_OOH |\
                        ATBIT_ZOOH | ATBIT_NOH | ATBIT_N_OH)

#define ATTYP_OS    (ATT_ACIDIC_CO | ATT_ACIDIC_S | ATT_OO | ATT_ZOO | ATT_NO /*| ATT_OTHER_NEG_O*/ | ATT_OTHER_ZO)
#define ATTYP_NP    (ATT_ATOM_N | ATT_ATOM_P)
#define ATTYP_N     (ATT_ATOM_N)
#define ATTYP_P     (ATT_ATOM_P)

/************* simple proton removal from acids **************************/
#define AR_ANY_OH       0  /* 1 => create unknown to be acidic anions */
#define AR_SIMPLE_STEPS 3
/* acidic groups for proton removal, step 1 */
#define AR_SIMPLE_MSK1  (ATBIT_COH | ATBIT_CSH | ATBIT_OOH | ATBIT_ZOOH | ATBIT_NOH | ATBIT_HalAcid)
#define AR_SIMPLE_TYP1  (ATT_ACIDIC_CO | ATT_ACIDIC_S | ATT_OO | ATT_ZOO | ATT_NO | ATT_HalAcid)
/* acidic groups for proton removal, step 2 */
#define AR_SIMPLE_MSK2  (AR_ANY_OH? (ATBIT_N_OH)  :0)
#define AR_SIMPLE_TYP2  (AR_ANY_OH? (ATT_N_O)     :0)
/* acidic groups for proton removal, step 3 */
#define AR_SIMPLE_MSK3  (AR_ANY_OH? (ATBIT_ZOH)   :0)
#define AR_SIMPLE_TYP3  (AR_ANY_OH? (ATT_OTHER_ZO):0)

/************* simple proton addition to acids **************************/
#define AA_ANY_O_Minus  0  /* 1 => neutralize unknown to be acidic anions */
#define AA_SIMPLE_STEPS 3
/* acidic groups for proton addition, step 1 */
#define AA_SIMPLE_MSK1  (ATBIT_CO_Minus | ATBIT_CS_Minus | ATBIT_OO_Minus | ATBIT_ZOO_Minus | ATBIT_NO_Minus | ATBIT_O_Minus | ATBIT_HalAnion)
#define AA_SIMPLE_TYP1  (ATT_ACIDIC_CO | ATT_ACIDIC_S | ATT_OO | ATT_ZOO | ATT_NO | ATT_OH_MINUS | ATT_HalAnion )
/* acidic groups for proton addition, step 2 */
#define AA_SIMPLE_MSK2  (AA_ANY_O_Minus? (ATBIT_N_O_Minus)               :0)
#define AA_SIMPLE_TYP2  (AA_ANY_O_Minus? (ATT_N_O)                       :0)
/* acidic groups for proton addition, step 3 */
#define AA_SIMPLE_MSK3  (AA_ANY_O_Minus? (ATBIT_ZO_Minus | ATBIT_O_Minus):0)
#define AA_SIMPLE_TYP3  (AA_ANY_O_Minus? (ATT_OTHER_ZO)                  :0)

#if ( FIX_NP_MINUS_BUG == 1 )
/* allow to add H(+) to =N(-) which previously was #N */
#undef AA_SIMPLE_STEPS
#define AA_SIMPLE_STEPS 4
#define AA_SIMPLE_MSK4  ATBIT_N_Minus
#define AA_SIMPLE_TYP4  ATT_NP_MINUS_V23
#endif
/************* hard proton removal from NP **************************/
/* (+) charge group for proton removal: mask & type */
#define PR_HARD_MSK_POS   ATBIT_MSK_NP
#define PR_HARD_TYP_POS   ATTYP_N
#define PR_HARD_TYP_POSP  ATTYP_P
/* (-) charge group for proton removal */
#define PR_HARD_MSK_NEG   (ATBIT_MSK_NP | ATBIT_MSK_OS)
#define PR_HARD_TYP_NEG   (ATTYP_N | ATTYP_OS)
/* H-group for proton removal */
#define PR_HARD_MSK_H     (ATBIT_MSK_NP | ATBIT_MSK_OS)
#define PR_HARD_TYP_H     (ATTYP_N | ATTYP_OS)

/************* hard proton removal from acids **************************/
/* (+) charge group for proton removal: mask & type */
#define AR_HARD_MSK_POS   ATBIT_MSK_NP
#define AR_HARD_TYP_POS   ATTYP_N
/* (-) charge group for proton removal */
#define AR_HARD_MSK_NEG   (ATBIT_MSK_NP | ATBIT_MSK_OS)
#define AR_HARD_TYP_NEG   (ATTYP_N | ATTYP_OS)
/* H-group acid for proton removal */
#define AR_HARD_MSK_HA    (ATBIT_CO | ATBIT_NO )
#define AR_HARD_TYP_HA    (ATT_ACIDIC_CO | ATT_NO)
/* H-group non-acid for proton removal */
#define AR_HARD_MSK_HN    ((ATBIT_MSK_NP | ATBIT_MSK_OS) & ~AR_HARD_MSK_HA)
#define AR_HARD_TYP_HN    ((ATTYP_N | ATTYP_OS) /*& ~AR_HARD_TYP_HA*/)

/************* hard proton addition to acids **************************/
/* (+) charge group for proton removal: mask & type */
#define AA_HARD_MSK_POS   ATBIT_MSK_NP
#define AA_HARD_TYP_POS   ATTYP_N
/* (-) charge group for negative charge removal */
#define AA_HARD_MSK_NEG   ((ATBIT_MSK_NP | ATBIT_MSK_OS) & ~(ATBIT_CO | ATBIT_NO ))
#define AA_HARD_TYP_NEG   (ATTYP_N | ATTYP_OS)
/* (-) charge group to accept negative charges  */
#define AA_HARD_MSK_CO    (ATBIT_CO | ATBIT_NO )
#define AA_HARD_TYP_CO    (ATT_ACIDIC_CO | ATT_NO)
/* H-group non-acid for proton removal */
#define AA_HARD_MSK_H     (ATBIT_MSK_NP | ATBIT_MSK_OS)
#define AA_HARD_TYP_H     (ATTYP_N | ATTYP_OS)


/*********************************************************************************/

#define BNS_MAX_NUM_FLOW_CHANGES (1+2*MAX_BOND_EDGE_CAP)

/* -- opiginal Pascal values --
#define NO_VERTEX     0
#define BLOSSOM_BASE -1
#define FIRST_INDX    1
*/


#define TREE_NOT_IN_M  0  /* not in T or T' */
#define TREE_IN_2      1  /* in T' and not s-reachable */
#define TREE_IN_2BLOSS 2  /* in T' and in a blossom, is s-reachable */
#define TREE_IN_1      3  /* in T and is s-reachable */

#define TREE_IS_S_REACHABLE(X) (Tree[X] >= TREE_IN_2BLOSS)
#define TREE_IS_ON_SCANQ TREE_IS_S_REACHABLE
/* #define TREE_IS_ON_SCANQ(X)    (Tree[X] != TREE_NOT_IN_M) */
#define TREE_MARK(X, MARK)           do{ if( Tree[X] < MARK ) Tree[X]=MARK; }while(0)

/**********************************************************************************
 *  store changes done to check whether an alternating path exists
 *  (see bSetBnsToCheckAltPath, bRestoreBnsAfterCheckAltPath)
 **********************************************************************************/
typedef struct tagAltPathChanges {
    /* caps changed in up to 2 vertices */
    VertexFlow nOldCapsVert[2][MAXVAL+1];
    Vertex     vOldVert[2];
    S_CHAR     bSetOldCapsVert[2]; /* number of caps to restore, including st-cap */
    /* save ids of the newly created temporary vertices */
    Vertex     vNewVertex[2];
    S_CHAR     bSetNew[2];         /* indicators whether to remove vertices */

} ALT_PATH_CHANGES;

/*************************************************************************/

/* local functions */

int RestoreRadicalsOnly( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at );
int bRadChangesAtomType( BN_STRUCT *pBNS, BN_DATA *pBD, Vertex v, Vertex v_1, Vertex v_2 );
int BnsAdjustFlowBondsRad( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms );
int SetAtomRadAndChemValFromVertexCapFlow( BN_STRUCT *pBNS, inp_ATOM *atom, int v1 );
int bNeedToTestTheFlow( int bond_type, int nTestFlow, int bTestForNonStereoBond );

int RestoreEdgeFlow( BNS_EDGE *edge, int delta, int bChangeFlow );
int SetAtomBondType( BNS_EDGE *edge, U_CHAR *bond_type12, U_CHAR *bond_type21, int delta, int bChangeFlow );
int RestoreBnStructFlow( BN_STRUCT *pBNS, int bChangeFlow );

int CompTGroupNumber( const void *tg1, const void *tg2 );
int CompCGroupNumber( const void *cg1, const void *cg2 );

/* Rings, Blocks, Non-stereo bonds */
int ReInitBnStructForAltBns( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bUnknAltAsNoStereo );
int MarkRingSystemsAltBns( BN_STRUCT* pBNS, int bUnknAltAsNoStereo );
int MarkNonStereoAltBns( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bUnknAltAsNoStereo );

/* called from BalancedNetworkSearch */

int GetVertexDegree( BN_STRUCT* pBNS, Vertex v );
/* Vertex Get2ndNeighbor1( BN_STRUCT* pBNS, Vertex u, EdgeIndex iedge ); not used */
Vertex Get2ndEdgeVertex( BN_STRUCT* pBNS, Edge uv );
Vertex GetVertexNeighbor( BN_STRUCT* pBNS, Vertex v, int neigh, EdgeIndex *iedge );
int GetEdgePointer( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv, BNS_EDGE **uv, S_CHAR *s_or_t );
int AugmentEdge( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv, int delta, S_CHAR bReverse, int bChangeFlow );
int rescap_mark( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv );
int rescap(  BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv );
Vertex FindBase( Vertex u, Vertex *BasePtr );
int FindPathToVertex_s( Vertex x, Edge *SwitchEdge, Vertex *BasePtr, Vertex *Path, int MaxPathLen );
Vertex MakeBlossom( BN_STRUCT* pBNS, Vertex *ScanQ, int *pQSize,
                    Vertex *Pu, Vertex *Pv, int max_len_Pu_Pv,
                    Edge *SwitchEdge, Vertex *BasePtr,
                    Vertex u, Vertex v, EdgeIndex iuv, Vertex b_u, Vertex b_v, S_CHAR *Tree );
int PullFlow( BN_STRUCT *pBNS, Edge *SwitchEdge, Vertex x, Vertex y, int delta, S_CHAR bReverse, int bChangeFlow );
int FindPathCap( BN_STRUCT* pBNS, Edge *SwitchEdge, Vertex x, Vertex y, int delta );


/*
int SetBondType( BNS_EDGE *edge, U_CHAR *bond_type12, U_CHAR *bond_type21, int delta, int bChangeFlow );
int SetBondsRestoreBnStructFlow( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bChangeFlow );
*/
int SetBondsFromBnStructFlow( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bChangeFlow0 );
int MarkAtomsAtTautGroups(  BN_STRUCT *pBNS, int num_atoms, BN_AATG *pAATG, int nEnd1, int nEnd2 );

int nMinFlow2Check( BN_STRUCT *pBNS, int iedge );
int nMaxFlow2Check( BN_STRUCT *pBNS, int iedge );
int nCurFlow2Check( BN_STRUCT *pBNS, int iedge );

/* bonds testing */
/*
int bRestoreFlowToCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd, int nTestFlow, inp_ATOM *at, int num_atoms, int bChangeFlow );
*/
int bSetFlowToCheckOneBond( BN_STRUCT *pBNS, int iedge, int flow, BNS_FLOW_CHANGES *fcd );
int bRestoreFlowAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd );
int bSetBondsAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd, int nTestFlow, inp_ATOM *at, int num_atoms, int bChangeFlow );
int BnsTestAndMarkAltBonds(  BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms, BNS_FLOW_CHANGES *fcd, int bChangeFlow, int nBondTypeToTest );
int bIsAltBond( int bond_type );
/* fix bonds */
int fix_special_bonds( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int edge_forbidden_mask );
int TempFix_NH_NH_Bonds( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms );
int CorrectFixing_NH_NH_Bonds( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms );

/* alt path testing */
int bSetBnsToCheckAltPath( BN_STRUCT *pBNS, int nVertDoubleBond, int nVertSingleBond, AT_NUMB type,
                           int path_type, ALT_PATH_CHANGES *apc, BNS_FLOW_CHANGES *fcd, int *nDots );
int bRestoreBnsAfterCheckAltPath( BN_STRUCT *pBNS, ALT_PATH_CHANGES *apc, int bChangeFlow );
Vertex GetGroupVertex(BN_STRUCT *pBNS, Vertex v1, AT_NUMB type);
BNS_IEDGE GetEdgeToGroupVertex( BN_STRUCT *pBNS, Vertex v1, AT_NUMB type);
int bAddNewVertex( BN_STRUCT *pBNS, int nVertDoubleBond, int nCap, int nFlow, int nMaxAdjEdges, int *nDots );
int AddNewEdge( BNS_VERTEX *p1, BNS_VERTEX *p2, BN_STRUCT *pBNS, int nEdgeCap, int nEdgeFlow );
int bAddStCapToAVertex( BN_STRUCT *pBNS, Vertex v1, Vertex v2, VertexFlow *nOldCapVertSingleBond, int *nDots, int bAdjacentDonors );

static void remove_alt_bond_marks(inp_ATOM *at, int num_atoms);
int bIsBnsEndpoint( BN_STRUCT *pBNS, int v );

/* protons removal, charge neutralization */
int is_acidic_CO( inp_ATOM *atom, int at_no );
int mark_at_type( inp_ATOM *atom, int num_atoms, int nAtTypeTotals[] );
int GetAtomChargeType( inp_ATOM *atom, int at_no, int nAtTypeTotals[], int *pMask, int bSubtract  );
int AddChangedAtHChargeBNS( inp_ATOM *at, int num_atoms, int nAtTypeTotals[], S_CHAR *mark );
int EliminatePlusMinusChargeAmbiguity( BN_STRUCT *pBNS, int num_atoms );
int AddOrRemoveExplOrImplH( int nDelta, inp_ATOM *at, int num_atoms, AT_NUMB at_no, T_GROUP_INFO *t_group_info );
int SubtractOrChangeAtHChargeBNS( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms,
                                 int nAtTypeTotals[], S_CHAR *mark, T_GROUP_INFO *t_group_info, int bSubtract );
int is_Z_atom( U_CHAR el_number );
int IsZOX( inp_ATOM *atom, int at_x, int ord );
int SimpleRemoveHplusNPO( inp_ATOM *at, int num_atoms, int nAtTypeTotals[], T_GROUP_INFO *t_group_info );
int CreateCGroupInBnStruct( inp_ATOM *at, int num_atoms,
                            BN_STRUCT *pBNS, int nType, int nMask, int nCharge );
int CreateTGroupInBnStruct( inp_ATOM *at, int num_atoms,
                            BN_STRUCT *pBNS, int nType, int nMask );
int RemoveLastGroupFromBnStruct( inp_ATOM *at, int num_atoms, int tg, BN_STRUCT *pBNS );
int SetInitCapFlowToCurrent( BN_STRUCT *pBNS );
int SimpleRemoveAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2remove );
int SimpleAddAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2add );
int HardRemoveAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2remove,
                              int *nNumCanceledCharges, BN_STRUCT *pBNS, BN_DATA *pBD );
int HardAddAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2add,
                              int *nNumCanceledCharges, BN_STRUCT *pBNS, BN_DATA *pBD );
int HardRemoveHplusNP( inp_ATOM *at, int num_atoms, int bCancelChargesAlways, int *nNumCanceledCharges,
                       BN_AATG *pAATG, BN_STRUCT *pBNS, BN_DATA *pBD  );
int RemoveNPProtonsAndAcidCharges( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, BN_STRUCT *pBNS, BN_DATA *pBD );
Vertex GetPrevVertex( BN_STRUCT* pBNS, Vertex y, Edge *SwitchEdge, EdgeIndex *iuv );
int bIgnoreVertexNonTACN_atom( BN_STRUCT* pBNS, Vertex u, Vertex v );
int bIgnoreVertexNonTACN_group( BN_STRUCT* pBNS, Vertex v, Vertex w, Edge *SwitchEdge );
int bIsRemovedHfromNHaion( BN_STRUCT* pBNS, Vertex u, Vertex v );
int bIsAggressiveDeprotonation( BN_STRUCT* pBNS, Vertex v, Vertex w, Edge *SwitchEdge );

int bIsAtomTypeHard( inp_ATOM *at, int endpoint, int nType, int nMask, int nCharge );
int bIsHDonorAccAtomType( inp_ATOM *at, int endpoint, int *cSubType );
int bIsNegAtomType( inp_ATOM *at, int i, int *cSubType );

#if ( BNS_RAD_SEARCH == 1 )
int RegisterRadEndpoint( BN_STRUCT *pBNS, BN_DATA *pBD, Vertex u);
int cmp_rad_endpoints( const void *a1, const void *a2 );
int cmp_endpoints_rad( const void *a1, const void *a2 );
#endif

int bHasChargedNeighbor( inp_ATOM *at, int iat );
/*********************************************************************************/

/**** prim(v) is v' *****/

#define prim(v) (Vertex)((v)^1)

/*********************************************************************************/

#define SwitchEdge_Vert1(u) SwitchEdge[u][0]
#define SwitchEdge_Vert2(u) Get2ndEdgeVertex( pBNS, SwitchEdge[u] )
#define SwitchEdge_IEdge(u) SwitchEdge[u][1]

/*********************************************************************************/
/*********************************************************************************/

/* returns value > 0 if a bond has been changed */
int RestoreEdgeFlow( BNS_EDGE *edge, int delta, int bChangeFlow )
{
    /*flow1 = edge->flow;*/ /* output from BNS */
    switch ( bChangeFlow & BNS_EF_CHNG_RSTR ) {
    case 0: /* the flow has not been permitted to change inside the BNS */
        /* nothing to do */
        /*flow1 = edge->flow;*/    /* output from BNS, the original flow value */
        /*flow2 = flow1 + delta;*/ /* the flow would be changed to this value by the BNS if permitted */
        break;
    case BNS_EF_CHNG_FLOW: /* the flow has been changed by the BNS; update flow0 */
        /*flow2 = edge->flow;*/    /* output from BNS, the changed value */
        /*flow1 = flow2 - delta;*/ /* the original flow value before the BNS */
        edge->flow0 = edge->flow;   /* SAVE NEW EDGE FLOW AS THE INITIAL FLOW FROM CHEM. BONDS */
        break;
    case BNS_EF_CHNG_RSTR: /* the flow has been changed by the BNS; requested to change it back */
        /*flow2 = edge->flow;*/    /* output from BNS, the changed value */
        /*flow1 = flow2 - delta;*/ /* the original flow value before the BNS */
        edge->flow = edge->flow-delta;    /* CHANGE EDGE FLOW BACK (RESTORE) */
        break;
    case BNS_EF_RSTR_FLOW: /* the flow has not been permitted to change inside the BNS */
        /* nothing to do */
        /*flow1 = edge->flow;*/    /* output from BNS, the original flow value */
        /*flow2 = flow1 + delta;*/ /* the flow would be changed to this value by the BNS if permitted */
        break;
    }
    return 0;
}
/*********************************************************************************/

/* returns value > 0 if a bond has been changed; do not change flow */
int SetAtomBondType( BNS_EDGE *edge, U_CHAR *bond_type12, U_CHAR *bond_type21, int delta, int bChangeFlow )
{
    int    flow1, flow2, tmp, ret = 0;
    int   bond_mark, bond_type, new_bond_type;

    if ( !edge->pass || !bond_type21 )
        return 0;

    switch ( bChangeFlow & BNS_EF_CHNG_RSTR ) {
    case 0:                    /* the flow has not been permitted to change inside the BNS: simulated in case of check one bond */
    case BNS_EF_RSTR_FLOW:     /* the flow has not been permitted to change inside the BNS: obsolete mode, unexpected bChangeFlow */
        flow1 = edge->flow0;   /* output from BNS, the original (old) flow value */
        flow2 = flow1 + delta; /* the flow would be changed to this value by the BNS if permitted */
        break;
    case BNS_EF_CHNG_FLOW:     /* the flow has been changed by the BNS */
    case BNS_EF_CHNG_RSTR:     /* the flow has been changed by the BNS; requested to change it back */
        flow2 = edge->flow;    /* output from BNS, the changed (new) value */
        flow1 = edge->flow0;   /* the original flow (old) value before the BNS */
        break;
    default:
        return 0; /* added 2006-03-21 */
    }
 
    if ( (bChangeFlow & BNS_EF_CHNG_BONDS) && (bChangeFlow & BNS_EF_ALTR_NS) !=BNS_EF_ALTR_NS ) {
        /* set new bond types according to the new flow values */
        new_bond_type = flow2+BOND_SINGLE;
        if ( *bond_type12 != new_bond_type ) {
            *bond_type12 = *bond_type21 = new_bond_type;
            ret ++;
        }
    } else
    if ( bChangeFlow & BNS_EF_ALTR_BONDS ) {
        if ( flow1 == flow2 ) {
            goto exit_function;
        }
        /* update alternating bond information */
        if ( flow1 > flow2 ) {
            /* make sure flow2 > flow1 */
            tmp   = flow1;
            flow1 = flow2;
            flow2 = tmp;
        }
        bond_mark = 0;
        switch( bond_type = (*bond_type12 & BOND_TYPE_MASK) ) {
        
        case BOND_SINGLE:
        case BOND_DOUBLE:
        case BOND_TRIPLE:
            /* assume that the input bond type fits either flow1 or flow2 */
            if ( flow1 == 0 && flow2 == 1 ) {
                if ( bChangeFlow & BNS_EF_SET_NOSTEREO ) {
                    bond_mark = BOND_MARK_ALT12NS;
                    bond_type = BOND_ALT12NS;
                } else {
                    bond_mark = BOND_MARK_ALT12;
                    bond_type = BOND_ALTERN;
                }
            } else
            if ( flow1 == 0 && flow2 == 2 ) {
                bond_mark = BOND_MARK_ALT13;
                bond_type = BOND_ALT_13;
            } else
            if ( flow1 == 1 && flow2 == 2 ) {
                bond_mark = BOND_MARK_ALT23;
                bond_type = BOND_ALT_23;
            } else {
                return BNS_BOND_ERR; /* error */
            }
            break;
        case BOND_TAUTOM:
            if ( flow1 == 0 && flow2 == 1 ) {
                bond_mark = BOND_MARK_ALT12NS;
            } else {
                return BNS_BOND_ERR; /* error */
            }
            break;

        default:
            new_bond_type = bond_type;
            bond_mark = (*bond_type12 & BOND_MARK_MASK);
            switch( bond_mark ) {
            case BOND_MARK_ALT12:
                if ( (bChangeFlow & BNS_EF_SET_NOSTEREO) && flow1 == 0 && flow2 == 1 ) {
                    bond_mark = BOND_MARK_ALT12NS;
                    new_bond_type = BOND_ALT12NS;
                    break;
                }
            case BOND_MARK_ALT12NS:
                if ( flow1 == 2 || flow2 == 2 ) {
                    bond_mark = BOND_MARK_ALT123;
                    new_bond_type = BOND_ALT_123;
                }
                break;
            case BOND_MARK_ALT13:
                if ( flow1 == 1 || flow2 == 1 ) {
                    bond_mark = BOND_MARK_ALT123;
                    new_bond_type = BOND_ALT_123;
                }
                break;
            case BOND_MARK_ALT23:
                if ( flow1 == 0 || flow2 == 0 ) {
                    bond_mark = BOND_MARK_ALT123;
                    new_bond_type = BOND_ALT_123;
                }
                break;
            case BOND_MARK_ALT123:
                break;
            
            case 0: /* special case: second alt bond testing */
                if ( flow1 == 0 && flow2 == 1 ) {
                    bond_mark = BOND_MARK_ALT12;
                } else
                if ( flow1 == 0 && flow2 == 2 ) {
                    bond_mark = BOND_MARK_ALT13;
                } else
                if ( flow1 == 1 && flow2 == 2 ) {
                    bond_mark = BOND_MARK_ALT23;
                } else {
                    return BNS_BOND_ERR; /* error */
                }
                break;


            default:
                return BNS_BOND_ERR; /* error */
            }
            switch( bond_type ) {
            case BOND_TAUTOM:
                break;
            case BOND_ALTERN:
            case BOND_ALT12NS:
            case BOND_ALT_123:
            case BOND_ALT_13:
            case BOND_ALT_23:
                bond_type = new_bond_type;
                break;
            default:
                return BNS_BOND_ERR; /* error */
            }
        }
        new_bond_type = bond_type | bond_mark;
        if ( new_bond_type != *bond_type12 ) {
            *bond_type12 = *bond_type21 = new_bond_type;
            ret ++;
        }
    }
exit_function:
    return ret;
}
/*********************************************************************************/
int RunBalancedNetworkSearch( BN_STRUCT *pBNS, BN_DATA *pBD, int bChangeFlow )
{
    /* Run BNS until no aug pass is found */
    int pass, delta=0, nSumDelta;
    nSumDelta = 0;
    for ( pass = 0; pass < pBNS->max_altp; pass ++ ) {
        pBNS->alt_path = pBNS->altp[pass];
        pBNS->bChangeFlow = 0;
        delta=BalancedNetworkSearch ( pBNS, pBD, bChangeFlow );
        ReInitBnData( pBD );
        if ( 0 < delta ) {
            pBNS->num_altp ++;
            nSumDelta += abs( delta );
        } else {
            break;
        }
    }
    if ( IS_BNS_ERROR(delta) )
        return delta;
    return nSumDelta; /* number of eliminated pairs of "dots"  */
}
/*********************************************************************************/
int SetAtomRadAndChemValFromVertexCapFlow( BN_STRUCT *pBNS, inp_ATOM *atom, int v1 )
{
    BNS_VERTEX *vert = pBNS->vert + v1;
    inp_ATOM   *at   = atom + v1;
    S_CHAR      cValue;
    int  nChanges = 0;
    /* set only on the 1st pass */
    if ( !vert->st_edge.pass ) {
        return 0;
    }
    /* adjust chem_bonds_valence */
    cValue = at->chem_bonds_valence - at->valence;
    if ( cValue  >= 0 && cValue != vert->st_edge.flow ) { 
        at->chem_bonds_valence = at->valence + vert->st_edge.flow;
        nChanges ++;
    }
    /* adjast radical */
    switch ( vert->st_edge.cap - vert->st_edge.flow ) {
    case 0:
        cValue = 0;
        break;
    case 1:
        cValue = RADICAL_DOUBLET;
        break;
    case 2:
        cValue = RADICAL_TRIPLET;
        break;
    default:
        return BNS_BOND_ERR;
    }
    if ( cValue != at->radical ) {
        at->radical = cValue;
        nChanges ++;
    }

    return nChanges;
}
/*********************************************************************************/
int AddChangedAtHChargeBNS( inp_ATOM *at, int num_atoms, int nAtTypeTotals[], S_CHAR *mark )
{
    int i, mask, num;
    for ( i = 0, num = 0; i < num_atoms; i ++ ) {
        if ( mark[i] ) {
            mark[i] = 0;
#if( FIX_NORM_BUG_ADD_ION_PAIR == 1 )
             /* add ignoring adjacent charges */
            at[i].at_type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, -2 );
#else
            at[i].at_type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 0 );
#endif
            num ++;
        }
    }
    return num;
}
/************************************************************************************/
/* eliminate neutral representation ambiguity: replace (+)--N==(-) with (+)==N--(-) */
/* here (+) is positive charge group, (-) is negative charge group, N is N or P     */
/* this reduces possibility of ion pair -OH => -O(+) + H(+) creation instead of     */
/* removing H(+) from N or P                                                        */
/*==== Call this function after alt path was found and new flows have been set. ====*/
/************************************************************************************/
int EliminatePlusMinusChargeAmbiguity( BN_STRUCT *pBNS, int num_atoms )
{
    int       pass, i, v0, v1, v2, ineigh1, /*ineigh0,*/ /*ineigh2,*/ vLast, n, delta, ret, err = 0;
    BNS_EDGE *edge;
    int       nFound, k;
    
    for ( pass = pBNS->num_altp-1, ret = 0; 0 <= pass; pass -- ) {
        
        pBNS->alt_path = pBNS->altp[pass];
        v1 = ALTP_START_ATOM(pBNS->alt_path);
        n  = ALTP_PATH_LEN(pBNS->alt_path);
        delta  = ALTP_DELTA(pBNS->alt_path);
        vLast = ALTP_END_ATOM(pBNS->alt_path);
        v0 = v2 = NO_VERTEX; /* negative number */

        for ( i = 0; i < n; i ++, delta = -delta, v0 = v1, v1 = v2 /*, ineigh0 = ineigh1*/ ) {
            ineigh1 = ALTP_THIS_ATOM_NEIGHBOR(pBNS->alt_path, i);  /* v1->v2 neighbor */
            /*ineigh2 = ALTP_NEXT_ATOM_NEIGHBOR(pBNS->alt_path, i);*/  /* v2->v1 neighbor */
            edge = pBNS->edge + pBNS->vert[v1].iedge[ineigh1];
            /* follow the BN Structure, not the inp_ATOM, to take care of swithching to
               t-groups, c-groups or other fictitious edges/vertices
            */
            v2 = edge->neighbor12 ^ v1;
            if ( v1 < num_atoms &&
                 ( v0 >= num_atoms && ( pBNS->vert[v0].type & BNS_VERT_TYPE_C_GROUP ) ||
                   v2 >= num_atoms && ( pBNS->vert[v2].type & BNS_VERT_TYPE_C_GROUP ) ) ) {
                int        cgPos, cgNeg;
                int        neighPos = -1, neighNeg = -1;
                BNS_EDGE *edgePos, *edgeNeg;
                nFound = 0;
                for ( k = pBNS->vert[v1].num_adj_edges-1; k >= 0 && (neighPos < 0 || neighNeg < 0); k -- ) {
                    BNS_EDGE   *next_edge = pBNS->edge + pBNS->vert[v1].iedge[k];
                    int         v         = next_edge->neighbor12 ^ v1;
                    if ( pBNS->vert[v].type & BNS_VERT_TYPE_C_GROUP ) {
                        if ( pBNS->vert[v].type & BNS_VERT_TYPE_C_NEGATIVE ) {
                            cgNeg    = v;
                            neighNeg = k;
                            nFound ++;
                        } else {
                            cgPos    = v;
                            neighPos = k;
                            nFound ++;
                        }
                    }
                }
                if ( 2 == nFound && neighPos >= 0 && neighNeg >= 0 ) {
                    /* both c-groups have been found */
                    edgePos = pBNS->edge + pBNS->vert[v1].iedge[neighPos];
                    edgeNeg = pBNS->edge + pBNS->vert[v1].iedge[neighNeg];
                    if ( edgePos->flow < edgeNeg->flow ) {
                        /* ambiguity found; replace (+cg)--N==(-cg) with (+cg)==N--(-cg) */
                        int dflow = edgeNeg->flow - edgePos->flow;
                        
                        edgePos->flow += dflow;
                        pBNS->vert[cgPos].st_edge.cap  += dflow;
                        pBNS->vert[cgPos].st_edge.flow += dflow;
                        
                        edgeNeg->flow -= dflow;
                        pBNS->vert[cgNeg].st_edge.cap  -= dflow;
                        pBNS->vert[cgNeg].st_edge.flow -= dflow;
                        ret ++;
                    }
                }
            }
        }
        
        if ( v2 != vLast ) {
            err = BNS_PROGRAM_ERR;
        }
    }
    return err? err : ret;
}
/**********************************************************************************************************
  nNum2Remove = number of H to remove
  num_H       = is number of H before the removal, including explicit H 
***********************************************************************************************************/
int AddOrRemoveExplOrImplH( int nDelta, inp_ATOM *at, int num_atoms, AT_NUMB at_no, T_GROUP_INFO *t_group_info )
{
    int       i, iso, nNumExplicit2Implicit;
    int       nNum2Remove, nNumRemovedExplicitH;
    S_CHAR    num_iso_H[NUM_H_ISOTOPES];
    int       tot_num_iso_H, num_H;
    inp_ATOM *at_H;

    if ( !nDelta ) {
        return 0;
    }
    /* add */
    if ( nDelta > 0 ) {
        at[at_no].num_H += nDelta;
        t_group_info->tni.nNumRemovedProtons --;
        return nDelta;
    }
    /* remove */
    nNum2Remove            = -nDelta;
    nNumRemovedExplicitH   = t_group_info->tni.nNumRemovedExplicitH; /* number of explicit H saved separately in
                                                                       at[num_atoms+i], i=0..nNumRemovedExplicitH-1 */
    tot_num_iso_H          = NUM_ISO_H(at,at_no);
    num_H                  = at[at_no].num_H;
    /*
    tot_num_iso_H          = NUM_ISO_H(at,at_no);
    num_H                  = at[at_no].num_H;
    nNumAtomExplicitH      = 0;
    nNumRemovedExplicitH   = t_group_info->tni.nNumRemovedExplicitH;
    tot_num_explicit_iso_H = 0;
    */
    at_H                   = at + num_atoms;
    memcpy( num_iso_H, at[at_no].num_iso_H, sizeof(num_iso_H)); 
    /*  Remove all explicit H, otherwise a false stereo can occur.
        Example: remove H(+) from the following substructure:

               H                                            H
       A      /                                      A     / 
        >X==N(+)    produces false stereogenic bond:  >X==N 
       B      \                                      B
               H
        
        To avoid this effect all explicit H atoms must be removed
    */
    nNumExplicit2Implicit = 0;
    for ( i = 0; i < nNumRemovedExplicitH; ) {
        if ( at_H[i].neighbor[0] == at_no ) {
            int m, k, orig_no = at_H[i].orig_at_number;
            nNumRemovedExplicitH --;
            nNumExplicit2Implicit ++;
            if ( nNumRemovedExplicitH > i ) {
                inp_ATOM at_i = at_H[i];
                memmove( at_H+i, at_H+i+1, sizeof(at_H[0])*(nNumRemovedExplicitH-i) );
                at_H[nNumRemovedExplicitH] = at_i; /* save removed H (for debugging purposes?) */
            }
            /* adjust 0D parities */
            if ( at[at_no].sb_parity[0] ) {
                for ( m = 0; m < MAX_NUM_STEREO_BONDS && at[at_no].sb_parity[m]; m ++ ) {
                    if ( at[at_no].sn_orig_at_num[m] == orig_no ) {
#ifdef _DEBUG
                        if ( at[at_no].sn_ord[m] >= 0 ) {
                            int stop = 1; /* sb maintenance error */
                        }
#endif
                        if ( at[at_no].valence >= MIN_NUM_STEREO_BOND_NEIGH ) {
                            at[at_no].sn_ord[m] = k = (at[at_no].sb_ord[m]==0);
                            at[at_no].sn_orig_at_num[m] = at[(int)at[at_no].neighbor[k]].orig_at_number;
                            if ( ATOM_PARITY_WELL_DEF( at[at_no].sb_parity[m] ) ) {
                                at[at_no].sb_parity[m] = 3 - at[at_no].sb_parity[m];
                            }
                        } else {
                            at[at_no].sn_ord[m] = -99; /* no sb neighbor exists anymore */
                            at[at_no].sn_orig_at_num[m] = 0;
                            if ( ATOM_PARITY_WELL_DEF( at[at_no].sb_parity[m] ) ) {
                                int pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
                                if ( 0 < get_opposite_sb_atom( at, at_no, at[at_no].sb_ord[m],
                                                  &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord ) ) {
                                    at[at_no].sb_parity[m] =
                                    at[pnxt_atom].sb_parity[pinxt_sb_parity_ord] = AB_PARITY_UNDF;
                                }
#ifdef _DEBUG
                                else {
                                    int stop = 1; /* sb maintenance error */
                                }
#endif
                            }
                        }
                    }
                }
            }
            /* do not increment i here: we have shifted next at_H[] element
               to the ith position and decremented nNumRemovedExplicitH */
        } else {
            i ++;
        }
    }

    for ( iso = -1; iso < NUM_H_ISOTOPES && 0 < nNum2Remove; iso ++ ) {
        /* each pass removes up to one H */
        if ( iso < 0 ) {
            /* try to remove non-isotopic */
            while ( tot_num_iso_H < num_H && 0 < nNum2Remove ) {
                /* non-isotopic H exists */
                num_H --;
                t_group_info->tni.nNumRemovedProtons ++;
                nNum2Remove --;
            }
        } else {
            /* remove isotopic */
            while ( num_iso_H[iso] && num_H && 0 < nNum2Remove ) {
                /* isotopic H exists */
                num_H --;
                num_iso_H[iso] --;
                t_group_info->tni.nNumRemovedProtonsIsotopic[iso] ++;
                t_group_info->tni.nNumRemovedProtons ++;
                nNum2Remove --;
            }
        }
    }
#if ( bRELEASE_VERSION != 1 )
    if ( nNum2Remove ) {
        int stop = 1; /* <BRKPT> Program error */
    }
#endif
    if ( nDelta + nNum2Remove < 0 ) {
        at[at_no].num_H = num_H;
        memcpy( at[at_no].num_iso_H, num_iso_H, sizeof(at[0].num_iso_H));
        t_group_info->tni.nNumRemovedExplicitH = nNumRemovedExplicitH;
    }
    return nDelta + nNum2Remove;
}

/*********************************************************************************/
int SubtractOrChangeAtHChargeBNS( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms,
                                 int nAtTypeTotals[], S_CHAR *mark, T_GROUP_INFO *t_group_info, int bSubtract )
{
    int       pass, i, v0, v1, v2, ineigh1, /*ineigh2,*/ vLast, n, delta, ret, err = 0;
    BNS_EDGE *edge;
    int       nDeltaH, nDeltaCharge;
    int       mask, type;
    
    for ( pass = pBNS->num_altp-1, ret = 0; 0 <= pass; pass -- ) {
        
        pBNS->alt_path = pBNS->altp[pass];
        v1 = ALTP_START_ATOM(pBNS->alt_path);
        n  = ALTP_PATH_LEN(pBNS->alt_path);
        delta  = ALTP_DELTA(pBNS->alt_path);
        vLast = ALTP_END_ATOM(pBNS->alt_path);
        v0 = v2 = NO_VERTEX;

        for ( i = 0; i < n; i ++, delta = -delta, v0 = v1, v1 = v2 ) {
            ineigh1 = ALTP_THIS_ATOM_NEIGHBOR(pBNS->alt_path, i);  /* v1->v2 neighbor */
            /*ineigh2 = ALTP_NEXT_ATOM_NEIGHBOR(pBNS->alt_path, i);*/  /* v2->v1 neighbor */
            edge = pBNS->edge + pBNS->vert[v1].iedge[ineigh1];
            /* follow the BN Structure, not the inp_ATOM, to take care of swithching to
               t-groups, c-groups or other fictitious edges/vertices
            */
            v2 = edge->neighbor12 ^ v1;
            if ( v1 < num_atoms && (v0 >= num_atoms || v2 >= num_atoms) ) {
                nDeltaH = nDeltaCharge = 0;
                if ( v0 >= num_atoms ) {
                    /* delta(v0-v1) = -delta(v1-v2) along the alternating path */
                    if ( pBNS->vert[v0].type & BNS_VERT_TYPE_TGROUP ) {
                        nDeltaH -= delta;
                    } else
                    if ( pBNS->vert[v0].type & BNS_VERT_TYPE_C_GROUP ) {
                        nDeltaCharge += delta;
                    }
                }
                if ( v2 >= num_atoms ) {
                    if ( pBNS->vert[v2].type & BNS_VERT_TYPE_TGROUP ) {
                        nDeltaH += delta;
                    } else
                    if ( pBNS->vert[v2].type & BNS_VERT_TYPE_C_GROUP ) {
                        nDeltaCharge -= delta;
                    }
                }
                if ( nDeltaH || nDeltaCharge ) {
                    if ( bSubtract ) {
                        if ( !mark[v1] ) {
                            /* first time the atom has been encountered: subtract */
#if( FIX_NORM_BUG_ADD_ION_PAIR == 1 )
                            type   = GetAtomChargeType( at, v1, nAtTypeTotals, &mask, 2 );
#else
                            type   = GetAtomChargeType( at, v1, nAtTypeTotals, &mask, 1 );
#endif
                            ret ++; /* number of changed atoms */
                            mark[v1] ++;
                        }
                    } else { /* Change */
                        at[v1].charge += nDeltaCharge;
                        if ( nDeltaH ) {
                            AddOrRemoveExplOrImplH( nDeltaH, at, num_atoms, (AT_NUMB)v1, t_group_info );
                        }
                        ret ++; /* number of changed atoms */
                    }
                } 
            }
        }
        
        if ( v2 != vLast ) {
            err = BNS_PROGRAM_ERR;
        }
    }
    return err? err : ret;
}
/*********************************************************************************/
int SetBondsFromBnStructFlow( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bChangeFlow0 )
{
    int       pass, i, v0, v1, v2, ineigh1, ineigh2, vLast, n, delta, ret, ret_val, err = 0;
    BNS_EDGE *edge;
    int       bMovingRad  = 0, bChangeFlowAdd;
    int       bChangeFlow = (bChangeFlow0 & ~BNS_EF_SET_NOSTEREO);
    /*
    bCheckMovingRad = (bChangeFlow & BNS_EF_ALTR_NS) == BNS_EF_ALTR_NS &&
                      pBNS->tot_st_cap > pBNS->tot_st_flow;
    */
    for ( pass = pBNS->num_altp-1, ret = 0; 0 <= pass; pass -- ) {
        pBNS->alt_path = pBNS->altp[pass];
        v1 = ALTP_START_ATOM(pBNS->alt_path);
        n  = ALTP_PATH_LEN(pBNS->alt_path);
        delta  = ALTP_DELTA(pBNS->alt_path);
        vLast = ALTP_END_ATOM(pBNS->alt_path);
        if ( (bChangeFlow0 & BNS_EF_SET_NOSTEREO) &&
             (pBNS->vert[v1].st_edge.cap0 > pBNS->vert[v1].st_edge.flow0 ||
              pBNS->vert[vLast].st_edge.cap0 > pBNS->vert[vLast].st_edge.flow0 ) ) {
            bMovingRad ++;
            bChangeFlowAdd = BNS_EF_SET_NOSTEREO;
            ret |= 2;
        } else {
            bChangeFlowAdd = 0;
        }
        /* start vertex */      
        if ( (bChangeFlow & BNS_EF_CHNG_RSTR) == BNS_EF_CHNG_RSTR) {
            /* restore s-v1 edge flow to the BNS this pass input value */
            ; /*pBNS->vert[v1].st_edge.flow -= delta;*/
        } else
        if ( (bChangeFlow & BNS_EF_SAVE_ALL) == BNS_EF_SAVE_ALL ) {
            if ( v1 < num_atoms ) {
                /* will produce wrong result if called for v1 next time? */
                ret_val = SetAtomRadAndChemValFromVertexCapFlow( pBNS, at, v1 );
                if ( ret_val < 0 ) {
                    err = BNS_PROGRAM_ERR;
                } else {
                    ret |= (ret_val > 0);
                }
            }
            /*pBNS->vert[v1].st_edge.flow0 = pBNS->vert[v1].st_edge.flow;*/
        }
        pBNS->vert[v1].st_edge.pass = 0;
        
        v0 = v2 = NO_VERTEX;
        for ( i = 0; i < n; i ++, delta = -delta, v0 = v1, v1 = v2 ) {
            ineigh1 = ALTP_THIS_ATOM_NEIGHBOR(pBNS->alt_path, i);  /* v1->v2 neighbor */
            ineigh2 = ALTP_NEXT_ATOM_NEIGHBOR(pBNS->alt_path, i);  /* v2->v1 neighbor */
            edge = pBNS->edge + pBNS->vert[v1].iedge[ineigh1];
            /* follow the BN Structure, not the inp_ATOM, to take care of swithching to
               t-groups, c-groups or other fictitious edges/vertices
            */

            v2 = edge->neighbor12 ^ v1;

            /* change at->chem_bonds_valence 2004-03-08 */
            if ( (bChangeFlow & BNS_EF_CHNG_BONDS) && v1 < num_atoms ) {
                if ( v0 >= num_atoms && v2 < num_atoms ) {
                    at[v1].chem_bonds_valence += delta; /* change in v1-v2 bond order */
                } else
                if ( v0 < num_atoms && v2 >= num_atoms && v0 != NO_VERTEX ) {
                    at[v1].chem_bonds_valence -= delta; /* change in v0-v1 bond order */
                }
            }
            
            if ( !edge->pass )
                continue;
    
            if ( v1 < num_atoms && ineigh1 < at[v1].valence &&
                 v2 < num_atoms && ineigh2 < at[v2].valence ) {
                if ( (bChangeFlow0 & BNS_EF_ALTR_NS )==BNS_EF_ALTR_NS &&
                     (bChangeFlow0 & BNS_EF_SAVE_ALL)==BNS_EF_SAVE_ALL ) {
                    /* 2004-07-02 special mode: save new ring bonds and mark as non-stereo non-ring bonds */
                    if ( at[v1].nRingSystem != at[v2].nRingSystem ) {
                        /* non-ring bond (bridge) */
                        bChangeFlowAdd = BNS_EF_ALTR_NS;
                    } else {
                        /* ring bond */
                        bChangeFlowAdd = 0;
                    }
                }
                /* change bonds on the first pass only: in this case all flow correspond to the BNS output */
                ret_val = SetAtomBondType( edge, &at[v1].bond_type[ineigh1], &at[v2].bond_type[ineigh2], delta, bChangeFlow | bChangeFlowAdd );
                if ( ret_val < 0 ) {
                    err = BNS_PROGRAM_ERR;
                } else {
                    ret |= (ret_val > 0);
                }
            }
            edge->pass = 0;
        }
        
        if ( v2 != vLast ) {
            err = BNS_PROGRAM_ERR;
        } else
        if ( (bChangeFlow & BNS_EF_CHNG_RSTR) == BNS_EF_CHNG_RSTR) {
            /* restore v2-t edge flow to the BNS this pass input value */
            /* "+=" instead of "-=" explanation: delta must have same sign as at the last edge */
            ; /*pBNS->vert[v2].st_edge.flow += delta; */
        } else
        if ( (bChangeFlow & BNS_EF_SAVE_ALL) == BNS_EF_SAVE_ALL ) {
            if ( v2 < num_atoms ) {
                ret_val = SetAtomRadAndChemValFromVertexCapFlow( pBNS, at, v2 );
                if ( ret_val < 0 ) {
                    err = BNS_PROGRAM_ERR;
                } else {
                    ret |= (ret_val > 0);
                }
            }
            /*pBNS->vert[v2].st_edge.flow0 = pBNS->vert[v2].st_edge.flow;*/
        }
        pBNS->vert[v2].st_edge.pass = 0;

    }
    return err? err : ret;
}
/*********************************************************************************/
int MarkAtomsAtTautGroups(  BN_STRUCT *pBNS, int num_atoms, BN_AATG *pAATG, int nEnd1, int nEnd2 )
{
    int       pass, i, j, v1, v2, ineigh1, ineigh2, vLast, vFirst, n, delta, err = 0;
    BNS_EDGE *edge;
    S_CHAR    cDelta[MAX_ALT_AATG_ARRAY_LEN];
    AT_NUMB   nVertex[MAX_ALT_AATG_ARRAY_LEN];
    int       nLenDelta = 0, last_i, nNumFound;

    for ( pass = pBNS->num_altp-1; 0 <= pass; pass -- ) {
        pBNS->alt_path = pBNS->altp[pass];
        vFirst =
        v1 = ALTP_START_ATOM(pBNS->alt_path);
        n  = ALTP_PATH_LEN(pBNS->alt_path);
        delta  = ALTP_DELTA(pBNS->alt_path);
        vLast  = ALTP_END_ATOM(pBNS->alt_path);
        v2     = NO_VERTEX;
        pAATG->nNumFound = 0; /* initialize */

        if ( nEnd1 != vFirst && nEnd1 != vLast ) {
            nEnd1 = -1; /* really not the end */
        }
        if ( nEnd2 != vFirst && nEnd2 != vLast ) {
            nEnd2 = -1; /* really not the end */
        }

        for ( i = 0; i < n; i ++, delta = -delta, v1 = v2 ) {
            ineigh1 = ALTP_THIS_ATOM_NEIGHBOR(pBNS->alt_path, i);  /* v1->v2 neighbor */
            ineigh2 = ALTP_NEXT_ATOM_NEIGHBOR(pBNS->alt_path, i);  /* v2->v1 neighbor */
            edge = pBNS->edge + pBNS->vert[v1].iedge[ineigh1];
            /* follow the BN Structure, not the inp_ATOM, to take care of swithching to
               t-groups, c-groups or other fictitious edges/vertices
            */
            v2 = edge->neighbor12 ^ v1;
            /*
            if ( v1 < num_atoms && v2 < num_atoms ) {
                continue;
            }
            */
            /* BNS increased edge flow by delta */
            if ( v1 >= num_atoms && 
                 ((pBNS->vert[v1].type & BNS_VERT_TYPE_TGROUP)||(pBNS->vert[v1].type & BNS_VERT_TYPE_TEMP)) &&
                 0 <= v2 && v2 <  num_atoms && (pBNS->vert[v2].type & BNS_VERT_TYPE_ATOM  ) ) {
                /*
                if ( !(pAATG->nMarkedAtom[v2] & AATG_MARK_IN_PATH) ) {
                    pAATG->nMarkedAtom[v2] |= AATG_MARK_IN_PATH;
                    pAATG->nNumFound ++;
                }
                */
                /* BNS increased bond order in  v1(t-group)-v2(atom) by delta: added delta attachments */
                if ( nLenDelta < MAX_ALT_AATG_ARRAY_LEN ) {
                    cDelta[nLenDelta]  = delta;
                    nVertex[nLenDelta] = v2;
                    nLenDelta ++;
                }
            } else
            if ( v2 >= num_atoms &&
                 ((pBNS->vert[v2].type & BNS_VERT_TYPE_TGROUP)||(pBNS->vert[v2].type & BNS_VERT_TYPE_TEMP)) &&
                 0 <= v1 && v1 <  num_atoms && (pBNS->vert[v1].type & BNS_VERT_TYPE_ATOM  ) ) {
                /*
                if ( !(pAATG->nMarkedAtom[v1] & AATG_MARK_IN_PATH) ) {
                    pAATG->nMarkedAtom[v1] |= AATG_MARK_IN_PATH;
                    pAATG->nNumFound ++;
                }
                */
                /* BNS increased bond order in  v1(atom)-v2(t-group) by delta: added delta attachments */
                if ( nLenDelta < MAX_ALT_AATG_ARRAY_LEN ) {
                    cDelta[nLenDelta]  = delta;
                    nVertex[nLenDelta] = v1;
                    nLenDelta ++;
                }
            } else
            /* special case when the testing 'dot' was placed on an atom (should be nEnd1 only) */
            if ( 0 <= v1 && v1 == nEnd1 || v1 == nEnd2 && 0 <= v2 && v2 < num_atoms ) {
                if ( nLenDelta < MAX_ALT_AATG_ARRAY_LEN ) {
                    cDelta[nLenDelta]  = -delta;
                    nVertex[nLenDelta] = v1;
                    nLenDelta ++;
                }
            } else
            if ( 0 <= v2 && v2 == nEnd1 || v2 == nEnd2 && 0 <= v1 && v1 < num_atoms ) {
                if ( nLenDelta < MAX_ALT_AATG_ARRAY_LEN ) {
                    cDelta[nLenDelta]  = -delta;
                    nVertex[nLenDelta] = v2;
                    nLenDelta ++;
                }
            }
        }
        
        if ( v2 != vLast ) {
            err = BNS_PROGRAM_ERR;
        } else {
            last_i = -1;
            nNumFound = 0;
            /* first run */
            for ( i = 1, j = 0; i < nLenDelta; j = i ++ ) {
                /* ignore sequences (-1,+1) and (+1,-1) in cDelta[] because they  */
                /* describe ordinary aug. paths of moving a single attachment     */
                /* we are looking for aug. paths describing movement of 2 or more */
                if ( cDelta[j] > 0 && cDelta[i] > 0 ||
                     cDelta[j] < 0 && cDelta[i] < 0 ) {
                    if ( j == last_i ) {
                        /* three attachments moved */
                        return 0;
                    }
                    v1 = nVertex[j];
                    if ( !(pAATG->nMarkedAtom[v1] & AATG_MARK_IN_PATH) ) {
                        nNumFound ++;
                    }
                    v2 = nVertex[i];
                    if ( !(pAATG->nMarkedAtom[v2] & AATG_MARK_IN_PATH) ) {
                        nNumFound ++;
                    }
                    last_i = i;
                }
            }
            if ( !nNumFound ) {
                return 0;
            }
            if ( nNumFound > 4 ) {
                return 0;
            }
            if ( nNumFound < 4 ) {
                return 0;
            }
            /* second run */
            for ( i = 1, j = 0; i < nLenDelta; j = i ++ ) {
                /* ignore sequences (-1,+1) and (+1,-1) in cDelta[] because they  */
                /* describe ordinary aug. paths of moving a single attachment     */
                /* we are looking for aug. paths describing movement of 2 or more */
                if ( cDelta[j] > 0 && cDelta[i] > 0 ||
                     cDelta[j] < 0 && cDelta[i] < 0 ) {
                    v1 = nVertex[i-1];
                    if ( !(pAATG->nMarkedAtom[v1] & AATG_MARK_IN_PATH) ) {
                        pAATG->nMarkedAtom[v1] |= AATG_MARK_IN_PATH;
                        pAATG->nNumFound ++;
                    }
                    v2 = nVertex[i];
                    if ( !(pAATG->nMarkedAtom[v2] & AATG_MARK_IN_PATH) ) {
                        pAATG->nMarkedAtom[v2] |= AATG_MARK_IN_PATH;
                        pAATG->nNumFound ++;
                    }
                }
            }
        }
    }
    return err? err : pAATG->nNumFound;
}
/*********************************************************************************/
int RestoreBnStructFlow( BN_STRUCT *pBNS, int bChangeFlow )
{
    int       pass, i, v1, v2, ineigh1, ineigh2, vLast, n, delta, ret, err = 0;
    BNS_EDGE *edge;
    
    for ( pass = pBNS->num_altp - 1, ret = 0; 0 <= pass; pass -- ) {
        pBNS->alt_path = pBNS->altp[pass];
        v1 = ALTP_START_ATOM(pBNS->alt_path);
        n  = ALTP_PATH_LEN(pBNS->alt_path);
        delta  = ALTP_DELTA(pBNS->alt_path);
        vLast  = ALTP_END_ATOM(pBNS->alt_path);
        v2     = NO_VERTEX;
        /* starting vertex */
        if ( (bChangeFlow & BNS_EF_CHNG_RSTR) == BNS_EF_CHNG_RSTR) {
            pBNS->vert[v1].st_edge.flow -= delta; /* restore s-v1 edge flow to the BNS input value */
        } else
        if ( (bChangeFlow & BNS_EF_SAVE_ALL) == BNS_EF_SAVE_ALL ) {
            pBNS->vert[v1].st_edge.flow0 = pBNS->vert[v1].st_edge.flow;
        }
        /* augmenting path edges */
        for ( i = 0; i < n; i ++, delta = -delta, v1 = v2 ) {
            ineigh1 = ALTP_THIS_ATOM_NEIGHBOR(pBNS->alt_path, i);  /* v1->v2 neighbor */
            ineigh2 = ALTP_NEXT_ATOM_NEIGHBOR(pBNS->alt_path, i);  /* v2->v1 neighbor */
            edge = pBNS->edge + pBNS->vert[v1].iedge[ineigh1];
            v2 = edge->neighbor12 ^ v1;
            RestoreEdgeFlow( edge, delta, bChangeFlow );
            edge->pass = 0;
        }
        /* ending vertex */
        if ( v2 != vLast ) {
            err = BNS_PROGRAM_ERR;
        } else
        if ( (bChangeFlow & BNS_EF_CHNG_RSTR) == BNS_EF_CHNG_RSTR) {
            /* restore v2-t edge flow to the original value */
            /* "+=" instead of "-=" explanation: delta must have same sign as at the last edge */
            pBNS->vert[v2].st_edge.flow += delta; 
        } else
        if ( (bChangeFlow & BNS_EF_SAVE_ALL) == BNS_EF_SAVE_ALL ) {
            pBNS->vert[v2].st_edge.flow0 = pBNS->vert[v2].st_edge.flow;
        }

    }
    return err? err : ret;
}
/***************************************************************************************/
int bNeedToTestTheFlow( int bond_type, int nTestFlow, int bTestForNonStereoBond )
{
    int nBondType   = ( BOND_TYPE_MASK & bond_type );
    int nBondAttrib = ( BOND_MARK_MASK & bond_type );

    if ( bTestForNonStereoBond ) {
        if ( nBondAttrib || nBondType == BOND_ALTERN || nBondType == BOND_ALT12NS ) {
            switch( nTestFlow ) {
            case 0:  /* single: can be 1 (single)?  */
                if ( nBondAttrib == BOND_MARK_ALT12NS||
                     nBondAttrib == BOND_MARK_ALT123 ||
                     nBondAttrib == BOND_MARK_ALT13  ) {
                    return 0; /* yes, already checked */
                }
                break;
             
            case 1:  /* double: can be 2 (double)? */
                if ( nBondAttrib == BOND_MARK_ALT12NS||
                     nBondAttrib == BOND_MARK_ALT123 ||
                     nBondAttrib == BOND_MARK_ALT23  ) {
                    return 0; /* yes, already checked */
                }
                break;
            case 2:  /* triple: can be 3 (triple)?  */
                if ( nBondAttrib == BOND_MARK_ALT13  ||
                     nBondAttrib == BOND_MARK_ALT123 ||
                     nBondAttrib == BOND_MARK_ALT23  ) {
                    return 0; /* yes, already checked */
                }
                break;
            }
        }
    } else {
        if ( nBondAttrib || nBondType == BOND_ALTERN || nBondType == BOND_ALT12NS ) {
            switch( nTestFlow ) {
            case 0:  /* single: can be 1 (single)?  */
                if ( nBondAttrib == BOND_MARK_ALT12  ||
                     nBondAttrib == BOND_MARK_ALT12NS||
                     nBondAttrib == BOND_MARK_ALT123 ||
                     nBondAttrib == BOND_MARK_ALT13  ) {
                    return 0;
                }
                break;
             
            case 1:  /* double: can be 2 (double)? */
                if ( nBondAttrib == BOND_MARK_ALT12  ||
                     nBondAttrib == BOND_MARK_ALT12NS||
                     nBondAttrib == BOND_MARK_ALT123 ||
                     nBondAttrib == BOND_MARK_ALT23  ) {
                    return 0; /* yes */
                }
                break;
            case 2:  /* triple: can be 3 (triple)?  */
                if ( nBondAttrib == BOND_MARK_ALT13  ||
                     nBondAttrib == BOND_MARK_ALT123 ||
                     nBondAttrib == BOND_MARK_ALT23  ) {
                    return 0;
                }
                break;
            }
        }
    }
    return 1;
}
/***********************************************************************************/
int nBondsValenceInpAt( const inp_ATOM *at, int *nNumAltBonds, int *nNumWrongBonds )
{
    int j, bond_type, nBondsValence = 0, nAltBonds = 0, nNumWrong = 0;
    for ( j = 0; j < at->valence; j ++ ) {
        bond_type = at->bond_type[j] & BOND_TYPE_MASK;
        switch( bond_type ) {
            case 0:  /* for structure from InChI reconstruction */
            case BOND_SINGLE:
            case BOND_DOUBLE:
            case BOND_TRIPLE:
                nBondsValence += bond_type;
                break;
            case BOND_ALTERN:
                nAltBonds ++;
                break;
            default:
                nNumWrong ++;
        }
    }
    switch ( nAltBonds ) {
    case 0:
        break;
    case 1:
        nBondsValence += 1; /* 1 or greater than 3 is wrong */
        nNumWrong ++;
        break;
    default:
        nBondsValence += nAltBonds+1;
        break;
    }
    if ( nNumAltBonds ) *nNumAltBonds = nAltBonds;
    if ( nNumWrongBonds ) *nNumWrongBonds = nNumWrong;
    return nBondsValence;
}
/***************************************************************************************/
/* if radical or has aromatic bonds then augment to the lowest "multiplicity" */
int BnsAdjustFlowBondsRad( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms )
{
    int  bError=0, nOrigDelta=0, ret, num_removed;

#if( CHECK_AROMBOND2ALT == 1 )    
    char *pcValMinusBondsVal = NULL;
    int  i, nValMinusBondsVal, nAltBonds, bIgnore;
    
    /* find valence excess (it may only be due to aromatic bonds) */
    for ( i = 0; i < num_atoms; i ++ ) {
        nValMinusBondsVal = (int)at[i].chem_bonds_valence - nBondsValenceInpAt( at+i, &nAltBonds, &bIgnore );
        bIgnore += (nAltBonds > 3);
        if ( !bIgnore && nValMinusBondsVal > 0 ) {
            if ( !pcValMinusBondsVal && 
                 !(pcValMinusBondsVal = (char *)inchi_calloc(num_atoms, sizeof(pcValMinusBondsVal[0])))) {
                bError = BNS_OUT_OF_RAM;
                goto exit_function;
            }
            /* mark atoms that have extra unsatisfied valence due to aromatic bonds */
            pcValMinusBondsVal[i] = nValMinusBondsVal + (at[i].radical == RADICAL_DOUBLET);
        }
    }
#endif /* CHECK_AROMBOND2ALT */

    /* match bonds to valences */
    do {
        num_removed = 0;
        ret = RunBalancedNetworkSearch( pBNS, pBD, BNS_EF_CHNG_FLOW );
        if ( IS_BNS_ERROR( ret ) ) {
            bError = ret;
        } else {
            nOrigDelta += ret;
            num_removed = pBNS->num_altp; /* number of augmenting paths */
            if ( ret > 0 ) {
                /* save new bonds in at[] and flows in pBNS and at[] */
                ret = SetBondsFromBnStructFlow( pBNS, at, num_atoms, BNS_EF_SAVE_ALL ); /* must include 1: 5=(4|1) */
                if ( IS_BNS_ERROR( ret ) ) {
                    bError = ret;
                }
                ret = RestoreBnStructFlow( pBNS, BNS_EF_SAVE_ALL ); /* must include 1: 5=(4|1) */
                if ( IS_BNS_ERROR( ret ) ) {
                    bError = ret;
                }
            }
            ReInitBnStructAltPaths( pBNS );
        }
    } while ( num_removed && num_removed == pBNS->max_altp && !bError );

#if( CHECK_AROMBOND2ALT == 1 )    
    /* check whether aromatic bonds have been replaces with alternating bonds */
    if ( !bError && pcValMinusBondsVal ) {
        for ( i = 0; i < num_atoms; i ++ ) {
            if ( !pcValMinusBondsVal[i] )
                continue;
            nValMinusBondsVal = (int)at[i].chem_bonds_valence - nBondsValenceInpAt( at+i, &nAltBonds, &bIgnore );
            if ( bIgnore ||
                 1 != (int)pcValMinusBondsVal[i] - (nValMinusBondsVal + (at[i].radical == RADICAL_DOUBLET)) ) {
                /* radical excess has not been reduced */
                bError = BNS_ALTBOND_ERR;
                break;
            }
        }
    }

exit_function:
    if ( pcValMinusBondsVal ) inchi_free( pcValMinusBondsVal );
#endif /* CHECK_AROMBOND2ALT */

    return bError? bError : nOrigDelta;
}

/***************************************************************************************/
int BnsTestAndMarkAltBonds(  BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms, BNS_FLOW_CHANGES *fcd, int bChangeFlow, int nBondTypeToTest )
{
    int ret, iat, ineigh, neigh;
    int nMinFlow, nMaxFlow, nTestFlow, nCurFlow;
    int iedge, bSuccess, bError, nDots, nChanges, bTestForNonStereoBond;
    /* Normalize bonds and find tautomeric groups */
    bError   = 0;
    nChanges = 0;
    bTestForNonStereoBond = pBNS->tot_st_cap > pBNS->tot_st_flow;
    for ( iat = 0; iat < num_atoms && !bError; iat ++ ) {
        for ( ineigh = 0; ineigh < at[iat].valence && !bError; ineigh ++ ) {
            neigh = at[iat].neighbor[ineigh];
            if ( neigh < iat )
                continue; /* we have already tested the bond */
            iedge = pBNS->vert[iat].iedge[ineigh];
            if ( IS_FORBIDDEN(pBNS->edge[iedge].forbidden, pBNS) )
                continue;
            if ( nBondTypeToTest && (at[iat].bond_type[ineigh] & BOND_TYPE_MASK) != nBondTypeToTest )
                continue;
            nMinFlow = nMinFlow2Check( pBNS, iedge );
            nMaxFlow = nMaxFlow2Check( pBNS, iedge );
            nCurFlow = nCurFlow2Check( pBNS, iedge );
            if ( nMinFlow == nMaxFlow ) {
                if ( nMaxFlow && bTestForNonStereoBond ) {
                    nTestFlow = nMaxFlow - (int)(pBNS->tot_st_cap - pBNS->tot_st_flow); /* temporary use of nTestFlow */
                    nMinFlow = inchi_max( 0, nTestFlow );
                } else {
                    continue;
                }
            }
            for ( nTestFlow = nMinFlow; nTestFlow <= nMaxFlow && !bError; nTestFlow ++ ) {
                if ( nTestFlow == nCurFlow )
                    continue;
                if ( !bNeedToTestTheFlow( at[iat].bond_type[ineigh], nTestFlow, bTestForNonStereoBond ) )
                    continue;
                bSuccess = 0;
                nDots = bSetFlowToCheckOneBond( pBNS, iedge, nTestFlow, fcd );
                if ( IS_BNS_ERROR(nDots) ) {
                    if ( nDots == BNS_CANT_SET_BOND ) {
                        ret = bRestoreFlowAfterCheckOneBond( pBNS, fcd );
                        if ( !IS_BNS_ERROR( ret ) ) {
                            continue;
                        }
                    }
                    bError = nDots;
                } else
                if ( nDots > 0 ) {
                    ret = RunBalancedNetworkSearch( pBNS, pBD, bChangeFlow );
                    if ( IS_BNS_ERROR( ret ) ) {
                        bError = ret;
                    } else
                    if ( ret > 0 ) {
                        if ( 2*ret == nDots ) {
                            ret = bSetBondsAfterCheckOneBond( pBNS, fcd, nTestFlow, at, num_atoms, bChangeFlow );
                            if ( IS_BNS_ERROR( ret ) ) {
                                bError = ret;
                            } else {
                                nChanges += (ret & 1);
                                ret = SetBondsFromBnStructFlow( pBNS, at, num_atoms, bChangeFlow );
                                if ( IS_BNS_ERROR( ret ) ) {
                                    bError = ret;
                                } else
                                if ( ret >= 0 ) {
                                    nChanges += (ret & 1);
                                    bSuccess = 1;
                                } else {
                                    bError = ret;
                                }
                            }
                        }
                        /* typically 2*ret < nDots; 2*ret > nDots should not happen. Check later */
                        ret = RestoreBnStructFlow( pBNS, bChangeFlow & BNS_EF_CHNG_RSTR);
                        if ( IS_BNS_ERROR( ret ) ) {
                            bError = ret;
                        }
                    }
                    /* --- reinitialize to repeat the calculations --- */
                    ReInitBnStructAltPaths( pBNS );
                } else
                if ( nDots == 0 ) {
                    ret = bSetBondsAfterCheckOneBond( pBNS, fcd, nTestFlow, at, num_atoms, bChangeFlow );
                    if ( IS_BNS_ERROR( ret ) ) {
                        bError = ret;
                    } else {
                        nChanges += (ret & 1);
                    }
                }
                ret = bRestoreFlowAfterCheckOneBond( pBNS, fcd );
                if ( IS_BNS_ERROR( ret ) ) {
                    bError = ret;
                }
            }
        }
    }
    return bError? bError : nChanges;
}
/************************************************************************/
static void remove_alt_bond_marks(inp_ATOM *at, int num_atoms)
{
    int i, j, val;
    for ( i = 0; i < num_atoms; i++ ) {
        for ( val = at[i].valence, j = 0; j < val; j ++ ) {
            at[i].bond_type[j] &= BOND_TYPE_MASK;
        }
    }
}
/***************************************************************************************/
int SetForbiddenEdges( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int forbidden_mask )
{
    static U_CHAR el_number_O;
    static U_CHAR el_number_C;
    static U_CHAR el_number_N;

    int i, j, neigh, num_found;
    BNS_IEDGE iedge;
    /*S_CHAR    edge_forbidden_mask = BNS_EDGE_FORBIDDEN_MASK;*/
    S_CHAR    edge_forbidden_mask = forbidden_mask;

    pBNS->edge_forbidden_mask |= forbidden_mask;

    if ( !el_number_C ) {
        el_number_O = (U_CHAR)get_periodic_table_number( "O" );
        el_number_C = (U_CHAR)get_periodic_table_number( "C" );
        el_number_N = (U_CHAR)get_periodic_table_number( "N" );
    }
    
    num_found = 0;

    for ( i = 0; i < num_atoms; i ++ ) {
        /* acetyl */
        if ( at[i].el_number == el_number_C && 3 == at[i].valence &&
             4 == at[i].chem_bonds_valence ) {
            int num_O = 0;
            int bond_to_O_val = 0;
            int forbidden_bond_pos = -1;
            int forbidden_bond_val = -1;
            for ( j = 0; j < at[i].valence; j ++ ) {
                neigh = at[i].neighbor[j];
                if ( at[neigh].el_number == el_number_O &&
                     at[neigh].valence   == 1  ) {
                    num_O ++;
                    bond_to_O_val += (at[i].bond_type[j] & BOND_TYPE_MASK);
                } else {
                    forbidden_bond_pos = j;
                    forbidden_bond_val = (at[i].bond_type[j] & BOND_TYPE_MASK);
                }
            }
            if ( 2 == num_O && 3 == bond_to_O_val && 1 == forbidden_bond_val ) {
                iedge = pBNS->vert[i].iedge[forbidden_bond_pos];
                pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                num_found ++;
            }
        } else
        /* nitro */
        if ( at[i].el_number == el_number_N && 3 == at[i].valence &&
             (4 == at[i].chem_bonds_valence || 5 == at[i].chem_bonds_valence) ) {
            int num_O = 0;
            int bond_to_O_val = 0;
            int forbidden_bond_pos = -1;
            int forbidden_bond_val = -1;
            for ( j = 0; j < at[i].valence; j ++ ) {
                neigh = at[i].neighbor[j];
                if ( at[neigh].el_number == el_number_O &&
                     at[neigh].valence   == 1  ) {
                    num_O ++;
                    bond_to_O_val += (at[i].bond_type[j] & BOND_TYPE_MASK);
                } else {
                    forbidden_bond_pos = j;
                    forbidden_bond_val = (at[i].bond_type[j] & BOND_TYPE_MASK);

                }
            }
            if ( 2 == num_O && (3 == bond_to_O_val || 4 == bond_to_O_val) && 1 == forbidden_bond_val ) {
                iedge = pBNS->vert[i].iedge[forbidden_bond_pos];
                pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                num_found ++;
            }
        }
    }
#if ( REMOVE_ION_PAIRS_FIX_BONDS == 1 )
    num_found += fix_special_bonds( pBNS, at, num_atoms, edge_forbidden_mask );
#endif
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
    num_found += TempFix_NH_NH_Bonds( pBNS, at, num_atoms );
#endif
    return num_found;
}
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
/************************************************************************/
int TempFix_NH_NH_Bonds( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms )
{
    static U_CHAR el_number_N;
    int i, j, neigh, num_found;
    BNS_IEDGE iedge;
    S_CHAR    edge_forbidden_mask = BNS_EDGE_FORBIDDEN_TEMP;
    if ( !el_number_N ) {
        el_number_N = (U_CHAR)get_periodic_table_number( "N" );
    }
    for ( i = 0, num_found = 0; i < num_atoms; i ++ ) {
        /* -NH-NH- or -NH-NH3 */
        if ( at[i].el_number == el_number_N && at[i].valence < 3 && at[i].num_H &&
             3 == at[i].chem_bonds_valence + at[i].num_H   &&
             at[i].chem_bonds_valence == at[i].valence &&
             !at[i].charge && !at[i].radical ) {
            for ( j = 0; j < at[i].valence; j ++ ) {
                neigh = at[i].neighbor[j];
                if ( neigh < i &&
                     at[neigh].el_number == el_number_N && at[neigh].valence < 3 && at[neigh].num_H &&
                     3 == at[neigh].chem_bonds_valence + at[neigh].num_H   &&
                     at[neigh].chem_bonds_valence == at[neigh].valence &&
                     !at[neigh].charge && !at[neigh].radical) {
                    iedge = pBNS->vert[i].iedge[j];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_found ++;
                }
            }
        }
    }
    return num_found;
}
/************************************************************************/
int CorrectFixing_NH_NH_Bonds( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms )
{
    static U_CHAR el_number_N;
    int i, j, neigh, num_found;
    BNS_IEDGE iedge;
    S_CHAR    edge_forbidden_mask = BNS_EDGE_FORBIDDEN_TEMP;
    if ( !el_number_N ) {
        el_number_N = (U_CHAR)get_periodic_table_number( "N" );
    }
    for ( i = 0, num_found = 0; i < num_atoms; i ++ ) {
        /* -NH-NH- or -NH-NH3 */
        if ( at[i].el_number == el_number_N && at[i].valence < 3 ) {
            for ( j = 0; j < at[i].valence; j ++ ) {
                neigh = at[i].neighbor[j];
                if ( neigh < i &&
                     at[neigh].el_number == el_number_N && at[neigh].valence < 3 ) {
                    if ( BOND_TYPE_SINGLE != (at[i].bond_type[j] & BOND_TYPE_MASK) ) {
                        iedge = pBNS->vert[i].iedge[j];
                        if ( pBNS->edge[iedge].forbidden & edge_forbidden_mask ) {
                            pBNS->edge[iedge].forbidden &= ~edge_forbidden_mask;
                            num_found ++;
                        }
                    }
                }
            }
        }
    }
    return num_found;
}
#endif
/************************************************************************/
/* fixes bonds set by remove_ion_pairs() in strutil.c */
int fix_special_bonds( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int forbidden_mask )
{   
    int num_changes = 0;

    /*                           0 1 2  3  4 5 6  7  8  9                   8  9  */
#if( FIX_REM_ION_PAIRS_Si_BUG == 1 )
    static const char    el[] = "N;P;As;Sb;O;S;Se;Te;C;Si;";  /* 8 elements + C, Si */
#else
    static const char    el[] = "N;P;As;Sb;O;S;Se;Te;C;Si";   /* 8 elements + C, Si */
#endif
    static char    en[12];         /* same number: 8 elements */
    static int     ne=0;           /* will be 8 and 10 */

#define ELEM_N_FST  0
#define ELEM_N_LEN  4
#define ELEM_O_FST  4
#define ELEM_O_LEN  4
#define ELEM_S_FST  (ELEM_O_FST+1)
#define ELEM_S_LEN  (ELEM_O_LEN-1)
#define ELEM_C_FST  8
#define ELEM_C_LEN  2

#define MAX_NEIGH 6

    int i, k, n1, n2, n3, n4, i1, i2, i3, i4, bond_type;
    inp_ATOM *a;
    char elname[ATOM_EL_LEN];
    int j[3], m[3], num_O, k_O, num_N, num_OH, num_OM, num_X, num_other, k_N;

    BNS_IEDGE iedge;
    /*S_CHAR    edge_forbidden_mask = BNS_EDGE_FORBIDDEN_MASK;*/
    S_CHAR    edge_forbidden_mask = forbidden_mask;

    pBNS->edge_forbidden_mask |= edge_forbidden_mask;

    if ( !ne ) { /* one time initialization */
        const char *b, *e;
        int  len;
        for ( b = el; e = strchr( b, ';'); b = e+1 ) {
            len = e-b;
            memcpy( elname, b, len );
            elname[len] = '\0';
            en[ne++] = get_periodic_table_number( elname );
        }
        en[ne]   = '\0';
        en[ne+1] = '\0';
    }
    for ( i = 0, a = at; i < num_atoms; i ++, a++ ) {
        if ( !a->charge && !a->radical && 
             2 <= a->chem_bonds_valence + NUMH(a,0) - get_el_valence( a->el_number, 0, 0 ) &&
             0 == num_of_H( at, i ) &&
             2 == nNoMetalBondsValence(at, i) + NUMH(a,0) - get_el_valence( a->el_number, 0, 0 ) &&
             NULL != memchr( en+ELEM_N_FST, at[i].el_number, ELEM_N_LEN) ) {
            /* found N(V), no H */
            if ( 2 == nNoMetalNumBonds(at, i) ) {
                /* #N= */
                /* fix bonds: double and triple: =N# so that bonds cannot be changed by the normalization */
#if( FIX_N_V_METAL_BONDS_GPF == 1 )
                if ( 0 > (i1 = nNoMetalNeighIndex( at, i )) ||
                     0 > (i2 = nNoMetalOtherNeighIndex( at, i,
                                 n1 = a->neighbor[i1]/* non-metal neighbor #1 */ ) ) ) {
                    /*num_err ++; */ /* do not count would-be original InChI v.1 buffer overflow GPF */
                    continue; /* v1 bug: 2 bonds to metal yield i1 < 0 and/or i2 < 0 => bounds violation */
                }
#else
                i1 = nNoMetalNeighIndex( at, i );
                n1 = a->neighbor[i1]; /* non-metal neighbor #1 */
                i2 = nNoMetalOtherNeighIndex( at, i, n1 );
#endif
                n2 = a->neighbor[i2]; /* non-metal neighbor #2 */
                /* forbid all edges to non-metals */
                iedge = pBNS->vert[i].iedge[i1];
                pBNS->edge[iedge].forbidden |= edge_forbidden_mask; /* fix bond to neighbor #1 */
                iedge = pBNS->vert[i].iedge[i2];
                pBNS->edge[iedge].forbidden |= edge_forbidden_mask; /* fix bond to neighbor #1 */
                num_changes ++;  /* added 11-15-2005 */
                /*                                                      i n3 */
                /* forbid single bond edge beyond the neighboring =N- as in #N=N- */
                if ( (at[i].bond_type[i1] & BOND_TYPE_MASK) == BOND_TYPE_DOUBLE ) {
                    i3 = i1;
                    n3 = n1;
                } else {
                    i3 = i2;
                    n3 = n2;
                }
                if ( 0 == NUMH(at, n3) && 2 == nNoMetalNumBonds( at, n3 ) &&
                     3 == nNoMetalBondsValence( at, n3 ) &&
                     NULL != memchr( en+ELEM_N_FST, at[n3].el_number, ELEM_N_LEN) &&
                     0 <= (k = nNoMetalOtherNeighIndex( at, n3, i )) ) {
                    /* found =N- ; forbid the edge*/
                    iedge = pBNS->vert[n3].iedge[k];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_changes ++;
                }

            } else
            if ( 3 == nNoMetalNumBonds(at, i) &&
                /*         |                      */
                /* found  =N=                     */
                /* locate all non-metal neighbors */
                 0 <= (j[0] = nNoMetalNeighIndex( at, i )) &&
                 0 <= (j[1] = nNoMetalOtherNeighIndex( at, i, m[0] = a->neighbor[j[0]] )) &&
                 0 <= (j[2] = nNoMetalOtherNeighIndex2( at, i, m[0], m[1] = a->neighbor[j[1]] )) ) {
                /* count specific neighbors: N(V)=N, N(V)-N N(V)=O, and N(V)-N */
                /* if there is a single neighbor connected by a double bond, namely
                   N(V)=N and/or N(V)=O, then fix the bond(s).
                   If N(V)=O was fixed then do not fix another bond */
                m[2] = a->neighbor[j[2]];
                num_O = num_N = 0;
                for ( k= 0; k < 3; k ++ ) {
                    n1 = m[k];
                    i1 = j[k];
                    if ( NULL != memchr( en+ELEM_N_FST, at[n1].el_number, ELEM_N_LEN) ) {
                        k_N = k;
                        num_N ++;
                    } else
                    if ( NULL != memchr( en+ELEM_O_FST, at[n1].el_number, ELEM_O_LEN) &&
                         1 == nNoMetalNumBonds( at, n1 ) ) {
                        k_O = k;
                        num_O ++;
                    }
                }
                num_other = 0;
                if ( 1 == num_O && 0 == at[n1=m[k_O]].charge && 0 == at[n1].radical &&
                     BOND_TYPE_DOUBLE == (at[i].bond_type[i1=j[k_O]] & BOND_TYPE_MASK) ) {
                    /* fix bond to neighbor =O */
                    iedge = pBNS->vert[i].iedge[i1];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_changes ++;
                    num_other ++; /* indicator: double to a terminal O has been fixed */
                }
                if ( !num_other && num_O <= 1 &&
                     1 == num_N && 0 == at[n1=m[k_N]].charge && 0 == at[n1].radical &&
                     BOND_TYPE_DOUBLE == (at[i].bond_type[i1=j[k_N]] & BOND_TYPE_MASK ) ) {
                    /* fix bond to neighbor =N */
                    iedge = pBNS->vert[i].iedge[i1];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_changes ++;
                }
                         
            } else
            if ( 4 == nNoMetalNumBonds(at, i) ) {
                /*         |                      */
                /* found  -N=N-                   */
                /*         |                      */
                /* locate non-metal neighbor connected by a double bond;
                 * if it is =N- then fix the double bond and the single bond beyond the neighbor
                 */
                num_N = 0;
                num_other = 0;
                for ( i1 = 0; i1 < at[i].valence; i1 ++ ) {
                    if ( BOND_TYPE_DOUBLE == (at[i].bond_type[i1] & BOND_TYPE_MASK) && 
                         !is_el_a_metal(at[n1=(int)at[i].neighbor[i1]].el_number) &&
                         NULL != memchr( en+ELEM_N_FST, at[n1].el_number, ELEM_N_LEN) ) {
                        num_N ++;
                        n2 = n1;
                        i2 = i1;
                    }
                }
                if ( 1 == num_N && 0 == NUMH(at, n2) &&
                     2 == nNoMetalNumBonds( at, n2 ) &&
                     3 == nNoMetalBondsValence(at, n2) &&
                     0 <= ( i3 = nNoMetalOtherNeighIndex( at, n2, i ) ) &&
                     BOND_TYPE_SINGLE == (at[n2].bond_type[i3] & BOND_TYPE_MASK) ) {
                    /* fix the single bond beyond the N(V) neighbor N(V)=N- */
                    iedge = pBNS->vert[n2].iedge[i3];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_changes ++;
                    /* fix the double bond */
                    iedge = pBNS->vert[i].iedge[i2];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_changes ++;
                }
            }
        } else
        if ( !a->charge && !a->radical && 
             2 <= a->chem_bonds_valence + NUMH(a,0) - get_el_valence( a->el_number, 0, 0 ) &&
             0 == num_of_H( at, i ) &&
             2 == nNoMetalBondsValence(at, i) + NUMH(a,0) - get_el_valence( a->el_number, 0, 0 ) &&
             NULL != memchr( en+ELEM_S_FST, a->el_number, ELEM_S_LEN) &&
             3 == nNoMetalNumBonds( at, i ) ) {
            /* found S(IV), no H, one double bond, total 3 bonds */
            /*       OH
                    /
              in O=S   (X != O) fix single bond O-X  (type 1)
                    \
                     X
             
                     X
                    /
              in Z=S    (X, Y != OH) fix double bond Z=S (type 2)
                    \
                     Y
            */
            num_N  = 0; /* number of non-metal neighbors connected by a double bond */
            num_OH = 0; /* number of neighbors OH connected by a single bond */
            num_OM = 0; /* number of charged neighbors O connected by a single bond */
            num_O  = 0; /* number of neighbors =O connected by a double bond */
            num_other = 0;
            for ( i1 = 0; i1 < a->valence; i1 ++ ) {
                n1=(int)a->neighbor[i1];
                if ( is_el_a_metal(at[n1].el_number) ) {
                    continue;
                }
                bond_type = (a->bond_type[i1] & BOND_TYPE_MASK);
                if ( BOND_TYPE_DOUBLE == bond_type ) {
                    num_N ++;
                    n2 = n1;
                    i2 = i1;
                    if ( NULL != memchr( en+ELEM_O_FST, at[n1].el_number, ELEM_O_LEN) ) {
                        num_O ++;
                    }
                } else
                if ( BOND_TYPE_SINGLE == bond_type &&
                     1 == nNoMetalNumBonds( at, n1 ) &&
                     1 == nNoMetalBondsValence(at, n1 ) &&
                     NULL != memchr( en+ELEM_O_FST, at[n1].el_number, ELEM_O_LEN) ) {
                    if ( 0 == at[n1].charge ) {
                        num_OH ++;
                        n3 = n1;
                        i3 = i1;
                    } else {
                        num_OM ++;
                    }
                } else {
                    num_other ++;
                    n4 = n1;
                    i4 = i1;
                }
            }
            if ( 1 == num_N && 1 == num_O && 1 == num_OH + num_OM ) {
                if ( 1 == num_other ) {
                    /* type 1: fix the single bond S-X */
                    iedge = pBNS->vert[i].iedge[i4];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_changes ++;
                }
            } else
            if ( 1 == num_N && !num_OH && !num_OM ) {
                int    bFound = 0; /* flag */
                int    bDoNotFixAnyBond = 0; /* flag */
                /* avoid case N=S-NH or N=S-N(-); N = N, P, As, Sb */
                if ( NULL != memchr( en+ELEM_N_FST, at[n2].el_number, ELEM_N_LEN) ) {
                    U_CHAR el_number = at[n2].el_number;
                    for ( i1 = 0; i1 < a->valence; i1 ++ ) {
                        n1=(int)a->neighbor[i1];
                        bond_type = (a->bond_type[i1] & BOND_TYPE_MASK);
                        if ( BOND_TYPE_SINGLE == bond_type &&
                             (NUMH(at, n1) || -1 == at[n1].charge) &&
                             el_number == at[n1].el_number ) {
                            i3 = i1;
                            n3 = n1;
                            bFound ++;
                        }
                    }
                }
                /* exception: check if Z==X and they belong to the same ring system */
                for ( i1 = 0; i1 < a->valence; i1 ++ ) {
                    if ( i1 != i2 ) {
                        n1=(int)a->neighbor[i1]; 
                        if ( at[n2].el_number   == at[n1].el_number &&
                             at[n2].nRingSystem == at[n1].nRingSystem ) {
                            bDoNotFixAnyBond ++;
                        }
                    }
                }

                if ( bDoNotFixAnyBond ) {
                    ; /* do nothing */
                } else
                if ( bFound ) {
                    if ( 1 == bFound && 
                         0 <= ( i4 = nNoMetalOtherNeighIndex2( at, i, n2, n3) ) ) {
                        /* fix bond i4 */
                        iedge = pBNS->vert[i].iedge[i4];
                        pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                        num_changes ++;
                    }
                } else {
                    /* fix the double bond >S=X */
                    iedge = pBNS->vert[i].iedge[i2];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_changes ++;
                    /* -- test later --
                    if ( 2 == nNoMetalNumBonds( at, n2 ) &&
                         0 <= ( i3 = nNoMetalOtherNeighIndex( at, n2, i ) ) ) {
                        iedge = pBNS->vert[n2].iedge[i3];
                        pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                        num_changes ++;
                    }
                    -------------------*/
                }
            }
        } else
        if ( !a->charge && !a->radical && 
             4 <= a->chem_bonds_valence + NUMH(a,0) - get_el_valence( a->el_number, 0, 0 ) &&
             0 == num_of_H( at, i ) &&
             4 == nNoMetalBondsValence(at, i) + NUMH(a,0) - get_el_valence( a->el_number, 0, 0 ) &&
             NULL != memchr( en+ELEM_S_FST, a->el_number, ELEM_S_LEN) &&
             4 == nNoMetalNumBonds( at, i ) ) {
            /* found S(VI), no H, two double bonds or one triple bond */
            /*     O 
                   || 
              in O=S--Y-   (X, Y -- non-terminal) fix single bonds S-X, S-Y  (type 1)
                    \
                     X--
             
                   O  
                   || 
              in O=S--O(-)    (X -- non-terminal) fix single bond S-X (type 2)
                    \
                     X--

                   O  
                   || 
              in O=S--OH     (X -- non-terminal) fix single bond S-X (type 3)
                    \
                     X--

            */
            int iN[4];  /* indexes of non-terminal neighbors connected by a single bond */
            num_N  = 0; /* number of non-metal neighbors connected by a double bond */
            num_OH = 0; /* number of neighbors OH connected by a single bond */
            num_OM = 0; /* number of non-terminal neighbors connected by a single bond */
            num_O  = 0; /* number of neighbors =O connected by a double bond */
            num_X  = 0; /* number of terminal atom X != O connected by a single bond */
            num_other = 0;
            for ( i1 = 0; i1 < a->valence; i1 ++ ) {
                n1=(int)a->neighbor[i1];
                if ( is_el_a_metal(at[n1].el_number) ) {
                    continue;
                }
                bond_type = (a->bond_type[i1] & BOND_TYPE_MASK);
                if ( BOND_TYPE_DOUBLE == bond_type ) {
                    num_N ++;
                    if ( (0 == at[n1].charge
#if( S_VI_O_PLUS_METAL_FIX_BOND == 1 )
                          || 1 == at[n1].charge && 2 == at[n1].valence
#endif
                         ) && 0 == at[n1].radical &&
                         0 == num_of_H( at, n1 ) &&
                         NULL != memchr( en+ELEM_O_FST, at[n1].el_number, ELEM_O_LEN) &&
                         1 == nNoMetalNumBonds( at, n1 ) ) {

                        num_O ++;
                    }
                } else
                if ( BOND_TYPE_SINGLE == bond_type &&
                     1 == nNoMetalNumBonds( at, n1 ) &&
                     NULL != memchr( en+ELEM_O_FST, at[n1].el_number, ELEM_O_LEN) &&
                     1 >= num_of_H( at, n1 ) &&
                     1 == (( 0 == at[n1].charge) && 1==num_of_H( at, n1 ))
                         +((-1 == at[n1].charge) && 0==num_of_H( at, n1 ))  ) {

                    num_OH ++; /* -OH or -O(-) */

                } else
                if ( BOND_TYPE_SINGLE == bond_type &&
                     1 < nNoMetalNumBonds( at, n1 ) ) {

                    iN[num_OM ++] = i1; /* non-terminal neighbor connected by a single bond */

                } else
                if ( BOND_TYPE_SINGLE == bond_type &&
                     1 == nNoMetalNumBonds( at, n1 ) ) {

                    num_X ++; /* other terminal neighbor connected by a single bond */

                } else {
                    num_other ++;
                }
            }
            if ( num_N == num_O && 2 == num_O && 2 == num_OH + num_OM + num_X && 0 == num_other ) {
                for ( i2 = 0; i2 < num_OM; i2 ++ ) {
                    i1 = iN[i2];
                    /* fix bond i1 */
                    iedge = pBNS->vert[i].iedge[i1];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_changes ++;
                }
            }
        } else
        if ( !a->charge && !a->radical && 
             6 <= a->chem_bonds_valence + NUMH(a,0) - get_el_valence( a->el_number, 0, 0 ) &&
             0 == num_of_H( at, i ) &&
             6 == nNoMetalBondsValence(at, i) + NUMH(a,0) - get_el_valence( a->el_number, 0, 0 ) &&
             NULL != memchr( en+ELEM_S_FST, a->el_number, ELEM_S_LEN) &&
             5 == nNoMetalNumBonds( at, i ) ) {
            /* found S(VIII), no H, three double bonds or two triple bond */
            /*

                   O  
                   || 
              in O=S--Y--     (X, Y -- non-terminal) fix single bond S-X, S-Y (type 4)
                  //\
                 O   X--

              note: this structure is a mistakenly drawn structure

                   O                  O        
                   ||                 ||       
                 O=S--O--Y--  or    O=S--Y--
                    \                  \       
                     X--                O--X--    


            */
            int iN[5];  /* indexes of non-terminal neighbors connected by a single bond */
            num_N  = 0; /* number of non-metal neighbors connected by a double bond */
            num_OH = 0; /* number of neighbors OH connected by a single bond */
            num_OM = 0; /* number of non-terminal neighbors connected by a single bond */
            num_O  = 0; /* number of neighbors =O connected by a double bond */
            num_X  = 0; /* number of terminal atom X != O connected by a single bond */
            num_other = 0;
            for ( i1 = 0; i1 < a->valence; i1 ++ ) {
                n1=(int)a->neighbor[i1];
                if ( is_el_a_metal(at[n1].el_number) ) {
                    continue;
                }
                bond_type = (a->bond_type[i1] & BOND_TYPE_MASK);
                if ( BOND_TYPE_DOUBLE == bond_type ) {
                    num_N ++;
                    if ( (0 == at[n1].charge
#if( S_VI_O_PLUS_METAL_FIX_BOND == 1 )
                          || 1 == at[n1].charge && 2 == at[n1].valence
#endif
                          )                        
                         && 0 == at[n1].radical &&
                         0 == num_of_H( at, n1 ) &&
                         NULL != memchr( en+ELEM_O_FST, at[n1].el_number, ELEM_O_LEN) &&
                         1 == nNoMetalNumBonds( at, n1 ) ) {

                        num_O ++;
                    }
                } else
                if ( BOND_TYPE_SINGLE == bond_type &&
                     1 == nNoMetalNumBonds( at, n1 ) &&
                     NULL != memchr( en+ELEM_O_FST, at[n1].el_number, ELEM_O_LEN) &&
                     1 >= num_of_H( at, n1 ) &&
                     1 == (( 0 == at[n1].charge) && 1==num_of_H( at, n1 ))
                         +((-1 == at[n1].charge) && 0==num_of_H( at, n1 ))  ) {

                    num_OH ++; /* -OH or -O(-) */

                } else
                if ( BOND_TYPE_SINGLE == bond_type &&
                     1 < nNoMetalNumBonds( at, n1 ) ) {

                    iN[num_OM ++] = i1; /* non-terminal neighbor connected by a single bond */

                } else
                if ( BOND_TYPE_SINGLE == bond_type &&
                     1 == nNoMetalNumBonds( at, n1 ) ) {

                    num_X ++; /* other terminal neighbor connected by a single bond */

                } else {
                    num_other ++;
                }
            }
            if ( num_N == num_O && 3 == num_O && 2 == num_OH + num_OM + num_X && 0 == num_other ) {
                for ( i2 = 0; i2 < num_OM; i2 ++ ) {
                    i1 = iN[i2];
                    /* fix bond i1 */
                    iedge = pBNS->vert[i].iedge[i1];
                    pBNS->edge[iedge].forbidden |= edge_forbidden_mask;
                    num_changes ++;
                }
            }
        }
    }
    return num_changes;
}

#define ALL_NONMETAL_Z  0
/***************************************************************************************/
int is_Z_atom( U_CHAR el_number )
{
    typedef enum tag_Z_elnumber {
        el_C ,
        el_N ,
        el_P ,
        el_As,
        el_Sb,
        el_S ,
        el_Se,
        el_Te,
        el_Cl,
        el_Br,
        el_I ,
#if ( ALL_NONMETAL_Z == 1 )
        el_B ,
        el_O ,
        el_Si,
        el_Ge,
        el_F ,
        el_At,
#endif
        el_len
    } Z_ELNUMBER;
    static U_CHAR el_numb[el_len];
/*
    return is_el_a_metal( (int)el_number );
*/
    if ( !el_numb[el_C] ) {
        el_numb[el_C ] = (U_CHAR)get_periodic_table_number( "C"  );
        el_numb[el_N ] = (U_CHAR)get_periodic_table_number( "N"  );
        el_numb[el_P ] = (U_CHAR)get_periodic_table_number( "P"  );
        el_numb[el_As] = (U_CHAR)get_periodic_table_number( "As" );
        el_numb[el_Sb] = (U_CHAR)get_periodic_table_number( "Sb" );
        el_numb[el_S ] = (U_CHAR)get_periodic_table_number( "S"  );
        el_numb[el_Se] = (U_CHAR)get_periodic_table_number( "Se" );
        el_numb[el_Te] = (U_CHAR)get_periodic_table_number( "Te" );
        el_numb[el_Cl] = (U_CHAR)get_periodic_table_number( "Cl" );
        el_numb[el_Br] = (U_CHAR)get_periodic_table_number( "Br" );
        el_numb[el_I ] = (U_CHAR)get_periodic_table_number( "I"  );
#if ( ALL_NONMETAL_Z == 1 )
        el_numb[el_B ] = (U_CHAR)get_periodic_table_number( "B"  );
        el_numb[el_O ] = (U_CHAR)get_periodic_table_number( "O"  );
        el_numb[el_Si] = (U_CHAR)get_periodic_table_number( "Si" );
        el_numb[el_Ge] = (U_CHAR)get_periodic_table_number( "Ge" );
        el_numb[el_F ] = (U_CHAR)get_periodic_table_number( "F"  );
        el_numb[el_At] = (U_CHAR)get_periodic_table_number( "At" );
#endif
    }
    if ( memchr( el_numb, el_number, el_len ) ) {
        return 1;
    } 
    return 0;

}
/***************************************************************************************/
int IsZOX( inp_ATOM *atom, int at_x, int ord )
{  /* detect O==Z--X, O=O,S,Se,Te */
    static U_CHAR el_number_O  = 0;
    static U_CHAR el_number_S  = 0;
    static U_CHAR el_number_Se = 0;
    static U_CHAR el_number_Te = 0;
    inp_ATOM *at_Z = atom + atom[at_x].neighbor[ord];

    int i, neigh, num_O;

    if ( !el_number_O ) {
        el_number_O  = (U_CHAR)get_periodic_table_number( "O" );
        el_number_S  = (U_CHAR)get_periodic_table_number( "S" );
        el_number_Se = (U_CHAR)get_periodic_table_number( "Se" );
        el_number_Te = (U_CHAR)get_periodic_table_number( "Te" );
    }
    for ( i = 0, num_O = 0; i < at_Z->valence; i ++ ) {
        neigh = at_Z->neighbor[i];
        if ( neigh == at_x ) {
            continue;
        }
        if ( atom[neigh].valence == 1 &&
             atom[neigh].chem_bonds_valence == 2 &&
             atom[neigh].charge  == 0 &&
             atom[neigh].radical == 0 &&
              (atom[neigh].el_number == el_number_O  ||
               atom[neigh].el_number == el_number_S  ||
               atom[neigh].el_number == el_number_Se ||
               atom[neigh].el_number == el_number_Te ) ) {
            num_O ++;
        }
    }
    return num_O;
}
/***************************************************************************************/
int GetAtomChargeType( inp_ATOM *atom, int at_no, int nAtTypeTotals[], int *pMask, int bSubtract  )
{
    static U_CHAR el_number_C  = 0;
    static U_CHAR el_number_O  = 0;
    static U_CHAR el_number_S  = 0;
    static U_CHAR el_number_Se = 0;
    static U_CHAR el_number_Te = 0;
    static U_CHAR el_number_P  = 0;
    static U_CHAR el_number_N  = 0;
    static U_CHAR el_number_H  = 0;

    static U_CHAR el_number_F  = 0;
    static U_CHAR el_number_Cl = 0;
    static U_CHAR el_number_Br = 0;
    static U_CHAR el_number_I  = 0;
    
    inp_ATOM *at = atom + at_no;
#if( FIX_NORM_BUG_ADD_ION_PAIR == 1 )
    int i, neigh, mask, bit, type, num_z, num_m, num_o, delta = bSubtract > 0 ? -1 : 1; /* 0 or -2 => add, 1 or 2 => subtract */
    int bNoAdjIon = (bSubtract==0 || bSubtract==1);
#else
    int i, neigh, mask, bit, type, num_z, num_m, num_o, delta = bSubtract? -1 : 1;
#endif
    int bUnsatNHasTerminalO = 0;
    if ( !el_number_C ) {
        el_number_C  = (U_CHAR)get_periodic_table_number( "C" );
        el_number_O  = (U_CHAR)get_periodic_table_number( "O" );
        el_number_S  = (U_CHAR)get_periodic_table_number( "S" );
        el_number_Se = (U_CHAR)get_periodic_table_number( "Se" );
        el_number_Te = (U_CHAR)get_periodic_table_number( "Te" );
        el_number_P  = (U_CHAR)get_periodic_table_number( "P" );
        el_number_N  = (U_CHAR)get_periodic_table_number( "N" );
        el_number_H  = (U_CHAR)get_periodic_table_number( "H" );
        el_number_F  = (U_CHAR)get_periodic_table_number( "F" );
        el_number_Cl = (U_CHAR)get_periodic_table_number( "Cl" );
        el_number_Br = (U_CHAR)get_periodic_table_number( "Br" );
        el_number_I  = (U_CHAR)get_periodic_table_number( "I" );
    }

    type = ATT_NONE;
    mask = 0;
    if ( at->radical && at->radical != RADICAL_SINGLET ) {
        goto exit_function;
    }
    if ( is_el_a_metal( at->el_number ) ) {
        goto exit_function; /* metal */
    }
    if ( at->charge < -1 || at->charge > 1 ) {
        goto exit_function;
    }
    if ( !at->valence && at->charge == 1 && !at->num_H && !at->radical && at->el_number == el_number_H ) {
        /* a proton (#1) */
        type = ATT_PROTON;
        mask = ATBIT_Proton;
        goto count_mask_bits;
    }
    if ( !at->valence && at->charge == -1 && !at->num_H && !at->radical &&
         ( at->el_number == el_number_F  ||
           at->el_number == el_number_Cl ||
           at->el_number == el_number_Br ||
           at->el_number == el_number_I )
         ) {
        /* a halogen anion (#2) */
        type = ATT_HalAnion;
        mask = ATBIT_HalAnion;
        goto count_mask_bits;
    }
#if ( HAL_ACID_H_XCHG == 1 )
    /* halogen/chalcogen acid */
    if ( !at->valence && at->charge == 0 && 1 == at->num_H && !at->radical &&
         ( at->el_number == el_number_F  ||
           at->el_number == el_number_Cl ||
           at->el_number == el_number_Br ||
           at->el_number == el_number_I ) ||
         !at->valence && at->charge == 0 && 2 == at->num_H && !at->radical &&
         ( at->el_number == el_number_O  ||
           at->el_number == el_number_S ||
           at->el_number == el_number_Se ||
           at->el_number == el_number_Te )
         ) {
        /* a halogen/chalcogen acid (#3) */
        type = ATT_HalAcid;
        mask = ATBIT_HalAcid;
        goto count_mask_bits;
    }
#endif    
    if ( detect_unusual_el_valence( at->el_number, at->charge, at->radical,
                                        at->chem_bonds_valence, at->num_H,
                                        at->valence ) ) {
        goto exit_function; /* unusual valence state */
    }
    /* check neighbors */
    for ( i = 0, num_z = 0, num_m = 0, num_o = 0; i < at->valence; i ++ ) {
        neigh = at->neighbor[i];
#if( FIX_NORM_BUG_ADD_ION_PAIR == 1 )
        if ( atom[neigh].charge < -1 || atom[neigh].charge > 1 ) {
            goto exit_function; /* neighboring charge */
        }
        if ( atom[neigh].charge && at->charge ) {
            if ( bNoAdjIon ) {
                goto exit_function; /* neighboring charge */
            }
            type = ATT_NONE;
            mask = 0;
            goto count_mask_bits;
        }
#else
        if ( atom[neigh].charge < -1 || atom[neigh].charge > 1 || atom[neigh].charge && at->charge ) {
            goto exit_function; /* neighboring charge */
        }
#endif
        if ( detect_unusual_el_valence( atom[neigh].el_number, atom[neigh].charge, atom[neigh].radical,
                                            atom[neigh].chem_bonds_valence, atom[neigh].num_H,
                                            atom[neigh].valence ) ) {
            goto exit_function; /* neighbor in unusual valence state */
        }
        if ( is_Z_atom( atom[neigh].el_number ) ) {
            num_z ++;
        }
        if ( is_el_a_metal( atom[neigh].el_number ) ) {
            num_m ++;
        }
        num_o += (atom[neigh].el_number == el_number_O);
        if ( at->el_number == el_number_N && at->valence == 2 && !at->charge &&
             /*at->valence < at->chem_bonds_valence &&*/
             atom[neigh].valence == 1 && atom[neigh].chem_bonds_valence == 2 &&
             (atom[neigh].el_number == el_number_O  ||
              atom[neigh].el_number == el_number_S  ||
              atom[neigh].el_number == el_number_Se ||
              atom[neigh].el_number == el_number_Te  )) {
            bUnsatNHasTerminalO ++;
        }
    }
    /* O, S, Se, Te */
    if ( at->el_number == el_number_O  ||
         at->el_number == el_number_S  ||
         at->el_number == el_number_Se ||
         at->el_number == el_number_Te ) {
        if ( at->charge == 1 ) {
            if ( at->num_H ) { /* #4 */
                type  = ATT_O_PLUS;
                mask |= ATBIT_OH_Plus;
            } else { /* #5 */
                type  = ATT_O_PLUS;
                mask |= ATBIT_O_Plus;
            }
        } else
        if ( at->valence > 1 ) {
            goto exit_function;  /* not a terminal atom #C1 */
        } else
        if ( at->valence && !(num_z || num_o) ) {
            if ( num_m == at->valence ) {
                goto exit_function; /* #C2 */
            }
            goto count_mask_bits;  /* #C3 count charges, no donor or acceptor found */
        } else
        /* here at->neigh[0] is one of: O, or Z={C, N, P, As, Sb, S, Se, Te, Cl, Br, I} */
        if ( at->valence ) {
            neigh = at->neighbor[0]; /* Z or O only */
            if ( !atom[neigh].charge && atom[neigh].el_number == el_number_C &&
                 atom[neigh].chem_bonds_valence > atom[neigh].valence ) {
                /* =C-OH, #C-OH, =C-O(-), #C-O(-), -C=O, =C=O; O = O, S, Se, Te */
                type = ATT_ACIDIC_CO;
                if ( at->num_H == 1 ) {
                    mask |= (ATBIT_COH);            /* #6: =C-OH, #C-OH; O=O,S,Se,Te */
                    /*nAtTypeTotals[ATTOT_NUM_COH] ++;*/
                } else
                if ( at->charge == -1 ) {
                    mask |= (ATBIT_CO_Minus);       /* #7: =C-O(-), #C-O(-); O=O,S,Se,Te */
                    /*nAtTypeTotals[ATTOT_NUM_CO_Minus] ++;*/
                } else
                if ( !at->num_H && !at->charge ) {
                    mask |= (ATBIT_CO);             /* #8 -C=O, =C=O; O=O,S,Se,Te */
                    /*nAtTypeTotals[ATTOT_NUM_CO] ++;*/
                } else {
                    mask |= (ATBIT_Errors);
                    /*nAtTypeTotals[ATTOT_NUM_Errors] ++;*/
                }
            } else
            if ( !atom[neigh].charge &&
                  (atom[neigh].el_number == el_number_O  ||
                   atom[neigh].el_number == el_number_S  ||
                   atom[neigh].el_number == el_number_Se ||
                   atom[neigh].el_number == el_number_Te )  &&
                  atom[neigh].chem_bonds_valence == atom[neigh].valence ) {
                /* -O-OH, -O-O(-); O = O, S, Se, Te */
                type = ATT_OO;
                if ( at->num_H == 1 ) {
                    mask |= (ATBIT_OOH);            /* #9 -O-OH */
                    /*nAtTypeTotals[ATTOT_NUM_OOH] ++;*/
                } else
                if ( at->charge == -1 ) {
                    mask |= (ATBIT_OO_Minus);       /* #10 -O-O(-) */
                    /*nAtTypeTotals[ATTOT_NUM_OO_Minus] ++;*/
                } else {
                    mask |= (ATBIT_Errors);
                    /*nAtTypeTotals[ATTOT_NUM_Errors] ++;*/
                }
            } else
            if ( !atom[neigh].charge &&
                 atom[neigh].chem_bonds_valence == atom[neigh].valence &&
                 atom[neigh].el_number == el_number_C &&
                 at->el_number         != el_number_O ) {
                /* >C-S(-), >C-SH; S = S, Se, Te  */
                type = ATT_ACIDIC_S;
                if ( at->num_H == 1 ) {
                    mask |= (ATBIT_CSH);            /* #11: >C-SH, >CH-SH, -CH2-SH; S = S, Se, Te */
                    /*nAtTypeTotals[ATTOT_NUM_CSH] ++;*/
                } else
                if ( at->charge == -1 ) {
                    mask |= (ATBIT_CS_Minus);       /* #12: >C-S(-), >CH-S(-), -CH2-S(-); S = S, Se, Te  */
                    /*nAtTypeTotals[ATTOT_NUM_CS_Minus] ++;*/
                } else {
                    mask |= (ATBIT_Errors);
                    /*nAtTypeTotals[ATTOT_NUM_Errors] ++;*/
                }
            } else
            if ( atom[neigh].el_number == el_number_N &&
                 atom[neigh].valence == 2 && (!atom[neigh].num_H || atom[neigh].num_H == 1 && atom[neigh].charge == 1) ) {
                 /* N or N(-) or NH(+) neighbor */
                type = ATT_NO; /* single bond only */
                if ( at->num_H == 1 ) {
                    mask |= (ATBIT_NOH);            /* #13: =N-OH, =NH(+)-OH, #N(+)-OH, -N(-)-OH; O = O, S, Se, Te */
                    /*nAtTypeTotals[ATTOT_NUM_NOH] ++;*/
                } else
                if ( at->charge == -1 ) {
                    mask |= (ATBIT_NO_Minus);       /* #14: =N-O(-); O = O, S, Se, Te */
                    /*nAtTypeTotals[ATTOT_NUM_NO_Minus] ++;*/
                } else
                if ( atom[neigh].charge == 1 || atom[neigh].charge == 0 ) {
                    mask |= (ATBIT_NO);             /* #15: =N(+)=O, -NH(+)=O -N=O */
                    /*nAtTypeTotals[ATTOT_NUM_NO] ++;*/
                } else {
                    mask |= (ATBIT_Errors);
                    /*nAtTypeTotals[ATTOT_NUM_Errors] ++;*/
                }
            } else
            if ( atom[neigh].el_number == el_number_N  ) {
                type = ATT_N_O; /* #16: single bond only */
                if ( at->num_H == 1 ) {
                    mask |= (ATBIT_N_OH);            /* #16: -NH-OH, >N-OH or >N(+)<OH; O = O, S, Se, Te */
                    /*nAtTypeTotals[ATTOT_NUM_NOH] ++;*/
                } else
                if ( at->charge == -1 ) {
                    mask |= (ATBIT_N_O_Minus);       /* #17: -NH-O(-), >N-O(-); O = O, S, Se, Te */
                    /*nAtTypeTotals[ATTOT_NUM_NO_Minus] ++;*/
                } else
                if ( atom[neigh].charge == 1 ) {
                    mask |= (ATBIT_N_O);             /* #18:  >N(+)=O */
                    /*nAtTypeTotals[ATTOT_NUM_NO] ++;*/
                } else {
                    mask |= (ATBIT_Errors);
                    /*nAtTypeTotals[ATTOT_NUM_Errors] ++;*/
                }
            } else
            if ( atom[neigh].el_number != el_number_C && atom[neigh].el_number != el_number_O &&
                 !is_el_a_metal( atom[neigh].el_number ) &&
                 atom[neigh].chem_bonds_valence > atom[neigh].valence ) {
                 /* =Z-OH, #Z-OH, =Z-O(-), #Z-O(-), -Z=O, =Z=O;
                    =Z(+)-OH, #Z(+)-OH, =Z-O(-), #Z-O(-), -Z(+)=O, =Z(+)=O; O = O, S, Se, Te */
                 /* neigh = Z\{N,C} = P, As, Sb, S, Se, Te, Cl, Br, I */
                if ( at->chem_bonds_valence == 1 && IsZOX( atom, at_no, 0 ) ) {
                    type = ATT_ZOO;
                    if ( at->num_H == 1 ) {
                        mask |= (ATBIT_ZOOH);            /* 18: O=Z-OH; O=O,S,Se,Te; Z may have charge */
                        /*nAtTypeTotals[ATTOT_NUM_ZOOH] ++;*/
                    } else
                    if ( at->charge == -1 ) {
                        mask |= (ATBIT_ZOO_Minus);       /* 19: O=Z-O(-); O = O, S, Se, Te */
                        /*nAtTypeTotals[ATTOT_NUM_ZOO_Minus] ++;*/
                    } else {
                        mask |= (ATBIT_Errors);
                        /*nAtTypeTotals[ATTOT_NUM_Errors] ++;*/
                    }
                } else {
                    type = ATT_OTHER_ZO;
                    if ( at->num_H == 1 ) {
                        mask |= (ATBIT_ZOH);            /* 20: =Z-OH, #Z-OH; O=O,S,Se,Te; Z may have charge */
                        /*nAtTypeTotals[ATTOT_NUM_ZOH] ++;*/
                    } else
                    if ( at->charge == -1 ) {
                        mask |= (ATBIT_ZO_Minus);       /* 21: =Z-O(-), #Z-O(-); O = O, S, Se, Te */
                        /*nAtTypeTotals[ATTOT_NUM_ZO_Minus] ++;*/
                    } else
                    if ( at->num_H == 0 ) {
                        mask |= (ATBIT_ZO);             /* 22: -Z=O, =Z=O; O=O,S,Se,Te; Z may have charge */
                        /*nAtTypeTotals[ATTOT_NUM_ZO] ++;*/
                    } else {
                        mask |= (ATBIT_Errors);
                        /*nAtTypeTotals[ATTOT_NUM_Errors] ++;*/
                    }
                }
            } else
            if ( at->charge == -1 && !is_el_a_metal( atom[neigh].el_number ) ) {
                /* >Z-O(-); O=O,S,Se,Te */
                type = ATT_OTHER_NEG_O;
                mask |= (ATBIT_O_Minus);            /* 23: -Z-O(-); O=O,S,Se,Te */
                /*nAtTypeTotals[ATTOT_NUM_ZO_Minus] ++;*/
            }
        } else
        if ( at->charge == -1 && at->num_H == 1 ) {
            type = ATT_OH_MINUS;
            mask |= (ATBIT_O_Minus);            /* 25: HO(-); O=O,S,Se,Te */
        }
    } else
    /* P, N, neutral valence = 3 (not 5) */
    if ( (at->el_number == el_number_N  ||
          at->el_number == el_number_P) &&
          0 <= at->valence && at->valence <= 3 &&
          at->chem_bonds_valence + at->num_H == 3 + at->charge ) {
        if ( at->valence && !(num_z /*|| num_o == at->valence*/) ) {
            if ( num_m == at->valence ) {
                goto exit_function;
            }
            goto count_mask_bits; /* N(III), N(-)(II), N(+)(IV) and same P that have only oxygen neighbors are ignored here */
        }
        type = (at->el_number == el_number_N)? ATT_ATOM_N : ATT_ATOM_P;
        switch ( at->charge ) {
        case -1:
            if (at->el_number == el_number_N) {
                mask |= (ATBIT_N_Minus); /* 26: -NH(-), =N(-), >N(-) */
                
                if ( at->num_H )
                    mask |= (ATBIT_NP_H); /* 27: -NH(-) */
#if ( FIX_NP_MINUS_BUG == 1 )
                else
                if ( at->valence == 1 && at->chem_bonds_valence >= 2 && (at->bond_type[0] & BOND_MARK_ALL) )
                    type |= ATT_NP_MINUS_V23; /* =N(-) created by normalization 2010-03-11 DT */
#endif
                
            }
            /*nAtTypeTotals[ATTOT_NUM_N_Minus] += (at->el_number == el_number_N);*/
            break;
        case 0:
            if ( at->num_H ) {
                mask |= (ATBIT_NP_H);          /* 28: -NH2, =NH, >NH */
                /*nAtTypeTotals[ATTOT_NUM_NP_H] ++;*/
            } else {
                if ( bUnsatNHasTerminalO == 1 ) {
                    mask |= (ATBIT_ON);            /* 29: -N=O,-N=OH(+) only, not =N-OH */
                } else {
                    mask |= (ATBIT_NP);            /* 30: -P=O,-P=OH(+), >N- =N- (incl. =N-OH) , #N */
                    /*nAtTypeTotals[ATTOT_NUM_NP] ++;*/
                }
            }
            break; /* ignore neutral N or P */
        case 1:
            if ( at->num_H ) {
                mask |= (ATBIT_NP_Proton);       /* 31: NH4(+), -NH3(+), =NH2(+), >NH2(+), =NH(+)-, >NH(+)-, #NH(+) */
                /*nAtTypeTotals[ATTOT_NUM_NP_Proton] ++;*/
            } else
            if ( at->chem_bonds_valence > at->valence ) {
                mask |= (ATBIT_NP_Plus);         /* =N(+)=, #N(+)-, =N(+)< */
                /*nAtTypeTotals[ATTOT_NUM_NP_Plus] ++;*/
            } else {
                type = 0; /* 32: ignore onium cations >N(+)< */
            }
            break;
        default:
            mask |= (1 << ATTOT_NUM_Errors);
            /*nAtTypeTotals[ATTOT_NUM_Errors] ++;*/
            break;
        }
    }
count_mask_bits:
    if ( nAtTypeTotals ) {
        if ( mask && !(mask & (ATBIT_Errors)) ) {
            for ( i = 0, bit = 1; i < ATTOT_ARRAY_LEN; i ++, bit <<= 1 ) {
                if ( bit & mask ) {
                    nAtTypeTotals[i] += delta;
                }
            }
        }
        /* count charges */
        if ( at->charge ) {
            nAtTypeTotals[ATTOT_TOT_CHARGE]  += delta * at->charge;
            nAtTypeTotals[ATTOT_NUM_CHARGES] += delta;
        }
    }
    if ( pMask ) {
        *pMask = mask;
    }
exit_function:
    if ( mask & (ATBIT_Errors) ) {
        type = 0;
        if ( nAtTypeTotals ) {
            nAtTypeTotals[ATTOT_NUM_Errors] += 1;
        }
    }
    return type;
}
/***************************************************************************************/
int SimpleRemoveHplusNPO( inp_ATOM *at, int num_atoms, int nAtTypeTotals[], T_GROUP_INFO *t_group_info )
{
    int i, mask, type, num_removed;
    for ( i = 0, num_removed = 0; i < num_atoms; i ++ ) {
        if ( (PR_SIMPLE_TYP & (type = GetAtomChargeType( at, i, NULL, &mask, 0 )) ) &&
             (PR_SIMPLE_MSK & mask ) ) {
#if ( bRELEASE_VERSION == 0 )
            if ( at[i].charge != 1 || at[i].num_H == 0 ) {
                return -1;  /* program error */
            }
#endif
            type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 1 ); /* subtract at[i] */
            at[i].charge = 0;
            AddOrRemoveExplOrImplH( -1, at, num_atoms, (AT_NUMB)i, t_group_info );
            /*at[i].num_H --;*/
            num_removed ++;
#if( FIX_NORM_BUG_ADD_ION_PAIR == 1 )
            type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 0 ); /* add changed at[i] */
#else
            type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 1 ); /* bug: subtract instead of add */
#endif
            /*
            if ( nAtTypeTotals ) {
                nAtTypeTotals[ATTOT_NUM_NP_Proton] --;
                if ( at[i].num_H ) {
                    nAtTypeTotals[ATTOT_NUM_NP_H] ++;
                } else {
                    nAtTypeTotals[ATTOT_NUM_NP] ++;
                }
                nAtTypeTotals[ATTOT_TOT_CHARGE] --;
                nAtTypeTotals[ATTOT_NUM_CHARGES] --;
            }
            */
        }
    }
    return num_removed;
}

/***************************************************************************************/
int bIsAtomTypeHard( inp_ATOM *at, int endpoint, int nType, int nMask, int nCharge )
{
    int        mask;    
    if ( (nType & GetAtomChargeType( at, endpoint, NULL, &mask, 0 )) && (mask & nMask)
#if( OPPOSITE_CHARGE_IN_CGROUP == 0 )
         && ( at[endpoint].charge == nCharge || !at[endpoint].charge )
#endif
    ) {
        return 1;
    }
    return 0;

}
/***************************************************************************************/
int bIsHDonorAccAtomType( inp_ATOM *at, int endpoint, int *cSubType )
{
    if ( bIsAtomTypeHard( at, endpoint, PR_HARD_TYP_H, PR_HARD_MSK_H, 0 ) ) {
        /* obtain donor/acceptor info */
        int neutral_valence = at[endpoint].chem_bonds_valence + at[endpoint].num_H - at[endpoint].charge;
        if ( neutral_valence != 2 /* O, S, Se, Te */ &&
             neutral_valence != 3 /* N, P */          ) {
            return -1; /* wrong endpoint neutral valence */
        } else {
            int edge_flow = at[endpoint].num_H;
            int num_bonds = at[endpoint].valence;
            int edge_cap  = neutral_valence - num_bonds; /* does not allow to reduce -NH3(+) to #N or -OH(+)- to -O- */
            edge_flow = inchi_min( edge_flow, edge_cap);
            /* what this means: */
            if ( edge_cap ) {
                if ( edge_cap > edge_flow )
                    *cSubType |= SALT_ACCEPTOR;
                if ( edge_flow )
                    *cSubType |= SALT_DONOR_H;
                return 4;
            }
        }
    }
    return -1;
}

/***************************************************************************************/
int bIsNegAtomType( inp_ATOM *at, int endpoint, int *cSubType )
{
    int sub_type = 0;
    if ( bIsAtomTypeHard( at, endpoint, PR_HARD_TYP_NEG, PR_HARD_MSK_NEG, -1 ) ) {
        /* obtain donor/acceptor info */
        int neutral_valence = at[endpoint].chem_bonds_valence + at[endpoint].num_H - at[endpoint].charge;
        if ( neutral_valence != 2 /* O, S, Se, Te */ &&
             neutral_valence != 3 /* N, P */          ) {
            return -1; /* wrong endpoint neutral valence */
        } else {
            int edge_flow = (at[endpoint].charge == -1);
            int num_bonds = at[endpoint].valence;
            int edge_cap  = neutral_valence - num_bonds - at[endpoint].num_H; /* does not allow to reduce -NH3(+) to #N or -OH(+)- to -O- */
            edge_flow = inchi_min( edge_flow, edge_cap);
            /* what this means: */
            if ( edge_cap ) {
                if ( edge_cap > edge_flow )
                    sub_type |= SALT_ACCEPTOR;
                if ( edge_flow ) {
                    sub_type |= SALT_DONOR_Neg;
                }
                if ( sub_type ) {
                    *cSubType |= sub_type;
                    return 4;
                }
            }
        }
    }
    return -1;
}
/***************************************************************************************/
int bIsHardRemHCandidate(  inp_ATOM *at, int i, int *cSubType )
{
    int ret1, ret2, ret;
    int sub_type = 0;
    ret1 = bIsHDonorAccAtomType( at, i, &sub_type );
    ret2 = bIsNegAtomType( at, i, &sub_type );
    ret = inchi_max(ret1, ret2);
    if ( ret > 0 && sub_type ) {
        *cSubType |= sub_type;
        return ret;
    }
    return -1;
}

/***************************************************************************************/
int CreateCGroupInBnStruct( inp_ATOM *at, int num_atoms,
                            BN_STRUCT *pBNS, int nType, int nMask, int nCharge )
{
    int         k, c_point, cg, centerpoint, fictpoint, type, ret = 0;
    int         num_cg       = 1;
    int         num_edges    = pBNS->num_edges;
    int         num_vertices = pBNS->num_vertices; /* new c-group bns-ID */
    BNS_VERTEX *vert_ficpoint, *ver_ficpont_prev;  /* fictitious vertex describing charge c-group */
    BNS_VERTEX *vertex_cpoint;
    BNS_EDGE   *edge;      /* edge between that vertex and the tautomeric c_point */
    int        mask, num_CPoints;

    /* Debug: check overflow */
    if ( num_vertices + num_cg >= pBNS->max_vertices ) {
        return BNS_VERT_EDGE_OVFL;
    }
    /* count new c-group edges */
    for ( c_point = 0, num_CPoints = 0; c_point < num_atoms; c_point ++ ) {
        if ( (nType & GetAtomChargeType( at, c_point, NULL, &mask, 0 )) && (mask & nMask)
#if( OPPOSITE_CHARGE_IN_CGROUP == 0 )
             && ( at[c_point].charge == nCharge || !at[c_point].charge )
#endif
           ) {
            num_CPoints ++;
        }
    }
    if ( !num_CPoints ) {
        return 0;
    }



    /* clear the new vertex */
    memset( pBNS->vert+num_vertices, 0, 1*sizeof(pBNS->vert[0]) );
    /* *old* Make sure the last t-group has the largest t-group ID:
       this is necessary to correctly add new edges and vertices for testing augmenting paths
    */
    /**************************************/
    /* initialize new fictitious vertex   */
    /* representing c-point group         */
    /**************************************/
    ver_ficpont_prev = pBNS->vert+num_vertices - 1;
    
    for ( cg = 0; cg < num_cg; cg ++, ver_ficpont_prev = vert_ficpoint ) {
        /*
          vert_ficpoint-1 is the last vertex;
          vert_ficpoint   is the being added vertex
          Note: nGroupNumber are not contiguous
        */
        vert_ficpoint                = pBNS->vert+num_vertices + cg;
        vert_ficpoint->iedge         = ver_ficpont_prev->iedge + ver_ficpont_prev->max_adj_edges;
        vert_ficpoint->max_adj_edges = num_CPoints+BNS_ADD_EDGES;
        vert_ficpoint->num_adj_edges = 0;
        vert_ficpoint->st_edge.flow  = vert_ficpoint->st_edge.flow0  = 0;
        vert_ficpoint->st_edge.cap   = vert_ficpoint->st_edge.cap0   = 0;
        vert_ficpoint->type          = BNS_VERT_TYPE_C_GROUP | ((nCharge<0)?BNS_VERT_TYPE_C_NEGATIVE:0);
    }

    /************************************************/
    /* connect c-points to the fictitious vertices  */
    /* representing c-point groups; set caps, flows */
    /************************************************/
    cg = 1;
    for ( c_point = 0; c_point < num_atoms; c_point ++ ) {
        if ( (nType & (type=GetAtomChargeType( at, c_point, NULL, &mask, 0 ))) && (mask & nMask)
#if( OPPOSITE_CHARGE_IN_CGROUP == 0 )
             && ( at[c_point].charge == nCharge || !at[c_point].charge)
#endif
        );
        else
            continue;
        fictpoint = cg + num_vertices - 1; /* c-group vertex index */
        vert_ficpoint = pBNS->vert + fictpoint; /* c-group vertex */
        vertex_cpoint = pBNS->vert + c_point;   /* c_point vertex */
        /* Debug: check overflow */
        if ( fictpoint >= pBNS->max_vertices ||
             num_edges >= pBNS->max_edges    ||
             vert_ficpoint->num_adj_edges >= vert_ficpoint->max_adj_edges ||
             vertex_cpoint->num_adj_edges >= vertex_cpoint->max_adj_edges ) {
            ret = BNS_VERT_EDGE_OVFL;
            break;
        }
        vertex_cpoint->type |= BNS_VERT_TYPE_C_POINT;
        if ( (KNOWN_ACIDIC_TYPE & type) && nCharge < 0 ) {
            vertex_cpoint->type |= pBNS->type_TACN;
        }
#if( FIX_CPOINT_BOND_CAP != 1 )  /* { */
        /* set capacity = 1 to the edges from the c_point to the centerpoint(s)     */
        /* if their current capacity is zero                                        */
        /* the centerpoint is any adjacent atom that is adjacent to a multiple bond */
        for ( k = 0; k < vertex_cpoint->num_adj_edges; k ++ ) {
            int iedge = vertex_cpoint->iedge[k];
            if ( !pBNS->edge[iedge].cap ) {
                /* single bond, possibly between c_point and centerpoint */
                centerpoint = (pBNS->edge[iedge].neighbor12 ^ c_point);
                if ( centerpoint < pBNS->num_atoms &&
                     pBNS->vert[centerpoint].st_edge.cap >= 1 ) {
                    int bond_type = (at[c_point].bond_type[k] & BOND_TYPE_MASK);
                    if ( bond_type == BOND_TAUTOM ||
                         bond_type == BOND_ALTERN ||
                         bond_type == BOND_SINGLE ) {
                        pBNS->edge[iedge].cap = 1;
                    }
                }
            }
        }
#endif /* } FIX_CPOINT_BOND_CAP */
        /* create a new edge connecting c_point to the new fictitious c-group vertex vert_ficpoint */
        edge = pBNS->edge + num_edges;
        edge->cap       = 1;
        edge->flow      = 0;
        edge->pass      = 0;
#if( RESET_EDGE_FORBIDDEN_MASK == 1 )
        edge->forbidden &= pBNS->edge_forbidden_mask; /* remove previous temporary ban */
#endif
        /* nCharge = +1: mark edge to c-point having no (+)-moveable charge with flow=1 */
        /* nCharge = -1: mark edge to c-point having -1 moveable charge with flow=1 */
        if ( nCharge==1 && at[c_point].charge != 1 || nCharge==-1 && at[c_point].charge == -1 )
        /*if ( !CHARGED_CPOINT(at,c_point) )*/
        {
            /* increment new edge flow, update st_edges of the adjacent vertices */
            edge->flow ++;
            /* increment c-group vertex st-flow & cap */
            vert_ficpoint->st_edge.flow ++;
            vert_ficpoint->st_edge.cap ++;
            /* increment c-point vertex st-flow & cap */
            vertex_cpoint->st_edge.flow ++;
            vertex_cpoint->st_edge.cap ++;
        }
#if( FIX_CPOINT_BOND_CAP == 1 ) /* { */
        /* set capacity = 1 to the edges from the c_point to the centerpoint(s)     */
        /* if their current capacity is zero                                        */
        /* the centerpoint is any adjacent atom that is adjacent to a multiple bond */
        for ( k = 0; k < vertex_cpoint->num_adj_edges; k ++ ) {
            int iedge = vertex_cpoint->iedge[k];
            VertexFlow  nNewCap = vertex_cpoint->st_edge.cap;
            centerpoint = (pBNS->edge[iedge].neighbor12 ^ c_point);
            if ( !pBNS->edge[iedge].cap ) {
                /* single bond, possibly between c_point and centerpoint */
                if ( centerpoint < pBNS->num_atoms &&
                     pBNS->vert[centerpoint].st_edge.cap >= 1 ) {
                    nNewCap = inchi_min( pBNS->vert[centerpoint].st_edge.cap, nNewCap );
                    nNewCap = inchi_min( nNewCap, MAX_BOND_EDGE_CAP );
                    pBNS->edge[iedge].cap = nNewCap;
                }
            }
#if( FIX_CPOINT_BOND_CAP2 == 1 ) /* multiple bond */
            else
            if ( centerpoint < pBNS->num_atoms &&
                 edge->flow && pBNS->edge[iedge].cap < MAX_BOND_EDGE_CAP ) {
                pBNS->edge[iedge].cap ++;
            }
#endif
        }
#endif  /* } FIX_CPOINT_BOND_CAP */
        /* connect edge to c_point and fictpoint and increment the counters of neighbors and edges */
        edge->neighbor1    = c_point; /* the smallest out of v1=endopoint and v2=num_vertices */
        edge->neighbor12   = c_point ^ fictpoint; /* v1 ^ v2 */
        vertex_cpoint->iedge[vertex_cpoint->num_adj_edges] = num_edges;
        vert_ficpoint->iedge[vert_ficpoint->num_adj_edges] = num_edges ++;
        edge->neigh_ord[0] = vertex_cpoint->num_adj_edges ++;
        edge->neigh_ord[1] = vert_ficpoint->num_adj_edges ++;
        edge->cap0  = edge->cap;
        edge->flow0 = edge->flow;

    }
    ret = pBNS->num_vertices; /* new c-group atom number */
    pBNS->num_edges     = num_edges;
    pBNS->num_vertices += num_cg;
    pBNS->num_c_groups += num_cg;

    return ret;
}

/*********************************************************************************/
int CreateTGroupInBnStruct( inp_ATOM *at, int num_atoms,
                            BN_STRUCT *pBNS, int nType, int nMask )
{
    int ret = 0;
    /* ret = ReInitBnStruct( pBNS ); */
    int         k, endpoint, tg, centerpoint, fictpoint;
    int         num_tg       = 1;
    int         num_edges    = pBNS->num_edges;
    int         num_vertices = pBNS->num_vertices;
    BNS_VERTEX *vert_ficpoint, *ver_ficpont_prev;  /* fictitious vertex describing t-group */
    BNS_VERTEX *vert_endpoint;
    BNS_EDGE   *edge;      /* edge between that vertex and the tautomeric endpoint */
     int        mask, num_endpoints, neutral_valence, edge_flow, edge_cap, num_bonds;

    /* Debug: check overflow */
    if ( num_vertices + num_tg >= pBNS->max_vertices ) {
        return BNS_VERT_EDGE_OVFL;
    }
    /* count new t-group edges */
    for ( endpoint = 0, num_endpoints = 0; endpoint < num_atoms; endpoint ++ ) {
        if ( (nType & GetAtomChargeType( at, endpoint, NULL, &mask, 0 )) && (mask & nMask)
           ) {
            num_endpoints ++;
        }
    }
    if ( !num_endpoints ) {
        return 0;
    }


    /* since t-group IDs may be not contiguous, clear all vertices that will be added.
       all-zeroes-vertex will be ignored by the BNS
    */
    memset( pBNS->vert+num_vertices, 0, num_tg*sizeof(pBNS->vert[0]) );
    /* *old* Make sure the last t-group has the largest t-group ID:
       this is necessary to correctly add new edges and vertices for testing augmenting paths
    */
    /**************************************/
    /* initialize new fictitious vertex   */
    /* representing t-point group         */
    /**************************************/
    ver_ficpont_prev = pBNS->vert+num_vertices - 1;

    for ( tg = 0; tg < num_tg; tg ++, ver_ficpont_prev = vert_ficpoint ) {
        /*
          vert_ficpoint-1 is the last vertex;
          vert_ficpoint   is the vertex that is being added
          Note: nGroupNumber are not contiguous
        */
        vert_ficpoint                = pBNS->vert+num_vertices + tg;
        vert_ficpoint->iedge         = ver_ficpont_prev->iedge + ver_ficpont_prev->max_adj_edges;
        vert_ficpoint->max_adj_edges = num_endpoints+BNS_ADD_EDGES+BNS_ADD_SUPER_TGROUP;
        vert_ficpoint->num_adj_edges = 0;
        vert_ficpoint->st_edge.flow  = vert_ficpoint->st_edge.flow0  = 0;
        vert_ficpoint->st_edge.cap   = vert_ficpoint->st_edge.cap0   = 0;
        vert_ficpoint->type         |= BNS_VERT_TYPE_TGROUP;
    }
    tg = 1;
    for ( endpoint = 0; endpoint < num_atoms; endpoint ++ ) {
        if ( (nType & GetAtomChargeType( at, endpoint, NULL, &mask, 0 )) && (mask & nMask));
        else
            continue;
        fictpoint = tg + num_vertices - 1;
        vert_ficpoint = pBNS->vert + fictpoint;
        vert_endpoint = pBNS->vert + endpoint;
        /* Debug: check overflow */
        if ( fictpoint >= pBNS->max_vertices ||
             num_edges >= pBNS->max_edges    ||
             vert_ficpoint->num_adj_edges >= vert_ficpoint->max_adj_edges ||
             vert_endpoint->num_adj_edges >= vert_endpoint->max_adj_edges ) {
            ret = BNS_VERT_EDGE_OVFL;
            break;
        }
        /* obtain donor/acceptor info */
        neutral_valence = at[endpoint].chem_bonds_valence + at[endpoint].num_H - at[endpoint].charge;
        if ( neutral_valence != 2 /* O, S, Se, Te */ &&
             neutral_valence != 3 /* N, P */          ) {
            ret = BNS_PROGRAM_ERR; /* wrong endpoint neutral valence */
            break;
        }
        edge_flow = at[endpoint].num_H;
        num_bonds = at[endpoint].valence;
        edge_cap  = neutral_valence - num_bonds; /* does not allow to reduce -NH3(+) to #N or -OH(+)- to -O- */
        if ( 3 == neutral_valence /* N or P */ && 1 < num_bonds ) {
            edge_cap ++; /* allow -NH2(+)- => -N=, >NH(+)- => >N- */
        }
        edge_flow = inchi_min( edge_flow, edge_cap);
        /*
        if ( !nGetEndpointInfo( at, endpoint, &eif ) ) {
            ret = BNS_BOND_ERR;
            break;
        }
        */
        vert_endpoint->type |= BNS_VERT_TYPE_ENDPOINT;
        /* create a new edge connecting endpoint to the new fictitious t-group vertex vert_ficpoint */
        edge = pBNS->edge + num_edges;
        edge->cap       = edge_cap;
        edge->flow      = edge_flow;
        edge->pass      = 0;
#if( RESET_EDGE_FORBIDDEN_MASK == 1 )
        edge->forbidden &= pBNS->edge_forbidden_mask;
#endif
        /* adjust st_flow and st_cap of the adjacent vertices */
        /* adjust t-group vertex st-flow & cap */
        vert_ficpoint->st_edge.flow += edge->flow;
        vert_ficpoint->st_edge.cap  += edge->flow;
        /* adjust endpoint vertex st-flow & cap */
        vert_endpoint->st_edge.flow += edge->flow;
        vert_endpoint->st_edge.cap  += edge->flow;

        /* adjust edge cap & flow according to the number of H and number of bonds */
        for ( k = 0; k < vert_endpoint->num_adj_edges; k ++ ) {
            int iedge = vert_endpoint->iedge[k];
            VertexFlow  nNewCap = vert_endpoint->st_edge.cap;
            if ( !pBNS->edge[iedge].cap ) {
                /* single bond, possibly between endpoint and centerpoint */
                centerpoint = (pBNS->edge[iedge].neighbor12 ^ endpoint);
                if ( centerpoint < pBNS->num_atoms &&
                     pBNS->vert[centerpoint].st_edge.cap >= 1 ) {
                    nNewCap = inchi_min( pBNS->vert[centerpoint].st_edge.cap, nNewCap );
                    nNewCap = inchi_min( nNewCap, MAX_BOND_EDGE_CAP );
                    pBNS->edge[iedge].cap = nNewCap;
                }
            }
        }

        /* connect edge to endpoint and fictpoint and increment the counters of neighbors and edges */
        edge->neighbor1    = endpoint; /* the smallest out of v1=endopoint and v2=num_vertices */
        edge->neighbor12   = endpoint ^ fictpoint; /* v1 ^ v2 */
        vert_endpoint->iedge[vert_endpoint->num_adj_edges] = num_edges;
        vert_ficpoint->iedge[vert_ficpoint->num_adj_edges] = num_edges ++;
        edge->neigh_ord[0] = vert_endpoint->num_adj_edges ++;
        edge->neigh_ord[1] = vert_ficpoint->num_adj_edges ++;
        edge->cap0  = edge->cap;
        edge->flow0 = edge->flow;
    }

    ret = pBNS->num_vertices; /* new t-group atom number */
    pBNS->num_edges     = num_edges;
    pBNS->num_vertices += num_tg;
    pBNS->num_t_groups += num_tg;

    return ret;
}
/*********************************************************************************/
int RemoveLastGroupFromBnStruct( inp_ATOM *at, int num_atoms, int tg, BN_STRUCT *pBNS )
{
    int ret = 0;
    /* ret = ReInitBnStruct( pBNS ); */
    int         k, endpoint, /*centerpoint, fictpoint,*/ iedge;
    int         num_edges    = pBNS->num_edges;
    int         num_vertices = pBNS->num_vertices;
    BNS_VERTEX *vert_ficpoint /*, *ver_ficpont_prev*/;  /* fictitious vertex describing t-group */
    BNS_VERTEX *vert_endpoint;
    BNS_EDGE   *edge;      /* edge between that vertex and the tautomeric endpoint */
    /*int        mask, num_endpoints, neutral_valence, edge_flow, edge_cap, num_bonds;*/
    int        is_t_group = 0, is_c_group = 0;

    /* Debug: check overflow */
    if ( pBNS->num_added_atoms + pBNS->num_c_groups + pBNS->num_t_groups + num_atoms >= pBNS->max_vertices ) {
        return BNS_VERT_EDGE_OVFL;
    }
    if ( tg + 1 != num_vertices ) {
        return BNS_VERT_EDGE_OVFL;
    }
    vert_ficpoint = pBNS->vert + tg;
    if ( vert_ficpoint->type & BNS_VERT_TYPE_TGROUP ) {
        is_t_group = 1;
    }
    if ( vert_ficpoint->type & BNS_VERT_TYPE_C_GROUP ) {
        is_c_group = 1;
        if ( vert_ficpoint->type & BNS_VERT_TYPE_C_NEGATIVE )
            is_c_group = 2;
    }
    for ( k = vert_ficpoint->num_adj_edges-1; 0 <= k; k -- ) {
        iedge         = vert_ficpoint->iedge[k];
        if ( iedge + 1 != num_edges ) {
            return BNS_VERT_EDGE_OVFL;
        }
        edge          = pBNS->edge + iedge;
        endpoint      = edge->neighbor12 ^ tg;
        vert_endpoint = pBNS->vert + endpoint;
        /* adjust st_flow, st_cap */
        vert_endpoint->st_edge.cap0  =
        vert_endpoint->st_edge.cap   -= edge->flow;
        vert_endpoint->st_edge.flow0 =
        vert_endpoint->st_edge.flow  -= edge->flow;
        if ( pBNS->type_TACN && (vert_endpoint->type & pBNS->type_TACN) == pBNS->type_TACN ) {
            vert_endpoint->type ^= pBNS->type_TACN;
        }
        if ( is_t_group ) {
            vert_endpoint->type ^= (vert_ficpoint->type & BNS_VERT_TYPE_ENDPOINT);
        }
        if ( is_c_group ) {
            vert_endpoint->type ^= (vert_ficpoint->type & BNS_VERT_TYPE_C_POINT);
        }
        /* remove edge */
        if ( edge->neigh_ord[0]+1 != vert_endpoint->num_adj_edges ) {
            return BNS_VERT_EDGE_OVFL;
        }
        vert_endpoint->num_adj_edges --;
        memset( edge, 0, sizeof(*edge) );
        num_edges --;
        if ( 1 == is_t_group && endpoint < num_atoms ) {
            at->endpoint = 0;
        }
        if ( 1 == is_c_group && endpoint < num_atoms ) {
            at->c_point = 0;
        }
    }
    memset( vert_ficpoint, 0, sizeof(*vert_ficpoint) );
    num_vertices --;

    pBNS->num_edges     = num_edges;
    pBNS->num_vertices  = num_vertices;
    if ( is_t_group )
        pBNS->num_t_groups --;
    if ( is_c_group )
        pBNS->num_c_groups --;

    return ret;
}
/******************************************************************************************/
int SetInitCapFlowToCurrent( BN_STRUCT *pBNS )
{
    int       i, j;
    BNS_EDGE *pEdge=NULL;
    for ( i = 0; i < pBNS->num_vertices; i ++ ) {
        pBNS->vert[i].st_edge.flow0 = pBNS->vert[i].st_edge.flow;
        pBNS->vert[i].st_edge.cap0 = pBNS->vert[i].st_edge.cap;
        for ( j = 0; j < pBNS->vert[i].num_adj_edges; j ++ ) {
            pEdge = pBNS->edge + pBNS->vert[i].iedge[j];
            pEdge->cap0  = pEdge->cap;
            pEdge->flow0 = pEdge->flow;
        }
    }
    return 0;
}
/******************************************************************************************/
int ArTypMask[] = { 
    AR_SIMPLE_TYP1, AR_SIMPLE_MSK1,
    AR_SIMPLE_TYP2, AR_SIMPLE_MSK2,
    AR_SIMPLE_TYP3, AR_SIMPLE_MSK3,
    0, 0 };
/******************************************************************************************/
int SimpleRemoveAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2remove )
{
    int i, j, max_j=-1, mask, type, num_removed;
    int num[AR_SIMPLE_STEPS+1], num_tot;

    for ( j = 0; ArTypMask[2*j]; j ++ ) {
        num[max_j = j] = 0;
    }

    for ( i = 0; i < num_atoms; i ++ ) {
        if ( !at[i].charge && at[i].num_H && (type = GetAtomChargeType( at, i, NULL, &mask, 0 )) ) {
            for ( j = 0; j <= max_j; j ++ ) {
                if ( (type & ArTypMask[2*j]) && (mask && ArTypMask[2*j+1]) ) {
                    num[j]  ++;
                    break;
                }
            }
        }
    }
    for ( j = 0, num_tot = 0; j <= max_j; j ++ ) {
        if ( (num_tot += num[j]) >= num2remove ) {
            max_j = j;
            break;
        }
    }
    if ( !num_tot ) {
        return 0;
    }
    for ( i = 0, num_removed = 0; i < num_atoms && num_removed < num2remove; i ++ ) {
        if ( !at[i].charge && at[i].num_H && (type = GetAtomChargeType( at, i, NULL, &mask, 0 )) ) {
            for ( j = 0; j <= max_j; j ++ ) {
                if ( num[j] && (type & ArTypMask[2*j]) && (mask && ArTypMask[2*j+1]) ) {
                    type = GetAtomChargeType( at, i, pAATG->nAtTypeTotals, &mask, 1 ); /* subtract at[i] */
                    num[j]  --;
                    at[i].charge --;
                    AddOrRemoveExplOrImplH( -1, at, num_atoms, (AT_NUMB)i, pAATG->t_group_info );
                    /*at[i].num_H  --;*/
                    num_removed ++;
                    type = GetAtomChargeType( at, i, pAATG->nAtTypeTotals, &mask, 0 ); /* add changed at[i] */
                    break;
                }
            }
        }
    }
    /*
    pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE]  -= num_removed;
    pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] += num_removed;
    */
    return num_removed;
}
/******************************************************************************************/
int bHasAcidicHydrogen( inp_ATOM *at, int i )
{
    int bFound = 0, j, type, mask;
    if ( !at[i].charge && at[i].num_H && (type = GetAtomChargeType( at, i, NULL, &mask, 0 )) ) {
        for ( j = 0; ArTypMask[2*j]; j ++ ) {
            if ( (type & ArTypMask[2*j]) && (mask & ArTypMask[2*j+1]) ) {
                bFound ++;
                break;
            }
        }
    }
    return bFound;
}
/******************************************************************************************/
int bHasOtherExchangableH ( inp_ATOM *at, int i )
{
    int bFound = 0, type, mask;
    if ( at[i].num_H && (type = GetAtomChargeType( at, i, NULL, &mask, 0 )) ) {
            if ( (type & ATT_ATOM_N) && (mask & ATBIT_NP_H) ) {
                bFound ++;
            }
    }
    return bFound;
}
/******************************************************************************************/
int AaTypMask[] = { 
    AA_SIMPLE_TYP1, AA_SIMPLE_MSK1,
#if ( FIX_NP_MINUS_BUG == 1 )
    AA_SIMPLE_TYP4, AA_SIMPLE_MSK4,   /* should not follow 0,0 pair */
#endif
    AA_SIMPLE_TYP2, AA_SIMPLE_MSK2,
    AA_SIMPLE_TYP3, AA_SIMPLE_MSK3,
    0, 0 };
/******************************************************************************************/
int SimpleAddAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2add )
{
    int i, j, max_j=-1, mask, type, num_added;
    int num[AR_SIMPLE_STEPS+1], num_tot;

    for ( j = 0; AaTypMask[2*j]; j ++ ) {
        num[max_j = j] = 0;
    }

    for ( i = 0; i < num_atoms; i ++ ) {
        if ( at[i].charge==-1 &&  (type = GetAtomChargeType( at, i, NULL, &mask, 0 )) ) {
            for ( j = 0; j <= max_j; j ++ ) {
                if ( (type & AaTypMask[2*j]) && (mask && AaTypMask[2*j+1]) ) {
                    num[j]  ++;
                    break;
                }
            }
        }
    }
    for ( j = 0, num_tot = 0; j <= max_j; j ++ ) {
        if ( (num_tot += num[j]) >= num2add ) {
            max_j = j;
            break;
        }
    }
    if ( !num_tot ) {
        return 0;
    }
    for ( i = 0, num_added = 0; i < num_atoms && num_added < num2add; i ++ ) {
        if ( at[i].charge==-1 && (type = GetAtomChargeType( at, i, NULL, &mask, 0 )) ) {
            for ( j = 0; j <= max_j; j ++ ) {
                if ( num[j] && (type & AaTypMask[2*j]) && (mask && AaTypMask[2*j+1]) ) {
                    type = GetAtomChargeType( at, i, pAATG->nAtTypeTotals, &mask, 1 ); /* subtract at[i] */
                    num[j]  --;
                    at[i].charge ++;
                    AddOrRemoveExplOrImplH( 1, at, num_atoms, (AT_NUMB)i, pAATG->t_group_info );
                    /*at[i].num_H  ++;*/
                    num_added ++;
                    type = GetAtomChargeType( at, i, pAATG->nAtTypeTotals, &mask, 0 ); /* add changed at[i] */
                    break;
                }
            }
        }
    }
    /*
    pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE]  += num_added;
    pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] -= num_added;
    */
    return num_added;
}
/******************************************************************************************/
int bHasAcidicMinus( inp_ATOM *at, int i )
{
    int bFound = 0, j, type, mask;
    if ( at[i].charge==-1 &&  (type = GetAtomChargeType( at, i, NULL, &mask, 0 )) ) {
        for ( j = 0; AaTypMask[2*j]; j ++ ) {
            if ( (type & AaTypMask[2*j]) && (mask & AaTypMask[2*j+1]) ) {
                bFound  ++;
                break;
            }
        }
    }
    return bFound;
}
/******************************************************************************************
Create 2 tautomeric groups: (1) for O on -C=O, (2) for the rest of the atoms.
Pull H from (2) to (1); remove later 
*******************************************************************************************/
int HardRemoveAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2remove, int *nNumCanceledCharges, BN_STRUCT *pBNS, BN_DATA *pBD )
{
    int cg_Plus      = 0;
    int cg_Minus     = 0;
    int tg_H_Other   = 0;
    int tg_H_Acid    = 0;

    int ret = 0, ret2;
    int nDelta, nNumChanges = 0, nNumMoved2AcidH = 0, nNumNeutralized = 0, nPrevNumCharges;

    int nPosCharges, nPosCharges2;
    int nNegCharges, nNegCharges2;
    /*
    int nNumNP_H, nNumNP_H2;
    int nNumOS_H, nNumOS_H2;
    */

    nPosCharges =  (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] + pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    nNegCharges =  (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] - pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    /*
    nNumNP_H    =  pAATG->nAtTypeTotals[ATTOT_NUM_NP_H] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_NP_Proton];
    nNumOS_H    =  pAATG->nAtTypeTotals[ATTOT_NUM_COH] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_CSH] + 
                   pAATG->nAtTypeTotals[ATTOT_NUM_ZOH];
    */

    /* prevent free exchange H <-> (-) */
    pBNS->type_CN   = (BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE);
    pBNS->type_T    = BNS_VERT_TYPE_TGROUP;
    pBNS->type_TACN = BNS_VERT_TYPE_ACID;
    /* create (+) charge group */
    cg_Plus = CreateCGroupInBnStruct( at, num_atoms, pBNS, AR_HARD_TYP_POS, AR_HARD_MSK_POS, 1 );
    /* create (-) charge group */
    /*
    if ( nAtTypeTotals[ATTOT_NUM_CO_Minus] + 
         nAtTypeTotals[ATTOT_NUM_CS_Minus] + 
         nAtTypeTotals[ATTOT_NUM_ZO_Minus] + 
         nAtTypeTotals[ATTOT_NUM_N_Minus] )
    */
    cg_Minus = CreateCGroupInBnStruct( at, num_atoms, pBNS, AR_HARD_TYP_NEG, AR_HARD_MSK_NEG, -1 );

    pBNS->type_CN   = (BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE);
    pBNS->type_T    = BNS_VERT_TYPE_TGROUP;
    pBNS->type_TACN = BNS_VERT_TYPE_ACID;
    
    /* create tautomeric group for non-acidic or negatively charged acidic O  */
    tg_H_Other = CreateTGroupInBnStruct( at, num_atoms, pBNS, AR_HARD_TYP_HN, AR_HARD_MSK_HN );

    /* create tautomeric group for possibly acidic O */
    tg_H_Acid = CreateTGroupInBnStruct( at, num_atoms, pBNS, AR_HARD_TYP_HA, AR_HARD_MSK_HA );
    if ( tg_H_Other >= num_atoms && tg_H_Acid >= num_atoms ) {
        /* find alt path to remove one proton */
        do {
            /* remove a proton */
            nPrevNumCharges = pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES];
            ret = bExistsAltPath( pBNS, pBD, pAATG, at, num_atoms,
                            tg_H_Other /*nVertDoubleBond*/, tg_H_Acid /*nVertSingleBond*/, ALT_PATH_MODE_REM_PROTON );
            if ( IS_BNS_ERROR( ret ) ) {
                return ret;
            }
            if ( ret & 1 ) {
                nDelta       = (ret & ~3) >> 2;
                nNumChanges += (0 != (ret & 2));
                if ( nDelta ) {
                    /* radical pair has disappeared */
                    ; /* goto quick_exit;*/
                }
                nNumMoved2AcidH ++;
                if ( nPrevNumCharges > pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] + 1 ) {
                    nNumNeutralized += (nPrevNumCharges - (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] - 1))/2;
                }
            }

        } while ( (ret & 1) && nNumMoved2AcidH < num2remove );

        /* neutralize: remove ion pairs like >N(+)=-O(-) => >N-=O */
        if ( (nNumMoved2AcidH /*|| bCancelChargesAlways*/) && cg_Minus >= num_atoms && cg_Plus >= num_atoms &&
             pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] > abs(pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE]) ) {
            do {
                nPrevNumCharges = pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES];
                ret = bExistsAltPath( pBNS, pBD, pAATG, at, num_atoms,
                                cg_Minus /*nVertDoubleBond*/, cg_Plus /*nVertSingleBond*/, ALT_PATH_MODE_REM_PROTON );
                if ( IS_BNS_ERROR( ret ) ) {
                    return ret;
                }
                if ( ret & 1 ) {
                    nDelta       = (ret & ~3) >> 2;
                    nNumChanges += (0 != (ret & 2));
                    if ( nDelta ) {
                        /* radical pair has disappeared */
                        ; /* goto quick_exit;*/
                    }
                    if ( nPrevNumCharges > pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] ) {
                        nNumNeutralized += (nPrevNumCharges - pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES])/2;
                    }
                }
            } while ( ret & 1 );
        }
    }

    ret = 0;
    if ( tg_H_Acid >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, tg_H_Acid, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
    if ( tg_H_Other >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, tg_H_Other, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
    if ( cg_Minus >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, cg_Minus, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
    if ( cg_Plus >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, cg_Plus, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
    
    pBNS->type_CN   = 0;
    pBNS->type_T    = 0;
    pBNS->type_TACN = 0;
    
    if ( ret ) {
        return ret;
    }
    if ( pAATG->nAtTypeTotals[ATTOT_NUM_CO_Minus] + pAATG->nAtTypeTotals[ATTOT_NUM_ZO_Minus] &&
         pAATG->nAtTypeTotals[ATTOT_NUM_N_Minus] ) {
    }
    
    nPosCharges2 = (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] + pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    nNegCharges2 = (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] - pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    /*
    nNumNP_H2    = pAATG->nAtTypeTotals[ATTOT_NUM_NP_H] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_NP_Proton];
    nNumOS_H2    = pAATG->nAtTypeTotals[ATTOT_NUM_COH] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_CSH] + 
                   pAATG->nAtTypeTotals[ATTOT_NUM_ZOH];
    */
    if ( (nPosCharges - nNegCharges) - (nPosCharges2 - nNegCharges2) != 0 ) {
        return BNS_PROGRAM_ERR;
    }

    if ( nNumCanceledCharges ) {
#if( FIX_CANCEL_CHARGE_COUNT_BUG == 1 )
        *nNumCanceledCharges += 2*nNumNeutralized;
#else
        *nNumCanceledCharges = 2*nNumNeutralized;
#endif
    }
    
    return nNumMoved2AcidH;
}
/******************************************************************************************/
int HardAddAcidicProtons( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2add, int *nNumCanceledCharges, BN_STRUCT *pBNS, BN_DATA *pBD )
{
    int cg_Plus        = 0;
    int cg_Minus_CO    = 0;
    int cg_Minus_Other = 0;
    int tg_H           = 0;

    int ret = 0, ret2;
    int nDelta, nNumChanges = 0, nNumMoved2AcidMinus = 0, nNumNeutralized = 0, nPrevNumCharges;

    int nPosCharges, nPosCharges2;
    int nNegCharges, nNegCharges2;
    /*
    int nNumNP_H, nNumNP_H2;
    int nNumOS_H, nNumOS_H2;
    */
    nPosCharges =  (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] + pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    nNegCharges =  (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] - pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    /*
    nNumNP_H    =  pAATG->nAtTypeTotals[ATTOT_NUM_NP_H] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_NP_Proton];
    nNumOS_H    =  pAATG->nAtTypeTotals[ATTOT_NUM_COH] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_CSH] + 
                   pAATG->nAtTypeTotals[ATTOT_NUM_ZOH];
    */
    /* prevent free exchange H <-> (-) */
    pBNS->type_CN   = (BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE);
    pBNS->type_T    = BNS_VERT_TYPE_TGROUP;
    pBNS->type_TACN = BNS_VERT_TYPE_ACID;
    /* create (+) charge group */
    cg_Plus = CreateCGroupInBnStruct( at, num_atoms, pBNS, AA_HARD_TYP_POS, AA_HARD_MSK_POS, 1 );
    /* create (-) charge group */
    /*
    if ( nAtTypeTotals[ATTOT_NUM_CO_Minus] + 
         nAtTypeTotals[ATTOT_NUM_CS_Minus] + 
         nAtTypeTotals[ATTOT_NUM_ZO_Minus] + 
         nAtTypeTotals[ATTOT_NUM_N_Minus] )
    */
    cg_Minus_CO    = CreateCGroupInBnStruct( at, num_atoms, pBNS, AA_HARD_TYP_CO, AA_HARD_MSK_CO, -1 );

    cg_Minus_Other = CreateCGroupInBnStruct( at, num_atoms, pBNS, AA_HARD_TYP_NEG, AA_HARD_MSK_NEG, -1 );

    pBNS->type_CN   = (BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE);
    pBNS->type_T    = BNS_VERT_TYPE_TGROUP;
    pBNS->type_TACN = BNS_VERT_TYPE_ACID;
    
    /* create tautomeric group for all H  */
    tg_H = CreateTGroupInBnStruct( at, num_atoms, pBNS, AA_HARD_TYP_H, AA_HARD_MSK_H );

    
    if ( cg_Minus_Other >= num_atoms && cg_Minus_CO >= num_atoms ) {
        /* find alt path to remove one proton */
        do {
            /* add a proton */
            nPrevNumCharges = pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES];
            ret = bExistsAltPath( pBNS, pBD, pAATG, at, num_atoms,
                            cg_Minus_Other /*nVertDoubleBond*/, cg_Minus_CO /*nVertSingleBond*/, ALT_PATH_MODE_REM_PROTON );
            if ( IS_BNS_ERROR( ret ) ) {
                return ret;
            }
            if ( ret & 1 ) {
                nDelta       = (ret & ~3) >> 2;
                nNumChanges += (0 != (ret & 2));
                if ( nDelta ) {
                    /* radical pair has disappeared */
                    ; /* goto quick_exit;*/
                }
                nNumMoved2AcidMinus ++;
                if ( nPrevNumCharges > pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] + 1 ) {
                    nNumNeutralized += (nPrevNumCharges - (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] - 1))/2;
                }
            }

        } while ( (ret & 1) && nNumMoved2AcidMinus < num2add );

        /* neutralize: remove ion pairs like >N(+)=-O(-) => >N-=O */
        if ( (nNumMoved2AcidMinus /*|| bCancelChargesAlways*/) && cg_Minus_Other >= num_atoms && cg_Plus >= num_atoms &&
             pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] > abs(pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE]) ) {
            do {
                nPrevNumCharges = pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES];
                ret = bExistsAltPath( pBNS, pBD, pAATG, at, num_atoms,
                                cg_Minus_Other /*nVertDoubleBond*/, cg_Plus /*nVertSingleBond*/, ALT_PATH_MODE_REM_PROTON );
                if ( IS_BNS_ERROR( ret ) ) {
                    return ret;
                }
                if ( ret & 1 ) {
                    nDelta       = (ret & ~3) >> 2;
                    nNumChanges += (0 != (ret & 2));
                    if ( nDelta ) {
                        /* radical pair has disappeared */
                        ; /* goto quick_exit;*/
                    }
                    if ( nPrevNumCharges > pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] ) {
                        nNumNeutralized += (nPrevNumCharges - pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES])/2;
                    }
                }
            } while ( ret & 1 );
        }
    }

    ret = 0;
    if ( tg_H >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, tg_H, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
    if ( cg_Minus_Other >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, cg_Minus_Other, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
    if ( cg_Minus_CO >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, cg_Minus_CO, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
    if ( cg_Plus >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, cg_Plus, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
    
    pBNS->type_CN   = 0;
    pBNS->type_T    = 0;
    pBNS->type_TACN = 0;

    if ( ret ) {
        return ret;
    }
    if ( pAATG->nAtTypeTotals[ATTOT_NUM_CO_Minus] + pAATG->nAtTypeTotals[ATTOT_NUM_ZO_Minus] &&
         pAATG->nAtTypeTotals[ATTOT_NUM_N_Minus] ) {
    }
    
    nPosCharges2 = (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] + pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    nNegCharges2 = (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] - pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    /*
    nNumNP_H2    = pAATG->nAtTypeTotals[ATTOT_NUM_NP_H] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_NP_Proton];
    nNumOS_H2    = pAATG->nAtTypeTotals[ATTOT_NUM_COH] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_CSH] + 
                   pAATG->nAtTypeTotals[ATTOT_NUM_ZOH];
    */
    if ( (nPosCharges - nNegCharges) - (nPosCharges2 - nNegCharges2) != 0 ) {
        return BNS_PROGRAM_ERR;
    }

    if ( nNumCanceledCharges ) {
#if( FIX_CANCEL_CHARGE_COUNT_BUG == 1 )
        *nNumCanceledCharges += 2*nNumNeutralized;
#else
        *nNumCanceledCharges = 2*nNumNeutralized;
#endif
    }
    
    return nNumMoved2AcidMinus;
}
/******************************************************************************************/
/* examples include removal of H from tautomeric O that belongs to the same t-group as N: */
/* >N(+)=-N=-OH =(taut.)=> >N(+)=-NH-=O =(+charge move)=> >N-=NH(+)-=O => >N-=N-=O + H(+) */
/******************************************************************************************/
int HardRemoveHplusNP( inp_ATOM *at, int num_atoms, int bCancelChargesAlways, int *nNumCanceledCharges,
                       BN_AATG *pAATG, BN_STRUCT *pBNS, BN_DATA *pBD  )
{

    int cg_Plus      = 0;
    int cg_Minus     = 0;
    int tg_H         = 0;
#if ( MOVE_PPLUS_TO_REMOVE_PROTONS == 1 )
    int cg_PlusP     = 0;
#endif
#if ( FIX_REM_PROTON_COUNT_BUG == 1 )
    int nPrevRemovedProtons, nCurrRemovedProtons;
#endif
    int ret = 0, ret2;
    int nDelta, nNumChanges = 0, nNumRemovedProtons = 0, nNumNeutralized = 0, nPrevNumCharges;

    int nPosCharges, nPosCharges2;
    int nNegCharges, nNegCharges2;
    /*
    int nNumNP_H, nNumNP_H2;
    int nNumOS_H, nNumOS_H2;
    */

    nPosCharges =  (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] + pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    nNegCharges =  (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] - pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    /*
    nNumNP_H    =  pAATG->nAtTypeTotals[ATTOT_NUM_NP_H] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_NP_Proton];
    nNumOS_H    =  pAATG->nAtTypeTotals[ATTOT_NUM_COH] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_CSH] + 
                   pAATG->nAtTypeTotals[ATTOT_NUM_ZOH];
    */
    /* prevent free exchange H <-> (-) */
    pBNS->type_CN   = (BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE);
    pBNS->type_T    = BNS_VERT_TYPE_TGROUP;
    pBNS->type_TACN = BNS_VERT_TYPE_ACID;
    /* create (+) charge group */
    cg_Plus = CreateCGroupInBnStruct( at, num_atoms, pBNS, PR_HARD_TYP_POS, PR_HARD_MSK_POS, 1 );
    /* create (-) charge group */
    /*
    if ( nAtTypeTotals[ATTOT_NUM_CO_Minus] + 
         nAtTypeTotals[ATTOT_NUM_CS_Minus] + 
         nAtTypeTotals[ATTOT_NUM_ZO_Minus] + 
         nAtTypeTotals[ATTOT_NUM_N_Minus] )
    */
#if ( MOVE_PPLUS_TO_REMOVE_PROTONS == 1 )
    cg_PlusP = CreateCGroupInBnStruct( at, num_atoms, pBNS, PR_HARD_TYP_POSP, PR_HARD_MSK_POS, 1 );
#endif
    cg_Minus = CreateCGroupInBnStruct( at, num_atoms, pBNS, PR_HARD_TYP_NEG, PR_HARD_MSK_NEG, -1 );
    
    /* create single tautomeric group */
    tg_H = CreateTGroupInBnStruct( at, num_atoms, pBNS, PR_HARD_TYP_H, PR_HARD_MSK_H );

    if ( tg_H >= num_atoms && cg_Plus >= num_atoms ) {

#if( FIX_N_MINUS_NORN_BUG == 1 )
        /* neutralize: remove ion pairs like >N(+)=-O(-) => >N-=O; >N(+)=-NH(-) => >N-=NH */
        if ( (nNumRemovedProtons || bCancelChargesAlways) && cg_Minus >= num_atoms && cg_Plus >= num_atoms &&
             pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] > abs(pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE]) ) {
            do {
                nPrevNumCharges = pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES];
#if ( FIX_REM_PROTON_COUNT_BUG == 1 )
                nPrevRemovedProtons = pAATG->t_group_info->tni.nNumRemovedProtons;
#endif
                ret = bExistsAltPath( pBNS, pBD, pAATG, at, num_atoms,
                                cg_Minus /*nVertDoubleBond*/, cg_Plus /*nVertSingleBond*/, ALT_PATH_MODE_REM_PROTON );
                if ( IS_BNS_ERROR( ret ) ) {
                    return ret;
                }
#if ( FIX_REM_PROTON_COUNT_BUG == 1 )
                nCurrRemovedProtons = pAATG->t_group_info->tni.nNumRemovedProtons;
                if ( nCurrRemovedProtons != nPrevRemovedProtons ) {
                    return BNS_RADICAL_ERR;
                }
#endif
                if ( ret & 1 ) {
                    nDelta       = (ret & ~3) >> 2;
                    nNumChanges += (0 != (ret & 2));
                    if ( nDelta ) {
                        /* radical pair has disappeared */
                        ; /* goto quick_exit;*/
                    }
                    if ( nPrevNumCharges > pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] ) {
                        nNumNeutralized += (nPrevNumCharges - pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES])/2;
                    }
                }
            } while ( ret & 1 );
        }
#endif
        /* find alt path to remove one proton */
        do {
            /* remove a proton */
            nPrevNumCharges = pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES];
#if ( FIX_REM_PROTON_COUNT_BUG == 1 )
            nPrevRemovedProtons = pAATG->t_group_info->tni.nNumRemovedProtons;
#endif
            ret = bExistsAltPath( pBNS, pBD, pAATG, at, num_atoms,
                            tg_H /*nVertDoubleBond*/, cg_Plus /*nVertSingleBond*/, ALT_PATH_MODE_REM_PROTON );
            if ( IS_BNS_ERROR( ret ) ) {
                return ret;
            }
#if ( FIX_REM_PROTON_COUNT_BUG == 1 )
            nCurrRemovedProtons = pAATG->t_group_info->tni.nNumRemovedProtons;
            if ( nCurrRemovedProtons != nPrevRemovedProtons + (ret & 1) ) {
                return BNS_RADICAL_ERR;
            }
#endif
            if ( ret & 1 ) {
                nDelta       = (ret & ~3) >> 2;
                nNumChanges += (0 != (ret & 2));
                if ( nDelta ) {
                    /* radical pair has disappeared */
                    ; /* goto quick_exit;*/
                }
                nNumRemovedProtons ++;
                if ( nPrevNumCharges > pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] + 1 ) {
                    nNumNeutralized += (nPrevNumCharges - (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] - 1))/2;
                }
            }

        } while ( ret & 1 );

        /* neutralize: remove ion pairs like >N(+)=-O(-) => >N-=O */
        if ( (nNumRemovedProtons || bCancelChargesAlways) && cg_Minus >= num_atoms && cg_Plus >= num_atoms &&
             pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] > abs(pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE]) ) {
            do {
                nPrevNumCharges = pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES];
#if ( FIX_REM_PROTON_COUNT_BUG == 1 )
                nPrevRemovedProtons = pAATG->t_group_info->tni.nNumRemovedProtons;
#endif
                ret = bExistsAltPath( pBNS, pBD, pAATG, at, num_atoms,
                                cg_Minus /*nVertDoubleBond*/, cg_Plus /*nVertSingleBond*/, ALT_PATH_MODE_REM_PROTON );
                if ( IS_BNS_ERROR( ret ) ) {
                    return ret;
                }
#if ( FIX_REM_PROTON_COUNT_BUG == 1 )
                nCurrRemovedProtons = pAATG->t_group_info->tni.nNumRemovedProtons;
                if ( nCurrRemovedProtons != nPrevRemovedProtons ) {
                    return BNS_RADICAL_ERR;
                }
#endif
                if ( ret & 1 ) {
                    nDelta       = (ret & ~3) >> 2;
                    nNumChanges += (0 != (ret & 2));
                    if ( nDelta ) {
                        /* radical pair has disappeared */
                        ; /* goto quick_exit;*/
                    }
                    if ( nPrevNumCharges > pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] ) {
                        nNumNeutralized += (nPrevNumCharges - pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES])/2;
                    }
                }
            } while ( ret & 1 );
        }
    }
    ret = 0;
    if ( tg_H >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, tg_H, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
    if ( cg_Minus >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, cg_Minus, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
#if ( MOVE_PPLUS_TO_REMOVE_PROTONS == 1 )
    if ( cg_PlusP >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, cg_PlusP, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }
#endif
    if ( cg_Plus >= num_atoms ) {
        ret2 = RemoveLastGroupFromBnStruct( at, num_atoms, cg_Plus, pBNS );
        if ( !ret && ret2 )
            ret = ret2;
    }

    pBNS->type_CN   = 0;
    pBNS->type_T    = 0;
    pBNS->type_TACN = 0;
    
    if ( ret ) {
        return ret;
    }
    if ( pAATG->nAtTypeTotals[ATTOT_NUM_CO_Minus] + pAATG->nAtTypeTotals[ATTOT_NUM_ZO_Minus] &&
         pAATG->nAtTypeTotals[ATTOT_NUM_N_Minus] ) {
    }
    
    nPosCharges2 = (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] + pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    nNegCharges2 = (pAATG->nAtTypeTotals[ATTOT_NUM_CHARGES] - pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE])/2;
    /*
    nNumNP_H2    = pAATG->nAtTypeTotals[ATTOT_NUM_NP_H] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_NP_Proton];
    nNumOS_H2    = pAATG->nAtTypeTotals[ATTOT_NUM_COH] +
                   pAATG->nAtTypeTotals[ATTOT_NUM_CSH] + 
                   pAATG->nAtTypeTotals[ATTOT_NUM_ZOH];
    */
    if ( (nPosCharges - nNegCharges) - (nPosCharges2 - nNegCharges2) != nNumRemovedProtons ) {
        return BNS_PROGRAM_ERR;
    }

    if ( nNumCanceledCharges ) {
#if( FIX_CANCEL_CHARGE_COUNT_BUG == 1 )
        *nNumCanceledCharges += 2*nNumNeutralized;
#else
        *nNumCanceledCharges = 2*nNumNeutralized;
#endif
    }
    
    return nNumRemovedProtons;
}

/***************************************************************************************/
int mark_at_type( inp_ATOM *atom, int num_atoms, int nAtTypeTotals[] )
{
    int i, max_num_ions, mask, type;
    /*int max_protons, max_O_Minus, num_H = 0, num_CO=0;*/
    if ( nAtTypeTotals ) {
        memset( nAtTypeTotals, 0, ATTOT_ARRAY_LEN * sizeof(nAtTypeTotals[0]) );
    }
    for ( i = 0; i < num_atoms; i++ ) {
        type = GetAtomChargeType( atom, i, nAtTypeTotals, &mask, 0 );
        atom[i].at_type = type;
        /*
        num_H  += ((type & PR_HARD_TYP_H)  && (mask & ATBIT_MSK_H));
        num_CO += ((type & AR_HARD_TYP_HA) && (mask & AR_HARD_MSK_HA));
        */
    }
    if ( nAtTypeTotals ) {
        /*
        max_protons = nAtTypeTotals[ATTOT_NUM_NP_Proton] + 
                      inchi_min(num_H, nAtTypeTotals[ATTOT_NUM_NP_Plus]);
        max_O_Minus = nAtTypeTotals[ATTOT_NUM_CO_Minus] + nAtTypeTotals[ATTOT_NUM_CS_Minus] + 
                      nAtTypeTotals[ATTOT_NUM_ZO_Minus] + nAtTypeTotals[ATTOT_NUM_OO_Minus] +
                      nAtTypeTotals[ATTOT_NUM_ZOO_Minus] + nAtTypeTotals[ATTOT_NUM_NO_Minus] +
                      nAtTypeTotals[ATTOT_NUM_O_Minus] +nAtTypeTotals[ATTOT_NUM_N_Minus];
                      ;
        max_num_ions = max_protons + max_O_Minus + nAtTypeTotals[ATTOT_NUM_CHARGES];
        */
        max_num_ions = nAtTypeTotals[ATTOT_NUM_CHARGES];
    } else {
        max_num_ions = 0;
    }
    return max_num_ions;
}
/***************************************************************************************/
int RemoveNPProtonsAndAcidCharges( inp_ATOM *at, int num_atoms, BN_AATG *pAATG, BN_STRUCT *pBNS, BN_DATA *pBD )
{

    /* prepare data structure */
    int num;
    int nNumCanceledCharges          = 0;
    int nNumHardRemovedProtons       = 0;
    int nNumHardRemovedAcidicProtons = 0;
    T_GROUP_INFO *t_group_info = pAATG->t_group_info;
    int ret=0, bError = 0;
    int bAllowHardRemove = (t_group_info->bTautFlags & TG_FLAG_TEST_TAUT__SALTS) &&
                           (t_group_info->bTautFlags & TG_FLAG_TEST_TAUT2_SALTS) &&
                           (t_group_info->bTautFlags & TG_FLAG_MOVE_POS_CHARGES ) &&
                           (t_group_info->bTautFlags & TG_FLAG_HARD_ADD_REM_PROTONS);
    if ( pAATG->nMarkedAtom && num_atoms < pAATG->nAllocLen ) {
        inchi_free( pAATG->nMarkedAtom );
        qzfree( pAATG->nEndPoint );
        memset( pAATG, 0, sizeof(*pAATG) );
    }
    if ( !pAATG->nMarkedAtom && (pAATG->nMarkedAtom = (S_CHAR *) inchi_malloc( num_atoms * sizeof(pAATG->nMarkedAtom[0]))) ) {
        pAATG->nAllocLen = num_atoms;
        pAATG->nNumFound = 0;
    }
    /* simple remove protons from N, P, and O,S,Se,Te */
    if ( num = pAATG->nAtTypeTotals[ATTOT_NUM_NP_Proton] + pAATG->nAtTypeTotals[ATTOT_NUM_OH_Plus] ) {
        ret = SimpleRemoveHplusNPO(at, num_atoms, pAATG->nAtTypeTotals, t_group_info);
        if ( ret != num ) {
            bError = BNS_PROGRAM_ERR;
            goto exit_function;
        }
        /*t_group_info->nNumRemovedProtons  += ret;*/
        t_group_info->tni.bNormalizationFlags |= (ret > 0)? FLAG_PROTON_NPO_SIMPLE_REMOVED : 0;
    }
    if ( (num = pAATG->nAtTypeTotals[ATTOT_NUM_NP_Plus]) && bAllowHardRemove ) {
        /* try to 'hard' remove more protons from N; charges may be canceled */
        ret = HardRemoveHplusNP(at, num_atoms, 1, &nNumCanceledCharges, pAATG, pBNS, pBD);
        if ( IS_BNS_ERROR( ret ) ) {
            bError = ret;
            goto exit_function;
        }
        nNumHardRemovedProtons            += ret;
        /*t_group_info->nNumRemovedProtons  += ret;*/
        t_group_info->tni.bNormalizationFlags |= (ret > 0)? FLAG_PROTON_NP_HARD_REMOVED : 0;
    }
    if ( pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE] > 0 ) {
        ret = SimpleRemoveAcidicProtons( at, num_atoms, pAATG, pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE] );
        if ( IS_BNS_ERROR( ret ) ) {
            bError = ret;
            goto exit_function;
        }
        /*t_group_info->nNumRemovedProtons  += ret;*/
        t_group_info->tni.bNormalizationFlags |= (ret > 0)? FLAG_PROTON_AC_SIMPLE_REMOVED : 0;
        if ( pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE] > 0 && bAllowHardRemove ) {
            ret = HardRemoveAcidicProtons( at, num_atoms, pAATG, pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE], &nNumCanceledCharges, pBNS, pBD );
            if ( IS_BNS_ERROR( ret ) ) {
                bError = ret;
                goto exit_function;
            }
            if ( ret > 0 ) {
                int ret2 = SimpleRemoveAcidicProtons( at, num_atoms, pAATG, ret );
                if ( ret2 != ret ) {
                    bError = BNS_PROGRAM_ERR;
                    goto exit_function;
                }
                /*t_group_info->nNumRemovedProtons  += ret;*/
                t_group_info->tni.bNormalizationFlags |= (ret > 0)? FLAG_PROTON_AC_HARD_REMOVED : 0;
                nNumHardRemovedAcidicProtons      += ret;
            }
        }
    } else
    if ( pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE] < 0 ) {
        ret = SimpleAddAcidicProtons( at, num_atoms, pAATG, -pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE] );
        if ( IS_BNS_ERROR( ret ) ) {
            bError = ret;
            goto exit_function;
        }
        /*t_group_info->nNumRemovedProtons  -= ret;*/
        /*
           CHECK_TACN == 1 prohibits replacing (-) on N with H unless H can be moved to N
           along an alternating path from another heteroatom (t-group will be detected).
        */
        t_group_info->tni.bNormalizationFlags |= (ret > 0)? FLAG_PROTON_AC_SIMPLE_ADDED : 0;
        if ( pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE] < 0 && bAllowHardRemove ) {
            ret = HardAddAcidicProtons( at, num_atoms, pAATG, -pAATG->nAtTypeTotals[ATTOT_TOT_CHARGE], &nNumCanceledCharges, pBNS, pBD );
            if ( IS_BNS_ERROR( ret ) ) {
                bError = ret;
                goto exit_function;
            }
            if ( ret > 0 ) {
                int ret2 = SimpleAddAcidicProtons( at, num_atoms, pAATG, ret );
                if ( ret2 != ret ) {
                    bError = BNS_PROGRAM_ERR;
                    goto exit_function;
                }
                /*t_group_info->nNumRemovedProtons  -= ret;*/
                t_group_info->tni.bNormalizationFlags |= (ret > 0)? FLAG_PROTON_AC_HARD_ADDED : 0;
                nNumHardRemovedAcidicProtons      -= ret;
            }
        }
    }
    t_group_info->tni.bNormalizationFlags |= nNumCanceledCharges? FLAG_PROTON_CHARGE_CANCEL : 0;

exit_function:
    if ( bError ) {
        ret = IS_BNS_ERROR(bError)? bError : BNS_PROGRAM_ERR;
    }
    return ret;
}
/***************************************************************************************/
/* main normalization procedure */
int mark_alt_bonds_and_taut_groups ( inp_ATOM *at, inp_ATOM *at_fixed_bonds_out, int num_atoms,
                                     T_GROUP_INFO *t_group_info, INCHI_MODE *inpbTautFlags, INCHI_MODE *inpbTautFlagsDone )


{
    BN_STRUCT *pBNS = NULL;
    BN_DATA   *pBD  = NULL;
    int bError, nChanges, nTotChanges, taut_found, salt_found, taut_pass, salt_pass, salt_step, ret, ret2, num;
    int  nOrigDelta, num_changed_bonds;
    int max_altp = BN_MAX_ALTP;
    int bChangeFlow = (BNS_EF_CHNG_RSTR | BNS_EF_ALTR_BONDS);
    BNS_FLOW_CHANGES fcd[BNS_MAX_NUM_FLOW_CHANGES+1];
    C_GROUP_INFO CGroupInfo;
    C_GROUP_INFO *c_group_info = &CGroupInfo;
    S_GROUP_INFO SGroupInfo;
    S_GROUP_INFO *s_group_info = &SGroupInfo;
    INCHI_MODE    *pbTautFlags     = t_group_info? &t_group_info->bTautFlags     : inpbTautFlags;
    INCHI_MODE    *pbTautFlagsDone = t_group_info? &t_group_info->bTautFlagsDone : inpbTautFlagsDone;

    int nAtTypeTotals[ATTOT_ARRAY_LEN];
    int nNumOrigTotAtoms;

    BN_AATG  aatg;
    BN_AATG *pAATG = &aatg;
    
    nChanges = 0;
    bError   = 0;

    memset( c_group_info, 0, sizeof(*c_group_info) );
    memset( s_group_info, 0, sizeof(*s_group_info) );
    memset( pAATG,        0, sizeof(*pAATG) );

    if ( (*pbTautFlags & TG_FLAG_MOVE_POS_CHARGES) && num_atoms > 1 ) {
        /* charge groups memory allocation */
        c_group_info->c_group     = (C_GROUP *)inchi_calloc(num_atoms/2, sizeof(c_group_info->c_group[0]));
        c_group_info->c_candidate = (C_CANDIDATE*)inchi_calloc(num_atoms, sizeof(c_group_info->c_candidate[0]));
        if (c_group_info->c_group && c_group_info->c_candidate) {
            c_group_info->max_num_c_groups   = num_atoms/2;
            c_group_info->max_num_candidates = num_atoms;
        } else {
            bError = BNS_OUT_OF_RAM; /* error: out of RAM */
            /*printf("BNS_OUT_OF_RAM-1: num_at=%d, c_gr=%lx c_can=%lx\n", num_atoms, c_group_info->c_group, c_group_info->c_candidate);*/
            goto exit_function;
        }
    }

    if ( *pbTautFlags & TG_FLAG_TEST_TAUT__SALTS ) {
        if ( t_group_info ) {
            /* salt groups memory allocation */
            s_group_info->s_candidate = (S_CANDIDATE*)inchi_calloc(num_atoms, sizeof(s_group_info->s_candidate[0]));
            if (s_group_info->s_candidate) {
                s_group_info->max_num_candidates = num_atoms;
            } else {
                bError = BNS_OUT_OF_RAM; /* error: out of RAM */
                /*printf("BNS_OUT_OF_RAM-2\n");*/
                goto exit_function;
            }
        }
    }
    if ( t_group_info ) {
        if ( t_group_info->tGroupNumber )
            inchi_free( t_group_info->tGroupNumber );
        t_group_info->tGroupNumber = (AT_NUMB *)inchi_calloc( 2*num_atoms+1, sizeof(t_group_info->tGroupNumber[0]) );
        if ( !t_group_info->tGroupNumber )  {
            /*printf("BNS_OUT_OF_RAM-9\n");*/
            bError = BNS_OUT_OF_RAM; /* error: out of RAM */
            goto exit_function;
        }
        num = t_group_info->tni.nNumRemovedExplicitH;
        memset ( &t_group_info->tni, 0, sizeof(t_group_info->tni) );
        t_group_info->tni.nNumRemovedExplicitH = num;
    }    
/*
again:
*/
    /* allocate Balanced Network Data Strucures; replace Alternating bonds with Single */
    if ( (pBNS = AllocateAndInitBnStruct( at, num_atoms, BNS_ADD_ATOMS, BNS_ADD_EDGES, max_altp, &num_changed_bonds )) &&
         (pBD  = AllocateAndInitBnData( pBNS->max_vertices )) ) {

        pBNS->pbTautFlags     = pbTautFlags;     /* carry through all functions */
        pBNS->pbTautFlagsDone = pbTautFlagsDone;  /* carry through all functions */

#if( BNS_PROTECT_FROM_TAUT == 1 )
        /* protect bonds to acetyl and nitro */
        SetForbiddenEdges( pBNS, at, num_atoms, BNS_EDGE_FORBIDDEN_MASK );
#endif
        /* set bonds in case of input "aromatic" bonds or multiple radicals */
        ret = BnsAdjustFlowBondsRad( pBNS, pBD, at, num_atoms );
        /* (here pair(s) of radicals could have disappeared from the atoms) */
        if ( IS_BNS_ERROR( ret ) ) {
            bError = ret;
            goto exit_function;
        }
        pBNS->tot_st_flow += 2*ret;
        /*return 0;*/ /* debug */
        nOrigDelta = ret;
        if ( pBNS->tot_st_cap > pBNS->tot_st_flow ) {
            /* has radical */
            bChangeFlow |= BNS_EF_SET_NOSTEREO;
        }
        /******************************************************************** 
         *  Remove protons from NH(+), but not PH(+)
         *  Add protons to COO(-) etc.
         *  or remove protons from COOH etc to make the organic part neutral
         *  Note: for now (-) from N(-) can be only canceled or moved to -C=O
         ********************************************************************/
        if ( ( *pbTautFlags & TG_FLAG_VARIABLE_PROTONS ) && t_group_info && 
             mark_at_type( at, num_atoms, nAtTypeTotals ) &&
             nAtTypeTotals[ATTOT_NUM_CHARGES] ) {
            /*
               the structure is simple to neutralize if it yields exactly 
                 num[H(+)]     = num[N,P H(+)]
                 num[N,S,O(-)] = num[=C-O(-)] + num[C-S(-)] + num[N(-)] + num[other O(-), S(-)]

                 and n(p) = num[H(+)] - num[N,S,O(-)] (no protons, no negative N,O,S condition)

               Additional check is needed if:
                 min{num[N,PH], num[N,P(+), not onium]} > 0 
                     => possibility to yield more H(+)

               min_charge = orig_charge(P,N,O,S) - n(p) - n(OH,SH)
               max_charge = orig_charge(P,N,O,S) - n(p) + n(O,S,N(-)) 
            */


            nNumOrigTotAtoms = t_group_info->tni.nNumRemovedExplicitH + num_atoms;
            pAATG->nAtTypeTotals = nAtTypeTotals;
            pAATG->t_group_info  = t_group_info;
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
            pBNS->edge_forbidden_mask |= BNS_EDGE_FORBIDDEN_TEMP;
#endif
            /***********************************************************/
            /*                                                         */
            /*          ( D E ) P R O T O N A T I O N                  */
            /*                                                         */
            /***********************************************************/
            ret = RemoveNPProtonsAndAcidCharges( at, num_atoms, pAATG, pBNS, pBD );
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
            pBNS->edge_forbidden_mask &= ~BNS_EDGE_FORBIDDEN_TEMP;
#endif
            if ( IS_BNS_ERROR( ret ) ) {
                bError = ret;
                goto exit_function;
            }
            if ( t_group_info->tni.bNormalizationFlags ) {
                SetInitCapFlowToCurrent( pBNS );
                if ( at_fixed_bonds_out ) {
                    /* copy modified initial tautomeric structure for displaying
                       Warning: implicit H counts in at_fixed_bonds_out include explicit Hs */
                    memcpy( at_fixed_bonds_out, at, nNumOrigTotAtoms*sizeof(at_fixed_bonds_out[0]) );
                    /* -- will be done in FillOutInputInfAtom() --
                    RemoveExcessiveImplicitH( num_atoms, t_group_info->tni.nNumRemovedExplicitH, at_fixed_bonds_out );
                    */
                }
            }
        }
        /****************** initial bonds normalization ***************/
        if ( *pbTautFlags & TG_FLAG_MOVE_POS_CHARGES ) {
            /******************* find moveable positive charges **********************/
            do {  /* cycling while ret>0 added 2004-06-04 */
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                    pBNS->edge_forbidden_mask |= BNS_EDGE_FORBIDDEN_TEMP;
                    CorrectFixing_NH_NH_Bonds( pBNS, at, num_atoms );
#endif
                    ret = MarkChargeGroups ( at, num_atoms, c_group_info, t_group_info, pBNS, pBD );
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                    pBNS->edge_forbidden_mask &= ~BNS_EDGE_FORBIDDEN_TEMP;
#endif
                if ( IS_BNS_ERROR( ret ) ) {
                    bError = ret;
                    goto exit_function;
                }
                if ( ret ) {
                    nChanges += ret;
                    ret2 = AddCGroups2BnStruct( pBNS, at, num_atoms, c_group_info );
                    if ( IS_BNS_ERROR( ret2 ) ) {
                        bError = ret2;
                        goto exit_function;
                    }
                    *pbTautFlagsDone |= TG_FLAG_MOVE_POS_CHARGES_DONE;
                }
            } while ( ret > 0 );
#if( BNS_RAD_SEARCH == 1 )
#else
            /* moveable charges may allow to cancel radicals -- check it */
            if ( pBNS->tot_st_cap > pBNS->tot_st_flow ) {
                ret = BnsAdjustFlowBondsRad( pBNS, pBD, at, num_atoms );
                if ( IS_BNS_ERROR( ret ) ) {
                    bError = ret;
                    goto exit_function;
                }
                if ( ret > 0 ) {
                    /*
                    pBNS->tot_st_flow += 2*ret;
                    ret = ReInitBnStruct( pBNS, at, num_atoms, 1 );
                    if ( IS_BNS_ERROR( ret ) ) {
                        bError = ret;
                        goto exit_function;
                    }
                    */
                    bError = BNS_RADICAL_ERR;
                    goto exit_function;
                }
            }
#endif
        }
        /************************************************************************/
        /********          test bonds for bond tautomerism         **************/
        /******** replace moveable bonds with "alternating" bonds  **************/
        /************************************************************************/
        ret  = BnsTestAndMarkAltBonds( pBNS, pBD, at, num_atoms, fcd, bChangeFlow, 0 );
        if ( IS_BNS_ERROR( ret ) ) {
            bError = ret;
            goto exit_function;
        }
        nChanges += ret;
        /*********************** end of initial bonds normalization *************/
        nTotChanges = 0;
        /* check for tautomerism */
        /* find new tautomer groups */
        salt_pass  = 0;
        salt_step  = 0;
        salt_found = 0;

        /*************************************************************/
        /*                                                           */
        /*           M A I N   C Y C L E   B E G I N                 */
        /*                                                           */
        /*************************************************************/

        do {
            nTotChanges += nChanges;
            nChanges     = 0;
            taut_pass    = 0;

            /**************** regular bond/H/(-)/positive charges tautomerism cycle begin **************/
            do {
                taut_pass ++;
                for ( taut_found = 0; 0 < (ret=MarkTautomerGroups( at, num_atoms, t_group_info, c_group_info, pBNS, pBD )); taut_found ++ )
                    ;
                if ( ret < 0 ) {
                    bError = ret;
                }
                if ( taut_found && !salt_pass ) {
                    *pbTautFlagsDone |= TG_FLAG_TEST_TAUT__ATOMS_DONE;
                }
                if ( taut_found || salt_found ) {
                    /****************** repeat bonds normalization ***************/
                    ret = ReInitBnStructAddGroups( pBNS, at, num_atoms, t_group_info, c_group_info );
                    if ( IS_BNS_ERROR( ret ) ) {
                        bError = ret;
                        goto exit_function;
                    }
#if( BNS_RAD_SEARCH == 1 )
#else
                    /* discovered moveable charges and H-atoms may allow to cancel radicals */
                    if ( pBNS->tot_st_cap > pBNS->tot_st_flow ) {
                        ret = BnsAdjustFlowBondsRad( pBNS, pBD, at, num_atoms );
                        if ( IS_BNS_ERROR( ret ) ) {
                            bError = ret;
                            goto exit_function;
                        }
                        if ( ret > 0 ) {
                            /*
                            pBNS->tot_st_flow += 2*ret;
                            ret = ReInitBnStruct( pBNS, at, num_atoms, 1 );
                            if ( IS_BNS_ERROR( ret ) ) {
                                bError = ret;
                                goto exit_function;
                            }
                            */
                            bError = BNS_RADICAL_ERR;
                            goto exit_function;
                        }
                    }
#endif
                    /****************** update bonds normalization ***************/
                    if ( *pbTautFlags & TG_FLAG_MOVE_POS_CHARGES ) {
                        /******************* find moveable charges ***************/
                        do {  /* cycling while ret>0 added 2004-06-04 */
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                            pBNS->edge_forbidden_mask |= BNS_EDGE_FORBIDDEN_TEMP;
                            CorrectFixing_NH_NH_Bonds( pBNS, at, num_atoms );
#endif
                            ret = MarkChargeGroups ( at, num_atoms, c_group_info, t_group_info, pBNS, pBD );
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                            pBNS->edge_forbidden_mask &= ~BNS_EDGE_FORBIDDEN_TEMP;
#endif
                            if ( IS_BNS_ERROR( ret ) ) {
                                bError = ret;
                                goto exit_function;
                            }
                            nChanges+= ret;
                            if ( ret > 0 ) {
                                ret2 = ReInitBnStructAddGroups( pBNS, at, num_atoms, t_group_info, c_group_info );
                                if ( IS_BNS_ERROR( ret2 ) ) {
                                    bError = ret2;
                                    goto exit_function;
                                }
                                *pbTautFlagsDone |= TG_FLAG_MOVE_POS_CHARGES_DONE;
                            }
                        } while ( ret > 0 );
                    }

                    /************************************************************************/
                    /********             find moveable bonds:                 **************/
                    /********          test bonds for bond tautomerism         **************/
                    /******** replace moveable bonds with "alternating" bonds  **************/
                    /************************************************************************/
                    ret   = BnsTestAndMarkAltBonds( pBNS, pBD, at, num_atoms, fcd, bChangeFlow, 0 );
                    if ( IS_BNS_ERROR( ret ) ) {
                        bError = ret;
                        goto exit_function;
                    }
                    nChanges+= ret;
                    /****************** end of update bonds normalization ***************/
                }
                salt_found = 0;

            } while( taut_found && !bError );
            /**************** regular bond/H/(-)/positive charges tautomerism cycle end **************/
            
            if ( bError ) {
                break;
            }

            /******************* 'salt' tautomerism permitted *************************/
            if ( *pbTautFlags & TG_FLAG_TEST_TAUT__SALTS ) {

                if ( *pbTautFlags & TG_FLAG_TEST_TAUT2_SALTS ) {
                    /*********** requested one or more "salt" attachement migrartion test ********/
                    if ( !nChanges && salt_pass && salt_step ) {
                        break;  /* done */
                    }
                    if ( !salt_step ) { /* salt step 0: process one attachment migrartion */
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                        pBNS->edge_forbidden_mask |= BNS_EDGE_FORBIDDEN_TEMP;
                        CorrectFixing_NH_NH_Bonds( pBNS, at, num_atoms );
#endif
                        salt_found = MarkSaltChargeGroups ( at, num_atoms, s_group_info,
                                                            t_group_info, c_group_info, pBNS, pBD );
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                        pBNS->edge_forbidden_mask &= ~BNS_EDGE_FORBIDDEN_TEMP;
#endif
                        if ( salt_found < 0 ) {
                            bError = salt_found;
                            break;
                        } else
                        if ( salt_found > 0 ) {
                            *pbTautFlagsDone |= TG_FLAG_TEST_TAUT__SALTS_DONE;
                        }
                        salt_step = !salt_found; 
                         /* if new 'salt' atoms have been found then repeat regular taut. search
                          *      MarkTautomerGroups() and do not perform salt step 1
                          * if new 'salt' atoms have NOT been found then switch to salt step 1
                          *      and never repeat salt step 0 for the current structure
                          */
                    }
                    if ( salt_step /*||
                         (t_group_info->tni.bNormalizationFlags & FLAG_NORM_CONSIDER_TAUT)*/ ) { 
                        /* salt step 1: process more than one attachment migration */
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                        pBNS->edge_forbidden_mask |= BNS_EDGE_FORBIDDEN_TEMP;
                        CorrectFixing_NH_NH_Bonds( pBNS, at, num_atoms );
#endif
                        salt_found = MarkSaltChargeGroups2 ( at, num_atoms, s_group_info,
                                                            t_group_info, c_group_info, pBNS, pBD );
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                        pBNS->edge_forbidden_mask &= ~BNS_EDGE_FORBIDDEN_TEMP;
#endif
                        if ( salt_found < 0 ) {
                            bError = salt_found;
                            break;
                        } else
                        if ( salt_found == 1 || salt_found == 5 ) {
                            *pbTautFlagsDone |= TG_FLAG_TEST_TAUT2_SALTS_DONE;
                            if ( salt_found == 5 ) {
                                *pbTautFlagsDone |= TG_FLAG_TEST_TAUT3_SALTS_DONE;
                            }
                            /* salt_found == 2 => only negative charges involved */
                        }
                    }

                    salt_pass ++;

                } else { /* !( *pbTautFlags & TG_FLAG_TEST_TAUT2_SALTS ) */
                   /*************** requested only one attachement migrartion test **********/
                    if ( !nChanges && salt_pass ) { /* one attachment migrartion */
                        break;
                    }            /* salt step 0: process one attachment migrartion */
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                        pBNS->edge_forbidden_mask |= BNS_EDGE_FORBIDDEN_TEMP;
                        CorrectFixing_NH_NH_Bonds( pBNS, at, num_atoms );
#endif
                        salt_found = MarkSaltChargeGroups ( at, num_atoms, s_group_info,
                                                            t_group_info, c_group_info, pBNS, pBD );
#if ( RESET_EDGE_FORBIDDEN_MASK == 0 )
                        pBNS->edge_forbidden_mask &= ~BNS_EDGE_FORBIDDEN_TEMP;
#endif
                    if ( salt_found < 0 ) {
                        bError = salt_found;
                        break;
                    } else
                    if ( salt_found > 0 ) {
                        *pbTautFlagsDone |= TG_FLAG_TEST_TAUT__SALTS_DONE;
                    }
                    salt_pass ++;
                } /* ( *pbTautFlags & TG_FLAG_TEST_TAUT2_SALTS ) */
            } /* ( *pbTautFlags & TG_FLAG_TEST_TAUT__SALTS ) */

        } while ( salt_found && !bError );
        /*************************************************************/
        /*                                                           */
        /*           M A I N   C Y C L E   E N D                     */
        /*                                                           */
        /*************************************************************/

        if ( *pbTautFlags & TG_FLAG_MERGE_TAUT_SALTS ) {
            if ( !bError && s_group_info /*&& s_group_info->num_candidates > 0*/ ) {
                ret = MergeSaltTautGroups( at, num_atoms, s_group_info,
                                     t_group_info, c_group_info, pBNS );
                if ( ret < 0 ) {
                    bError = ret;
                } else
                if ( ret > 0 ) {
                    *pbTautFlagsDone |= TG_FLAG_MERGE_TAUT_SALTS_DONE;
                }
            }
        }
        if ( !bError && t_group_info &&
             (t_group_info->bTautFlags & TG_FLAG_VARIABLE_PROTONS) &&
             (t_group_info->bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE|TG_FLAG_FOUND_ISOTOPIC_H_DONE)) ) {
            ret = MakeIsotopicHGroup( at, num_atoms, s_group_info, t_group_info );
            if ( ret < 0 ) {
                bError = ret;
            }
        }
        /* success */
        remove_alt_bond_marks( at, num_atoms);
        
        /************************************************
         *  Temporarily ignore all non-alternating bonds
         *  and mark non-ring alt bonds non-stereogenic
         ************************************************/

        ReInitBnStructForAltBns( pBNS, at, num_atoms, 0 );
        MarkRingSystemsAltBns( pBNS, 0 );
        MarkNonStereoAltBns( pBNS, at, num_atoms, 0 );
#if( FIX_EITHER_DB_AS_NONSTEREO == 1 )
        /* second time unknown ("Either") alternating bonds are treated as non-stereogenic */
        /* stereobonds bonds that lost stereo get "Either" stereo_type */
        ReInitBnStructForAltBns( pBNS, at, num_atoms, 1 );
        MarkRingSystemsAltBns( pBNS, 1 );
        MarkNonStereoAltBns( pBNS, at, num_atoms, 1 );
#endif
    } else {
        bError = BNS_OUT_OF_RAM;
        /*printf("BNS_OUT_OF_RAM-3\n");*/
    }

exit_function:    
    pBNS = DeAllocateBnStruct( pBNS );
    pBD  = DeAllocateBnData( pBD );
/*#if( MOVE_CHARGES == 1 )*/
    if ( c_group_info )  {
        if ( c_group_info->c_group ) {
            inchi_free( c_group_info->c_group );
        }
        if ( c_group_info->c_candidate ) {
            inchi_free( c_group_info->c_candidate );
        }
    }
/*#endif*/
    if ( s_group_info && s_group_info->s_candidate ) {
        inchi_free( s_group_info->s_candidate );
    }
    if ( pAATG && pAATG->nMarkedAtom ) {
        inchi_free( pAATG->nMarkedAtom );
        qzfree( pAATG->nEndPoint );
        /*qzfree( pAATG->nAtTypeTotals );*/ /* nAtTypeTotals is a stack array */
    }
    if ( t_group_info && t_group_info->tGroupNumber ) {
        inchi_free( t_group_info->tGroupNumber );
        t_group_info->tGroupNumber = NULL;
    }

    if ( !bError && num_atoms == 1 && at[0].at_type == ATT_PROTON && t_group_info && !t_group_info->tni.nNumRemovedExplicitH ) {
        /* remove single isolated proton */
        t_group_info->tni.nNumRemovedProtons   = 1;
        t_group_info->tni.bNormalizationFlags |= FLAG_PROTON_SINGLE_REMOVED;
        if ( at[0].iso_atw_diff ) {
            t_group_info->tni.nNumRemovedProtonsIsotopic[at[0].iso_atw_diff-1] ++;
        }
        if ( at_fixed_bonds_out ) {
            memcpy( at_fixed_bonds_out, at, num_atoms*sizeof(at_fixed_bonds_out[0]) );
        }
        /*num_atoms --;*/
    }


    /*
       Additional currently unused info:

       nOrigDelta > 0: original structure has been changed
                       due to fiund augmenting path(s)
       nChanges   > 0: either alt. bonds or taut. groups have been found 
    */

    return bError? bError : num_atoms;  /* ret = 0 => success, any other => error */

}
/*********************************************************************************/
int nMaxFlow2Check( BN_STRUCT *pBNS, int iedge )
{
    BNS_EDGE *pEdge    = pBNS->edge + iedge;
    int       nMaxFlow = (pEdge->cap & EDGE_FLOW_MASK); /* edge cap */

    if ( nMaxFlow > MAX_BOND_EDGE_CAP ) {
        nMaxFlow = MAX_BOND_EDGE_CAP;
    }
    return nMaxFlow;
}
/*********************************************************************************/
int nCurFlow2Check( BN_STRUCT *pBNS, int iedge )
{
    BNS_EDGE *pEdge    = pBNS->edge + iedge;
    int       nCurFlow = (pEdge->flow & EDGE_FLOW_MASK); /* edge flow */
    return    nCurFlow;
}
/*********************************************************************************/
int nMinFlow2Check( BN_STRUCT *pBNS, int iedge )
{
    BNS_EDGE *pEdge    = pBNS->edge + iedge;
    Vertex    v1       = pEdge->neighbor1;
    Vertex    v2       = v1 ^ pEdge->neighbor12;
    int       f12      = (pEdge->flow & EDGE_FLOW_MASK);
    int       rescap1, rescap2, rescap12, i, iedge_i;

    if ( f12 > 0 ) {
        for ( i = 0, rescap1 = 0; i < pBNS->vert[v1].num_adj_edges; i ++ ) {
            iedge_i = pBNS->vert[v1].iedge[i];
            if ( iedge_i == iedge )
                continue;
            rescap1 += (pBNS->edge[iedge_i].cap & EDGE_FLOW_MASK) - (pBNS->edge[iedge_i].flow & EDGE_FLOW_MASK);
        }
        for ( i = 0, rescap2 = 0; i < pBNS->vert[v2].num_adj_edges; i ++ ) {
            iedge_i = pBNS->vert[v2].iedge[i];
            if ( iedge_i == iedge )
                continue;
            rescap2 += (pBNS->edge[iedge_i].cap & EDGE_FLOW_MASK) - (pBNS->edge[iedge_i].flow & EDGE_FLOW_MASK);
        }
        rescap12 = inchi_min( rescap1, rescap2 );
        rescap12 = inchi_min( rescap12, f12 );
        return f12-rescap12;
    }
    return 0;
}
/**********************************************************************************/
int bSetBondsAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd, int nTestFlow, inp_ATOM *at, int num_atoms, int bChangeFlow0 )
{
    int    ifcd, iedge, new_flow, ret_val, nChanges = 0, bError=0;
    int    bChangeFlow;
    Vertex v1, v2;
    int    ineigh1, ineigh2;
    BNS_EDGE *pEdge;
    
    bChangeFlow0 &= ~BNS_EF_CHNG_RSTR;  /* do not change pEdge flow in SetBondType */
    if ( !bChangeFlow0 )
        return 0;

    bChangeFlow = (bChangeFlow0 & ~BNS_EF_SET_NOSTEREO);
    /* find the next to the last changed */
    if ( bChangeFlow0 & BNS_EF_SET_NOSTEREO ) {

        for ( ifcd = 0; NO_VERTEX != (iedge = fcd[ifcd].iedge); ifcd ++ ) {
            iedge = fcd[ifcd].iedge;
            pEdge = pBNS->edge + iedge;
            if ( !pEdge->pass ) {
                continue;
            }
        
            if ( !ifcd && nTestFlow>=0 ) {
                new_flow = nTestFlow;
            } else {
                new_flow = (int)pEdge->flow;
            }

            v1 = pEdge->neighbor1;
            v2 = pEdge->neighbor12 ^ v1;
            if ( v1 < num_atoms && v2 < num_atoms && new_flow != pEdge->flow0 ) {
                if ( (pBNS->vert[v1].st_edge.cap0 == pBNS->vert[v1].st_edge.flow0) !=
                     (pBNS->vert[v1].st_edge.cap  == pBNS->vert[v1].st_edge.flow ) ||
                     (pBNS->vert[v2].st_edge.cap0 == pBNS->vert[v2].st_edge.flow0) !=
                     (pBNS->vert[v2].st_edge.cap  == pBNS->vert[v2].st_edge.flow )) {
                    bChangeFlow |= BNS_EF_SET_NOSTEREO;
                    nChanges    |= BNS_EF_SET_NOSTEREO;
                }
            }
        }

    } else {
        for ( ifcd = 0; NO_VERTEX != (iedge = fcd[ifcd].iedge); ifcd ++ )
            ;
    }

    /* restore in reversed order to correctly handle vertex changed more than once */
    for ( ifcd -= 1; 0 <= ifcd; ifcd -- ) {

        iedge = fcd[ifcd].iedge;
        pEdge = pBNS->edge + iedge;
        if ( !pEdge->pass ) {
            continue;
        }
        
        if ( !ifcd && nTestFlow>=0 ) {
            new_flow = nTestFlow;
        } else {
            new_flow = (int)pEdge->flow;
        }

        v1 = pEdge->neighbor1;
        v2 = pEdge->neighbor12 ^ v1;
        if ( v1 < num_atoms && v2 < num_atoms && bChangeFlow && new_flow != pEdge->flow0 ) {
            ineigh1 = pEdge->neigh_ord[0];
            ineigh2 = pEdge->neigh_ord[1];
            ret_val = SetAtomBondType( pEdge, &at[v1].bond_type[ineigh1], &at[v2].bond_type[ineigh2], new_flow-pEdge->flow0, bChangeFlow );
            if ( !IS_BNS_ERROR( ret_val ) ) {
                nChanges |= (ret_val > 0);
            } else {
                bError = ret_val;
            }
        }
        pEdge->pass = 0;
    }
    return bError? bError : nChanges;
}
/**********************************************************************************/
int bRestoreFlowAfterCheckOneBond( BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd )
{
    int ifcd, iedge;
    Vertex v1, v2;
    BNS_EDGE *pEdge;
    
    /* find the next to the last changed */
    for ( ifcd = 0; NO_VERTEX != (iedge = fcd[ifcd].iedge); ifcd ++ )
        ;

    /* restore in reversed order to correctly handle vertex changed more than once */
    for ( ifcd -= 1; 0 <= ifcd; ifcd -- ) {

        /* restore edge flow & cap */
        iedge = fcd[ifcd].iedge;
        pEdge = pBNS->edge + iedge;
        pEdge->flow = fcd[ifcd].flow;
        pEdge->cap  = fcd[ifcd].cap;
        pEdge->pass = 0;

        /* restore st-flow, cap */
        if ( NO_VERTEX != (v1 = fcd[ifcd].v1) ) {
            pBNS->vert[v1].st_edge.flow = fcd[ifcd].flow_st1;
            pBNS->vert[v1].st_edge.cap  = fcd[ifcd].cap_st1;
            pBNS->vert[v1].st_edge.pass = 0;
        }
        if ( NO_VERTEX != (v2 = fcd[ifcd].v2) ) {
            pBNS->vert[v2].st_edge.flow = fcd[ifcd].flow_st2;
            pBNS->vert[v2].st_edge.cap  = fcd[ifcd].cap_st2;
            pBNS->vert[v2].st_edge.pass = 0;
        }
    }
    return 0;
}
/**********************************************************************************/
int bSetFlowToCheckOneBond( BN_STRUCT *pBNS, int iedge, int flow, BNS_FLOW_CHANGES *fcd )
{
    BNS_EDGE *pEdge    = pBNS->edge + iedge;
    int       f12      = (pEdge->flow & EDGE_FLOW_MASK); /* the original flow */
    int       ifcd     = 0;
    int       nDots    = 0;
    int       i, iedge_i;

    fcd[ifcd].iedge    = NO_VERTEX;

    if ( f12 < flow ) {
        /* Increase edge flow: Grab flow from the neighbors and delete it: set flow12=cap12 = 0 */
        /************************************************************************************/
        /* For example, simulate a new fixed double bond in place of a single bond and      */
        /* creates ONE or NONE (in case of a radical on adjacent atom) augmenting paths and */
        /* makes it impossible for the BNS to set same flow as it originally was            */
        /************************************************************************************/
        Vertex    v1       = pEdge->neighbor1;
        Vertex    v2       = v1 ^ pEdge->neighbor12;
        Vertex    v_i;        /* neighbor of v1 or v2 */
        BNS_EDGE *pEdge_i;
        int       delta1, delta2, f, st_edge_rescap;

        if ( (pBNS->vert[v1].st_edge.cap  & EDGE_FLOW_ST_MASK) < flow ||
             (pBNS->vert[v2].st_edge.cap  & EDGE_FLOW_ST_MASK) < flow ) {
            return BNS_CANT_SET_BOND;
        }
        if ( (pBNS->vert[v1].st_edge.flow & EDGE_FLOW_ST_MASK) < f12 ||
             (pBNS->vert[v2].st_edge.flow & EDGE_FLOW_ST_MASK) < f12 ) {
            return BNS_CAP_FLOW_ERR;
        }

        
        fcd[ifcd].iedge    = iedge;
        fcd[ifcd].flow     = pEdge->flow;
        fcd[ifcd].cap      = pEdge->cap;
        
        fcd[ifcd].v1       = v1;
        fcd[ifcd].flow_st1 = pBNS->vert[v1].st_edge.flow;
        fcd[ifcd].cap_st1  = pBNS->vert[v1].st_edge.cap;
        
        fcd[ifcd].v2       = v2;
        fcd[ifcd].flow_st2 = pBNS->vert[v2].st_edge.flow;
        fcd[ifcd].cap_st2  = pBNS->vert[v2].st_edge.cap;

        fcd[++ifcd].iedge    = NO_VERTEX; /* mark the end of the fcd[] data */
        pEdge->pass |= 64;

        delta1 = delta2 = flow - f12;

        if ( f12 > 0 ) {
            /* remove old edge flow from the flow and cap of the adjacent vertices' st-edges */
            pBNS->vert[v1].st_edge.cap  = ((pBNS->vert[v1].st_edge.cap  & EDGE_FLOW_ST_MASK)-f12) | (pBNS->vert[v1].st_edge.cap  & ~EDGE_FLOW_ST_MASK);
            pBNS->vert[v2].st_edge.cap  = ((pBNS->vert[v2].st_edge.cap  & EDGE_FLOW_ST_MASK)-f12) | (pBNS->vert[v2].st_edge.cap  & ~EDGE_FLOW_ST_MASK);
            pBNS->vert[v1].st_edge.flow = ((pBNS->vert[v1].st_edge.flow & EDGE_FLOW_ST_MASK)-f12) | (pBNS->vert[v1].st_edge.flow & ~EDGE_FLOW_ST_MASK);
            pBNS->vert[v2].st_edge.flow = ((pBNS->vert[v2].st_edge.flow & EDGE_FLOW_ST_MASK)-f12) | (pBNS->vert[v2].st_edge.flow & ~EDGE_FLOW_ST_MASK);
            /* delete current edge flow and capacity */
            pEdge->flow = (pEdge->flow & ~EDGE_FLOW_MASK);
        }
        pEdge->cap  = (pEdge->cap  & ~EDGE_FLOW_MASK);

        /* grab the adjacent vertex1 radical (st_edge_rescap) if it exists */
        st_edge_rescap = (pBNS->vert[v1].st_edge.cap & EDGE_FLOW_ST_MASK) - (pBNS->vert[v1].st_edge.flow & EDGE_FLOW_ST_MASK);
        while ( st_edge_rescap && delta1 ) { 
            st_edge_rescap --;  /* grab the radical */
            delta1 --;
            pBNS->vert[v1].st_edge.cap = ((pBNS->vert[v1].st_edge.cap & EDGE_FLOW_ST_MASK)-1) | (pBNS->vert[v1].st_edge.cap & ~EDGE_FLOW_ST_MASK);
            nDots --;
        }
        /* grab the adjacent vertex2 radical (st_edge_rescap) if it exists */
        st_edge_rescap = (pBNS->vert[v2].st_edge.cap & EDGE_FLOW_ST_MASK) - (pBNS->vert[v2].st_edge.flow & EDGE_FLOW_ST_MASK);
        while ( st_edge_rescap && delta2 ) {
            st_edge_rescap --;  /* grab the radical */
            delta2 --;
            pBNS->vert[v2].st_edge.cap = ((pBNS->vert[v2].st_edge.cap & EDGE_FLOW_ST_MASK)-1) | (pBNS->vert[v2].st_edge.cap & ~EDGE_FLOW_ST_MASK);
            nDots --;
        }
        /* grab flows from v1 neighbors */
        for ( i = 0; delta1 && i < pBNS->vert[v1].num_adj_edges; i ++ ) {
            iedge_i = pBNS->vert[v1].iedge[i];
            if ( iedge_i == iedge )
                continue;
            pEdge_i = pBNS->edge + iedge_i;
            if ( IS_FORBIDDEN(pEdge_i->forbidden, pBNS) )
                continue;
            f = (pEdge_i->flow & EDGE_FLOW_MASK);
            if ( f ) {
                v_i     = pEdge_i->neighbor12 ^ v1;

                fcd[ifcd].iedge    = iedge_i;
                fcd[ifcd].flow     = pEdge_i->flow;
                fcd[ifcd].cap      = pEdge_i->cap;
                
                fcd[ifcd].v1       = v_i;
                fcd[ifcd].flow_st1 = pBNS->vert[v_i].st_edge.flow;
                fcd[ifcd].cap_st1  = pBNS->vert[v_i].st_edge.cap;
                
                fcd[ifcd].v2       = NO_VERTEX;
                fcd[ifcd].flow_st2 = 0;
                fcd[ifcd].cap_st2  = 0;

                fcd[++ifcd].iedge    = NO_VERTEX; /* mark the end of the fcd[] data */
                pEdge_i->pass |= 64;

                while ( f && delta1 ) {
                    f --;
                    delta1 --;
                    pEdge_i->flow = ((pEdge_i->flow & EDGE_FLOW_MASK) - 1) | (pEdge_i->flow & ~EDGE_FLOW_MASK);
                    pBNS->vert[v_i].st_edge.flow = ((pBNS->vert[v_i].st_edge.flow & EDGE_FLOW_ST_MASK)-1) | (pBNS->vert[v_i].st_edge.flow & ~EDGE_FLOW_ST_MASK);
                    /* next 2 lines added 01-22-2002 */
                    pBNS->vert[v1].st_edge.cap  = ((pBNS->vert[v1].st_edge.cap  & EDGE_FLOW_ST_MASK)-1) | (pBNS->vert[v1].st_edge.cap  & ~EDGE_FLOW_ST_MASK);
                    pBNS->vert[v1].st_edge.flow = ((pBNS->vert[v1].st_edge.flow & EDGE_FLOW_ST_MASK)-1) | (pBNS->vert[v1].st_edge.flow & ~EDGE_FLOW_ST_MASK);

                    nDots ++;
                }
                    
            }
        }
        /* grab flows from v2 neighbors */
        for ( i = 0; delta2 && i < pBNS->vert[v2].num_adj_edges; i ++ ) {
            iedge_i = pBNS->vert[v2].iedge[i];
            if ( iedge_i == iedge )
                continue;
            pEdge_i = pBNS->edge + iedge_i;
            if ( IS_FORBIDDEN(pEdge_i->forbidden, pBNS) )
                continue;
            f = (pEdge_i->flow & EDGE_FLOW_MASK);
            if ( f ) {
                v_i     = pEdge_i->neighbor12 ^ v2;

                fcd[ifcd].iedge    = iedge_i;
                fcd[ifcd].flow     = pEdge_i->flow;
                fcd[ifcd].cap      = pEdge_i->cap;
                
                fcd[ifcd].v1       = v_i;
                fcd[ifcd].flow_st1 = pBNS->vert[v_i].st_edge.flow;
                fcd[ifcd].cap_st1  = pBNS->vert[v_i].st_edge.cap;
                
                fcd[ifcd].v2       = NO_VERTEX;
                fcd[ifcd].flow_st2 = 0;
                fcd[ifcd].cap_st2  = 0;

                fcd[++ifcd].iedge    = NO_VERTEX; /* mark the end of the fcd[] data */
                pEdge_i->pass |= 64;

                while ( f && delta2 ) {
                    f --;
                    delta2 --;
                    pEdge_i->flow = ((pEdge_i->flow & EDGE_FLOW_MASK) - 1) | (pEdge_i->flow & ~EDGE_FLOW_MASK);
                    pBNS->vert[v_i].st_edge.flow = ((pBNS->vert[v_i].st_edge.flow & EDGE_FLOW_ST_MASK)-1) | (pBNS->vert[v_i].st_edge.flow & ~EDGE_FLOW_ST_MASK);
                    /* next 2 lines added 01-22-2002 */
                    pBNS->vert[v2].st_edge.cap  = ((pBNS->vert[v2].st_edge.cap  & EDGE_FLOW_ST_MASK)-1) | (pBNS->vert[v2].st_edge.cap  & ~EDGE_FLOW_ST_MASK);
                    pBNS->vert[v2].st_edge.flow = ((pBNS->vert[v2].st_edge.flow & EDGE_FLOW_ST_MASK)-1) | (pBNS->vert[v2].st_edge.flow & ~EDGE_FLOW_ST_MASK);

                    nDots ++;
                }
                    
            }
        }
        if ( delta1 || delta2 ) {
            return BNS_CANT_SET_BOND;
        }
    }

    if ( f12 >= flow ) {
        /* Decrease edge flow: Redirect flow to the neighbors and delete it on the edge: set flow12=cap12 = 0 */
        /* f12==flow fixes flow through the edge so that BNS cannot change it */
        /**********************************************************************************************/
        /* For example, simulate a removal of a double bond and create ONE or NONE augmenting path    */
        /* Make it impossible for BNS to set same flow as it originally was                           */
        /**********************************************************************************************/
        Vertex    v1       = pEdge->neighbor1;
        Vertex    v2       = (v1 ^ pEdge->neighbor12);
        int       delta;
        /* if NOT (st-cap >= st-flow >= f12 >= flow) then error in the BN structure */
        if ( (pBNS->vert[v1].st_edge.flow & EDGE_FLOW_ST_MASK) < f12 ||
             (pBNS->vert[v2].st_edge.flow & EDGE_FLOW_ST_MASK) < f12 ||
             (pBNS->vert[v1].st_edge.cap  & EDGE_FLOW_ST_MASK) < flow ||
             (pBNS->vert[v2].st_edge.cap  & EDGE_FLOW_ST_MASK) < flow ) {
            return BNS_CAP_FLOW_ERR;
        }
        fcd[ifcd].iedge    = iedge;
        fcd[ifcd].flow     = pEdge->flow;
        fcd[ifcd].cap      = pEdge->cap;
        
        fcd[ifcd].v1       = v1;
        fcd[ifcd].flow_st1 = pBNS->vert[v1].st_edge.flow;
        fcd[ifcd].cap_st1  = pBNS->vert[v1].st_edge.cap;
        
        fcd[ifcd].v2       = v2;
        fcd[ifcd].flow_st2 = pBNS->vert[v2].st_edge.flow;
        fcd[ifcd].cap_st2  = pBNS->vert[v2].st_edge.cap;

        fcd[++ifcd].iedge    = NO_VERTEX; /* mark the end of the fcd[] data */
        pEdge->pass |= 64;

        delta = f12 - flow;
        /* remove current edge flow from st-edges */
        /* -- seem to be a bug --
        pBNS->vert[v1].st_edge.flow = ((pBNS->vert[v1].st_edge.flow & EDGE_FLOW_ST_MASK)-delta) | (pBNS->vert[v1].st_edge.flow & ~EDGE_FLOW_ST_MASK);
        pBNS->vert[v2].st_edge.flow = ((pBNS->vert[v2].st_edge.flow & EDGE_FLOW_ST_MASK)-delta) | (pBNS->vert[v2].st_edge.flow & ~EDGE_FLOW_ST_MASK);
        */

        /* replacement to the above 2 lines 01-16-2002 */
        /* remove old edge flow from the flow of the adjacent vertices' st-edges */

        pBNS->vert[v1].st_edge.flow = ((pBNS->vert[v1].st_edge.flow & EDGE_FLOW_ST_MASK)-f12) | (pBNS->vert[v1].st_edge.flow & ~EDGE_FLOW_ST_MASK);
        pBNS->vert[v2].st_edge.flow = ((pBNS->vert[v2].st_edge.flow & EDGE_FLOW_ST_MASK)-f12) | (pBNS->vert[v2].st_edge.flow & ~EDGE_FLOW_ST_MASK);

        /* added 01-16-2002: reduce st-cap if new flow > 0  */
        /* remove new edge flow from the cap of the adjacent vertices' st-edges */
        pBNS->vert[v1].st_edge.cap = ((pBNS->vert[v1].st_edge.cap & EDGE_FLOW_ST_MASK)-flow) | (pBNS->vert[v1].st_edge.cap & ~EDGE_FLOW_ST_MASK);
        pBNS->vert[v2].st_edge.cap = ((pBNS->vert[v2].st_edge.cap & EDGE_FLOW_ST_MASK)-flow) | (pBNS->vert[v2].st_edge.cap & ~EDGE_FLOW_ST_MASK);
        
        /* delete current edge flow and capacity */
        pEdge->flow = (pEdge->flow & ~EDGE_FLOW_MASK);
        pEdge->cap  = (pEdge->cap  & ~EDGE_FLOW_MASK);
        nDots = 2*delta;
    }
    return nDots;
}
/**********************************************************************************/
/* Connect new (fictitious, temporary) vertex to to nVertDoubleBond by a new edge */
/* Add radical (set st-cap=1) to the new vertex, set cap=1 to the new edge        */
/* Add radical (set st-cap=1) to nVertSingleBond                                  */
/* Find augmenting path connecting new vertex to nVertSingleBond                  */
/* This corresponds to moving H-atom from nVertSingleBond to nVertDoubleBond      */
/**********************************************************************************/
int bAddNewVertex( BN_STRUCT *pBNS, int nVertDoubleBond, int nCap, int nFlow, int nMaxAdjEdges, int *nDots )
{
    Vertex      vlast    = pBNS->num_vertices - 1;
    Vertex      vnew     = pBNS->num_vertices;
    Vertex      v2       = nVertDoubleBond;
    BNS_VERTEX *pVert2   = pBNS->vert + v2;   /* pointer to an old vertex */
    BNS_VERTEX *pNewVert = pBNS->vert + vnew; /* pointer to a new vertex */
    
    EdgeIndex   iedge    = pBNS->num_edges;
    BNS_EDGE   *pEdge    = pBNS->edge + iedge; /* pointer to a new edge */

    if ( iedge >= pBNS->max_edges || vnew >= pBNS->max_vertices ) {
        return BNS_VERT_EDGE_OVFL; /* edges or vertices overflow */
    }
    if ( (pBNS->vert[vlast].iedge - pBNS->iedge) + pBNS->vert[vlast].max_adj_edges + nMaxAdjEdges >= pBNS->max_iedges ) {
        return BNS_VERT_EDGE_OVFL; /* iedges overflow */
    }
    if ( pVert2->num_adj_edges >= pVert2->max_adj_edges || nMaxAdjEdges <= 0 ) {
        return BNS_VERT_EDGE_OVFL; /* neighbors overflow */
    }
    /* fill out the new edge, set its cap and flow, connect */
    /* memset( pEdge, 0, sizeof(*pEdge) ); */
    pEdge->cap  = pEdge->cap0  = nCap;
    pEdge->flow = pEdge->flow0 = nFlow;
    pEdge->pass = 0;
    pEdge->neighbor1  = v2;
    pEdge->neighbor12 = v2 ^ vnew;
    pEdge->forbidden  = 0;
    /* fill out the new vertex */
    /* memset( pNewVert, 0, sizeof(*pNewVert) ); */
    pNewVert->max_adj_edges = nMaxAdjEdges;
    pNewVert->num_adj_edges = 0;
    pNewVert->st_edge.cap0  = pNewVert->st_edge.cap   = nCap;
    pNewVert->st_edge.flow0 = pNewVert->st_edge.flow  = nFlow;
    pNewVert->st_edge.pass  = 0; /* add initialization; added 2006-03-25 */
    pNewVert->iedge         = pBNS->vert[vlast].iedge + pBNS->vert[vlast].max_adj_edges;
    pNewVert->type          = BNS_VERT_TYPE_TEMP;
    *nDots += nCap - nFlow;

    pEdge->neigh_ord[v2>vnew] = pVert2->num_adj_edges;
    pEdge->neigh_ord[v2<vnew] = pNewVert->num_adj_edges;

    /* connect new edge to v2 */
    pVert2->iedge[pVert2->num_adj_edges ++] = iedge;
    /* connect new edge to vnew */
    pNewVert->iedge[pNewVert->num_adj_edges ++] = iedge;

    /* fix v2 flow and cap */
    *nDots -= (int)pVert2->st_edge.cap - (int)pVert2->st_edge.flow;
    pVert2->st_edge.flow += nFlow;
    if ( pVert2->st_edge.cap < pVert2->st_edge.flow ) {
        pVert2->st_edge.cap = pVert2->st_edge.flow;
    }
    *nDots += (int)pVert2->st_edge.cap - (int)pVert2->st_edge.flow;

    pBNS->num_edges ++;
    pBNS->num_vertices ++;

    return vnew;
}

/*****************************************************************************************************/
int AddNewEdge( BNS_VERTEX *p1, BNS_VERTEX *p2, BN_STRUCT *pBNS, int nEdgeCap, int nEdgeFlow )
{
    int ip1 = p1 - pBNS->vert;
    int ip2 = p2 - pBNS->vert;
    int ie  = pBNS->num_edges;
    BNS_EDGE *e = pBNS->edge + ie;
    /* debug: check bounds */
    if ( ip1 >= pBNS->max_vertices || ip1 < 0 ||
         ip2 >= pBNS->max_vertices || ip2 < 0 ||
         ie  >= pBNS->max_edges    || ie  < 0 ||
         (p1->iedge - pBNS->iedge) < 0 ||
         (p1->iedge - pBNS->iedge) + p1->max_adj_edges > pBNS->max_iedges ||
         (p2->iedge - pBNS->iedge) < 0 ||
         (p2->iedge - pBNS->iedge) + p2->max_adj_edges > pBNS->max_iedges ||
         p1->num_adj_edges >= p1->max_adj_edges ||
         p2->num_adj_edges >= p2->max_adj_edges  ) {
        return BNS_VERT_EDGE_OVFL;
    }
    /* clear the edge */
    memset( e, 0, sizeof(*e) );
    /* connect */
    e->neighbor1  = inchi_min( ip1, ip2 );
    e->neighbor12 = ip1 ^ ip2;
    p1->iedge[p1->num_adj_edges] = ie;
    p2->iedge[p2->num_adj_edges] = ie;
    e->neigh_ord[ip1 > ip2] = p1->num_adj_edges ++;
    e->neigh_ord[ip1 < ip2] = p2->num_adj_edges ++;
    e->cap = e->cap0 = nEdgeCap;
    e->flow = e->flow0 = nEdgeFlow;
    p1->st_edge.flow += nEdgeFlow;
    p2->st_edge.flow += nEdgeFlow;
    if ( p1->st_edge.cap < p1->st_edge.flow ) {
        p1->st_edge.cap = p1->st_edge.flow;
    }
    if ( p2->st_edge.cap < p2->st_edge.flow ) {
        p2->st_edge.cap = p2->st_edge.flow;
    }
    pBNS->num_edges ++;
    return ie;
}
/**********************************************************************************/
BNS_IEDGE GetEdgeToGroupVertex( BN_STRUCT *pBNS, Vertex v1, AT_NUMB type)
{
    if ( v1 < pBNS->num_atoms ) {
        Vertex      v2;
        BNS_EDGE   *pEdge1;
        BNS_VERTEX *pVert1 = pBNS->vert+v1;
        int         i      = pVert1->num_adj_edges-1;

        while( 0 <= i ) {
            pEdge1 = pBNS->edge + pVert1->iedge[i];
            v2 = pEdge1->neighbor12 ^ v1;
            if ( pBNS->vert[v2].type == type ) {
                return IS_FORBIDDEN(pEdge1->forbidden, pBNS)? NO_VERTEX : pVert1->iedge[i];
            }
            i --;
        }
        return NO_VERTEX; /* not found t-group */
    } else
    if ( v1 < pBNS->num_vertices ) {
        return NO_VERTEX;
    }
    return BNS_VERT_EDGE_OVFL;
}
/**********************************************************************************/
Vertex GetGroupVertex(BN_STRUCT *pBNS, Vertex v1, AT_NUMB type)
{
    if ( v1 < pBNS->num_atoms ) {
        Vertex      v2;
        BNS_EDGE   *pEdge1;
        BNS_VERTEX *pVert1 = pBNS->vert+v1;
        int         i      = pVert1->num_adj_edges-1;

        AT_NUMB type2;

        if ( type == BNS_VERT_TYPE_ENDPOINT )
            type2 = BNS_VERT_TYPE_TGROUP;
        else
        if ( type == BNS_VERT_TYPE_C_POINT )
            type2 = BNS_VERT_TYPE_C_GROUP;
        else
            type2 = 0;

        if ( (pVert1->type & type) == type ) {
            while( 0 <= i ) {
                pEdge1 = pBNS->edge + pVert1->iedge[i];
                v2 = pEdge1->neighbor12 ^ v1;
                if ( pBNS->vert[v2].type == type2 ) {
                    if ( IS_FORBIDDEN(pEdge1->forbidden, pBNS) ) {
                        return NO_VERTEX;
                    }
                    return v2;
                }
                i --;
            }
        }
        return BNS_BOND_ERR; /* not found t-group */

    } else
    if ( v1 < pBNS->num_vertices ) {
        return NO_VERTEX;
    }
    return BNS_VERT_EDGE_OVFL;
}
/**********************************************************************************/
int bAddStCapToAVertex( BN_STRUCT *pBNS, Vertex v1, Vertex v2, VertexFlow *nOldCapVertSingleBond, int *nDots, int bAdjacentDonors )
{
        BNS_VERTEX *pVert1   = pBNS->vert + v1;
        BNS_VERTEX *pVert;
        BNS_EDGE   *pEdge;
        Vertex      v;
        int         i, n;
        VertexFlow  nNewCap;
        /* Change v1: increment its st-cap */
        n = 0;
        nOldCapVertSingleBond[n++] = pVert1->st_edge.cap;
        /*if ( pVert1->st_edge.cap == pVert1->st_edge.flow ) {*/
            pVert1->st_edge.cap ++;
            *nDots += 1;
        /*}*/
        /* increment caps of adjacent edges if 
           (1) the neighbor has st-cap != 0 and 
           (2) (edge cap==0) OR (nSumEdgeCap < pVert1->st_edge.cap && pVert->st_edge.flow > pVert1->st_edge.cap)
        */
        if ( !(pVert1->type & BNS_VERT_TYPE_ANY_GROUP) ) {
            /*
            AT_NUMB nSumEdgeCap = 0;
            for ( i = 0; i < pVert1->num_adj_edges; i ++ ) {
                pEdge = pBNS->edge + pVert1->iedge[i];
                nSumEdgeCap += pEdge->cap;
            }
            */
            /* do not increment caps of t-group or c-group edges */
            for ( i = 0; i < pVert1->num_adj_edges; i ++ ) {
                pEdge = pBNS->edge + pVert1->iedge[i];
                nOldCapVertSingleBond[n++] = pEdge->cap; /* save edge cap */
                v = pEdge->neighbor12 ^ v1;
                if ( v == v2 && !bAdjacentDonors ) {
                    continue;
                }
                pVert = pBNS->vert + v;
                if ( pVert->type & BNS_VERT_TYPE_ANY_GROUP )
                    continue;
                nNewCap = inchi_min(pVert->st_edge.cap, pVert1->st_edge.cap);
                nNewCap = inchi_min(nNewCap, MAX_BOND_EDGE_CAP);
                pEdge->cap = nNewCap; /* change edge cap */
                /*
                if ( pVert->st_edge.cap > 0 && !pEdge->cap ) {
                    pEdge->cap ++;
                } else
                if ( pVert->st_edge.flow > pVert1->st_edge.cap &&
                     pEdge->cap < MAX_BOND_EDGE_CAP &&
                     nSumEdgeCap < pVert1->st_edge.cap ) {
                    pEdge->cap ++;
                }
                */
            }
        }

        return n; /* number of elements in nOldCapVertSingleBond[*] */
}
/**********************************************************************************/

#define BNS_CHK_ALTP_NO_ALTPATH  0
#define BNS_CHK_ALTP_SAME_TGROUP 1
#define BNS_CHK_ALTP_SAME_VERTEX 2
#define BNS_CHK_ALTP_SET_SUCCESS 4
/**********************************************************************************/
int bSetBnsToCheckAltPath( BN_STRUCT *pBNS, int nVertDoubleBond, int nVertSingleBond, AT_NUMB type, 
                           int path_type, ALT_PATH_CHANGES *apc, BNS_FLOW_CHANGES *fcd, int *nDots )
{
    
    if ( !pBNS->vert[nVertDoubleBond].st_edge.flow &&
        
         !( path_type == ALT_PATH_MODE_REM2H_CHG ||
            path_type == ALT_PATH_MODE_ADD2H_CHG ||
            path_type == ALT_PATH_MODE_REM2H_TST ||
            path_type == ALT_PATH_MODE_ADD2H_TST )
        ) {
        return BNS_CHK_ALTP_NO_ALTPATH;
    } else {
        
        
        Vertex      vNew;
        Vertex      v1       = nVertSingleBond;
        Vertex      v2       = nVertDoubleBond;

        BNS_VERTEX *pVert1   = pBNS->vert + v1;
        BNS_VERTEX *pVert2   = pBNS->vert + v2;
        int n, bAdjacentDonors = 0;
        int ifcd = 0;
        
        Vertex     t1=NO_VERTEX;
        Vertex     t2=NO_VERTEX;
        int        iapc;

/*#if( TEST_REMOVE_S_ATOMS == 1 )*/ /* && ALT_PATH_MODE_4_SALT == path_type */
        if ( ( *pBNS->pbTautFlags & TG_FLAG_TEST_TAUT2_SALTS ) &&

             ALT_PATH_MODE_4_SALT2   == path_type &&
             (BNS_VERT_TYPE_ENDPOINT & type)      ) {

/*
---------------------------------------------------------
     \   action |  DB action (v2)   |   SB action (v1)  |
vertex \        | accept H @ vertex | donate H @ vertex |
type     \      | nVertDoubleBond   | nVertSingleBond   |
----------------+-------------------+-------------------+                 
    -ZH (v1)    |  error            |   -ZH(.)          |
(cap>0 on edge  |                   |   increment       |
 except v1-v2)  |                   |   st-cap on Z     |
----------------+-------------------+-------------------+                 
    =Z  (v2)    |  =Z-(.)           |   error           |
  (st-flow>0)   |  add fict vertex  |                   |
                |  with st-cap=1    |                   |
----------------+-------------------+-------------------+                 
  endpoint      |  T(.)             |   T-(.)           |
  of t-group    |  increment        |   add fict vertex |
  represented   |  st-cap on T      |   with st-cap=1   |
  by fictitious |                   |                   |
  vertex T      |                   |                   |
---------------------------------------------------------              
*/

            int         bSet_v1;  /* indicator: v1 has been set */
            int         bSet_v2;  /* indicator: v2 has been set */
            int         i;

            Vertex      v1t      = NO_VERTEX;
            Vertex      v2t      = NO_VERTEX;
            Vertex      v1Act, v2Act;
            Vertex      v;

            memset( apc, 0, sizeof(*apc) );
            fcd[ifcd].iedge = NO_VERTEX;
            *nDots = 0;

            if ( v1 == v2 ) {
                return BNS_CHK_ALTP_SAME_VERTEX;
            }

            /* check whether v1 has neighbors adjacent to 
               multiple bonds
            */
            for ( i = 0, n = 0; i < pVert1->num_adj_edges; i ++ ) {
                v = (pBNS->edge + pVert1->iedge[i])->neighbor12 ^ v1; /* v is adjacent to v1 */
                if ( v == v2 )
                    continue; /* ignore connection to v2 */
                n += (pBNS->vert[v].st_edge.cap > 0);
            }
            if ( !n ) {
                return BNS_CHK_ALTP_NO_ALTPATH; /* the vertex cannot have flow */
            }

            v1Act = v1;
            v2Act = v2;
        
            /* find t-group that contains v1 */
            if ( (pVert1->type & type) == type ) {
                v1t = GetGroupVertex(pBNS, v1, type);
                if ( IS_BNS_ERROR( v1t ) ) {
                    return v1t;
                }
                if ( v1t != NO_VERTEX ) {
                    v1Act = v1t;
                }
            }
            /* find t-group that contains v2 */
            if ( (pVert2->type & type) == type ) {
                v2t = GetGroupVertex(pBNS, v2, type);
                if ( IS_BNS_ERROR( v2t ) ) {
                    return v2t;
                }
                if ( v2t != NO_VERTEX ) {
                    v2Act = v2t;
                }
            }
            if ( v1t != NO_VERTEX && v1t == v2t ) {
                return BNS_CHK_ALTP_SAME_TGROUP;
            }
        
            bSet_v1 = bSet_v2 = 0;
            /* create new edges adjacent to v1t or v2 */
            iapc = 0;
            if ( v1t != NO_VERTEX ) {
                /* create new edge and vertex, connect to v1t */
                vNew = bAddNewVertex( pBNS, v1t, 1, 0, 1, nDots );
                if ( IS_BNS_ERROR(vNew) ) {
                    return vNew;
                }
                apc->vNewVertex[iapc] = vNew;
                apc->bSetNew[iapc]    = 1;
                bSet_v1 = 1;
                iapc ++;
            }
            if ( v2t == NO_VERTEX ) {
                /* create new edge and vertex, connect to v2 */
                vNew = bAddNewVertex( pBNS, v2, 1, 0, 1, nDots );
                if ( IS_BNS_ERROR(vNew) ) {
                    return vNew;
                }
                apc->vNewVertex[iapc] = vNew;
                apc->bSetNew[iapc]    = 1;
                bSet_v2 = 1;
                iapc ++;
            }

            /* add st-cap to v1 and/or v2t */
            iapc = 0;
            if ( !bSet_v1 ) {
                /* add st-cap to v1 */
                if ( v1t != NO_VERTEX ) {
                    return BNS_BOND_ERR;
                }
                n = bAddStCapToAVertex( pBNS, v1, v2Act, apc->nOldCapsVert[iapc], nDots, 0 );
                apc->bSetOldCapsVert[iapc] = n;
                apc->vOldVert[iapc]        = v1;
                iapc ++;
            }
            if ( !bSet_v2 ) {
                /* add st-cap to v2t */
                if ( v2t == NO_VERTEX ) {
                    return BNS_BOND_ERR;
                }
                n = bAddStCapToAVertex( pBNS, v2t, v1Act, apc->nOldCapsVert[iapc], nDots, 0 );
                apc->bSetOldCapsVert[iapc] = n;
                apc->vOldVert[iapc]        = v2t;
                iapc ++;
            }
            if ( *nDots < 0 || *nDots %2 ) {
                return BNS_SET_ALTP_ERR;
            }
            return BNS_CHK_ALTP_SET_SUCCESS;

        }
        /* ( *pBNS->pbTautFlags & TG_FLAG_TEST_TAUT2_SALTS ) */
/*#endif*/ /*  ( TEST_REMOVE_S_ATOMS == 1 && ALT_PATH_MODE_4_SALT == path_type ) */
        /*****************************************************************/
        if ( path_type == ALT_PATH_MODE_REM2H_CHG ||
             path_type == ALT_PATH_MODE_ADD2H_CHG ||
             path_type == ALT_PATH_MODE_REM2H_TST ||
             path_type == ALT_PATH_MODE_ADD2H_TST ) { /* added 2004-03-18 */
            
            int bDonors = (path_type == ALT_PATH_MODE_REM2H_CHG) || (path_type == ALT_PATH_MODE_REM2H_TST);

            int         bSet_v1;  /* indicator: v1 has been set */
            int         bSet_v2;  /* indicator: v2 has been set */
            int         i, cap = 1;

            Vertex      v1t      = NO_VERTEX;
            Vertex      v2t      = NO_VERTEX;
            Vertex      v1Act, v2Act;
            Vertex      v;

            memset( apc, 0, sizeof(*apc) );
            fcd[ifcd].iedge = NO_VERTEX;
            *nDots = 0;
            /*
            if ( v1 == v2 ) {
                return BNS_CHK_ALTP_SAME_VERTEX;
            }
            */
            /* check whether v1 and v2 have proper neighbors  */
            for ( i = 0, n = bAdjacentDonors = 0; i < pVert1->num_adj_edges; i ++ ) {
                v = (pBNS->edge + pVert1->iedge[i])->neighbor12 ^ v1; /* v is adjacent to v1 */
                /* do not ignore connection to v2
                if ( v == v2 )
                    continue;
                */
                n += bDonors ? (pBNS->vert[v].st_edge.cap > 0) : ((pBNS->edge + pVert1->iedge[i])->flow > 0);
                bAdjacentDonors += bDonors ? (v == v2) && ((pBNS->edge + pVert1->iedge[i])->flow < MAX_BOND_EDGE_CAP) : 0;
                     /* two donors connected by a single or double bond */ 
            }
            if ( !n && !bAdjacentDonors ) {
                return BNS_CHK_ALTP_NO_ALTPATH; /* the vertex cannot have flow */
            }
            for ( i = 0, n = bAdjacentDonors = 0; i < pVert2->num_adj_edges; i ++ ) {
                v = (pBNS->edge + pVert2->iedge[i])->neighbor12 ^ v2; /* v is adjacent to v2 */
                /* do not ignore connection to v1
                if ( v == v1 )
                    continue;
                 */
                n += bDonors ? (pBNS->vert[v].st_edge.cap > 0) : ((pBNS->edge + pVert2->iedge[i])->flow > 0);
                bAdjacentDonors += bDonors ? (v == v1) && ((pBNS->edge + pVert2->iedge[i])->flow < MAX_BOND_EDGE_CAP ) : 0;
                     /* two donors connected by a single or double bond */
            }
            if ( !n && !bAdjacentDonors ) {
                return BNS_CHK_ALTP_NO_ALTPATH; /* the vertex cannot have flow */
            }

            v1Act = v1;
            v2Act = v2;
        
            /* find t-group that contains v1 */
            if ( (pVert1->type & type) == type ) {
                v1t = GetGroupVertex(pBNS, v1, type);
                if ( BNS_BOND_ERR == v1t ) {
                    v1t = NO_VERTEX;
                } else
                if ( IS_BNS_ERROR( v1t ) ) {
                    return v1t;
                } else
                if ( v1t != NO_VERTEX ) {
                    v1Act = v1t;
                }
            }
            /* find t-group that contains v2 */
            if ( (pVert2->type & type) == type ) {
                v2t = GetGroupVertex(pBNS, v2, type);
                if ( BNS_BOND_ERR == v2t ) {
                    v2t = NO_VERTEX;
                } else
                if ( IS_BNS_ERROR( v2t ) ) {
                    return v2t;
                } else
                if ( v2t != NO_VERTEX ) {
                    v2Act = v2t;
                }
            }
            
            if ( v1t != NO_VERTEX && v1t == v2t ) {
                cap = 2; /* same t-group */
            }
            
            /*  bAddNewVertex: (bDonors != 0) == (vit != NO_VERTEX), i=1,2 */
            bSet_v1 = bSet_v2 = 0;
            /* create new edges adjacent to v1t or v2 */
            iapc = 0;
            if ( (bDonors != 0) == (v1t != NO_VERTEX) ) {
                /* create new edge and vertex, connect to v1Act */
                vNew = bAddNewVertex( pBNS, v1Act, cap, 0, 1, nDots );
                if ( IS_BNS_ERROR(vNew) ) {
                    return vNew;
                }
                apc->vNewVertex[iapc] = vNew;
                apc->bSetNew[iapc]    = 1;
                bSet_v1 = 1;
                iapc ++;
            }
            if (  (bDonors != 0) == (v2t != NO_VERTEX) && cap == 1 ) {
                /* create new edge and vertex, connect to v2Act; do not do it if cap==2 */
                vNew = bAddNewVertex( pBNS, v2Act, cap, 0, 1, nDots );
                if ( IS_BNS_ERROR(vNew) ) {
                    return vNew;
                }
                apc->vNewVertex[iapc] = vNew;
                apc->bSetNew[iapc]    = 1;
                bSet_v2 = 1;
                iapc ++;
            } else
            if (  (bDonors != 0) == (v2t != NO_VERTEX) ) {
                bSet_v2 = 1;
            }

            /* add st-cap to v1 and/or v2t */
            iapc = 0;
            /* if cap=2 then just increment st_cap 2 times */
            if ( !bSet_v1 ) {
                /* add st-cap to v1 */
                if ( (bDonors != 0) == (v1t != NO_VERTEX) ) {
                    return BNS_BOND_ERR;
                }
                n = bAddStCapToAVertex( pBNS, v1Act, v2Act, apc->nOldCapsVert[iapc], nDots, bAdjacentDonors );
                apc->bSetOldCapsVert[iapc] = n;
                apc->vOldVert[iapc]        = v1Act;
                iapc ++;
            }
            if ( !bSet_v2 ) {
                /* add st-cap to v2t */
                if ( (bDonors != 0) == (v2t != NO_VERTEX) ) {
                    return BNS_BOND_ERR;
                }
                n = bAddStCapToAVertex( pBNS, v2Act, v1Act, apc->nOldCapsVert[iapc], nDots, bAdjacentDonors );
                apc->bSetOldCapsVert[iapc] = n;
                apc->vOldVert[iapc]        = v2Act;
                iapc ++;
            }
            if ( *nDots < 0 || *nDots %2 ) {
                return BNS_SET_ALTP_ERR;
            }
            return BNS_CHK_ALTP_SET_SUCCESS;

        }
        /**************************************************************************/
        if ( path_type == ALT_PATH_MODE_REM_PROTON ) { /* added 2004-03-05 */
            if ( v1 >= 0 && v2 >= 0 && 
                 (pVert1->type & BNS_VERT_TYPE_ANY_GROUP) &&
                 (pVert2->type & BNS_VERT_TYPE_ANY_GROUP) ) {
                /* create new edge and vertex, connect to v2 */
                if ( (pBNS->vert[v1].type & BNS_VERT_TYPE_C_GROUP) &&
                     (pBNS->vert[v1].st_edge.flow == 2*pBNS->vert[v1].num_adj_edges ) ) {
                    /* so far in a charge group max edge flow = 1 2004-03-08 */
                    return BNS_CHK_ALTP_NO_ALTPATH;
                }
                memset( apc, 0, sizeof(*apc) );
                fcd[ifcd].iedge = NO_VERTEX;
                *nDots = 0;
                iapc = 0;

                vNew = bAddNewVertex( pBNS, v2, 1, 0, 1, nDots );
                if ( IS_BNS_ERROR(vNew) ) {
                    return vNew;
                }
                apc->vNewVertex[iapc] = vNew;
                apc->bSetNew[iapc]    = 1;
                /*iapc ++;*/
                /* add st-cap (dot) to v1 */
                n = bAddStCapToAVertex( pBNS, v1, v2, apc->nOldCapsVert[iapc], nDots, 0 );
                apc->bSetOldCapsVert[iapc] = n;
                apc->vOldVert[iapc]        = v1;
                iapc ++;
                return BNS_CHK_ALTP_SET_SUCCESS;
            }
        }

#if ( NEUTRALIZE_ENDPOINTS == 1 ) /* { */


        *nDots = 0;
        memset( apc, 0, sizeof(*apc) );
        fcd[ifcd].iedge = NO_VERTEX;

        if ( type & BNS_VERT_TYPE_ENDPOINT ) {
            BNS_IEDGE  iedge;
            AT_NUMB    type2;
            int ret2;
            /* prohibit charge movement */
            type2 =  BNS_VERT_TYPE_C_GROUP;
            iedge = GetEdgeToGroupVertex( pBNS, v1, type2 );
            if (iedge != NO_VERTEX ) {
                /*  set flow=1 on an edge to a c-group vertex to make sure there is no positive charge
                 *  when moving tautomeric H-atoms
                 */
                ret2 = bSetFlowToCheckOneBond( pBNS, iedge, 1, fcd+ifcd );
                if ( IS_BNS_ERROR(ret2) ) {
                    return ret2;
                }
                *nDots += ret2;
                while ( fcd[ifcd].iedge != NO_VERTEX ) {
                    ifcd ++;
                }
            }
            iedge = GetEdgeToGroupVertex( pBNS, v2, type2 );
            if (iedge != NO_VERTEX ) {
                /*  set flow=1 on an edge to a c-group vertex to make sure there is no positive charge
                 *  when moving tautomeric H-atoms
                 */
                ret2 = bSetFlowToCheckOneBond( pBNS, iedge, 1, fcd+ifcd );
                if ( IS_BNS_ERROR(ret2) ) {
                    return ret2;
                }
                *nDots += ret2;
                while ( fcd[ifcd].iedge != NO_VERTEX ) {
                    ifcd ++;
                }
            }
            /* set hydrogen counts */
            type2 = BNS_VERT_TYPE_TGROUP;
            iedge = GetEdgeToGroupVertex( pBNS, v1, type2 );
            if (iedge != NO_VERTEX ) {
                /*  set flow=1 on an edge to a t-group vertex to make sure there is
                 *  a moveable hydrogen atom or (-) on v1 when moving tautomeric H-atoms
                 */
#if( FIX_H_CHECKING_TAUT == 1 )
                ret2 = bSetFlowToCheckOneBond( pBNS, iedge, 1, fcd+ifcd );
                if ( IS_BNS_ERROR(ret2) ) {
                    return ret2;
                }
                *nDots += ret2;
                while ( fcd[ifcd].iedge != NO_VERTEX ) {
                    ifcd ++;
                }
#else
                t1 = pBNS->edge[iedge].neighbor12 ^ v1;
#endif
            }
            iedge = GetEdgeToGroupVertex( pBNS, v2, type2 );
            if (iedge != NO_VERTEX ) {
                /*  set flow=0 on an edge to a t-group vertex to make sure there is
                 *  no moveable hydrogen atom or (-) on v2 when moving tautomeric H-atoms
                 */
#if( FIX_H_CHECKING_TAUT == 1 )
                ret2 = bSetFlowToCheckOneBond( pBNS, iedge, 0, fcd+ifcd );
                if ( IS_BNS_ERROR(ret2) ) {
                    return ret2;
                }
                *nDots += ret2;
                while ( fcd[ifcd].iedge != NO_VERTEX ) {
                    ifcd ++;
                }
#else
                t2 = pBNS->edge[iedge].neighbor12 ^ v2;
#endif
            }

#if( FIX_H_CHECKING_TAUT == 1 )
#else
            if ( t1 == t2 && t1 != NO_VERTEX ) {
                return BNS_CHK_ALTP_SAME_TGROUP;
            }
#endif
            iapc = 0;
            /* create new edge and vertex with cap=1 at v2 and/or t1 */
            if ( t1 != NO_VERTEX ) {
                /* create new edge and vertex, connect to t1 */
                vNew = bAddNewVertex( pBNS, t1, 1/*cap*/, 0/*flow*/, 1/*max_adj_edges*/, nDots );
                if ( IS_BNS_ERROR(vNew) ) {
                    return vNew;
                }
                apc->vNewVertex[iapc] = vNew;
                apc->bSetNew[iapc]    = 1;
                iapc ++;
            }
            if ( t2 == NO_VERTEX ) {
                /* create new edge and vertex, connect to v2 */
                vNew = bAddNewVertex( pBNS, v2, 1/*cap*/, 0/*flow*/, 1/*max_adj_edges*/, nDots );
                if ( IS_BNS_ERROR(vNew) ) {
                    return vNew;
                }
                apc->vNewVertex[iapc] = vNew;
                apc->bSetNew[iapc]    = 1;
                iapc ++;
            }

            /* add st-cap to v1 and/or v2t */
            iapc = 0;
            if ( t1 == NO_VERTEX ) {
                /* add st-cap to v1 */
                n = bAddStCapToAVertex( pBNS, v1, (Vertex)(t2 == NO_VERTEX? v2:t2), apc->nOldCapsVert[iapc], nDots, 0 );
                apc->bSetOldCapsVert[iapc] = n;
                apc->vOldVert[iapc]        = v1;
                iapc ++;
            }
            if ( t2 != NO_VERTEX ) {
                /* add st-cap to t2 */
                n = bAddStCapToAVertex( pBNS, t2, (Vertex)(t1 == NO_VERTEX? v1:t1), apc->nOldCapsVert[iapc], nDots, 0 );
                apc->bSetOldCapsVert[iapc] = n;
                apc->vOldVert[iapc]        = t2;
                iapc ++;
            }
        } else {
            /* create new edge and vertex, connect to v2 */
            vNew = bAddNewVertex( pBNS, v2, 1 /* cap*/, 0 /* flow */, 1 /* max_adj_edges */, nDots );
            if ( IS_BNS_ERROR(vNew) ) {
                return vNew;
            }
            apc->vNewVertex[0] = vNew;
            apc->bSetNew[0]    = 1;

            /* add st-cap to v1 */
            n = bAddStCapToAVertex( pBNS, v1, v2, apc->nOldCapsVert[0], nDots, 0 );
            apc->bSetOldCapsVert[0] = n;
            apc->vOldVert[0]        = v1;
        }
#else  /* } NEUTRALIZE_ENDPOINTS == 0 {*/

        *nDots = 0;
        memset( apc, 0, sizeof(*apc) );
        fcd[ifcd].iedge = NO_VERTEX;

        /* create new edge and vertex, connect to v2 */
        vNew = bAddNewVertex( pBNS, v2, 1 /* cap*/, 0 /* flow */, 1 /* max_adj_edges */, nDots, 0 );
        if ( IS_BNS_ERROR(vNew) ) {
            return vNew;
        }
        apc->vNewVertex[0] = vNew;
        apc->bSetNew[0]    = 1;

        /* add st-cap to v1 */
        n = bAddStCapToAVertex( pBNS, v1, v2, apc->nOldCapsVert[0], nDots );
        apc->bSetOldCapsVert[0] = n;
        apc->vOldVert[0]        = v1;
#endif /* } NEUTRALIZE_ENDPOINTS */

        if ( *nDots < 0 || *nDots %2 ) {
            return BNS_SET_ALTP_ERR;
        }
        return BNS_CHK_ALTP_SET_SUCCESS;
    }
    /*return BNS_CHK_ALTP_NO_ALTPATH;*/
}
/**********************************************************************************/
int bRestoreBnsAfterCheckAltPath( BN_STRUCT *pBNS, ALT_PATH_CHANGES *apc, int bChangeFlow )
/* int nVertDoubleBond, int nVertSingleBond, int nNewVertex, AT_NUMB *nOldCapVertSingleBond */
{
    BNS_EDGE   *pEdge;
    Vertex      vNew;
    Vertex      vOld;
    BNS_VERTEX *pOldVert;
    BNS_VERTEX *pNewVert;
    int i, j, n, ret;

    ret = 0;
    if ( bChangeFlow & BNS_EF_UPD_H_CHARGE ) {
        /* remove new temp. vertices and edges connectong them to the structure */
        for ( i = sizeof(apc->bSetNew)/sizeof(apc->bSetNew[0])-1; 0 <= i; i -- ) {
            if ( apc->bSetNew[i] ) {
                vNew = apc->vNewVertex[i];
                pNewVert = pBNS->vert + vNew;
                for ( j = 0; j < pNewVert->num_adj_edges; j ++ ) {
                    pEdge    = pBNS->edge+pNewVert->iedge[j];
                    vOld     = pEdge->neighbor12 ^ vNew;
                    pOldVert = pBNS->vert + vOld;
                    pOldVert->st_edge.flow -= pEdge->flow;
                    pOldVert->st_edge.cap  -= pEdge->flow;
                    /* disconnect new edge from pOldVert */
                    pOldVert->iedge[--pOldVert->num_adj_edges] = 0;
                    /* clear the new edge */
                    memset( pEdge, 0, sizeof(*pEdge) );
                    /* and decrement the total number of edges */
                    pBNS->num_edges --;
                }
                /* clear the new vertex */
                memset( pNewVert, 0, sizeof( pNewVert ) );
                /* and decrement the total number of vertices (new vertice ids are contiguous */
                pBNS->num_vertices --;
                ret ++;
            }
        }
        /* Restore changed caps of old vertices */
        for ( i = sizeof(apc->bSetOldCapsVert)/sizeof(apc->bSetOldCapsVert[0])-1; 0 <= i; i -- ) {
            if ( n = apc->bSetOldCapsVert[i] ) {
                pOldVert   = pBNS->vert + apc->vOldVert[i];
                if ( pOldVert->st_edge.flow <= apc->nOldCapsVert[i][0] ) { 
                    pOldVert->st_edge.cap = apc->nOldCapsVert[i][0];
                    n --;
                    ret ++;
                    for ( j = 0; j < n && j < pOldVert->num_adj_edges; j ++ ) {
                        pEdge = pBNS->edge + pOldVert->iedge[j];
                        pEdge->cap = apc->nOldCapsVert[i][j+1];
                    }
                }
            }
        }
     } else {
        /* Restore changed caps of old vertices */
        for ( i = sizeof(apc->bSetOldCapsVert)/sizeof(apc->bSetOldCapsVert[0])-1; 0 <= i; i -- ) {
            if ( n = apc->bSetOldCapsVert[i] ) {
                pOldVert   = pBNS->vert + apc->vOldVert[i];
                pOldVert->st_edge.cap = apc->nOldCapsVert[i][0];
                n --;
                ret ++;
                for ( j = 0; j < n && j < pOldVert->num_adj_edges; j ++ ) {
                    pEdge = pBNS->edge + pOldVert->iedge[j];
                    pEdge->cap = apc->nOldCapsVert[i][j+1];
                }
            }
        }

        /* remove new temp. vertices and edges connectong them to the structure */
        for ( i = sizeof(apc->bSetNew)/sizeof(apc->bSetNew[0])-1; 0 <= i; i -- ) {
            if ( apc->bSetNew[i] ) {
                vNew = apc->vNewVertex[i];
                pNewVert = pBNS->vert + vNew;
                for ( j = 0; j < pNewVert->num_adj_edges; j ++ ) {
                    pEdge    = pBNS->edge+pNewVert->iedge[j];
                    vOld     = pEdge->neighbor12 ^ vNew;
                    pOldVert   = pBNS->vert + vOld;
                    /* disconnect new edge from pOldVert */
                    pOldVert->iedge[--pOldVert->num_adj_edges] = 0;
                    /* clear the new edge */
                    memset( pEdge, 0, sizeof(*pEdge) );
                    /* and decrement the total number of edges */
                    pBNS->num_edges --;
                }
                /* clear the new vertex */
                memset( pNewVert, 0, sizeof( pNewVert ) );
                /* and decrement the total number of vertices (new vertice ids are contiguous */
                pBNS->num_vertices --;
                ret ++;
            }
        }
    }
    return 0;
}
/**********************************************************************************/
int bExistsAnyAltPath( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms,
                      int nVert2, int nVert1, int path_type )
{
    int nRet1, nRet2;
    nRet1 = bExistsAltPath( pBNS, pBD, NULL, at, num_atoms, nVert2, nVert1, path_type );
    if ( nRet1 > 0 )
        return nRet1;
    nRet2 = bExistsAltPath( pBNS, pBD, NULL, at, num_atoms, nVert1, nVert2, path_type );
    if ( nRet2 > 0 )
        return nRet2;
    if ( IS_BNS_ERROR( nRet1 ) )
        return nRet1;
    if ( IS_BNS_ERROR( nRet2 ) )
        return nRet2;
    return 0;
}

#define ALT_PATH_TAUTOM 1
#define ALT_PATH_CHARGE 2
#define ALT_PATH_4_SALT 3



/**********************************************************************************/
int bIsBnsEndpoint( BN_STRUCT *pBNS, int v )
{
    int i, vt;
    BNS_VERTEX    *pVert;  /* vertices */
    BNS_EDGE      *pEdge;  /* edges */

    if ( 0 <= v && v < pBNS->num_atoms && (pVert = pBNS->vert+v) && (pVert->type & BNS_VERT_TYPE_ENDPOINT) ) {
        for ( i = pVert->num_adj_edges - 1; 0 <= i; i -- ) {
            pEdge = pBNS->edge + pVert->iedge[i];
            vt = pEdge->neighbor12 ^ v;
            if ( pBNS->vert[vt].type & BNS_VERT_TYPE_TGROUP ) {
                return !IS_FORBIDDEN(pEdge->forbidden, pBNS);
            }
        }
    }
    return 0;
}
/**********************************************************************************/
#if ( BNS_RAD_SEARCH == 1 )
/**********************************************************************************/
int bRadChangesAtomType( BN_STRUCT *pBNS, BN_DATA *pBD, Vertex v, Vertex v_1, Vertex v_2 )
{
    
    EdgeIndex iuv;
    Vertex v_O, v_ChgOrH;
    /* the previous atom along the path: should be a terminal atom */
    if ( v_1 == NO_VERTEX ) {
        v_1 = GetPrevVertex( pBNS, v, pBD->SwitchEdge, &iuv );
    }
    v_O = v_1 / 2 - 1;
    if ( v_O < 0 || v_O >= pBNS->num_atoms ) {
        return 0;
    }
    /* make sure v_O is a terminal atom: its second neighbor is not an atom */
    if ( pBNS->vert[pBNS->edge[pBNS->vert[v_O].iedge[1]].neighbor12 ^ v_O].type & BNS_VERT_TYPE_ATOM ) {
        return 0;
    }
    /* the next to previous vertex vertex along the path: should be a Charge or Taut group vertex */
    if ( v_2 == NO_VERTEX ) {
        v_2 = GetPrevVertex( pBNS, v_1, pBD->SwitchEdge, &iuv );
    }
    v_ChgOrH = v_2 / 2 - 1;
    if ( v_ChgOrH < pBNS->num_atoms ) {
        return 0;
    }
    /* make sure v_ChgOrH is a charge or taut_group */
    if ( pBNS->vert[v_ChgOrH].type & (BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_C_GROUP) )
        return 1;
    return 0;
}
/**********************************************************************************/
int RegisterRadEndpoint( BN_STRUCT *pBNS, BN_DATA *pBD, Vertex u)
{
    EdgeIndex iuv;
    int       i, num_found;
    Vertex    v, w;
    Vertex    u_last, v2;
    switch( pBD->bRadSrchMode ) {
    case  RAD_SRCH_NORM:
        /* go backwards along alt path and stop at the 1st found atom (not a fictitious vertex) */
        /* we need only vertices where a radical may be moved, therefore exclude u%2=1 (odd) vertices */
        /* atom number = u/2-1; u = 0 or 1 is 's' or 't' vertices, respectively, they are not atoms  */
        num_found = 0;
        while ( u > Vertex_t && (u % 2 || u/2 > pBNS->num_atoms ) ) {
            u = GetPrevVertex( pBNS, u, pBD->SwitchEdge, &iuv );
        }
        w = u/2 - 1; /* Check whether u is a radical endpoint */
        if ( Vertex_t < u && w < pBNS->num_atoms &&
             pBNS->vert[w].st_edge.cap == (pBNS->vert[w].st_edge.flow & EDGE_FLOW_ST_MASK) ) {
            /* u is an atom; it is not a radical atom */
            /* now search for the starting radical atom by following the path back from u */
            v = u_last = u;
            while( v > Vertex_t ) {
                u = v;
                v = GetPrevVertex( pBNS, u, pBD->SwitchEdge, &iuv ); /* Radical endpoint */
            }
            /* check whether u is a radical atom */
            if ( !(u%2) && Vertex_t < u &&
                 (u = u/2 - 1) < pBNS->num_atoms &&
                 pBNS->vert[u].st_edge.cap > (pBNS->vert[u].st_edge.flow & EDGE_FLOW_ST_MASK) ) {
                /* at pBNS->vert[u] we have found the radical that originated the path */
                /* pBD->RadEndpoints[2k] is the radical, pBD->RadEndpoints[2k+1] is the farthest atom */
                /* to which the radical may be moved (farthest reachable atom) */

                /* add *all* atoms that may receive radical from u_rad */
                /* exception: at2 in: ==(+/-/H)---at1==at2(possible rad endpoint) if pBNS->type_TACN */
                
                for ( v = u_last; v > Vertex_t; v = GetPrevVertex( pBNS, v, pBD->SwitchEdge, &iuv ) ) {
                    if ( !(v%2) && (v2 = v/2 - 1) < pBNS->num_atoms &&
                         pBNS->vert[v2].st_edge.cap == (pBNS->vert[v2].st_edge.flow & EDGE_FLOW_ST_MASK) ) {
                        /* check exception */
                        if ( pBNS->type_TACN &&
                             bRadChangesAtomType( pBNS, pBD, v, NO_VERTEX, NO_VERTEX ) ) {
                            continue;
                        }
                        /* add */
                        for ( i = 0; i < pBD->nNumRadEndpoints; i += 2 ) {
                            /* check whether this pair, (u,w), has already been saved */
                            if ( u  == pBD->RadEndpoints[i] &&
                                 v2 == pBD->RadEndpoints[i+1] ) {
                                break;
                            }
                        }
                        if ( i >= pBD->nNumRadEndpoints ) {
                            /* add new (u,w) pair */
                            if ( pBD->nNumRadEndpoints+2 <= pBD->max_num_vertices  ) {
                                /* add */
                                pBD->RadEndpoints[pBD->nNumRadEndpoints ++] = u; /* radical */
                                pBD->RadEndpoints[pBD->nNumRadEndpoints ++] = v2; /* endpoint */
                                num_found ++;
                                /*return 1;*/ /* registered */
                            } else {
                                return BNS_VERT_EDGE_OVFL;
                            }
                        }
                    }
                }
                if ( num_found ) {
                    return 1;
                }
            }
        }
        break;

    case RAD_SRCH_FROM_FICT:
        /* find nearest atom accessible from a fictitious vertex */
        /* go backwards along alt path and stop at the 1st found atom (not a fictitious vertex) */
        v = u;
        w = NO_VERTEX; /* the nearest atom -- radical-endpoint */
        u = NO_VERTEX; /* fictitious vertex carrying a radical */
        while ( v > Vertex_t ) {
            u = v;
            if ( !(v % 2) && v/2 <= pBNS->num_atoms && 
                 pBNS->vert[v/2-1].st_edge.cap - pBNS->vert[v/2-1].st_edge.flow < 2 ) {
                w = v; /* vertex w is atom that may be singlet or doublet but not triplet */
            }
            v = GetPrevVertex( pBNS, u, pBD->SwitchEdge, &iuv );
        }
        v = u/2 - 1; /* vertex u may be the radical from which the path originated; w is the nearest atom */
        if ( w == NO_VERTEX || u == NO_VERTEX || w % 2 || u == w || v < pBNS->num_atoms ||
             pBNS->vert[v].st_edge.cap == pBNS->vert[v].st_edge.flow ||
             (w = w/2 - 1) >= pBNS->num_atoms ) {
            break; /* reject */
        }
        u = v;
        /* at pBNS->vert[u] we have found the radical that originated the path, w is the nearest atom */
        for ( i = 0; i < pBD->nNumRadEndpoints; i += 2 ) {
            if ( u == pBD->RadEndpoints[i] &&
                 w == pBD->RadEndpoints[i+1] ) {
                break; /* this pair has already been stored */
            }
        }
        if ( i >= pBD->nNumRadEndpoints ) {
            /* a new pair has been found */
            if ( pBD->nNumRadEndpoints+2 <= pBD->max_num_vertices  ) {
                /* add */
                pBD->RadEndpoints[pBD->nNumRadEndpoints ++] = u; /* radical */
                pBD->RadEndpoints[pBD->nNumRadEndpoints ++] = w; /* endpoint */
                return 1; /* registered */
            } else {
                return BNS_VERT_EDGE_OVFL;
            }
        }
        break;
    }

    return 0; /* rejected */
}
/**********************************************************************************/
int cmp_rad_endpoints( const void *a1, const void *a2 )
{
    /* Vertex radical_vertex, radical_endpoint */
    const Vertex *p1 = (const Vertex *)a1;
    const Vertex *p2 = (const Vertex *)a2;
    if ( p1[0] < p2[0] )
        return -1;
    if ( p1[0] > p2[0] )
        return 1;
    if ( p1[1] < p2[1] )
        return -1;
    if ( p1[1] > p2[1] )
        return 1;
    return 0;
}
/**********************************************************************************/
int RemoveRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at )
{
    BNS_EDGE   *e;
    EdgeIndex   ie;
    BNS_VERTEX *p1, *p2;
    Vertex      v1, v2;
    int         i, delta, rad;
    for ( i = pBD->nNumRadEdges-1; 0 <= i; i -- ) {
        ie = pBD->RadEdges[i];
        if ( ie < 0 || ie >= pBNS->num_edges ) {
            goto error_exit;
        }
        e = pBNS->edge + ie;
        v1 = e->neighbor1;
        v2 = e->neighbor12 ^ v1;   /* v2 > v1 <=> v2 was added later */
        if ( ie + 1 != pBNS->num_edges || 
             v1 < 0 || v1 >= pBNS->num_vertices ||
             v2 < 0 || v2 >= pBNS->num_vertices ) {
            goto error_exit;
        }
        p1 = pBNS->vert + v1;
        p2 = pBNS->vert + v2;

        if ( p2->iedge[p2->num_adj_edges-1] != ie ||
             p1->iedge[p1->num_adj_edges-1] != ie ) {
            goto error_exit;
        }
        p2->num_adj_edges --;
        p1->num_adj_edges --;
        p2->iedge[p2->num_adj_edges] = 0;
        p1->iedge[p1->num_adj_edges] = 0;
        p2->st_edge.flow -= e->flow;
        p1->st_edge.flow -= e->flow;

        if ( !p2->num_adj_edges && v2 >= pBNS->num_atoms ) {
            if ( v2+1 != pBNS->num_vertices ) {
                goto error_exit;
            }
            memset( p2, 0, sizeof(*p2) );
            pBNS->num_vertices --;
        }
        if ( !p1->num_adj_edges && v1 >= pBNS->num_atoms ) {
            if ( v1+1 != pBNS->num_vertices ) {
                goto error_exit;
            }
            memset( p1, 0, sizeof(*p1) );
            pBNS->num_vertices --;
        }
        if ( at && v1 < pBNS->num_atoms ) {
            delta = p1->st_edge.cap - p1->st_edge.flow;
            rad   = at[v1].radical;
            switch( delta ) {
            case 0:
                if ( rad == RADICAL_DOUBLET )
                    rad = 0;
                break;
            case 1:
                if ( rad != RADICAL_DOUBLET )
                    rad = RADICAL_DOUBLET;
            }
            at[v1].radical = rad;
        }
        memset( e, 0, sizeof(*e) );
        pBNS->num_edges --;
    }
    pBD->nNumRadEdges = 0;
    pBD->nNumRadicals = 0;
    pBD->bRadSrchMode = RAD_SRCH_NORM;
    return 0;
error_exit:
    return BNS_PROGRAM_ERR;
}
/**********************************************************************************/
int RestoreRadicalsOnly( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at )
{
    BNS_EDGE   *e;
    EdgeIndex   ie;
    BNS_VERTEX *p1, *p2;
    Vertex      v1, v2;
    int         i, delta, rad;
    int         p1_num_adj_edges, p2_num_adj_edges;

    for ( i = pBD->nNumRadEdges-1; 0 <= i; i -- ) {
        ie = pBD->RadEdges[i];
        if ( ie < 0 || ie >= pBNS->num_edges ) {
            goto error_exit;
        }
        e = pBNS->edge + ie;
        v1 = e->neighbor1;         /* atom */
        v2 = e->neighbor12 ^ v1;   /* v2 > v1 <=> v2 was added later */
        if ( v1 < 0 || v1 >= pBNS->num_atoms ||
             v2 < pBNS->num_atoms || v2 >= pBNS->num_vertices ) {
            goto error_exit;
        }
        p1 = pBNS->vert + v1;
        p2 = pBNS->vert + v2;

        p1_num_adj_edges = e->neigh_ord[0];
        p2_num_adj_edges = e->neigh_ord[1];

        if ( p2->iedge[p2_num_adj_edges] != ie ||
             p1->iedge[p1_num_adj_edges] != ie ) {
            goto error_exit;
        }

        if ( at && v1 < pBNS->num_atoms ) {
            delta = p1->st_edge.cap - p1->st_edge.flow + e->flow;
            rad   = at[v1].radical;
            switch( delta ) {
            case 0:
                if ( rad == RADICAL_DOUBLET )
                    rad = 0;
                break;
            case 1:
                if ( rad != RADICAL_DOUBLET )
                    rad = RADICAL_DOUBLET;
            }
            at[v1].radical = rad;
        }
    }
    return 0;
error_exit:
    return BNS_PROGRAM_ERR;
}
/**********************************************************************************/
int SetRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode )
{
    int ret, i, j, k, num_new_edges, delta;
    BNS_VERTEX *pRad, *pEndp;
    Vertex     wRad, vRad, vEndp, nNumRadicals;
    int        nDots=0 /* added initialization, 2006-03 */, nNumEdges;
    if ( pBNS->tot_st_cap <= pBNS->tot_st_flow ) {
        return 0;
    }
    pBD->nNumRadEndpoints = 0;
    pBD->nNumRadEdges     = 0;
    pBD->bRadSrchMode     = bRadSrchMode;
    pBNS->alt_path = pBNS->altp[0];
    pBNS->bChangeFlow = 0;
    ret = BalancedNetworkSearch( pBNS, pBD, BNS_EF_RAD_SRCH );
    ReInitBnData( pBD );
    ReInitBnStructAltPaths( pBNS );
    if ( !ret && pBD->nNumRadEndpoints >= 2 ) {
        /* sort by radical locations */
        qsort( pBD->RadEndpoints, pBD->nNumRadEndpoints/2, 2*sizeof(pBD->RadEndpoints[0]), cmp_rad_endpoints );
        num_new_edges = 0;
        nNumRadicals  = 0;
        /* create new vertices (type=BNS_VERT_TYPE_TEMP) and edges with flow=cap=1 */
        /* connecting the new vertices radical vertices */
        for ( i = 0; i < pBD->nNumRadEndpoints; i = j ) {
            wRad = pBD->RadEndpoints[i];
            pRad = pBNS->vert + wRad;
            delta = pRad->st_edge.cap - (pRad->st_edge.flow & EDGE_FLOW_ST_MASK);
            if ( delta <= 0 ) {
                delta = 1;
            }
            nNumEdges = 0;
            for ( j = i; j < pBD->nNumRadEndpoints && wRad == pBD->RadEndpoints[j] ; j += 2 ) {
                nNumEdges ++;
            }
            /* add new aux vertex to the radical atom/vertex */
            vRad = bAddNewVertex( pBNS, wRad, delta, delta, nNumEdges+1, &nDots );
            if ( IS_BNS_ERROR( vRad ) ) {
                ret = vRad;
                goto error_exit;
            }
            pRad     = pBNS->vert + vRad;
            pBD->RadEdges[pBD->nNumRadEdges ++] = pRad->iedge[pRad->num_adj_edges-1];
            /* replace references to vertex wRad with vRad */
            for ( k = i, nNumEdges = 0; k < j; k += 2 ) {
                pBD->RadEndpoints[k] = vRad;
            }
            nNumRadicals ++;
        }
        /* all vRad vertex indices should be in the range vFirstNewVertex...vFirstNewVertex+nNumRadicals-1 */
        /* connect new vertices to the radical endpoints thus replacing radicals with even-length alternating cycles */
        for ( i = 0; i < pBD->nNumRadEndpoints; i = j ) {
            vRad = pBD->RadEndpoints[i];
            pRad = pBNS->vert + vRad;
            for ( j = i; j < pBD->nNumRadEndpoints && vRad == pBD->RadEndpoints[j] ; j += 2 ) {
                /* connect vew vertex pRad to radical endpoints */
                vEndp    = pBD->RadEndpoints[j+1];
                pEndp    = pBNS->vert + vEndp;
                ret = AddNewEdge( pRad, pEndp, pBNS, 1, 0 );
                if ( IS_BNS_ERROR( ret ) ) {
                    goto error_exit;
                }
                pBD->RadEdges[pBD->nNumRadEdges ++] = ret;
            }
        }
        pBD->nNumRadicals = nNumRadicals;
        return nNumRadicals; /* done */
    }
    return 0; /* nothing to do */

error_exit:
    RemoveRadEndpoints( pBNS, pBD, NULL );
    return ret;

}
/**********************************************************************************/
#define MAX_NUM_RAD  256
/***************************************************************************/
int SetRadEndpoints2( BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode )
{
    int ret = 0, i, j, k, n, num_new_edges, delta = 1;
    BNS_VERTEX *pRad, *pEndp;
    Vertex     wRad, vRad, vEndp, nNumRadicals;
    Vertex     vRadList[MAX_NUM_RAD], vRadEqul[MAX_NUM_RAD];
    int        nNumRad = 0;
    int        edge_flow;
    int        nDots=0 /* added initialization, 2006-03 */, nNumEdges;
    NodeSet    VertSet;
    if ( pBNS->tot_st_cap <= pBNS->tot_st_flow ) {
        return 0;
    }
    /* find all radicals: their vertices have st_cap-st_flow=delta */
    /* save radical atom numbers in vRadList[] and remove radical by making st_cap=st_flow */
    for ( i = 0; i < pBNS->num_atoms; i ++ ) {
        if ( pBNS->vert[i].st_edge.cap - delta == (pBNS->vert[i].st_edge.flow & EDGE_FLOW_ST_MASK) ) {
            if ( nNumRad < MAX_NUM_RAD ) {
                pBNS->vert[i].st_edge.cap -= delta;
                pBNS->tot_st_cap          -= delta;
                vRadList[nNumRad] = i;       /* radical position; i > j <=> vRadList[i] > vRadList[j]  */
                vRadEqul[nNumRad] = nNumRad; /* the smallest radical atom that has reachable
                                              * atoms in common with this radical atom
                                              * always keep vRadEqul[nNumRad] <= nNumRad */
                nNumRad ++;
            }
        }
    }
    if ( pBNS->tot_st_cap - pBNS->tot_st_flow > nNumRad ) {
        return BNS_CAP_FLOW_ERR; /* extra st_cap on non-atoms or program error */
    }
    memset( &VertSet, 0, sizeof(VertSet) );
    /* find reachable atoms by enabling each radical separately */
    for ( j = 0; j < nNumRad; j ++ ) {
        i  = vRadList[j];
        pBD->nNumRadEndpoints = 0;
        pBD->nNumRadEdges     = 0;
        pBD->bRadSrchMode     = bRadSrchMode;
        pBNS->alt_path = pBNS->altp[0];
        pBNS->bChangeFlow     = 0;
        pBNS->vert[i].st_edge.cap += delta; /* enable single radical */
        pBNS->tot_st_cap          += delta;
        ret = BalancedNetworkSearch( pBNS, pBD, BNS_EF_RAD_SRCH ); /* find reachable atoms */
        ReInitBnData( pBD );
        ReInitBnStructAltPaths( pBNS );
        pBD->bRadSrchMode     = RAD_SRCH_NORM;
        pBNS->vert[i].st_edge.cap -= delta; /* disable single radical */
        pBNS->tot_st_cap          -= delta;
        if ( IS_BNS_ERROR( ret ) ) {
            goto error_exit;
        } else
        if ( ret ) {
            ret = BNS_RADICAL_ERR; /* found augmenting path: should not happen since only one radical was enabled */
            goto error_exit;
        }
        if ( !ret && pBD->nNumRadEndpoints >= 2 ) {
            /* sort by: primary_key=radical locations, secondary_key=radical endoint */
            qsort( pBD->RadEndpoints, pBD->nNumRadEndpoints/2, 2*sizeof(pBD->RadEndpoints[0]), cmp_rad_endpoints );
            if ( pBD->RadEndpoints[0] != i || pBD->RadEndpoints[pBD->nNumRadEndpoints-2] != i ) {
                ret = BNS_RADICAL_ERR; /* more than one radical vertex */
                goto error_exit;
            }
            if ( nNumRad > 1 ) {
                /* if more than one radical then save reachable atoms in bitmaps to allow */
                /* faster finding whether same atoms are reachable by two or more radicals */
                /* Later merge such sets */
                if ( NULL == VertSet.bitword ) {
                    SetBitCreate( );
                    if ( !NodeSetCreate( &VertSet, pBNS->num_atoms, nNumRad ) ) {
                        ret = BNS_OUT_OF_RAM; /* out of RAM */
                        goto error_exit;
                    }
                }
                NodeSetFromRadEndpoints( &VertSet, j, pBD->RadEndpoints, pBD->nNumRadEndpoints);
                /* do not allow any radical center be treated as a reachable atom: */
                RemoveFromNodeSet( &VertSet, j, vRadList, nNumRad );
            }
        }
    }
    /* restore radical st_cap so that st_cap-st_flow=delta */
    for ( j = 0; j < nNumRad; j ++ ) {
        i  = vRadList[j];
        pBNS->vert[i].st_edge.cap += delta;
        pBNS->tot_st_cap          += delta;
    }
    /* merge lists that have common radical endpoints */
    /* defect: if vertex sets i and j do not intersect they will be compared 2 times */
    /* total up to nNumRad*(nNumRad-1)/2 calls to DoNodeSetsIntersect() */
    if ( nNumRad > 1 ) {
        for ( i = 0; i < nNumRad; i ++ ) {
            if ( vRadEqul[i] != i )
                continue;
            do {
                n = 0;
                for ( j = i+1; j < nNumRad; j ++ ) {
                    if ( vRadEqul[j] != j )
                        continue;
                    if ( DoNodeSetsIntersect( &VertSet, i, j) ) {
                        AddNodeSet2ToNodeSet1( &VertSet, i, j);
                        vRadEqul[j] = i; /* Set j was copied to set i; i < j */
                        n ++;
                    }
                }
            } while( n );
        }
        /* fill out pBD->RadEndpoints[] */
        for ( i = 0, n = 0; i < nNumRad; i ++ ) {
            if ( i == vRadEqul[i] ) {
                if ( !IsNodeSetEmpty( &VertSet, i) ) {
                    /* store equivalent radicals */
                    for ( j = i+1; j < nNumRad; j ++ ) {
                        if (i == vRadEqul[j] ) {
                            pBD->RadEndpoints[n++] =  vRadList[i];
                            pBD->RadEndpoints[n++] = -vRadList[j]-2; /* equivalent radical, alvays not zero */
                        }
                    }
                    /* store endpoints */
                    n = AddNodesToRadEndpoints( &VertSet, i, pBD->RadEndpoints, vRadList[i], n, pBD->max_len_Pu_Pv );
                    if ( n < 0 ) {
                        ret = BNS_RADICAL_ERR; /* pBD->RadEndpoints overflow */
                        goto error_exit;
                    }
                } else {
                    pBD->RadEndpoints[n++] =  vRadList[i];
                    pBD->RadEndpoints[n++] =  -1; /* immobile radical, only one edge to add */
                }
            }
        }
        pBD->nNumRadEndpoints = n;
        NodeSetFree( &VertSet );
    } else
    if ( nNumRad == 1 && !pBD->nNumRadEndpoints ) {
        /* 2006-07-30: a single radical; no possible endpoint found */
        for ( i = 0, n = 0; i < nNumRad; i ++ ) {
            pBD->RadEndpoints[n++] =  vRadList[i];
            pBD->RadEndpoints[n++] =  -1; /* immobile radical, only one edge to add */
        }
        pBD->nNumRadEndpoints = n;
    }

    if ( !ret && pBD->nNumRadEndpoints >= 2 ) {
        /* already sorted by radical locations */
        num_new_edges = 0;
        nNumRadicals  = 0;
        /**************************************************************************
         * create new vertices (type=BNS_VERT_TYPE_TEMP) and edges with flow=cap=1
         * connecting the new vertices radical vertices
         *
         *  
         * Original structure:    atom A is a radical center    A==B--C*--D==E
         *   A*--B==C--D==E       atoms C and E are reachable:  A==B--C===D--E*
         *
         * Resultant temporary structure:
         *   A---B==C--D==E                     
         *  ||     /     /                      
         *  ||    /    /          The additional new vertex (*) and its
         *  ||   /   /            3 edges replace the radical with alternating
         *  ||  /  /              circuits that allow same bond changes
         *  || / /                as moving the radical to atoms C or E.
         *  ||//                  "Double bonds" here have edge cap=1, flow=1
         *  (*)                   "Single bonds" have edge cap=1, flow=0
         *                        
         *   The "equivalent radical centers" (which have at least one reachable atom
         *   in common) are connected to (*) with "double bonds" (edge cap=1, flow=1).
         *   Reachable non-radical atoms are connected by edges with cap=1, flow=0
         *   After running BNS to find alt.path a "double bond" from (*) may move
         *   to another atom thus muving the radical.
         *
         *   Number of additional (*) vertices = number of sets of
         *   "equivalent radical centers".
         *   Each such a set may include one or more radical centers.
         *
         *   The radicals will be re-created in RemoveRadEndpoints()
         ***************************************************************************/
        for ( i = 0; i < pBD->nNumRadEndpoints; i = j ) {
            wRad = pBD->RadEndpoints[i];
            pRad = pBNS->vert + wRad;
            delta = pRad->st_edge.cap - (pRad->st_edge.flow & EDGE_FLOW_ST_MASK);
            if ( delta <= 0 ) {
                delta = 1;
            }
            nNumEdges = 0;
            for ( j = i; j < pBD->nNumRadEndpoints && wRad == pBD->RadEndpoints[j] ; j += 2 ) {
                nNumEdges += (pBD->RadEndpoints[j+1] != -1); /* immobile radicals have one edge only */
            }
            /* add new aux vertex to the radical atom/vertex making st_cap-st_flow=0 */
            /* in case of immobile radical there will be no additional eddges since nNumEdges=0 */
            vRad = bAddNewVertex( pBNS, wRad, delta, delta, nNumEdges+1, &nDots );
            if ( IS_BNS_ERROR( vRad ) ) {
                ret = vRad;
                goto error_exit;
            }
            pRad     = pBNS->vert + vRad;
            pBD->RadEdges[pBD->nNumRadEdges ++] = pRad->iedge[pRad->num_adj_edges-1];
            /* replace references to vertex wRad with vRad */
            for ( k = i, nNumEdges = 0; k < j; k += 2 ) {
                pBD->RadEndpoints[k] = vRad;
            }
            nNumRadicals ++;
        }
        /* all vRad vertex indices should be in the range vFirstNewVertex...vFirstNewVertex+nNumRadicals-1 */
        /* connect new vertices to the radical endpoints thus replacing radicals with even-length alternating cycles */
        for ( i = 0; i < pBD->nNumRadEndpoints; i = j ) {
            vRad = pBD->RadEndpoints[i];
            pRad = pBNS->vert + vRad;
            for ( j = i; j < pBD->nNumRadEndpoints && vRad == pBD->RadEndpoints[j] ; j += 2 ) {
                /* connect vew vertex pRad to radical endpoints */
                vEndp     = pBD->RadEndpoints[j+1];
                if ( vEndp == -1 )
                    continue;
                if ( vEndp < 0 ) {
                    edge_flow = 1;
                    vEndp = -vEndp - 2; /* equivalent radical centers */
                } else {
                    edge_flow = 0;
                }
                pEndp    = pBNS->vert + vEndp;
                ret = AddNewEdge( pRad, pEndp, pBNS, 1, edge_flow );
                if ( IS_BNS_ERROR( ret ) ) {
                    goto error_exit;
                }
                pBD->RadEdges[pBD->nNumRadEdges ++] = ret;
            }
        }
        pBD->nNumRadicals = nNumRadicals;
        return nNumRadicals; /* done */
    }
    return 0; /* nothing to do */

error_exit:
    RemoveRadEndpoints( pBNS, pBD, NULL );
    NodeSetFree( &VertSet );
    return ret;

}


#else
/**********************************************************************************/
int SetRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode )
{
    return 0;
}
int RemoveRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at )
{
    return 0;
}
int SetRadEndpoints2( BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode )
{
    return 0;
}
#endif
/**********************************************************************************/
/* Return value ret bits if not IS_BNS_ERROR(ret):

    ret & 1         => Success
    ret & 2         => Bonds changed to Alt
    (ret & ~3) >> 2 => nDelta: number of removed dots
*/
int bExistsAltPath( BN_STRUCT *pBNS, BN_DATA *pBD, BN_AATG *pAATG, inp_ATOM *at, int num_atoms,
                    int nVertDoubleBond, int nVertSingleBond, int path_type )
{
    ALT_PATH_CHANGES apc;
    int ret, ret_val, bError, bSuccess, bChangeFlow=0, nDots, nDelta, bDoMarkChangedBonds = 1;
    int bAdjustRadicals = 0;
    AT_NUMB          type;
    BNS_FLOW_CHANGES fcd[4*BNS_MAX_NUM_FLOW_CHANGES+1];
    ENDPOINT_INFO    eif;
#if( KETO_ENOL_TAUT == 1 )
    ENDPOINT_INFO    eif2;
#endif

    /* initialize */
    switch( path_type ) {
    case ALT_PATH_MODE_TAUTOM:
        /* Check for alt path allowing to move H and (-). Purpose: confirm possible tautomerism */
        type        = BNS_VERT_TYPE_ENDPOINT;
        bChangeFlow = BNS_EF_CHNG_RSTR;
        if ( !at[nVertSingleBond].endpoint &&
             (!nGetEndpointInfo( at, nVertSingleBond, &eif ) || !eif.cDonor ) )
            return 0;
        if ( !at[nVertDoubleBond].endpoint &&
             (!nGetEndpointInfo( at, nVertDoubleBond, &eif ) || !eif.cAcceptor ) )
            return 0;
        break;

#if( KETO_ENOL_TAUT == 1 )
    case ALT_PATH_MODE_TAUTOM_KET:
        /* Check for alt path allowing to move H and (-). Purpose: confirm possible tautomerism */
        type        = BNS_VERT_TYPE_ENDPOINT;
        bChangeFlow = BNS_EF_CHNG_RSTR;

        if ( !at[nVertSingleBond].endpoint &&
             (!nGetEndpointInfo_KET( at, nVertSingleBond, &eif ) || !eif.cDonor ) )
            return 0;
        if ( !at[nVertDoubleBond].endpoint &&
             (!nGetEndpointInfo_KET( at, nVertDoubleBond, &eif2 ) || !eif2.cAcceptor ) )
            return 0;
        /*
        if ( eif.cKetoEnolCode + eif2.cKetoEnolCode != 3 )
            return 0;
        */
        break;

#endif
    case ALT_PATH_MODE_CHARGE:
        /* Find alt path allowing to move (+). Purpose: establish "charge groups",
           mark alt. bonds due to (+) charge movement */
        type        = BNS_VERT_TYPE_C_POINT;
        bChangeFlow = (BNS_EF_CHNG_RSTR | BNS_EF_ALTR_BONDS);
        break;

    case ALT_PATH_MODE_4_SALT:
    case ALT_PATH_MODE_4_SALT2:
        /* Find alt paths allowing to move (-) and H between "acidic oxygen atoms".
           Purpose: mark alt bonds due to this "long range" tautomerism. */
        type        = BNS_VERT_TYPE_ENDPOINT;
        bChangeFlow = (BNS_EF_CHNG_RSTR | BNS_EF_ALTR_BONDS);
        if ( !bIsBnsEndpoint( pBNS, nVertSingleBond ) /* !at[nVertSingleBond].endpoint*/ &&
             (!nGetEndpointInfo( at, nVertSingleBond, &eif ) || !eif.cDonor ) )
            return 0;
        if ( !bIsBnsEndpoint( pBNS, nVertDoubleBond ) /* !at[nVertDoubleBond].endpoint*/ &&
             (!nGetEndpointInfo( at, nVertDoubleBond, &eif ) || !eif.cAcceptor ) )
            return 0;
        memset( &apc, 0, sizeof(apc) );
        break;

    case ALT_PATH_MODE_REM2H_CHG:
        bChangeFlow |= BNS_EF_ALTR_BONDS; /* fall through */
    case ALT_PATH_MODE_REM2H_TST:
        bChangeFlow |= BNS_EF_CHNG_RSTR;
        type        = BNS_VERT_TYPE_ENDPOINT;
        /* allow non-tautomeric donors or any tautomeric atom */
        if ( !bIsBnsEndpoint( pBNS, nVertSingleBond ) /* not linked to a t-group or the edge forbidden */&&
             (!nGetEndpointInfo( at, nVertSingleBond, &eif ) || !eif.cDonor ) ) /* not a donor */
            return 0;
        if ( !bIsBnsEndpoint( pBNS, nVertDoubleBond ) /* not connected to a t-group */ &&
             (!nGetEndpointInfo( at, nVertDoubleBond, &eif ) || !eif.cDonor ) )
            return 0;
        memset( &apc, 0, sizeof(apc) );
        break;

    case ALT_PATH_MODE_ADD2H_CHG:
        bChangeFlow |= BNS_EF_ALTR_BONDS; /* fall through */
    case ALT_PATH_MODE_ADD2H_TST:
        bChangeFlow |= BNS_EF_CHNG_RSTR;
        type        = BNS_VERT_TYPE_ENDPOINT;
        /* allow non-tautomeric acceptors or any tautomeric atom */
        if ( !bIsBnsEndpoint( pBNS, nVertSingleBond ) /* !at[nVertSingleBond].endpoint*/ &&
             (!nGetEndpointInfo( at, nVertSingleBond, &eif ) || !eif.cAcceptor ) )
            return 0;
        if ( !bIsBnsEndpoint( pBNS, nVertDoubleBond ) /* !at[nVertSingleBond].endpoint*/ &&
             (!nGetEndpointInfo( at, nVertDoubleBond, &eif ) || !eif.cAcceptor ) )
            return 0;
        break;

    case ALT_PATH_MODE_REM_PROTON:
        /* alt path is between the t-group (nVertDoubleBond) and
           the (+)-charge group (nVertSingleBond) */
        type                = 0;
        /*bDoMarkChangedBonds = 0;*/
        bChangeFlow         = (BNS_EF_SAVE_ALL | BNS_EF_UPD_H_CHARGE) | BNS_EF_ALTR_NS; /* added BNS_EF_ALTR_NS: set non-stereo altern non-ring bonds 2004-07-02*/
        break;
    default:
        type = 0;
        bChangeFlow = BNS_EF_CHNG_RSTR;
        break;
    }
    
    bError      = 0;
    bSuccess    = 0;
    nDelta      = 0;

    ret = SetRadEndpoints2( pBNS, pBD, RAD_SRCH_NORM );
    if ( IS_BNS_ERROR( ret ) ) {
        return ret;
    }

    /* set BNS to check alt path */
    ret = bSetBnsToCheckAltPath( pBNS, nVertDoubleBond, nVertSingleBond, type, path_type, &apc, fcd, &nDots );
    switch( ret ) {
    case  BNS_CHK_ALTP_NO_ALTPATH:
        ret = RemoveRadEndpoints( pBNS, pBD, NULL );
        return ret;
    case BNS_CHK_ALTP_SAME_TGROUP:
        bSuccess = 1;
        goto reinit_BNS;
    case BNS_CHK_ALTP_SAME_VERTEX:
        ret = RemoveRadEndpoints( pBNS, pBD, NULL );
        return ret? ret : 1; /* very strange ... set a breakpoint here */
    case BNS_CHK_ALTP_SET_SUCCESS:
        break;  /* actually check the existence of the altpath */
    case BNS_CANT_SET_BOND:
        goto reinit_BNS;
    default:
        ret_val = RemoveRadEndpoints( pBNS, pBD, NULL );
        if ( IS_BNS_ERROR( ret ) ) {
            return ret;
        }
        return BNS_PROGRAM_ERR;
    }

    bAdjustRadicals = ( (bChangeFlow & BNS_EF_UPD_RAD_ORI) && !(bChangeFlow & BNS_EF_RSTR_FLOW) );

    /*****************************************************************
     * nDots = 2 for ALT_PATH_CHARGE (checking moveable positive charges)
     * Now nDots for ALT_PATH_TAUTOM or ALT_PATH_4_SALT can be greater
     * because some of the bonds are effectively removed and dots
     * (vertex st-caps) may be added
     * -- to make sure there is no (+) charge on a tautomeric endpoint
     * -- to fix positions of moveable tautomeric attachements
     *    (H and (-)-charges) at the ends of an alt path
     */

    /* run BNS */

    ret = RunBalancedNetworkSearch( pBNS, pBD, bChangeFlow );
    if ( IS_BNS_ERROR( ret ) ) {
        bError = ret;
    } else
    if ( ret > 0 ) {
        if ( 2*ret >= nDots ) {
            nDelta = 2*ret - nDots;  /* non-zero means augmentation created another alt. path -- between radicals */
            if ( pAATG && pAATG->nMarkedAtom ) {
                if ( pAATG->nAtTypeTotals && (bChangeFlow & BNS_EF_UPD_H_CHARGE) ) {
                    memset( pAATG->nMarkedAtom, 0, num_atoms*sizeof(pAATG->nMarkedAtom[0]) );
                    /* mark atoms that have charge or H changed, check their input types (that is, before changes),
                       and subtract their input charge/H from nAtTypeTotals */
                    SubtractOrChangeAtHChargeBNS( pBNS, at, num_atoms, pAATG->nAtTypeTotals, pAATG->nMarkedAtom, NULL, 1 );
                    /* ZChange charges and/or H, update t_group_info, do not check types or change nAtTypeTotals */
                    /* Atom types will be checked and nAtTypeTotals will be changed in
                       AddChangedAtHChargeBNS() later */
                    SubtractOrChangeAtHChargeBNS( pBNS, at, num_atoms, NULL, NULL, pAATG->t_group_info, 0 );
                } else
                if ( !pAATG->nAtTypeTotals ){
                    bDoMarkChangedBonds = MarkAtomsAtTautGroups(  pBNS, num_atoms, pAATG, nVertSingleBond, nVertDoubleBond );
                    if ( bDoMarkChangedBonds < 0 ) {
                        bError = bDoMarkChangedBonds;
                        bDoMarkChangedBonds = 0;
                    }
                }
            }
            if ( bDoMarkChangedBonds ) {
                /* mark bonds that were changed to configure bond testing */
                ret_val = bSetBondsAfterCheckOneBond( pBNS, fcd, -1, at, num_atoms, bChangeFlow );
                if ( IS_BNS_ERROR( ret_val ) ) {
                    bError = ret_val;
                }
                /*ret = SetBondsRestoreBnStructFlow( pBNS, at, num_atoms, bChangeFlow );*/
                /* mark all other changed bonds */
                ret = SetBondsFromBnStructFlow( pBNS, at, num_atoms, bChangeFlow );
                if ( IS_BNS_ERROR( ret ) ) {
                    bError = ret;
                } else
                if ( !(ret & 1) && !(ret_val & 1) ) {
                    bSuccess = 1;
                } else
                if ( ((ret & 1) || (ret_val & 1)) &&
                     (bChangeFlow & BNS_EF_ALTR_BONDS) || (bChangeFlow & BNS_EF_UPD_H_CHARGE) ) {
                    /* some bonds have been changed to alternating */
                    bSuccess = 3;
                } else {
                    bError = BNS_BOND_ERR;
                }
                if ( !bError && pAATG && pAATG->nMarkedAtom && (bChangeFlow & BNS_EF_UPD_H_CHARGE) ) {
                    /* Update radicals to avoid errors in atom type check in AddChangedAtHChargeBNS() */
                    if ( bAdjustRadicals ) {
                        ret_val = RestoreRadicalsOnly( pBNS, pBD, at );
                        if ( IS_BNS_ERROR( ret_val ) ) {
                            bError = ret_val;
                        }
                    }
                    /* Check atom types of marked atoms and add charge/H changes to nAtTypeTotals */
                    /* Changing atoms were marked in the 1st call to SubtractOrChangeAtHChargeBNS(..., 1) above */
                    AddChangedAtHChargeBNS( at, num_atoms, pAATG->nAtTypeTotals, pAATG->nMarkedAtom );
                    if ( bChangeFlow & BNS_EF_CHNG_FLOW ) {
                        /* eliminate ambiguities in already changed flow:
                           replace (+)--N==(-) with (+)==N--(-) (both represent neutral N) */
                        EliminatePlusMinusChargeAmbiguity( pBNS, num_atoms );
                    }
                }
            }
        }
        ret = RestoreBnStructFlow( pBNS, bChangeFlow & BNS_EF_CHNG_RSTR);
        if ( IS_BNS_ERROR( ret ) ) {
            bError = ret;
        }
    }
reinit_BNS:
    /* --- reinitialize to repeat the calculations --- */
    bRestoreBnsAfterCheckAltPath( pBNS, &apc, bChangeFlow & BNS_EF_UPD_H_CHARGE );
    bRestoreFlowAfterCheckOneBond( pBNS, fcd );
    ret_val = RemoveRadEndpoints( pBNS, pBD, bAdjustRadicals? at : NULL );
    ReInitBnStructAltPaths( pBNS );
    return bError? bError : ret_val? ret_val : (bSuccess + 4*nDelta);
}


/*********************************************************************************/
BN_STRUCT* AllocateAndInitBnStruct( inp_ATOM *at, int num_atoms, int nMaxAddAtoms, int nMaxAddEdges, int max_altp, int *pNum_changed_bonds )
{
    BN_STRUCT   *pBNS         = NULL;
    BNS_VERTEX  *vert;

    int    neigh, num_changed_bonds=0;
    U_CHAR bond_type, bond_mark;

    int i, j, k, n_edges, num_bonds, num_edges, f1, f2, edge_cap, edge_flow, st_cap, st_flow, flag_alt_bond;
    int tot_st_cap, tot_st_flow;
    int max_tg, max_edges, max_vertices, len_alt_path, max_iedges, num_altp;
#if( BNS_RAD_SEARCH == 1 )
    int num_rad = 0;

    nMaxAddEdges += 1;
#endif
#if( FIX_NUM_TG == 1 )
    max_tg = inchi_max( num_atoms / 2, 5);
#else
    max_tg = num_atoms;
#endif
    num_changed_bonds = 0;
    
    for ( i = 0, num_bonds = 0; i < num_atoms; i ++ ) {
        num_bonds += at[i].valence;
#if( BNS_RAD_SEARCH == 1 )
        num_rad   += (at[i].radical == RADICAL_DOUBLET);
#endif
    }
    /* each atom has enough edges to belong to a tautomeric group + nMaxAddEdges */
    /* number of atoms is large enough to accommodate max. possible number of t-groups + nMaxAddAtoms */
    /* max_altp cannot be larger than BN_MAX_ALTP = 16 */
    num_edges    = (num_bonds /= 2);
    /* +1 for a super-tautomeric group */
    max_vertices = num_atoms + nMaxAddAtoms + max_tg + 1;
    /* +max_tg for edges between t-groups and super-tautomeric group */
    max_edges    = num_edges + (nMaxAddEdges + NUM_KINDS_OF_GROUPS)*max_vertices + max_tg;
#if( BNS_RAD_SEARCH == 1 )
    if ( num_rad ) {
        max_vertices *= 2;
        max_edges    *= 2;
    }
#endif
    max_iedges   = 2*max_edges;
    len_alt_path = max_vertices+iALTP_HDR_LEN+1; /* may overflow if an edge is traversed in 2 directions */

    if ( !( pBNS           = (BN_STRUCT   *)inchi_calloc( 1,           sizeof(BN_STRUCT)) )  ||
         !( pBNS->edge     = (BNS_EDGE    *)inchi_calloc( max_edges,   sizeof(BNS_EDGE)) )   ||
         !( pBNS->vert     = (BNS_VERTEX  *)inchi_calloc( max_vertices,sizeof(BNS_VERTEX)) ) ||
         !( pBNS->iedge    = (BNS_IEDGE   *)inchi_calloc( max_iedges,  sizeof(BNS_IEDGE)) ) ) { 
        return DeAllocateBnStruct( pBNS );
    }
    /* alt path init */
    for ( num_altp = 0; num_altp < max_altp && num_altp < BN_MAX_ALTP; num_altp ++ ) {
        if ( !( pBNS->altp[num_altp] = (BNS_ALT_PATH*)inchi_calloc( len_alt_path,sizeof(BNS_ALT_PATH))) ) {
            return DeAllocateBnStruct( pBNS );
        }
        ALTP_ALLOCATED_LEN(pBNS->altp[num_altp]) = len_alt_path;
        pBNS->len_alt_path                 = len_alt_path;  /* ??? duplication ??? */
        /* re-init */
        ALTP_DELTA(pBNS->altp[num_altp])         = 0;
        ALTP_START_ATOM(pBNS->altp[num_altp])    = NO_VERTEX;
        ALTP_END_ATOM(pBNS->altp[num_altp])      = NO_VERTEX;
        ALTP_PATH_LEN(pBNS->altp[num_altp])      = 0;
    }
    pBNS->alt_path = NULL;
    pBNS->num_altp = 0;
    pBNS->max_altp = num_altp;


    /* fill vertices (no connectivity) */
    pBNS->vert[0].iedge = pBNS->iedge;
    for ( i = 0; i < num_atoms; i ++ ) {
        k = pBNS->vert[i].max_adj_edges = at[i].valence + (nMaxAddEdges + NUM_KINDS_OF_GROUPS);
        pBNS->vert[i+1].iedge = pBNS->vert[i].iedge + k;
    }
    pBNS->num_atoms       = num_atoms;      /* number of real atoms */
    pBNS->num_added_atoms = 0;
    pBNS->num_t_groups    = 0;              /* number of added t-groups */
    pBNS->num_c_groups    = 0;
    pBNS->nMaxAddAtoms    = nMaxAddAtoms;
    pBNS->nMaxAddEdges    = nMaxAddEdges;

    pBNS->num_vertices    = num_atoms;      /* current number of vertices, a sum of
                                               pBNS->num_atoms
                                               pBNS->num_t_groups
                                               pBNS->num_added_atoms
                                            */
    pBNS->max_vertices    = max_vertices;
    

    pBNS->num_bonds       = num_bonds;      /* number of real edges (bonds) */
    pBNS->max_edges       = max_edges;
    pBNS->max_iedges      = max_iedges;

    /* 
       To remove t-groups and added atoms:
       In atoms i = 0..pBNS->num_atoms-1
            pBNS->vert[i].num_adj_edges = pBNS->vert[i].max_adj_edges - pBNS->nMaxAddEdges - NUM_KINDS_OF_GROUPS;
       pBNS->num_vertices    = pBNS->num_atoms;
       pBNS->num_edges       = pBNS->num_bonds;
       pBNS->num_added_atoms = 0;
       pBNS->num_t_groups    = 0;
       pBNS->num_added_edges = 0;

       ALTP_DELTA(pBNS->alt_path)      = 0;
       ALTP_START_ATOM(pBNS->alt_path) = NO_VERTEX;
       ALTP_END_ATOM(pBNS->alt_path)   = NO_VERTEX;
       ALTP_PATH_LEN(pBNS->alt_path)   = 0;

    */


    /* fill edges and connectivity */
    tot_st_cap = tot_st_flow = 0;
    for ( i = 0, n_edges = 0; i < num_atoms; i ++ ) {
        vert    = &pBNS->vert[i];
        st_cap  = 0;
        st_flow = 0;
        flag_alt_bond = 0;
        for ( j = 0; j < at[i].valence; j ++ ) {
            neigh = at[i].neighbor[j];
            /* find this bond at the neighbor */
            for ( k = 0; k < at[neigh].valence; k ++ ) {
                if ( at[neigh].neighbor[k] == i ) {
                    break;
                }
            }
            bond_type = (at[i].bond_type[j] & BOND_TYPE_MASK);
            bond_mark = (at[i].bond_type[j] & ~BOND_TYPE_MASK);
            if ( bond_type != BOND_SINGLE && bond_type != BOND_DOUBLE &&
                 bond_type != BOND_TRIPLE /*&& bond_type != BOND_ALTERN*/ ) {
                /* make Unknown or Alternating bonds single */
                bond_type = 1;
                at[i].bond_type[j] = bond_mark | bond_type;
                num_changed_bonds ++;
            }
            if ( neigh > i ) {
                /* this is the first time we encounter this bond */
                f1 = MAX_AT_FLOW(at[i]);
                f2 = MAX_AT_FLOW(at[neigh]);
                edge_flow = bond_type-1;
                if ( edge_flow > MAX_BOND_EDGE_CAP ) {
                    flag_alt_bond ++;
                    edge_flow = 0;  /* BNS will determine flows (that is, bonds) */
                    edge_cap  = AROM_BOND_EDGE_CAP;
                } else {
#if ( 0 && KETO_ENOL_TAUT == 1 )  /* ????? */
                    edge_cap  = inchi_max(f1, f2);
#else
                    edge_cap  = inchi_min(f1, f2);
#endif
                    edge_cap  = inchi_min(edge_cap, MAX_BOND_EDGE_CAP); /* max capacity = 2 means up to triple bond */
                }
                pBNS->edge[n_edges].neighbor1    = (AT_NUMB)i;
                pBNS->edge[n_edges].neighbor12   = (AT_NUMB)(i ^ neigh);
                pBNS->edge[n_edges].flow =
                pBNS->edge[n_edges].flow0        = edge_flow;
                pBNS->edge[n_edges].cap  =
                pBNS->edge[n_edges].cap0         = edge_cap;
                pBNS->edge[n_edges].neigh_ord[0] = j;
                pBNS->edge[n_edges].neigh_ord[1] = k;
                pBNS->edge[n_edges].pass         = 0;
                pBNS->edge[n_edges].forbidden    = 0;

                vert->iedge[j] = pBNS->vert[neigh].iedge[k] = n_edges ++;
            } else {
                /* this is the second time we encounter this bond. It was stored at */
                int  iedge = pBNS->vert[neigh].iedge[k];
                edge_cap   = pBNS->edge[iedge].cap;
                edge_flow  = pBNS->edge[iedge].flow;
            }
            st_flow += edge_flow;
            st_cap  += edge_cap;
        }
        vert->num_adj_edges = j;
        vert->st_edge.cap   =
        vert->st_edge.cap0  = MAX_AT_FLOW(at[i]);
        vert->st_edge.flow  =
        vert->st_edge.flow0 = st_flow;
        vert->type          = BNS_VERT_TYPE_ATOM;
        tot_st_cap  += vert->st_edge.cap;
        tot_st_flow += vert->st_edge.flow;
    }
    *pNum_changed_bonds = num_changed_bonds/2;

    pBNS->num_edges       =  n_edges;   /* number of edges */
    pBNS->num_added_edges = 0;

    pBNS->tot_st_cap  = tot_st_cap;
    pBNS->tot_st_flow = tot_st_flow;

    return pBNS;
}
/*********************************************************************************/
BN_STRUCT* DeAllocateBnStruct( BN_STRUCT *pBNS )
{
    int i;
    if ( pBNS ) {
        if ( pBNS->edge ) {
            inchi_free( pBNS->edge );
        }
        for ( i = 0; i < pBNS->max_altp && i < BN_MAX_ALTP ; i ++ ) {
            if ( pBNS->altp[i] ) {
                inchi_free( pBNS->altp[i] );
            }
        }
        if ( pBNS->vert ) {
            if ( pBNS->vert[0].iedge ) {
                inchi_free( pBNS->vert[0].iedge );
            }
            inchi_free( pBNS->vert );
        }
        inchi_free( pBNS );
    }
    return NULL;
}
/*********************************************************************************/
int ReInitBnStructAltPaths( BN_STRUCT *pBNS )
{
    int i;
    for ( i = 0; i < pBNS->max_altp && i < BN_MAX_ALTP; i ++ ) {
        if ( pBNS->altp[i] ) {
            ALTP_DELTA(pBNS->altp[i])      = 0;
            ALTP_PATH_LEN(pBNS->altp[i])   = 0;
            ALTP_START_ATOM(pBNS->altp[i]) = NO_VERTEX;
            ALTP_END_ATOM(pBNS->altp[i])   = NO_VERTEX;
        }
    }
    pBNS->alt_path = NULL;
    pBNS->num_altp = 0;
    return i;
}
/*********************************************************************************/
int ReInitBnStructAddGroups( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, T_GROUP_INFO *tgi, C_GROUP_INFO *cgi )
{
    int ret;
    /* strip all t-groups and c-groups */
    ret = ReInitBnStruct( pBNS, at, num_atoms, 0 );
    if ( ret ) {
        ret = BNS_REINIT_ERR;
        goto exit_function;
    }
/*#if( MOVE_CHARGES == 1 )*/
    if ( *pBNS->pbTautFlags & TG_FLAG_MOVE_POS_CHARGES ) {
    /* add c-groups */
    ret = AddCGroups2BnStruct( pBNS, at, num_atoms, cgi );
    if ( IS_BNS_ERROR( ret ) ) {
        goto exit_function;
    }
    }
/*#endif*/
    /* add t-groups */
    ret = AddTGroups2BnStruct( pBNS, at, num_atoms, tgi );
    if ( IS_BNS_ERROR( ret ) ) {
        goto exit_function;
    }
exit_function:
    return ret;
}
/*********************************************************************************/
int ReInitBnStruct( BN_STRUCT *pBNS, inp_ATOM *at, int num_at, int bRemoveGroupsFromAtoms )
{
    int i, vfict, kfict, iedgefict, endpoint, centerpoint, iedge, k;
    int ret = 0;
    if ( pBNS ) {
        if ( pBNS->vert && pBNS->edge ) {
            /* debug */
            for ( k = 0, i = 0; k < pBNS->num_edges; k ++ ) {
                if ( pBNS->edge[k].pass ) {
                    i ++;
                }
            }
            ret += i * 100;
            /* restore flow and cap on edges to vertices connected to fictitious atoms */
            for ( vfict = pBNS->num_atoms; vfict < pBNS->num_vertices; vfict ++ ) {
                for ( kfict = 0; kfict < pBNS->vert[vfict].num_adj_edges; kfict ++ ) {
                    iedgefict = pBNS->vert[vfict].iedge[kfict]; /* fictitious edge to the endpoint */
                    endpoint  = pBNS->edge[iedgefict].neighbor12 ^ vfict;  /* the endpoint */
                    /* to simlify restore cap and flow in ALL edges to the endpoint */
                    if ( bRemoveGroupsFromAtoms && endpoint < num_at ) {
                        at[endpoint].c_point  = 0;
                        at[endpoint].endpoint = 0;
                    }
                    for ( k = 0; k < pBNS->vert[endpoint].num_adj_edges; k ++ ) {
                        iedge = pBNS->vert[endpoint].iedge[k]; /* edge to endpoint */
                        centerpoint = pBNS->edge[iedge].neighbor12 ^ endpoint;
                        pBNS->edge[iedge].cap       = pBNS->edge[iedge].cap0;
                        pBNS->edge[iedge].flow      = pBNS->edge[iedge].flow0;
                        pBNS->edge[iedge].pass      = 0;
#if( RESET_EDGE_FORBIDDEN_MASK == 1 )
                        pBNS->edge[iedge].forbidden &= pBNS->edge_forbidden_mask;
#endif
                        pBNS->vert[centerpoint].st_edge.cap  = pBNS->vert[centerpoint].st_edge.cap0;
                        pBNS->vert[centerpoint].st_edge.flow = pBNS->vert[centerpoint].st_edge.flow0;
                    }
                    pBNS->vert[endpoint].st_edge.cap  = pBNS->vert[endpoint].st_edge.cap0;
                    pBNS->vert[endpoint].st_edge.flow = pBNS->vert[endpoint].st_edge.flow0;
                    pBNS->vert[endpoint].type &= BNS_VERT_TYPE_ATOM;
                }
            }
            /* reset number of neighbors */
            if ( pBNS->num_edges > pBNS->num_bonds ) {
                for ( i = 0; i < pBNS->num_atoms; i ++ ) {
                    pBNS->vert[i].num_adj_edges =
                    pBNS->vert[i].max_adj_edges - pBNS->nMaxAddEdges - NUM_KINDS_OF_GROUPS;
                }
            }
        } else {
            ret += 2;
        }
        if ( !pBNS->edge ) {
            ret += 4;
        }
        if ( !pBNS->iedge ) {
            ret += 8;
        }

        ReInitBnStructAltPaths( pBNS );

        pBNS->num_vertices    = pBNS->num_atoms;
        pBNS->num_edges       = pBNS->num_bonds;
        pBNS->num_added_atoms = 0;
        pBNS->num_t_groups    = 0;
        pBNS->num_c_groups    = 0;
        pBNS->num_added_edges = 0;

    } else {
        ret += 1;
    }

    return ret;
}
/*********************************************************************************/
int CompTGroupNumber( const void *tg1, const void *tg2 )
{
    return (int)((const T_GROUP *)tg1)->nGroupNumber - (int)((const T_GROUP *)tg2)->nGroupNumber;
}
/*********************************************************************************/
int CompCGroupNumber( const void *cg1, const void *cg2 )
{
    return (int)((const C_GROUP *)cg1)->nGroupNumber - (int)((const C_GROUP *)cg2)->nGroupNumber;
}
/*********************************************************************************/
int AddTGroups2BnStruct( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, T_GROUP_INFO *tgi  )
{
    int ret = 0;
    /* ret = ReInitBnStruct( pBNS ); */
    if ( tgi && tgi->num_t_groups && tgi->t_group ) {
        int         i, k, endpoint, centerpoint, fictpoint;
        int         num_tg       = tgi->num_t_groups;
        int         num_edges    = pBNS->num_edges;
        int         num_vertices = pBNS->num_vertices;
        BNS_VERTEX *vert_ficpoint, *ver_ficpont_prev;  /* fictitious vertex describing t-group */
        BNS_VERTEX *vert_endpoint;
        BNS_EDGE   *edge;      /* edge between that vertex and the tautomeric endpoint */
        int        nMaxTGroupNumber = 0;
        ENDPOINT_INFO eif;

        /* Debug: check overflow */
        if ( num_vertices + num_tg >= pBNS->max_vertices ) {
            return BNS_VERT_EDGE_OVFL;
        }
        /* find the largest t-group ID */
        for ( i = 0; i < num_tg; i ++ ) {
            if ( tgi->t_group[i].nGroupNumber > nMaxTGroupNumber ) {
                nMaxTGroupNumber = tgi->t_group[i].nGroupNumber;
            }
        }
        /* since t-group IDs may be not contiguous, clear all vertices that will be added.
           all-zeroes-vertex will be ignored by the BNS
        */
        memset( pBNS->vert+num_vertices, 0, nMaxTGroupNumber*sizeof(pBNS->vert[0]) );
        /* Make sure the last t-group has the largest t-group ID:
           this is necessary to correctly add new edges and vertices for testing augmenting paths
        */
#if( bRELEASE_VERSION != 1 )
        insertions_sort( tgi->t_group, num_tg, sizeof(tgi->t_group[0]), CompTGroupNumber );
        for ( i = 1; i < num_tg; i ++ ) {
            if ( 1 != tgi->t_group[i].nGroupNumber - tgi->t_group[i-1].nGroupNumber ) {
                return BNS_BOND_ERR;
            }
        }
#else
        if ( nMaxTGroupNumber != tgi->t_group[num_tg-1].nGroupNumber ) {
            insertions_sort( tgi->t_group, num_tg, sizeof(tgi->t_group[0]), CompTGroupNumber );
        }
#endif
        /* initialize new fictitious vertices */
        ver_ficpont_prev = pBNS->vert+num_vertices - 1;

        for ( i = 0; i < num_tg; i ++, ver_ficpont_prev = vert_ficpoint ) {
            /*
              vert_ficpoint-1 is the last vertex;
              vert_ficpoint   is the vertex that is being added
              Note: nGroupNumber are not contiguous
            */
            vert_ficpoint                = pBNS->vert+num_vertices + tgi->t_group[i].nGroupNumber - 1;
            vert_ficpoint->iedge         = ver_ficpont_prev->iedge + ver_ficpont_prev->max_adj_edges;
            vert_ficpoint->max_adj_edges = tgi->t_group[i].nNumEndpoints+BNS_ADD_EDGES+BNS_ADD_SUPER_TGROUP;
            vert_ficpoint->num_adj_edges = 0;
            vert_ficpoint->st_edge.flow  = vert_ficpoint->st_edge.flow0  = 0;
            vert_ficpoint->st_edge.cap   = vert_ficpoint->st_edge.cap0   = 0;
            vert_ficpoint->type          = BNS_VERT_TYPE_TGROUP;
        }

        for ( endpoint = 0; endpoint < num_atoms; endpoint ++ ) {
            if ( !at[endpoint].endpoint )
                continue;
            fictpoint = at[endpoint].endpoint + num_vertices - 1;
            vert_ficpoint = pBNS->vert + fictpoint;
            vert_endpoint = pBNS->vert + endpoint;
            /* Debug: check overflow */
            if ( fictpoint >= pBNS->max_vertices ||
                 num_edges >= pBNS->max_edges    ||
                 vert_ficpoint->num_adj_edges >= vert_ficpoint->max_adj_edges ||
                 vert_endpoint->num_adj_edges >= vert_endpoint->max_adj_edges ) {
                ret = BNS_VERT_EDGE_OVFL;
                break;
            }
            /* obtain donor/acceptor info */
            if ( !nGetEndpointInfo( at, endpoint, &eif ) ) {
#if( KETO_ENOL_TAUT == 1 )
                if ( !((tgi->bTautFlags & TG_FLAG_KETO_ENOL_TAUT ) &&
                        nGetEndpointInfo_KET( at, endpoint, &eif )) )
#endif
                {
                    ret = BNS_BOND_ERR;
                    break;
                }
            }

            vert_endpoint->type |= BNS_VERT_TYPE_ENDPOINT;

            /* set capacity = 1 to the edges from the endpoint to the centerpoint(s) */
            for ( k = 0; k < vert_endpoint->num_adj_edges; k ++ ) {
                int iedge = vert_endpoint->iedge[k];
                if ( !pBNS->edge[iedge].cap ) {
                    /* single bond, possibly between endpoint and centerpoint */
                    centerpoint = (pBNS->edge[iedge].neighbor12 ^ endpoint);
                    if ( centerpoint < pBNS->num_atoms &&
                         pBNS->vert[centerpoint].st_edge.cap >= 1 ) {
                        int bond_type = (at[endpoint].bond_type[k] & BOND_TYPE_MASK);
                        if (bond_type == BOND_TAUTOM  ||
                            bond_type == BOND_ALTERN  ||
                            bond_type == BOND_ALT12NS ||
                            bond_type == BOND_SINGLE ) {
                            pBNS->edge[iedge].cap = 1;
                        }
                    }
                }
            }
            /* create a new edge connecting endpoint to the new fictitious t-group vertex vert_ficpoint */
            edge = pBNS->edge + num_edges;
            edge->cap       = 1;
            edge->flow      = 0;
            edge->pass      = 0;
#if( RESET_EDGE_FORBIDDEN_MASK == 1 )
            edge->forbidden &= pBNS->edge_forbidden_mask;
#endif
            /* later include case when the charge change allows the endpoint to become tautomeric */
            /* mark endoint having moveable H atom with flow=1 */

            /* -- old "no charges" version -- */
            /* if (at[endpoint].chem_bonds_valence == at[endpoint].valence) */
            /* -- the following line takes charges into account -- */
            if ( eif.cDonor ) /* means the endpoint has an H-atom to donate */
            {
                /* increment edge flow */
                edge->flow ++;
                /* increment one vertex st-flow & cap */
                vert_ficpoint->st_edge.flow ++;
                vert_ficpoint->st_edge.cap ++;
                /* increment another vertex st-flow & cap */
                vert_endpoint->st_edge.flow ++;
                vert_endpoint->st_edge.cap ++;
            }
            /* connect edge to endpoint and fictpoint and increment the counters of neighbors and edges */
            edge->neighbor1    = endpoint; /* the smallest out of v1=endopoint and v2=num_vertices */
            edge->neighbor12   = endpoint ^ fictpoint; /* v1 ^ v2 */
            vert_endpoint->iedge[vert_endpoint->num_adj_edges] = num_edges;
            vert_ficpoint->iedge[vert_ficpoint->num_adj_edges] = num_edges ++;
            edge->neigh_ord[0] = vert_endpoint->num_adj_edges ++;
            edge->neigh_ord[1] = vert_ficpoint->num_adj_edges ++;
            edge->cap0  = edge->cap;
            edge->flow0 = edge->flow;
        }

        pBNS->num_edges     = num_edges;
        pBNS->num_vertices += nMaxTGroupNumber;
        pBNS->num_t_groups  = num_tg;

    }
    return ret;
}
/******************************************************************************************/
/*#if( MOVE_CHARGES == 1 )*/ /* { */
/*********************************************************************************/
int AddCGroups2BnStruct( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi  )
{
    int ret = 0;
    /* ret = ReInitBnStruct( pBNS ); */
    if ( cgi && cgi->num_c_groups && cgi->c_group ) {
        int         i, k, c_point, centerpoint, fictpoint;
        int         num_cg       = cgi->num_c_groups;
        int         num_edges    = pBNS->num_edges;
        int         num_vertices = pBNS->num_vertices;
        BNS_VERTEX *vert_ficpoint, *ver_ficpont_prev;  /* fictitious vertex describing charge c-group */
        BNS_VERTEX *vertex_cpoint;
        BNS_EDGE   *edge;      /* edge between that vertex and the tautomeric c_point */
        int        nMaxCGroupNumber = 0;

        /* Debug: check overflow */
        if ( num_vertices + num_cg >= pBNS->max_vertices ) {
            return BNS_VERT_EDGE_OVFL;
        }
        /* find the largest t-group ID */
        for ( i = 0; i < num_cg; i ++ ) {
            if ( cgi->c_group[i].nGroupNumber > nMaxCGroupNumber ) {
                nMaxCGroupNumber = cgi->c_group[i].nGroupNumber;
            }
        }
        /* since t-group IDs may be not contiguous, clear all vertices that will be added.
           all-zeroes-vertex will be ignored by the BNS
        */
        memset( pBNS->vert+num_vertices, 0, nMaxCGroupNumber*sizeof(pBNS->vert[0]) );
        /* Make sure the last t-group has the largest t-group ID:
           this is necessary to correctly add new edges and vertices for testing augmenting paths
        */
#if( bRELEASE_VERSION != 1 )
        insertions_sort( cgi->c_group, num_cg, sizeof(cgi->c_group[0]), CompCGroupNumber );
        for ( i = 1; i < num_cg; i ++ ) {
            if ( 1 != cgi->c_group[i].nGroupNumber - cgi->c_group[i-1].nGroupNumber ) {
                return BNS_BOND_ERR;
            }
        }
#else
        if ( nMaxCGroupNumber != cgi->c_group[num_cg-1].nGroupNumber ) {
            insertions_sort( cgi->c_group, num_cg, sizeof(cgi->c_group[0]), CompCGroupNumber );
        }
#endif
        /**************************************/
        /* initialize new fictitious vertices */
        /* representing c-point groups        */
        /**************************************/
        ver_ficpont_prev = pBNS->vert+num_vertices - 1;
        
        for ( i = 0; i < num_cg; i ++, ver_ficpont_prev = vert_ficpoint ) {
            /*
              vert_ficpoint-1 is the last vertex;
              vert_ficpoint   is the being added vertex
              Note: nGroupNumber are not contiguous
            */
            vert_ficpoint                = pBNS->vert+num_vertices + cgi->c_group[i].nGroupNumber - 1;
            vert_ficpoint->iedge         = ver_ficpont_prev->iedge + ver_ficpont_prev->max_adj_edges;
            vert_ficpoint->max_adj_edges = cgi->c_group[i].num_CPoints+BNS_ADD_EDGES;
            vert_ficpoint->num_adj_edges = 0;
            vert_ficpoint->st_edge.flow  = vert_ficpoint->st_edge.flow0  = 0;
            vert_ficpoint->st_edge.cap   = vert_ficpoint->st_edge.cap0   = 0;
            vert_ficpoint->type          = BNS_VERT_TYPE_C_GROUP;
        }

        /************************************************/
        /* connect c-points to the fictitious vertices  */
        /* representing c-point groups; set caps, flows */
        /************************************************/
        for ( c_point = 0; c_point < num_atoms; c_point ++ ) {
            if ( !at[c_point].c_point )
                continue;
            fictpoint = at[c_point].c_point + num_vertices - 1; /* c-group vertex index */
            vert_ficpoint = pBNS->vert + fictpoint; /* c-group vertex */
            vertex_cpoint = pBNS->vert + c_point;   /* c_point vertex */
            /* Debug: check overflow */
            if ( fictpoint >= pBNS->max_vertices ||
                 num_edges >= pBNS->max_edges    ||
                 vert_ficpoint->num_adj_edges >= vert_ficpoint->max_adj_edges ||
                 vertex_cpoint->num_adj_edges >= vertex_cpoint->max_adj_edges ) {
                ret = BNS_VERT_EDGE_OVFL;
                break;
            }
            vertex_cpoint->type |= BNS_VERT_TYPE_C_POINT;
#if( FIX_CPOINT_BOND_CAP != 1 )  /* { */
            /* set capacity = 1 to the edges from the c_point to the centerpoint(s)     */
            /* if their current capacity is zero                                        */
            /* the centerpoint is any adjacent atom that is adjacent to a multiple bond */
            for ( k = 0; k < vertex_cpoint->num_adj_edges; k ++ ) {
                int iedge = vertex_cpoint->iedge[k];
                if ( !pBNS->edge[iedge].cap ) {
                    /* single bond, possibly between c_point and centerpoint */
                    centerpoint = (pBNS->edge[iedge].neighbor12 ^ c_point);
                    if ( centerpoint < pBNS->num_atoms &&
                         pBNS->vert[centerpoint].st_edge.cap >= 1 ) {
                        int bond_type = (at[c_point].bond_type[k] & BOND_TYPE_MASK);
                        if ( bond_type == BOND_TAUTOM ||
                             bond_type == BOND_ALTERN ||
                             bond_type == BOND_SINGLE ) {
                            pBNS->edge[iedge].cap = 1;
                        }
                    }
                }
            }
#endif /* } FIX_CPOINT_BOND_CAP */
            /* create a new edge connecting c_point to the new fictitious c-group vertex vert_ficpoint */
            edge = pBNS->edge + num_edges;
            edge->cap       = 1;
            edge->flow      = 0;
            edge->pass      = 0;
#if( RESET_EDGE_FORBIDDEN_MASK == 1 )
            edge->forbidden &= pBNS->edge_forbidden_mask;
#endif
            /* mark edge to c-point having NO moveable charge with flow=1 */
            if ( !CHARGED_CPOINT(at,c_point) ) {
                /* increment edge flow */
                edge->flow ++;
                /* increment c-group vertex st-flow & cap */
                vert_ficpoint->st_edge.flow ++;
                vert_ficpoint->st_edge.cap ++;
                /* increment c-point vertex st-flow & cap */
                vertex_cpoint->st_edge.flow ++;
                vertex_cpoint->st_edge.cap ++;
            }
#if( FIX_CPOINT_BOND_CAP == 1 ) /* { */
            /* set capacity = 1 to the edges from the c_point to the centerpoint(s)     */
            /* if their current capacity is zero                                        */
            /* the centerpoint is any adjacent atom that is adjacent to a multiple bond */
            for ( k = 0; k < vertex_cpoint->num_adj_edges; k ++ ) {
                int iedge = vertex_cpoint->iedge[k];
                VertexFlow  nNewCap = vertex_cpoint->st_edge.cap;
                centerpoint = (pBNS->edge[iedge].neighbor12 ^ c_point);
                if ( !pBNS->edge[iedge].cap ) {
                    /* single bond, possibly between c_point and centerpoint */
                    if ( centerpoint < pBNS->num_atoms &&
                         pBNS->vert[centerpoint].st_edge.cap >= 1 ) {
                        nNewCap = inchi_min( pBNS->vert[centerpoint].st_edge.cap, nNewCap );
                        nNewCap = inchi_min( nNewCap, MAX_BOND_EDGE_CAP );
                        pBNS->edge[iedge].cap = nNewCap;
                    }
#if( FIX_CPOINT_BOND_CAP2 == 1 ) /* multiple bond */
                    else
                    if ( centerpoint < pBNS->num_atoms &&
                         edge->flow && pBNS->edge[iedge].cap < MAX_BOND_EDGE_CAP ) {
                        pBNS->edge[iedge].cap ++;
                    }
#endif
                }
            }
#endif  /* } FIX_CPOINT_BOND_CAP */
            /* connect edge to c_point and fictpoint and increment the counters of neighbors and edges */
            edge->neighbor1    = c_point; /* the smallest out of v1=endopoint and v2=num_vertices */
            edge->neighbor12   = c_point ^ fictpoint; /* v1 ^ v2 */
            vertex_cpoint->iedge[vertex_cpoint->num_adj_edges] = num_edges;
            vert_ficpoint->iedge[vert_ficpoint->num_adj_edges] = num_edges ++;
            edge->neigh_ord[0] = vertex_cpoint->num_adj_edges ++;
            edge->neigh_ord[1] = vert_ficpoint->num_adj_edges ++;
            edge->cap0  = edge->cap;
            edge->flow0 = edge->flow;

        }

        pBNS->num_edges     = num_edges;
        pBNS->num_vertices += nMaxCGroupNumber;
        pBNS->num_c_groups  = num_cg;

    }
    return ret;
}
/*#endif*/ /* } MOVE_CHARGES == 1 */

/*********************************************************************************/
void ClearAllBnDataVertices( Vertex *v, Vertex value, int size )
{
    int i;
    for ( i = 0; i < size; i ++ ) {
        v[i] = value;
    }
}

/*********************************************************************************/
void ClearAllBnDataEdges( Edge *e, Vertex value, int size )
{
    int i;
    for ( i = 0; i < size; i ++ ) {
        e[i][0] = value;
    }
}
/*********************************************************************************/
BN_DATA *DeAllocateBnData( BN_DATA *pBD )
{
    if ( pBD ) {
        if ( pBD->BasePtr )
            inchi_free( pBD->BasePtr );
        if ( pBD->SwitchEdge )
            inchi_free( pBD->SwitchEdge );
        if ( pBD->Tree )
            inchi_free( pBD->Tree );
        if ( pBD->ScanQ )
            inchi_free( pBD->ScanQ );
        if ( pBD->Pu )
            inchi_free( pBD->Pu );
        if ( pBD->Pv )
            inchi_free( pBD->Pv );
#if( BNS_RAD_SEARCH == 1 )
        if ( pBD->RadEndpoints ) {
            inchi_free( pBD->RadEndpoints );
        }
        if ( pBD->RadEdges ) {
            inchi_free( pBD->RadEdges );
        }
#endif
        inchi_free( pBD );
    }
    return NULL;
}

/*********************************************************************************/
BN_DATA *AllocateAndInitBnData( int max_num_vertices ) 
{
    BN_DATA *pBD = NULL;
    int      max_len_Pu_Pv;
    max_num_vertices = 2*max_num_vertices+2;
    max_len_Pu_Pv    = max_num_vertices/2+1;
    max_len_Pu_Pv   += max_len_Pu_Pv % 2; /* even length */
    if ( !(pBD             = (BN_DATA *) inchi_calloc( 1, sizeof(BN_DATA) ) ) ||
         !(pBD->BasePtr    = (Vertex *)  inchi_calloc( max_num_vertices, sizeof(Vertex) ) ) ||
         !(pBD->SwitchEdge = (Edge   *)  inchi_calloc( max_num_vertices, sizeof(Edge  ) ) ) ||
         !(pBD->Tree       = (S_CHAR *)  inchi_calloc( max_num_vertices, sizeof(S_CHAR) ) ) ||
         !(pBD->ScanQ      = (Vertex *)  inchi_calloc( max_num_vertices, sizeof(Vertex) ) ) ||
         !(pBD->Pu         = (Vertex *)  inchi_calloc( max_len_Pu_Pv,    sizeof(Vertex) ) ) ||
#if( BNS_RAD_SEARCH == 1 )
         !(pBD->RadEndpoints = (Vertex *)   inchi_calloc( max_len_Pu_Pv, sizeof(Vertex) ) ) ||
         !(pBD->RadEdges     = (EdgeIndex*) inchi_calloc( max_len_Pu_Pv, sizeof(EdgeIndex) ) ) ||
#endif
         !(pBD->Pv         = (Vertex *)  inchi_calloc( max_len_Pu_Pv,    sizeof(Vertex) ) ) 
       ) {
        pBD = DeAllocateBnData( pBD );
    } else {
        /* Initialize data */
        ClearAllBnDataEdges(pBD->SwitchEdge, NO_VERTEX, max_num_vertices);
        ClearAllBnDataVertices(pBD->BasePtr, NO_VERTEX, max_num_vertices);
        memset(pBD->Tree, TREE_NOT_IN_M, max_num_vertices);
        pBD->QSize = -1;
        pBD->max_len_Pu_Pv    = max_len_Pu_Pv;
        pBD->max_num_vertices = max_num_vertices;
#if( BNS_RAD_SEARCH == 1 )
        pBD->nNumRadEndpoints = 0;
#endif
    }
    return pBD;
}
/*********************************************************************************/
int ReInitBnData( BN_DATA *pBD )
{
    int i, ret = 0;
    Vertex u, v;
    if ( pBD ) {
        if ( !pBD->ScanQ ) {
            ret += 2;
        }
        if ( !pBD->BasePtr ) {
            ret += 4;
        }
        if ( !pBD->SwitchEdge ) {
            ret += 8;
        }
        if ( !pBD->Tree ) {
            ret += 16;
        }
        if ( !ret ) {
            for ( i = 0; i <= pBD->QSize; i ++ ) {
                u = pBD->ScanQ[i];
                v = prim(u);
                pBD->BasePtr[u]     =
                pBD->BasePtr[v]     =
                pBD->SwitchEdge_Vert1(u) =
                pBD->SwitchEdge_Vert1(v) = NO_VERTEX;
                pBD->Tree[u]        =
                pBD->Tree[v]        = TREE_NOT_IN_M;
            }
        }
        pBD->QSize = -1;
        if ( !pBD->Pu ) {
            ret += 32;
        }
        if ( !pBD->Pv ) {
            ret += 64;
        }
    } else {
        ret += 1;
    }
    
    return ret;
}


/*******************************************************/
int GetVertexDegree( BN_STRUCT* pBNS, Vertex v )
{
    int i = v / 2 - 1;
    if ( i >= 0 ) {
        if ( pBNS->vert[i].st_edge.cap > 0 ) {
            return pBNS->vert[i].num_adj_edges+1; /* add 1 neighbor for s or t */
        } else {
            return 0;   /* since the edge s-v has zero capacity, we ignore vertex v */
        }
    } else {
        return pBNS->num_vertices;
    }
}
/*********************************************************************************/

Vertex Get2ndEdgeVertex( BN_STRUCT* pBNS, Edge uv )
{
    /*
    Vertex ret;
    */
    if ( uv[1] >= 0 ) {
        /* -- debug --
        if ( uv[1] > pBNS->num_edges || uv[0] > 2*pBNS->num_vertices+3 ) {
            int stop = 1;
        }
        ret = ((uv[0]-2) ^ (2*pBNS->edge[uv[1]].neighbor12+1)) + 2;
        if ( ret > 2*pBNS->num_vertices+3 ) {
            int stop = 1;
        }
        return ret;
          -- end debug -- */
        return ((uv[0]-2) ^ (2*pBNS->edge[uv[1]].neighbor12+1)) + 2;
        /*short u = uv[0]-FIRST_INDX; */
        /*short t = 2*(((u / 2 - 1) ^ pBNS->edge[uv[1]].neighbor12) + 1) + ((u+1) & 1) + FIRST_INDX; */
        /*return t; */
    }
    if ( uv[0] <= 1 )
        return -(1 + uv[1]); /* vertex1 is s or t, return x or y */
    else
        return uv[0] % 2; /* vertex1 is x or y, return s or t; never happens? -- NSC 3737, 7634,... */
}    
/*********************************************************************************/
Vertex GetVertexNeighbor( BN_STRUCT* pBNS, Vertex v, int neigh, EdgeIndex *iedge )
{
    /* neigh = 0 => the neighbor is s or t except case when v is s or t. */
    /* v= FIRST_INDX or FIRST_INDX+1: v is s or t respectively */
    int i, neighbor;
    if ( (i = v - 2) >= 0 ) {
        /* neighbor of x or y */
        if ( neigh ) {
            neigh --;
            /* x or y */
            *iedge = pBNS->vert[i/2].iedge[neigh];
            if ( !(pBNS->edge[*iedge].cap & EDGE_FLOW_MASK) || IS_FORBIDDEN(pBNS->edge[*iedge].forbidden, pBNS) ) {
                return NO_VERTEX;
            }
            neighbor = (i ^ (2 * pBNS->edge[*iedge].neighbor12 + 1)) + 2; /* parity opposite to v parity */
        } else {
            /* neighbor of x or y is s or t */
            neighbor = (v & 1); /* s or t, same parity as v */
            *iedge   = -( neighbor + 1 );
        }
    } else {
        /* neighbor of s or t: x or y, same parity as v */
        if ( !(pBNS->vert[neigh].st_edge.cap & EDGE_FLOW_ST_MASK) ) {
            return NO_VERTEX;
        }
        neighbor = 2*neigh + 2 + (v & 1); /* parity same as the parity of v */
        *iedge   = -( neighbor + 1 );
    }
    return neighbor;
}
/*********************************************************************************/
int GetEdgePointer( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv, BNS_EDGE **uv, S_CHAR *s_or_t )
{
    int i = u / 2 - 1;
    int j = v / 2 - 1;
    int bBackward = BNS_WRONG_PARMS;
    *uv = NULL;
    if ( i >= 0 ) {
        /* u is an atom */
        if ( j >= 0 ) {
            /* v is an atom */
            if ( (u+v)%2 ) {
                *uv = pBNS->edge+iuv;
                bBackward = ( u & 1 );
                *s_or_t = 0;
            }
        } else
        /* v is s or t */
        if ( v >= 0 && !((u+v)%2) ) {
            *uv = (BNS_EDGE*)&pBNS->vert[i].st_edge;
            bBackward = !(v & 1);
            *s_or_t = v+3; /* 3=> v=s, 4=> v=t */
        }
    } else 
    if ( j >= 0 ) {
        /* u is s or t */
        if ( u >= 0 && !((u+v)%2) ) {
            /* v is an atom */
            *uv = (BNS_EDGE*)&pBNS->vert[j].st_edge;
            bBackward = (u & 1 );
            *s_or_t = u+1; /* 1=> u=s, 2=> u=t */
        }
    }
    return bBackward;
}
/*********************************************************************************/
int AugmentEdge( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv, int delta, S_CHAR bReverse, int bChangeFlow )
{
    int f, flow, ret=0;

    BNS_ST_EDGE    *pst_edge;
    BNS_EDGE       *pedge;
    S_CHAR          s_or_t;
    int bBackward = GetEdgePointer( pBNS, u, v, iuv, &pedge, &s_or_t );
    if ( !IS_BNS_ERROR(bBackward) ) {
        if ( bBackward ) {
            delta = -delta;
        }
        
        if ( s_or_t ) {
            pst_edge = ( BNS_ST_EDGE *) pedge;
            flow = pst_edge->flow;
            f    =  (flow & EDGE_FLOW_ST_MASK) + delta; /* new flow */
            if ( !delta ) {
                /*((BNS_ST_EDGE *)pedge)->flow = pst_edge->flow & ~EDGE_FLOW_ST_PATH;*/
                pst_edge->flow = pst_edge->flow & ~EDGE_FLOW_ST_PATH;
            } else {
                int cap = pst_edge->cap;
                if ( f < 0 || f > cap ) {
                    ret = BNS_WRONG_PARMS;
                } else {
                    if ( !(bChangeFlow & BNS_EF_CHNG_FLOW) ) {
                        f -= delta; /* do not actually change the flow, only find the augmenting path */
                    } else
                    if ( delta ) {
                        /*((BNS_ST_EDGE *)pedge)->pass ++;*/
                        pst_edge->pass ++;
                    }
                    flow = (flow & ~(EDGE_FLOW_ST_PATH | EDGE_FLOW_ST_MASK)) + f;
                    /*((BNS_ST_EDGE *)pedge)->flow = flow;*/
                    pst_edge->flow = flow;
                    /*((BNS_ST_EDGE *)pedge)->delta += delta; */
                    if ( bReverse ) {
                        /* u <- v; Note: in case of bReverse s_or_t has actually been determined
                                         for the u' <- v' pair; therefore s and t should be switched
                                         in order to correctly determine the 1st or the last atom
                                         on the augmenting path.
                        */
                        switch( s_or_t ) {
                        case 1: /* u = t: t<-v, v is the last vertex */
                            ALTP_END_ATOM(pBNS->alt_path) = v / 2 - 1;
                            break;
                        case 2: /* u = s: s<-v, error */
                            ret = BNS_WRONG_PARMS;
                            break;
                        case 3: /* v = t: u<-t, error */
                            ret = BNS_WRONG_PARMS;
                            break;
                        case 4: /* v = s: u<-s, u is the first vertex */
                            ALTP_START_ATOM(pBNS->alt_path) = u / 2 - 1;
                            ALTP_DELTA(pBNS->alt_path) = delta;
                            break;
                        default:
                            ret = BNS_WRONG_PARMS;
                            break;
                        }
                    } else {
                        /* u -> v */
                        switch( s_or_t ) {
                        case 1: /* u = s: s->v, v is the first vertex */
                            ALTP_START_ATOM(pBNS->alt_path) = v / 2 - 1;
                            ALTP_DELTA(pBNS->alt_path) = delta;
                            break;
                        case 2: /* u = t: t->v, error */
                            ret = BNS_WRONG_PARMS;
                            break;
                        case 3: /* v = s: u->s, error */
                            ret = BNS_WRONG_PARMS;
                            break;
                        case 4: /* v = t: u->t, u is the last vertex */
                            ALTP_END_ATOM(pBNS->alt_path) = u / 2 - 1;
                            break;
                        default:
                            ret = BNS_WRONG_PARMS;
                            break;
                        }
                    }
                }
            }
        } else {
            f = (pedge->flow & EDGE_FLOW_MASK) + delta; 
            if ( !delta ) {
                pedge->flow &= ~EDGE_FLOW_PATH;
            } else {
                if ( f < 0 || f > pedge->cap ) {
                    ret = BNS_WRONG_PARMS;
                } else {
                    AT_NUMB   iu = u / 2 - 1;
                    AT_NUMB   iv = v / 2 - 1;
                    int       indx;
                    if ( !(bChangeFlow & BNS_EF_CHNG_FLOW) ) {
                        f -= delta; /* do not actually change the flow, only find the augmenting path */
                    } else
                    if ( delta ) {
                        pedge->pass ++;
                    }
                    pedge->flow = (pedge->flow & ~(EDGE_FLOW_PATH | EDGE_FLOW_MASK)) | f;
                    if ( ALTP_MAY_ADD(pBNS->alt_path) ) {
                        /* bReverse? u <- v : u -> v */
                        indx = bReverse? (pedge->neighbor1 == iv) : (pedge->neighbor1 == iu);
                        ALTP_CUR_THIS_ATOM_NEIGHBOR(pBNS->alt_path) = pedge->neigh_ord[1-indx];
                        ALTP_CUR_NEXT_ATOM_NEIGHBOR(pBNS->alt_path) = pedge->neigh_ord[indx];
                        ALTP_NEXT(pBNS->alt_path);
                    } else {
                        ALTP_OVERFLOW(pBNS->alt_path) = 1;
                        ret = BNS_ALTPATH_OVFL;
                    }
                }
            }
        }
        return ret? ret : f;

    }
    return bBackward;
}
/*********************************************************************************/
/* find residual capacity and mark the edge as belonging to the augmenting path  */
/*********************************************************************************/
int rescap_mark( BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv )
{
    BNS_ST_EDGE    *pst_edge;
    BNS_EDGE       *pedge;

    int    f, flow;
    S_CHAR s_or_t;
    int    bBackward = GetEdgePointer( pBNS, u, v, iuv, &pedge, &s_or_t );
    
    if ( !IS_BNS_ERROR( bBackward ) ) {
        
        if ( s_or_t ) {
            pst_edge = (BNS_ST_EDGE *)pedge;
            flow = pst_edge->flow;
            f    = (flow & EDGE_FLOW_ST_MASK);
            if ( !bBackward ) {
                f = (int)pst_edge->cap - f;
            }
            if ( flow & EDGE_FLOW_ST_PATH ) {
                pBNS->bNotASimplePath ++;
                f /= 2;   /* this is the second time we pass the same edge: reduce flow by a factor of 2 */
            } else {
                pst_edge->flow |= EDGE_FLOW_ST_PATH; /* mark the edge */
            }
        } else {
            flow = pedge->flow;
            f    = flow & EDGE_FLOW_MASK;
            if ( !bBackward ) {
                f = (int)pedge->cap - f;
            }
            if ( flow & EDGE_FLOW_PATH ) {
                f /= 2;    /* this is the second time we pass the same edge: reduce flow by a factor of 2 */
                pBNS->bNotASimplePath ++;
            } else {
                pedge->flow |= EDGE_FLOW_PATH;  /* mark the edge */
            }
        }
        return f;

    }
    return bBackward;
}
/********************************************************************************
Get previous vertex in the searched path
z is SwitchEdge_Vert2(y) != y. Go backward from z to y
*********************************************************************************/
Vertex GetPrevVertex( BN_STRUCT* pBNS, Vertex y, Edge *SwitchEdge, EdgeIndex *iuv )
{
    Vertex w, z, x2, y2, n;
    EdgeIndex iwy;

    w   = SwitchEdge_Vert1(y);
    z   = SwitchEdge_Vert2(y);
    iwy = SwitchEdge_IEdge(y);
    if ( z == y ) {
        *iuv = iwy;
        return w;
    }
    x2 = prim(y);
    y2 = prim(z);
    n = 0;
    while ( y2 != NO_VERTEX ) {
        w   = SwitchEdge_Vert1(y2);
        z   = SwitchEdge_Vert2(y2);
        iwy = SwitchEdge_IEdge(y2);
        if ( w == x2 ) {
            *iuv = iwy;
            /*return z; */
            return (y + z)%2? z : prim(z);
        }
        n ++;
#ifdef _DEBUG
        if ( n ) {
            int stop = 1;
        }
#endif
        if ( w == y2 )
            return NO_VERTEX;
        y2 = w;
    }
    return y2;
}
#define CHECK_TACN  1
/*********************************************************************************
   the purpose is to avoid paths
     (H-group)[u]---atom[v]---((-)-cgroup)[w],  where
   atom[v] is not acidic and (-) and H are not interchangeable without
   explicit bond tautomerism.
   It is important that acidic atoms are only O,S,Se,Te and should have
   only one chemical bond. Only because of this an early rejection of
   the vertex v (before it gets on SCANQ) is possible.

   CHECK_TACN == 1 prohibits replacing (-) on N with H unless H can be moved to N
   along an alternating path from another heteroatom (t-group will be detected).

**********************************************************************************/
#if( FIX_TACN_POSSIBLE_BUG == 1 ) /* { */
/*********************************************************************************/
int bIgnoreVertexNonTACN_atom( BN_STRUCT* pBNS, Vertex u, Vertex v )
{
#define TYPE_T   1   /* t-group [also called H-group] */
#define TYPE_CN  2   /* (-)c-group */
    int    i, degree, ret, u_is_taut=0, w_is_taut, num_allowed=0, num_found_groups=0;
    Vertex w;
    EdgeIndex ivw;
    if ( !pBNS->type_TACN || u <= 1 || v <= 1 ||
         (pBNS->vert[v/2-1].type & pBNS->type_TACN) ) {
        return 0; /* add/remove H(+) is allowed for acidic atoms */
    }
    if ( !pBNS->type_T || !pBNS->type_CN )
        return 0; /* should not happen */
    u_is_taut = ((pBNS->vert[u/2-1].type & pBNS->type_T)  == pBNS->type_T )? TYPE_T  :
                ((pBNS->vert[u/2-1].type & pBNS->type_CN) == pBNS->type_CN)? TYPE_CN : 0;
    if ( u_is_taut ) {
        /* u is either t-group vertex or (-) c-group */
        degree = GetVertexDegree( pBNS, v );
        for ( i = 0; i < degree; i ++ ) {
            /* v = vert[u].neighbor[i]; */
            w = GetVertexNeighbor( pBNS, v, i, &ivw );
            if ( w == NO_VERTEX || w <= 1 ) {
                continue; /* the atom has only single bonds or it is s or t, ignore it */
            }
            if ( w != u && (ret = rescap(pBNS, v, w, ivw)) > 0 ) {
                num_allowed ++;
                w_is_taut = ((pBNS->vert[w/2-1].type & pBNS->type_CN) == pBNS->type_CN)? TYPE_CN :
                            ((pBNS->vert[w/2-1].type & pBNS->type_T)  == pBNS->type_T )? TYPE_T  : 0;
                if ( (u_is_taut | w_is_taut) == (TYPE_T | TYPE_CN) ) {
                    num_found_groups ++;
                }
            }
        }
        if ( num_found_groups && num_allowed == 1 ) {
            return 1; /* reject */
        }
    }
    return 0;
#undef TYPE_T
#undef TYPE_CN
}
#else  /* } FIX_TACN_POSSIBLE_BUG { */
/*********************************************************************************/
int bIgnoreVertexNonTACN_atom( BN_STRUCT* pBNS, Vertex u, Vertex v )
{
    int    i, degree, ret, u_is_taut=0, num_allowed=0, num_found_groups=0;
    Vertex w;
    EdgeIndex ivw;
    if ( !pBNS->type_TACN || u <= 1 || v <= 1 ||
         (pBNS->vert[v/2-1].type & pBNS->type_TACN) ) {
        return 0; /* add/remove H(+) is allowed for acidic atoms */
    }
    if ( !pBNS->type_T || !pBNS->type_CN )
        return 0; /* should not happen */
    if ( (u_is_taut  = (pBNS->vert[u/2-1].type & pBNS->type_T)  == pBNS->type_T) ||
         (             (pBNS->vert[u/2-1].type & pBNS->type_CN) == pBNS->type_CN)  ) {
        /* u is either t-group vertex or (-) c-group */
        degree = GetVertexDegree( pBNS, v );
        for ( i = 0; i < degree; i ++ ) {
            /* v = vert[u].neighbor[i]; */
            w = GetVertexNeighbor( pBNS, v, i, &ivw );
            if ( w == NO_VERTEX || w <= 1 ) {
                continue; /* the atom has only single bonds or it is s or t, ignore it */
            }
            if ( w != u && (ret = rescap(pBNS, v, w, ivw)) > 0 ) {
                num_allowed ++;
                if ( (u_is_taut? ((pBNS->vert[w/2-1].type & pBNS->type_CN) == pBNS->type_CN) :
                                 ((pBNS->vert[w/2-1].type & pBNS->type_T)  == pBNS->type_T ) ) ) {
                    num_found_groups ++;
                }
            }
        }
        if ( num_found_groups && num_allowed == 1 ) {
            return 1; /* reject */
        }

    }
    return 0;
}
#endif /* } FIX_TACN_POSSIBLE_BUG */
/*********************************************************************************
   the purpose is to avoid paths
     (H-group)[u]---atom[v]---((-)-cgroup)[w],  where
   atom[v] is not acidic and (-) and H are not interchangeable without
   explicit bond tautomerism.
   It is important that acidic atoms are only O,S,Se,Te and should have
   only one chemical bond. Only because of this an early rejection of
   the vertex v (before it gets on SCANQ) is possible.
**********************************************************************************/
#if( FIX_TACN_POSSIBLE_BUG == 1 ) /* { */
/*********************************************************************************/
int bIgnoreVertexNonTACN_group( BN_STRUCT* pBNS, Vertex v, Vertex w, Edge *SwitchEdge )
{
#define TYPE_T   1   /* t-group [also called H-group] */
#define TYPE_CN  2   /* (-)c-group */
    int    u_is_taut=0, w_is_taut=0;
    Vertex u;
    EdgeIndex iuv;
    if ( v <= 1 || w <= 1 )
        return 0;
#if ( CHECK_TACN == 1 )
    if ( !pBNS->type_TACN ||
         (pBNS->vert[v/2-1].type & pBNS->type_TACN) ) {
        return 0;
    }
    if ( !pBNS->type_T || !pBNS->type_CN )
        return 0; /* should not happen */
#endif
    u = GetPrevVertex( pBNS, v, SwitchEdge, &iuv );
    /*
    u   = SwitchEdge_Vert1(v);
    iuv = SwitchEdge_IEdge(v);
    */
    if ( u == NO_VERTEX || iuv < 0 )
        return 0; /* should not happen */
    /* check edge adjacency */
    if ( pBNS->edge[iuv].neighbor1  != (u/2-1)  && pBNS->edge[iuv].neighbor1 != v/2-1 ||
         (pBNS->edge[iuv].neighbor12 ^ (u/2-1)) != (v/2-1) ) {
        return 0; /* !!! should not happen !!! */
    }
         
#if ( CHECK_TACN == 1 )
    u_is_taut = ((pBNS->vert[u/2-1].type & pBNS->type_T)  == pBNS->type_T )? TYPE_T :
                ((pBNS->vert[u/2-1].type & pBNS->type_CN) == pBNS->type_CN)? TYPE_CN : 0;
    w_is_taut = ((pBNS->vert[w/2-1].type & pBNS->type_T)  == pBNS->type_T )? TYPE_T :
                ((pBNS->vert[w/2-1].type & pBNS->type_CN) == pBNS->type_CN)? TYPE_CN : 0;
    if ( (u_is_taut | w_is_taut) == (TYPE_T | TYPE_CN ) ) {
        /* rescap must have already been checked */
        return 1;
    }
#endif

    return 0;
#undef TYPE_T
#undef TYPE_CN
}
#else  /* } FIX_TACN_POSSIBLE_BUG { */
/*********************************************************************************/
int bIgnoreVertexNonTACN_group( BN_STRUCT* pBNS, Vertex v, Vertex w, Edge *SwitchEdge )
{
    int    u_is_taut=0, w_is_taut=0;
    Vertex u;
    EdgeIndex iuv;
    if ( v <= 1 || w <= 1 )
        return 0;
#if ( CHECK_TACN == 1 )
    if ( !pBNS->type_TACN ||
         (pBNS->vert[v/2-1].type & pBNS->type_TACN) ) {
        return 0;
    }
    if ( !pBNS->type_T || !pBNS->type_CN )
        return 0; /* should not happen */
#endif
    u = GetPrevVertex( pBNS, v, SwitchEdge, &iuv );
    /*
    u   = SwitchEdge_Vert1(v);
    iuv = SwitchEdge_IEdge(v);
    */
    if ( u == NO_VERTEX || iuv < 0 )
        return 0; /* should not happen */
    /* check edge adjacency */
    if ( pBNS->edge[iuv].neighbor1  != (u/2-1)  && pBNS->edge[iuv].neighbor1 != v/2-1 ||
         (pBNS->edge[iuv].neighbor12 ^ (u/2-1)) != (v/2-1) ) {
        return 0; /* !!! should not happen !!! */
    }
         
#if ( CHECK_TACN == 1 )
    if ( ((u_is_taut  = (pBNS->vert[u/2-1].type & pBNS->type_T)  == pBNS->type_T) ||
          (             (pBNS->vert[u/2-1].type & pBNS->type_CN) == pBNS->type_CN)) &&
         ((w_is_taut  = (pBNS->vert[w/2-1].type & pBNS->type_T)  == pBNS->type_T) ||
          (             (pBNS->vert[w/2-1].type & pBNS->type_CN) == pBNS->type_CN)) &&
         u_is_taut + w_is_taut == 1 ) {
        /* rescap must have already been checked */
        return 1;
    }
#endif

    return 0;
}
#endif /* } FIX_TACN_POSSIBLE_BUG { */

#if( FIX_KEEP_H_ON_NH_ANION == 1 )
/*********************************************************************************/
/* detect an attempt to remove H from -NH(-) to make =N(-);                      */
/* all taut atoma except N are 'acidic'                                          */
/*********************************************************************************/
int bIsRemovedHfromNHaion( BN_STRUCT* pBNS, Vertex u, Vertex v )
{
    int    i, u2, v2, vat2;
    Vertex vtg, vat;
    BNS_VERTEX *pvAT, *pvCN;
    BNS_EDGE   *pEdge;
    if ( !pBNS->type_TACN || u <= 1 || v <= 1 ||
         u%2 || !(v%2) /* the edge flow may only increase */ ) {
        return 0;
    }
    if ((pBNS->vert[u2 = u/2-1].type & pBNS->type_TACN) ||
        (pBNS->vert[v2 = v/2-1].type & pBNS->type_TACN) ) {
        return 0; /* add/remove H is allowed for acidic atoms */
    }
    if ( !pBNS->type_T || !pBNS->type_CN )
        return 0; /* should not happen */
    /* find which of u, v vertices is N and which is t-group */
    if ( ((pBNS->vert[u2].type & pBNS->type_T)  == pBNS->type_T ) && v2 < pBNS->num_atoms ) {
        vtg = u;
        vat = v;
    } else
    if ( ((pBNS->vert[v2].type & pBNS->type_T)  == pBNS->type_T ) && u2 < pBNS->num_atoms ) {
        vtg = v;
        vat = u;
    } else {
        return 0;
    }
    vat2 = vat/2-1;
    pvAT = pBNS->vert + vat2;  /* atom */
    for ( i = pvAT->num_adj_edges-1; 0 <= i; i -- ) {
        pEdge = pBNS->edge + pvAT->iedge[i];
        pvCN  = pBNS->vert + (pEdge->neighbor12 ^ vat2);
        if ( ((pvCN->type & pBNS->type_CN) == pBNS->type_CN) &&  pEdge->flow > 0 ) {
            return 1; /* detected */
        }
    }
    return 0;
}
#endif
#if ( FIX_AVOID_ADP == 1 )
/************************************************************************
   Detect  (tg)-N=A-A=A-A=N-(tg)
                        u v  w
     k =    5   4 3 2 1 0 1  2
            ^
          odd number means ADP
*************************************************************************/
int bIsAggressiveDeprotonation( BN_STRUCT* pBNS, Vertex v, Vertex w, Edge *SwitchEdge )
{
#define TYPE_T   1   /* t-group [also called H-group] */
#define TYPE_CN  2   /* (-)c-group */
#define TYPE_AT  4
    int    k, v2, u2, w2, u2_next, type0, type1, type2, type;
    Vertex u, u_next;
    EdgeIndex iuv;
    if ( v <= 1 || w <= 1 )
        return 0;

    if ( !pBNS->type_TACN || !pBNS->type_T || !pBNS->type_CN )
        return 0; /* should not happen */
    v2 = v/2 - 1;
    w2 = w/2 - 1;
    if ( v2 >= pBNS->num_atoms || w2 < pBNS->num_atoms )
        goto cross_edge;

    if ( !((pBNS->vert[w2].type & pBNS->type_T)  == pBNS->type_T ) &&
         !((pBNS->vert[w2].type & pBNS->type_CN) == pBNS->type_CN) )
        goto cross_edge;
    /* v ia an atom, w is a t-group, v != w' */
    for ( k = 0, u = v; 1 < (u_next = u, u = GetPrevVertex( pBNS, u, SwitchEdge, &iuv )); k ++ ) {
        u2 = u/2 - 1;
        if ( u2 >= pBNS->num_atoms ) {
            /* moving backward along the alt path we have found a vertex
               that is not an atom. Possibly it is a t- or (-)c-group */
            if ( !( k % 2 ) ) {
                return 0; /* even vertex -- always okay */
            }
            if ( !((pBNS->vert[u2].type & pBNS->type_T)  == pBNS->type_T ) &&
                 !((pBNS->vert[u2].type & pBNS->type_CN) == pBNS->type_CN) ) {
                /* not a t- or (-)c-group */
                return 0;
            }
            u2_next = u_next/2 - 1;
            if ( !(pBNS->vert[v2     ].type & pBNS->type_TACN) &&
                 !(pBNS->vert[u2_next].type & pBNS->type_TACN)  ) {
                /* none of the atoms at the ends are N */
                return 0;
            }
            return 1;
        }
    }
    return 0;
cross_edge:
    /*****************************************************************************
     * v and w (v=w') are same vertex reached with opposite "phases".
     * w cannot be (t) because this would have been detected earlier -- ???
     *   (t)-A=A-A=A-A=A-(t)
     *           v 
     *    3  2 1 0 1 2 3  4
     *   kv               kw
     *   (kv + kw)%2 == 1  <==> aggressive deprotonation
     *****************************************************************************/
    if ( v == prim(w) ) {
        type0 = 0;
        if ( v2 >= pBNS->num_atoms ) {
            type0 = ((pBNS->vert[v2].type & pBNS->type_T)  == pBNS->type_T )? TYPE_T :
                    ((pBNS->vert[v2].type & pBNS->type_CN) == pBNS->type_CN)? TYPE_CN : 0;
        }

        

    }
          
    
    return 0;
}
#endif

/*********************************************************************************/
int rescap(  BN_STRUCT* pBNS, Vertex u, Vertex v, EdgeIndex iuv )
{
    BNS_ST_EDGE    *pst_edge;
    BNS_EDGE       *pedge;

    int f;
    S_CHAR s_or_t;
    int    bBackward = GetEdgePointer( pBNS, u, v, iuv, &pedge, &s_or_t );
    if ( !IS_BNS_ERROR(bBackward) ) {
        
        if ( s_or_t ) {
            pst_edge = (BNS_ST_EDGE *)pedge;
            f    = (pst_edge->flow & EDGE_FLOW_ST_MASK);
            if ( !bBackward ) {
                f = (int)pst_edge->cap-f;
            }
        } else {
            f    = (pedge->flow & EDGE_FLOW_MASK);
            if ( !bBackward ) {
                f = (int)pedge->cap-f;
            }
        }
        return f;

    }
    return bBackward; /* error */
}
/*********************************************************************************
    W.Kocay, D.Stone,
    "An Algorithm for Balanced Flows",
    The Journal of Combinatorial Mathematics and Combinatorial Computing,
    vol. 19 (1995) pp. 3-31
    
    W.Kocay, D.Stone,
    "Balanced network flows",
    Bulletin of the Institute of Combinatorics and its Applications,
     vol. 7 (1993), pp. 17--32

    N = a balanced network (bipartite directed graph) of:
        n=2*V+2 vertices (incl. s (source) and t (target);
          each other vertex i is included 2 times:
          set X (x[i]) of V vertices (v) and a set Y (y[j]) of
          V complementary vertices (v') so that x[i]' = y[i], x[i]''=x[i], and
        m=2*E+2*V edges (each original undirected edge i-j is represented as
          2 directed edges: x[i]->y[j] and x[j]->y[i]; plus
          V edges s->x[i] and V edges y[j]->t)
    v'   = a complement vertex to v
    v'u' = (uv)' = a complement edge to uv
    rescap(uv) = cap(uv)-f(uv) if uv is a forward edge
               = f(uv)         if uv is a backward edge
    (i)   0 <= f(uv) <= cap(uv)
    (ii)  f+(u) = f-(u) for all u in X, Y where f+(u) is a total flow out of u,
                  and f-(u) is a total flow into u
    (iii) f(uv) = f((uv)') (balanced flow condition)

    S     = a set of all s-searchable vertices
    S-    = all other vertices
    if t in S, then N contains a valid augmenting path P, so that flow can be
            augmented on both P and P'
    if t not in S, then let K=[S,S-], the set of all edges directed from S to S-.
            K is an edge-cut that has a special structure.
    Let
    A  = {x[i], y[i] |x[i] not in S, y[i]     in S}
    B  = {x[i], y[i] |x[i]     in S, y[i] not in S}
    C  = {x[i], y[i] |x[i]     in S, y[i]     in S}
    D  = {x[i], y[i] |x[i] not in S, y[i] not in S}
    N[C]  = subgraph of N induced by C;
            it consists of of a number of connected components C[i]
    K[i]  = those edges of K with one endpoint in C[i].

    If t is in S- then

    i)   Each C[i] = C[i]'
    ii)  There are no edges between C and D
    iii) Each K[i] has odd capacity, it is called a balanced edge-cut.

    balcap(K) = cap(K) - odd(K), where odd(K) is the number of connected components in N[C].
    Name "odd(K)" is because each cap(K[i]) is odd.

    Max-Balanced-Flow-Min-Balanced-Cut Theorem:

    Let f be a balanced flow in N, and let K be any balanced edge-cut.
    The value of a maximum balanced flow equals the capacity of a minimum
    edge-cut, that is, val(f) = balcap(K) when f is maximum and K is minimum.

*********************************************************************************/

/*******************************************************/
/*                                                     */
/*        VERTEX NUMBERING                             */
/*                                                     */
/*   Total number of atoms    = n                      */
/*   Total number of vertices = 2*n+2                  */
/*******************************************************/
/*                                                     */
/* atom numbering starts from 0:                       */
/*                                                     */
/* atoms    s t x0 y0 x1 y1 ... xi   yi   ...xn   yn   */
/* vertices 0 1 2  3  4  5  ... 2i-2 2i-1 ...2n-2 2n-1 */
/*                                                     */
/* atom = vertex/2-1; if negative then s or t          */
/*                                                     */
/* vertex = (atom + 1) * 2 + i; i=0 for x, i=1 for y   */
/*                                                     */
/*******************************************************/

/*********************************************************************************/
/* v' variable is called prim(v) for now */
int BalancedNetworkSearch ( BN_STRUCT* pBNS, BN_DATA *pBD, int bChangeFlow )
{   /* N has source s, target t, the mirror network M is constructed */
    /* the tree T contains a valid sv-path for each v in T. Simultaneously the complementary
       tree T' is built as indicated in comments. The trees T and T' must have no edges or
       vertices in common. Initially T will be built as in breadth-first-search, and T' will
       be the complementary tree, it will contain the complementary valid v't-path.
    */

    Vertex          *BasePtr =       pBD->BasePtr;
    Edge            *SwitchEdge =    pBD->SwitchEdge;
    S_CHAR          *Tree =          pBD->Tree;
    Vertex          *ScanQ =         pBD->ScanQ;     
    int              QSize =         pBD->QSize;
    Vertex          *Pu =            pBD->Pu;
    Vertex          *Pv =            pBD->Pv;
    int              max_len_Pu_Pv=  pBD->max_len_Pu_Pv;

    /* added to translate into C */
    int i, k, degree, delta, ret = 0;
    Vertex u, b_u, v, b_v, w;
    EdgeIndex iuv;
#if( BNS_RAD_SEARCH == 1 )
    int              n, bRadSearch   = (BNS_EF_RAD_SRCH & bChangeFlow) && pBD->RadEndpoints;
    BRS_MODE         bRadSrchMode    = RAD_SRCH_NORM;
    int              bRadSearchPrelim = 0;
    if ( bRadSearch ) {
        pBD->nNumRadEndpoints = 0;
        bRadSrchMode          = pBD->bRadSrchMode;
        bRadSearchPrelim      = pBNS->type_TACN && bRadSrchMode == RAD_SRCH_NORM;
    }
#endif

/*  -- Always --
    Vertex_s = FIRST_INDX;
    Vertex_t = Vertex_s+1;
*/
    QSize = k = 0;     /* put s on ScanQ = set S */
    ScanQ[QSize] = Vertex_s;
    BasePtr[Vertex_t] = Vertex_s;
    BasePtr[Vertex_s] = BLOSSOM_BASE; /* create initial blossom C(Vertex_s) with base s */
    Tree[Vertex_s]    = TREE_IN_1;

    do {
        u = ScanQ[k]; /* select u from the head of ScanQ */
        /* since u is on the queue, it has a blossom C(U) with base b_u */
        b_u = FindBase( u, BasePtr );
        degree = GetVertexDegree( pBNS, u );
#if( BNS_RAD_SEARCH == 1 )
        n = 0;
#endif
        for ( i = 0; i < degree; i ++ ) {
            /* v = vert[u].neighbor[i]; */
            v = GetVertexNeighbor( pBNS, u, i, &iuv );
            if ( v == NO_VERTEX ) {
                continue; /* the atom has only single bonds, ignore it */
            }
#if( BNS_RAD_SEARCH == 1 )
            if ( !k && bRadSrchMode == RAD_SRCH_FROM_FICT && v/2 <= pBNS->num_atoms ) {
                continue; /* start from fict. vertices only */
            }
            if ( bRadSearchPrelim && v/2 > pBNS->num_atoms ) {
                continue; /* during initial add/remove H allow radical movement only through real atoms */
            }
#endif
            if ( /* PrevPt[u] != v ** avoid edges of T */
                 ( SwitchEdge_Vert1(u) != v || SwitchEdge_Vert2(u) != u )  /* avoid edges of T */
                 && (ret = rescap(pBNS, u, v, iuv)) > 0 ) {
                /* special treatment to prevent H<->(-) replacement on non-acidic atoms */
                /*----------------------------------------------------------------------*/
                if ( pBNS->type_TACN ) {
                    if ( bIgnoreVertexNonTACN_atom( pBNS, u, v ) ) {
                        continue;
                    }
                    if ( bIgnoreVertexNonTACN_group( pBNS, u, v, SwitchEdge ) ) {
                        continue;
                    }
#if( FIX_KEEP_H_ON_NH_ANION == 1 )
                    if ( bIsRemovedHfromNHaion( pBNS, u, v ) ) {
                        continue;
                    }
#endif
#if ( FIX_AVOID_ADP == 1 )
                    if ( bIsAggressiveDeprotonation( pBNS, u, v, SwitchEdge ) ) {
                        continue;
                    }
#endif

                }
                /*----------------------------------------------------------------------*/
                b_v = FindBase(v, BasePtr); /* Notation: b_x is a base of x */
                
                if ( b_v == NO_VERTEX ) {  /* originally 0 instead of NO_VERTEX */
                    /* Important note: following "A. Note of Implementing the Algorithm" from the
                       article by Kocay and Stone, all references to PrevPt[a] have been
                       replaced with SwitchEdge[a][0]=SwitchEdge_Vert1(a); to be on a safe side
                       the check whether (SwitchEdge_Vert2(a)==a) has been added.
                    */

                    /* v is not yet in T or T' -- add it to T and M */
                    /*PrevPt[v] = u; ** this effectively adds v to T and v' to T' */
                    QSize ++;
                    ScanQ[QSize]  = v;
                    TREE_MARK(v, TREE_IN_1); /* mark v s-reachable (in T) */
                    TREE_MARK(prim(v), TREE_IN_2); /* mark v' in T' */

                    /* Switch Edges: If u in T then Pu (a valid su-path) can be constructed
                       by successfully executing u=PrevPt[u] until u=s.
                       For vertices in T' the situation is different.
                       Suppose uv and v'u' are added to a mirror network M creating a blossom:

                                    s            T    (Note: in the code v' is prim(v),
                     T             / \                       u' is prim(u), etc.)
                                  /   ...
                                 w=x1
                                / \          u in T,   v in T'
                            u=y2   v'=y3
                                \ /      <--- two added edges uv and (uv)'=v'u'
                     --------    X   ------------intersection of the edges ------------------------
                                / \      <--- (the edges intersection in the picture is shown as X)
                           u'=x2   v=x3      u' it T', v' in T
                                \ /
                    T'           w'=y1
                                  \   ...
                                   \ /
                                    t           T'

                       Vertices v and u' now become s-reachable;
                       The valid paths to v and u' must use the edges uv and v'u' respectively.

                       For each vertex z in S we define a switch-edge that allows a valid sz-path
                       to be constructed, SwitchEdge[v]=uv and SwitchEdge[u']=v'u'. (We don't
                       know the direction of the edge uv, it may be (u,v) or (v,u). In either case,
                       the complementary edge is indicated by v'u').

                       Vertex w' also becomes s-reachable when uv is added to M, and a valid sw'-path
                       must use one of uv and v'u'. Therefore we choose one of them, say uv (see below
                       the rule of choosing the switch-edge), and define SwitchEdge[w'] = uv.

                       When the addition of an edge uv to M causes a vertex z to become s-reachable
                       (where z was previously non-reachable), z is placed on the ScanQ, that is, into S.
                       The edge uv is said to be a switch-edge for z.

                       Rule: We choose the the order of the vertices uv to be such that the valid sz-path
                       consists of a valid su-path, followed by edge uv, followed by a valid vz-path.

                       For vertices z in T we can take SwitchEdge[z]=yz where y=PrevPt[z] since
                       it is the edge yz that allows z to be s-reachable.
                       For vertices z not in S we take SwitchEdge[z]=NONE.

                    */

                    /* v is not yet in T or T' -- add it to T and M */
                    SwitchEdge_Vert1(v) = u; /* this effectively adds uv and v'u' to M */
                    SwitchEdge_IEdge(v) = iuv;

                    BasePtr[prim(v)] = v;
                    BasePtr[v]  = BLOSSOM_BASE; /* create a trivial blossom C(v) with base v */
#if( BNS_RAD_SEARCH == 1 )
                    n ++;
#endif
                } else
                if ( TREE_IS_S_REACHABLE(prim(v)) /*Is_s_Reachable(prim(v)*/
                     /* if v' is reachable, an st-path is given by P(u)-uv-P'(v') */
                     /*&& PrevPt[prim(u)] != prim(v) ** avoid edges of T' */
                     && (SwitchEdge_Vert1(prim(u)) != prim(v) || SwitchEdge_Vert2(prim(u)) != prim(u)) /* avoid edges of T' */
                     && b_u != b_v
                     && !(pBNS->type_TACN && bIgnoreVertexNonTACN_group( pBNS, prim(v), u, SwitchEdge ))
#if( FIX_KEEP_H_ON_NH_ANION == 1 )
                     && !(pBNS->type_TACN && bIsRemovedHfromNHaion( pBNS, prim(v), u ))
#endif
                     ) {
#if( BNS_RAD_SEARCH == 1 )
                    n ++;
#endif
                     /* there is now a valid sv-path via u avoiding b_v (unless v==b_v)
                       => u, v, u', and v' now all become part of the same connected component of M[C] */
                    w = MakeBlossom( pBNS, ScanQ, &QSize, Pu, Pv, max_len_Pu_Pv, SwitchEdge, BasePtr, u, v, iuv, b_u, b_v, Tree );
                    /* this constructed the new blossom and returned its base */
                    if ( IS_BNS_ERROR( w ) ) {
                        pBD->QSize = QSize;
                        return w; /* error */
                    }
                    b_u = w; /* the new base of C(u) */
                    if ( prim(w) == Vertex_t ) {
                        /* t is now s-reachable, a valid augmenting path P exists in M */
                        delta = FindPathCap( pBNS, SwitchEdge, Vertex_s, Vertex_t, 10000 ); /* compute the residual capacity of P + P' */
                        if ( IS_BNS_ERROR( delta ) ) {
                            pBD->QSize = QSize;
                            return delta; /* error */
                        }
#if( ALLOW_ONLY_SIMPLE_ALT_PATH == 1 )
                        if ( pBNS->bNotASimplePath || abs(delta) > 1 ) {
                            delta = 0;
                        }
#endif
                        if ( delta ) {
                            pBNS->bChangeFlow |= (bChangeFlow & BNS_EF_CHNG_FLOW);
                        }
                        ret = PullFlow( pBNS, SwitchEdge, Vertex_s, Vertex_t, delta, 0, bChangeFlow ); /* augment on a pair of valid st-paths */
                        pBD->QSize = QSize;
                        return ( IS_BNS_ERROR(ret)? ret : delta );
                    }
                }
            } else
            if ( IS_BNS_ERROR( ret ) ) {
                pBD->QSize = QSize;
                return ret; /* error */
            }
        }
#if( BNS_RAD_SEARCH == 1 )
        if ( bRadSearch && !n ) {
            /* the BNS stopped at u */
            n = RegisterRadEndpoint( pBNS, pBD, u);
            if ( IS_BNS_ERROR( n ) ) {
                pBD->QSize = QSize;
                return n;
            }
        }
#endif
        k ++; /* advance ScanQ */
    } while( k <= QSize );
    /* if this point is reached, no valid augmenting path exists, ScanQ contains
       the set S of all s-reachable vertices and K=[S,S-] is a minimum balanced edge-cut */
    /* ClearFlowMarks( vert, num_vert); */
    pBD->QSize = QSize;
    return 0;
}
/********************************************************************
Blossoms.

  The vertices of a mirror network M consist T U T'. Intersection T ^ T' is empty.

  The edges of M consist of switch-edges and their complements because edges
  are added in complementary pairs, one of which is always a switch-edge.

  The base of every blossom is in T.
  Let C(i) be a blossom with base b_i. Since C(i)=C(i)', C(i) contains vertices of T and T'.
  Since every valid sv-path to v in C(i) contains b_i, b_i is the first s-reachable vertex of C(i).

  Suppose the mirror network M contains a valid sz-path P(z) to all vertices z in ScanQ.
  Every vertex of P(z) is s-reachable therefore its vertices are all in blossoms or
  trivial blossoms.

  Let z be an s-reachable vertex and P(z) be a valid path in M.
  Then every valid sz-path in M contains exactly the same sequence of blossom bases as P(z).

*********************************************************************/

/***********************************************************************
    BasePtr[u] = -2  NO_VERTEX       u is not a blossom
                 -1  BLOSSOM_BASE    u is the base of its blossom
                  v                  a vertex closer to the base
************************************************************************/
Vertex FindBase( Vertex u, Vertex *BasePtr )
{
    if ( BasePtr[u] == NO_VERTEX ) {
        return NO_VERTEX;
    } else
    if ( BasePtr[u] == BLOSSOM_BASE ) {
        return u;
    } else {
        Vertex b;
        b = FindBase(BasePtr[u], BasePtr );
        BasePtr[u] = b; /* path compression */
        return b;
    }
}

/*********************************************************************/
/* returns index of the last path element and the path               */
/*********************************************************************/
int FindPathToVertex_s( Vertex x, Edge *SwitchEdge, Vertex *BasePtr, Vertex *Path, int MaxPathLen )
{
    /* x is the base of a blossom, construct a valid Path of blossom bases to s */
    int i = 0;
    Path[i] = x;
    while ( x != Vertex_s ) {
        x = FindBase(SwitchEdge_Vert1(x), BasePtr);
        if ( ++i < MaxPathLen ) {
            Path[i] = x;
        } else {
            return BNS_WRONG_PARMS;
        }
    }
    return i;
}
/*********************************************************************/
/* make a blossom                                                    */
/*********************************************************************/
Vertex MakeBlossom( BN_STRUCT* pBNS, Vertex *ScanQ, int *pQSize,
                    Vertex *Pu, Vertex *Pv, int max_len_Pu_Pv,
                    Edge *SwitchEdge, Vertex *BasePtr,
                    Vertex u, Vertex v, EdgeIndex iuv, Vertex b_u, Vertex b_v, S_CHAR *Tree )
{
    /* In order to find the base of the new blossom, the paths
       P(u) and P(v') are constructed and compared in order to
       find the last blossom base they have in common which
       is reachable on a valid path.

       Edge uv connects two blossoms, their bases are b_u and b_v.
    */
    Vertex w, z;
    int len_Pu, len_Pv;
    int i, j;
    EdgeIndex izw;

    len_Pu = FindPathToVertex_s( b_u, SwitchEdge, BasePtr, Pu, max_len_Pu_Pv );
    if ( IS_BNS_ERROR( len_Pu ) ) {
        return len_Pu;
    }
    len_Pv = FindPathToVertex_s( b_v, SwitchEdge, BasePtr, Pv, max_len_Pu_Pv );
    if ( IS_BNS_ERROR( len_Pv ) ) {
        return len_Pv;
    }
    i = len_Pu;
    j = len_Pv;
    /* initially Pu[i] and Pv[j] both equal to s, but their first elements are different */
    /* find the last blossom base common to Pu and Pv */
    while ( i >= 0 && j >= 0 && Pu[i] == Pv[j]  ) {
        /* was (Pu[i]==Pv[j] && i>=0 && j>=0) => tested Pu[-1], Pv[-1] <- pointed by W.Ihlenfeldt 08-26-2004*/
        i --;
        j --;
    }
    i ++;
    w   = Pu[i]; /* w is the last common vertex */
    z   = SwitchEdge_Vert1(w);
    izw = SwitchEdge_IEdge(w);
    /* now extend the blossom if rescap(zw) >= 2 */
    while ( w != Vertex_s && rescap(pBNS, z, w, izw) >= 2 )
    {
        i++;
        w = Pu[i];
        z = SwitchEdge_Vert1(w);
        izw = SwitchEdge_IEdge(w);
    }
    /* w is the base of the new blossom */
    /* first follow the path Pu from w to b_u */
    for ( i = i-1; i >= 0; i -- ) {
        z = Pu[i];  /* z is the base of the blossom */
        BasePtr[z] = w;
        BasePtr[prim(z)] = w; /* w is the new base of the blossom */
        /* z and z' may already be part of a blossom that is being
           swallowed into a larger blossom.
           We don't want to change the switch edge in that case.
        */
        
        if ( !TREE_IS_ON_SCANQ(prim(z)) /*!IsInScanQ(prim(z)) */) 
        {
            SwitchEdge_Vert1(prim(z)) = prim(v);     /* set the switch edge of z' */
            /* SwitchEdge[prim(z)][1] = prim(u);  */
            SwitchEdge_IEdge(prim(z)) = iuv;
            (*pQSize) ++;
            ScanQ[*pQSize] = prim(z);    /* add z' to ScanQ */
            TREE_MARK(prim(z), TREE_IN_2BLOSS); /* mark z' s-reachable */
        }
    }
    /* now follow the path Pv */
    for ( j = j; j >= 0; j -- ) {
        z = Pv[j]; /* z is the base of the blossom */
        BasePtr[z] = w;
        BasePtr[prim(z)] = w; /* w is the new base of the blossom */
        /* z and z' may already be part of a blossom that is being
           swallowed into a larger blossom.
           We don't want to change the switch edge in that case.
        */

        if ( !TREE_IS_ON_SCANQ(prim(z)) /*!IsInScanQ(prim(z)) */ )
        {
            SwitchEdge_Vert1(prim(z)) = u;     /* set the switch edge of z' */
            /* SwitchEdge[prim(z)][1] = v;  */
            SwitchEdge_IEdge(prim(z)) = iuv;
            (*pQSize) ++;
            ScanQ[*pQSize] = prim(z);    /* add z' to ScanQ */
            TREE_MARK(prim(z), TREE_IN_2BLOSS); /* mark z' s-reachable */
        }
    }

    if ( !TREE_IS_ON_SCANQ(prim(w)) /* !IsInScanQ(prim(w))*/ ) 
    {   /* add w' to the blossom */
        SwitchEdge_Vert1(prim(w)) = u;     /* set the switch edge of w' */
        /* SwitchEdge[prim(w)][1] = v;  */
        SwitchEdge_IEdge(prim(w)) = iuv;
        (*pQSize) ++;
        ScanQ[*pQSize] = prim(w); /* add w' to ScanQ */
        TREE_MARK(prim(w), TREE_IN_2BLOSS);  /* mark w' s-reachable */
    }
    return w;
}
/*****************************************************************************
    When t is found to be s-reachable, a valid st-path P is known to exist.
    Its complementary path P' is also valid. Once the residual capacity
    delta(P) is known, the flow is augmented by calling PullFlow(s,t,delta).
    It constructs the path P by using the switch-edges.
    Let uv=SwitchEdge[t].
    Then P is given by a valid su-path, followed by the edge uv, followed by
    a valid vt-path.
    PullFlow is a recursive procedure that constructs the path and its complement.

    Let wz=SwitchEdge[y]. PullFlow(x, y, delta) uses the xw- and zy-portions of P
    (see below). Since it must also augment on P' simultaneously, the zy-portion
    is replaced by the y'z'-portion.

         x                  y'
         |                  |            P:   x--w--z--y
       P |                  |  P'        P':  y'-z'-w'-x'
         |                  o
         o                   \           
        /   w'          z     \          
       /   o----\   /----o    /
       \  /      \ /      \  /
        \/        X        \/
        /\       / \       /\
       /  \ w   /   \  z' /  \
       \   o----     ----o   /           Using a switch-edge wz and w'z'
        \                   /            to construct P and P'          
         o                 o  
         |                 | 
         |                 |
         x'                y
    
      
 ******************************************************************************/

int PullFlow( BN_STRUCT *pBNS, Edge *SwitchEdge, Vertex x, Vertex y, int delta, S_CHAR bReverse, int bChangeFlow )
{ /* 
     Augment the flow by delta on all edges on a path P
     between x and y in the order of the path;
     AugmentEdge( pBNS, w, z, iwz, delta, 0 ) means the path is in w->z direction
     AugmentEdge( pBNS, w, z, iwz, delta, 1 ) means the path is in w<-z direction

     Unlike PullFlow in the paper by Kocay & Stone, here the augmentation
     always starts at "s", proceeds sequentially through the path end terminates at "t".
     Since we do not really need the complement path, PullFlow ignores it.

  */
  
    Vertex w, z;
    EdgeIndex iwz;
    int ret = 0;

    w   = SwitchEdge_Vert1(y);
    z   = SwitchEdge_Vert2(y);
    iwz = SwitchEdge_IEdge(y);
    if ( bReverse ) {
        /* P consists of a path from x to w, then wz, then a path from z to y.  */
        /* z may equal y, in which case z is just PrevPt[y] */
        if ( z != y ) {
            ret = PullFlow( pBNS, SwitchEdge, prim(y), prim(z), delta, (S_CHAR)(1-bReverse), bChangeFlow ); /* augment between z and y */
        }
        if ( !IS_BNS_ERROR(ret) ) {
            ret = AugmentEdge( pBNS, w, z, iwz, delta, bReverse, bChangeFlow);
        }
        /* Do not augment the complementary path: AugmentEdge( prim(z), prim(w), vert, delta); */
        /* w may equal x, in which case there is no need to call PullFlow(x, w) */
        if ( w != x && !IS_BNS_ERROR(ret) ) {
            ret = PullFlow( pBNS, SwitchEdge, x, w, delta, bReverse, bChangeFlow ); /* augment between x and w */
        }
    } else {
        /* P consists of a path from x to w, then wz, then a path from z to y.  */
        /* w may equal x, in which case there is no need to call PullFlow(x, w) */
        if ( w != x && !IS_BNS_ERROR(ret) ) {
            ret = PullFlow( pBNS, SwitchEdge, x, w, delta, bReverse, bChangeFlow ); /* augment between x and w */
        }
        if ( !IS_BNS_ERROR(ret) ) {
            ret = AugmentEdge( pBNS, w, z, iwz, delta, bReverse, bChangeFlow);
        }
        /* z may equal y, in which case z is just PrevPt[y] */
        if ( z != y && !IS_BNS_ERROR(ret) ) {
            ret = PullFlow( pBNS, SwitchEdge, prim(y), prim(z), delta, (S_CHAR)(1-bReverse), bChangeFlow ); /* augment between z and y */
        }
    }
    return ret;
}

/********************************************************************************
Before augmenting on the two paths, it is necessary to find delta(P).
This can be done by following the paths and computing the minimum
residual capacity of all edges on P. An edge on both P and P' counts
for only half of its actual residual capacity, since augmentng on P by
delta will simutaneously reduce its capacity on P' by delta.
The path P can only be followed by using the switch-edges, as in PullFlow(...).
FindPathCap( x, y, delta ) is a recursive procedure that finds the residual
capacity on the portion of P between x and y. delta is the minimum capacity
found so far along the path.
********************************************************************************/
int FindPathCap( BN_STRUCT* pBNS, Edge *SwitchEdge, Vertex x, Vertex y, int delta )
{ /* find the minimum residual capacity of all edges
     between x and y in a valid st-path P.
     delta is the minimum found so far
     the vertices occur in order s,...,x,...,y,...,t along P
     the vertices occur in order s,...,y',...,x',...,t along P'
  */
    Vertex w, z, iwz;
    int    cap, delta2;
    static int level;

    if ( level ++ > 50 ) {
#ifdef _DEBUG
        int stop = 1;
#else
    ;
#endif
    }


    w   = SwitchEdge_Vert1(y);
    z   = SwitchEdge_Vert2(y); /* wz is on the path P */
    iwz = SwitchEdge_IEdge(y); /* edge index */

    /* rescap_mark() detects edges passed 2 times and reduces rescap */
    cap = rescap_mark( pBNS, w, z, iwz );

    if ( IS_BNS_ERROR( cap ) ) {
        level --;
        return cap;
    }
    if ( cap < delta ) {
        delta = cap;
    }
    /* P consists of a path from x to w, then wz, then a path from z to y */
    if ( w != x ) {
        delta2 = FindPathCap( pBNS, SwitchEdge, x, w, delta );
        delta = inchi_min( delta2, delta );
    }
    if ( z != y ) {
        delta2 = FindPathCap( pBNS, SwitchEdge, prim(y), prim(z), delta );
        delta = inchi_min( delta2, delta );
    }
    level --;
    return delta;
}


/* BT = bond types */
#define BT_ALTERN_BOND           1      /* 1-2, possibly stereo */
#define BT_OTHER_ALTERN_BOND     2      /* 1-3, 2-3, 1-2-3 alternating non-stereo non-taut bonds */

#define BT_ALTERN_NS_BOND        4

#define BT_TAUTOM_BOND           8

#define BT_ALTERN_UNKN_BOND     16

#define BT_IGNORE_BOND           0

#define BT_NONSTEREO_MASK        (BT_TAUTOM_BOND|BT_ALTERN_NS_BOND)

#define BT_ALT_BOND_MASK         (BT_ALTERN_BOND|BT_OTHER_ALTERN_BOND)

#define BT_NONTAUT_BOND_MASK     (BT_ALTERN_BOND|BT_OTHER_ALTERN_BOND|BT_ALTERN_NS_BOND)

/* BNS members redefinitions for finding non-stereo bonds */
/* BNS_EDGE */
#define nBlockNumberAltBns   flow   /* internal variable of the DFS traversal: mark traversed bonds */
#define nNumAtInBlockAltBns  cap
#define nBondTypeInpAltBns   pass    /* 0=>cannot be stereo at all, 1=>alt or taut non-stereo, 2=>can be stereo */
#define nBondNonStereoAltBns cap     /* 1=>found to be non-stereogenic although BondTypeInp=2; 0=>as in BondTypeInp */

#if( BNS_MARK_ONLY_BLOCKS == 1 )   /* { */
/* BNS_VERTEX */
#define bCutVertexAltBns         st_edge.cap0  /* cut-vertex flag */
#define nRingSystemAltBns        st_edge.cap   /* ordering number of a ring system */
#define nNumAtInRingSystemAltBns st_edge.flow0 /* number of vertices in a ring system */
#define nBlockSystemAltBns       st_edge.flow  /* result of the DFS traversal: even cirquit must be within one block */

#endif  /* } */

#define valenceAltBns            num_adj_edges



/********************************************************************************/
int MarkRingSystemsAltBns( BN_STRUCT* pBNS, int bUnknAltAsNoStereo )
{
    AT_NUMB   *nStackAtom = NULL;
    int        nTopStackAtom;
    AT_NUMB   *nRingStack = NULL;
    int        nTopRingStack; /* was AT_NUMB */
    AT_NUMB   *nBondStack = NULL;
    int        nTopBondStack;
    AT_NUMB   *nDfsNumber = NULL;
    AT_NUMB   *nLowNumber = NULL;
    S_CHAR    *cNeighNumb = NULL;
    AT_NUMB    nDfs;
    AT_NUMB    nNumAtInRingSystem;
    int        i, j, u, w, start, nNumRingSystems, nNumStartChildren;
    BNS_VERTEX *at       = pBNS->vert;
    BNS_EDGE   *bond     = pBNS->edge;
    int        num_atoms = pBNS->num_atoms;
    int        num_edges = pBNS->num_bonds;

    /*  allocate arrays */
    nStackAtom = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nStackAtom[0]));
    nRingStack = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nRingStack[0]));
    nDfsNumber = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nDfsNumber[0]));
    nLowNumber = (AT_NUMB *)inchi_malloc(num_atoms*sizeof(nLowNumber[0]));
    nBondStack = (AT_NUMB *)(num_edges? inchi_malloc(num_edges*sizeof(nBondStack[0])):NULL); /* special case: no bonds 2006-03 */
    cNeighNumb = (S_CHAR  *)inchi_malloc(num_atoms*sizeof(cNeighNumb[0]));
    /*  check allocation */
    if ( !nStackAtom || !nRingStack || !nDfsNumber || !nLowNumber || !nBondStack && num_edges || !cNeighNumb 
        ) {
        nNumRingSystems = CT_OUT_OF_RAM;  /*  program error */ /*   <BRKPT> */
        goto exit_function;
    }

    /********************************************
     *
     * Find Cut-vertices & Blocks
     *
     *             1\      /5   has 3 blocks (maximal subgraphs that
     *  Example:   | >3--4< |   are nonseparable by deleting a single vertex):
     *             2/      \6   (1,2,3, has 3 bonds), (3,4, has 1 bond), and (4,5,6, has 3 bonds)
     *
     *                          Cut-vertices or articulation points are 
     *                          intersections of the blocks: points 3 and 4.
     ********************************************/

    /********************************************************

       RingSystemAlt are atoms connected by alternating bonds
       (as must be indicated in bIsAltBond()):

       BOND_ALTERN
       BOND_ALT_123
       BOND_ALT_13 
       BOND_ALT_23
       
      Since other bonds may be present, we possibly need
      to restart to move to another component 
    *********************************************************/

    nNumRingSystems = 0;
    memset( nDfsNumber, 0, num_atoms*sizeof(nDfsNumber[0]));

    for ( start = 0; start < num_atoms; start ++ ) {
        if ( nDfsNumber[start] )
            continue;
        for ( i = 0; i < at[start].valenceAltBns; i ++ ) {
            if ( bond[at[start].iedge[i]].nBondTypeInpAltBns & BT_ALTERN_BOND )
                goto found_alt;
        }
        continue;

found_alt:


        /*  initiation */
        u               = start; /*  start atom */
        nDfs            = 0;
        nTopStackAtom   =-1;
        nTopRingStack   =-1;
        nTopBondStack   =-1;
        memset( cNeighNumb, 0, num_atoms*sizeof(cNeighNumb[0]));
        /*  push the start atom on the stack */
        nLowNumber[u] = nDfsNumber[u] = ++nDfs;
        nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
        nRingStack[++nTopRingStack] = (AT_NUMB)u;

        nNumStartChildren = 0;

        do {
            /* advance */
            /*while ( (int)at[i=nStackAtom[nTopStackAtom]].valenceAltBns > (j = (int)cNeighNumb[i]) )*/
            /* replaced due to missing sequence point */
            while ( i=(int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i], (int)at[i].valenceAltBns > j )
            {
                cNeighNumb[i] ++;
                if ( !(bond[w=at[i].iedge[j]].nBondTypeInpAltBns & BT_ALT_BOND_MASK) ) {
                    continue;
                }
                u = (int)(bond[at[i].iedge[j]].neighbor12 ^ i);
                if ( !nDfsNumber[u] ) {
                    /* tree edge, 1st visit -- advance */
                    nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
                    nRingStack[++nTopRingStack] = (AT_NUMB)u;
                    nBondStack[++nTopBondStack] = (AT_NUMB)w;
                    nLowNumber[u] = nDfsNumber[u] = ++nDfs;
                    nNumStartChildren += (i == start);
                } else
                if ( !nTopStackAtom || u != (int)nStackAtom[nTopStackAtom-1] ) { /*  may comment out ? */
                    /* back edge: u is not a predecessor of i */
                    if ( nDfsNumber[u] < nDfsNumber[i] ) {
                        /* Back edge, 1st visit: u is ancestor of i. Save and compare */
                        nBondStack[++nTopBondStack] = (AT_NUMB)w;
                        if ( nLowNumber[i] > nDfsNumber[u] ) {
                            nLowNumber[i] = nDfsNumber[u];
                        }
                    }
                }
            }
            cNeighNumb[i] = 0;

            /* back up */
            if ( i != start ) {
                u = (int)nStackAtom[nTopStackAtom-1]; /* predecessor of i */
                if ( nLowNumber[i] >= nDfsNumber[u] ) {
                    /* output the block; the block was entered through its first bond u->i */
                    nNumRingSystems ++;
                    /*at[u].nBlockSystemAltBns = nNumRingSystems;*/ /* mark the atom */
                    nNumAtInRingSystem = 1;
                    /*
                    if ( u != start || nNumStartChildren > 1 ) {
                        at[u].bCutVertexAltBns += 1;  // mark cut-vertex (articulation point)
                    }
                    */
                    while ( nTopRingStack >= 0 ) {
                        j = nRingStack[nTopRingStack--];
                        /*at[j].nBlockSystemAltBns = nNumRingSystems;*/ /*  mark the atom */
                        nNumAtInRingSystem ++;
                        if ( i == j ) {
                            break;
                        }
                    }
                    while ( nTopBondStack >= 0 ) {
                        w = nBondStack[nTopBondStack--];
                        bond[w].nBlockNumberAltBns  = nNumRingSystems; /*  mark the bond */
                        bond[w].nNumAtInBlockAltBns = nNumAtInRingSystem;
                        if ( i == bond[w].neighbor1 && u == (i ^ bond[w].neighbor12) ||
                             u == bond[w].neighbor1 && i == (u ^ bond[w].neighbor12)) {
                            break;
                        }
                    }
                } else
                if ( nLowNumber[u] > nLowNumber[i] ) {
                    /* inherit */
                    nLowNumber[u] = nLowNumber[i];
                }
            }
        } while ( --nTopStackAtom >= 0 );
    }

#if( BNS_MARK_ONLY_BLOCKS != 1 )  /* { */

    /********************************************
     *
     * Find Ring Systems
     * Including chain atoms X: A-X-B, where the bonds (of any kind) are bridges.
     *
     ********************************************/

    /*  initiation */

    nNumRingSystems = 0;

    for ( start = 0; start < num_atoms; start ++ ) {
        if ( at[start].nRingSystemAltBns )
            continue;
        for ( i = 0; i < at[start].valenceAltBns; i ++ ) {
            if ( bond[at[start].iedge[i]].nBondTypeInpAltBns & BT_ALT_BOND_MASK )
                goto found_alt2;
        }
        continue;

found_alt2:

        u               = start; /*  start atom */
        nDfs            = 0;
        nTopStackAtom   =-1;
        nTopRingStack   =-1;
        memset( nDfsNumber, 0, num_atoms*sizeof(nDfsNumber[0]));
        memset( cNeighNumb, 0, num_atoms*sizeof(cNeighNumb[0]));
        /*  push the start atom on the stack */
        nLowNumber[u] = nDfsNumber[u] = ++nDfs;
        nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
        nRingStack[++nTopRingStack] = (AT_NUMB)u;

        do {
            /* advance */
advance_ring:
            /*if ( (int)at[i=nStackAtom[nTopStackAtom]].valenceAltBns > (j = (int)cNeighNumb[i]) )*/
            /* replaced due to missing sequence point */
            if ( i=(int)nStackAtom[nTopStackAtom], j = (int)cNeighNumb[i], (int)at[i].valenceAltBns > j )
            {
                cNeighNumb[i] ++;
                if ( !(bond[at[i].iedge[j]].nBondTypeInpAltBns & BT_ALTERN_BOND) ) {
                    goto advance_ring;
                }
                u = (int)(bond[at[i].iedge[j]].neighbor12 ^ i);
                if ( !nDfsNumber[u] ) {
                    /* tree edge, 1st visit -- advance */
                    nStackAtom[++nTopStackAtom] = (AT_NUMB)u;
                    nRingStack[++nTopRingStack] = (AT_NUMB)u;
                    nLowNumber[u] = nDfsNumber[u] = ++nDfs;
                } else
                if ( !nTopStackAtom || u != (int)nStackAtom[nTopStackAtom-1] ) {
                    /* back edge: u is not a predecessor of i */
                    if ( nDfsNumber[u] < nDfsNumber[i] ) {
                        /* Back edge, 1st visit: u is ancestor of i. Compare */
                        if ( nLowNumber[i] > nDfsNumber[u] ) {
                            nLowNumber[i] = nDfsNumber[u];
                        }
                    }
                }
                goto advance_ring;
            } else {
                cNeighNumb[i] = 0;
            }

            /* back up */
            if ( nDfsNumber[i] == nLowNumber[i] ) {
                /*  found a ring system */
                nNumRingSystems ++;

                /*  unwind nRingStack[] down to i */

                /*  count atoms in a ring system */
                for ( nNumAtInRingSystem = 0, j =  nTopRingStack; 0 <= j; j -- ) {
                    nNumAtInRingSystem ++;
                    if ( i == (int)nRingStack[j] ) {
                        break;
                    }
                }
                while ( nTopRingStack >= 0 ) {
                    j = (int)nRingStack[nTopRingStack--];
                    at[j].nRingSystemAltBns        = (AT_NUMB)nNumRingSystems; /*  ring system id */
                    at[j].nNumAtInRingSystemAltBns = nNumAtInRingSystem;
                    if ( i == j ) {
                        /*  reached atom on the top of nStackAtom[] stack  */
                        break;
                    }
                }
            } else
            if ( nTopStackAtom > 0 ) {
                j = (int)nStackAtom[nTopStackAtom-1];
                /* inherit nLowNumber */
                if ( nLowNumber[j] > nLowNumber[i] ) {
                    nLowNumber[j] = nLowNumber[i];
                }
            }
        } while ( --nTopStackAtom >= 0 );
    }

#endif /* }  BNS_MARK_ONLY_BLOCKS != 1 */

exit_function:
    if ( nStackAtom )
        inchi_free( nStackAtom );
    if ( nRingStack )
        inchi_free( nRingStack );
    if ( nDfsNumber )
        inchi_free( nDfsNumber );
    if ( nLowNumber )
        inchi_free( nLowNumber );
    if ( nBondStack )
        inchi_free( nBondStack );
    if ( cNeighNumb )
        inchi_free( cNeighNumb );
    return nNumRingSystems;
}

/*********************************************************************************/
int ReInitBnStructForAltBns( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bUnknAltAsNoStereo )
{
    Vertex v, v2;
    int ret, bond_type, num_to_test, j;
    BNS_EDGE   *pBond;
    BNS_VERTEX *pAtom;
    /* strip all t-groups and c-groups */
    num_to_test = 0;
    if ( bUnknAltAsNoStereo ) {
        for ( j = 0; j < pBNS->num_edges; j ++ ) {
            pBNS->edge[j].pass = 0;
        }
    }
    ret = ReInitBnStruct( pBNS, at, num_atoms, 0 );
    if ( ret || pBNS->num_atoms != num_atoms || pBNS->num_vertices != num_atoms || pBNS->num_bonds != pBNS->num_edges ) {
        ret = BNS_REINIT_ERR;
        goto exit_function;
    }
    /* eliminate bonds and fix st-caps */
    for ( v = 0; v < num_atoms; v ++ ) {
        pAtom = pBNS->vert + v;
        for ( j = 0; j < pAtom->valenceAltBns; j ++ ) {
            pBond = pBNS->edge + pAtom->iedge[j];
            if ( pBond->neighbor1 == v ) {
                bond_type = (at[v].bond_type[j] & BOND_TYPE_MASK);
                v2 = pBond->neighbor12 ^ v;
                if ( at[v].endpoint || at[v2].endpoint ) {
                    bond_type = 0; /* any bond to an endpoint considered non-stereogenic */
                }
#if( FIX_EITHER_DB_AS_NONSTEREO == 1 )
                if ( bUnknAltAsNoStereo ) {
                    if ( bond_type == BOND_ALTERN && at[v].bond_stereo[j] == STEREO_DBLE_EITHER ) {
                        bond_type = 0; /* treat unknown (Either) ALT bond as non-stereo */
                    }
                }
#endif
                switch ( bond_type ) {
                
                case BOND_ALTERN :
                    pBond->nBondTypeInpAltBns = BT_ALTERN_BOND;
                    num_to_test ++;
                    break;
                
                case BOND_ALT_123:
                case BOND_ALT_13 :
                case BOND_ALT_23 :
                    pBond->nBondTypeInpAltBns = BT_OTHER_ALTERN_BOND;
                    break;

                case BOND_TAUTOM :
                    pBond->nBondTypeInpAltBns = BT_TAUTOM_BOND;
                    break;
                
                case BOND_ALT12NS:
                    pBond->nBondTypeInpAltBns = BT_ALTERN_NS_BOND;
                    break;

                case 0:
                case BOND_SINGLE :
                case BOND_DOUBLE :
                case BOND_TRIPLE :
                    pBond->nBondTypeInpAltBns = BT_IGNORE_BOND;
                    break;

                default:
                    pBond->nBondTypeInpAltBns = BT_IGNORE_BOND;
                    break;

                }
                pBond->nBondNonStereoAltBns =
                pBond->nBlockNumberAltBns   =
                pBond->nNumAtInBlockAltBns  = 0;

#if( RESET_EDGE_FORBIDDEN_MASK == 1 )
                pBond->forbidden &= pBNS->edge_forbidden_mask;
#endif
            }
        }
        pAtom->bCutVertexAltBns         =
        pAtom->nRingSystemAltBns        =
        pAtom->nNumAtInRingSystemAltBns =
        pAtom->nBlockSystemAltBns       = 0;
    }

    return num_to_test;
exit_function:
    return ret;
}

/************************************************************************************/
int MarkNonStereoAltBns( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bUnknAltAsNoStereo )
{
    int       num_bonds = pBNS->num_bonds;
    int       ret;
    int       ibond, ib1, ib2;
    BNS_EDGE *pBond;
    Vertex    iat1, iat2;

    ret = 0;

    if ( pBNS->num_atoms != num_atoms || pBNS->num_vertices != num_atoms || pBNS->num_bonds != pBNS->num_edges ) {
        ret = BNS_REINIT_ERR;
        goto exit_function;
    }
    if ( bUnknAltAsNoStereo ) {
        for ( ibond=0; ibond < num_bonds; ibond ++ ) {
            pBond = pBNS->edge + ibond;
            if ( pBond->nBondTypeInpAltBns != BT_ALTERN_BOND && pBond->nBondTypeInpAltBns != BT_IGNORE_BOND ) {
                continue;
            }
            iat1 = pBond->neighbor1;
            iat2 = pBond->neighbor12 ^ iat1;
            ib1  = pBond->neigh_ord[0];
            ib2  = pBond->neigh_ord[1];
            if ( /* alt bond non-adjacent to a taut. endpoint: */
                  (pBond->nBondTypeInpAltBns == BT_ALTERN_BOND &&
                  pBond->nNumAtInBlockAltBns <= 3 )  /* non-ring bond */ ||
                 /* alt bond adjacent to a taut. endpoint: */
                  (pBond->nBondTypeInpAltBns == BT_IGNORE_BOND &&
                  (at[iat1].bond_type[ib1] & BOND_TYPE_MASK) == BOND_ALTERN )
                ) {
                if ( (at[iat1].bond_type[ib1] & BOND_TYPE_MASK) == BOND_ALTERN ) {
                    /* bond_type = BOND_ALT12NS; */
                    at[iat1].bond_stereo[ib1] =
                    at[iat2].bond_stereo[ib2] =STEREO_DBLE_EITHER;
                    ret ++;
                }
            } 
        }
    } else {
        for ( ibond=0; ibond < num_bonds; ibond ++ ) {
            pBond = pBNS->edge + ibond;
            if ( pBond->nBondTypeInpAltBns != BT_ALTERN_BOND && pBond->nBondTypeInpAltBns != BT_IGNORE_BOND ) {
                continue;
            }
            iat1 = pBond->neighbor1;
            iat2 = pBond->neighbor12 ^ iat1;
            ib1  = pBond->neigh_ord[0];
            ib2  = pBond->neigh_ord[1];
            if ( /* alt bond non-adjacent to a taut. endpoint: */
                 (pBond->nBondTypeInpAltBns == BT_ALTERN_BOND &&
                  pBond->nNumAtInBlockAltBns <= 3 ) /* non-ring bond */ ||
                 /* alt bond adjacent to a taut. endpoint: */
                 (pBond->nBondTypeInpAltBns == BT_IGNORE_BOND &&
                  (at[iat1].bond_type[ib1] & BOND_TYPE_MASK) == BOND_ALTERN ) 
                ) {
                at[iat1].bond_type[ib1] =
                at[iat2].bond_type[ib2] =BOND_ALT12NS;
                ret ++;
            } 
        }
    }

exit_function:

    return ret;
}

#if( READ_INCHI_STRING == 1 )
/*****************************************************************/
#ifndef RI_ERR_ALLOC
/* from ichirvrs.h */
#define RI_ERR_ALLOC   (-1)
#define RI_ERR_SYNTAX  (-2)
#define RI_ERR_PROGR   (-3)
#endif
/*****************************************************************/
int bHasChargedNeighbor( inp_ATOM *at, int iat )
{
    int i;
    for( i = 0; i < at[iat].valence; i ++ ) {
        if ( at[(int)at[iat].neighbor[i]].charge )
            return 1;
    }
    return 0;
}
/***********************************************************************************************
    *num_protons_to_add = nToBeRemovedByNormFromRevrs

    nToBeRemovedByNormFromRevrs > 0: less protons should be allowed to be
                                     added by the Normalization of the Reconstructed Structure
    nToBeRemovedByNormFromRevrs < 0: prepare more H(+) to be removed by
                                     the InChI Normalization of the Reconstructed Structure

    OrigStruct -> NormOrig + n(orig)*H(+) 
    RevrStruct -> NormRevr + n(revr)*H(+) 
    nToBeRemovedByNormFromRevrs = n(orig) - n(revr)  [each may be negative]
    
    n(orig) > n(revr) or nToBeRemovedByNormFromRevrs > 0 means:
    -----------------------------------------------------------
    - Too many protons were added by the Normalization to the Reconstructed Structure
      (a) n(revr) < 0 => protons were added while they should not have been added;
          Solution: "neutralize" (-) charged proton acceptors by moving charges to other atoms
                     on the condition ADP cannot add in another way;
      (b) n(orig) > n(revr) => 0  => too few protons were removed
          Solution: (the easiest) attach H(+) to =O or -N< or -N=
          Solution: move (+) from N or OH to an atom adjacent to (-) charge or to
                    an atom that is not N.
    
    n(orig) < n(revr) or nToBeRemovedByNormFromRevrs < 0 means:
    -----------------------------------------------------------
    - Too few protons were added by the Normalization to the Reconstructed Stucture
      (a) n(orig) < 0 => protons were not added while they should have been added;
          Solution: move (-) to O by replacing =O with -O(-)
      (b) 0 <= n(orig) < n(revr) => too many protons were removed

   Note: it is critically important to takr into account cumbersome Normalization
     Total Charge: if it is >= 0 then no H(+) may be removed from -OH or by ADP
     However, if N(+) is present then ADP will always try to remove a proton
*********************************************************************************************/
int AddRemoveProtonsRestr( inp_ATOM *at, int num_atoms, int *num_protons_to_add,
                           int nNumProtAddedByRestr, INCHI_MODE bNormalizationFlags,
                           int num_tg, int nChargeRevrs, int nChargeInChI )
{
    int i, j, ret = 0;
    int nAtTypeTotals[ATTOT_ARRAY_LEN];
    int   num_prot = *num_protons_to_add;
    int   type, mask, bSuccess, nTotCharge, nNumSuccess = 0;
    int max_j_Aa=-1, max_j_Ar=-1;

/* for the reference:

#define FLAG_NORM_CONSIDER_TAUT      ( FLAG_PROTON_NPO_SIMPLE_REMOVED | \
                                       FLAG_PROTON_NP_HARD_REMOVED    | \
                                       FLAG_PROTON_AC_SIMPLE_ADDED    | \
                                       FLAG_PROTON_AC_SIMPLE_REMOVED  | \
                                       FLAG_PROTON_AC_HARD_REMOVED    | \
                                       FLAG_PROTON_AC_HARD_ADDED      | \
                                       FLAG_PROTON_SINGLE_REMOVED     | \
                                       FLAG_PROTON_CHARGE_CANCEL    )

#define FLAG_FORCE_SALT_TAUT         ( FLAG_PROTON_NP_HARD_REMOVED  | \
                                       FLAG_PROTON_AC_HARD_REMOVED  | \
                                       FLAG_PROTON_AC_HARD_ADDED    )

*/
    /* if ChargeRevrs > nChargeInChI then we should prevent proton addition or facilitate proton removal
          a typical case is (=) on N or O instead of C(-)

       if ChargeRevrs < nChargeInChI then we should prevent proton removal or facilitate proton addition 
    */
    
    mark_at_type( at, num_atoms, nAtTypeTotals );
    for ( i = nTotCharge = 0; i < num_atoms; i ++ ) {
        nTotCharge += at[i].charge;
    }
    /* size for SimpleAddAcidicProtons() */
    for ( max_j_Aa = 0; AaTypMask[2*max_j_Aa]; max_j_Aa ++ )
        ;
    /* size for SimpleRemoveAcidicProtons */
    for ( max_j_Ar = 0; ArTypMask[2*max_j_Ar]; max_j_Ar ++ )
        ;
    if ( num_prot < 0 && nAtTypeTotals[ATTOT_TOT_CHARGE]-nNumProtAddedByRestr <= 0 ) {
        /* remove proton(s) */
        /* use test from SimpleAddAcidicProtons() to test whether removal of H(+) from =C-OH, etc. is correct */
        for ( i = 0; i < num_atoms && num_prot; i ++ ) {
            /* choose an atom */
            if ( at[i].sb_parity[0] || at[i].p_parity || at[i].charge ||
                 !at[i].num_H || at[i].radical || bHasChargedNeighbor( at, i ) ) {
                continue;
            }
            /* try to remove a proton and check whether InChI would add it back */
            at[i].charge --;
            at[i].num_H  --;
            type = GetAtomChargeType( at, i, NULL, &mask, 0 );
            at[i].charge ++;
            at[i].num_H  ++;

            if ( type ) {
                for ( bSuccess = 0, j = 0; j < max_j_Aa; j ++ ) {
                    if ( bSuccess = (type & AaTypMask[2*j]) && (mask && AaTypMask[2*j+1]) ) {
                        break; /* the proton may be added to this atom */
                    }
                }
                if ( bSuccess ) {
                    type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 1 ); /* subtract at[i] */
                    at[i].charge --;
                    at[i].num_H  --;
                    type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 0 ); /* add changed at[i] */
                    num_prot ++; /* success */
                    nNumSuccess ++;
                }
            }
        }
    }
    if ( num_prot < 0 && num_tg && nAtTypeTotals[ATTOT_TOT_CHARGE]-nNumProtAddedByRestr <= 0 ) {
        /* alternative proton removal: O=C-NH => (-)O-C=N, O and N are taut. endpoints */
        int endp2, centp, k, i0, k0;
        for ( i = 0; i < num_atoms; i ++ ) {
            /* choose an atom */
            if ( !at[i].endpoint || at[i].sb_parity[0] || at[i].p_parity ||
                 at[i].radical || at[i].charge || bHasChargedNeighbor( at, i ) ) {
                continue;
            }
            /* looking for tautomeric =O */
            if ( 1 != at[i].valence || BOND_TYPE_DOUBLE != at[i].bond_type[0] || at[i].num_H ||
                 2 != get_endpoint_valence( at[i].el_number ) ) {
                continue;
            }
            centp = at[i].neighbor[0];
            if ( at[centp].sb_parity[0] || at[centp].p_parity || !is_centerpoint_elem( at[centp].el_number ) ) {
                continue;
            }
            /* found a possible centerpoint, looking for -NH endpoint */
            for ( k = 0; k < at[centp].valence; k ++ ) {
                if ( at[centp].bond_type[k] != BOND_TYPE_SINGLE ) {
                    continue;
                }
                endp2 = at[centp].neighbor[k];
                if ( at[endp2].endpoint != at[i].endpoint ||
                     !at[endp2].num_H || at[endp2].charge ||
                     at[endp2].sb_parity[0] || at[endp2].p_parity ||
                     at[endp2].valence != at[endp2].chem_bonds_valence ||
                     3 != at[endp2].chem_bonds_valence + at[endp2].num_H ||
                     3 != get_endpoint_valence( at[endp2].el_number ) ) {
                    continue;
                }
                /* find bonds in reciprocal ajacency lists */
                for ( i0 = 0; i0 < at[centp].valence && i != at[centp].neighbor[i0]; i0 ++ )
                    ;
                for ( k0 = 0; k0 < at[endp2].valence && centp != at[endp2].neighbor[k0]; k0 ++ )
                    ;
                if ( i0 == at[centp].valence || k0 == at[endp2].valence ) {
                    return RI_ERR_PROGR;
                }
                /* -NH has been found */
                type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 1 ); /* subtract at[i] */
                type = GetAtomChargeType( at, endp2, nAtTypeTotals, &mask, 1 ); /* subtract at[endp2] */
                
                at[i].bond_type[0] --;
                at[centp].bond_type[i0] --;
                at[i].chem_bonds_valence --;
                at[i].charge --;

                at[endp2].bond_type[k0] ++;
                at[centp].bond_type[k] ++;
                at[endp2].chem_bonds_valence ++;
                at[endp2].num_H --;

                num_prot ++;
                nNumSuccess ++;
                
                type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 0 ); /* add at[i] */
                type = GetAtomChargeType( at, endp2, nAtTypeTotals, &mask, 0 ); /* add at[endp2] */
            }
        }
    }
    if ( num_prot > 0 ) {
        /* add protons */
        /* 1. Use test from SimpleRemoveAcidicProtons() to test whether addition of H(+) to =C-O(-), etc. is correct */
        for ( i = 0; i < num_atoms && num_prot && nAtTypeTotals[ATTOT_TOT_CHARGE]-nNumProtAddedByRestr >= 0; i ++ ) {
            /* choose an atom */
            if ( at[i].sb_parity[0] || at[i].p_parity || at[i].num_H ||
                 at[i].charge != -1 || at[i].radical || bHasChargedNeighbor( at, i ) ) {
                continue;
            }
            /* try to add a proton and check whether InChI would remove it back */
            at[i].charge ++;
            at[i].num_H ++;
            type = GetAtomChargeType( at, i, NULL, &mask, 0 );
            at[i].charge --;
            at[i].num_H --;
            
            if ( type ) {
                for ( bSuccess = 0, j = 0; j < max_j_Ar; j ++ ) {
                    if ( bSuccess = (type & ArTypMask[2*j]) && (mask && ArTypMask[2*j+1]) ) {
                        break;
                    }
                }
                if ( bSuccess ) {
                    type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 1 ); /* subtract at[i] */
                    at[i].charge ++;
                    at[i].num_H  ++;
                    type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 0 ); /* add changed at[i] */
                    num_prot --; /* success */
                    nNumSuccess ++;
                }
            }
        }
        /* 2. Use test from SimpleRemoveHplusNPO() */
        for ( i = 0; i < num_atoms && num_prot; i ++ ) {
            /* choose an atom */
            if ( at[i].sb_parity[0] || at[i].p_parity ||
                 at[i].charge || at[i].radical || bHasChargedNeighbor( at, i ) ) {
                continue;
            }
            /* try to add a proton and check whether InChI would remove it back */
            at[i].num_H ++;
            at[i].charge ++;
            bSuccess = (PR_SIMPLE_TYP & (type = GetAtomChargeType( at, i, NULL, &mask, 0 )) ) &&
                       (PR_SIMPLE_MSK & mask );
            at[i].num_H --;  /* failed */
            at[i].charge --;
            if ( bSuccess ) {
                type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 1 ); /* subtract at[i] */
                at[i].num_H ++;
                at[i].charge ++;
                type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 0 ); /* add changed at[i] */
                num_prot --;     /* succeeded */
                nNumSuccess ++;
            }
        }
    }

    if ( num_prot < 0 && (bNormalizationFlags & FLAG_PROTON_AC_HARD_ADDED) && 1 == num_tg &&
         nAtTypeTotals[ATTOT_TOT_CHARGE]-nNumProtAddedByRestr <= 0 ) {
        /* try to remove protons from tautomeric N (specific ADP must be present) */
        int nNumAcceptors_DB_O=0, nNumDonors_SB_NH=0, num_max, num_success;
        for ( i = 0; i < num_atoms; i ++ ) {
            /* choose an atom */
            if ( !at[i].endpoint || at[i].radical ||
                 at[i].sb_parity[0] || at[i].p_parity || bHasChargedNeighbor( at, i ) ) {
                continue;
            }
            type = GetAtomChargeType( at, i, NULL, &mask, 0 );
            if ( (type & AA_HARD_TYP_CO) && (mask & AA_HARD_MSK_CO) ) {
                nNumAcceptors_DB_O ++;
            } else
            if ( (type == ATT_ATOM_N ) && (mask == ATBIT_NP_H) && !at[i].charge &&
                 at[i].valence == at[i].chem_bonds_valence ) {
                nNumDonors_SB_NH ++;
            }
        }
        num_max = inchi_min( nNumAcceptors_DB_O, nNumDonors_SB_NH );
        for ( i = 0, num_success = 0; i < num_atoms && num_success < num_max && num_prot < 0; i ++ ) {
            /* choose an atom */
            if ( !at[i].endpoint|| at[i].radical || at[i].sb_parity[0] ||
                 at[i].p_parity || bHasChargedNeighbor( at, i ) ) {
                continue;
            }
            type = GetAtomChargeType( at, i, NULL, &mask, 0 );
            if ( (type == ATT_ATOM_N ) && (mask == ATBIT_NP_H) && !at[i].charge &&
                 at[i].valence == at[i].chem_bonds_valence ) {
                type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 1 ); /* subtract at[i] */
                at[i].num_H --;
                at[i].charge --;
                type = GetAtomChargeType( at, i, nAtTypeTotals, &mask, 0 ); /* add changed at[i] */
                num_prot ++;
                num_success ++;
                nNumSuccess ++;
            }
        }
    }
/*exit_function:*/
    *num_protons_to_add = num_prot;
    return ret<0? ret : nNumSuccess;    
}
/*****************************************************************/
int AddRemoveIsoProtonsRestr( inp_ATOM *at, int num_atoms, NUM_H num_protons_to_add[], int num_tg )
{
    int i, j, k, n, ret = 0;
    int   nNumSuccess = 0, min_at, max_at, num_H, num_iso_H, num_expl_H, num_expl_iso_H;
    int   iCurIso; /* 0=> 1H, 1=> D, 2=> T */
    int   iCurMode, iCurMode1, iCurMode2; /* 0=> Not Endpoints, 1=> Endpoints */
    static U_CHAR el_number_H = 0;

    /* distribute isotopes from  heaviest to lightest; pick up atoms in order 1. Not endpoints; 2. Endpoints */
    iCurMode1 = 0;
    iCurMode2 = num_tg ? 1 : 0;
    if ( !el_number_H ) {
        el_number_H = (U_CHAR) get_periodic_table_number( "H" );
    }
    for ( iCurMode = iCurMode1; iCurMode <= iCurMode2; iCurMode ++ ) {
        for ( iCurIso = 2; 0 <= iCurIso; iCurIso -- ) {
            /* check for isotopic H to add */
            if ( !num_protons_to_add[iCurIso] ) {
                continue;
            }
            if ( 0 > num_protons_to_add[iCurIso] ) {
                ret = RI_ERR_PROGR;
                goto exit_function;
            }
            /* limits for atom scanning */
            min_at = 0;
            max_at = num_atoms;
            /* cycle withio the limits */
            for ( i = min_at; i < max_at && 0 < num_protons_to_add[iCurIso]; i ++ ) {
                /* pick an atom */
                if ( iCurMode ) {
                    if ( at[i].endpoint )
                        j = i;  /* atom number */
                    else
                        continue;
                } else
                if ( !at[i].endpoint &&
                     1 == bHeteroAtomMayHaveXchgIsoH( at, i ) ) { /* atom number */
                    j = i;
                } else
                if ( at[i].el_number == el_number_H && at[i].charge == 1 &&
                    !at[i].valence && !at[i].radical && !at[i].iso_atw_diff ) {
                    /* proton, not isotopic; make it isotopic */
                    at[i].iso_atw_diff = 1 + iCurIso;
                    num_protons_to_add[iCurIso] --;
                    nNumSuccess ++;
                    continue;
                } else {
                    continue;
                }
                /* j is the atom number */
                /* count implicit H */
                num_H      = at[j].num_H;
                num_iso_H  = NUM_ISO_H(at,j);
                while ( num_H > 0 && num_protons_to_add[iCurIso] > 0  ) {
                    /* substitute one implicit H with an isotopic atom H */
                    at[j].num_iso_H[iCurIso] ++;
                    at[j].num_H --;
                    num_protons_to_add[iCurIso] --;
                    num_H --;
                    num_iso_H ++;
                    nNumSuccess ++;
                }
                /* count explicit H */
                num_expl_H = num_expl_iso_H = 0;
                for ( k = 0; k < at[j].valence && num_atoms <= (n=at[j].neighbor[k]); k ++ ) {
                    num_expl_H     += (0 == at[n].iso_atw_diff);
                    num_expl_iso_H += (0 != at[n].iso_atw_diff);
                }
                while ( num_expl_H > 0 && num_protons_to_add[iCurIso] > 0  ) {
                    /* substitute one explicit H with an isotopic atom H */
                    n = at[j].neighbor[num_expl_H];
                    if ( at[n].iso_atw_diff ) {
                        ret = RI_ERR_PROGR;
                        goto exit_function;
                    }
                    at[n].iso_atw_diff = 1 + iCurIso;
                    num_expl_H --;
                    num_expl_iso_H ++;
                    num_protons_to_add[iCurIso] --;
                    nNumSuccess ++;
                }
            }
        }
    }
exit_function:
    return ret<0? ret : nNumSuccess;    
}

#endif

