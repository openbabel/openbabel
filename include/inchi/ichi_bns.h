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


#ifndef _INCHI_BNS_H_
#define _INCHI_BNS_H_


#include "incomdef.h"
#include "inpdef.h"


/*#define FIX_SRU_CYCLIZING_PS_BONDS_IN_BNS 1*/

#define BN_MAX_ALTP  16
/*#define MAX_VERTEX 1024*/ /* including s; if vert[] has num_vert then MAX_VERTEX has (2*num_vert+2+FIRST_INDX) elements */

/* forward declarations */

struct BalancedNetworkStructure;
struct BalancedNetworkData;
struct tagTautomerGroupsInfo;
struct tagChargeGroupsInfo;
struct BN_AtomsAtTautGroup;
struct tagSaltChargeCandidate;

/* define BNS types */

typedef int  Vertex;
typedef int  EdgeIndex;
typedef int  Edge[2];         /* Edge[0] = vertex1, Edge[1] = iedge or -(1+vertex1) if vertex2 = s or t */
typedef int  BNS_IEDGE;
typedef int  EdgeFlow;
typedef int  VertexFlow;


#define BNS_EDGE_FORBIDDEN_MASK  1
#define BNS_EDGE_FORBIDDEN_TEMP  2
#define BNS_EDGE_FORBIDDEN_TEST  4

/* BNS vertex types */

#define BNS_VERT_TYPE_ATOM          0x0001
#define BNS_VERT_TYPE_ENDPOINT      0x0002  /* attribute */
#define BNS_VERT_TYPE_TGROUP        0x0004
#define BNS_VERT_TYPE_C_POINT       0x0008
#define BNS_VERT_TYPE_C_GROUP       0x0010
#define BNS_VERT_TYPE_SUPER_TGROUP  0x0020
#define BNS_VERT_TYPE_TEMP          0x0040

#define BNS_VERT_TYPE__AUX          0x0080  /* vertex added to build charge substructures */
#define BNS_VERT_TYPE_C_NEGATIVE    0x0100  /* negative charge group; attribute, should be used with BNS_VERT_TYPE_C_GROUP */
#define BNS_VERT_TYPE_ACID          0x0200  /* only for this type are allowed paths: t_group-atom-c_group_neg (path_TACN) */
#define BNS_VERT_TYPE_CARBON_GR     0x0400  /* charge of carbon atom; should be used with BNS_VT_C_POS, BNS_VT_C_NEG */
#define BNS_VERT_TYPE_METAL_GR      0x0800  /* metal atom group; may be used alone or with BNS_VT_M_POS, BNS_VT_M_NEG */

#define BNS_VERT_TYPE_ANY_GROUP    (BNS_VERT_TYPE_TGROUP | BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_SUPER_TGROUP)

/* InChI->Structure */

#define BNS_VT_C_POS     BNS_VERT_TYPE_C_GROUP                               /* positive charge group, heteroat */
#define BNS_VT_C_NEG     (BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE)  /* negative charge group, heteroat */
#define BNS_VT_C_POS_C   (BNS_VT_C_POS | BNS_VERT_TYPE_CARBON_GR)            /* positive charge group, C, Si, Ge, Sn */
#define BNS_VT_C_NEG_C   (BNS_VT_C_NEG | BNS_VERT_TYPE_CARBON_GR)            /* negative charge group, C, Si, Ge, Sn */
#define BNS_VT_C_POS_M   (BNS_VT_C_POS | BNS_VERT_TYPE_METAL_GR)             /* positive charge group, metal */
#define BNS_VT_C_NEG_M   (BNS_VT_C_NEG | BNS_VERT_TYPE_METAL_GR)             /* negative charge group, metal */
#define BNS_VT_M_GROUP   BNS_VERT_TYPE_METAL_GR                              /* metal-group, flower vertex */

#define BNS_VT_C_POS_ALL  (BNS_VERT_TYPE_SUPER_TGROUP | BNS_VERT_TYPE_C_GROUP)    /* supergroup (+) */
#define BNS_VT_C_NEG_ALL  (BNS_VT_C_POS_ALL | BNS_VERT_TYPE_C_NEGATIVE) /* supergroup (-) */

#define BNS_VT_CHRG_STRUCT  (BNS_VERT_TYPE__AUX | BNS_VERT_TYPE_TEMP)          /* ChargeStruct vertex */
#define BNS_VT_YVCONNECTOR  BNS_VERT_TYPE__AUX                                 /* group connection */

#define IS_BNS_VT_C_OR_CSUPER_GR(X) ((X) & BNS_VT_C_POS)
#define IS_BNS_VT_C_GR(X)           (((X) & BNS_VT_C_POS_ALL) == BNS_VERT_TYPE_C_GROUP)
#define IS_BNS_VT_CM_GR(X)          (((X) & BNS_VT_C_POS_M) == BNS_VT_C_POS_M) /* metal charge group */
#define IS_BNS_VT_M_GR(X)           ((X) == BNS_VERT_TYPE_METAL_GR )  /* metal flower base or vertices */
#define IS_BNS_VT_YVCONNECTOR(X)    (((X) & BNS_VERT_TYPE__AUX) && !((X) & BNS_VERT_TYPE_TEMP))
#define IS_BNS_VT_CHRG_STRUCT(X)    (((X) & BNS_VERT_TYPE__AUX) &&  ((X) & BNS_VERT_TYPE_TEMP))
#define IS_BNS_VT_ATOM(X)           ((X) & BNS_VERT_TYPE_ATOM)

#define BNS_ADD_SUPER_TGROUP 1  /* reserve one more edge for a t-group to connect to a single super-t-group */
#define NUM_KINDS_OF_GROUPS  2  /* 1 accounts for t-group kind, one more 1 accounts for c-group kind */

#define BNS_ADD_ATOMS        2  /* max. number of fictitious atoms to add (except t-gtoups) */
#define BNS_ADD_EDGES        1  /* max. number of edges to add to each atom (except edges to a t-group or c-group) */

typedef enum tagAltPathConst {
    iALTP_MAX_LEN,    /* 0 */
    iALTP_FLOW,       /* 1 */
    iALTP_PATH_LEN,   /* 2 */
    iALTP_START_ATOM, /* 3 */
    iALTP_END_ATOM,   /* 4 */
    iALTP_NEIGHBOR,   /* 5 */
    iALTP_HDR_LEN = iALTP_NEIGHBOR
} ALT_CONST;

#define ALTP_PATH_LEN(altp)             (altp)[iALTP_PATH_LEN].number  /* number of bonds = number of atoms-1*/
#define ALTP_END_ATOM(altp)             (altp)[iALTP_END_ATOM].number
#define ALTP_START_ATOM(altp)           (altp)[iALTP_START_ATOM].number
#define ALTP_THIS_ATOM_NEIGHBOR(altp,X) (altp)[iALTP_NEIGHBOR+(X)].ineigh[0]  /* 0 <= X < path_len */
#define ALTP_NEXT_ATOM_NEIGHBOR(altp,X) (altp)[iALTP_NEIGHBOR+(X)].ineigh[1]
#define ALTP_CUR_THIS_ATOM_NEIGHBOR(altp) (altp)[iALTP_NEIGHBOR+ALTP_PATH_LEN(altp)].ineigh[0]  /* 0 <= X < path_len */
#define ALTP_CUR_NEXT_ATOM_NEIGHBOR(altp) (altp)[iALTP_NEIGHBOR+ALTP_PATH_LEN(altp)].ineigh[1]
#define ALTP_NEXT(altp)                 (++ALTP_PATH_LEN(altp))
#define ALTP_PREV(altp)                 (--ALTP_PATH_LEN(altp))
#define ALTP_MAY_ADD(altp)              (iALTP_NEIGHBOR + (altp)[iALTP_PATH_LEN].number < (altp)[iALTP_MAX_LEN].number)
#define ALTP_ALLOCATED_LEN(altp)        (altp)[iALTP_MAX_LEN].number
#define ALTP_DELTA(altp)                 (altp)[iALTP_FLOW].flow[0]
#define ALTP_OVERFLOW(altp)             (altp)[iALTP_FLOW].flow[1]

#define Vertex_s 0
#define Vertex_t 1

#define NO_VERTEX    -2
#define BLOSSOM_BASE -1

#define ADD_CAPACITY_RADICAL        1   /* add capacity to radical */

#define MAX_BOND_EDGE_CAP           2  /* triple bond */
#define AROM_BOND_EDGE_CAP          1
#define MAX_TGROUP_EDGE_CAP         2  /* -NH2 provides max. capacity */

/* edge to s or t */
#define EDGE_FLOW_ST_MASK       0x3fff  /* mask for flow */
#define EDGE_FLOW_ST_PATH       0x4000  /* mark: the edge belongs to the augmenting path */

/* edges between other vertices */
/* EdgeFlow WAS defined as S_SHORT; change from S_CHAR made 9-23-2005 */
#define EDGE_FLOW_MASK          0x3fff  /* mask for flow */
#define EDGE_FLOW_PATH          0x4000  /* mark: the edge belongs to the augmenting path */

/*********************************************************************************/
#if ( ADD_CAPACITY_RADICAL == 1 )  /* { */
/*  -- do not treat triplets as moving dots -- 2004-02-18 --
#define MAX_AT_FLOW(X) (((X).chem_bonds_valence - (X).valence)+\
                       ((is_centerpoint_elem((X).el_number)||get_endpoint_valence((X).el_number))?\
                           (((X).radical==RADICAL_DOUBLET)+2*((X).radical==RADICAL_TRIPLET)):0))
*/
#define MAX_AT_FLOW(X) (((X).chem_bonds_valence - (X).valence)+\
                       ((is_centerpoint_elem((X).el_number)||get_endpoint_valence((X).el_number))?\
                           (((X).radical==RADICAL_DOUBLET)/*+2*((X).radical==RADICAL_TRIPLET)*/):0))


#else /* } ADD_CAPACITY_RADICAL { */

#define MAX_AT_FLOW(X) (((X).chem_bonds_valence - (X).valence)

#endif  /* } ADD_CAPACITY_RADICAL */

/**************************** BNS_EDGE ************************************/
typedef struct BnsEdge {
    AT_NUMB   neighbor1;                      /* the smaller neighbor */
    AT_NUMB   neighbor12;                     /* neighbor1 ^ neighbor2 */
    AT_NUMB   neigh_ord[2];                   /* ordering number of the neighbor: [0]: at<neighbor, [1]: at>neighbor */
    EdgeFlow  cap;                            /* Edge capacity */
    EdgeFlow  cap0;                           /* Initial edge capacity */
    EdgeFlow  flow;                           /* Edge flow */
    EdgeFlow  flow0;                          /* Initial flow */
    /*S_CHAR    delta; */
    S_CHAR    pass;                           /* number of times changed in AugmentEdge() */
    S_CHAR    forbidden;
} BNS_EDGE;

/**************************** BNS_ST_EDGE ************************************/
typedef struct BnsStEdge {
    VertexFlow cap;                            /* Edge capacity */
    VertexFlow cap0;                           /* Initial edge capacity */
    VertexFlow flow;                           /* Edge flow */
    VertexFlow flow0;                          /* Initial edge flow */
    S_CHAR     pass;                           /* number of times changed in AugmentEdge() */
    /*S_CHAR     delta; */
} BNS_ST_EDGE;

/**************************** BNS_VERTEX ************************************/
typedef struct BnsVertex {

    BNS_ST_EDGE st_edge;                     /* 0,1 capacity and flow of the edge to s or t */
    AT_NUMB     type;                        /* 2, atom, t-group, or added atom: BNS_VERT_TYPE_TGROUP, etc. */
    AT_NUMB     num_adj_edges;               /* 3, actual number of neighbors incl. t-groups, excl. s or t */
    AT_NUMB     max_adj_edges;               /* 4, including reserved */
    /*S_CHAR      path_neigh[2];*/               /* 5 found path information */
    /* indexes of Edges */
    BNS_IEDGE  *iedge;                       /* 6 a pointer to the array of edge indexes adjacent to this vertex */
}BNS_VERTEX;

/**************************** BNS_ALT_PATH ************************************/
typedef union BnsAltPath {
    VertexFlow   flow[2];
    Vertex       number;
    AT_NUMB      ineigh[2];
} BNS_ALT_PATH;

/**************************** BN_STRUCT ************************************/
typedef struct BalancedNetworkStructure {

    int num_atoms;        /* number of real atoms */
    /*int len_atoms; */       /* size of filled with real atoms portion of the BNS_VERTEX data */
    int num_added_atoms;  /* number of added fictitious atoms */
    int nMaxAddAtoms;     /* max. number of atoms to add (not including t-groups) */
    int num_c_groups;     /* number of added c-groups */
    int num_t_groups;     /* number of added t-groups */
    int num_vertices;     /* total number currently in effect; includes t-groups and added atoms */
    /*int len_vertices; */    /* allocation size for BNS_VERTEX data */
    int num_bonds;        /* number of real bonds/2 = number of edges between real atoms */
    int num_edges;        /* number of currently in effect */
    int num_iedges;       /* added 9-16-2005; used only in InChI Reversing */
    int num_added_edges;  /* number of added edges (not including edges to t-groups) */
    int nMaxAddEdges;     /* max. number edges of add to each atom (not including edges to t-groups) */

    int max_vertices;     /* allocation size for BNS_VERTEX structures */
    int max_edges;        /* allocation size for edge[]; iedge has length 2*max_edges */
    int max_iedges;       /* allocation size for iedge */

    int tot_st_cap;       /* not used, only calculated */
    int tot_st_flow;      /* not used, only calculated */

    int len_alt_path;     /* length of alt_path[] */

    int bNotASimplePath;  /* alternating path traversed same bond 2 times in opposite directions */
    int bChangeFlow;      /* actually change flow */

    BNS_VERTEX    *vert;  /* vertices */
    BNS_EDGE      *edge;  /* edges */
    BNS_IEDGE     *iedge;
    BNS_ALT_PATH  *alt_path;           /* current altp[] element */
    BNS_ALT_PATH  *altp[BN_MAX_ALTP];  /* keep alt. paths */

    int            max_altp;
    int            num_altp;

    INCHI_MODE     *pbTautFlags;     /* carry it through all functions; never NULL */
    INCHI_MODE     *pbTautFlagsDone; /* carry it through all functions; never NULL */
    AT_NUMB        type_TACN; /* BNS_VERT_TYPE_ACID: if non-zero than only for it path type_T-type_TACN-type_CN allowed */
    AT_NUMB        type_T;    /* BNS_VERT_TYPE_TGROUP */
    AT_NUMB        type_CN;   /* BNS_VERT_TYPE_C_GROUP | BNS_VERT_TYPE_C_NEGATIVE */
    S_CHAR         edge_forbidden_mask;
    /* v. 1.05 */
    struct tagINCHI_CLOCK *ic;
    struct tagInchiTime *ulTimeOutTime;
} BN_STRUCT;

/********************* BN_DATA *******************************************/
typedef enum tagBnsRadSrchMode {
    RAD_SRCH_NORM = 0,   /* normal search for normalization */
    RAD_SRCH_FROM_FICT = 1    /* search from fict. vertices to atoms */
} BRS_MODE;
typedef struct BalancedNetworkData {
    Vertex          *BasePtr;    /*[MAX_VERTEX];  pointer toward the base of C(v) */
    Edge            *SwitchEdge; /*[MAX_VERTEX];  a pair of vertices and an edge, implemented here as [*][2] array */
    S_CHAR          *Tree;       /*[MAX_VERTEX];  indicates presence in ScanQ, T, T', s-reachability */
    Vertex          *ScanQ;      /*[MAX_VERTEX];  contains the set S of s-reachable vertices */
    int             QSize;       /* index of the  last element added to ScanQ */
    Vertex          *Pu;         /*[MAX_VERTEX/2+1] */
    Vertex          *Pv;         /*[MAX_VERTEX/2+1] */
    int             max_num_vertices; /* allocation size of all except Pu, Pv */
    int             max_len_Pu_Pv;    /* allocation size of Pu and Pv */
#if ( BNS_RAD_SEARCH == 1 )
    Vertex         *RadEndpoints; /*[MAX_VERTEX*/
    int             nNumRadEndpoints;
    EdgeIndex      *RadEdges;
    int             nNumRadEdges;
    int             nNumRadicals;
    BRS_MODE        bRadSrchMode; /* 1 => connect fict. vertices-radicals to the accessible atoms */
#endif
} BN_DATA;

/* internal array size */
#define MAX_ALT_AATG_ARRAY_LEN 127
/* detected endpoint markings */
#define AATG_MARK_IN_PATH      1  /* atom in path detected by the BNS */
#define AATG_MARK_WAS_IN_PATH  2  /* found to be in path before next level */
/* output */
#define AATG_MARK_MAIN_TYPE   4  /* atom O-"salt" */
#define AATG_MARK_OTHER_TYPE  8  /* other atom to be tested */

struct tagTautomerGroupsInfo;   /* forward declaration */

/******************** atoms in alt path through taut group ****************/
typedef struct BN_AtomsAtTautGroup {
    int      nAllocLen;
    int      nNumFound;
    int      nNumMainAdj2Tgroup;
    int      nNumOthersAdj2Tgroup;
    AT_NUMB *nEndPoint;            /* original t-group number */
    S_CHAR  *nMarkedAtom;          /* atom mark, see AATG_MARK_* */
    int     *nAtTypeTotals;
    struct tagTautomerGroupsInfo *t_group_info;
} BN_AATG;


/************ store changes in flow and capacity to test a bond ****************/

typedef struct tagBNS_FLOW_CHANGES {
    BNS_IEDGE  iedge;
    EdgeFlow   flow;
    EdgeFlow   cap;
    Vertex     v1;
    VertexFlow cap_st1;
    VertexFlow flow_st1;
    Vertex     v2;
    VertexFlow cap_st2;
    VertexFlow flow_st2;
} BNS_FLOW_CHANGES;


#define ALT_PATH_MODE_TAUTOM     1
#define ALT_PATH_MODE_CHARGE     2
#define ALT_PATH_MODE_4_SALT     3   /* mark alt bonds along the path */
#define ALT_PATH_MODE_4_SALT2    4   /* mark alt bonds along the path, path to taut. group fict. vertex if exists */
#define ALT_PATH_MODE_REM2H_CHG  5   /* remove 2 H along alt. path AH-=-BH => A=-=B and change bonds to alternating */
#define ALT_PATH_MODE_ADD2H_CHG  6   /* add 2 H along alt. path A=-=B => AH-=-BH    and change bonds to alternating */
#define ALT_PATH_MODE_REM2H_TST  7   /* test-remove 2 H along alt. path AH-=-BH => A=-=B; restore changed bonds */
#define ALT_PATH_MODE_ADD2H_TST  8   /* test-add 2 H along alt. path A=-=B => AH-=-BH; restore changed bonds */
#define ALT_PATH_MODE_REM_PROTON 9   /* remove proton, adjust bonds, charges, H-counts 2004-03-05 */
#if ( KETO_ENOL_TAUT == 1 )
#define ALT_PATH_MODE_TAUTOM_KET 10  /* same as ALT_PATH_MODE_TAUTOM, applies to C=-OH or CH-=O; H may be (-) */
#endif

#if ( TAUT_PT_22_00 == 1 )
#define ALT_PATH_MODE_TAUTOM_PT_22_00 11
#endif
#if ( TAUT_PT_16_00 == 1 )
#define ALT_PATH_MODE_TAUTOM_PT_16_00 12
#endif
#if ( TAUT_PT_06_00 == 1 )
#define ALT_PATH_MODE_TAUTOM_PT_06_00 13
#endif
#if ( TAUT_PT_39_00 == 1 )
#define ALT_PATH_MODE_TAUTOM_PT_39_00 14
#endif
#if ( TAUT_PT_13_00 == 1 )
#define ALT_PATH_MODE_TAUTOM_PT_13_00 15
#endif
#if ( TAUT_PT_18_00 == 1 )
#define ALT_PATH_MODE_TAUTOM_PT_18_00 16
#endif

typedef U_SHORT  bitWord;
#define BIT_WORD_MASK  ((bitWord)~0)

typedef struct tagNodeSet {
    bitWord **bitword;
    int num_set; /* number of sets */
    int len_set; /* number of bitWords in each set */
} NodeSet;


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


/*********************************************************************************
  bChangeFlow:
      1 => change flow inside the BNS search
      3 => change flow inside the BNS search and undo the flow change in the BNS structure here
      4 => change bonds in the structure according to the flow
      8 => make altern. bonds in the structure

  Note: (bChangeFlow & 1) == 1 is needed for multiple runs
**********************************************************************************/

/* "EF" = "Edge Flow" */
#define BNS_EF_CHNG_FLOW      1  /* change Balanced Network (BN) flow inside the BNS search */
#define BNS_EF_RSTR_FLOW      2  /* undo BN flow changes after BNS */
#define BNS_EF_CHNG_RSTR      (BNS_EF_CHNG_FLOW | BNS_EF_RSTR_FLOW)
#define BNS_EF_CHNG_BONDS     4  /* change bonds in the structure according to the BN flow */
#define BNS_EF_ALTR_BONDS     8  /* make altern. bonds in the structure if the flow has changed */
#define BNS_EF_UPD_RAD_ORI   16  /* update BN flow0 & Atom radical values:
                                    flow0 := flow, radical:=st_cap - st_flow */
#define BNS_EF_SET_NOSTEREO  32  /* in combination with BNS_EF_ALTR_BONDS only:
                                    ALT12 bond cannot be stereogenic */
#define BNS_EF_UPD_H_CHARGE  64  /* update charges and H-counts according to change flow to c- and t-group vertices */

#define BNS_EF_SAVE_ALL     (BNS_EF_CHNG_FLOW | BNS_EF_CHNG_BONDS | BNS_EF_UPD_RAD_ORI)
#define BNS_EF_ALTR_NS      (BNS_EF_ALTR_BONDS | BNS_EF_SET_NOSTEREO)

#define BNS_EF_RAD_SRCH     128  /* search for rafical paths closures */


    struct tagCANON_GLOBALS;

    int  NodeSetCreate( struct tagCANON_GLOBALS *pCG, NodeSet *pSet, int n, int L );
    void NodeSetFree( struct tagCANON_GLOBALS *pCG, NodeSet *pSet );

    int  IsNodeSetEmpty( NodeSet *cur_nodes, int k );
    int  DoNodeSetsIntersect( NodeSet *cur_nodes, int k1, int k2 );
    void AddNodeSet2ToNodeSet1( NodeSet *cur_nodes, int k1, int k2 );
    void NodeSetFromRadEndpoints( struct tagCANON_GLOBALS *pCG, NodeSet *cur_nodes, int k, /*Node *v*/ Vertex RadEndpoints[], int num_v );
    void RemoveFromNodeSet( struct tagCANON_GLOBALS *pCG, NodeSet *cur_nodes, int k, Vertex v[], int num_v );
    int  AddNodesToRadEndpoints( struct tagCANON_GLOBALS *pCG, NodeSet *cur_nodes, int k, Vertex RadEndpoints[], Vertex vRad, int nStart, int nLen );


    int nExists2AtMoveAltPath( struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD,
                               struct BN_AtomsAtTautGroup *pAATG, inp_ATOM *at, int num_atoms,
                               int jj2, int jj1, struct tagSaltChargeCandidate *s_candidate, int nNumCandidates,
                               AT_NUMB *nForbiddenAtom, int nNumForbiddenAtoms );

    int bExistsAltPath( struct tagCANON_GLOBALS *pCG,
                        struct BalancedNetworkStructure *pBNS,
                        struct BalancedNetworkData *pBD,
                        struct BN_AtomsAtTautGroup *pAATG,
                        inp_ATOM *at,
                        int num_atoms,
                        int nVertDoubleBond,
                        int nVertSingleBond,
                        int path_type );
    int bExistsAnyAltPath( struct tagCANON_GLOBALS *pCG, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD,
                           inp_ATOM *at, int num_atoms, int nVertDoubleBond, int nVertSingleBond, int path_type );
    int AddTGroups2BnStruct( struct tagCANON_GLOBALS *pCG, struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_atoms,
                             struct tagTautomerGroupsInfo *tgi );
    int AddSuperTGroup2BnStruct( struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_atoms,
                                 struct tagTautomerGroupsInfo *tgi );
    int AddCGroups2BnStruct( struct tagCANON_GLOBALS *pCG, struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_atoms,
                             struct tagChargeGroupsInfo *cgi );

    int ReInitBnStruct( struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_at, int bRemoveGroupsFromAtoms );
    int ReInitBnStructAddGroups( struct tagCANON_GLOBALS *pCG, struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_atoms,
                                 struct tagTautomerGroupsInfo *tgi, struct tagChargeGroupsInfo *cgi );


    int DisconnectTestAtomFromTGroup( struct BalancedNetworkStructure *pBNS, int v1, int *pv2, BNS_FLOW_CHANGES *fcd );
    int DisconnectTGroupFromSuperTGroup( struct BalancedNetworkStructure *pBNS, int v1, int *pv1, int *pv2,
                                         BNS_FLOW_CHANGES *fcd );
    int ReconnectTestAtomToTGroup( struct BalancedNetworkStructure *pBNS, int v1, int v2, int ie, BNS_FLOW_CHANGES *fcd );

    int bIsHardRemHCandidate( inp_ATOM *at, int i, int *cSubType );

    /* moved from ichi_bns.c 2005-08-23 */
    int RunBalancedNetworkSearch( BN_STRUCT *pBNS, BN_DATA *pBD, int bChangeFlow );
    BN_STRUCT* AllocateAndInitBnStruct( inp_ATOM *at, int num_atoms, int nMaxAddAtoms, int nMaxAddEdges, int max_altp, int *num_changed_bonds );
    BN_STRUCT* DeAllocateBnStruct( BN_STRUCT *pBNS );
    int ReInitBnStructAltPaths( BN_STRUCT *pBNS );
    int ReInitBnStructForMoveableAltBondTest( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms );
    void ClearAllBnDataVertices( Vertex *v, Vertex value, int size );
    void ClearAllBnDataEdges( Edge *e, Vertex value, int size );
    BN_DATA *DeAllocateBnData( BN_DATA *pBD );
    BN_DATA *AllocateAndInitBnData( int max_num_vertices );
    int ReInitBnData( BN_DATA *pBD );
    int SetForbiddenEdges( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int edge_forbidden_mask,
                           int nebend, int *ebend );
    /* main function: find augmenting path */
    int BalancedNetworkSearch( BN_STRUCT* pBNS, BN_DATA *pBD, int bChangeFlow );

    int SetRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode );
    int SetRadEndpoints2( struct tagCANON_GLOBALS *pCG, BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode );
    int RemoveRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at );

    int AddRemoveProtonsRestr( inp_ATOM *at, int num_atoms, int *num_protons_to_add,
                               int nNumProtAddedByRestr, INCHI_MODE bNormalizationFlags,
                               int num_tg, int nChargeRevrs, int nChargeInChI );
    int AddRemoveIsoProtonsRestr( inp_ATOM *at, int num_atoms, NUM_H num_protons_to_add[], int num_tg );


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif



#endif    /* _INCHI_BNS_H_ */
