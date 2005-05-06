/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.00
 * April 13, 2005
 * Developed at NIST
 */

#ifndef __INCHI_BNS_H___
#define __INCHI_BNS_H___

#define BN_MAX_ALTP  16
#define MAX_VERTEX 1024 /* including s; if vert[] has num_vert then MAX_VERTEX has (2*num_vert+2+FIRST_INDX) elements */

/* forward declarations */

struct BalancedNetworkStructure;
struct BalancedNetworkData;
struct tagTautomerGroupsInfo;
struct tagChargeGroupsInfo;
struct BN_AtomsAtTautGroup;
struct tagSaltChargeCandidate;

/* define BNS types */

typedef S_SHORT  Vertex;
typedef S_SHORT  EdgeIndex;
typedef S_SHORT  Edge[2];         /* Edge[0] = vertex1, Edge[1] = iedge or -(1+vertex1) if vertex2 = s or t */
typedef S_SHORT BNS_IEDGE;
typedef S_SHORT EdgeFlow;
typedef S_SHORT VertexFlow;


#define BNS_EDGE_FORBIDDEN_MASK  1
#define BNS_EDGE_FORBIDDEN_TEMP  2

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

} BN_STRUCT;

/********************* BN_DATA *******************************************/
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

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


int nExists2AtMoveAltPath( struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD,
                           struct BN_AtomsAtTautGroup *pAATG, inp_ATOM *at, int num_atoms,
                           int jj2, int jj1, struct tagSaltChargeCandidate *s_candidate, int nNumCandidates,
                           AT_NUMB *nForbiddenAtom, int nNumForbiddenAtoms);
int bExistsAltPath( struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD,
                    struct BN_AtomsAtTautGroup *pAATG, inp_ATOM *at, int num_atoms, int nVertDoubleBond, int nVertSingleBond, int path_type );
int bExistsAnyAltPath( struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD,
                       inp_ATOM *at, int num_atoms, int nVertDoubleBond, int nVertSingleBond, int path_type );
int AddTGroups2BnStruct( struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_atoms,
                         struct tagTautomerGroupsInfo *tgi  );
int AddSuperTGroup2BnStruct( struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_atoms,
                             struct tagTautomerGroupsInfo *tgi  );
int AddCGroups2BnStruct( struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_atoms,
                         struct tagChargeGroupsInfo *cgi );

int ReInitBnStruct( struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_at, int bRemoveGroupsFromAtoms );
int ReInitBnStructAddGroups( struct BalancedNetworkStructure *pBNS, inp_ATOM *at, int num_atoms,
                             struct tagTautomerGroupsInfo *tgi, struct tagChargeGroupsInfo *cgi );


int DisconnectTestAtomFromTGroup( struct BalancedNetworkStructure *pBNS, int v1, int *pv2, BNS_FLOW_CHANGES *fcd );
int DisconnectTGroupFromSuperTGroup( struct BalancedNetworkStructure *pBNS, int v1, int *pv1, int *pv2,
                                     BNS_FLOW_CHANGES *fcd );
int ReconnectTestAtomToTGroup( struct BalancedNetworkStructure *pBNS, int v1, int v2, int ie, BNS_FLOW_CHANGES *fcd );

int bIsHardRemHCandidate(  inp_ATOM *at, int i, int *cSubType );

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __INCHI_BNS_H___ */
