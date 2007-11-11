/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02-beta
 * August 23, 2007
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


/* input/output format */
#ifndef __INPDEF_H__
#define __INPDEF_H__

/*^^^ */
#include "mode.h"
#include "incomdef.h" 
#include "ichidrp.h"
/*^^^ */

#define bDrawingLabelLeftShift endpoint    /* for drawing only */
typedef S_SHORT ST_CAP_FLOW;

/* inp_ATOM::at_type */
#define ATT_NONE         0x0000
#define ATT_ACIDIC_CO    0x0001
#define ATT_ACIDIC_S     0x0002
#define ATT_OO           0x0004
#define ATT_ZOO          0x0008
#define ATT_NO           0x0010
#define ATT_N_O          0x0020
#define ATT_ATOM_N       0x0040
#define ATT_ATOM_P       0x0080
#define ATT_OTHER_NEG_O  0x0100
#define ATT_OTHER_ZO     0x0200   /* -Z=O or =Z=O */
#define ATT_OH_MINUS     0x0400   /* OH(-), O=O,S,Se,Te */
#define ATT_O_PLUS       0x0800   /* -OH2(+), =OH(+), -OH(+)-, OH3(+), =O(+)-, etc; O=O,S,Se,Te */
#define ATT_PROTON       0x1000
#define ATT_HalAnion     0x2000
#define ATT_HalAcid      0x4000

#define AT_FLAG_ISO_H_POINT 0x01  /* may have isotopic H */

#define PERIODIC_NUMBER_H  1

#ifndef NUMH
#define NUM_ISO_H(AT,N) (AT[N].num_iso_H[0]+AT[N].num_iso_H[1]+AT[N].num_iso_H[2])
#define NUMH(AT,N)     (AT[N].num_H+NUM_ISO_H(AT,N))
#endif

#define FlagSC_0D  1  /* bUsed0DParity */
#define FlagSB_0D  2  /* bUsed0DParity */

#define SB_PARITY_FLAG  0x38 /* mask for disconnected metal parity if it is different */
#define SB_PARITY_SHFT  3    /* number of right shift bits to get disconnected metal parity */
#define SB_PARITY_MASK  0x07
#define SB_PARITY_1(X) (X & SB_PARITY_MASK)  /* refers to connected structure */
#define SB_PARITY_2(X) (((X) >> SB_PARITY_SHFT) & SB_PARITY_MASK) /* refers to connected structure */



typedef struct tagInputAtom {
    char          elname[ATOM_EL_LEN]; /* chem. element name */
    U_CHAR        el_number;               /* number of the element in the Periodic Table */
    AT_NUMB       neighbor[MAXVAL];        /* positions (from 0) of the neighbors in the inp_ATOM array */
    AT_NUMB       orig_at_number;          /* original atom number */
    AT_NUMB       orig_compt_at_numb;      /* atom number within the component before terminal H removal */
    S_CHAR        bond_stereo[MAXVAL];     /* 1=Up,4=Either,6=Down; this atom is at the pointing wedge,
                                              negative => on the opposite side; 3=Either double bond  */
    U_CHAR        bond_type[MAXVAL];       /* 1..4; 4="aromatic", should be discouraged on input */

    S_CHAR        valence;                 /* number of bonds = number of neighbors */
    S_CHAR        chem_bonds_valence;      /* sum of bond types (type 4 needs special treatment) */
    S_CHAR        num_H;                   /* number of implicit hydrogens including D and T    */
    S_CHAR        num_iso_H[NUM_H_ISOTOPES]; /* number of implicit 1H, 2H(D), 3H(T) < 16 */
    S_CHAR        iso_atw_diff;            /* =0 => natural isotopic abundances  */
                                           /* >0 => (mass) - (mass of the most abundant isotope) + 1 */
                                           /* <0 => (mass) - (mass of the most abundant isotope) */
    S_CHAR        charge;                  /* charge */
    S_CHAR        radical;                 /* RADICAL_SINGLET, RADICAL_DOUBLET, or RADICAL_TRIPLET */
    S_CHAR        bAmbiguousStereo;
    S_CHAR        cFlags;                  /* AT_FLAG_ISO_H_POINT */
    AT_NUMB       at_type;                 /* ATT_NONE, ATT_ACIDIC */
    AT_NUMB       component;               /* number of the structure component > 0 */
    AT_NUMB       endpoint;                /* id of a tautomeric group */
    AT_NUMB       c_point;                 /* id of a positive charge group */
    double        x;
    double        y;
    double        z;
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

#if( FIND_RING_SYSTEMS == 1 )
    S_CHAR  bCutVertex;
    AT_NUMB nRingSystem;
    AT_NUMB nNumAtInRingSystem;
    AT_NUMB nBlockSystem;

#if( FIND_RINS_SYSTEMS_DISTANCES == 1 )
    AT_NUMB nDistanceFromTerminal;       /* terminal atom or ring system has 1, next has 2, etc. */
#endif

#endif
} inp_ATOM;

typedef struct tagOrigAtom {
    /* initially filled out by MolfileToOrigAtom */
    /* may be changed by disconnecting salts and disconnecting metals */
    inp_ATOM *at;
    int num_dimensions;
    int num_inp_bonds;
    int num_inp_atoms;
    /* may be changed by disconnecting salts and disconnecting metals */
    int num_components;    /* set by MarkDisconnectedComponents() and disconnecting metals */
    int bDisconnectSalts;  /* whether salt disconnection is possible */
    int bDisconnectCoord;  /* 0 if no disconnection needed else (Num Implicit H to disconnect)+1 */
#if( bRELEASE_VERSION == 0 )
    int bExtract;
#endif
    AT_NUMB *nCurAtLen;      /* has max_num_components elements */
    AT_NUMB *nOldCompNumber; /* 0 or component number in previous numbering */
    int      nNumEquSets;  /* number of found component equivalence sets */
    AT_NUMB *nEquLabels; /* num_inp_atoms elements, value>0 marks atoms in the set #value  */
    AT_NUMB *nSortedOrder; /* num_components elements, values = 1..num_components; only if num_components > 1  */
    int bSavedInINCHI_LIB[INCHI_NUM];
    int bPreprocessed[INCHI_NUM];
    MOL_COORD *szCoord;
} ORIG_ATOM_DATA;

typedef struct tagOriginalStruct {
    int num_atoms;
    char *szAtoms;
    char *szBonds;
    char *szCoord;
} ORIG_STRUCT;

typedef struct tagAtomParmsForDrawing {
    char      at_string[ATOM_INFO_LEN];
    int       DrawingLabelLeftShift;
    int       DrawingLabelLength;
    AT_NUMB   nCanonNbr;               /* if zero then do not use all data for the atom */
    AT_NUMB   nCanonEquNbr;
    AT_NUMB   nTautGroupCanonNbr;
    AT_NUMB   nTautGroupEquNbr;
    S_CHAR    cFlags;                  /* AT_FLAG_ISO_H_POINT */
#ifdef DISPLAY_DEBUG_DATA
    int       nDebugData;
#endif
    S_CHAR    cHighlightTheAtom;
    S_CHAR    cStereoCenterParity;
    S_CHAR    cStereoBondParity[MAX_STEREO_BONDS];
    S_CHAR    cStereoBondWarning[MAX_STEREO_BONDS];
    S_CHAR    cStereoBondNumber[MAX_STEREO_BONDS];
} inf_ATOM;


#define INF_STEREO_ABS         0x0001
#define INF_STEREO_REL         0x0002
#define INF_STEREO_RAC         0x0004
#define INF_STEREO_NORM        0x0008
#define INF_STEREO_INV         0x0010
#define INF_STEREO             0x0020
#define INF_STEREO_ABS_REL_RAC (INF_STEREO_ABS | INF_STEREO_REL | INF_STEREO_RAC)
#define INF_STEREO_NORM_INV    (INF_STEREO_NORM | INF_STEREO_INV)

#define MAX_LEN_REMOVED_PROTONS 128

typedef struct tagInfoAtomData {
    inf_ATOM  *at;
    int        num_at;
    AT_NUMB    StereoFlags;
    AT_NUMB    num_components;
    AT_NUMB    *pStereoFlags;

    int        nNumRemovedProtons;
    int        num_removed_iso_H; /* number of exchangable isotopic H */
    NUM_H      num_iso_H[NUM_H_ISOTOPES]; /* number of exchangable isotopic H */
    char       szRemovedProtons[MAX_LEN_REMOVED_PROTONS];
} INF_ATOM_DATA;

typedef struct tagInputAtomData {
    inp_ATOM *at;
    inp_ATOM *at_fixed_bonds; /* tautomeric case, added or removed H */
    int       num_at;
    int       num_removed_H;
    int       num_bonds;
    int       num_isotopic;
    int       bExists;
    int       bDeleted;
    int       bHasIsotopicLayer;
    int       bTautomeric;
    int       bTautPreprocessed;
    int       nNumRemovedProtons;
    NUM_H     nNumRemovedProtonsIsotopic[NUM_H_ISOTOPES]; /* isotopic composition of removed protons, not included in num_iso_H[] */
    NUM_H     num_iso_H[NUM_H_ISOTOPES]; /* isotopic H on tautomeric atoms and those in nIsotopicEndpointAtomNumber */
    INCHI_MODE  bTautFlags;
    INCHI_MODE  bTautFlagsDone;
    INCHI_MODE  bNormalizationFlags;
} INP_ATOM_DATA;
typedef INP_ATOM_DATA INP_ATOM_DATA2[TAUT_NUM];

typedef struct tagNormCanonFlags {
    INCHI_MODE  bTautFlags[INCHI_NUM][TAUT_NUM];
    INCHI_MODE  bTautFlagsDone[INCHI_NUM][TAUT_NUM];
    INCHI_MODE  bNormalizationFlags[INCHI_NUM][TAUT_NUM];
    int        nCanonFlags[INCHI_NUM][TAUT_NUM];
} NORM_CANON_FLAGS;

typedef struct tagCompositeAtomData {
    inp_ATOM *at;
    int       num_at;
    int       num_removed_H;
    int       num_bonds;
    int       num_isotopic;
    int       bExists;
    int       bDeleted;    /* unused */
    int       bHasIsotopicLayer;
    int       bTautomeric;
    int       nNumRemovedProtons;
    NUM_H     nNumRemovedProtonsIsotopic[NUM_H_ISOTOPES]; /* isotopic composition of removed protons, not included in num_iso_H[] */
    NUM_H     num_iso_H[NUM_H_ISOTOPES]; /* isotopic H on tautomeric atoms and those in nIsotopicEndpointAtomNumber */

    AT_NUMB   *nOffsetAtAndH;
    int       num_components;
} COMP_ATOM_DATA;
/*
typedef COMP_ATOM_DATA COMP_ATOM_DATA3[TAUT_NUM+1];
*/
#define ADD_LEN_STRUCT_FPTRS 100  /* allocation increments */
typedef long INCHI_FPTR;
typedef struct tagStructFptrs {
    INCHI_FPTR *fptr;      /* input:  fptr[cur_fptr]   = file pointer to the structure to read */
                          /* output: fptr[cur_fptr+1] = file pointer to the next structure or EOF */
    int        len_fptr;  /* allocated length of fptr */
    int        cur_fptr;  /* input: k-1 to read the kth struct, k = 1, 2, 3,...; left unchanged; struct number := cur_fptr+1 */
    int        max_fptr;  /* length of the filled out portion of fptr */
} STRUCT_FPTRS;

#define FLAG_INP_AT_CHIRAL         1
#define FLAG_INP_AT_NONCHIRAL      2
#define FLAG_SET_INP_AT_CHIRAL     4
#define FLAG_SET_INP_AT_NONCHIRAL  8

/* BILLY 8/6/04 */
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

int MolfileToInpAtom( FILE *inp_molfile, int bDoNotAddH, inp_ATOM **at, MOL_COORD **szCoord, int max_num_at,
                      int *num_dimensions, int *num_bonds, const char *pSdfLabel, char *pSdfValue,
                      long *Id, long *lMolfileNumber, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr );
int MolfileToOrigAtom( FILE *inp_molfile, ORIG_ATOM_DATA *orig_at_data, int bMergeAllInputStructures,
                       int bGetOrigCoord, int bDoNotAddH,
                       const char *pSdfLabel, char *pSdfValue, long *lSdfId, long *lMolfileNumber,
                       INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr );
int INChIToOrigAtom( FILE *inp_molfile, ORIG_ATOM_DATA *orig_at_data, int bMergeAllInputStructures,
                       int bGetOrigCoord, int bDoNotAddH, INPUT_TYPE nInputType,
                       char *pSdfLabel, char *pSdfValue, long *lSdfId,
                       INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr );

int MarkDisconnectedComponents( ORIG_ATOM_DATA *orig_at_data, int bProcessOldCompNumbers );
int DisconnectSalts( ORIG_ATOM_DATA *orig_inp_data, int bDisconnect );
int DisconnectMetals( ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone );
int bMayDisconnectMetals( ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone );
int bHasMetalAtom( ORIG_ATOM_DATA *orig_inp_data );
int FixAdjacentRadicals( int num_inp_atoms, inp_ATOM *at ); /* FIX_ADJ_RAD == 1 */
int fix_odd_things( int num_atoms, inp_ATOM *at, int bFixBug );
int post_fix_odd_things( int num_atoms, inp_ATOM *at );
int remove_ion_pairs( int num_atoms, inp_ATOM *at );

int bFoundFeature( inp_ATOM *at, int num_atoms );
int CopyMOLfile(FILE *inp_file, long fPtrStart, long fPtrEnd, 
                FILE *prb_file, /*^^^ was: INCHI_FILE */
                long nNumb);

void FreeInpAtom( inp_ATOM **at );
void FreeInfAtom( inf_ATOM **at );
void FreeOrigAtData( ORIG_ATOM_DATA *orig_at_data );
void FreeInpAtomData( INP_ATOM_DATA *inp_at_data );
void FreeCompAtomData( COMP_ATOM_DATA *inp_at_data );
void FreeInfoAtomData( INF_ATOM_DATA *inf_at_data );

int FixUnkn0DStereoBonds(inp_ATOM *at, int num_at);

inf_ATOM *CreateInfAtom( int num_atoms );
inp_ATOM *CreateInpAtom( int num_atoms );

int CreateInfoAtomData( INF_ATOM_DATA *inf_at_data, int num_atoms, int num_components );
int AllocateInfoAtomData( INF_ATOM_DATA *inf_at_data, int num_atoms, int num_components );
int DuplicateInfoAtomData( INF_ATOM_DATA *inf_at_data_to, const INF_ATOM_DATA *inf_at_data_from);
int CreateInpAtomData( INP_ATOM_DATA *inp_at_data, int num_atoms, int create_at_fixed_bonds );
int CreateCompAtomData( COMP_ATOM_DATA *inp_at_data, int num_atoms, int num_components, int bIntermediateTaut );
#ifndef INCHI_ANSI_ONLY
int DisplayInputStructure( char *szOutputString, inp_ATOM  *at, INF_ATOM_DATA *inf_at_data, int num_at, DRAW_PARMS *dp );
#endif
void PrintFileName( const char *fmt, 
                   FILE *output_file, /*^^^ was: INCHI_FILE */
                   const char *szFname );
void MySleep( unsigned long ms );

#ifndef __ICHITIME_H__
struct tagInchiTime;
int bInchiTimeIsOver( struct tagInchiTime *TickEnd );
#endif

int get_endpoint_valence( U_CHAR el_number );
#if( KETO_ENOL_TAUT == 1 )
int get_endpoint_valence_KET( U_CHAR el_number );
#endif

#if( TEST_RENUMB_ATOMS == 1 )  /*  { */
int CopyInpAtomData( INP_ATOM_DATA *dest_inp_at_data, INP_ATOM_DATA *src_inp_at_data );
void RenumbInpAtomData( INP_ATOM_DATA *dest_inp_at_data, INP_ATOM_DATA *src_inp_at_data, AT_RANK *new_ord );
void MakeNewOrd( int num_atoms, AT_RANK *new_ord );
#endif

int ReconcileAllCmlBondParities( inp_ATOM *at, int num_atoms, int bDisconnected );


/* BILLY 8/6/04 */
#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif  /* __INPDEF_H__ */
