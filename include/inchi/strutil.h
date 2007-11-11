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


#ifndef __STRUTIL_H__
#define __STRUTIL_H__

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/* forward declaration */
struct tagTautomerGroupsInfo;

int ExtractConnectedComponent(  inp_ATOM *at, int num_at, int component_number, inp_ATOM *component_at );
int SetConnectedComponentNumber( inp_ATOM *at, int num_at, int component_number );
INChI     *Alloc_INChI( inp_ATOM *at, int num_at, int *found_num_bonds, int *found_num_isotopic, int nAllocMode );
int       Free_INChI(INChI **ppINChI);
int       Free_INChI_Members(INChI *pINChI);
int       Free_INChI_Stereo( INChI_Stereo *pINChI_Stereo );

INChI_Aux *Alloc_INChI_Aux( int num_at, int num_isotopic_atoms, int nAllocMode, int bOrigData );
int       Free_INChI_Aux( INChI_Aux **ppINChI_Aux );
int  Create_INChI( INChI **ppINChI, INChI_Aux **ppINChI_Aux, ORIG_ATOM_DATA *orig_inp_data,
                  inp_ATOM *inp_at, INP_ATOM_DATA *inp_norm_data[2],
                  int num_inp_at, INCHI_MODE nUserMode,
                  INCHI_MODE *pbTautFlags, INCHI_MODE *pbTautFlagsDone,
                  struct tagInchiTime *ulMaxTime, struct tagTautomerGroupsInfo *ti_out, char *pStrErrStruct);
int FillOutInfAtom(inp_ATOM *norm_at, INF_ATOM_DATA *inf_norm_at_data, int init_num_at,
                   int num_removed_H, int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H *nNumRemovedProtonsIsotopic, int bIsotopic,
                   INChI *pINChI, INChI_Aux *pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode );
int FillOutCompositeCanonInfAtom(COMP_ATOM_DATA *composite_norm_data, INF_ATOM_DATA *inf_norm_at_data,
                                 int bIsotopic, int bTautomeric,
                                 PINChI2 *pINChI2, PINChI_Aux2 *pINChI_Aux2, int bAbcNumbers, INCHI_MODE nMode);
#if( FIX_DALKE_BUGS == 1 )
char *AllocateAndFillHillFormula( INChI *pINChI );
#endif

typedef enum tagInchiDiffBits {
    IDIF_PROBLEM       = 0x00000001, /* severe: at least one InChI does not exist */
    IDIF_NUM_AT        = 0x00000001, /* severe: different number of atoms in the skeleton */
    IDIF_ATOMS         = 0x00000001, /* severe: diiferent types of skeleton atoms */
    IDIF_NUM_EL        = 0x00000001, /* severe: formulas differ in another element */
    IDIF_CON_LEN       = 0x00000001, /* severe: different connection table lengths */
    IDIF_CON_TBL       = 0x00000001, /* severe: different connection tables */
    IDIF_POSITION_H    = 0x00000002, /* difference in non-taut (Mobile-H) or all H (Fixed-H) location/number */
    IDIF_MORE_FH       = 0x00000004, /* extra fixed H */
    IDIF_LESS_FH       = 0x00000008, /* missing fixed H */
    IDIF_MORE_H        = 0x00000010, /* formulas differ in number of H */
    IDIF_LESS_H        = 0x00000020, /* formulas differ in number of H */
    /*IDIF_TAUT_LEN      = 0x00000008,*/ /* different lengths of tautomer lists */
    IDIF_NO_TAUT       = 0x00000040, /* restored structure has no taut groups while the original InChI has some */
    IDIF_WRONG_TAUT    = 0x00000080, /* restored has tautomerism while the original does not have it */
    IDIF_SINGLE_TG     = 0x00000100, /* restored has 1 taut. group while the original InChI has multiple tg */
    IDIF_MULTIPLE_TG   = 0x00000200, /* restored has multiple tg while the original InChI has only one tg */
    IDIF_NUM_TG        = 0x00000400, /* different number of tautomeric groups */
    /*IDIF_LESS_TG_ENDP  = 0x00000200,*/ /* restores structure has less taut. endpoints */
    /*IDIF_MORE_TG_ENDP  = 0x00000400,*/ /* restores structure has more taut. endpoints */
    IDIF_EXTRA_TG_ENDP = 0x00000800, /* extra tautomeric endpoint(s) in restored structure */
    IDIF_MISS_TG_ENDP  = 0x00001000, /* one or more tg endpoint is not in the restored structure */
    IDIF_DIFF_TG_ENDP  = 0x00002000, /* lists of tg endpoints are different */
    IDIF_TG            = 0x00004000, /* different tautomeric groups */
    IDIF_NUM_ISO_AT    = 0x00008000, /* ?severe: restored struct. has different number of isotopic atoms */
    IDIF_ISO_AT        = 0x00010000, /* ?severe: restored struct. has different locations/isotopes of isotopic atoms */
    IDIF_CHARGE        = 0x00020000, /* restored structure has different charge */
    IDIF_REM_PROT      = 0x00040000, /* proton(s) removed/added from the restored structure */
    IDIF_REM_ISO_H     = 0x00080000, /* isotopic H removed */
    IDIF_SC_INV        = 0x00100000, /* restores structure has different inversion stereocenter mark */
    IDIF_SC_PARITY     = 0x00200000, /* restored structure has stereoatoms or allenes with different parity */
    IDIF_SC_EXTRA_UNDF = 0x00400000, /* restored structure has extra undefined stereocenter(s) */
    IDIF_SC_EXTRA      = 0x00800000, /* restored structure has extra stereocenter(s) */
    IDIF_SC_MISS_UNDF  = 0x01000000,  /* restored structure has not some undefined stereocenter(s) */
    IDIF_SC_MISS       = 0x02000000, /* restored structure has not some stereocenters that are not undefined */
    IDIF_SB_PARITY     = 0x04000000, /* restored structure has stereobonds or cumulenes with different parity */
    IDIF_SB_EXTRA_UNDF = 0x08000000, /* restored structure has extra undefined stereobond(s) */
    IDIF_SB_EXTRA      = 0x10000000, /* restored structure has extra stereobond(s) */
    IDIF_SB_MISS_UNDF  = 0x20000000, /* restored structure has not some undefined stereocenters */
    IDIF_SB_MISS       = 0x40000000  /* restored structure has not some stereobonds that are not undefined */
} IDIF;


#define IDIFF_SB (IDIF_SB_PARITY | IDIF_SB_EXTRA_UNDF | IDIF_SB_EXTRA | IDIF_SB_MISS_UNDF | IDIF_SB_MISS)
#define IDIFF_SC (IDIF_SC_PARITY | IDIF_SC_EXTRA_UNDF | IDIF_SC_EXTRA | IDIF_SC_MISS_UNDF | IDIF_SC_MISS)

#define IDIFF_CONSTIT (IDIF_POSITION_H | IDIF_MORE_FH | IDIF_LESS_FH | IDIF_MORE_H | IDIF_LESS_H |\
                       IDIF_NO_TAUT | IDIF_WRONG_TAUT | IDIF_SINGLE_TG | IDIF_MULTIPLE_TG | \
                       IDIF_NUM_TG | IDIF_EXTRA_TG_ENDP | IDIF_MISS_TG_ENDP | IDIF_TG | \
                       IDIF_NUM_ISO_AT | IDIF_ISO_AT | IDIF_CHARGE | IDIF_REM_PROT | IDIF_REM_ISO_H |\
                       IDIF_DIFF_TG_ENDP)
#define IDIFF_STEREO  (IDIF_SC_INV | IDIF_SC_PARITY | IDIF_SC_EXTRA_UNDF | IDIF_SC_EXTRA | \
                       IDIF_SC_MISS_UNDF | IDIF_SC_MISS | IDIF_SB_PARITY | IDIF_SB_EXTRA_UNDF |\
                       IDIF_SB_EXTRA | IDIF_SB_MISS_UNDF | IDIF_SB_MISS)


/*************************************************************************************/
#define ICR_MAX_ENDP_IN1_ONLY 32
#define ICR_MAX_ENDP_IN2_ONLY 32
#define ICR_MAX_DIFF_FIXED_H  32
#define ICR_MAX_SB_IN1_ONLY   32
#define ICR_MAX_SB_IN2_ONLY   32
#define ICR_MAX_SC_IN1_ONLY   32
#define ICR_MAX_SC_IN2_ONLY   32
#define ICR_MAX_SB_UNDF       32
#define ICR_MAX_SC_UNDF       32
typedef struct tagInChICompareResults {

    INCHI_MODE flags;

    int tot_num_H1;
    int tot_num_H2;
    int num_taut_H1;
    int num_taut_H2;
    int num_taut_M1;
    int num_taut_M2;
    
    /* 1 => InChI from reversed struct. 2 => input InChI */

    AT_NUMB endp_in1_only[ICR_MAX_ENDP_IN1_ONLY]; /* endpoint canonical number = index+1 */
    int     num_endp_in1_only;

    AT_NUMB endp_in2_only[ICR_MAX_ENDP_IN2_ONLY]; /* endpoint canonical number = index+1 */
    int     num_endp_in2_only;

    AT_NUMB diff_pos_H_at[ICR_MAX_DIFF_FIXED_H];  /* non-tautomeric H */
    S_CHAR  diff_pos_H_nH[ICR_MAX_DIFF_FIXED_H];
    int     num_diff_pos_H;

    AT_NUMB fixed_H_at1_more[ICR_MAX_DIFF_FIXED_H];  /* extra fixed_H */
    S_CHAR  fixed_H_nH1_more[ICR_MAX_DIFF_FIXED_H];
    int     num_fixed_H1_more;

    AT_NUMB fixed_H_at2_more[ICR_MAX_DIFF_FIXED_H];  /* missed fixed_H */
    S_CHAR  fixed_H_nH2_more[ICR_MAX_DIFF_FIXED_H];
    int     num_fixed_H2_more;

    AT_NUMB sc_in1_only[ICR_MAX_SC_IN1_ONLY];
    int     num_sc_in1_only;
    AT_NUMB sc_in2_only[ICR_MAX_SC_IN2_ONLY];
    int     num_sc_in2_only;

    AT_NUMB sb_in1_only[ICR_MAX_SB_IN1_ONLY];
    int     num_sb_in1_only;
    AT_NUMB sb_in2_only[ICR_MAX_SB_IN2_ONLY];
    int     num_sb_in2_only;

    AT_NUMB sb_undef_in1_only[ICR_MAX_SC_UNDF];
    int     num_sb_undef_in1_only;
    AT_NUMB sb_undef_in2_only[ICR_MAX_SC_UNDF];
    int     num_sb_undef_in2_only;
    
    AT_NUMB sc_undef_in1_only[ICR_MAX_SB_UNDF];
    int     num_sc_undef_in1_only;
    AT_NUMB sc_undef_in2_only[ICR_MAX_SB_UNDF];
    int     num_sc_undef_in2_only;

} ICR;


INCHI_MODE CompareReversedINChI2( INChI *i1 /* InChI from reversed struct */, INChI *i2 /* input InChI */, INChI_Aux *a1, INChI_Aux *a2, ICR *picr, int *err );
int CompareIcr( ICR *picr1, ICR *picr2, INCHI_MODE *pin1, INCHI_MODE *pin2, INCHI_MODE mask );

int CompareReversedINChI( INChI *i1, INChI *i2, INChI_Aux *a1, INChI_Aux *a2 );
const char *CompareReversedInchiMsg( int code );

#define EQL_EXISTS   1
#define EQL_SP3      2
#define EQL_SP3_INV  4
#define EQL_SP2      8
int Eql_INChI_Stereo( INChI_Stereo *s1, int eql1, INChI_Stereo *s2, int eql2, int bRelRac );

int Eql_INChI_Isotopic( INChI *i1, INChI *i2 );

#define EQL_EQU     0
#define EQL_EQU_TG  1
#define EQL_EQU_ISO 2
int Eql_INChI_Aux_Equ( INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2 );

#define EQL_NUM      0
#define EQL_NUM_INV  1
#define EQL_NUM_ISO  2
int Eql_INChI_Aux_Num( INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2 );

int bHasEquString( AT_NUMB *LinearCT, int nLenCT );

int CompINChINonTaut2(const void *p1, const void *p2);
int CompINChITaut2(const void *p1, const void *p2);
int CompINChI2(const INCHI_SORT *p1, const INCHI_SORT *p2, int bTaut, int bCompareIsotopic);
int CompINChITautVsNonTaut(const INCHI_SORT *p1, const INCHI_SORT *p2, int bCompareIsotopic);


typedef enum tagDiffINChISegments { /* r = repetitive, n = non-repetitive */
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

typedef enum tagDiffINChILayers {
    DIFL_M,  /* 0 main layer */
    DIFL_MI, /* 1 main isotopic */
    DIFL_F,  /* 2 fixed-H */
    DIFL_FI, /* 3 fixed-H isotopic */
    DIFL_LENGTH /* number of layers */
} DIF_LAYERS;

/* Value meaning */
typedef enum tagMarkDiff {
    DIFV_BOTH_EMPTY = 0, /* both this and the component in the preceding namesake segment are empty  */
    DIFV_EQL2PRECED = 1, /* equal to the component in the preceding namesake segment                 */
    DIFV_NEQ2PRECED = 2, /* different from the component in the preceding namesake segment           */
    DIFV_IS_EMPTY   = 4, /* is empty while the preceding namesake segment is not empty               */
    DIFV_FI_EQ_MI   = 8, /* FI stereo component is equal to the component in the MI namesake segment */
                         /* while M and F components are empty                                       */
    /* decision_F = bitmask: bits that should not be present */
    /* decision_T = bitmask: at least one of the bits should be present */
    /* decision = true if( !( BITS & decision_F ) && ( BITS & decision_F ) ) */
    DIFV_OUTPUT_EMPTY_T = (DIFV_IS_EMPTY), /* bits present for empty segment output */
    DIFV_OUTPUT_EMPTY_F = (DIFV_EQL2PRECED | DIFV_NEQ2PRECED | DIFV_FI_EQ_MI), /* bits NOT present */

    DIFV_OUTPUT_OMIT_F  = (DIFV_NEQ2PRECED | DIFV_IS_EMPTY), /* bits NOT present for omitting */

    DIFV_OUTPUT_FILL_T  = (DIFV_EQL2PRECED | DIFV_NEQ2PRECED | DIFV_FI_EQ_MI)

} DIF_VALUES;

typedef enum tagINChISegmAction {
    INCHI_SEGM_OMIT  = 0,
    INCHI_SEGM_FILL  = 1, /* the value is used in str_LineEnd() */
    INCHI_SEGM_EMPTY = 2  /* the value is used in str_LineEnd() */
} INCHI_SEGM_ACTION;

int CompINChILayers(const INCHI_SORT *p1, const INCHI_SORT *p2, char sDifSegs[][DIFS_LENGTH] );
int MarkUnusedAndEmptyLayers( char sDifSegs[][DIFS_LENGTH] );
int INChI_SegmentAction( char cDifSegs );

#define FLAG_SORT_PRINT_TRANSPOS_BAS   1 /* transposition in the main InChI layer */
#define FLAG_SORT_PRINT_TRANSPOS_REC   2 /* transposition in the reconnected InChI layer */
#define FLAG_SORT_PRINT_NO_NFIX_H_BAS  4 /* no fixed H non-isotopic in the main InChI layer */
#define FLAG_SORT_PRINT_NO_NFIX_H_REC  8 /* no fixed H non-isotopic in the reconnected InChI layer */
#define FLAG_SORT_PRINT_NO_IFIX_H_BAS 16 /* no fixed H isotopic in the main InChI layer */
#define FLAG_SORT_PRINT_NO_IFIX_H_REC 32 /* no fixed H isotopic in the the reconnected InChI layer */
#define FLAG_SORT_PRINT_ReChI_PREFIX  64 /* Output ReChI instead of InChI */

int OutputINChI1( char *pStr, int nStrLen, INCHI_SORT *pINChISortTautAndNonTaut[][TAUT_NUM], int iINChI,
                ORIG_STRUCT *pOrigStruct,
                int bDisconnectedCoord, int bOutputType, int bINChIOutputOptions, int bXml, int bAbcNumbers,
                int bCtPredecessors, int bNoStructLabels,
                int num_components[], int num_non_taut[], int num_taut[],
                INCHI_FILE *output_file, INCHI_FILE *log_file, int num_input_struct,
                const char *szSdfLabel, const char *szSdfValue, long lSdfId, int *pSortPrintINChIFlags );
int OutputINChI2( char *pStr, int nStrLen, INCHI_SORT *pINChISortTautAndNonTaut[][TAUT_NUM], int iINChI,
                ORIG_STRUCT *pOrigStruct,
                int bDisconnectedCoord, int bOutputType, int bINChIOutputOptions, int bXml, int bAbcNumbers,
                int bCtPredecessors, int bNoStructLabels,
                int num_components[], int num_non_taut[], int num_taut[],
                INCHI_FILE *output_file, INCHI_FILE *log_file, int num_input_struct,
                const char *szSdfLabel, const char *szSdfValue, long lSdfId, int *pSortPrintINChIFlags );
int SaveEquComponentsInfoAndSortOrder ( int iINChI, INCHI_SORT *pINChISort[TAUT_NUM], int *num_components,
                                        ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data,
                                        COMP_ATOM_DATA composite_norm_data[TAUT_NUM],
                                        int bCompareComponents );
int OutputINChIXmlRootStartTag( INCHI_FILE *output_file );
int OutputINChIXmlRootEndTag( INCHI_FILE *output_file );
int OutputINChIXmlError( INCHI_FILE *output_file, char *pStr, int nStrLen, int ind,
                        /*int nErrorNumber,*/ char *szErrorText, int bError );
int OutputINChIPlainError( INCHI_FILE *output_file, char *pStr, int nStrLen,
                           char *pErrorText, int bError );
int OutputINChIXmlStructStartTag( INCHI_FILE *output_file, char *pStr, int ind /* indent*/, int nStrLen, int bNoStructLabels,
                                 int num_input_struct, const char *szSdfLabel, const char *szSdfValue );
int OutputINChIXmlStructEndTag( INCHI_FILE *output_file, char *pStr, int nStrLen, int ind );

int GetInpStructErrorType(INPUT_PARMS *ip, int err, char *pStrErrStruct, int num_inp_atoms );
int ProcessStructError( INCHI_FILE *output_file, INCHI_FILE *log_file, /*int err,*/ char *pStrErrStruct, int nErrorType,
                         int *bXmlStructStarted, long num_inp, INPUT_PARMS *ip, char *pStr, int nStrLen );

int bNumHeterAtomHasIsotopicH( inp_ATOM *atom, int num_atoms );

int WriteToSDfile( const INP_ATOM_DATA *inp_at_data, INCHI_FILE* fcb, const char* name, const char* comment,
                   const char *szLabel, const char *szValue );
int WriteOrigAtomDataToSDfile( const ORIG_ATOM_DATA *inp_at_data, INCHI_FILE* fcb, const char* name, const char* comment,
                   int bChiral, int bAtomsDT, const char *szLabel, const char *szValue);
int bIsMetalSalt( inp_ATOM *at, int i );


extern const char gsMissing[];
extern const char gsEmpty[];
extern const char gsSpace[];
extern const char gsEqual[];
/*
#define gsMissing  "is missing"
#define gsEmpty    ""
#define gsSpace    " "
#define gsEqual    "="
*/
/* format string for SDF_LBL_VAL(L,V): %s%s%s%s (four strings) */
/*#define SDF_LBL_VAL(L,V)  ((L)&&(L)[0])?gsSpace:gsEmpty, ((L)&&(L)[0])?L:gsEmpty, ((L)&&(L)[0])? (((V)&&(V)[0])?gsEqual:gsSpace):gsEmpty, ((L)&&(L)[0])?((V)&&(V)[0]?V:gsMissing):gsEmpty*/
#define SDF_LBL_VAL(L,V)  ((L)&&(L)[0])?gsSpace:gsEmpty, ((L)&&(L)[0])?L:gsEmpty, ((L)&&(L)[0])? (((V)&&(V)[0])?gsEqual:gsSpace):gsEmpty, ((V)&&(V)[0])?V:((L)&&(L)[0])?gsMissing:gsEmpty

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __STRUTIL_H__ */

