/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.00
 * April 13, 2005
 * Developed at NIST
 */

#ifndef __STRUTIL_H__
#define __STRUTIL_H__


#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

int ExtractConnectedComponent(  inp_ATOM *at, int num_at, int component_number, inp_ATOM *component_at );
int SetConnectedComponentNumber( inp_ATOM *at, int num_at, int component_number );
INChI     *Alloc_INChI( inp_ATOM *at, int num_at, int *found_num_bonds, int *found_num_isotopic, int nAllocMode );
int       Free_INChI(INChI **ppINChI);
INChI_Aux *Alloc_INChI_Aux( int num_at, int num_isotopic_atoms, int nAllocMode, int bOrigData );
int       Free_INChI_Aux( INChI_Aux **ppINChI_Aux );
int  Create_INChI( INChI **ppINChI, INChI_Aux **ppINChI_Aux, ORIG_ATOM_DATA *orig_inp_data,
                  inp_ATOM *inp_at, INP_ATOM_DATA *inp_norm_data[2],
                  int num_inp_at, INCHI_MODE nUserMode,
                  INCHI_MODE *pbTautFlags, INCHI_MODE *pbTautFlagsDone,
                  struct tagInchiTime *ulMaxTime, char *pStrErrStruct);
int FillOutInfAtom(inp_ATOM *norm_at, INF_ATOM_DATA *inf_norm_at_data, int init_num_at,
                   int num_removed_H, int nNumRemovedProtons, NUM_H *nNumRemovedProtonsIsotopic, int bIsotopic,
                   INChI *pINChI, INChI_Aux *pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode );
int FillOutCompositeCanonInfAtom(COMP_ATOM_DATA *composite_norm_data, INF_ATOM_DATA *inf_norm_at_data,
                                 int bIsotopic, int bTautomeric,
                                 PINChI2 *pINChI2, PINChI_Aux2 *pINChI_Aux2, int bAbcNumbers, INCHI_MODE nMode);

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
                         int *bXmlStructStarted, int num_inp, INPUT_PARMS *ip, char *pStr, int nStrLen );

int bNumHeterAtomHasIsotopicH( inp_ATOM *atom, int num_atoms );

int WriteToSDfile( const INP_ATOM_DATA *inp_at_data, INCHI_FILE* fcb, const char* name, const char* comment,
                   const char *szLabel, const char *szValue );
int WriteOrigAtomDataToSDfile( const ORIG_ATOM_DATA *inp_at_data, INCHI_FILE* fcb, const char* name, const char* comment,
                   int bChiral, const char *szLabel, const char *szValue);


extern char gsMissing[];
extern char gsEmpty[];
extern char gsSpace[];
extern char gsEqual[];
/* format string for SDF_LBL_VAL(L,V): %s%s%s%s (four strings) */
/*#define SDF_LBL_VAL(L,V)  ((L)&&(L)[0])?gsSpace:gsEmpty, ((L)&&(L)[0])?L:gsEmpty, ((L)&&(L)[0])? (((V)&&(V)[0])?gsEqual:gsSpace):gsEmpty, ((L)&&(L)[0])?((V)&&(V)[0]?V:gsMissing):gsEmpty*/
#define SDF_LBL_VAL(L,V)  ((L)&&(L)[0])?gsSpace:gsEmpty, ((L)&&(L)[0])?L:gsEmpty, ((L)&&(L)[0])? (((V)&&(V)[0])?gsEqual:gsSpace):gsEmpty, ((V)&&(V)[0])?V:((L)&&(L)[0])?gsMissing:gsEmpty

#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif /* __STRUTIL_H__ */

